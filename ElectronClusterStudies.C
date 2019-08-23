#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "THStack.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TColor.h"
#include <vector>

using namespace std;

float DeltaPhiDir(float Phi1,float Phi2){
  float deltaphi=Phi1-Phi2;
  if(deltaphi>M_PI){
    deltaphi=deltaphi-2*M_PI;
  }
  if(deltaphi<(-M_PI)){
    deltaphi=2*M_PI+deltaphi;
  }
  return deltaphi;
}

float DeltaPhi(float Phi1,float Phi2){
  float deltaphi=fabs(Phi1-Phi2);
  if(deltaphi>M_PI){
    deltaphi=2*M_PI-deltaphi;
  }
  return deltaphi;
}

void CalculatePerformance(const TH1F *const pTH1F, float &resolution, float &resolutionError)
{
  static const float FLOAT_MAX(std::numeric_limits<float>::max());
  
  if (NULL == pTH1F)
    return;
  
  if (5 > pTH1F->GetEntries())
    {
      std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
      return;
    }
  
  // Calculate raw properties of distribution
  float sum = 0., total = 0.;
  double sx = 0., sxx = 0.;
  const unsigned int nbins(pTH1F->GetNbinsX());
  
  for (unsigned int i = 0; i <= nbins; ++i)
    {
      const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
      const float yi(pTH1F->GetBinContent(i));
      sx += yi * binx;
      sxx += yi * binx * binx;
      total += yi;
    }
  
  const float rawMean(sx / total);
  const float rawMeanSqua(sxx / total);
  const float rawRms(std::sqrt(rawMeanSqua - rawMean * rawMean));
  
  sum = 0.;
  unsigned int is0 = 0;
  
  for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
      sum += pTH1F->GetBinContent(i);
      is0 = i;
    }
  
  // Calculate truncated properties
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX), rms(FLOAT_MAX);
  float high(0.f);
  
  for (unsigned int istart = 0; istart <= is0; ++istart)
    {
      double sumn = 0.;
      double csum = 0.;
      double sumx = 0.;
      double sumxx = 0.;
      unsigned int iend = 0;
      
      for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
	  const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
	  const float yi(pTH1F->GetBinContent(i));
	  csum += yi;
	  
	  if (sumn < 0.9 * total)
            {
	      sumn += yi;
	      sumx += yi * binx;
	      sumxx+= yi * binx * binx;
	      iend = i;
            }
        }
      
        const float localMean(sumx / sumn);
        const float localMeanSqua(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSqua - localMean * localMean));
	
        if (localRms < rmsmin)
	  {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;
	    frac = rms / mean * std::sqrt(2) * 100.;
	    efrac = frac / std::sqrt(total);
	  }
	resolution = frac;
	resolutionError = efrac;
    }
  std::cout<<"resolution/error "<<resolution<<"/"<<resolutionError<<" mean "<<std::endl;
}


void CalculatePerformanceNoMean(const TH1F *const pTH1F, float &resolution, float &resolutionError)
{
  static const float FLOAT_MAX(std::numeric_limits<float>::max());
  
  if (NULL == pTH1F)
    return;
  
  if (5 > pTH1F->GetEntries())
    {
      std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
      return;
    }
  
  // Calculate raw properties of distribution
  float sum = 0., total = 0.;
  double sx = 0., sxx = 0.;
  const unsigned int nbins(pTH1F->GetNbinsX());
  
  for (unsigned int i = 0; i <= nbins; ++i)
    {
      const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
      const float yi(pTH1F->GetBinContent(i));
      sx += yi * binx;
      sxx += yi * binx * binx;
      total += yi;
    }
  
  const float rawMean(sx / total);
  const float rawMeanSqua(sxx / total);
  const float rawRms(std::sqrt(rawMeanSqua - rawMean * rawMean));

  sum = 0.;
  unsigned int is0 = 0;
  
  for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
      sum += pTH1F->GetBinContent(i);
      is0 = i;
    }
  
  // Calculate truncated properties
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX), rms(FLOAT_MAX);
  float high(0.f);
  
  for (unsigned int istart = 0; istart <= is0; ++istart)
    {
      double sumn = 0.;
      double csum = 0.;
      double sumx = 0.;
      double sumxx = 0.;
      unsigned int iend = 0;
      
      for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
	  const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
	  const float yi(pTH1F->GetBinContent(i));
	  csum += yi;
	  
	  if (sumn < 0.9 * total)
            {
	      sumn += yi;
	      sumx += yi * binx;
	      sumxx+= yi * binx * binx;
	      iend = i;
            }
        }
      
        const float localMean(sumx / sumn);
        const float localMeanSqua(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSqua - localMean * localMean));
	
        if (localRms < rmsmin)
	  {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;
	    frac = rms* 100.;
	    efrac = frac / std::sqrt(total);
	  }
	resolution = frac;
	resolutionError = efrac;
    }

}



TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void CLICdpStyle()
{
  gROOT->SetStyle("Plain"); /*Default white background for all plots*/
  /* set bkg color of all to kWhite: white, but not 0*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetStatColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(kWhite);
  
  
   /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
   gStyle->SetPaperSize(20, 26); 
   /* No yellow border around histogram*/
   gStyle->SetDrawBorder(0);
   /* remove border of canvas*/
   gStyle->SetCanvasBorderMode(0);
   /* remove border of pads*/
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetLegendBorderSize(0);
  
   /* default text size*/
   gStyle->SetTextSize(0.05);
   gStyle->SetTitleSize(0.06,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   gStyle->SetLabelOffset(0.015,"xyz");
   gStyle->SetTitleOffset(1.2,"yz"); //equivalent to: gStyle->SetTitleYOffset(1.2);
   gStyle->SetTitleOffset(1.17,"x");



   /* Use visible font for all text*/
   int font = 42; 
   gStyle->SetTitleFont(font);
   gStyle->SetTitleFontSize(0.06);
   gStyle->SetStatFont(font);
   gStyle->SetStatFontSize(0.07);
   gStyle->SetTextFont(font);
   gStyle->SetLabelFont(font,"xyz");
   gStyle->SetTitleFont(font,"xyz");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetStatBorderSize(1);
   //ROSA
   //gStyle->SetLegendFont(font);

   /* big marker points*/
   gStyle->SetMarkerStyle(1);
   gStyle->SetLineWidth(2);  
   gStyle->SetMarkerSize(1.2);
   /*set palette in 2d histogram to nice and colorful one*/
   gStyle->SetPalette(1,0); 

   /*No title for histograms*/
   gStyle->SetOptTitle(0);
   /* show the errors on the stat box */
   gStyle->SetOptStat(0); 
   /* show errors on fitted parameters*/
   gStyle->SetOptFit(0); 
   /* number of decimals used for errors*/
   gStyle->SetEndErrorSize(5);   

   /* set line width to 2 by default so that histograms are visible when printed small
      idea: emphasize the data, not the frame around*/
   gStyle->SetHistLineWidth(2);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetFuncWidth(2);
   gStyle->SetHistLineColor(kBlack);
   gStyle->SetFuncColor(kBlack);
   gStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadTopMargin(0.11);
   gStyle->SetPadRightMargin(0.08);
   gStyle->SetPadLeftMargin(0.17);
   
   //set the number of divisions to show
   gStyle->SetNdivisions(506, "xy");
   
   //turn off xy grids
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   
   //set the tick mark style
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   gStyle->SetCanvasDefW(800);
   gStyle->SetCanvasDefH(700);

   gROOT->ForceStyle();
}


TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,0,50,800,300);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fill_singleParticle_histograms(TFile* file, std::vector<TH1F*> h_hist_vec, std::vector<TH2F*> h2_hist_vec, int signal_particle_ID){

  TTree* tree = (TTree*)file->Get("showerData");
  

  vector<int> *track_nHits=0;
  vector<int> *track_ndf=0;
  vector<float> *track_chi2=0;
  vector<float> *track_sigmaPOverP=0;

  vector<float> *track_x_atCalo=0;
  vector<float> *track_y_atCalo=0;
  vector<float> *track_z_atCalo=0;

  vector<float> *track_px_atCalo=0;
  vector<float> *track_py_atCalo=0;
  vector<float> *track_pz_atCalo=0;
  vector<float> *track_p_atCalo=0;

  vector<float> *track_p_atIP=0;
  vector<float> *track_pt_atIP=0;
  vector<float> *track_phi_atIP=0;
  vector<float> *track_theta_atIP=0;

  vector<float> *track_p_innermostHit=0;
  vector<float> *track_pt_innermostHit=0;
  vector<float> *track_x_innermostHit=0;
  vector<float> *track_y_innermostHit=0;
  vector<float> *track_z_innermostHit=0;

  vector<float> *track_p_outermostRHit=0;
  vector<float> *track_pt_outermostRHit=0;
  vector<float> *track_x_outermostRHit=0;
  vector<float> *track_y_outermostRHit=0;
  vector<float> *track_z_outermostRHit=0;

  vector<float> *cluster_energy=0;

  vector<float> *cluster_hit_x=0;
  vector<float> *cluster_hit_y=0;
  vector<float> *cluster_hit_z=0;
  vector<float> *cluster_hit_E=0;
  vector<int> *cluster_hit_type=0;
  vector<int> *cluster_hit_index=0;

  vector<float> *no_cluster_hit_E=0;
  vector<float> *no_cluster_hit_x=0;
  vector<float> *no_cluster_hit_y=0;
  vector<float> *no_cluster_hit_z=0;

  vector<float> *reco_Px=0;
  vector<float> *reco_Py=0;
  vector<float> *reco_Pz=0;
  vector<float> *reco_E=0;
  vector<float> *reco_CosTheta=0;
  vector<int> *reco_PDGID=0;
 vector<float> *reco_clusters_energy=0;

  vector<float> *true_x=0;
  vector<float> *true_y=0;
  vector<float> *true_z=0;

  vector<float> *true_Px=0;
  vector<float> *true_Py=0;
  vector<float> *true_Pz=0;
  vector<float> *true_E=0;
  vector<float> *true_CosTheta=0;
  vector<int> *true_PDGID=0;
  vector<int> *true_GenStatus=0; //0 for irradiated photons


  tree->SetBranchAddress("true_Px", &true_Px);
  tree->SetBranchAddress("true_Py", &true_Py);
  tree->SetBranchAddress("true_Pz", &true_Pz);
  tree->SetBranchAddress("true_Energy", &true_E);
  tree->SetBranchAddress("true_CosTheta", &true_CosTheta);
  tree->SetBranchAddress("true_PDGID", &true_PDGID);
  tree->SetBranchAddress("true_GenStatus", &true_GenStatus);

  tree->SetBranchAddress("reco_Px", &reco_Px);
  tree->SetBranchAddress("reco_Py", &reco_Py); 
  tree->SetBranchAddress("reco_Pz", &reco_Pz);
  tree->SetBranchAddress("reco_Energy", &reco_E);
  tree->SetBranchAddress("reco_CosTheta", &reco_CosTheta);
  tree->SetBranchAddress("reco_PDGID", &reco_PDGID);
  tree->SetBranchAddress("reco_clusters_energy", &reco_clusters_energy);

  tree->SetBranchAddress("track_nHits", &track_nHits);
  tree->SetBranchAddress("track_ndf", &track_ndf);
  tree->SetBranchAddress("track_chi2",&track_chi2);
  tree->SetBranchAddress("track_sigmaPOverP", &track_sigmaPOverP);

  tree->SetBranchAddress("track_x_atCalo", &track_x_atCalo);
  tree->SetBranchAddress("track_y_atCalo", &track_y_atCalo);
  tree->SetBranchAddress("track_z_atCalo", &track_z_atCalo);

  tree->SetBranchAddress("track_px_atCalo", &track_px_atCalo);
  tree->SetBranchAddress("track_py_atCalo", &track_py_atCalo);
  tree->SetBranchAddress("track_pz_atCalo", &track_pz_atCalo);
  tree->SetBranchAddress("track_p_atCalo", &track_p_atCalo);

  tree->SetBranchAddress("track_pt_atIP", &track_pt_atIP);
  tree->SetBranchAddress("track_Phi_atIP", &track_phi_atIP);
  tree->SetBranchAddress("track_Theta_atIP", &track_theta_atIP);
  tree->SetBranchAddress("track_p_atIP", &track_p_atIP);

  tree->SetBranchAddress("track_x_innermostHit",&track_x_innermostHit);
  tree->SetBranchAddress("track_y_innermostHit",&track_y_innermostHit);
  tree->SetBranchAddress("track_z_innermostHit",&track_z_innermostHit);

  tree->SetBranchAddress("track_x_outermostRHit",&track_x_outermostRHit);
  tree->SetBranchAddress("track_y_outermostRHit",&track_y_outermostRHit);
  tree->SetBranchAddress("track_z_outermostRHit",&track_z_outermostRHit);

  tree->SetBranchAddress("track_p_innermostHit",&track_p_innermostHit);
  tree->SetBranchAddress("track_pt_innermostHit",&track_pt_innermostHit);
  tree->SetBranchAddress("track_p_outermostRHit",&track_p_outermostRHit);
  tree->SetBranchAddress("track_pt_outermostRHit",&track_pt_outermostRHit);

  tree->SetBranchAddress("cluster_energy", &cluster_energy);

  tree->SetBranchAddress("cluster_hit_x", &cluster_hit_x);
  tree->SetBranchAddress("cluster_hit_y", &cluster_hit_y);
  tree->SetBranchAddress("cluster_hit_z", &cluster_hit_z);
  tree->SetBranchAddress("cluster_hit_index", &cluster_hit_index);
  tree->SetBranchAddress("cluster_hit_type", &cluster_hit_type);
  tree->SetBranchAddress("cluster_hit_E", &cluster_hit_E);

  tree->SetBranchAddress("no_cluster_hit_x", &no_cluster_hit_x);
  tree->SetBranchAddress("no_cluster_hit_y", &no_cluster_hit_y);
  tree->SetBranchAddress("no_cluster_hit_z", &no_cluster_hit_z);

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    tree->GetEntry(i_entry);
    if(i_entry%10000==0){
      std::cout<<"in entry "<<i_entry<<std::endl;
    }
    bool found_true_quark=true;
    TLorentzVector trueV;

    float true_ph_rad_largest=0;
    float true_ph_rad_all=0;

    for(unsigned int i=0;i<true_Px->size();i++){
      if((*true_GenStatus)[i]==1 && abs((*true_PDGID)[i])==signal_particle_ID){
	trueV.SetPxPyPzE((*true_Px)[i],(*true_Py)[i],(*true_Pz)[i],(*true_E)[i]);
	//break;
      }
      if((*true_GenStatus)[i]==0 && abs((*true_PDGID)[i])==22){
	true_ph_rad_all+=(*true_E)[i];
	if((*true_E)[i]>true_ph_rad_largest){
	  true_ph_rad_largest=(*true_E)[i];
	}
      }
    }
    if(trueV.CosTheta()>0.70 || trueV.CosTheta()<0.20 ){
      //std::cout<<"event "<<i_entry<<" should be skipped "<<std::endl;
      //continue;
    }

    int index_lead_signal=-1;
    float E_reco_lead_signal=-1;
    float E_sum_all_ph=0;
    for(unsigned int i=0;i<reco_Px->size();i++){
      if(abs((*reco_PDGID)[i])==signal_particle_ID){
	if((*reco_E)[i]> E_reco_lead_signal){
	  E_reco_lead_signal=(*reco_E)[i];
	  index_lead_signal=i;
	}
      }
      if(abs((*reco_PDGID)[i])==22){
	E_sum_all_ph+=(*reco_E)[i];
      }
    }
    
    float E_sum_all_ph_1_deg_veto=0;
    float E_sum_all_ph_1_deg_veto_10deg=0;
    float E_sum_all_ph_1_deg_veto_5deg=0;
    float E_sum_all_ph_1_deg_veto_2_5deg=0;
    float E_sum_all_ph_0_5_deg_veto_10deg=0;
    float E_sum_all_ph_0_5_deg_veto_5deg=0;
    float E_sum_all_ph_0_5_deg_veto_2_5deg=0;
    float E_sum_all_ph_0_2_deg_veto_10deg=0;
    float E_sum_all_ph_0_2_deg_veto_5deg=0;
    float E_sum_all_ph_0_2_deg_veto_2_5deg=0;
    float E_sum_all_ph_0_1_deg_veto_10deg=0;
    float E_sum_all_ph_0_1_deg_veto_5deg=0;
    float E_sum_all_ph_0_1_deg_veto_2_5deg=0;
    float E_sum_all_NH_0_1_deg_veto_10deg=0;
    float E_sum_all_NH_0_1_deg_veto_5deg=0;
    float E_sum_all_NH_0_1_deg_veto_2_5deg=0;

    float E_sum_all_ph_no_veto_10deg=0;
    float E_sum_all_ph_no_veto_5deg=0;
    float E_sum_all_ph_no_veto_2_5deg=0;
    float E_sum_all_NH_no_veto_10deg=0;
    float E_sum_all_NH_no_veto_5deg=0;
    float E_sum_all_NH_no_veto_2_5deg=0;

    float E_sum_all_El_no_veto_10deg=0;
    float E_sum_all_El_no_veto_5deg=0;
    
    for(unsigned int i=0;i<reco_Px->size();i++){
      if(abs((*reco_PDGID)[i])==signal_particle_ID){
	if((*reco_E)[i]> E_reco_lead_signal){
	  E_reco_lead_signal=(*reco_E)[i];
	  index_lead_signal=i;
	}
      }
      if(abs((*reco_PDGID)[i])==22){
	E_sum_all_ph+=(*reco_E)[i];
      }
    }

    TLorentzVector recoV;
 
    if(E_reco_lead_signal>0){
      recoV.SetPxPyPzE((*reco_Px)[index_lead_signal],(*reco_Py)[index_lead_signal],(*reco_Pz)[index_lead_signal],(*reco_E)[index_lead_signal]);

      int index_track1=-1;
      int index_track0=-1;

      int index_track0_min=-1;

      float track_p_min=400;
      float track_pt_min=400;
      float track_phi_min=400;

      float pt_track1=-1;

      if(track_p_atIP->size()>1){
	for(unsigned int i=0;i<track_p_atIP->size();i++){
	  //std::cout<<" track "<<i <<" pt/p/N/chi2ndf/sigmapoverp "<<(*track_pt_atIP)[i]<<"/"<<(*track_p_atIP)[i]<<"/"<<(*track_nHits)[i]<<"/"<<(*track_chi2)[i]/(float)(*track_ndf)[i]<<"/"<<(*track_sigmaPOverP)[i]<<std::endl;
	  if( (fabs((*track_p_atIP)[i]-recoV.P())/recoV.P())<track_p_min && (fabs((*track_pt_atIP)[i]-recoV.Pt())/recoV.Pt())<track_pt_min && fabs((*track_phi_atIP)[i]-recoV.Phi())<track_phi_min){//proven to be sufficient to match to the track in question
	    track_pt_min=fabs((*track_pt_atIP)[i]-recoV.Pt())/recoV.Pt();
	    track_p_min=fabs((*track_p_atIP)[i]-recoV.P())/recoV.P();
	    track_phi_min=fabs((*track_phi_atIP)[i]-recoV.Phi());
	    if(index_track0>-1){
	      if((*track_pt_atIP)[i]>pt_track1){
		pt_track1=(*track_pt_atIP)[i];
		index_track1=index_track0;
	      }
	    }
	    index_track0=i;
	    //std::cout<<"px/py/pz of track "<<i<< " and reco are the same "<<fabs((*track_p_atIP)[i]-recoV.P())/recoV.P()<<" "<< fabs((*track_pt_atIP)[i]-recoV.Pt())/recoV.Pt()<<" "<<fabs((*track_phi_atIP)[i]-recoV.Phi())<<std::endl;
	  }else{
	    if((*track_pt_atIP)[i]>pt_track1){
	      pt_track1=(*track_pt_atIP)[i];
	      index_track1=i;
	    }
	    //std::cout<<"far off track px/py/pz of track "<<i<< " and reco are the same "<<fabs((*track_p_atIP)[i]-recoV.P())/recoV.P()<<" "<< fabs((*track_pt_atIP)[i]-recoV.Pt())/recoV.Pt()<<" "<<fabs((*track_phi_atIP)[i]-recoV.Phi())<<std::endl;
	  }
	}
      }else{
	index_track0=0;
      }
      if(track_p_atIP->size()>1){
	index_track1=1;
      }
      //std::cout<<"end of track "<<index_track0<<" p/pt/phi min "<<track_p_min<<"/"<<track_pt_min<<"/"<<track_phi_min<<std::endl;
      if(index_track0==-1){
	std::cout<<"WTF no track found"<<std::endl;
      }
      for(unsigned int i=0;i<reco_Px->size();i++){
	if(i==index_lead_signal){
	  continue;
	}
	if(abs((*reco_PDGID)[i])==11){
	  TLorentzVector tmp_recoEl;
	  tmp_recoEl.SetPxPyPzE((*reco_Px)[i],(*reco_Py)[i],(*reco_Pz)[i],(*reco_E)[i]);
	  
	  if((tmp_recoEl.Angle(recoV.Vect())*TMath::RadToDeg())<10.0){
	    E_sum_all_El_no_veto_10deg+=(*reco_E)[i];
	    if((tmp_recoEl.Angle(recoV.Vect())*TMath::RadToDeg())<5.0){
	      E_sum_all_El_no_veto_5deg+=(*reco_E)[i];
	    }
	  }/*
	  if((fabs(tmp_recoEl.Theta()-recoV.Theta())*TMath::RadToDeg())>0.2){
	    continue;
	  }
	  //cut on dphi angle
	  if((DeltaPhi(tmp_recoEl.Phi(),recoV.Phi())*TMath::RadToDeg())<10.0){
	    E_sum_all_El_no_veto_10deg+=(*reco_E)[i];
	    if((DeltaPhi(tmp_recoEl.Phi(),recoV.Phi())*TMath::RadToDeg())<5.0){
	      E_sum_all_El_no_veto_5deg+=(*reco_E)[i];
	    }
	    }*/
	}

	if(abs((*reco_PDGID)[i])==2112){
	  TLorentzVector tmp_recoNH;
	  tmp_recoNH.SetPxPyPzE((*reco_Px)[i],(*reco_Py)[i],(*reco_Pz)[i],(*reco_E)[i]);
	  
	  if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())<10.0){
	    E_sum_all_NH_no_veto_10deg+=(*reco_E)[i];
	    if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
	      E_sum_all_NH_0_1_deg_veto_10deg+=(*reco_E)[i];
	    }
	    if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())<5.0){
	      E_sum_all_NH_no_veto_5deg+=(*reco_E)[i];
	      if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
		E_sum_all_NH_0_1_deg_veto_5deg+=(*reco_E)[i];
	      }
	      if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())<2.5){
	        E_sum_all_NH_no_veto_2_5deg+=(*reco_E)[i];
		if((tmp_recoNH.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
		  E_sum_all_NH_0_1_deg_veto_2_5deg+=(*reco_E)[i];
		}
	      }
	    }
	  }/*
	  if((fabs(tmp_recoNH.Theta()-recoV.Theta())*TMath::RadToDeg())>0.2){
	    continue;
	  }
	  //cut on dphi angle
	  if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())<10.0){
	    E_sum_all_NH_no_veto_10deg+=(*reco_E)[i];
	    if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
	      E_sum_all_NH_0_1_deg_veto_10deg+=(*reco_E)[i];
	    }
	    if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())<5.0){
	      E_sum_all_NH_no_veto_5deg+=(*reco_E)[i];
	      if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
		E_sum_all_NH_0_1_deg_veto_5deg+=(*reco_E)[i];
	      }
	      if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())<2.5){
		E_sum_all_NH_no_veto_2_5deg+=(*reco_E)[i];
		if((DeltaPhi(tmp_recoNH.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
		  E_sum_all_NH_0_1_deg_veto_2_5deg+=(*reco_E)[i];
		}
	    }
	    }
	    }*/
	}
	if(abs((*reco_PDGID)[i])==22){
	  TLorentzVector tmp_reco;
	  tmp_reco.SetPxPyPzE((*reco_Px)[i],(*reco_Py)[i],(*reco_Pz)[i],(*reco_E)[i]);
	  //cut on 3D angle
	  
	  if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>1.0){
	    E_sum_all_ph_1_deg_veto=(*reco_E)[i];
	  }
	  if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())<10.0){
	    E_sum_all_ph_no_veto_10deg+=(*reco_E)[i];
	    if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
	      E_sum_all_ph_0_1_deg_veto_10deg+=(*reco_E)[i];
	      if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.2){
		E_sum_all_ph_0_2_deg_veto_10deg+=(*reco_E)[i];
		if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.5){
		  E_sum_all_ph_0_5_deg_veto_10deg+=(*reco_E)[i];
		  if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>1.0){
		    E_sum_all_ph_1_deg_veto_10deg+=(*reco_E)[i];
		  }
		}
	      }
	    }

	    if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())<5.0){
	      E_sum_all_ph_no_veto_5deg+=(*reco_E)[i];
	      if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
		E_sum_all_ph_0_1_deg_veto_5deg+=(*reco_E)[i];
		if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.2){
		  E_sum_all_ph_0_2_deg_veto_5deg+=(*reco_E)[i];
		  if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.5){
		    E_sum_all_ph_0_5_deg_veto_5deg+=(*reco_E)[i];
		    if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>1.0){
		      E_sum_all_ph_1_deg_veto_5deg+=(*reco_E)[i];
		    }
		  }
		}
	      }
	      if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())<2.5){
	        E_sum_all_ph_no_veto_2_5deg+=(*reco_E)[i];
		if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.01){
		  E_sum_all_ph_0_1_deg_veto_2_5deg+=(*reco_E)[i];
		  if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.2){
		    E_sum_all_ph_0_2_deg_veto_2_5deg+=(*reco_E)[i];
		    if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>0.5){
		      E_sum_all_ph_0_5_deg_veto_2_5deg+=(*reco_E)[i];
		      if((tmp_reco.Angle(recoV.Vect())*TMath::RadToDeg())>1.0){
			E_sum_all_ph_1_deg_veto_2_5deg+=(*reco_E)[i];
		      }
		    }
		  }
		}
	      }
	    }
	  }//cut on 3D angle
	  /*
	  if((fabs(tmp_reco.Theta()-recoV.Theta())*TMath::RadToDeg())>0.2){
	    continue;
	  }
	  //cut on dphi angle
	  if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>1.0){
	    E_sum_all_ph_1_deg_veto=(*reco_E)[i];
	  }
	  if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())<10.0){
	    E_sum_all_ph_no_veto_10deg+=(*reco_E)[i];
	    if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
	      E_sum_all_ph_0_1_deg_veto_10deg+=(*reco_E)[i];
	      if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.2){
		E_sum_all_ph_0_2_deg_veto_10deg+=(*reco_E)[i];
		if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.5){
		  E_sum_all_ph_0_5_deg_veto_10deg+=(*reco_E)[i];
		  if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>1.0){
		    E_sum_all_ph_1_deg_veto_10deg+=(*reco_E)[i];
		  }
		}
	      }
	    }

	    if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())<5.0){
	      E_sum_all_ph_no_veto_5deg+=(*reco_E)[i];
	      if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
		E_sum_all_ph_0_1_deg_veto_5deg+=(*reco_E)[i];
		if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.2){
		  E_sum_all_ph_0_2_deg_veto_5deg+=(*reco_E)[i];
		  if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.5){
		    E_sum_all_ph_0_5_deg_veto_5deg+=(*reco_E)[i];
		    if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>1.0){
		      E_sum_all_ph_1_deg_veto_5deg+=(*reco_E)[i];
		    }
		  }
		}
	      }
	      if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())<2.5){
		E_sum_all_ph_no_veto_2_5deg+=(*reco_E)[i];
		if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.01){
		  E_sum_all_ph_0_1_deg_veto_2_5deg+=(*reco_E)[i];
		  if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.2){
		    E_sum_all_ph_0_2_deg_veto_2_5deg+=(*reco_E)[i];
		    if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>0.5){
		      E_sum_all_ph_0_5_deg_veto_2_5deg+=(*reco_E)[i];
		      if((DeltaPhi(tmp_reco.Phi(),recoV.Phi())*TMath::RadToDeg())>1.0){
			E_sum_all_ph_1_deg_veto_2_5deg+=(*reco_E)[i];
		      }
		    }
		  }
		}
	      }
	    }
	  }//cut on dphi
	  */
	}
      }

      //std::cout<<"hist_vec size"<<h_hist_vec.size()<<std::endl;

      h_hist_vec[0]->Fill(recoV.E()/trueV.E());
      h_hist_vec[1]->Fill((recoV.E()+E_sum_all_ph)/trueV.E());
      if((true_ph_rad_largest/trueV.E())>0.10){
	h_hist_vec[12]->Fill(recoV.E()/trueV.E());
	h_hist_vec[13]->Fill((recoV.E()+E_sum_all_ph)/trueV.E());
      }else{
	h_hist_vec[14]->Fill(recoV.E()/trueV.E());
	h_hist_vec[15]->Fill((recoV.E()+E_sum_all_ph)/trueV.E());
      }
      if(true_ph_rad_all==0){
	h_hist_vec[16]->Fill(recoV.E()/trueV.E());
	h_hist_vec[17]->Fill((recoV.E()+E_sum_all_ph)/trueV.E());
      }
      if((true_ph_rad_all/trueV.E())<0.10){
	h_hist_vec[18]->Fill(recoV.E()/trueV.E());
	h_hist_vec[19]->Fill((recoV.E()+E_sum_all_ph)/trueV.E());

	h_hist_vec[20]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto)/trueV.E());
	h_hist_vec[22]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_10deg)/trueV.E());
	h_hist_vec[24]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[26]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[28]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_10deg)/trueV.E());
	h_hist_vec[30]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_5deg)/trueV.E());
	h_hist_vec[32]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[34]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_10deg)/trueV.E());
	h_hist_vec[36]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_5deg)/trueV.E());
	h_hist_vec[38]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[40]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_10deg)/trueV.E());
	h_hist_vec[42]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[44]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_2_5deg)/trueV.E());

	h_hist_vec[71]->Fill(track_pt_atIP->size());
	h_hist_vec[72]->Fill((*track_nHits)[index_track0]);
	if(fabs(recoV.CosTheta())<0.7){
	  h_hist_vec[73]->Fill(sqrt(pow((*track_x_outermostRHit)[index_track0],2)+pow((*track_y_outermostRHit)[index_track0],2)));
	  h_hist_vec[74]->Fill(sqrt(pow((*track_x_innermostHit)[index_track0],2)+pow((*track_y_innermostHit)[index_track0],2)));
	}else if(fabs(recoV.CosTheta())>0.85){
	  h_hist_vec[75]->Fill(fabs((*track_z_outermostRHit)[index_track0]));
	  h_hist_vec[76]->Fill(fabs((*track_z_innermostHit)[index_track0]));
	}
	h_hist_vec[77]->Fill((*track_pt_outermostRHit)[index_track0]/(*track_pt_innermostHit)[index_track0]);
	h_hist_vec[78]->Fill((*track_p_outermostRHit)[index_track0]/(*track_p_innermostHit)[index_track0]);
	h_hist_vec[79]->Fill((*track_sigmaPOverP)[index_track0]);
	h_hist_vec[80]->Fill((*track_chi2)[index_track0]/(float)(*track_ndf)[index_track0]);
	h_hist_vec[81]->Fill((*track_p_atIP)[index_track0]/recoV.E());
	h_hist_vec[82]->Fill((*track_p_atIP)[index_track0]/(*reco_clusters_energy)[index_lead_signal]);
	h_hist_vec[96]->Fill((*reco_clusters_energy)[index_lead_signal]/trueV.E());
	if( index_track1 >-1){
	  h_hist_vec[83]->Fill((*track_nHits)[index_track1]);
	  if(fabs(recoV.CosTheta())<0.7){
	    h_hist_vec[84]->Fill(sqrt(pow((*track_x_outermostRHit)[index_track1],2)+pow((*track_y_outermostRHit)[index_track1],2)));
	    h_hist_vec[85]->Fill(sqrt(pow((*track_x_innermostHit)[index_track1],2)+pow((*track_y_innermostHit)[index_track1],2)));
	  }else if(fabs(recoV.CosTheta())>0.85){
	    h_hist_vec[86]->Fill(fabs((*track_z_outermostRHit)[index_track1]));
	    h_hist_vec[87]->Fill(fabs((*track_z_innermostHit)[index_track1]));
	  }
	  h_hist_vec[88]->Fill((*track_pt_outermostRHit)[index_track1]/(*track_pt_innermostHit)[index_track1]);
	  h_hist_vec[89]->Fill((*track_p_outermostRHit)[index_track1]/(*track_p_innermostHit)[index_track1]);
	  h_hist_vec[90]->Fill((*track_sigmaPOverP)[index_track1]);
	  h_hist_vec[91]->Fill((*track_chi2)[index_track1]/(float)(*track_ndf)[index_track1]);
	  h_hist_vec[92]->Fill((*track_p_atIP)[index_track1]/recoV.E());
	  h_hist_vec[93]->Fill((*track_p_atIP)[index_track1]/(*reco_clusters_energy)[index_lead_signal]);
	  h_hist_vec[94]->Fill((*track_p_atIP)[index_track1]/(*track_p_atIP)[index_track0]);
	  //if(track_p_atIP->size()==2){
	  //std::cout<<" a lot distance here "<<sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_innermostHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_innermostHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_innermostHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_outermostRHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_outermostRHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_outermostRHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2))<<std::endl;
	  //}
	  h_hist_vec[95]->Fill(sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2)));
	  //}
	}
	h_hist_vec[98]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[100]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());
	h_hist_vec[102]->Fill((recoV.E()+E_sum_all_ph_no_veto_2_5deg+E_sum_all_NH_no_veto_2_5deg)/trueV.E());

	h_hist_vec[104]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_10deg+E_sum_all_NH_0_1_deg_veto_10deg)/trueV.E());
 	h_hist_vec[106]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_5deg+E_sum_all_NH_0_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[108]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_2_5deg+E_sum_all_NH_0_1_deg_veto_2_5deg)/trueV.E());

  	h_hist_vec[110]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg)/trueV.E());
 	h_hist_vec[112]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg)/trueV.E());
 
	h_hist_vec[114]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[116]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

	h_hist_vec[118]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg)/trueV.E());
 	h_hist_vec[120]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg)/trueV.E());

	h_hist_vec[122]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[124]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

	h_hist_vec[126]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg)/trueV.E());
 	h_hist_vec[128]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg)/trueV.E());

	h_hist_vec[130]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[132]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

      }else{
	//std::cout<<"high RI/RO/ZI/ZO/sPOverP/pOverClu/POverE "<<sqrt(pow((*track_x_innermostHit)[index_track0],2)+pow((*track_y_innermostHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_outermostRHit)[index_track0],2)+pow((*track_y_outermostRHit)[index_track0],2))<<"/"<<(*track_z_innermostHit)[index_track0]<<"/"<<(*track_z_outermostRHit)[index_track0]<<"/"<<(*track_sigmaPOverP)[index_track0]<<"/"<<(*track_p_atIP)[index_track0]/(*reco_clusters_energy)[index_lead_signal]<<"/"<<(*track_p_atIP)[index_track0]/recoV.E()<<std::endl;
	h_hist_vec[21]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto)/trueV.E());
	h_hist_vec[23]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_10deg)/trueV.E());
	h_hist_vec[25]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[27]->Fill((recoV.E()+E_sum_all_ph_1_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[29]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_10deg)/trueV.E());
	h_hist_vec[31]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_5deg)/trueV.E());
	h_hist_vec[33]->Fill((recoV.E()+E_sum_all_ph_0_5_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[35]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_10deg)/trueV.E());
	h_hist_vec[37]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_5deg)/trueV.E());
	h_hist_vec[39]->Fill((recoV.E()+E_sum_all_ph_0_2_deg_veto_2_5deg)/trueV.E());
	h_hist_vec[41]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_10deg)/trueV.E());
	h_hist_vec[43]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[45]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_2_5deg)/trueV.E());

	h_hist_vec[46]->Fill(track_pt_atIP->size());
	h_hist_vec[47]->Fill((*track_nHits)[index_track0]);
	//std::cout<<"RI/RO/ZI/ZO/sPOverP/pOverClu/POverE "<<sqrt(pow((*track_x_innermostHit)[index_track0],2)+pow((*track_y_innermostHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_outermostRHit)[index_track0],2)+pow((*track_y_outermostRHit)[index_track0],2))<<"/"<<(*track_z_innermostHit)[index_track0]<<"/"<<(*track_z_outermostRHit)[index_track0]<<"/"<<(*track_sigmaPOverP)[index_track0]<<"/"<<(*track_p_atIP)[index_track0]/(*reco_clusters_energy)[index_lead_signal]<<"/"<<(*track_p_atIP)[index_track0]/recoV.E()<<std::endl;
	if(fabs(recoV.CosTheta())<0.7){
	  h_hist_vec[48]->Fill(sqrt(pow((*track_x_outermostRHit)[index_track0],2)+pow((*track_y_outermostRHit)[index_track0],2)));
	  h_hist_vec[49]->Fill(sqrt(pow((*track_x_innermostHit)[index_track0],2)+pow((*track_y_innermostHit)[index_track0],2)));
	}else if(fabs(recoV.CosTheta())>0.85){
	  h_hist_vec[50]->Fill(fabs((*track_z_outermostRHit)[index_track0]));
	  h_hist_vec[51]->Fill(fabs((*track_z_innermostHit)[index_track0]));
	}
	h_hist_vec[52]->Fill((*track_pt_outermostRHit)[index_track0]/(*track_pt_innermostHit)[index_track0]);
	h_hist_vec[53]->Fill((*track_p_outermostRHit)[index_track0]/(*track_p_innermostHit)[index_track0]);
	h_hist_vec[54]->Fill((*track_sigmaPOverP)[index_track0]);
	h_hist_vec[55]->Fill((*track_chi2)[index_track0]/(float)(*track_ndf)[index_track0]);
	h_hist_vec[56]->Fill((*track_p_atIP)[index_track0]/recoV.E());
	h_hist_vec[57]->Fill((*track_p_atIP)[index_track0]/(*reco_clusters_energy)[index_lead_signal]);
	h_hist_vec[97]->Fill((*reco_clusters_energy)[index_lead_signal]/trueV.E());
	if( index_track1 >-1){
	  h_hist_vec[58]->Fill((*track_nHits)[index_track1]);
	  if(fabs(recoV.CosTheta())<0.7){
	    h_hist_vec[59]->Fill(sqrt(pow((*track_x_outermostRHit)[index_track1],2)+pow((*track_y_outermostRHit)[index_track1],2)));
	    h_hist_vec[60]->Fill(sqrt(pow((*track_x_innermostHit)[index_track1],2)+pow((*track_y_innermostHit)[index_track1],2)));
	  }else if(fabs(recoV.CosTheta())>0.85){
	    h_hist_vec[61]->Fill(fabs((*track_z_outermostRHit)[index_track1]));
	    h_hist_vec[62]->Fill(fabs((*track_z_innermostHit)[index_track1]));
	  }
	  h_hist_vec[63]->Fill((*track_pt_outermostRHit)[index_track1]/(*track_pt_innermostHit)[index_track1]);
	  h_hist_vec[64]->Fill((*track_p_outermostRHit)[index_track1]/(*track_p_innermostHit)[index_track1]);
	  h_hist_vec[65]->Fill((*track_sigmaPOverP)[index_track1]);
	  h_hist_vec[66]->Fill((*track_chi2)[index_track1]/(float)(*track_ndf)[index_track1]);
	  h_hist_vec[67]->Fill((*track_p_atIP)[index_track1]/recoV.E());
	  h_hist_vec[68]->Fill((*track_p_atIP)[index_track1]/(*reco_clusters_energy)[index_lead_signal]);
	  h_hist_vec[69]->Fill((*track_p_atIP)[index_track1]/(*track_p_atIP)[index_track0]);
	  h_hist_vec[70]->Fill(sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2)));
	  //if(track_p_atIP->size()==2){
	  //std::cout<<" few try distance here "<<sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_innermostHit)[index_track1]-(*track_x_innermostHit)[index_track0],2) + pow((*track_y_innermostHit)[index_track1]-(*track_y_innermostHit)[index_track0],2) + pow((*track_z_innermostHit)[index_track1]-(*track_z_innermostHit)[index_track0],2))<<"/"<<sqrt(pow((*track_x_outermostRHit)[index_track1]-(*track_x_outermostRHit)[index_track0],2) + pow((*track_y_outermostRHit)[index_track1]-(*track_y_outermostRHit)[index_track0],2) + pow((*track_z_outermostRHit)[index_track1]-(*track_z_outermostRHit)[index_track0],2))<<std::endl;

	}
	h_hist_vec[99]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[101]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());
 	h_hist_vec[103]->Fill((recoV.E()+E_sum_all_ph_no_veto_2_5deg+E_sum_all_NH_no_veto_2_5deg)/trueV.E());

	h_hist_vec[105]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_10deg+E_sum_all_NH_0_1_deg_veto_10deg)/trueV.E());
 	h_hist_vec[107]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_5deg+E_sum_all_NH_0_1_deg_veto_5deg)/trueV.E());
	h_hist_vec[109]->Fill((recoV.E()+E_sum_all_ph_0_1_deg_veto_2_5deg+E_sum_all_NH_0_1_deg_veto_2_5deg)/trueV.E());

  	h_hist_vec[111]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg)/trueV.E());
 	h_hist_vec[113]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg)/trueV.E());
 
	h_hist_vec[115]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[117]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

	h_hist_vec[119]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg)/trueV.E());
 	h_hist_vec[121]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg)/trueV.E());

	h_hist_vec[123]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[125]->Fill(((*reco_clusters_energy)[index_lead_signal]+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

	h_hist_vec[127]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg)/trueV.E());
 	h_hist_vec[129]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg)/trueV.E());

	h_hist_vec[131]->Fill((recoV.E()+E_sum_all_ph_no_veto_10deg+E_sum_all_El_no_veto_10deg+E_sum_all_NH_no_veto_10deg)/trueV.E());
 	h_hist_vec[133]->Fill((recoV.E()+E_sum_all_ph_no_veto_5deg+E_sum_all_El_no_veto_5deg+E_sum_all_NH_no_veto_5deg)/trueV.E());

      }
    }

    int ind_clus_0=-1;
    int ind_clus_1=-1;
    float clu0_E=-1;
    float clu1_E=-1;
    
    //check for leading and subleading cluster index
    for(unsigned int cl=0; cl<cluster_energy->size();cl++){
      if((*cluster_energy)[cl]>clu0_E){
	clu1_E=clu0_E;
	ind_clus_1=ind_clus_0;
	clu0_E=(*cluster_energy)[cl];
	ind_clus_0=cl;
      }else if((*cluster_energy)[cl]>clu1_E){
	clu1_E=(*cluster_energy)[cl];
	ind_clus_1=cl;
      }
    }
    if(E_reco_lead_signal>0){
      if(clu0_E>0){
	h_hist_vec[2]->Fill(recoV.E()+clu0_E);
      }else{
	h_hist_vec[2]->Fill(recoV.E());
      }
    }

    if(track_x_atCalo->size()==1){
      //TVectorLorentz trackstateV;
      //trackstateV.SetPtEtaPhiM((*track_px_atCalo)[0],(*track_py_atCalo)[0],(*track_pz_atCalo)[0]);
      //calculate the deltaR between the track state at CALO and the photon candidates

      int ind_recs_0=-1;
      int ind_recs_1=-1;
      float rec0_E=-1;
      float rec1_E=-1;
      for(unsigned int cl=0; cl<reco_E->size();cl++){
	if(abs((*reco_PDGID)[cl])!=22){
	  continue;
	}
	if((*reco_E)[cl]>rec0_E){
	  rec1_E=rec0_E;
	  ind_recs_1=ind_recs_0;
	  rec0_E=(*reco_E)[cl];
	  ind_recs_0=cl;
	}else if((*reco_E)[cl]>rec1_E){
	  rec1_E=(*reco_E)[cl];
	  ind_recs_1=cl;
	}
      }
      //cosTheta=pz/p

      float theta_trackstate=acos((*track_pz_atCalo)[0]/(*track_p_atCalo)[0]);
      /*
      if(ind_recs_0>-1){
	h2_hist_vec[10]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_0],(*reco_Px)[ind_recs_0]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_0])-theta_trackstate)*TMath::RadToDeg());
	h2_hist_vec[11]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_0],(*reco_Px)[ind_recs_0]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_0])-theta_trackstate)*TMath::RadToDeg(),rec0_E);
	h2_hist_vec[16]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_0],(*reco_Px)[ind_recs_0]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_0])-recoV.Theta())*TMath::RadToDeg());
	h2_hist_vec[17]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_0],(*reco_Px)[ind_recs_0]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_0])-recoV.Theta())*TMath::RadToDeg(),(*reco_E)[ind_recs_0]);
      }
      if(ind_recs_1>-1){
	h2_hist_vec[12]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_1],(*reco_Px)[ind_recs_1]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_1])-theta_trackstate)*TMath::RadToDeg());
	h2_hist_vec[13]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_1],(*reco_Px)[ind_recs_1]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_1])-theta_trackstate)*TMath::RadToDeg(),rec1_E);
	h2_hist_vec[18]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_1],(*reco_Px)[ind_recs_1]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_1])-recoV.Theta())*TMath::RadToDeg());
	h2_hist_vec[19]->Fill(DeltaPhiDir(atan2((*reco_Py)[ind_recs_1],(*reco_Px)[ind_recs_1]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[ind_recs_1])-recoV.Theta())*TMath::RadToDeg(),(*reco_E)[ind_recs_1]);
      }


      for(unsigned int i=0;i<reco_Px->size();i++){
	if(abs((*reco_PDGID)[i])==22){
	h2_hist_vec[8]->Fill(DeltaPhiDir(atan2((*reco_Py)[i],(*reco_Px)[i]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[i])-theta_trackstate)*TMath::RadToDeg());
	h2_hist_vec[9]->Fill(DeltaPhiDir(atan2((*reco_Py)[i],(*reco_Px)[i]),atan2((*track_py_atCalo)[0],(*track_px_atCalo)[0]))*TMath::RadToDeg(), (acos((*reco_CosTheta)[i])-theta_trackstate)*TMath::RadToDeg(),(*reco_E)[i]);

	h2_hist_vec[14]->Fill(DeltaPhiDir(atan2((*reco_Py)[i],(*reco_Px)[i]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[i])-recoV.Theta())*TMath::RadToDeg());
	h2_hist_vec[15]->Fill(DeltaPhiDir(atan2((*reco_Py)[i],(*reco_Px)[i]),recoV.Phi())*TMath::RadToDeg(), (acos((*reco_CosTheta)[i])-recoV.Theta())*TMath::RadToDeg(),(*reco_E)[i]);
	}
	}*/
    }
    
    if(track_x_atCalo->size()==1){
      //exactly one track in the event
      /*
      for(unsigned int c=0;c<cluster_hit_x->size();c++){
	//TH2 vec [0]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for ALL clusters,i.e. all cluster hits --> only ECAL hits
	//TH2 vec [1]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for ALL clusters,i.e. all cluster hits, weight entries by the hit energies --> check for ECAL hits only
	//TH2 vec [2]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for cluster 0,i.e. all cluster hits --> only ECAL hits
	//TH2 vec [3]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for cluster 0,i.e. all cluster hits, energy weighted --> only ECAL hits
	//TH2 vec [4]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for most energetic cluster,i.e. all cluster hits --> only ECAL hits
	//TH2 vec [5]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for most energetic cluster,i.e. all cluster hits, energy weighted --> only ECAL hits
	//TH2 vec [6]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for second energetic cluster ,i.e. all cluster hits --> only ECAL hits
	//TH2 vec [7]//2D plot r-z plane corrected for track state impact point at CaloFace -->check for second energetic cluster ,i.e. all cluster hits, energy weighted --> only ECAL hits
	if((*cluster_hit_type)[c]==1){//consider ONLY ECAL hits
	  h2_hist_vec[0]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0]);
	  h2_hist_vec[1]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0],(*cluster_hit_E)[c]);
	  if((*cluster_hit_index)[c]==0){
	    h2_hist_vec[2]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0]);
	    h2_hist_vec[3]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0],(*cluster_hit_E)[c]);
	  }
	  if((*cluster_hit_index)[c]==ind_clus_0){
	    h2_hist_vec[4]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0]);
	    h2_hist_vec[5]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0],(*cluster_hit_E)[c]);
	  }
	  if((*cluster_hit_index)[c]==ind_clus_1){
	    h2_hist_vec[6]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0]);
	    h2_hist_vec[7]->Fill(sqrt(pow((*cluster_hit_x)[c]-(*track_x_atCalo)[0],2)+ pow((*cluster_hit_y)[c]-(*track_y_atCalo)[0],2)), (*cluster_hit_z)[c]-(*track_z_atCalo)[0],(*cluster_hit_E)[c]);
	  }
	}
	}*/
    }
  }
}


void ElectronClusterPlotsFull(){

  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

  //const char* final_histo_name="/eos/user/w/weberma2/data/validation171212/electron_cluster_shape_em10_em100_em1000_dress_ph_el_nh_noClusterHits_dphiSelectionOnly_dtheta_max_0_2degrees_lowest_dphiVeto_0_01_wTracksInfo.root";
  const char* final_histo_name="/eos/user/w/weberma2/data/validation171212/electron_cluster_shape_em10_em100_em1000_dress_ph_el_nh_noClusterHits_dAngleSelection_lowest_angleVeto_0_01_wTracksInfo.root";

  int signal_particleID=11;
  bool is_electron=false;
  if(is_electron){
    signal_particleID=11;
  }

  string label_legend= "e^{-}, E_{true}=10 GeV";

  unsigned int n_bins100=100;
  float lim_E_rel_low=0.00;
  float lim_E_rel_high=2.00;

 
  TFile* file_em10_pt1=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em10_fullClusterTrackInfo_noHits_Full_pt1.root");
  TFile* file_em10_pt2=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em10_fullClusterTrackInfo_noHits_Full_pt2.root");
  TFile* file_em10_pt3=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em10_fullClusterTrackInfo_noHits_Full_pt3.root");
  TFile* file_em10_pt4=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em10_fullClusterTrackInfo_noHits_Full_pt4.root");

  TFile* file_em100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em100_fullClusterTrackInfo_noHits_Full.root");
  
  TFile* file_em1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em1000_fullClusterTrackInfo_noHits_Full.root");

  //TFile* file_em10_2=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_test_electronoutput_em10_fullClusterHitTrackInfo_pt_AllHits_pt5.root");

  TFile* file_histogram=new TFile(final_histo_name,"recreate");



  TH1F* h_em100_reco_el_over_gen_el = new TH1F("h_em100_reco_el_over_gen_el","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_all_ph = new TH1F("h_em100_reco_el_over_gen_el_dress_all_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_with_em100clu0 = new TH1F("h_em100_reco_el_over_gen_el_dress_with_em100clu0","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_with_em100clu0_DR50 = new TH1F("h_em100_reco_el_over_gen_el_dress_with_em100clu0_DR50","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_track_pt_last_track_pt_first = new TH1F("h_em100_reco_el_track_pt_last_track_pt_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_track_pt_atCalo_track_pt_atIP = new TH1F("h_em100_reco_el_track_pt_atCalo_track_pt_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_track_p_last_track_p_first = new TH1F("h_em100_reco_el_track_p_last_track_p_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_track_p_atCalo_track_p_atIP = new TH1F("h_em100_reco_el_track_p_atCalo_track_p_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_no_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_no_true_ph = new TH1F("h_em100_reco_el_over_gen_el_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph = new TH1F("h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_nTracks_true_ph_sum_0_10_E_orig = new TH1F("h_em100_nTracks_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em100_TrackNHits_true_ph_sum_0_10_E_orig = new TH1F("h_em100_TrackNHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em100_TrackOuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_TrackOuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em100_TrackInnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_TrackInnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em100_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_TrackOuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_TrackInnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em100_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em100_sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em100_Track_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em100_chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em100_Track_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em100_Track1NHits_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1NHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em100_Track1OuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_Track1OuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em100_Track1InnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_Track1InnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,1500);
  TH1F* h_em100_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_Track1OuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_Track1InnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em100_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em100_Track1_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);


  TH1F* h_em100_nTracks_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_nTracks_no_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em100_TrackNHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_TrackNHits_no_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em100_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_TrackOuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em100_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_TrackInnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em100_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_TrackOuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_TrackInnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em100_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em100_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em100_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em100_Track1NHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1NHits_no_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em100_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_Track1OuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em100_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_barrel_Track1InnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,1500);
  TH1F* h_em100_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_Track1OuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_endcap_Track1InnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em100_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em100_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em100_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em100_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em100_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);

  TH1F* h_em100_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
 
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_em100;
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_all_ph);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_with_em100clu0);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_with_em100clu0_DR50);
  hist_vec_em100.push_back(h_em100_reco_el_track_pt_last_track_pt_first);
  hist_vec_em100.push_back(h_em100_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_track_pt_atCalo_track_pt_atIP);
  hist_vec_em100.push_back(h_em100_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_track_p_last_track_p_first);
  hist_vec_em100.push_back(h_em100_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_track_p_atCalo_track_p_atIP);
  hist_vec_em100.push_back(h_em100_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_no_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_no_true_ph);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_nTracks_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackNHits_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackOuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackInnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_p_OverECluster_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_Track1NHits_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1OuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1InnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_OverECluster_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_nTracks_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackNHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_Track1NHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
 
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig);
  
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
 
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig); 
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig); 
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
 
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig); 
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
 
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig);
 
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em100.push_back(h_em100_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig);


  for(unsigned int i=0;i<hist_vec_em100.size();i++){
    hist_vec_em100[i]->Sumw2();
    hist_vec_em100[i]->SetLineWidth(2);
    hist_vec_em100[i]->SetLineColor(kBlack);
    hist_vec_em100[i]->SetLineStyle(2);
  }
  /*
  TH2F* h2_em100_cluster_deltaR_vs_deltaZ=new TH2F("h2_em100_cluster_deltaR_vs_deltaZ", "h2_em100_cluster_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_cluster_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_cluster_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em100_cluster_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em100_cluster_deltaR_vs_deltaZ_EWeight", "h2_em100_cluster_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_cluster_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_cluster_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em100_cluster0_deltaR_vs_deltaZ=new TH2F("h2_em100_cluster0_deltaR_vs_deltaZ", "h2_em100_cluster0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_cluster0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_cluster0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em100_cluster0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em100_cluster0_deltaR_vs_deltaZ_EWeight", "h2_em100_cluster0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_cluster0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_cluster0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em100_clusterE0_deltaR_vs_deltaZ=new TH2F("h2_em100_clusterE0_deltaR_vs_deltaZ", "h2_em100_clusterE0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_clusterE0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_clusterE0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em100_clusterE0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em100_clusterE0_deltaR_vs_deltaZ_EWeight", "h2_em100_clusterE0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_clusterE0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_clusterE0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em100_clusterE1_deltaR_vs_deltaZ=new TH2F("h2_em100_clusterE1_deltaR_vs_deltaZ", "h2_em100_clusterE1_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_clusterE1_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_clusterE1_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em100_clusterE1_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em100_clusterE1_deltaR_vs_deltaZ_EWeight", "h2_em100_clusterE1_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em100_clusterE1_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em100_clusterE1_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  float delta_phi_min = -50;
  float delta_phi_max = 30;

  float delta_theta_min = -60;
  float delta_theta_max = 60;


  int nbins_angle2D=20;

  TH2F* h2_em100_ph_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em100_ph_deltaPhi_vs_deltaTheta_track", "h2_em100_ph_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_ph_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_ph_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em100_phE0_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em100_phE0_deltaPhi_vs_deltaTheta_track", "h2_em100_phE0_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_phE0_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_phE0_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em100_phE1_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em100_phE1_deltaPhi_vs_deltaTheta_track", "h2_em100_phE1_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_phE1_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_phE1_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  float delta_phi_min_reco = -4;
  float delta_phi_max_reco = 15;

  float delta_theta_min_reco = -1;
  float delta_theta_max_reco = 1;


  TH2F* h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el", "h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el", "h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el", "h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  */
  
  std::vector<TH2F*> hist2D_vec_em100;
  /*
  hist2D_vec_em100.push_back(h2_em100_cluster_deltaR_vs_deltaZ);
  hist2D_vec_em100.push_back(h2_em100_cluster_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em100.push_back(h2_em100_cluster0_deltaR_vs_deltaZ);
  hist2D_vec_em100.push_back(h2_em100_cluster0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em100.push_back(h2_em100_clusterE0_deltaR_vs_deltaZ);
  hist2D_vec_em100.push_back(h2_em100_clusterE0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em100.push_back(h2_em100_clusterE1_deltaR_vs_deltaZ);
  hist2D_vec_em100.push_back(h2_em100_clusterE1_deltaR_vs_deltaZ_EWeight);

  hist2D_vec_em100.push_back(h2_em100_ph_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em100.push_back(h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em100.push_back(h2_em100_phE0_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em100.push_back(h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em100.push_back(h2_em100_phE1_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em100.push_back(h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_track);

  hist2D_vec_em100.push_back(h2_em100_ph_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em100.push_back(h2_em100_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em100.push_back(h2_em100_phE0_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em100.push_back(h2_em100_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em100.push_back(h2_em100_phE1_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em100.push_back(h2_em100_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  */
  for(unsigned int i=0;i<hist2D_vec_em100.size();i++){
    hist2D_vec_em100[i]->Sumw2();
    //hist2D_vec_em100[i]->
  }

  fill_singleParticle_histograms(file_em100,hist_vec_em100,hist2D_vec_em100,signal_particleID);
  //fill_singleParticle_histograms(file_em100_2,hist_vec_em100,hist2D_vec_em100,signal_particleID);


  TH1F* h_em1000_reco_el_over_gen_el = new TH1F("h_em1000_reco_el_over_gen_el","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_all_ph = new TH1F("h_em1000_reco_el_over_gen_el_dress_all_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_with_em1000clu0 = new TH1F("h_em1000_reco_el_over_gen_el_dress_with_em1000clu0","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_with_em1000clu0_DR50 = new TH1F("h_em1000_reco_el_over_gen_el_dress_with_em1000clu0_DR50","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_track_pt_last_track_pt_first = new TH1F("h_em1000_reco_el_track_pt_last_track_pt_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_track_pt_atCalo_track_pt_atIP = new TH1F("h_em1000_reco_el_track_pt_atCalo_track_pt_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_track_p_last_track_p_first = new TH1F("h_em1000_reco_el_track_p_last_track_p_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_track_p_atCalo_track_p_atIP = new TH1F("h_em1000_reco_el_track_p_atCalo_track_p_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_no_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_no_true_ph = new TH1F("h_em1000_reco_el_over_gen_el_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph = new TH1F("h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_nTracks_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_nTracks_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em1000_TrackNHits_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_TrackNHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em1000_TrackOuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_TrackOuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em1000_TrackInnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_TrackInnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em1000_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_TrackOuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em1000_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_TrackInnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em1000_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em1000_Track_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em1000_Track_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em1000_Track1NHits_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1NHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em1000_Track1OuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_Track1OuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em1000_Track1InnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_Track1InnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em1000_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_Track1OuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em1000_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_Track1InnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em1000_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em1000_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em1000_Track1_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);

  TH1F* h_em1000_nTracks_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_nTracks_no_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em1000_TrackNHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_TrackNHits_no_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em1000_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_TrackOuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em1000_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_TrackInnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em1000_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_TrackOuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em1000_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_TrackInnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em1000_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em1000_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em1000_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em1000_Track1NHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1NHits_no_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em1000_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_Track1OuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em1000_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_barrel_Track1InnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,1500);
  TH1F* h_em1000_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_Track1OuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em1000_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_endcap_Track1InnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em1000_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em1000_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em1000_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em1000_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em1000_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);

  TH1F* h_em1000_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_em1000;
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_all_ph);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_with_em1000clu0);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_with_em1000clu0_DR50);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_pt_last_track_pt_first);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_pt_atCalo_track_pt_atIP);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_p_last_track_p_first);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_p_atCalo_track_p_atIP);
  hist_vec_em1000.push_back(h_em1000_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_no_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_no_true_ph);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_nTracks_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackNHits_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackOuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackInnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_p_OverECluster_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_Track1NHits_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1OuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1InnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_OverECluster_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_nTracks_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackNHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_Track1NHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
 
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig);
  
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);  
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);  
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig);
 
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em1000.push_back(h_em1000_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig);



  for(unsigned int i=0;i<hist_vec_em1000.size();i++){
    hist_vec_em1000[i]->Sumw2();
    hist_vec_em1000[i]->SetLineWidth(2);
    hist_vec_em1000[i]->SetLineColor(kBlack);
    hist_vec_em1000[i]->SetLineStyle(2);
  }
  /*
  TH2F* h2_em1000_cluster_deltaR_vs_deltaZ=new TH2F("h2_em1000_cluster_deltaR_vs_deltaZ", "h2_em1000_cluster_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_cluster_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_cluster_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em1000_cluster_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em1000_cluster_deltaR_vs_deltaZ_EWeight", "h2_em1000_cluster_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_cluster_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_cluster_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em1000_cluster0_deltaR_vs_deltaZ=new TH2F("h2_em1000_cluster0_deltaR_vs_deltaZ", "h2_em1000_cluster0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_cluster0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_cluster0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em1000_cluster0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em1000_cluster0_deltaR_vs_deltaZ_EWeight", "h2_em1000_cluster0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_cluster0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_cluster0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em1000_clusterE0_deltaR_vs_deltaZ=new TH2F("h2_em1000_clusterE0_deltaR_vs_deltaZ", "h2_em1000_clusterE0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_clusterE0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_clusterE0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em1000_clusterE0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em1000_clusterE0_deltaR_vs_deltaZ_EWeight", "h2_em1000_clusterE0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_clusterE0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_clusterE0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em1000_clusterE1_deltaR_vs_deltaZ=new TH2F("h2_em1000_clusterE1_deltaR_vs_deltaZ", "h2_em1000_clusterE1_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_clusterE1_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_clusterE1_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em1000_clusterE1_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em1000_clusterE1_deltaR_vs_deltaZ_EWeight", "h2_em1000_clusterE1_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em1000_clusterE1_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em1000_clusterE1_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em1000_ph_0deltaPhi_vs_deltaTheta_track=new TH2F("h2_em1000_ph_0deltaPhi_vs_deltaTheta_track", "h2_em1000_ph_0deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track", "h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em1000_phE0_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em1000_phE0_deltaPhi_vs_deltaTheta_track", "h2_em1000_phE0_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em1000_phE1_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em1000_phE1_deltaPhi_vs_deltaTheta_track", "h2_em1000_phE1_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");


  TH2F* h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el", "h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el", "h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el", "h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  */
  
  std::vector<TH2F*> hist2D_vec_em1000;
  /*
  hist2D_vec_em1000.push_back(h2_em1000_cluster_deltaR_vs_deltaZ);
  hist2D_vec_em1000.push_back(h2_em1000_cluster_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em1000.push_back(h2_em1000_cluster0_deltaR_vs_deltaZ);
  hist2D_vec_em1000.push_back(h2_em1000_cluster0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em1000.push_back(h2_em1000_clusterE0_deltaR_vs_deltaZ);
  hist2D_vec_em1000.push_back(h2_em1000_clusterE0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em1000.push_back(h2_em1000_clusterE1_deltaR_vs_deltaZ);
  hist2D_vec_em1000.push_back(h2_em1000_clusterE1_deltaR_vs_deltaZ_EWeight);

  hist2D_vec_em1000.push_back(h2_em1000_ph_0deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em1000.push_back(h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em1000.push_back(h2_em1000_phE0_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em1000.push_back(h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em1000.push_back(h2_em1000_phE1_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em1000.push_back(h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_track);

  hist2D_vec_em1000.push_back(h2_em1000_ph_0deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em1000.push_back(h2_em1000_ph_0deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em1000.push_back(h2_em1000_phE0_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em1000.push_back(h2_em1000_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em1000.push_back(h2_em1000_phE1_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em1000.push_back(h2_em1000_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  */
  for(unsigned int i=0;i<hist2D_vec_em1000.size();i++){
    hist2D_vec_em1000[i]->Sumw2();
    //hist2D_vec_em1000[i]->
  }

  fill_singleParticle_histograms(file_em1000,hist_vec_em1000,hist2D_vec_em1000,signal_particleID);
  //fill_singleParticle_histograms(file_em1000_2,hist_vec_em1000,hist2D_vec_em1000,signal_particleID);





  TH1F* h_em10_reco_el_over_gen_el = new TH1F("h_em10_reco_el_over_gen_el","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_all_ph = new TH1F("h_em10_reco_el_over_gen_el_dress_all_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_with_em10clu0 = new TH1F("h_em10_reco_el_over_gen_el_dress_with_em10clu0","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_with_em10clu0_DR50 = new TH1F("h_em10_reco_el_over_gen_el_dress_with_em10clu0_DR50","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_track_pt_last_track_pt_first = new TH1F("h_em10_reco_el_track_pt_last_track_pt_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_track_pt_atCalo_track_pt_atIP = new TH1F("h_em10_reco_el_track_pt_atCalo_track_pt_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_track_pt_atCalo_track_pt_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_track_p_last_track_p_first = new TH1F("h_em10_reco_el_track_p_last_track_p_first","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_track_p_atCalo_track_p_atIP = new TH1F("h_em10_reco_el_track_p_atCalo_track_p_atIP","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_track_p_atCalo_track_p_atIP_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_no_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_no_true_ph = new TH1F("h_em10_reco_el_over_gen_el_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph = new TH1F("h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_nTracks_true_ph_sum_0_10_E_orig = new TH1F("h_em10_nTracks_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em10_TrackNHits_true_ph_sum_0_10_E_orig = new TH1F("h_em10_TrackNHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em10_TrackOuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_TrackOuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em10_TrackInnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_TrackInnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em10_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_TrackOuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em10_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_TrackInnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em10_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em10_sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em10_Track_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em10_chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em10_Track_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em10_Track1NHits_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1NHits_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em10_Track1OuterR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_Track1OuterR_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em10_Track1InnerR_barrel_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_Track1InnerR_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em10_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_Track1OuterZ_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em10_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_Track1InnerZ_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em10_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig = new TH1F("h_em10__Track1sigmaPOverP_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em10_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig = new TH1F("h_em10__Track1chi2OverNDF_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em10_Track1_p_OverEReco_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_OverEReco_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_p_OverECluster_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_OverECluster_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);


  TH1F* h_em10_nTracks_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_nTracks_no_true_ph_sum_0_10_E_orig","",10,-0.5,9.5);
  TH1F* h_em10_TrackNHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_TrackNHits_no_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em10_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_TrackOuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em10_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_TrackInnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,1500);
  TH1F* h_em10_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_TrackOuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em10_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_TrackInnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em10_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em10_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em10_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);

  TH1F* h_em10_Track1NHits_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1NHits_no_true_ph_sum_0_10_E_orig","",20,-0.5,19.5);
  TH1F* h_em10_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_Track1OuterR_no_true_ph_sum_0_10_E_orig","",n_bins100,50,1500);
  TH1F* h_em10_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_barrel_Track1InnerR_no_true_ph_sum_0_10_E_orig","",n_bins100,30,200);
  TH1F* h_em10_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_Track1OuterZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2200);
  TH1F* h_em10_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_endcap_Track1InnerZ_no_true_ph_sum_0_10_E_orig","",n_bins100,0,400);
  TH1F* h_em10_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_ori","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig","",n_bins100,0.75,1.25);
  TH1F* h_em10_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig","",n_bins100,0,0.05);
  TH1F* h_em10_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig","",n_bins100,0,5.0);
  TH1F* h_em10_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig","",n_bins100,0,2.0);
  TH1F* h_em10_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig","",n_bins100,0,3000);

  TH1F* h_em10_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);
  TH1F* h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig = new TH1F("h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig","",n_bins100,lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_em10;
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_all_ph);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_with_em10clu0);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_with_em10clu0_DR50);
  hist_vec_em10.push_back(h_em10_reco_el_track_pt_last_track_pt_first);
  hist_vec_em10.push_back(h_em10_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_track_pt_atCalo_track_pt_atIP);
  hist_vec_em10.push_back(h_em10_reco_el_track_pt_last_track_pt_first_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_track_p_last_track_p_first);
  hist_vec_em10.push_back(h_em10_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_track_p_atCalo_track_p_atIP);
  hist_vec_em10.push_back(h_em10_reco_el_track_p_last_track_p_first_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_all_ph_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_no_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_no_true_ph);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_all_ph_no_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_5deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_2deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);
 
  hist_vec_em10.push_back(h_em10_nTracks_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackNHits_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackOuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackInnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackOuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackInnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_p_OverECluster_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_Track1NHits_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1OuterR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1InnerR_barrel_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1OuterZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1InnerZ_endcap_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_ptLastOverFirst_Hit_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_pLastOverFirstHit_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_sigmaPOverP_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_chi2OverNDF_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_OverEReco_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_OverECluster_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_Over_Track0_p_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_3D_dist_1st_track0_last_hit_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_nTracks_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackNHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackOuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackInnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackOuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_TrackInnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track_p_OverECluster_no_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_Track1NHits_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1OuterR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1InnerR_barrel_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1OuterZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1InnerZ_endcap_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_ptLastOverFirst_Hit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_pLastOverFirstHit_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_sigmaPOverP_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_chi2OverNDF_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_OverEReco_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_OverECluster_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_p_Over_Track0_p_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_Track1_3D_dist_1st_track0_last_hit_no_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_clustE_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clustE_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
 
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_no_veto_2_5deg_true_ph_sum_0_10_E_orig);
  
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_5deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wNH_0_1deg_veto_2_5deg_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);  
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);  
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wNH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wNH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);  
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_NH_10deg_over_gen_el_E_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_clust_E_dress_ph_wEl_NH_5deg_over_gen_el_E_true_ph_sum_0_10_E_orig);

  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_no_veto_5deg_true_ph_sum_0_10_E_orig);
 
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_no_true_ph_sum_0_10_E_orig );
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_10deg_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_no_true_ph_sum_0_10_E_orig);
  hist_vec_em10.push_back(h_em10_reco_el_over_gen_el_dress_wEl_NH_no_veto_5deg_true_ph_sum_0_10_E_orig);

  for(unsigned int i=0;i<hist_vec_em10.size();i++){
    hist_vec_em10[i]->Sumw2();
    hist_vec_em10[i]->SetLineWidth(2);
    hist_vec_em10[i]->SetLineColor(kBlack);
    hist_vec_em10[i]->SetLineStyle(2);
  }
  /*
  TH2F* h2_em10_cluster_deltaR_vs_deltaZ=new TH2F("h2_em10_cluster_deltaR_vs_deltaZ", "h2_em10_cluster_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_cluster_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_cluster_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em10_cluster_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em10_cluster_deltaR_vs_deltaZ_EWeight", "h2_em10_cluster_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_cluster_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_cluster_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em10_cluster0_deltaR_vs_deltaZ=new TH2F("h2_em10_cluster0_deltaR_vs_deltaZ", "h2_em10_cluster0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_cluster0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_cluster0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em10_cluster0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em10_cluster0_deltaR_vs_deltaZ_EWeight", "h2_em10_cluster0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_cluster0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_cluster0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em10_clusterE0_deltaR_vs_deltaZ=new TH2F("h2_em10_clusterE0_deltaR_vs_deltaZ", "h2_em10_clusterE0_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_clusterE0_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_clusterE0_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em10_clusterE0_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em10_clusterE0_deltaR_vs_deltaZ_EWeight", "h2_em10_clusterE0_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_clusterE0_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_clusterE0_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em10_clusterE1_deltaR_vs_deltaZ=new TH2F("h2_em10_clusterE1_deltaR_vs_deltaZ", "h2_em10_clusterE1_deltaR vs deltaZ",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_clusterE1_deltaR_vs_deltaZ->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_clusterE1_deltaR_vs_deltaZ->SetYTitle("cluster_hit_z-trackAtCalo_z");
  TH2F* h2_em10_clusterE1_deltaR_vs_deltaZ_EWeight=new TH2F("h2_em10_clusterE1_deltaR_vs_deltaZ_EWeight", "h2_em10_clusterE1_deltaR vs deltaZ_EWeight",n_bins100,-10.,150.,n_bins100,-100.,100.);
  h2_em10_clusterE1_deltaR_vs_deltaZ_EWeight->SetXTitle("cluster_hit_R-trackAtCalo_R");
  h2_em10_clusterE1_deltaR_vs_deltaZ_EWeight->SetYTitle("cluster_hit_z-trackAtCalo_z");

  TH2F* h2_em10_ph_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em10_ph_deltaPhi_vs_deltaTheta_track", "h2_em10_ph_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_ph_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_ph_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em10_phE0_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em10_phE0_deltaPhi_vs_deltaTheta_track", "h2_em10_phE0_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_phE0_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_phE0_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");

  TH2F* h2_em10_phE1_deltaPhi_vs_deltaTheta_track=new TH2F("h2_em10_phE1_deltaPhi_vs_deltaTheta_track", "h2_em10_phE1_deltaPhi_vs_deltaTheta_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_phE1_deltaPhi_vs_deltaTheta_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_phE1_deltaPhi_vs_deltaTheta_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");
  TH2F* h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track=new TH2F("h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track", "h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track",nbins_angle2D,delta_phi_min,delta_phi_max,nbins_angle2D,delta_theta_min,delta_theta_max);
  h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetXTitle("reco_ph_phi-trackAtCalo_phi");
  h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track->SetYTitle("reco_ph_theta-trackAtCalo_theta");


  TH2F* h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el", "h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el", "h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");

  TH2F* h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el=new TH2F("h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el", "h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  TH2F* h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el=new TH2F("h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el", "h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el",nbins_angle2D,delta_phi_min_reco,delta_phi_max_reco,nbins_angle2D,delta_theta_min_reco,delta_theta_max_reco);
  h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetXTitle("reco_ph_phi-reco_el_phi");
  h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el->SetYTitle("reco_ph_theta-reco_el_theta");
  */
  
  std::vector<TH2F*> hist2D_vec_em10;
  /*
  hist2D_vec_em10.push_back(h2_em10_cluster_deltaR_vs_deltaZ);
  hist2D_vec_em10.push_back(h2_em10_cluster_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em10.push_back(h2_em10_cluster0_deltaR_vs_deltaZ);
  hist2D_vec_em10.push_back(h2_em10_cluster0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em10.push_back(h2_em10_clusterE0_deltaR_vs_deltaZ);
  hist2D_vec_em10.push_back(h2_em10_clusterE0_deltaR_vs_deltaZ_EWeight);
  hist2D_vec_em10.push_back(h2_em10_clusterE1_deltaR_vs_deltaZ);
  hist2D_vec_em10.push_back(h2_em10_clusterE1_deltaR_vs_deltaZ_EWeight);

  hist2D_vec_em10.push_back(h2_em10_ph_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em10.push_back(h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em10.push_back(h2_em10_phE0_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em10.push_back(h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_track);
  hist2D_vec_em10.push_back(h2_em10_phE1_deltaPhi_vs_deltaTheta_track);
  hist2D_vec_em10.push_back(h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_track);

  hist2D_vec_em10.push_back(h2_em10_ph_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em10.push_back(h2_em10_ph_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em10.push_back(h2_em10_phE0_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em10.push_back(h2_em10_phE0_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  hist2D_vec_em10.push_back(h2_em10_phE1_deltaPhi_vs_deltaTheta_reco_el);
  hist2D_vec_em10.push_back(h2_em10_phE1_deltaPhi_vs_deltaTheta_EWeight_reco_el);
  */
  for(unsigned int i=0;i<hist2D_vec_em10.size();i++){
    hist2D_vec_em10[i]->Sumw2();
    //hist2D_vec_em10[i]->
  }

  fill_singleParticle_histograms(file_em10_pt1,hist_vec_em10,hist2D_vec_em10,signal_particleID);
  fill_singleParticle_histograms(file_em10_pt2,hist_vec_em10,hist2D_vec_em10,signal_particleID);
  fill_singleParticle_histograms(file_em10_pt3,hist_vec_em10,hist2D_vec_em10,signal_particleID);
  fill_singleParticle_histograms(file_em10_pt4,hist_vec_em10,hist2D_vec_em10,signal_particleID);


  file_histogram->Write();
  file_histogram->Close();

}
