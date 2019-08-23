#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1F.h"
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
#include "TProfile.h"
#include "TColor.h"
#include <vector>

using namespace std;


double fnc_dscb(double*xx,double*pp)
{
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];
  
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

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

void CalculatePerformance(const TH1F *const pTH1F, float &resolution, float &resolutionError, float &mean, float &meanError)
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
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX)/*, mean(FLOAT_MAX)*/, low(FLOAT_MAX), rms(FLOAT_MAX);
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
	    meanError = mean / std::sqrt(total); 
	  }
	resolution = frac;
	resolutionError = efrac;
    }
  std::cout<<"resolution/error "<<resolution<<"/"<<resolutionError<<" mean "<<std::endl;
}


void CalculatePerformanceNoMean(const TH1F *const pTH1F, float &resolution, float &resolutionError, float &mean, float &meanError)
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
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX)/*, mean(FLOAT_MAX)*/, low(FLOAT_MAX), rms(FLOAT_MAX);
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
	    meanError = mean / std::sqrt(total);
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


void fill_background_missingPT_histograms(TFile* file, std::vector<TH1F*> h_hist_vec){

  TTree* tree = (TTree*)file->Get("showerData");
  
  int eventNumber=0;
  float Z_mcE=0;

  float d1_cosTheta=0;

  float E_totPFO=0;
  float E_totPi=0;
  float E_totPh=0;
  float E_totE=0;
  float E_totMu=0;
  float E_totN=0;

  float px_totPFO=0;
  float px_totPi=0;
  float px_totPh=0;
  float px_totE=0;
  float px_totMu=0;
  float px_totN=0;

  float py_totPFO=0;
  float py_totPi=0;
  float py_totPh=0;
  float py_totE=0;
  float py_totMu=0;
  float py_totN=0;

  float E_trueNeut=0;
  float E_trueZ2=0;

  float px_trueNeut=0;
  float px_trueZ2=0;

  float py_trueNeut=0;
  float py_trueZ2=0;

  float pz_trueNeut=0;
  float pz_trueZ2=0;

  vector<float> *jetExclPx=0;
  vector<float> *jetExclPy=0;
  vector<float> *jetExclPz=0;
  vector<float> *jetExclE=0;
  vector<float> *jetExclCosTheta=0;

  vector<float> *jetInclPx=0;
  vector<float> *jetInclPy=0;
  vector<float> *jetInclPz=0;
  vector<float> *jetInclE=0;
  vector<float> *jetInclCosTheta=0;



  vector<float> *trueMEPx=0;
  vector<float> *trueMEPy=0;
  vector<float> *trueMEPz=0;
  vector<float> *trueMEE=0;
  vector<int> *trueMEPDGID=0;

  tree->SetBranchAddress("trueME_Px", &trueMEPx);
  tree->SetBranchAddress("trueME_Py", &trueMEPy);
  tree->SetBranchAddress("trueME_Pz", &trueMEPz);
  tree->SetBranchAddress("trueME_E", &trueMEE);
  tree->SetBranchAddress("trueME_PDGID", &trueMEPDGID);

  tree->SetBranchAddress("d1_mcCosTheta",&d1_cosTheta);


  tree->SetBranchAddress("eventNumber",&eventNumber);
  tree->SetBranchAddress("Z_mcE",&Z_mcE);

  tree->SetBranchAddress("E_totPFO", &E_totPFO);
  tree->SetBranchAddress("E_totPi", &E_totPi);
  tree->SetBranchAddress("E_totPh", &E_totPh);
  tree->SetBranchAddress("E_totE", &E_totE);
  tree->SetBranchAddress("E_totMu", &E_totMu);
  tree->SetBranchAddress("E_totN", &E_totN);

  tree->SetBranchAddress("px_totPFO", &px_totPFO);
  tree->SetBranchAddress("px_totPi", &px_totPi);
  tree->SetBranchAddress("px_totPh", &px_totPh);
  tree->SetBranchAddress("px_totE", &px_totE);
  tree->SetBranchAddress("px_totMu", &px_totMu);
  tree->SetBranchAddress("px_totN", &px_totN);

  tree->SetBranchAddress("py_totPFO", &py_totPFO);
  tree->SetBranchAddress("py_totPi", &py_totPi);
  tree->SetBranchAddress("py_totPh", &py_totPh);
  tree->SetBranchAddress("py_totE", &py_totE);
  tree->SetBranchAddress("py_totMu", &py_totMu);
  tree->SetBranchAddress("py_totN", &py_totN);

  tree->SetBranchAddress("E_trueNeut", &E_trueNeut);
  tree->SetBranchAddress("E_trueZ2", &E_trueZ2);
  //Z2 is the hadronically decaying Z,  
  tree->SetBranchAddress("px_trueNeut", &px_trueNeut);
  tree->SetBranchAddress("px_trueZ2", &px_trueZ2);

  tree->SetBranchAddress("py_trueNeut", &py_trueNeut);
  tree->SetBranchAddress("py_trueZ2", &py_trueZ2);

  tree->SetBranchAddress("pz_trueNeut", &pz_trueNeut);
  tree->SetBranchAddress("pz_trueZ2", &pz_trueZ2);


  tree->SetBranchAddress("jetExclPx", &jetExclPx);
  tree->SetBranchAddress("jetExclPy", &jetExclPy);
  tree->SetBranchAddress("jetExclPz", &jetExclPz);
  tree->SetBranchAddress("jetExclE", &jetExclE);
  tree->SetBranchAddress("jetExclCosTheta", &jetExclCosTheta);

  tree->SetBranchAddress("jetInclPx", &jetInclPx);
  tree->SetBranchAddress("jetInclPy", &jetInclPy);
  tree->SetBranchAddress("jetInclPz", &jetInclPz);
  tree->SetBranchAddress("jetInclE", &jetInclE);
  tree->SetBranchAddress("jetInclCosTheta", &jetInclCosTheta);

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    tree->GetEntry(i_entry);
    bool forward_quarks=false;
    if(trueMEPDGID!=NULL){
      for(unsigned int i=0;i<trueMEPDGID->size();i++){
	if(abs((*trueMEPDGID)[i])<7){
	  TLorentzVector tmp;
	  tmp.SetPxPyPzE((*trueMEPx)[i],(*trueMEPy)[i],(*trueMEPz)[i],(*trueMEE)[i]);
	  //if( (tmp.Theta()*TMath::RadToDeg())<60. || (tmp.Theta()*TMath::RadToDeg())>135.){
	  if(fabs(tmp.CosTheta())>0.70){
	    forward_quarks=true;
	  }
	}
      }
    }
    if(forward_quarks){
      //std::cout<<"fw quarks "<<std::endl;
      continue;
    }

       //if(fabs(sqrt(E_trueNeut*E_trueNeut-px_trueNeut*px_trueNeut-py_trueNeut*py_trueNeut-pz_trueNeut*pz_trueNeut))>300.){
       //continue;
       //}
    h_hist_vec[0]->Fill(-px_totPFO);
    h_hist_vec[1]->Fill(sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut));
    if(jetExclCosTheta->size()==2){
      if(fabs((*jetExclCosTheta)[0])<0.9 && fabs((*jetExclCosTheta)[1])<0.9){
	h_hist_vec[2]->Fill(sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO));
	h_hist_vec[3]->Fill(E_totPFO);
	h_hist_vec[4]->Fill((*jetExclE)[0]+(*jetExclE)[1]);
	h_hist_vec[5]->Fill(TMath::RadToDeg()*atan2(-py_totPFO,-px_totPFO));
	h_hist_vec[6]->Fill(sqrt(pow((*jetExclPx)[0]+(*jetExclPx)[1],2)+pow((*jetExclPy)[0]+(*jetExclPy)[1],2)));
	h_hist_vec[7]->Fill(TMath::RadToDeg()*DeltaPhi(atan2((*jetExclPy)[0],(*jetExclPx)[0]),atan2((*jetExclPy)[1],(*jetExclPx)[1])));
      }
    }
    h_hist_vec[8] ->Fill(sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut));
    h_hist_vec[9]->Fill(TMath::RadToDeg()*atan2(py_trueNeut,px_trueNeut));
    h_hist_vec[10]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    //0,25,50,75,100,140,180,220,260,340,420,500,750,1500
    float trueMET=sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut);
    //std::cout<<"delta MET/true/reco MET/ dphi "<<sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut)<<"/"<<sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut)<<"/"<<sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)<<"/"<<TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO))<<std::endl;
    if(trueMET<25){
      h_hist_vec[11]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[12]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<50){
      h_hist_vec[13]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[14]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<75){
      h_hist_vec[15]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[16]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<100){
      h_hist_vec[17]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[18]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<140){
      h_hist_vec[19]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[20]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<180){
      h_hist_vec[21]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[22]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<220){
      h_hist_vec[23]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[24]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<260){
      h_hist_vec[25]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[26]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<340){
      h_hist_vec[27]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[28]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<420){
      h_hist_vec[29]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[30]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<500){
      h_hist_vec[31]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[32]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else if(trueMET<750){
      h_hist_vec[33]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[34]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }else{
      h_hist_vec[35]->Fill((sqrt(px_totPFO*px_totPFO+py_totPFO*py_totPFO)  - sqrt(px_trueNeut*px_trueNeut+py_trueNeut*py_trueNeut))/trueMET);
      h_hist_vec[36]->Fill(TMath::RadToDeg()*DeltaPhiDir(atan2(py_trueNeut,px_trueNeut),atan2(-py_totPFO,-px_totPFO)));
    }
  }
  delete tree;
}

  
void fill_resolution_histograms(TFile* file, std::vector<TH1F*> h_hist_vec){

  TTree* tree = (TTree*)file->Get("showerData");
  
  int eventNumber=0;
  float Z_mcE=0;
  
  float d1_cosTheta=0;
  
  float E_totPFO=0;
  float E_totPi=0;
  float E_totPh=0;
  float E_totE=0;
  float E_totMu=0;
  float E_totN=0;

  float px_totPFO=0;
  float px_totPi=0;
  float px_totPh=0;
  float px_totE=0;
  float px_totMu=0;
  float px_totN=0;

  float py_totPFO=0;
  float py_totPi=0;
  float py_totPh=0;
  float py_totE=0;
  float py_totMu=0;
  float py_totN=0;


  float E_true_totAll=0;
  float E_true_totPi=0;
  float E_true_totPh=0;
  float E_true_totE=0;
  float E_true_totMu=0;
  float E_true_totN=0;
  float E_true_totP=0;
  float E_true_totK=0;

  float px_true_totAll=0;
  float px_true_totPi=0;
  float px_true_totPh=0;
  float px_true_totE=0;
  float px_true_totMu=0;
  float px_true_totN=0;
  float px_true_totP=0;
  float px_true_totK=0;

  float py_true_totAll=0;
  float py_true_totPi=0;
  float py_true_totPh=0;
  float py_true_totE=0;
  float py_true_totMu=0;
  float py_true_totN=0;
  float py_true_totP=0;
  float py_true_totK=0;

  vector<float> *jetExclPx=0;
  vector<float> *jetExclPy=0;
  vector<float> *jetExclPz=0;
  vector<float> *jetExclE=0;
  vector<float> *jetExclCosTheta=0;

  vector<float> *jetInclPx=0;
  vector<float> *jetInclPy=0;
  vector<float> *jetInclPz=0;
  vector<float> *jetInclE=0;
  vector<float> *jetInclCosTheta=0;


  tree->SetBranchAddress("jetExclPx", &jetExclPx);
  tree->SetBranchAddress("jetExclPy", &jetExclPy);
  tree->SetBranchAddress("jetExclPz", &jetExclPz);
  tree->SetBranchAddress("jetExclE", &jetExclE);
  tree->SetBranchAddress("jetExclCosTheta", &jetExclCosTheta);

  tree->SetBranchAddress("jetInclPx", &jetInclPx);
  tree->SetBranchAddress("jetInclPy", &jetInclPy);
  tree->SetBranchAddress("jetInclPz", &jetInclPz);
  tree->SetBranchAddress("jetInclE", &jetInclE);
  tree->SetBranchAddress("jetInclCosTheta", &jetInclCosTheta);


  tree->SetBranchAddress("d1_mcCosTheta",&d1_cosTheta);

  tree->SetBranchAddress("eventNumber",&eventNumber);
  tree->SetBranchAddress("Z_mcE",&Z_mcE);

  tree->SetBranchAddress("E_totPFO", &E_totPFO);
  tree->SetBranchAddress("E_totPi", &E_totPi);
  tree->SetBranchAddress("E_totPh", &E_totPh);
  tree->SetBranchAddress("E_totE", &E_totE);
  tree->SetBranchAddress("E_totMu", &E_totMu);
  tree->SetBranchAddress("E_totN", &E_totN);


  tree->SetBranchAddress("px_totPFO", &px_totPFO);
  tree->SetBranchAddress("px_totPi", &px_totPi);
  tree->SetBranchAddress("px_totPh", &px_totPh);
  tree->SetBranchAddress("px_totE", &px_totE);
  tree->SetBranchAddress("px_totMu", &px_totMu);
  tree->SetBranchAddress("px_totN", &px_totN);

  tree->SetBranchAddress("py_totPFO", &py_totPFO);
  tree->SetBranchAddress("py_totPi", &py_totPi);
  tree->SetBranchAddress("py_totPh", &py_totPh);
  tree->SetBranchAddress("py_totE", &py_totE);
  tree->SetBranchAddress("py_totMu", &py_totMu);
  tree->SetBranchAddress("py_totN", &py_totN);

  tree->SetBranchAddress("E_true_totAll", &E_true_totAll);
  tree->SetBranchAddress("E_true_totPi", &E_true_totPi);
  tree->SetBranchAddress("E_true_totPh", &E_true_totPh);
  tree->SetBranchAddress("E_true_totE", &E_true_totE);
  tree->SetBranchAddress("E_true_totMu", &E_true_totMu);
  tree->SetBranchAddress("E_true_totN", &E_true_totN);
  tree->SetBranchAddress("E_true_totK", &E_true_totK);
  tree->SetBranchAddress("E_true_totP", &E_true_totP);

  tree->SetBranchAddress("px_true_totAll", &px_true_totAll);
  tree->SetBranchAddress("px_true_totPi", &px_true_totPi);
  tree->SetBranchAddress("px_true_totPh", &px_true_totPh);
  tree->SetBranchAddress("px_true_totE", &px_true_totE);
  tree->SetBranchAddress("px_true_totMu", &px_true_totMu);
  tree->SetBranchAddress("px_true_totN", &px_true_totN);
  tree->SetBranchAddress("px_true_totK", &px_true_totK);
  tree->SetBranchAddress("px_true_totP", &px_true_totP);

  tree->SetBranchAddress("py_true_totAll", &py_true_totAll);
  tree->SetBranchAddress("py_true_totPi", &py_true_totPi);
  tree->SetBranchAddress("py_true_totPh", &py_true_totPh);
  tree->SetBranchAddress("py_true_totE", &py_true_totE);
  tree->SetBranchAddress("py_true_totMu", &py_true_totMu);
  tree->SetBranchAddress("py_true_totN", &py_true_totN);
  tree->SetBranchAddress("py_true_totK", &py_true_totK);
  tree->SetBranchAddress("py_true_totP", &py_true_totP);

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);
    if(fabs(d1_cosTheta)<0.1){
      h_hist_vec[0]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.2){
      h_hist_vec[1]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.3){
      h_hist_vec[2]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.4){
      h_hist_vec[3]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.5){
      h_hist_vec[4]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.6){
      h_hist_vec[5]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.7){
      h_hist_vec[6]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.8){
      h_hist_vec[7]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.9){
      h_hist_vec[8]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.925){
      h_hist_vec[9]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.950){
      h_hist_vec[10]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.975){
      h_hist_vec[11]->Fill(E_totPFO/E_true_totAll);
    }else if(fabs(d1_cosTheta)<0.985){
      h_hist_vec[12]->Fill(E_totPFO/E_true_totAll);
    }
    if(fabs(d1_cosTheta)<0.7){
      h_hist_vec[13]->Fill(E_totPFO/E_true_totAll);
    }
   if(jetExclCosTheta->size()==2){
     if(fabs((*jetExclCosTheta)[0])<0.7 && fabs((*jetExclCosTheta)[1])<0.7){
       h_hist_vec[14]->Fill(E_totPFO/E_true_totAll);
     }
   }
  }
  delete tree;
}

void fill_B_and_D_histograms(TFile* file, std::vector<TH1F*> h_hist_vec){

  TTree* tree = (TTree*)file->Get("showerData");

  vector<float> *true_Eptx=0;
  vector<float> *true_Epty=0;
  vector<float> *true_Eptz=0;
  vector<float> *true_Eptr=0;
  vector<int> *true_PDGID=0;

  vector<float> *trueTauDaughterVtxx=0;
  vector<float> *trueTauDaughterVtxy=0;
  vector<float> *trueTauDaughterVtxz=0;
  vector<float> *trueTauDaughterVtxr=0;
  vector<int> *trueTauDaughterMotherPDGID=0;



  tree->SetBranchAddress("true_Eptx", &true_Eptx);
  tree->SetBranchAddress("true_Epty", &true_Epty);
  tree->SetBranchAddress("true_Eptz", &true_Eptz);
  tree->SetBranchAddress("true_Eptr", &true_Eptr);
  tree->SetBranchAddress("true_PDGID", &true_PDGID);

  tree->SetBranchAddress("trueTauDaughterVtxx", &trueTauDaughterVtxx);
  tree->SetBranchAddress("trueTauDaughterVtxy", &trueTauDaughterVtxy);
  tree->SetBranchAddress("trueTauDaughterVtxz", &trueTauDaughterVtxz);
  tree->SetBranchAddress("trueTauDaughterVtxr", &trueTauDaughterVtxr);
  tree->SetBranchAddress("trueTauDaughterMotherPDGID", &trueTauDaughterMotherPDGID);


  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);
    for(unsigned int i =0;i<true_Eptx->size();i++){
      //B hadrons
      if( (abs((*true_PDGID)[i])>500 && abs((*true_PDGID)[i])<550) || (abs((*true_PDGID)[i])>5100 && abs((*true_PDGID)[i])<5550)){
	h_hist_vec[0]->Fill((*true_Eptr)[i]);
	h_hist_vec[1]->Fill((*true_Eptz)[i]);
	//B mesons
	if(abs((*true_PDGID)[i])>500 && abs((*true_PDGID)[i])<550){
	  h_hist_vec[2]->Fill((*true_Eptr)[i]);
	  h_hist_vec[3]->Fill((*true_Eptz)[i]);
	}
	//B baryons
	if(abs((*true_PDGID)[i])>5100 && abs((*true_PDGID)[i])<5550){
	  h_hist_vec[4]->Fill((*true_Eptr)[i]);
	  h_hist_vec[5]->Fill((*true_Eptz)[i]);
	}
      }
      //C hadrons
      if( (abs((*true_PDGID)[i])>400 && abs((*true_PDGID)[i])<450) || (abs((*true_PDGID)[i])>4100 && abs((*true_PDGID)[i])<4550)){
	h_hist_vec[6]->Fill((*true_Eptr)[i]);
	h_hist_vec[7]->Fill((*true_Eptz)[i]);
	//C mesons
	if(abs((*true_PDGID)[i])>400 && abs((*true_PDGID)[i])<450){
	  h_hist_vec[8]->Fill((*true_Eptr)[i]);
	  h_hist_vec[9]->Fill((*true_Eptz)[i]);
	}
	//C baryons
	if(abs((*true_PDGID)[i])>4100 && abs((*true_PDGID)[i])<4550){
	  h_hist_vec[10]->Fill((*true_Eptr)[i]);
	  h_hist_vec[11]->Fill((*true_Eptz)[i]);
	}
      }
      for(unsigned int i =0;i<trueTauDaughterVtxx->size();i++){
      //B hadrons
      if( (abs((*trueTauDaughterMotherPDGID)[i])>500 && abs((*trueTauDaughterMotherPDGID)[i])<550) || (abs((*trueTauDaughterMotherPDGID)[i])>5100 && abs((*trueTauDaughterMotherPDGID)[i])<5550)){
	h_hist_vec[12]->Fill((*trueTauDaughterVtxr)[i]);
	h_hist_vec[13]->Fill((*trueTauDaughterVtxz)[i]);
	//B mesons
	if(abs((*trueTauDaughterMotherPDGID)[i])>500 && abs((*trueTauDaughterMotherPDGID)[i])<550){
	  h_hist_vec[14]->Fill((*trueTauDaughterVtxr)[i]);
	  h_hist_vec[15]->Fill((*trueTauDaughterVtxz)[i]);
	}
	//B baryons
	if(abs((*trueTauDaughterMotherPDGID)[i])>5100 && abs((*trueTauDaughterMotherPDGID)[i])<5550){
	  h_hist_vec[16]->Fill((*trueTauDaughterVtxr)[i]);
	  h_hist_vec[17]->Fill((*trueTauDaughterVtxz)[i]);
	}
      }
      //C hadrons
      if( (abs((*trueTauDaughterMotherPDGID)[i])>400 && abs((*trueTauDaughterMotherPDGID)[i])<450) || (abs((*trueTauDaughterMotherPDGID)[i])>4100 && abs((*trueTauDaughterMotherPDGID)[i])<4550)){
	h_hist_vec[18]->Fill((*trueTauDaughterVtxr)[i]);
	h_hist_vec[19]->Fill((*trueTauDaughterVtxz)[i]);
	//C mesons
	if(abs((*trueTauDaughterMotherPDGID)[i])>400 && abs((*trueTauDaughterMotherPDGID)[i])<450){
	  h_hist_vec[20]->Fill((*trueTauDaughterVtxr)[i]);
	  h_hist_vec[21]->Fill((*trueTauDaughterVtxz)[i]);
	}
	//C baryons
	if(abs((*trueTauDaughterMotherPDGID)[i])>4100 && abs((*trueTauDaughterMotherPDGID)[i])<4550){
	  h_hist_vec[22]->Fill((*trueTauDaughterVtxr)[i]);
	  h_hist_vec[23]->Fill((*trueTauDaughterVtxz)[i]);
	}
      }
      }
    }
  }
  delete tree;
}




void JERPlotsFull(){

  CLICdpStyle();
  gROOT->ProcessLine("#include <vector>");

  const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles181123/Zuds_CT_181011_histoFitRange_0_00_to_2_00_RedXRange_centralQuark_cosTheta_0_70_emptyLastCosThetaBin.root";


  string label_legend= "#gamma, E_{true}=1 GeV";


  unsigned int n_bins100=100;
  float lim_E_rel_low=0.0;
  float lim_E_rel_high=2.00;

  const unsigned int nRegionBins(14);
  float pRegionBinEdges[nRegionBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 0.985,1.025};

    


  float resolutionTH1(0.f), resolutionErrorTH1(0.f),meanTH1(0.f),meanErrorTH1(0.f);
  resolutionTH1 = 0.0;
  resolutionErrorTH1 = 0.0;

  meanTH1 = 0.0;
  meanErrorTH1 = 0.0;
  
  /*
    TFile* file_CLIC_100_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/180724_gcc62/photonStudy_Zuds100_10812_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_200_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/180724_gcc62/photonStudy_Zuds200_10816_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_380_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/180724_gcc62/photonStudy_Zuds380_10820_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_500_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles//180724_gcc62/photonStudy_Zuds500_10824_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_1500_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/180724_gcc62/photonStudy_Zuds1500_10828_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_3000_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/180724_gcc62/photonStudy_Zuds3000_10832_no3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  */

  TFile* file_CLIC_100_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds100_11621_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_200_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds200_11629_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_380_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds380_11637_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_500_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds500_11483_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_1500_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds1500_11495_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_3000_CT=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11645_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");

  TFile* file_CLIC_100_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds100_11503_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_200_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds200_11507_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_380_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds380_11511_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_500_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds500_11515_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_1500_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds1500_11527_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_CLIC_3000_CT_noSC=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11251_noOverlay_noSC_CLIC_o3_v14_CT_PandoraPFOs.root");


  TFile* file_noOverlay      =TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_z_n1n1_3000_11614_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_OverlaySelected=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_z_n1n1_3000_11615_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  TFile* file_OverlayTight=   TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_z_n1n1_3000_11615_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  TFile* file_OverlayAll=     TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_z_n1n1_3000_11615_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_OverlayLoose=   TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_z_n1n1_3000_11615_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");

  TFile* file_noZuds3000_Overlay      =TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11645_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_Zuds3000_OverlaySelected=TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11646_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  TFile* file_Zuds3000_OverlayTight=   TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11646_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  TFile* file_Zuds3000_OverlayAll=     TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11646_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  TFile* file_Zuds3000_OverlayLoose=   TFile::Open("/eos/user/w/weberma2/data/photonMacroFiles/181011_gcc62/photonStudy_Zuds3000_11646_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");

  TFile* file_bb_GENInfo     = TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/180518_gcc62_CT/photonStudy_bb3000_10225_DR200_noOverlay_CLIC_o3_v14_SWC_CLIC_GEN_PandoraPFO_FullMC_GENEvents.root");
  TFile* file_tt_GENInfo     = TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/180518_gcc62_CT/photonStudy_tt3000_10229_DR200_noOverlay_CLIC_o3_v14_SWC_CLIC_GEN_PandoraPFO_FullMC_GENEventsForBAndCHadrons.root");

  TFile* file_histogram=new TFile(final_histo_name,"recreate");


  int n_bins1000 = 100;
  float lim_Vtx_low = -500;
  float lim_Vtx_high = 500;
  /*
  TH1F* h_bbar3000_GEN_b_Hadron_EptR =  new TH1F("h_bbar3000_GEN_b_Hadron_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Hadron_EptZ =  new TH1F("h_bbar3000_GEN_b_Hadron_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Meson_EptR =  new TH1F("h_bbar3000_GEN_b_Meson_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Meson_EptZ =  new TH1F("h_bbar3000_GEN_b_Meson_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Baryon_EptR =  new TH1F("h_bbar3000_GEN_b_Baryon_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Baryon_EptZ =  new TH1F("h_bbar3000_GEN_b_Baryon_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_bbar3000_GEN_b_Hadron_DaughtVtxR =  new TH1F("h_bbar3000_GEN_b_Hadron_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Hadron_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_b_Hadron_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Meson_DaughtVtxR =  new TH1F("h_bbar3000_GEN_b_Meson_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Meson_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_b_Meson_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Baryon_DaughtVtxR =  new TH1F("h_bbar3000_GEN_b_Baryon_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_b_Baryon_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_b_Baryon_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_bbar3000_GEN_c_Hadron_EptR =  new TH1F("h_bbar3000_GEN_c_Hadron_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Hadron_EptZ =  new TH1F("h_bbar3000_GEN_c_Hadron_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Meson_EptR =  new TH1F("h_bbar3000_GEN_c_Meson_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Meson_EptZ =  new TH1F("h_bbar3000_GEN_c_Meson_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Baryon_EptR =  new TH1F("h_bbar3000_GEN_c_Baryon_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Baryon_EptZ =  new TH1F("h_bbar3000_GEN_c_Baryon_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_bbar3000_GEN_c_Hadron_DaughtVtxR =  new TH1F("h_bbar3000_GEN_c_Hadron_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Hadron_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_c_Hadron_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Meson_DaughtVtxR =  new TH1F("h_bbar3000_GEN_c_Meson_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Meson_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_c_Meson_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Baryon_DaughtVtxR =  new TH1F("h_bbar3000_GEN_c_Baryon_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_bbar3000_GEN_c_Baryon_DaughtVtxZ =  new TH1F("h_bbar3000_GEN_c_Baryon_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  
  std::vector<TH1F*> hist_vec_3TeV_bb;
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Hadron_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Hadron_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Meson_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Meson_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Baryon_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Baryon_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Hadron_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Hadron_DaughtVtxZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Meson_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Meson_DaughtVtxZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Baryon_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_b_Baryon_DaughtVtxZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Hadron_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Hadron_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Meson_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Meson_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Baryon_EptR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Baryon_EptZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Hadron_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Hadron_DaughtVtxZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Meson_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Meson_DaughtVtxZ);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Baryon_DaughtVtxR);
  hist_vec_3TeV_bb.push_back(h_bbar3000_GEN_c_Baryon_DaughtVtxZ);

  for(unsigned int i=0;i<hist_vec_3TeV_bb.size();i++){
    hist_vec_3TeV_bb[i]->Sumw2();
    hist_vec_3TeV_bb[i]->SetLineWidth(2);
    hist_vec_3TeV_bb[i]->SetLineColor(kBlack);
  }

  //fill_B_and_D_histograms(file_bb_GENInfo,hist_vec_3TeV_bb);
  
  TH1F* h_ttar3000_GEN_b_Hadron_EptR =  new TH1F("h_ttar3000_GEN_b_Hadron_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Hadron_EptZ =  new TH1F("h_ttar3000_GEN_b_Hadron_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Meson_EptR =  new TH1F("h_ttar3000_GEN_b_Meson_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Meson_EptZ =  new TH1F("h_ttar3000_GEN_b_Meson_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Baryon_EptR =  new TH1F("h_ttar3000_GEN_b_Baryon_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Baryon_EptZ =  new TH1F("h_ttar3000_GEN_b_Baryon_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_ttar3000_GEN_b_Hadron_DaughtVtxR =  new TH1F("h_ttar3000_GEN_b_Hadron_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Hadron_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_b_Hadron_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Meson_DaughtVtxR =  new TH1F("h_ttar3000_GEN_b_Meson_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Meson_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_b_Meson_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Baryon_DaughtVtxR =  new TH1F("h_ttar3000_GEN_b_Baryon_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_b_Baryon_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_b_Baryon_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_ttar3000_GEN_c_Hadron_EptR =  new TH1F("h_ttar3000_GEN_c_Hadron_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Hadron_EptZ =  new TH1F("h_ttar3000_GEN_c_Hadron_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Meson_EptR =  new TH1F("h_ttar3000_GEN_c_Meson_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Meson_EptZ =  new TH1F("h_ttar3000_GEN_c_Meson_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Baryon_EptR =  new TH1F("h_ttar3000_GEN_c_Baryon_EptR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Baryon_EptZ =  new TH1F("h_ttar3000_GEN_c_Baryon_EptZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);

  TH1F* h_ttar3000_GEN_c_Hadron_DaughtVtxR =  new TH1F("h_ttar3000_GEN_c_Hadron_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Hadron_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_c_Hadron_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Meson_DaughtVtxR =  new TH1F("h_ttar3000_GEN_c_Meson_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Meson_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_c_Meson_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Baryon_DaughtVtxR =  new TH1F("h_ttar3000_GEN_c_Baryon_DaughtVtxR","", n_bins1000, 0,lim_Vtx_high);
  TH1F* h_ttar3000_GEN_c_Baryon_DaughtVtxZ =  new TH1F("h_ttar3000_GEN_c_Baryon_DaughtVtxZ","", n_bins1000, lim_Vtx_low,lim_Vtx_high);
  
  std::vector<TH1F*> hist_vec_3TeV_tt;
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Hadron_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Hadron_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Meson_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Meson_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Baryon_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Baryon_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Hadron_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Hadron_DaughtVtxZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Meson_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Meson_DaughtVtxZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Baryon_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_b_Baryon_DaughtVtxZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Hadron_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Hadron_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Meson_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Meson_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Baryon_EptR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Baryon_EptZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Hadron_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Hadron_DaughtVtxZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Meson_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Meson_DaughtVtxZ);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Baryon_DaughtVtxR);
  hist_vec_3TeV_tt.push_back(h_ttar3000_GEN_c_Baryon_DaughtVtxZ);

  for(unsigned int i=0;i<hist_vec_3TeV_tt.size();i++){
    hist_vec_3TeV_tt[i]->Sumw2();
    hist_vec_3TeV_tt[i]->SetLineWidth(2);
    hist_vec_3TeV_tt[i]->SetLineColor(kBlack);
  }

  fill_B_and_D_histograms(file_tt_GENInfo,hist_vec_3TeV_tt);
  */

  TGraphErrors *gre_DR07_RMS90_JER_0_70_wSC = new TGraphErrors(6);
  gre_DR07_RMS90_JER_0_70_wSC ->SetName("gre_DR07_RMS90_JER_0_70_wSC");
  gre_DR07_RMS90_JER_0_70_wSC ->SetTitle("");
  gre_DR07_RMS90_JER_0_70_wSC ->SetFillColor(1);
  gre_DR07_RMS90_JER_0_70_wSC ->SetLineColor(kBlack);
  gre_DR07_RMS90_JER_0_70_wSC ->SetMarkerColor(kBlack);
  gre_DR07_RMS90_JER_0_70_wSC ->SetMarkerStyle(kFullCircle);

  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_100_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_100_CT;
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_100_CT.push_back( h_100_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_100_CT.size();i++){
    hist_vec_100_CT[i]->Sumw2();
    hist_vec_100_CT[i]->SetLineWidth(2);
    hist_vec_100_CT[i]->SetLineColor(kCyan+1);
    hist_vec_100_CT[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_100_CT,hist_vec_100_CT);

  TH1F* h_100_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_100", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_100_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_100_CT_E_SigmaCB_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_100_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_100", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_100_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_100_CT_E_RMS90_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_100_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_100", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_100_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_100_CT_E_RMS_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_100_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_100", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_100_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_100_CT_E_JER_Mean_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_100_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_100_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_100_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_100_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_100_CT[i]->GetRMS()*100. );
    h_100_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_100_CT[i]->GetRMSError()*100.);
    h_100_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_100_CT[i]->GetMean());
    h_100_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_100_CT[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_100_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_100_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_100_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_100_CT[i]->GetFunction("fdscb_tmp");
    h_100_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_100_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  

  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_100_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(0,50,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(0,0,resolutionErrorTH1);

  std::cout<<"100 done"<<std::endl;
 

  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_200_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_200_CT;
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_200_CT.push_back( h_200_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_200_CT.size();i++){
    hist_vec_200_CT[i]->Sumw2();
    hist_vec_200_CT[i]->SetLineWidth(2);
    hist_vec_200_CT[i]->SetLineColor(kYellow+1);
    hist_vec_200_CT[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_200_CT,hist_vec_200_CT);
  TH1F* h_200_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_200", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_200_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_200_CT_E_SigmaCB_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_200_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_200", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_200_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_200_CT_E_RMS90_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_200_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_200", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_200_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_200_CT_E_RMS_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_200_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_200", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_200_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_200_CT_E_JER_Mean_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_200_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_200_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_200_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_200_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_200_CT[i]->GetRMS()*100. );
    h_200_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_200_CT[i]->GetRMSError()*100.);
    h_200_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_200_CT[i]->GetMean());
    h_200_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_200_CT[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_200_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_200_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_200_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_200_CT[i]->GetFunction("fdscb_tmp");
    h_200_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_200_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  

  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_200_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(1,100,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(1,0,resolutionErrorTH1);

  std::cout<<"200 done"<<std::endl;
  //200 done here, 380 starts now

  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_380_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_380_CT;
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_380_CT.push_back( h_380_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_380_CT.size();i++){
    hist_vec_380_CT[i]->Sumw2();
    hist_vec_380_CT[i]->SetLineWidth(2);
    hist_vec_380_CT[i]->SetLineColor(kOrange);
    hist_vec_380_CT[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_380_CT,hist_vec_380_CT);

  TH1F* h_380_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_380", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_380_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_380_CT_E_SigmaCB_summary_histogram->SetLineColor(kOrange);
  h_380_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_380_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_380", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_380_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_380_CT_E_RMS90_summary_histogram->SetLineColor(kOrange);
  h_380_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_380_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_380", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_380_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_380_CT_E_RMS_summary_histogram->SetLineColor(kOrange);
  h_380_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_380_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_380", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_380_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_380_CT_E_JER_Mean_summary_histogram->SetLineColor(kOrange);
  h_380_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_380_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_380_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_380_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_380_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_380_CT[i]->GetRMS()*100. );
    h_380_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_380_CT[i]->GetRMSError()*100.);
    h_380_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_380_CT[i]->GetMean());
    h_380_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_380_CT[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_380_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_380_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_380_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_380_CT[i]->GetFunction("fdscb_tmp");
    h_380_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_380_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  

  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_380_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(2,190,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(2,0,resolutionErrorTH1);
  std::cout<<"380 done"<<std::endl;

  //380 done here, 380 starts now

  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_500_CT;
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_500_CT.push_back( h_500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_500_CT.size();i++){
    hist_vec_500_CT[i]->Sumw2();
    hist_vec_500_CT[i]->SetLineWidth(2);
    hist_vec_500_CT[i]->SetLineColor(kGreen-2);
    hist_vec_500_CT[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_500_CT,hist_vec_500_CT);

  TH1F* h_500_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_500", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_500_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_500_CT_E_SigmaCB_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_500_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_500", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_500_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_500_CT_E_RMS90_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_500_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_500", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_500_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_500_CT_E_RMS_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_500_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_500", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_500_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_500_CT_E_JER_Mean_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_500_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_500_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_500_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_500_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_500_CT[i]->GetRMS()*100. );
    h_500_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_500_CT[i]->GetRMSError()*100.);
    h_500_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_500_CT[i]->GetMean());
    h_500_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_500_CT[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_500_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_500_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_500_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_500_CT[i]->GetFunction("fdscb_tmp");
    h_500_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_500_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_500_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(3,250,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(3,0,resolutionErrorTH1);

  std::cout<<"500 done"<<std::endl;
  //500 done, go now to 1500

  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_1500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_1500_CT;
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_1500_CT.push_back( h_1500_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_1500_CT.size();i++){
    hist_vec_1500_CT[i]->Sumw2();
    hist_vec_1500_CT[i]->SetLineWidth(2);
    hist_vec_1500_CT[i]->SetLineColor(kBlue);
    hist_vec_1500_CT[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_1500_CT,hist_vec_1500_CT);

  TH1F* h_1500_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_1500", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_1500_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_1500_CT_E_SigmaCB_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_1500_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_1500", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_1500_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_1500_CT_E_RMS90_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_1500_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_1500", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_1500_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_1500_CT_E_RMS_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_1500_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_1500", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_1500_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_1500_CT_E_JER_Mean_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_1500_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_1500_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_1500_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_1500_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_1500_CT[i]->GetRMS()*100. );
    h_1500_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_1500_CT[i]->GetRMSError()*100.);
    h_1500_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_1500_CT[i]->GetMean());
    h_1500_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_1500_CT[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_1500_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_1500_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_1500_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_1500_CT[i]->GetFunction("fdscb_tmp");
    h_1500_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_1500_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_1500_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(4,750,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(4,0,resolutionErrorTH1);

  std::cout<<"1500 done"<<std::endl;
  //1500 done, go now to 3000

  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_3000_CT_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_3000_CT;
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_3000_CT.push_back( h_3000_CT_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_3000_CT.size();i++){
    hist_vec_3000_CT[i]->Sumw2();
    hist_vec_3000_CT[i]->SetLineWidth(2);
    hist_vec_3000_CT[i]->SetLineColor(kRed);
    hist_vec_3000_CT[i]->SetLineStyle(1);
  }
  fill_resolution_histograms(file_CLIC_3000_CT,hist_vec_3000_CT);

  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->SetLineColor(kBlack);
  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->SetLineColor(kBlue);
  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->Scale(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->Integral());
  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->SetLineColor(kMagenta);
  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->Scale(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->Integral());
  h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->Scale(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->Integral());

  TCanvas *resolutionGraphCanvas_JER_totE_3000_RMS90_CT_fancy = setUpperCanvas("resolutionGraphCanvas_JER_totE_3000_RMS90_CT_fancy");
  //resolutionGraphCanvas_JER_totE_3000_RMS90_CT_fancy->cd();
  //TLegend *leg_JER_totE_3000_RMS90_CT_FullSummary = resolutionGraphCanvas_JER_totE_3000_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_totE_3000_RMS90_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetBorderSize(0);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetTextAlign(12);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetTextSize(0.050);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetTextFont(42);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetMargin(0.15);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetLineColor(1);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetLineStyle(1);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetLineWidth(1);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetFillColor(0);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetFillStyle(1001);
  leg_JER_totE_3000_RMS90_CT_FullSummary->SetHeader("1500 GeV Jets");
  leg_JER_totE_3000_RMS90_CT_FullSummary->AddEntry(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->DrawCopy("h,e"),"0.80<|cos#theta|<0.90");
  leg_JER_totE_3000_RMS90_CT_FullSummary->AddEntry(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->DrawCopy("h,e,same"),"0.900<|cos#theta|<0.925");
  leg_JER_totE_3000_RMS90_CT_FullSummary->AddEntry(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->DrawCopy("h,e,same"),"0.925<|cos#theta|<0.95");
  leg_JER_totE_3000_RMS90_CT_FullSummary->AddEntry(h_3000_CT_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->DrawCopy("h,e,same"),"0.975<|cos#theta|<1.00");
  leg_JER_totE_3000_RMS90_CT_FullSummary->Draw();

  //l->DrawLatex(x,y,label.c_str());

  TH1F* h_3000_CT_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_3000", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_3000_CT_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_3000_CT_E_SigmaCB_summary_histogram->SetLineColor(kRed);
  h_3000_CT_E_SigmaCB_summary_histogram->SetLineStyle(2);

  TH1F* h_3000_CT_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_3000", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_3000_CT_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_E_RMS90_summary_histogram->SetLineWidth(2);
  h_3000_CT_E_RMS90_summary_histogram->SetLineColor(kRed);
  h_3000_CT_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_3000_CT_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_3000", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_3000_CT_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_E_RMS_summary_histogram->SetLineWidth(2);
  h_3000_CT_E_RMS_summary_histogram->SetLineColor(kRed);
  h_3000_CT_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_3000_CT_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_3000", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_3000_CT_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_3000_CT_E_JER_Mean_summary_histogram->SetLineColor(kRed);
  h_3000_CT_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_3000_CT[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_3000_CT_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_3000_CT_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_3000_CT_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_3000_CT[i]->GetRMS()*100. );
    h_3000_CT_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_3000_CT[i]->GetRMSError()*100.);
    h_3000_CT_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_3000_CT[i]->GetMean());
    h_3000_CT_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_3000_CT[i]->GetMeanError());
    //double sided crystal ball
    /*TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_3000_CT[i]->Fit("fit_gaus_tmp","R");
    hist_vec_3000_CT[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_3000_CT[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_3000_CT[i]->GetFunction("fdscb_tmp");
    h_3000_CT_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_3000_CT_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
    */
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_3000_CT[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPoint(5,1500,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_wSC->SetPointError(5,0,resolutionErrorTH1);
  //try the noSC histograms now

  TGraphErrors *gre_DR07_RMS90_JER_0_70_noSC = new TGraphErrors(6);
  gre_DR07_RMS90_JER_0_70_noSC ->SetName("gre_DR07_RMS90_JER_0_70_noSC");
  gre_DR07_RMS90_JER_0_70_noSC ->SetTitle("");
  gre_DR07_RMS90_JER_0_70_noSC ->SetFillColor(1);
  gre_DR07_RMS90_JER_0_70_noSC ->SetLineColor(kRed);
  gre_DR07_RMS90_JER_0_70_noSC ->SetMarkerColor(kRed);
  gre_DR07_RMS90_JER_0_70_noSC ->SetMarkerStyle(kFullSquare);

  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_100_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_100_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_100_CT_noSC;
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_100_CT_noSC.push_back( h_100_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_100_CT_noSC.size();i++){
    hist_vec_100_CT_noSC[i]->Sumw2();
    hist_vec_100_CT_noSC[i]->SetLineWidth(2);
    hist_vec_100_CT_noSC[i]->SetLineColor(kCyan+1);
    hist_vec_100_CT_noSC[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_100_CT_noSC,hist_vec_100_CT_noSC);

  TH1F* h_100_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_100", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_100_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_100_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_100_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_100", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_100_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_100_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_100_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_100", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_100_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_100_CT_noSC_E_RMS_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_100_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_100", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_100_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_100_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_100_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_100_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kCyan+1);
  h_100_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_100_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_100_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_100_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_100_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_100_CT_noSC[i]->GetRMS()*100. );
    h_100_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_100_CT_noSC[i]->GetRMSError()*100.);
    h_100_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_100_CT_noSC[i]->GetMean());
    h_100_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_100_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_100_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_100_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_100_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_100_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_100_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_100_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_100_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(0,50,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(0,0,resolutionErrorTH1);

  std::cout<<"100 done"<<std::endl;


  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_200_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_200_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_200_CT_noSC;
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_200_CT_noSC.push_back( h_200_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_200_CT_noSC.size();i++){
    hist_vec_200_CT_noSC[i]->Sumw2();
    hist_vec_200_CT_noSC[i]->SetLineWidth(2);
    hist_vec_200_CT_noSC[i]->SetLineColor(kYellow+1);
    hist_vec_200_CT_noSC[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_200_CT_noSC,hist_vec_200_CT_noSC);
  TH1F* h_200_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_200", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_200_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_200_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_200_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_200", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_200_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_200_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_200_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_200", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_200_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_200_CT_noSC_E_RMS_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_200_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_200", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_200_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_200_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_200_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_200_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kYellow+1);
  h_200_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_200_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_200_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_200_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_200_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_200_CT_noSC[i]->GetRMS()*100. );
    h_200_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_200_CT_noSC[i]->GetRMSError()*100.);
    h_200_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_200_CT_noSC[i]->GetMean());
    h_200_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_200_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_200_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_200_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_200_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_200_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_200_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_200_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_200_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(1,100,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(1,0,resolutionErrorTH1);

  std::cout<<"200 done"<<std::endl;
  //200 done here, 380 starts now

  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_380_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_380_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_380_CT_noSC;
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_380_CT_noSC.push_back( h_380_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_380_CT_noSC.size();i++){
    hist_vec_380_CT_noSC[i]->Sumw2();
    hist_vec_380_CT_noSC[i]->SetLineWidth(2);
    hist_vec_380_CT_noSC[i]->SetLineColor(kOrange);
    hist_vec_380_CT_noSC[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_380_CT_noSC,hist_vec_380_CT_noSC);

  TH1F* h_380_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_380", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_380_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_380_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kOrange);
  h_380_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_380_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_380", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_380_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_380_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kOrange);
  h_380_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_380_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_380", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_380_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_380_CT_noSC_E_RMS_summary_histogram->SetLineColor(kOrange);
  h_380_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_380_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_380", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_380_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_380_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_380_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_380_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kOrange);
  h_380_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_380_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_380_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_380_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_380_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_380_CT_noSC[i]->GetRMS()*100. );
    h_380_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_380_CT_noSC[i]->GetRMSError()*100.);
    h_380_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_380_CT_noSC[i]->GetMean());
    h_380_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_380_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_380_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_380_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_380_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_380_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_380_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_380_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_380_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(2,190,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(2,0,resolutionErrorTH1);
  std::cout<<"380 done"<<std::endl;
  //380 done here, 380 starts now

  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_500_CT_noSC;
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_500_CT_noSC.push_back( h_500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_500_CT_noSC.size();i++){
    hist_vec_500_CT_noSC[i]->Sumw2();
    hist_vec_500_CT_noSC[i]->SetLineWidth(2);
    hist_vec_500_CT_noSC[i]->SetLineColor(kGreen-2);
    hist_vec_500_CT_noSC[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_500_CT_noSC,hist_vec_500_CT_noSC);

  TH1F* h_500_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_500", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_500_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_500_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_500_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_500", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_500_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_500_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_500_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_500", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_500_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_500_CT_noSC_E_RMS_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_500_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_500", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_500_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_500_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_500_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_500_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kGreen-2);
  h_500_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_500_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_500_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_500_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_500_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_500_CT_noSC[i]->GetRMS()*100. );
    h_500_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_500_CT_noSC[i]->GetRMSError()*100.);
    h_500_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_500_CT_noSC[i]->GetMean());
    h_500_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_500_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_500_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_500_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_500_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_500_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_500_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_500_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_500_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(3,250,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(3,0,resolutionErrorTH1);
  std::cout<<"500 done"<<std::endl;
  //500 done, go now to 1500

  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_1500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_1500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_1500_CT_noSC;
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_1500_CT_noSC.push_back( h_1500_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_1500_CT_noSC.size();i++){
    hist_vec_1500_CT_noSC[i]->Sumw2();
    hist_vec_1500_CT_noSC[i]->SetLineWidth(2);
    hist_vec_1500_CT_noSC[i]->SetLineColor(kBlue);
    hist_vec_1500_CT_noSC[i]->SetLineStyle(3);
  }
  fill_resolution_histograms(file_CLIC_1500_CT_noSC,hist_vec_1500_CT_noSC);

  TH1F* h_1500_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_1500", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);
  TH1F* h_1500_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_1500", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_1500_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_1500_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_1500_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_1500", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_1500_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_1500_CT_noSC_E_RMS_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_1500_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_1500", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kBlue);
  h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_1500_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_1500_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_1500_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_1500_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_1500_CT_noSC[i]->GetRMS()*100. );
    h_1500_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_1500_CT_noSC[i]->GetRMSError()*100.);
    h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_1500_CT_noSC[i]->GetMean());
    h_1500_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_1500_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_1500_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_1500_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_1500_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_1500_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_1500_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_1500_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(4,750,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(4,0,resolutionErrorTH1);
  std::cout<<"1500 done"<<std::endl;
  //1500 done, go now to 3000

  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);
  TH1F* h_3000_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70 =  new TH1F("h_3000_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70","", n_bins100, lim_E_rel_low,lim_E_rel_high);

  std::vector<TH1F*> hist_vec_3000_CT_noSC;
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_to_0_10);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_10_to_0_20);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_20_to_0_30);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_30_to_0_40);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_40_to_0_50);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_50_to_0_60);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_60_to_0_70);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70_to_0_80);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_950_to_0_975);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_70);
  hist_vec_3000_CT_noSC.push_back( h_3000_CT_noSC_E_rel_totVis_RMS_j1j2_cosT_0_70);
  for(unsigned int i=0;i<hist_vec_3000_CT_noSC.size();i++){
    hist_vec_3000_CT_noSC[i]->Sumw2();
    hist_vec_3000_CT_noSC[i]->SetLineWidth(2);
    hist_vec_3000_CT_noSC[i]->SetLineColor(kRed);
    hist_vec_3000_CT_noSC[i]->SetLineStyle(1);
  }
  fill_resolution_histograms(file_CLIC_3000_CT_noSC,hist_vec_3000_CT_noSC);

  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->SetLineColor(kBlack);
  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->SetLineColor(kBlue);
  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->Scale(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->Integral());
  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->SetLineColor(kMagenta);
  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->Scale(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->Integral());
  h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->Scale(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->Integral()/h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->Integral());

  TCanvas *resolutionGraphCanvas_JER_totE_3000_RMS90_CT_noSC_fancy = setUpperCanvas("resolutionGraphCanvas_JER_totE_3000_RMS90_CT_noSC_fancy");
  //resolutionGraphCanvas_JER_totE_3000_RMS90_CT_noSC_fancy->cd();
  //TLegend *leg_JER_totE_3000_RMS90_CT_noSC_FullSummary = resolutionGraphCanvas_JER_totE_3000_RMS90_CT_noSC_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_totE_3000_RMS90_CT_noSC_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetBorderSize(0);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetTextAlign(12);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetTextSize(0.050);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetTextFont(42);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetMargin(0.15);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetLineColor(1);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetLineStyle(1);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetLineWidth(1);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetFillColor(0);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetFillStyle(1001);
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->SetHeader("1500 GeV Jets");
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->AddEntry(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_80_to_0_90->DrawCopy("h,e"),"0.80<|cos#theta|<0.90");
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->AddEntry(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_900_to_0_925->DrawCopy("h,e,same"),"0.900<|cos#theta|<0.925");
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->AddEntry(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950->DrawCopy("h,e,same"),"0.925<|cos#theta|<0.95");
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->AddEntry(h_3000_CT_noSC_E_rel_totVis_RMS_d1_cosT_0_975_to_0_985->DrawCopy("h,e,same"),"0.975<|cos#theta|<1.00");
  leg_JER_totE_3000_RMS90_CT_noSC_FullSummary->Draw();

  //l->DrawLatex(x,y,label.c_str());

  TH1F* h_3000_CT_noSC_E_SigmaCB_summary_histogram=new TH1F("JER_SigmaCBVsCosTheta_CT_noSC_3000", "#sigma (E_{j}^{reco} / E_{j}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_noSC_E_SigmaCB_summary_histogram ->SetYTitle("#sigma(E_{j}) / E_{j} [%]");
  h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetLineWidth(2);
  h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetLineColor(kRed);
  h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetLineStyle(2);

  TH1F* h_3000_CT_noSC_E_RMS90_summary_histogram=new TH1F("JER_RMS90VsCosTheta_CT_noSC_3000", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_noSC_E_RMS90_summary_histogram ->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
  h_3000_CT_noSC_E_RMS90_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_noSC_E_RMS90_summary_histogram->SetLineWidth(2);
  h_3000_CT_noSC_E_RMS90_summary_histogram->SetLineColor(kRed);
  h_3000_CT_noSC_E_RMS90_summary_histogram->SetLineStyle(3);

  TH1F* h_3000_CT_noSC_E_RMS_summary_histogram=new TH1F("JER_RMSVsCosTheta_CT_noSC_3000", "RMS(E_{j}) / Mean(E_{j}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_noSC_E_RMS_summary_histogram ->SetYTitle("RMS (E_{j}) / Mean(E_{j}) [%]");
  h_3000_CT_noSC_E_RMS_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_noSC_E_RMS_summary_histogram->SetLineWidth(2);
  h_3000_CT_noSC_E_RMS_summary_histogram->SetLineColor(kRed);
  h_3000_CT_noSC_E_RMS_summary_histogram->SetLineStyle(3);

  TH1F* h_3000_CT_noSC_E_JER_Mean_summary_histogram=new TH1F("JER_MeanVsCosTheta_CT_noSC_3000", "Mean(E_{tot}^{RECO}/E_{tot}^{true}) vs |cos#theta|", nRegionBins, pRegionBinEdges);
  h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetYTitle("Mean(E_{tot}^{RECO}/E_{tot}^{true})");
  h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetXTitle("|cos#theta|");
  h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetLineWidth(2);
  h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetLineColor(kRed);
  h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetLineStyle(3);

  for(unsigned int i=0;i<13;i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformance( hist_vec_3000_CT_noSC[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    h_3000_CT_noSC_E_RMS90_summary_histogram->SetBinContent(i+1,resolutionTH1 );
    h_3000_CT_noSC_E_RMS90_summary_histogram ->SetBinError(i+1, resolutionErrorTH1);
    h_3000_CT_noSC_E_RMS_summary_histogram->SetBinContent(i+1,hist_vec_3000_CT_noSC[i]->GetRMS()*100. );
    h_3000_CT_noSC_E_RMS_summary_histogram ->SetBinError(i+1, hist_vec_3000_CT_noSC[i]->GetRMSError()*100.);
    h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetBinContent(i+1,hist_vec_3000_CT_noSC[i]->GetMean());
    h_3000_CT_noSC_E_JER_Mean_summary_histogram->SetBinError(i+1,hist_vec_3000_CT_noSC[i]->GetMeanError());
    //double sided crystal ball
    /*TF1* fit_gaus_tmp = new TF1("fit_gaus_tmp", "gaus",0.5,1.5);
    hist_vec_3000_CT_noSC[i]->Fit("fit_gaus_tmp","R");
    hist_vec_3000_CT_noSC[i]->GetFunction("fit_gaus_tmp");
    TF1* fdscb_tmp = new TF1("fdscb_tmp",fnc_dscb,0.,3.,7);
    fdscb_tmp->SetParameter(0,fit_gaus_tmp->GetParameter(0));
    //norm
    fdscb_tmp->SetParLimits (0,0.5*fit_gaus_tmp->GetParameter(0),2.0*fit_gaus_tmp->GetParameter(0));
    fdscb_tmp->SetParameter(1,fit_gaus_tmp->GetParameter(1));
    //gaussian mean
    fdscb_tmp->SetParLimits (1,0.8,1.2);
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParLimits(2,fit_gaus_tmp->GetParameter(2)*0.25,4.*fit_gaus_tmp->GetParameter(2));
    //-a1: parameter, from when the gaussian core should count
    fdscb_tmp->SetParameter(2,fit_gaus_tmp->GetParameter(2));
    fdscb_tmp->SetParameter (3,2.0);
    fdscb_tmp->SetParLimits (3,0.25,15.0);
    //a2 parameter, up to when the gaussian core should count
    fdscb_tmp->SetParameter (5,2.0);
    fdscb_tmp->SetParLimits (5,0.25,15.0);
    fdscb_tmp->SetParameter(4,5.);
    fdscb_tmp->SetParLimits(4,0.,25.);
    fdscb_tmp->SetParameter(6,5.);
    fdscb_tmp->SetParLimits(6,0.,25.);
    int fitstatus(0); 
    fitstatus = hist_vec_3000_CT_noSC[i]->Fit(fdscb_tmp,"R");
    fdscb_tmp = hist_vec_3000_CT_noSC[i]->GetFunction("fdscb_tmp");
    h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetBinContent(i+1,sqrt(2.)*fdscb_tmp->GetParameter(2)*100. );
    h_3000_CT_noSC_E_SigmaCB_summary_histogram->SetBinError(i+1,sqrt(2.)*fdscb_tmp->GetParError(2)*100. );
    */
  }  
  resolutionTH1=0; 
  resolutionErrorTH1=0;
  meanTH1=0;
  meanErrorTH1=0;
  CalculatePerformance( hist_vec_3000_CT_noSC[13], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPoint(5,1500,resolutionTH1);
  gre_DR07_RMS90_JER_0_70_noSC->SetPointError(5,0,resolutionErrorTH1);




  int n_bins100_MET=1000;

  float lim_E_low=2500;
  float lim_E_low_jets=2000;
  float lim_E_high_jets=4000;
  float lim_E_high=6000;
  float lim_MET_rel_low=-5.0;
  float lim_MET_rel_high=5.00;

  TH1F* h_OverlayLoose_px_miss_totPFO =  new TH1F("h_OverlayLoose_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_OverlayLoose_delta_MET_totPFO_MET_true =  new TH1F("h_OverlayLoose_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_OverlayLoose_MET_totPFO =  new TH1F("h_OverlayLoose_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayLoose_E_totPFO =  new TH1F("h_OverlayLoose_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_OverlayLoose_E_jetSum =  new TH1F("h_OverlayLoose_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_OverlayLoose_METPhi_totPFO =  new TH1F("h_OverlayLoose_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_OverlayLoose_MHT_totPFO =  new TH1F("h_OverlayLoose_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayLoose_dphi_j1j2_totPFO =  new TH1F("h_OverlayLoose_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_OverlayLoose_MET_true =  new TH1F("h_OverlayLoose_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayLoose_METPhi_true =  new TH1F("h_OverlayLoose_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_OverlayLoose_dphi_MET_true_reco =  new TH1F("h_OverlayLoose_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_0_25 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_0_25 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_25_50 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_25_50 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_50_75 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_50_75 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_75_100 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_75_100 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_100_140 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_100_140 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_140_180 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_140_180 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_180_220 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_180_220 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_220_260 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_220_260 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_260_340 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_260_340 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_340_420 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_340_420 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_420_500 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_420_500 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_500_750 =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_500_750 =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayLoose_deltaMET_over_trueMET_750_Inf =  new TH1F("h_OverlayLoose_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayLoose_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_OverlayLoose_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_OverlayLoose;
  h_vec_OverlayLoose.push_back(h_OverlayLoose_px_miss_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_delta_MET_totPFO_MET_true);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_MET_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_E_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_E_jetSum);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_METPhi_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_MHT_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_j1j2_totPFO);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_MET_true);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_METPhi_true);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_MET_true_reco);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_0_25);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_0_25);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_25_50);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_25_50);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_50_75);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_50_75);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_75_100);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_75_100);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_100_140);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_100_140);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_140_180);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_140_180);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_180_220);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_180_220);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_220_260);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_220_260);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_260_340);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_260_340);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_340_420);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_340_420);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_420_500);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_420_500);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_500_750);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_500_750);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_deltaMET_over_trueMET_750_Inf);
  h_vec_OverlayLoose.push_back(h_OverlayLoose_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_OverlayLoose.size();i++){
    h_vec_OverlayLoose[i]->SetLineWidth(3);
    h_vec_OverlayLoose[i]->SetLineColor(kBlue);
    h_vec_OverlayLoose[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_OverlayLoose,h_vec_OverlayLoose);

  const unsigned int nRegionBins_MET(13);
  //float pRegionBinEdges_MET[nRegionBins_MET + 1] = {0, 25., 50., 75., 100., 140., 180., 220., 260., 340., 420., 500., 750., 1500.};
  float pRegionBinEdges_MET[nRegionBins_MET + 1] = {0, 25., 50., 75., 100., 140., 180., 220., 260., 340., 420., 500., 750., 1010.};

  TH1F* h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlue);
  h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlue);
  h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  int ibin=1;
  
  for(unsigned int i=11;i<h_vec_OverlayLoose.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayLoose[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayLoose[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }
  

  TH1F* h_OverlaySelected_px_miss_totPFO =  new TH1F("h_OverlaySelected_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_OverlaySelected_delta_MET_totPFO_MET_true =  new TH1F("h_OverlaySelected_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_OverlaySelected_MET_totPFO =  new TH1F("h_OverlaySelected_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlaySelected_E_totPFO =  new TH1F("h_OverlaySelected_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_OverlaySelected_E_jetSum =  new TH1F("h_OverlaySelected_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_OverlaySelected_METPhi_totPFO =  new TH1F("h_OverlaySelected_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_OverlaySelected_MHT_totPFO =  new TH1F("h_OverlaySelected_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlaySelected_dphi_j1j2_totPFO =  new TH1F("h_OverlaySelected_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_OverlaySelected_MET_true =  new TH1F("h_OverlaySelected_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_OverlaySelected_METPhi_true =  new TH1F("h_OverlaySelected_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_OverlaySelected_dphi_MET_true_reco =  new TH1F("h_OverlaySelected_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_0_25 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_0_25 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_25_50 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_25_50 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_50_75 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_50_75 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_75_100 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_75_100 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_100_140 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_100_140 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_140_180 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_140_180 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_180_220 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_180_220 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_220_260 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_220_260 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_260_340 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_260_340 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_340_420 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_340_420 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_420_500 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_420_500 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_500_750 =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_500_750 =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlaySelected_deltaMET_over_trueMET_750_Inf =  new TH1F("h_OverlaySelected_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlaySelected_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_OverlaySelected_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_OverlaySelected;
  h_vec_OverlaySelected.push_back(h_OverlaySelected_px_miss_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_delta_MET_totPFO_MET_true);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_MET_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_E_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_E_jetSum);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_METPhi_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_MHT_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_j1j2_totPFO);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_MET_true);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_METPhi_true);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_MET_true_reco);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_0_25);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_0_25);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_25_50);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_25_50);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_50_75);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_50_75);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_75_100);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_75_100);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_100_140);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_100_140);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_140_180);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_140_180);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_180_220);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_180_220);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_220_260);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_220_260);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_260_340);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_260_340);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_340_420);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_340_420);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_420_500);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_420_500);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_500_750);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_500_750);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_deltaMET_over_trueMET_750_Inf);
  h_vec_OverlaySelected.push_back(h_OverlaySelected_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_OverlaySelected.size();i++){
    h_vec_OverlaySelected[i]->SetLineWidth(3);
    h_vec_OverlaySelected[i]->SetLineColor(kRed-7);
    h_vec_OverlaySelected[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_OverlaySelected,h_vec_OverlaySelected);

  TH1F* h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kRed-7);
  h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kRed-7);
  h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_OverlaySelected.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlaySelected[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlaySelected[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }

  TH1F* h_OverlayTight_px_miss_totPFO =  new TH1F("h_OverlayTight_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_OverlayTight_delta_MET_totPFO_MET_true =  new TH1F("h_OverlayTight_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_OverlayTight_MET_totPFO =  new TH1F("h_OverlayTight_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayTight_E_totPFO =  new TH1F("h_OverlayTight_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_OverlayTight_E_jetSum =  new TH1F("h_OverlayTight_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_OverlayTight_METPhi_totPFO =  new TH1F("h_OverlayTight_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_OverlayTight_MHT_totPFO =  new TH1F("h_OverlayTight_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayTight_dphi_j1j2_totPFO =  new TH1F("h_OverlayTight_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_OverlayTight_MET_true =  new TH1F("h_OverlayTight_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayTight_METPhi_true =  new TH1F("h_OverlayTight_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_OverlayTight_dphi_MET_true_reco =  new TH1F("h_OverlayTight_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_0_25 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_0_25 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_25_50 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_25_50 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_50_75 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_50_75 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_75_100 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_75_100 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_100_140 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_100_140 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_140_180 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_140_180 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_180_220 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_180_220 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_220_260 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_220_260 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_260_340 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_260_340 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_340_420 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_340_420 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_420_500 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_420_500 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_500_750 =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_500_750 =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayTight_deltaMET_over_trueMET_750_Inf =  new TH1F("h_OverlayTight_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayTight_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_OverlayTight_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_OverlayTight;
  h_vec_OverlayTight.push_back(h_OverlayTight_px_miss_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_delta_MET_totPFO_MET_true);
  h_vec_OverlayTight.push_back(h_OverlayTight_MET_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_E_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_E_jetSum);
  h_vec_OverlayTight.push_back(h_OverlayTight_METPhi_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_MHT_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_j1j2_totPFO);
  h_vec_OverlayTight.push_back(h_OverlayTight_MET_true);
  h_vec_OverlayTight.push_back(h_OverlayTight_METPhi_true);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_MET_true_reco);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_0_25);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_0_25);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_25_50);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_25_50);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_50_75);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_50_75);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_75_100);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_75_100);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_100_140);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_100_140);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_140_180);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_140_180);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_180_220);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_180_220);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_220_260);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_220_260);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_260_340);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_260_340);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_340_420);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_340_420);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_420_500);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_420_500);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_500_750);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_500_750);
  h_vec_OverlayTight.push_back(h_OverlayTight_deltaMET_over_trueMET_750_Inf);
  h_vec_OverlayTight.push_back(h_OverlayTight_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_OverlayTight.size();i++){
    h_vec_OverlayTight[i]->SetLineWidth(3);
    h_vec_OverlayTight[i]->SetLineColor(kGreen+2);
    h_vec_OverlayTight[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_OverlayTight,h_vec_OverlayTight);

  TH1F* h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kGreen+2);
  h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kGreen+2);
  h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_OverlayTight.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayTight[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayTight[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }


  TH1F* h_noOverlay_px_miss_totPFO =  new TH1F("h_noOverlay_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_noOverlay_delta_MET_totPFO_MET_true =  new TH1F("h_noOverlay_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_noOverlay_MET_totPFO =  new TH1F("h_noOverlay_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_noOverlay_E_totPFO =  new TH1F("h_noOverlay_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_noOverlay_E_jetSum =  new TH1F("h_noOverlay_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_noOverlay_METPhi_totPFO =  new TH1F("h_noOverlay_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_noOverlay_MHT_totPFO =  new TH1F("h_noOverlay_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_noOverlay_dphi_j1j2_totPFO =  new TH1F("h_noOverlay_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_noOverlay_MET_true =  new TH1F("h_noOverlay_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_noOverlay_METPhi_true =  new TH1F("h_noOverlay_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_noOverlay_dphi_MET_true_reco =  new TH1F("h_noOverlay_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_0_25 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_0_25 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_25_50 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_25_50 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_50_75 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_50_75 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_75_100 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_75_100 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_100_140 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_100_140 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_140_180 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_140_180 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_180_220 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_180_220 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_220_260 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_220_260 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_260_340 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_260_340 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_340_420 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_340_420 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_420_500 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_420_500 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_500_750 =  new TH1F("h_noOverlay_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_500_750 =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_noOverlay_deltaMET_over_trueMET_750_Inf =  new TH1F("h_noOverlay_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noOverlay_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_noOverlay_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_noOverlay;
  h_vec_noOverlay.push_back(h_noOverlay_px_miss_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_delta_MET_totPFO_MET_true);
  h_vec_noOverlay.push_back(h_noOverlay_MET_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_E_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_E_jetSum);
  h_vec_noOverlay.push_back(h_noOverlay_METPhi_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_MHT_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_j1j2_totPFO);
  h_vec_noOverlay.push_back(h_noOverlay_MET_true);
  h_vec_noOverlay.push_back(h_noOverlay_METPhi_true);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_MET_true_reco);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_0_25);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_0_25);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_25_50);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_25_50);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_50_75);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_50_75);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_75_100);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_75_100);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_100_140);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_100_140);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_140_180);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_140_180);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_180_220);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_180_220);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_220_260);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_220_260);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_260_340);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_260_340);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_340_420);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_340_420);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_420_500);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_420_500);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_500_750);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_500_750);
  h_vec_noOverlay.push_back(h_noOverlay_deltaMET_over_trueMET_750_Inf);
  h_vec_noOverlay.push_back(h_noOverlay_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_noOverlay.size();i++){
    h_vec_noOverlay[i]->SetLineWidth(3);
    h_vec_noOverlay[i]->SetLineColor(kBlack);
    h_vec_noOverlay[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_noOverlay,h_vec_noOverlay);

  TH1F* h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlack);
  h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlack);
  h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_noOverlay.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_noOverlay[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_noOverlay_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_noOverlay[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_noOverlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }



  TH1F* h_OverlayAll_px_miss_totPFO =  new TH1F("h_OverlayAll_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_OverlayAll_delta_MET_totPFO_MET_true =  new TH1F("h_OverlayAll_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_OverlayAll_MET_totPFO =  new TH1F("h_OverlayAll_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayAll_E_totPFO =  new TH1F("h_OverlayAll_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_OverlayAll_E_jetSum =  new TH1F("h_OverlayAll_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_OverlayAll_METPhi_totPFO =  new TH1F("h_OverlayAll_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_OverlayAll_MHT_totPFO =  new TH1F("h_OverlayAll_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayAll_dphi_j1j2_totPFO =  new TH1F("h_OverlayAll_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_OverlayAll_MET_true =  new TH1F("h_OverlayAll_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_OverlayAll_METPhi_true =  new TH1F("h_OverlayAll_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_OverlayAll_dphi_MET_true_reco =  new TH1F("h_OverlayAll_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_0_25 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_0_25 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_25_50 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_25_50 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_50_75 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_50_75 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_75_100 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_75_100 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_100_140 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_100_140 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_140_180 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_140_180 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_180_220 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_180_220 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_220_260 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_220_260 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_260_340 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_260_340 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_340_420 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_340_420 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_420_500 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_420_500 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_500_750 =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_500_750 =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_OverlayAll_deltaMET_over_trueMET_750_Inf =  new TH1F("h_OverlayAll_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_OverlayAll_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_OverlayAll_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_OverlayAll;
  h_vec_OverlayAll.push_back(h_OverlayAll_px_miss_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_delta_MET_totPFO_MET_true);
  h_vec_OverlayAll.push_back(h_OverlayAll_MET_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_E_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_E_jetSum);
  h_vec_OverlayAll.push_back(h_OverlayAll_METPhi_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_MHT_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_j1j2_totPFO);
  h_vec_OverlayAll.push_back(h_OverlayAll_MET_true);
  h_vec_OverlayAll.push_back(h_OverlayAll_METPhi_true);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_MET_true_reco);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_0_25);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_0_25);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_25_50);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_25_50);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_50_75);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_50_75);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_75_100);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_75_100);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_100_140);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_100_140);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_140_180);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_140_180);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_180_220);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_180_220);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_220_260);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_220_260);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_260_340);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_260_340);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_340_420);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_340_420);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_420_500);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_420_500);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_500_750);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_500_750);
  h_vec_OverlayAll.push_back(h_OverlayAll_deltaMET_over_trueMET_750_Inf);
  h_vec_OverlayAll.push_back(h_OverlayAll_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_OverlayAll.size();i++){
    h_vec_OverlayAll[i]->SetLineWidth(3);
    h_vec_OverlayAll[i]->SetLineColor(kMagenta);
    h_vec_OverlayAll[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_OverlayAll,h_vec_OverlayAll);

  TH1F* h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kMagenta);
  h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kMagenta);
  h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  
  for(unsigned int i=11;i<h_vec_OverlayAll.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayAll[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_OverlayAll[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }
  n_bins100_MET=50;

  //3000 overlay
  TH1F* h_Zuds3000_OverlayLoose_px_miss_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_Zuds3000_OverlayLoose_delta_MET_totPFO_MET_true =  new TH1F("h_Zuds3000_OverlayLoose_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_Zuds3000_OverlayLoose_MET_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayLoose_E_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_Zuds3000_OverlayLoose_E_jetSum =  new TH1F("h_Zuds3000_OverlayLoose_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_Zuds3000_OverlayLoose_METPhi_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_Zuds3000_OverlayLoose_MHT_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayLoose_dphi_j1j2_totPFO =  new TH1F("h_Zuds3000_OverlayLoose_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_Zuds3000_OverlayLoose_MET_true =  new TH1F("h_Zuds3000_OverlayLoose_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayLoose_METPhi_true =  new TH1F("h_Zuds3000_OverlayLoose_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_Zuds3000_OverlayLoose_dphi_MET_true_reco =  new TH1F("h_Zuds3000_OverlayLoose_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_Zuds3000_OverlayLoose;
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_px_miss_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_delta_MET_totPFO_MET_true);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_MET_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_E_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_E_jetSum);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_METPhi_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_MHT_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_j1j2_totPFO);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_MET_true);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_METPhi_true);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_MET_true_reco);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_0_25);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_0_25);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_25_50);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_25_50);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_50_75);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_50_75);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_75_100);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_75_100);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_100_140);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_100_140);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_140_180);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_140_180);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_180_220);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_180_220);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_220_260);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_220_260);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_260_340);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_260_340);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_340_420);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_340_420);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_420_500);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_420_500);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_500_750);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_500_750);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_deltaMET_over_trueMET_750_Inf);
  h_vec_Zuds3000_OverlayLoose.push_back(h_Zuds3000_OverlayLoose_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_Zuds3000_OverlayLoose.size();i++){
    h_vec_Zuds3000_OverlayLoose[i]->SetLineWidth(3);
    h_vec_Zuds3000_OverlayLoose[i]->SetLineColor(kBlue);
    h_vec_Zuds3000_OverlayLoose[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_Zuds3000_OverlayLoose,h_vec_Zuds3000_OverlayLoose);

  TH1F* h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlue);
  h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlue);
  h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  
  for(unsigned int i=11;i<h_vec_Zuds3000_OverlayLoose.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayLoose[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_Zuds3000_OverlayLoose_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayLoose[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_Zuds3000_OverlayLoose_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }
  

  TH1F* h_Zuds3000_OverlaySelected_px_miss_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_Zuds3000_OverlaySelected_delta_MET_totPFO_MET_true =  new TH1F("h_Zuds3000_OverlaySelected_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_Zuds3000_OverlaySelected_MET_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlaySelected_E_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_Zuds3000_OverlaySelected_E_jetSum =  new TH1F("h_Zuds3000_OverlaySelected_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_Zuds3000_OverlaySelected_METPhi_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_Zuds3000_OverlaySelected_MHT_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlaySelected_dphi_j1j2_totPFO =  new TH1F("h_Zuds3000_OverlaySelected_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_Zuds3000_OverlaySelected_MET_true =  new TH1F("h_Zuds3000_OverlaySelected_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlaySelected_METPhi_true =  new TH1F("h_Zuds3000_OverlaySelected_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_Zuds3000_OverlaySelected_dphi_MET_true_reco =  new TH1F("h_Zuds3000_OverlaySelected_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_Zuds3000_OverlaySelected;
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_px_miss_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_delta_MET_totPFO_MET_true);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_MET_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_E_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_E_jetSum);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_METPhi_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_MHT_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_j1j2_totPFO);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_MET_true);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_METPhi_true);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_MET_true_reco);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_0_25);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_0_25);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_25_50);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_25_50);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_50_75);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_50_75);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_75_100);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_75_100);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_100_140);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_100_140);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_140_180);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_140_180);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_180_220);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_180_220);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_220_260);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_220_260);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_260_340);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_260_340);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_340_420);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_340_420);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_420_500);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_420_500);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_500_750);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_500_750);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_deltaMET_over_trueMET_750_Inf);
  h_vec_Zuds3000_OverlaySelected.push_back(h_Zuds3000_OverlaySelected_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_Zuds3000_OverlaySelected.size();i++){
    h_vec_Zuds3000_OverlaySelected[i]->SetLineWidth(3);
    h_vec_Zuds3000_OverlaySelected[i]->SetLineColor(kRed-7);
    h_vec_Zuds3000_OverlaySelected[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_Zuds3000_OverlaySelected,h_vec_Zuds3000_OverlaySelected);

  TH1F* h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kRed-7);
  h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kRed-7);
  h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_Zuds3000_OverlaySelected.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlaySelected[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_Zuds3000_OverlaySelected_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlaySelected[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_Zuds3000_OverlaySelected_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }

  TH1F* h_Zuds3000_OverlayTight_px_miss_totPFO =  new TH1F("h_Zuds3000_OverlayTight_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_Zuds3000_OverlayTight_delta_MET_totPFO_MET_true =  new TH1F("h_Zuds3000_OverlayTight_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_Zuds3000_OverlayTight_MET_totPFO =  new TH1F("h_Zuds3000_OverlayTight_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayTight_E_totPFO =  new TH1F("h_Zuds3000_OverlayTight_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_Zuds3000_OverlayTight_E_jetSum =  new TH1F("h_Zuds3000_OverlayTight_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_Zuds3000_OverlayTight_METPhi_totPFO =  new TH1F("h_Zuds3000_OverlayTight_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_Zuds3000_OverlayTight_MHT_totPFO =  new TH1F("h_Zuds3000_OverlayTight_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayTight_dphi_j1j2_totPFO =  new TH1F("h_Zuds3000_OverlayTight_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_Zuds3000_OverlayTight_MET_true =  new TH1F("h_Zuds3000_OverlayTight_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayTight_METPhi_true =  new TH1F("h_Zuds3000_OverlayTight_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_Zuds3000_OverlayTight_dphi_MET_true_reco =  new TH1F("h_Zuds3000_OverlayTight_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayTight_deltaMET_over_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayTight_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_Zuds3000_OverlayTight;
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_px_miss_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_delta_MET_totPFO_MET_true);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_MET_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_E_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_E_jetSum);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_METPhi_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_MHT_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_j1j2_totPFO);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_MET_true);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_METPhi_true);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_MET_true_reco);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_0_25);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_0_25);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_25_50);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_25_50);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_50_75);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_50_75);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_75_100);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_75_100);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_100_140);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_100_140);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_140_180);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_140_180);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_180_220);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_180_220);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_220_260);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_220_260);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_260_340);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_260_340);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_340_420);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_340_420);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_420_500);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_420_500);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_500_750);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_500_750);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_deltaMET_over_trueMET_750_Inf);
  h_vec_Zuds3000_OverlayTight.push_back(h_Zuds3000_OverlayTight_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_Zuds3000_OverlayTight.size();i++){
    h_vec_Zuds3000_OverlayTight[i]->SetLineWidth(3);
    h_vec_Zuds3000_OverlayTight[i]->SetLineColor(kGreen+2);
    h_vec_Zuds3000_OverlayTight[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_Zuds3000_OverlayTight,h_vec_Zuds3000_OverlayTight);

  TH1F* h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kGreen+2);
  h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kGreen+2);
  h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_Zuds3000_OverlayTight.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayTight[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_Zuds3000_OverlayTight_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayTight[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_Zuds3000_OverlayTight_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }


  TH1F* h_noZuds3000_Overlay_px_miss_totPFO =  new TH1F("h_noZuds3000_Overlay_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_noZuds3000_Overlay_delta_MET_totPFO_MET_true =  new TH1F("h_noZuds3000_Overlay_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_noZuds3000_Overlay_MET_totPFO =  new TH1F("h_noZuds3000_Overlay_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_noZuds3000_Overlay_E_totPFO =  new TH1F("h_noZuds3000_Overlay_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_noZuds3000_Overlay_E_jetSum =  new TH1F("h_noZuds3000_Overlay_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_noZuds3000_Overlay_METPhi_totPFO =  new TH1F("h_noZuds3000_Overlay_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_noZuds3000_Overlay_MHT_totPFO =  new TH1F("h_noZuds3000_Overlay_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_noZuds3000_Overlay_dphi_j1j2_totPFO =  new TH1F("h_noZuds3000_Overlay_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_noZuds3000_Overlay_MET_true =  new TH1F("h_noZuds3000_Overlay_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_noZuds3000_Overlay_METPhi_true =  new TH1F("h_noZuds3000_Overlay_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_noZuds3000_Overlay_dphi_MET_true_reco =  new TH1F("h_noZuds3000_Overlay_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_0_25 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_0_25 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_25_50 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_25_50 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_50_75 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_50_75 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_75_100 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_75_100 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_100_140 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_100_140 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_140_180 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_140_180 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_180_220 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_180_220 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_220_260 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_220_260 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_260_340 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_260_340 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_340_420 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_340_420 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_420_500 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_420_500 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_500_750 =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_500_750 =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_noZuds3000_Overlay_deltaMET_over_trueMET_750_Inf =  new TH1F("h_noZuds3000_Overlay_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_noZuds3000_Overlay_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_noZuds3000_Overlay_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_noZuds3000_Overlay;
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_px_miss_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_delta_MET_totPFO_MET_true);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_MET_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_E_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_E_jetSum);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_METPhi_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_MHT_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_j1j2_totPFO);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_MET_true);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_METPhi_true);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_MET_true_reco);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_0_25);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_0_25);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_25_50);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_25_50);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_50_75);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_50_75);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_75_100);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_75_100);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_100_140);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_100_140);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_140_180);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_140_180);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_180_220);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_180_220);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_220_260);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_220_260);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_260_340);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_260_340);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_340_420);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_340_420);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_420_500);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_420_500);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_500_750);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_500_750);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_deltaMET_over_trueMET_750_Inf);
  h_vec_noZuds3000_Overlay.push_back(h_noZuds3000_Overlay_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_noZuds3000_Overlay.size();i++){
    h_vec_noZuds3000_Overlay[i]->SetLineWidth(3);
    h_vec_noZuds3000_Overlay[i]->SetLineColor(kBlack);
    h_vec_noZuds3000_Overlay[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_noZuds3000_Overlay,h_vec_noZuds3000_Overlay);

  TH1F* h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlack);
  h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kBlack);
  h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  for(unsigned int i=11;i<h_vec_noZuds3000_Overlay.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_noZuds3000_Overlay[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_noZuds3000_Overlay_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_noZuds3000_Overlay[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_noZuds3000_Overlay_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }



  TH1F* h_Zuds3000_OverlayAll_px_miss_totPFO =  new TH1F("h_Zuds3000_OverlayAll_px_miss_totPFO","", n_bins100_MET, -250.,250.);
  TH1F* h_Zuds3000_OverlayAll_delta_MET_totPFO_MET_true =  new TH1F("h_Zuds3000_OverlayAll_delta_MET_totPFO","", n_bins100_MET, -100,100);
  TH1F* h_Zuds3000_OverlayAll_MET_totPFO =  new TH1F("h_Zuds3000_OverlayAll_MET_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayAll_E_totPFO =  new TH1F("h_Zuds3000_OverlayAll_E_totPFO","", n_bins100_MET, lim_E_low,lim_E_high);
  TH1F* h_Zuds3000_OverlayAll_E_jetSum =  new TH1F("h_Zuds3000_OverlayAll_E_jetSum","", n_bins100_MET, lim_E_low_jets,lim_E_high_jets);
  TH1F* h_Zuds3000_OverlayAll_METPhi_totPFO =  new TH1F("h_Zuds3000_OverlayAll_METPhi_totPFO","", n_bins100_MET, -180.0,180.);
  TH1F* h_Zuds3000_OverlayAll_MHT_totPFO =  new TH1F("h_Zuds3000_OverlayAll_MHT_totPFO","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayAll_dphi_j1j2_totPFO =  new TH1F("h_Zuds3000_OverlayAll_dphi_j1j2_totPFO","", n_bins100_MET, 160,180);
  TH1F* h_Zuds3000_OverlayAll_MET_true =  new TH1F("h_Zuds3000_OverlayAll_MET_true","", n_bins100_MET, 0,1000);
  TH1F* h_Zuds3000_OverlayAll_METPhi_true =  new TH1F("h_Zuds3000_OverlayAll_METPhi_true","", n_bins100_MET, -180,180);
  TH1F* h_Zuds3000_OverlayAll_dphi_MET_true_reco =  new TH1F("h_Zuds3000_OverlayAll_dphi_MET_true_reco","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_0_25","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_0_25 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_0_25","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_25_50","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_25_50 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_25_50","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_50_75","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_50_75 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_50_75","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_75_100","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_75_100 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_75_100","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_100_140","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_100_140 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_100_140","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_140_180","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_140_180 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_140_180","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_180_220","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_180_220 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_180_220","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_220_260","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_220_260 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_220_260","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_260_340","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_260_340 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_260_340","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_340_420","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_340_420 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_340_420","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_420_500","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_420_500 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_420_500","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_500_750","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_500_750 =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_500_750","", n_bins100_MET, -90.,90.);
  TH1F* h_Zuds3000_OverlayAll_deltaMET_over_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayAll_deltaMET_over_trueMET_750_Inf","", n_bins100_MET, lim_MET_rel_low,lim_MET_rel_high);
  TH1F* h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_750_Inf =  new TH1F("h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_750_Inf","", n_bins100_MET, -90.,90.);

  std::vector<TH1F*>h_vec_Zuds3000_OverlayAll;
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_px_miss_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_delta_MET_totPFO_MET_true);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_MET_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_E_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_E_jetSum);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_METPhi_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_MHT_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_j1j2_totPFO);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_MET_true);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_METPhi_true);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_MET_true_reco);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_0_25);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_0_25);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_25_50);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_25_50);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_50_75);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_50_75);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_75_100);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_75_100);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_100_140);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_100_140);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_140_180);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_140_180);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_180_220);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_180_220);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_220_260);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_220_260);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_260_340);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_260_340);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_340_420);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_340_420);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_420_500);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_420_500);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_500_750);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_500_750);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_deltaMET_over_trueMET_750_Inf);
  h_vec_Zuds3000_OverlayAll.push_back(h_Zuds3000_OverlayAll_dphi_recoMET_trueMET_750_Inf);

  for(unsigned int i=0;i<h_vec_Zuds3000_OverlayAll.size();i++){
    h_vec_Zuds3000_OverlayAll[i]->SetLineWidth(3);
    h_vec_Zuds3000_OverlayAll[i]->SetLineColor(kMagenta);
    h_vec_Zuds3000_OverlayAll[i]->Sumw2();
  }

  fill_background_missingPT_histograms(file_Zuds3000_OverlayAll,h_vec_Zuds3000_OverlayAll);

  TH1F* h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}(#Delta#slash{p}_{T}/#slash{p}_{T,true}) [%]");
  h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kMagenta);
  h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  TH1F* h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram=new TH1F("h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins_MET, pRegionBinEdges_MET);
  h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetYTitle("RMS_{90}#Delta#phi(#vec{#slash{p}}_{T,reco},#vec{#slash{p}}_{T,true})");
  h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetXTitle("#slash{p}_{T,true} [GeV]");
  h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineWidth(3);
  h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineColor(kMagenta);
  h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetLineStyle(1);

  ibin=1;
  
  for(unsigned int i=11;i<h_vec_Zuds3000_OverlayAll.size();i++){
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayAll[i], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    std::cout<<"bin "<<i<<" res/error"<<resolutionTH1<<"/"<<resolutionErrorTH1<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1);
    h_Zuds3000_OverlayAll_METRes_RMS90VsMET_0_1500_summary_histogram ->SetBinError(ibin,resolutionErrorTH1);
    resolutionTH1=0; 
    resolutionErrorTH1=0;
    meanTH1=0;
    meanErrorTH1=0;
    CalculatePerformanceNoMean( h_vec_Zuds3000_OverlayAll[i+1], resolutionTH1, resolutionErrorTH1,meanTH1,meanErrorTH1);
    //RMS is calculated in %, to get the real RMS out divide by 100
    std::cout<<"bin "<<i+1<<" res/error"<<resolutionTH1/100.<<"/"<<resolutionErrorTH1/100.<<" fill "<<ibin<<std::endl;
    h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinContent(ibin,resolutionTH1/100.);
    h_Zuds3000_OverlayAll_PhiRes_RMS90VsMET_0_1500_summary_histogram->SetBinError(ibin,resolutionErrorTH1/100.);
    i++;
    //effectively we always consider a pair of histograms
    ibin+=1;
  }

 
  file_histogram->cd();
  gre_DR07_RMS90_JER_0_70_wSC->Write();
  gre_DR07_RMS90_JER_0_70_noSC->Write();

  file_histogram->Write();
  file_histogram->Close();

}
