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
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TColor.h"
#include "TVector3.h"
#include <vector>
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

using namespace std;

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
    gStyle->SetTitleSize(0.07,"xyz");
    gStyle->SetLabelSize(0.06,"xyz");
    /* title offset: distance between given text and axis, here x,y,z*/
    gStyle->SetLabelOffset(0.015,"xyz");
    gStyle->SetTitleOffset(1.2,"yz"); //equivalent to: gStyle->SetTitleYOffset(1.2);
    gStyle->SetTitleOffset(1.0,"x");
    
    
    
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
    gStyle->SetMarkerSize(1.5);
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
    gStyle->SetFuncColor(kRed);
    gStyle->SetLabelColor(kBlack,"xyz");
    
    //set the margins
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadTopMargin(0.08);
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


float DeltaPhi(float Phi1,float Phi2){
  float deltaphi=fabs(Phi1-Phi2);
  if(deltaphi>M_PI){
    deltaphi=2*M_PI-deltaphi;
  }
  return deltaphi;
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


TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,100,50,690,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}

TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,100,50,690,250);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void AddUnderflowTH1(TH1* hist){
  hist->SetBinContent(1,hist->GetBinContent(1)+hist->GetBinContent(0));
  hist->SetBinError(1,sqrt(pow(hist->GetBinError(1),2)+pow(hist->GetBinError(0),2)));
  hist->SetBinContent(0,0.);
  hist->SetBinError(0,0.);
}
void AddOverUnderflowTH1(TH1* hist){
  hist->SetBinContent(1,hist->GetBinContent(1)+hist->GetBinContent(0));
  hist->SetBinError(1,sqrt(pow(hist->GetBinError(1),2)+pow(hist->GetBinError(0),2)));
  hist->SetBinContent(0,0.);
  hist->SetBinError(0,0.);
  hist->SetBinContent(hist->GetNbinsX(),hist->GetBinContent(hist->GetNbinsX())+hist->GetBinContent(hist->GetNbinsX()+1));
  hist->SetBinError(hist->GetNbinsX(),sqrt(pow(hist->GetBinError(hist->GetNbinsX()),2)+pow(hist->GetBinError(hist->GetNbinsX()+1),2)));
  hist->SetBinContent(hist->GetNbinsX()+1,0.);
  hist->SetBinError(hist->GetNbinsX()+1,0.);
}

//Double_t line_function(Double_t *x, Double_t *par_corr){
//100= par_corr[0]*x[0]+par_corr[1]*x[1];
//}/eos/user/w/weberma2/data/pionMacroFiles/171221//pionStudy


int plotReducedShower(){
  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

 
   TFile* file_photons=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/171221/pionStudy_em5_ILC171221_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune.root");

   //now resolution curves for photons
   //for em 5/10/15/20/25/30/40/50/100/250/500 --> 1000,1500 later
  TFile* file_photons_fitECAL=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/171109/pionStudy_ph1500_ILC171109_CT_FitFW_wRefit_CLIC_o3_v13_SWC_CLIC_K0L_newTune.root");

  //for photon events 1/5/10/15/30/50/100/200/500/1000/1500
  const char* final_histo_name="/afs/cern.ch/user/w/weberma2/validation171212/pionStudy_photon_neutron_ILC171109_CT_FitFW_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_default_finalhistos.root";


  bool is_electron=true;
  bool is_muon=false;
  bool is_photon=false;
  bool is_neutron=false;

  bool run_cluster_hit_part=false;


  float lim_energy_low = 3;
  float lim_energy_high = 7;

  float lim_energy_true_low = 3;
  float lim_energy_true_high = 7;

  int test_signal_particle = 211;
  string label_legend= "e^{-}, E_{true}=5.0 GeV";
  const char*  header_label="e^{-}, E=5.0 GeV";
  const char*  signal_label="pion identified";
  if(is_electron){
    test_signal_particle=11;
    signal_label="electrons";
  }

 if(is_muon){
    test_signal_particle=13;
    signal_label="muons";
  }

 if(is_photon){
   test_signal_particle=22;
   signal_label="photons";
 }
 if(is_neutron){
   test_signal_particle=2112;
   signal_label="neutrons";
 }


  TTree* tree_photons = (TTree*)file_photons->Get("showerData");

  int eventNumber=0;
    
  vector<float> *trueE=0;
  vector<float> *truePx=0;
  vector<float> *truePy=0;
  vector<float> *truePz=0;
  vector<float> *truePhi=0;
  vector<float> *trueCosTheta=0;
  //vector<float> *trueTheta=0;
  vector<int> *truePDGID=0;
  vector<int> *trueGenStatus=0;
  vector<int> *trueNumDaughters=0;
  vector<int> *trueDecayTrackerCalo=0;//1 if decayed in tracker, 2 if decayed in tracker, 3 if coming from backscatter, 4 if left the detector
  vector<int> *trueMotherDecayTrackerCalo=0;//1 if decayed in tracker, 2 if decayed in calo

  vector<float> *recoE=0;
  vector<float> *recoPx=0;
  vector<float> *recoPy=0;
  vector<float> *recoPz=0;
  vector<float> *recoPhi=0;
  vector<float> *recoCosTheta=0;
  //vector<float> *recoTheta=0;
  vector<float> *recoPhi_logERW=0;
  vector<float> *recoCosTheta_logERW=0;
  vector<float> *recoTheta_logERW=0;
  vector<int> *recoCharge=0;
  vector<int> *recoPDGID=0;
  vector<float> *recotrack0p=0;
  vector<float> *recotrack0pt=0;
  vector<float> *recotrack0Chi2OverNdof=0;
  vector<float> *recotrack0nHits=0;
  vector<float> *recoClusterEnergy=0;
  vector<float> *recoEEB=0;
  vector<float> *recoEEE=0;
  vector<float> *recoEEP=0;
  vector<float> *recoEHB=0;
  vector<float> *recoEHE=0;
  vector<float> *recoEHP=0;
  vector<float> *recoEMB=0;
  vector<int> *recoNHitsMB=0;

  vector<int> *recoNhitsEB=0;
  vector<int> *recoNhitsEE=0;
  vector<int> *recoNhitsEP=0;
  vector<int> *recoNhitsHB=0;
  vector<int> *recoNhitsHE=0;
  vector<int> *recoNhitsHP=0;
  vector<int> *recoNtracks=0;
  vector<int> *recoNclusters=0;

  tree_photons->SetBranchAddress("eventNumber",&eventNumber);
  tree_photons->SetBranchAddress("true_Energy",&trueE);
  tree_photons->SetBranchAddress("true_Px",&truePx);
  tree_photons->SetBranchAddress("true_Py",&truePy);
  tree_photons->SetBranchAddress("true_Pz",&truePz);
  tree_photons->SetBranchAddress("true_Phi",&truePhi);
  tree_photons->SetBranchAddress("true_GenStatus",&trueGenStatus);
  tree_photons->SetBranchAddress("true_decayTrackerCalo",&trueDecayTrackerCalo);
  tree_photons->SetBranchAddress("true_motherDecayTrackerCalo",&trueMotherDecayTrackerCalo);
  tree_photons->SetBranchAddress("true_CosTheta",&trueCosTheta);
  //tree_photons->SetBranchAddress("true_Theta",&trueTheta);
  tree_photons->SetBranchAddress("true_PDGID",&truePDGID);
  tree_photons->SetBranchAddress("true_numDaughters",&trueNumDaughters);


  tree_photons->SetBranchAddress("reco_Energy",&recoE);
  tree_photons->SetBranchAddress("reco_Px",&recoPx);
  tree_photons->SetBranchAddress("reco_Py",&recoPy);
  tree_photons->SetBranchAddress("reco_Pz",&recoPz);
  tree_photons->SetBranchAddress("reco_Phi",&recoPhi);
  tree_photons->SetBranchAddress("reco_CosTheta",&recoCosTheta);
  //tree_photons->SetBranchAddress("reco_Theta",&recoTheta);
  tree_photons->SetBranchAddress("reco_Phi_logERW",&recoPhi_logERW);
  tree_photons->SetBranchAddress("reco_CosTheta_logERW",&recoCosTheta_logERW);
  tree_photons->SetBranchAddress("reco_Theta_logERW",&recoTheta_logERW);
  tree_photons->SetBranchAddress("reco_PDGID",&recoPDGID);
  tree_photons->SetBranchAddress("reco_Charge",&recoCharge);
  tree_photons->SetBranchAddress("reco_track0_p",&recotrack0p);
  tree_photons->SetBranchAddress("reco_track0_pt",&recotrack0pt);
  tree_photons->SetBranchAddress("reco_track0_chi2OverNdof",&recotrack0Chi2OverNdof);
  tree_photons->SetBranchAddress("reco_track0_nHits",&recotrack0nHits);
  tree_photons->SetBranchAddress("reco_clusters_energy",&recoClusterEnergy);
  tree_photons->SetBranchAddress("reco_nClusters",&recoNclusters);
  tree_photons->SetBranchAddress("reco_nTracks",&recoNtracks);

  tree_photons->SetBranchAddress("reco_E_EB",&recoEEB);
  tree_photons->SetBranchAddress("reco_E_EE",&recoEEE);
  tree_photons->SetBranchAddress("reco_E_EP",&recoEEP);
  tree_photons->SetBranchAddress("reco_E_HB",&recoEHB);
  tree_photons->SetBranchAddress("reco_E_HE",&recoEHE);
  tree_photons->SetBranchAddress("reco_E_HP",&recoEHP);
  tree_photons->SetBranchAddress("reco_E_MB",&recoEMB);
  tree_photons->SetBranchAddress("reco_nhitsEB",&recoNhitsEB);
  tree_photons->SetBranchAddress("reco_nhitsEE",&recoNhitsEE);
  tree_photons->SetBranchAddress("reco_nhitsEP",&recoNhitsEP);
  tree_photons->SetBranchAddress("reco_nhitsHB",&recoNhitsHB);
  tree_photons->SetBranchAddress("reco_nhitsHE",&recoNhitsHE);
  tree_photons->SetBranchAddress("reco_nhitsHP",&recoNhitsHP);
  tree_photons->SetBranchAddress("reco_nhitsMB",&recoNHitsMB);

  int n_bins100 = 49;
  float lim_phi_low = -TMath::Pi();
  float lim_phi_high = TMath::Pi();
  //file_histogram->Open();



  TH1F* h_phi_pions =  new TH1F("h_phi_pions","", n_bins100, lim_phi_low,lim_phi_high);
  h_phi_pions->Sumw2();
  h_phi_pions->SetLineWidth(2);

  TH1F* h_phi_else =  new TH1F("h_phi_else","", n_bins100, lim_phi_low,lim_phi_high);
  h_phi_else->Sumw2();
  h_phi_else->SetLineWidth(2);
  h_phi_else->SetLineColor(2);

  float lim_cosTheta_low = -1.0;
  float lim_cosTheta_high = 1.0;
    
    float lim_Theta_low = 0;
    float lim_Theta_high = 25;

  TH1F* h_nHitsMB_pionInEvent =  new TH1F("h_nHitsMB_pionInEvent","", 11, -0.5,10.5);
  h_nHitsMB_pionInEvent->Sumw2();
  h_nHitsMB_pionInEvent->SetLineWidth(2);

  TH1F* h_nHitsMB_noPionInEvent =  new TH1F("h_nHitsMB_noPionInEvent","", 11, -0.5,10.5);
  h_nHitsMB_noPionInEvent->Sumw2();
  h_nHitsMB_noPionInEvent->SetLineWidth(2);
  h_nHitsMB_noPionInEvent->SetLineColor(2);

  TH1F* h_E_MB_pionInEvent =  new TH1F("h_E_MB_pionInEvent","", 25, -0,2.2);
  h_E_MB_pionInEvent->Sumw2();
  h_E_MB_pionInEvent->SetLineWidth(2);

  TH1F* h_E_MB_noPionInEvent =  new TH1F("h_E_MB_noPionInEvent","", 25, -0,2.2);
  h_E_MB_noPionInEvent->Sumw2();
  h_E_MB_noPionInEvent->SetLineWidth(2);
  h_E_MB_noPionInEvent->SetLineColor(2);

  TH1F* h_cosTheta_pions =  new TH1F("h_cosTheta_pions","", n_bins100, lim_cosTheta_low,lim_cosTheta_high);
  h_cosTheta_pions->Sumw2();
  h_cosTheta_pions->SetLineWidth(2);

    TH1F* h_Theta_pions =  new TH1F("h_Theta_pions","", n_bins100, lim_Theta_low,lim_Theta_high);
    h_Theta_pions->Sumw2();
    h_Theta_pions->SetLineWidth(2);

  TH1F* h_cosTheta_else =  new TH1F("h_cosTheta_else","", n_bins100, lim_cosTheta_low,lim_cosTheta_high);
  h_cosTheta_else->Sumw2();
  h_cosTheta_else->SetLineWidth(2);
  h_cosTheta_else->SetLineColor(2);
    
    TH1F* h_Theta_else =  new TH1F("h_Theta_else","", n_bins100, lim_Theta_low,lim_Theta_high);
    h_Theta_else->Sumw2();
    h_Theta_else->SetLineWidth(2);
    h_Theta_else->SetLineColor(2);


  TH1F* h_cluster1_E_pionInEvent =  new TH1F("h_cluster1_E_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_pionInEvent->Sumw2();
  h_cluster1_E_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster1_E_noPionInEvent =  new TH1F("h_cluster1_E_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_noPionInEvent->Sumw2();
  h_cluster1_E_noPionInEvent->SetLineWidth(2);
  h_cluster1_E_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster1_E_EB_pionInEvent =  new TH1F("h_cluster1_E_EB_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_EB_pionInEvent->Sumw2();
  h_cluster1_E_EB_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster1_E_EB_noPionInEvent =  new TH1F("h_cluster1_E_EB_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_EB_noPionInEvent->Sumw2();
  h_cluster1_E_EB_noPionInEvent->SetLineWidth(2);
  h_cluster1_E_EB_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster1_E_EE_pionInEvent =  new TH1F("h_cluster1_E_EE_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_EE_pionInEvent->Sumw2();
  h_cluster1_E_EE_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster1_E_EE_noPionInEvent =  new TH1F("h_cluster1_E_EE_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_EE_noPionInEvent->Sumw2();
  h_cluster1_E_EE_noPionInEvent->SetLineWidth(2);
  h_cluster1_E_EE_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster1_E_HB_pionInEvent =  new TH1F("h_cluster1_E_HB_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_HB_pionInEvent->Sumw2();
  h_cluster1_E_HB_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster1_E_HB_noPionInEvent =  new TH1F("h_cluster1_E_HB_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_HB_noPionInEvent->Sumw2();
  h_cluster1_E_HB_noPionInEvent->SetLineWidth(2);
  h_cluster1_E_HB_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster1_E_HE_pionInEvent =  new TH1F("h_cluster1_E_HE_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_HE_pionInEvent->Sumw2();
  h_cluster1_E_HE_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster1_E_HE_noPionInEvent =  new TH1F("h_cluster1_E_HE_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster1_E_HE_noPionInEvent->Sumw2();
  h_cluster1_E_HE_noPionInEvent->SetLineWidth(2);
  h_cluster1_E_HE_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster2_E_pionInEvent =  new TH1F("h_cluster2_E_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_pionInEvent->Sumw2();
  h_cluster2_E_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster2_E_noPionInEvent =  new TH1F("h_cluster2_E_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_noPionInEvent->Sumw2();
  h_cluster2_E_noPionInEvent->SetLineWidth(2);
  h_cluster2_E_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster2_E_EB_pionInEvent =  new TH1F("h_cluster2_E_EB_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_EB_pionInEvent->Sumw2();
  h_cluster2_E_EB_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster2_E_EB_noPionInEvent =  new TH1F("h_cluster2_E_EB_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_EB_noPionInEvent->Sumw2();
  h_cluster2_E_EB_noPionInEvent->SetLineWidth(2);
  h_cluster2_E_EB_noPionInEvent->SetLineColor(2);

  TH2F* h_all_hits_Z_vs_R =  new TH2F("h_all_hits_Z_vs_R","caloHits, no signal identified", 140, -250,250, 100,4400,6450);
  h_all_hits_Z_vs_R->Sumw2();
  h_all_hits_Z_vs_R->GetXaxis()->SetTitle("|Z_{hit}|");
  h_all_hits_Z_vs_R->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_all_hits_Z_vs_R_noPionInEvent =  new TH2F("h_all_hits_Z_vs_R_noPionInEvent ","caloHits, no signal identified", 140, -250,250, 100,4400,6450);
  h_all_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_all_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_all_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_all_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_all_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_all_hits_Z_vs_R_pionInEvent =  new TH2F("h_all_hits_Z_vs_R_pionInEvent ","calohits, signal identified", 140, -250,250, 100,4400,6450);
  h_all_hits_Z_vs_R_pionInEvent->Sumw2();
  h_all_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_all_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_all_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_all_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_no_cluster_hits_Z_vs_R_noPionInEvent =  new TH2F("h_no_cluster_hits_Z_vs_R_noPionInEvent ","caloHits, no signal identified", 140, -250,250, 100,4400,6450);
  h_no_cluster_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_no_cluster_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_no_cluster_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_no_cluster_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_no_cluster_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_no_cluster_hits_Z_vs_R_pionInEvent =  new TH2F("h_no_cluster_hits_Z_vs_R_pionInEvent ","calohits, signal identified", 140, -250,250, 100,4400,6450);
  h_no_cluster_hits_Z_vs_R_pionInEvent->Sumw2();
  h_no_cluster_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_no_cluster_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_no_cluster_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_no_cluster_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_cluster_hits_Z_vs_R_noPionInEvent =  new TH2F("h_cluster_hits_Z_vs_R_noPionInEvent ","caloHits, no signal identified", 140, -250,250, 100,4400,6450);
  h_cluster_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_cluster_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_cluster_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_cluster_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_cluster_hits_Z_vs_R_pionInEvent =  new TH2F("h_cluster_hits_Z_vs_R_pionInEvent ","calohits, signal identified", 140, -250,250, 100,4400,6450);
  h_cluster_hits_Z_vs_R_pionInEvent->Sumw2();
  h_cluster_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_cluster_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_cluster_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_cluster1_hits_Z_vs_R_pionInEvent =  new TH2F("h_cluster1_hits_Z_vs_R_pionInEvent ","cluster1 calohits, signal identified", 140, -250,250, 100,4400,6450);
  h_cluster1_hits_Z_vs_R_pionInEvent->Sumw2();
  h_cluster1_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_cluster1_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_cluster1_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster1_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_cluster1_hits_Z_vs_R_noPionInEvent =  new TH2F("h_cluster1_hits_Z_vs_R_noPionInEvent ","cluster1 calohits, no signal identified", 140, -250,250, 100,4400,6450);
  h_cluster1_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_cluster1_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_cluster1_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_cluster1_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster1_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_RECO_cluster_hits_Z_vs_R_pionInEvent =  new TH2F("h_RECO_cluster_hits_Z_vs_R_pionInEvent ","calohits of signal particle", 140, -250,250, 100,4400,6450);
  h_RECO_cluster_hits_Z_vs_R_pionInEvent->Sumw2();
  h_RECO_cluster_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_RECO_cluster_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_RECO_cluster_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_RECO_cluster_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_RECO_cluster_hits_Z_vs_R_noPionInEvent =  new TH2F("h_RECO_cluster_hits_Z_vs_R_noPionInEvent ","", 140, -250,250, 100,4400,6450);
  h_RECO_cluster_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_RECO_cluster_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_RECO_cluster_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_RECO_cluster_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_RECO_cluster_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");


  TH2F* h_cluster2_hits_Z_vs_R_pionInEvent =  new TH2F("h_cluster2_hits_Z_vs_R_pionInEvent ","cluster2 calohits, signal identified", 140, 2050,2250, 100,1400,1680);
  h_cluster2_hits_Z_vs_R_pionInEvent->Sumw2();
  h_cluster2_hits_Z_vs_R_pionInEvent->SetLineWidth(2);
  h_cluster2_hits_Z_vs_R_pionInEvent->SetLineColor(2);
  h_cluster2_hits_Z_vs_R_pionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster2_hits_Z_vs_R_pionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH2F* h_cluster2_hits_Z_vs_R_noPionInEvent =  new TH2F("h_cluster2_hits_Z_vs_R_noPionInEvent ","cluster 2 calohits, no signal identified", 140, 2050,2250, 100,1400,1680);
  h_cluster2_hits_Z_vs_R_noPionInEvent->Sumw2();
  h_cluster2_hits_Z_vs_R_noPionInEvent->SetLineWidth(2);
  h_cluster2_hits_Z_vs_R_noPionInEvent->SetLineColor(2);
  h_cluster2_hits_Z_vs_R_noPionInEvent->GetXaxis()->SetTitle("|Z_{hit}|");
  h_cluster2_hits_Z_vs_R_noPionInEvent->GetYaxis()->SetTitle("R_{hit}");

  TH1F* h_cluster2_E_EE_pionInEvent =  new TH1F("h_cluster2_E_EE_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_EE_pionInEvent->Sumw2();
  h_cluster2_E_EE_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster2_E_EE_noPionInEvent =  new TH1F("h_cluster2_E_EE_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_EE_noPionInEvent->Sumw2();
  h_cluster2_E_EE_noPionInEvent->SetLineWidth(2);
  h_cluster2_E_EE_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster2_E_HB_pionInEvent =  new TH1F("h_cluster2_E_HB_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_HB_pionInEvent->Sumw2();
  h_cluster2_E_HB_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster2_E_HB_noPionInEvent =  new TH1F("h_cluster2_E_HB_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_HB_noPionInEvent->Sumw2();
  h_cluster2_E_HB_noPionInEvent->SetLineWidth(2);
  h_cluster2_E_HB_noPionInEvent->SetLineColor(2);

  TH1F* h_cluster2_E_HE_pionInEvent =  new TH1F("h_cluster2_E_HE_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_HE_pionInEvent->Sumw2();
  h_cluster2_E_HE_pionInEvent->SetLineWidth(2);

  TH1F* h_cluster2_E_HE_noPionInEvent =  new TH1F("h_cluster2_E_HE_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_cluster2_E_HE_noPionInEvent->Sumw2();
  h_cluster2_E_HE_noPionInEvent->SetLineWidth(2);
  h_cluster2_E_HE_noPionInEvent->SetLineColor(2);


  TH1F* h_energy_pions =  new TH1F("h_energy_pions","", n_bins100, lim_energy_low,lim_energy_high);
  h_energy_pions->Sumw2();
  h_energy_pions->SetLineWidth(2);

  TH1F* h_energy_else_pionInEvent =  new TH1F("h_energy_else_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_energy_else_pionInEvent->Sumw2();
  h_energy_else_pionInEvent->SetLineWidth(2);
  h_energy_else_pionInEvent->SetLineColor(2);

  TH1F* h_energy_true_ph_pionInEvent =  new TH1F("h_energy_true_ph_pionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_ph_pionInEvent->Sumw2();
  h_energy_true_ph_pionInEvent->SetLineWidth(2);

  TH1F* h_energy_true_e_pionInEvent =  new TH1F("h_energy_true_e_pionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_e_pionInEvent->Sumw2();
  h_energy_true_e_pionInEvent->SetLineWidth(2);
  h_energy_true_e_pionInEvent->SetLineColor(2);

  TH1F* h_energy_true_sum_pionInEvent =  new TH1F("h_energy_true_sum_pionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_sum_pionInEvent->Sumw2();
  h_energy_true_sum_pionInEvent->SetLineWidth(2);
  h_energy_true_sum_pionInEvent->SetLineColor(kGreen+2);

  TH1F* h_energy_else_noPionInEvent =  new TH1F("h_energy_else_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_energy_else_noPionInEvent->Sumw2();
  h_energy_else_noPionInEvent->SetLineWidth(2);
  h_energy_else_noPionInEvent->SetLineColor(2);

  TH1F* h_energy_true_ph_noPionInEvent =  new TH1F("h_energy_true_ph_noPionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_ph_noPionInEvent->Sumw2();
  h_energy_true_ph_noPionInEvent->SetLineWidth(2);

  TH1F* h_energy_true_e_noPionInEvent =  new TH1F("h_energy_true_e_noPionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_e_noPionInEvent->Sumw2();
  h_energy_true_e_noPionInEvent->SetLineWidth(2);
  h_energy_true_e_noPionInEvent->SetLineColor(2);

  TH1F* h_energy_true_sum_noPionInEvent =  new TH1F("h_energy_true_sum_noPionInEvent","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_sum_noPionInEvent->Sumw2();
  h_energy_true_sum_noPionInEvent->SetLineWidth(2);
  h_energy_true_sum_noPionInEvent->SetLineColor(kGreen+2);

  TH1F* h_energy_else =  new TH1F("h_energy_else","", n_bins100, lim_energy_low,lim_energy_high);
  h_energy_else->Sumw2();
  h_energy_else->SetLineWidth(2);
  h_energy_else->SetLineColor(2);

  TH1F* h_energy_true_ph =  new TH1F("h_energy_true_ph","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_ph->Sumw2();
  h_energy_true_ph->SetLineWidth(2);

  TH1F* h_energy_true_e =  new TH1F("h_energy_true_e","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_e->Sumw2();
  h_energy_true_e->SetLineWidth(2);
  h_energy_true_e->SetLineColor(2);

  TH1F* h_energy_true_sum =  new TH1F("h_energy_true_sum","", n_bins100, lim_energy_true_low,lim_energy_true_high);
  h_energy_true_sum->Sumw2();
  h_energy_true_sum->SetLineWidth(2);
  h_energy_true_sum->SetLineColor(kGreen+2);

  TH1F* h_energyTot_pionInEvent =  new TH1F("h_energyTot_pionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_energyTot_pionInEvent->Sumw2();
  h_energyTot_pionInEvent->SetLineWidth(2);
 
  TH1F* h_energyTot_noPionInEvent =  new TH1F("h_energyTot_noPionInEvent","", n_bins100, lim_energy_low,lim_energy_high);
  h_energyTot_noPionInEvent->Sumw2();
  h_energyTot_noPionInEvent->SetLineWidth(2);
  h_energyTot_noPionInEvent->SetLineColor(2);

  TH1F* h_clusterESum_pionInEvent =  new TH1F("h_clusterEsum_pionInEvent","", n_bins100, -0.1,lim_energy_high);
  h_clusterESum_pionInEvent->Sumw2();
  h_clusterESum_pionInEvent->SetLineWidth(2);
 
  TH1F* h_clusterESum_noPionInEvent =  new TH1F("h_clusterESum_noPionInEvent","", n_bins100, -0.1,lim_energy_high);
  h_clusterESum_noPionInEvent->Sumw2();
  h_clusterESum_noPionInEvent->SetLineWidth(2);
  h_clusterESum_noPionInEvent->SetLineColor(2);

  TH1F* h_nTot_pionInEvent =  new TH1F("h_nPFOs_pionInEvent","", 10, 0.5,10.5);
  h_nTot_pionInEvent->Sumw2();
  h_nTot_pionInEvent->SetLineWidth(2);

  TH1F* h_nTot_noPionInEvent =  new TH1F("h_nPFOs_noPionInEvent","", 10, 0.5,10.5);
  h_nTot_noPionInEvent->Sumw2();
  h_nTot_noPionInEvent->SetLineWidth(2);
  h_nTot_noPionInEvent->SetLineColor(2);


  TH1F* h_nClusters_pionInEvent =  new TH1F("h_nClusters_pionInEvent","", 11, -0.5,10.5);
  h_nClusters_pionInEvent->Sumw2();
  h_nClusters_pionInEvent->SetLineWidth(2);
 
  TH1F* h_nClusters_noPionInEvent =  new TH1F("h_nClusters_noPionInEvent","", 11, -0.5,10.5);
  h_nClusters_noPionInEvent->Sumw2();
  h_nClusters_noPionInEvent->SetLineWidth(2);
  h_nClusters_noPionInEvent->SetLineColor(2);

  

  TH1F* h_truePhi_pionInEvent =  new TH1F("h_truePhi_pionInEvent","", n_bins100, lim_phi_low,lim_phi_high);
  h_truePhi_pionInEvent->Sumw2();
  h_truePhi_pionInEvent->SetLineWidth(2);
 
  TH1F* h_truePhi_noPionInEvent =  new TH1F("h_truePhi_noPionInEvent","", n_bins100, lim_phi_low,lim_phi_high);
  h_truePhi_noPionInEvent->Sumw2();
  h_truePhi_noPionInEvent->SetLineWidth(2);
  h_truePhi_noPionInEvent->SetLineColor(2);

  TH1F* h_truePhi_all =  new TH1F("h_truePhi_all","", n_bins100, lim_phi_low,lim_phi_high);
  h_truePhi_all->Sumw2();
  h_truePhi_all->SetLineWidth(2);
  h_truePhi_all->SetLineColor(2);

  TH1F* h_trueCosTheta_pionInEvent =  new TH1F("h_trueCosTheta_pionInEvent","", n_bins100, lim_cosTheta_low,lim_cosTheta_high);
  h_trueCosTheta_pionInEvent->Sumw2();
  h_trueCosTheta_pionInEvent->SetLineWidth(2);
 
  TH1F* h_trueCosTheta_noPionInEvent =  new TH1F("h_trueCosTheta_noPionInEvent","", n_bins100, lim_cosTheta_low,lim_cosTheta_high);
  h_trueCosTheta_noPionInEvent->Sumw2();
  h_trueCosTheta_noPionInEvent->SetLineWidth(2);
  h_trueCosTheta_noPionInEvent->SetLineColor(2);

  TH1F* h_trueCosTheta_all =  new TH1F("h_trueCosTheta_all","", n_bins100, lim_cosTheta_low,lim_cosTheta_high);
  h_trueCosTheta_all->Sumw2();
  h_trueCosTheta_all->SetLineWidth(2);
  h_trueCosTheta_all->SetLineColor(2);
  
    TH1F* h_trueTheta_pionInEvent =  new TH1F("h_trueTheta_pionInEvent","", n_bins100, lim_Theta_low,lim_Theta_high);
    h_trueTheta_pionInEvent->Sumw2();
    h_trueTheta_pionInEvent->SetLineWidth(2);
    
    TH1F* h_trueTheta_noPionInEvent =  new TH1F("h_trueTheta_noPionInEvent","", n_bins100, lim_Theta_low,lim_Theta_high);
    h_trueTheta_noPionInEvent->Sumw2();
    h_trueTheta_noPionInEvent->SetLineWidth(2);
    h_trueTheta_noPionInEvent->SetLineColor(2);
    
    TH1F* h_trueTheta_all =  new TH1F("h_trueTheta_all","", n_bins100, lim_Theta_low,lim_Theta_high);
    h_trueTheta_all->Sumw2();
    h_trueTheta_all->SetLineWidth(2);
    h_trueTheta_all->SetLineColor(2);
    
  TH1F* h_dEPFO12_pionInEvent =  new TH1F("h_dEPFO12_pionInEvent","", n_bins100, 0.5,lim_energy_high);
  h_dEPFO12_pionInEvent->Sumw2();
  h_dEPFO12_pionInEvent->SetLineWidth(2);

  TH1F* h_dEPFO12_noPionInEvent =  new TH1F("h_dEPFO12_noPionInEvent","", n_bins100,0.5,lim_energy_high);
  h_dEPFO12_noPionInEvent->Sumw2();
  h_dEPFO12_noPionInEvent->SetLineWidth(2);
  h_dEPFO12_noPionInEvent->SetLineColor(2);
  h_dEPFO12_noPionInEvent->SetLineStyle(2);

  TH1F* h_dphiPFO12_pionInEvent =  new TH1F("h_dphiPFO12_pionInEvent","", n_bins100, -0.7,0.5);
  h_dphiPFO12_pionInEvent->Sumw2();
  h_dphiPFO12_pionInEvent->SetLineWidth(2);

  TH1F* h_dphiPFO12_noPionInEvent =  new TH1F("h_dphiPFO12_noPionInEvent","", n_bins100,-0.7,0.5);
  h_dphiPFO12_noPionInEvent->Sumw2();
  h_dphiPFO12_noPionInEvent->SetLineWidth(2);
  h_dphiPFO12_noPionInEvent->SetLineColor(2);
  h_dphiPFO12_noPionInEvent->SetLineStyle(2);

  TH1F* h_dthetaPFO12_pionInEvent =  new TH1F("h_dthetaPFO12_pionInEvent","", n_bins100, -0.5,0.5);
  h_dthetaPFO12_pionInEvent->Sumw2();
  h_dthetaPFO12_pionInEvent->SetLineWidth(2);

  TH1F* h_dthetaPFO12_noPionInEvent =  new TH1F("h_dthetaPFO12_noPionInEvent","", n_bins100,-0.5,0.5);
  h_dthetaPFO12_noPionInEvent->Sumw2();
  h_dthetaPFO12_noPionInEvent->SetLineWidth(2);
  h_dthetaPFO12_noPionInEvent->SetLineColor(2);
  h_dthetaPFO12_noPionInEvent->SetLineStyle(2);


  TH1F* h_PFOPDGID_pionInEvent =  new TH1F("h_PFOPDGID_pionInEvent","", 3000, -249.5,2750.5);
  h_PFOPDGID_pionInEvent->Sumw2();
  h_PFOPDGID_pionInEvent->SetLineWidth(2);

  TH1F* h_PFOPDGID_noPionInEvent =  new TH1F("h_PFOPDGID_noPionInEvent","", 3000, -249.5,2750.5);
  h_PFOPDGID_noPionInEvent->Sumw2();
  h_PFOPDGID_noPionInEvent->SetLineWidth(2);
  h_PFOPDGID_noPionInEvent->SetLineColor(2);
  h_PFOPDGID_noPionInEvent->SetLineStyle(2);
  
  float lim_dtheta_low=-0.004;
  float lim_dtheta_high=0.004;

  TH1F* h_PFOSignal_deltaTheta_Barrel =  new TH1F("h_PFOSignal_deltaTheta_Barrel","", n_bins100, lim_dtheta_low,lim_dtheta_high);
  h_PFOSignal_deltaTheta_Barrel->Sumw2();
  h_PFOSignal_deltaTheta_Barrel->SetLineWidth(2);

  TH1F* h_PFOSignal_deltaTheta_Endcap =  new TH1F("h_PFOSignal_deltaTheta_Endcap","", n_bins100, lim_dtheta_low,lim_dtheta_high);
  h_PFOSignal_deltaTheta_Endcap->Sumw2();
  h_PFOSignal_deltaTheta_Endcap->SetLineWidth(2);
  h_PFOSignal_deltaTheta_Endcap->SetLineColor(2);
  h_PFOSignal_deltaTheta_Endcap->SetLineStyle(2);

  TH1F* h_PFOSignal_deltaPhi_Barrel =  new TH1F("h_PFOSignal_deltaPhi_Barrel","", n_bins100, lim_dtheta_low,lim_dtheta_high);
  h_PFOSignal_deltaPhi_Barrel->Sumw2();
  h_PFOSignal_deltaPhi_Barrel->SetLineWidth(2);

 TH1F* h_PFOSignal_deltaPhi_Endcap =  new TH1F("h_PFOSignal_deltaPhi_Endcap","", n_bins100, lim_dtheta_low,lim_dtheta_high);
  h_PFOSignal_deltaPhi_Endcap->Sumw2();
  h_PFOSignal_deltaPhi_Endcap->SetLineWidth(2);
  h_PFOSignal_deltaPhi_Endcap->SetLineColor(2);
  h_PFOSignal_deltaPhi_Endcap->SetLineStyle(2);


    TEfficiency* SignalIDEfficiencyVsTrueTheta = new TEfficiency("SignalIDEffVsMCTrueTheta","",75,0.,25);
    SignalIDEfficiencyVsTrueTheta->SetTitle("signal identification efficiency, e^{-}, E=3.0 GeV;#theta_{true};");
    
  TEfficiency* SignalIDEfficiencyVsTrueCosTheta = new TEfficiency("SignalIDEffVsMCTrueCosTheta","",75,-1.0,1.0);
  SignalIDEfficiencyVsTrueCosTheta->SetTitle("signal identification efficiency, e^{-}, E=3.0 GeV;cos#theta_{true};");

  TEfficiency* muonEfficiencyVsTrueCosTheta_noPionInEvent = new TEfficiency("muonEffVsMCTrueCosTheta_noSignal","",75,-1.0,1.0);
  muonEfficiencyVsTrueCosTheta_noPionInEvent->SetTitle("muon efficiency on missidentified pion events, e^{-}, E=3.0 GeV;cos#theta_{true};");
  TEfficiency* muonMuonEfficiencyVsTrueCosTheta_noPionInEvent = new TEfficiency("muonMuonEffVsMCTrueCosTheta_noSignal","",75,-1.0,1.0);
  muonMuonEfficiencyVsTrueCosTheta_noPionInEvent->SetTitle("muon or muon ID efficiency on missidentified pion events, e^{-}, E=3.0 GeV;cos#theta_{true};");

 TEfficiency* PionVetoEfficiencyVsTrueCosTheta_noPionInEvent = new TEfficiency("PionVetoVsMCTrueCosTheta_noSignal","",75,-1.0,1.0);
 PionVetoEfficiencyVsTrueCosTheta_noPionInEvent->SetTitle("pion Veto Efficiency on missidentified muon events, e^{-}, E=3.0 GeV;cos#theta_{true};");

 TEfficiency* PionEfficiencyVsTrueCosTheta_noPionInEvent = new TEfficiency("PionVsMCTrueCosTheta_noSignal","",75,-1.0,1.0);
 PionEfficiencyVsTrueCosTheta_noPionInEvent->SetTitle("pion Efficiency on missidentified muon events, e^{-}, E=3.0 GeV;cos#theta_{true};");

TEfficiency* PionEfficiencyVsTrueTheta_noPionInEvent = new TEfficiency("PionVsMCTrueTheta_noSignal","",75,0.,25);
PionEfficiencyVsTrueTheta_noPionInEvent->SetTitle("pion Efficiency on missidentified muon events, e^{-}, E=3.0 GeV;#theta_{true};");

  unsigned int signal_events=0;
  unsigned int no_signal_events=0;

  //for(unsigned int i_entry=0;i_entry<750;i_entry++){
  for(unsigned int i_entry=0;i_entry<tree_photons->GetEntries();i_entry++){
    tree_photons->GetEntry(i_entry);

    if(i_entry%250==0){
      std::cout<<"in entry "<<i_entry<<" test for "<<test_signal_particle<<std::endl;
    }

    bool veto_transition=false;
    bool fill_true_phi_values=false;
    bool fill_reco_phi_values=false;
    bool fill_true_CosTheta_values=false;
    bool fill_reco_CosTheta_values=false;
    for(unsigned int i=0;i<trueGenStatus->size();i++){
      if((*trueGenStatus)[i]==1 && ( (acos((*trueCosTheta)[i])*TMath::RadToDeg())<10 || ((acos((*trueCosTheta)[i])*TMath::RadToDeg())>20 && (acos((*trueCosTheta)[i])*TMath::RadToDeg())<160) || (acos((*trueCosTheta)[i])*TMath::RadToDeg())>170)){//TR default 80 %
	veto_transition=true;
	//break;
      }
    }
    if(veto_transition){
      continue;
    }

    //distance hit implement in pion shower macro
    bool signal_event=false;
    float energy_tot=0;
    unsigned int nClusters_tot=0;
    bool pass_highE_extras=false;
    float extras_E=0;    
    TLorentzVector TVtrue;
    bool found_truth=false;
    
    float E_true_else=0;
    float E_true_signal=0;
    
    bool decayintracker=false;

    bool foundCluster=false;//found one cluster

    bool reachesecal=true;


    unsigned int index_true=0;
    
    for(unsigned int m=0;m<trueGenStatus->size();m++){
      if((*trueGenStatus)[m]==1 && (*trueDecayTrackerCalo)[m]==1 ){
	//std::cout<<"decay in tracker "<<eventNumber<<std::endl;
	decayintracker=true;
      }
      //if((*truePDGID)[m]==22 && (*trueE)[m]>1){
      //pass_highE_photon=true;
      //}
      //std::cout<<"test values PDG/GenStatus/decayTC/motDecTC "<<(*truePDGID)[m]<<"/"<<(*trueDecayTrackerCalo)[m]<<"/"<<(*trueMotherDecayTrackerCalo)[m]<<std::endl;
      if(abs((*truePDGID)[m])==test_signal_particle && (*trueGenStatus)[m]==0){
	//particles which are created in the detector, but also don't decay inside the tracker volume
	E_true_signal+=(*trueE)[m];
      }else if((*trueGenStatus)[m]==0){//particles which are created in the detector, but also don't decay inside the tracker volume, no signal
	E_true_else+=(*trueE)[m];
      }
      if(abs((*truePDGID)[m])==test_signal_particle && (*trueGenStatus)[m]==1){//the status gen 1 particle
	index_true=m;
	TVtrue.SetPxPyPzE((*truePx)[m],(*truePy)[m],(*truePz)[m],(*trueE)[m]);
	h_trueCosTheta_all->Fill((*trueCosTheta)[m]);
    h_trueTheta_all->Fill(acos(fabs((*trueCosTheta)[m]))*TMath::RadToDeg());
	h_truePhi_all->Fill((*truePhi)[m]);
	TVector3 TV3true((*truePx)[m],(*truePy)[m],(*truePz)[m]); 
	found_truth=true;
      }
      if((*trueGenStatus)[m]!=1){
	extras_E+=(*trueE)[m];
      }
    }//true stuff done
    h_energy_true_ph->Fill(E_true_else);
    h_energy_true_e->Fill(E_true_signal);
    h_energy_true_sum->Fill(E_true_else+E_true_signal);
    if(found_truth==false){
      std::cout<<"WTF< where is truth"<<std::endl;
    }
    if(truePhi->size()==1){
    }
    TLorentzVector TVreco;

    int index0=-1;
    int index1=-1;
    
    float E_index0=-1;
    float E_index1=-1;
    bool isElectronEvent=false;
    bool isMuonEvent=false;
    bool isPionEvent=false;
    bool isPhotonEvent=false;
    int index_signal=-1;
    for(unsigned int s=0;s<recoE->size();s++){
      if((*recoE)[s]>E_index0){
          E_index1=E_index0;
          E_index0=(*recoE)[s];
          index1=index0;
          index0=s;
      }else if ((*recoE)[s]>E_index1){
          E_index1=(*recoE)[s];
          index1=s;
      }
      nClusters_tot+=(*recoNclusters)[s];
      if(abs((*recoPDGID)[s])==211){
          isPionEvent=true;
	//std::cout<<"pion found "<<std::endl;
      }
      if(abs((*recoPDGID)[s])==22){
          isPhotonEvent=true;
      }
      if(abs((*recoPDGID)[s])==11){
          isElectronEvent=true;
	//std::cout<<"muon found "<<std::endl;
      }
      if(abs((*recoPDGID)[s])==13){
          isMuonEvent=true;
	//std::cout<<"muon found "<<std::endl;
      }
      if(abs((*recoPDGID)[s])==test_signal_particle){
          TVreco.SetPxPyPzE((*recoPx)[s],(*recoPy)[s],(*recoPz)[s],(*recoE)[s]);
          if(index_signal==-1){
              if(fabs((*recoE)[s]-(*trueE)[index_true])<fabs((*recoE)[index_signal]-(*trueE)[index_true])){
                  index_signal=s;
              }
          }
          if(fabs((*trueCosTheta)[index_true])<0.7){
              h_PFOSignal_deltaTheta_Barrel->Fill(acos((*recoCosTheta)[s])-acos((*trueCosTheta)[index_true]));
              h_PFOSignal_deltaPhi_Barrel->Fill(DeltaPhiDir((*recoPhi)[s],(*truePhi)[index_true]));
          }else{
              h_PFOSignal_deltaTheta_Endcap->Fill(acos((*recoCosTheta)[s])-acos((*trueCosTheta)[index_true]));
              h_PFOSignal_deltaPhi_Endcap->Fill(DeltaPhiDir((*recoPhi)[s],(*truePhi)[index_true]));
          }
          double dotproduct3D=(TVreco.Px()*TVtrue.Px()+TVreco.Py()*TVtrue.Py()+TVreco.Pz()*TVtrue.Pz())/(TVreco.P()*TVtrue.P());
	  //0.2 corresponds to about 11.46 degrees
          if((acos(dotproduct3D)*TMath::RadToDeg())<10.5 /*&& (fabs(TVtrue.Energy()-TVreco.Energy())<(100.*0.15*sqrt(TVtrue.Energy())))*/){
              //std::cout<<"rat1/rat2 "<<fabs(TVtrue.Energy()-TVreco.Energy())/sqrt(TVtrue.Energy())<<"/"<< fabs(TVtrue.Energy()-TVreco.Energy())/TVtrue.Energy()<<std::endl;
              signal_event=true;
          }
          h_energy_pions->Fill((*recoE)[s]);
          h_phi_pions->Fill((*recoPhi)[s]);
          h_cosTheta_pions->Fill((*recoCosTheta)[s]);
          h_Theta_pions->Fill(acos(fabs((*recoCosTheta)[s]))*TMath::RadToDeg());
          }else{
              h_energy_else->Fill((*recoE)[s]);
              h_phi_else->Fill((*recoPhi)[s]);
              h_cosTheta_else->Fill((*recoCosTheta)[s]);
              h_Theta_else->Fill(acos(fabs((*recoCosTheta)[s]))*TMath::RadToDeg());
          }
        energy_tot+=(*recoE)[s];
        }//reco energy done
    SignalIDEfficiencyVsTrueCosTheta->Fill(signal_event, (*trueCosTheta)[index_true]);
    SignalIDEfficiencyVsTrueTheta->Fill(signal_event, acos(fabs((*trueCosTheta)[index_true]))*TMath::RadToDeg());
    if(signal_event){
      h_nHitsMB_pionInEvent->Fill((*recoNHitsMB)[0]);
      h_E_MB_pionInEvent->Fill((*recoEMB)[0]);
      signal_events+=1;
      if(index1>-1){
	h_dphiPFO12_pionInEvent->Fill(DeltaPhiDir((*recoPhi)[index0],(*recoPhi)[index1]));
	h_dthetaPFO12_pionInEvent->Fill(acos((*recoCosTheta)[index0])-acos((*recoCosTheta)[index1]));
	h_dEPFO12_pionInEvent->Fill((*recoE)[index0]-(*recoE)[index1]);
      }
      h_energy_true_ph_pionInEvent->Fill(E_true_else);
      h_energy_true_e_pionInEvent->Fill(E_true_signal);
      h_energy_true_sum_pionInEvent->Fill(E_true_else+E_true_signal);
      h_energyTot_pionInEvent->Fill(energy_tot);
      h_nTot_pionInEvent->Fill(recoE->size());
      for(unsigned int s=0;s<recoE->size();s++){
	h_PFOPDGID_pionInEvent->Fill((*recoPDGID)[s]);
      }      
      for(unsigned int s=0;s<trueE->size();s++){
	if(abs((*truePDGID)[s])==test_signal_particle && (*trueGenStatus)[s]==1){
	  h_trueCosTheta_pionInEvent->Fill((*trueCosTheta)[s]);
      h_trueTheta_pionInEvent->Fill(acos(fabs((*trueCosTheta)[s]))*TMath::RadToDeg());
	  h_truePhi_pionInEvent->Fill((*truePhi)[s]);
	}
      }
    }else{
      h_nHitsMB_noPionInEvent->Fill((*recoNHitsMB)[0]);
      h_E_MB_noPionInEvent->Fill((*recoEMB)[0]);
      no_signal_events+=1;
      muonEfficiencyVsTrueCosTheta_noPionInEvent->Fill((isElectronEvent),(*trueCosTheta)[index_true]);
      muonMuonEfficiencyVsTrueCosTheta_noPionInEvent->Fill((isElectronEvent || isMuonEvent),(*trueCosTheta)[index_true]);
      PionVetoEfficiencyVsTrueCosTheta_noPionInEvent->Fill((!isPionEvent),(*trueCosTheta)[index_true]);
      PionEfficiencyVsTrueCosTheta_noPionInEvent->Fill((isPionEvent),(*trueCosTheta)[index_true]);
        PionEfficiencyVsTrueTheta_noPionInEvent->Fill((isPionEvent),acos(fabs((*trueCosTheta)[index_true]))*TMath::RadToDeg());
      if(index1>-1){
	h_dphiPFO12_noPionInEvent->Fill(DeltaPhiDir((*recoPhi)[index0],(*recoPhi)[index1]));
	h_dthetaPFO12_noPionInEvent->Fill(acos((*recoCosTheta)[index0])- acos((*recoCosTheta)[index1]));
	h_dEPFO12_noPionInEvent->Fill((*recoE)[index0]-(*recoE)[index1]);
      }
      h_energy_true_ph_noPionInEvent->Fill(E_true_else);
      h_energy_true_e_noPionInEvent->Fill(E_true_signal);
      h_energy_true_sum_noPionInEvent->Fill(E_true_else+E_true_signal);
      h_energyTot_noPionInEvent->Fill(energy_tot);
      h_nTot_noPionInEvent->Fill(recoE->size());
      for(unsigned int s=0;s<trueE->size();s++){
	if(abs((*truePDGID)[s])==test_signal_particle && (*trueGenStatus)[s]==1){
	  h_trueCosTheta_noPionInEvent->Fill((*trueCosTheta)[s]);
      h_trueTheta_noPionInEvent->Fill(acos(fabs((*trueCosTheta)[s]))*TMath::RadToDeg());
	  h_truePhi_noPionInEvent->Fill((*truePhi)[s]);
	}
      }
      for(unsigned int s=0;s<recoE->size();s++){
	h_PFOPDGID_noPionInEvent->Fill((*recoPDGID)[s]);
      }
    }
  }
  
  TCanvas* can_E_particles=setUpperCanvas("can_E_particles");
  can_E_particles->cd();

  TLegend* leg_E_particles=new TLegend(0.45,0.65,0.89,0.89);
  leg_E_particles->SetBorderSize(0);
  leg_E_particles->SetFillStyle(0);
  leg_E_particles->SetHeader(header_label);
  leg_E_particles->AddEntry(h_energy_pions->DrawCopy("hist,e"),signal_label);
  leg_E_particles->AddEntry(h_energy_else->DrawCopy("hist,e,same"),"else");
  leg_E_particles->Draw();
    
    //h_energy_pions->Fit("gaus");
    //h_energy_pions->Draw("p,e,same");
    
    
    
    //TF1* fit_energy_resolution =  (TF1*)h_energy_pions->GetFunction("gaus");
    //std::cout<<"fit sigma over fit mean "<<fit_energy_resolution->GetParameter(2)<<"/"<<fit_energy_resolution->GetParameter(1)<<"/"<<fit_energy_resolution->GetParameter(2)/fit_energy_resolution->GetParameter(1)<<std::endl;
  
 std::cout<<"particle E signal"<<h_energy_pions->GetMean()/10.<<"/"<<h_energy_pions->GetMeanError()/10.<<"/"<<h_energy_pions->GetRMS()/10.<<"/"<<h_energy_pions->GetRMSError()/10.<<std::endl;

  TCanvas* can_E_tot_event=setUpperCanvas("can_E_tot_event");
  can_E_tot_event->cd();

  TLegend* leg_E_tot_event=new TLegend(0.45,0.65,0.89,0.89);
  leg_E_tot_event->SetBorderSize(0);
  leg_E_tot_event->SetFillStyle(0);
  leg_E_tot_event->SetHeader(header_label);
  leg_E_tot_event->AddEntry(h_energyTot_pionInEvent->DrawCopy("hist,e"),signal_label);
  //leg_E_tot_event->AddEntry(h_energyTot_noPionInEvent->DrawCopy("hist,e,same"),"no signal particle in event");
  leg_E_tot_event->Draw();
  std::cout<<"particle E tot"<<h_energyTot_pionInEvent->GetMean()/10.<<"/"<<h_energyTot_pionInEvent->GetMeanError()/10.<<"/"<<h_energyTot_pionInEvent->GetRMS()/10.<<"/"<<h_energyTot_pionInEvent->GetRMSError()/10.<<std::endl;

  
  
  TCanvas* can_phi_particle=setUpperCanvas("can_phi_particle");
  can_phi_particle->cd();

  TLegend* leg_phi_particles=new TLegend(0.45,0.65,0.89,0.89);
  leg_phi_particles->SetBorderSize(0);
  leg_phi_particles->SetFillStyle(0);
  leg_phi_particles->SetHeader(header_label);
  leg_phi_particles->AddEntry(h_phi_pions->DrawCopy("hist,e"),signal_label);
  leg_phi_particles->AddEntry(h_phi_else->DrawCopy("hist,e,same"),"else");
  leg_phi_particles->Draw();
  
  TCanvas* can_cosTheta_particle=setUpperCanvas("can_cosTheta_particle");
  can_cosTheta_particle->cd();

  TLegend* leg_cosTheta_particles=new TLegend(0.45,0.65,0.89,0.89);
  leg_cosTheta_particles->SetBorderSize(0);
  leg_cosTheta_particles->SetFillStyle(0);
  leg_cosTheta_particles->SetHeader(header_label);
  leg_cosTheta_particles->AddEntry(h_cosTheta_pions->DrawCopy("hist,e"),signal_label);
  leg_cosTheta_particles->AddEntry(h_cosTheta_else->DrawCopy("hist,e,same"),"else");
  leg_cosTheta_particles->Draw();
  
  TCanvas* can_nPFOs_particle=setUpperCanvas("can_nPFOs_event");
  can_nPFOs_particle->cd();

  TLegend* leg_nPFOs_particles=new TLegend(0.45,0.65,0.89,0.89);
  leg_nPFOs_particles->SetBorderSize(0);
  leg_nPFOs_particles->SetFillStyle(0);
  leg_nPFOs_particles->SetHeader(header_label);
  leg_nPFOs_particles->AddEntry(h_nTot_pionInEvent->DrawCopy("hist,e"),signal_label);
  leg_nPFOs_particles->AddEntry(h_nTot_noPionInEvent->DrawCopy("hist,e,same"),"else");
  leg_nPFOs_particles->Draw();

  TCanvas* can_nClusters_particle=setUpperCanvas("can_nClusters_event");
  can_nClusters_particle->cd();

  TLegend* leg_nClusters_particles=new TLegend(0.45,0.65,0.89,0.89);
  leg_nClusters_particles->SetBorderSize(0);
  leg_nClusters_particles->SetFillStyle(0);
  leg_nClusters_particles->AddEntry(h_nClusters_pionInEvent->DrawCopy("hist,e"),signal_label);
  leg_nClusters_particles->AddEntry(h_nClusters_noPionInEvent->DrawCopy("hist,e,same"),"else");
  leg_nClusters_particles->Draw();

  
  TCanvas* can_PFOPDGID=setUpperCanvas("can_PFOPDGID_event");
  can_PFOPDGID->cd();
  
  TLegend* leg_PFOPDGIDs=new TLegend(0.45,0.65,0.89,0.89);
  leg_PFOPDGIDs->SetBorderSize(0);
  leg_PFOPDGIDs->SetFillStyle(0);
  leg_PFOPDGIDs->SetHeader(header_label);
  leg_PFOPDGIDs->AddEntry(h_PFOPDGID_pionInEvent->DrawCopy("hist,e"),"PFO-PDGIDs, signal found in event");
  leg_PFOPDGIDs->AddEntry(h_PFOPDGID_noPionInEvent->DrawCopy("hist,e,same"),"PFO-PDGIDs, no signal in event");
  leg_PFOPDGIDs->Draw();

  TCanvas* can_truePhi_event=setUpperCanvas("can_truePhi_event");
  can_truePhi_event->cd();

  TLegend* leg_truePhi_events=new TLegend(0.45,0.65,0.89,0.89);
  leg_truePhi_events->SetBorderSize(0);
  leg_truePhi_events->SetFillStyle(0);
  leg_truePhi_events->AddEntry(h_truePhi_pionInEvent->DrawCopy("hist,e"),signal_label);
  leg_truePhi_events->AddEntry(h_truePhi_noPionInEvent->DrawCopy("hist,e,same"),"else");
  leg_truePhi_events->Draw();

  TCanvas* can_trueCosTheta_event=setUpperCanvas("can_trueCosTheta_event");
  can_trueCosTheta_event->cd();

  TLegend* leg_trueCosTheta_events=new TLegend(0.45,0.65,0.89,0.89);
  leg_trueCosTheta_events->SetBorderSize(0);
  leg_trueCosTheta_events->SetFillStyle(0);
  leg_trueCosTheta_events->SetHeader(header_label);
  leg_trueCosTheta_events->AddEntry(h_trueCosTheta_pionInEvent->DrawCopy("hist,e"),signal_label);
  leg_trueCosTheta_events->AddEntry(h_trueCosTheta_noPionInEvent->DrawCopy("hist,e,same"),"else");
  leg_trueCosTheta_events->Draw();
    
    TCanvas* can_trueTheta_event=setUpperCanvas("can_trueTheta_event");
    can_trueTheta_event->cd();
    
    TLegend* leg_trueTheta_events=new TLegend(0.45,0.65,0.89,0.89);
    leg_trueTheta_events->SetBorderSize(0);
    leg_trueTheta_events->SetFillStyle(0);
    leg_trueTheta_events->SetHeader(header_label);
    leg_trueTheta_events->AddEntry(h_trueTheta_pionInEvent->DrawCopy("hist,e"),signal_label);
    leg_trueTheta_events->AddEntry(h_trueTheta_noPionInEvent->DrawCopy("hist,e,same"),"else");
    leg_trueTheta_events->Draw();

  
  TCanvas* can_deltaPhi_Signal=setUpperCanvas("can_deltaPhi_Signal");
  can_deltaPhi_Signal->cd();

  TLegend* leg_deltaPhi_Signals=new TLegend(0.45,0.65,0.89,0.89);
  leg_deltaPhi_Signals->SetBorderSize(0);
  leg_deltaPhi_Signals->SetFillStyle(0);
  leg_deltaPhi_Signals->SetHeader("e^{-}, E=3.0 GeV");
  leg_deltaPhi_Signals->AddEntry(h_PFOSignal_deltaPhi_Barrel->DrawCopy("hist,e"),"|cos(#theta_{true})|<0.7");
  leg_deltaPhi_Signals->AddEntry(h_PFOSignal_deltaPhi_Endcap->DrawCopy("hist,e,same"),"|cos(#theta_{true})|#geq0.7");
  leg_deltaPhi_Signals->Draw();

  TCanvas* can_deltaTheta_Signal=setUpperCanvas("can_deltaTheta_Signal");
  can_deltaTheta_Signal->cd();

  TLegend* leg_deltaTheta_Signals=new TLegend(0.45,0.65,0.89,0.89);
  leg_deltaTheta_Signals->SetBorderSize(0);
  leg_deltaTheta_Signals->SetFillStyle(0);
  leg_deltaTheta_Signals->SetHeader("e^{-}, E=3.0 GeV");
  leg_deltaTheta_Signals->AddEntry(h_PFOSignal_deltaTheta_Barrel->DrawCopy("hist,e"),"|cos(#theta_{true})|<0.7");
  leg_deltaTheta_Signals->AddEntry(h_PFOSignal_deltaTheta_Endcap->DrawCopy("hist,e,same"),"|cos(#theta_{true})|#geq0.7");
  leg_deltaTheta_Signals->Draw();

  
  TEfficiency* tEff_SignalPion_CosTheta = new TEfficiency(*h_trueCosTheta_pionInEvent,*h_trueCosTheta_all);
  tEff_SignalPion_CosTheta->SetTitle("pion identification efficiency, E=30 GeV; cos#theta_{true};");
    
    TEfficiency* tEff_SignalPion_Theta = new TEfficiency(*h_trueTheta_pionInEvent,*h_trueTheta_all);
    tEff_SignalPion_Theta->SetTitle("pion identification efficiency, E=30 GeV; #theta_{true};");

  TEfficiency* tEff_SignalPion_Phi = new TEfficiency(*h_truePhi_pionInEvent,*h_truePhi_all);
  tEff_SignalPion_Phi->SetTitle("signal particle identification efficiency; #Phi_{true}; ");

 


  TCanvas* can_eff_CosTheta_Signal=setUpperCanvas("can_eff_CosTheta_Signal");
  can_eff_CosTheta_Signal->cd();

  tEff_SignalPion_CosTheta->Draw();

  TLegend* leg_eff_CosTheta_Signal=new TLegend(0.45,0.65,0.89,0.89);
  leg_eff_CosTheta_Signal->SetBorderSize(0);
  leg_eff_CosTheta_Signal->SetFillStyle(0);
  leg_eff_CosTheta_Signal->SetHeader("muon identification efficiency, E=3.0 GeV");
  leg_eff_CosTheta_Signal->Draw();
    
    TCanvas* can_eff_Theta_Signal=setUpperCanvas("can_eff_Theta_Signal");
    can_eff_Theta_Signal->cd();
    
    tEff_SignalPion_Theta->Draw();
    
    TLegend* leg_eff_Theta_Signal=new TLegend(0.45,0.65,0.89,0.89);
    leg_eff_Theta_Signal->SetBorderSize(0);
    leg_eff_Theta_Signal->SetFillStyle(0);
    //leg_eff_Theta_Signal->SetHeader("pion identification efficiency, ILCSoft-17-08-23, E=5 GeV");
    leg_eff_Theta_Signal->Draw();
    /*
  
  TCanvas* can_eff_phi_Signal=setUpperCanvas("can_eff_phi_Signal");
  can_eff_phi_Signal->cd();

  tEff_SignalPion_Phi->Draw();

  TLegend* leg_eff_phi_Signal=new TLegend(0.45,0.65,0.89,0.89);
  leg_eff_phi_Signal->SetBorderSize(0);
  leg_eff_phi_Signal->SetFillStyle(0);
  leg_eff_phi_Signal->SetHeader("Signal identification efficiency, TruthTrack FitFW");
  leg_eff_phi_Signal->Draw();
  

  TCanvas* can_trueE_sums=setUpperCanvas("can_trueE_sums");
  can_trueE_sums->cd();

  TLegend* leg_trueE_sums=new TLegend(0.45,0.65,0.89,0.89);
  leg_trueE_sums->SetBorderSize(0);
  leg_trueE_sums->SetFillStyle(0);
  leg_trueE_sums->SetHeader("PFA e^{-},E=5 GeV, True Particle Energy distribution");
  leg_trueE_sums->AddEntry(h_energy_true_e->DrawCopy("hist,e"),signal_label);
  leg_trueE_sums->AddEntry(h_energy_true_ph->DrawCopy("hist,e,same"),"other particles");
  leg_trueE_sums->AddEntry(h_energy_true_sum->DrawCopy("hist,e,same"),"all particles");
  leg_trueE_sums->Draw();


  TCanvas* can_trueE_sums_noPionInEvent=setUpperCanvas("can_trueE_sums_noPionInEvent");
  can_trueE_sums_noPionInEvent->cd();

  TLegend* leg_trueE_sums_noPionInEvent=new TLegend(0.45,0.65,0.89,0.89);
  leg_trueE_sums_noPionInEvent->SetBorderSize(0);
  leg_trueE_sums_noPionInEvent->SetFillStyle(0);
  leg_trueE_sums_noPionInEvent->SetHeader("PFA e^{-},E=5 GeV, True Particle Energy (no signal identified)","h");
  leg_trueE_sums_noPionInEvent->AddEntry(h_energy_true_e_noPionInEvent->DrawCopy("hist,e"),signal_label);
  leg_trueE_sums_noPionInEvent->AddEntry(h_energy_true_ph_noPionInEvent->DrawCopy("hist,e,same"),"other particles");
  leg_trueE_sums_noPionInEvent->AddEntry(h_energy_true_sum_noPionInEvent->DrawCopy("hist,e,same"),"all particles");
  leg_trueE_sums_noPionInEvent->Draw();


  TCanvas* can_trueE_sums_pionInEvent=setUpperCanvas("can_trueE_sums_pionInEvent");
  can_trueE_sums_pionInEvent->cd();

  TLegend* leg_trueE_sums_pionInEvent=new TLegend(0.45,0.65,0.89,0.89);
  leg_trueE_sums_pionInEvent->SetBorderSize(0);
  leg_trueE_sums_pionInEvent->SetFillStyle(0);
  leg_trueE_sums_pionInEvent->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_trueE_sums_pionInEvent->AddEntry(h_energy_true_e_pionInEvent->DrawCopy("hist,e"),signal_label);
  leg_trueE_sums_pionInEvent->AddEntry(h_energy_true_ph_pionInEvent->DrawCopy("hist,e,same"),"other particles");
  leg_trueE_sums_pionInEvent->AddEntry(h_energy_true_sum_pionInEvent->DrawCopy("hist,e,same"),"all particles");
  leg_trueE_sums_pionInEvent->Draw();




  TCanvas* can_dphiPFO12=setUpperCanvas("can_dphiPFO12");
  can_dphiPFO12->cd();

  TLegend* leg_dphiPFO12=new TLegend(0.45,0.65,0.89,0.89);
  leg_dphiPFO12->SetBorderSize(0);
  leg_dphiPFO12->SetFillStyle(0);
  leg_dphiPFO12->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_dphiPFO12->AddEntry(h_dphiPFO12_pionInEvent->DrawCopy("hist,e"),"signal identified");
  leg_dphiPFO12->AddEntry(h_dphiPFO12_noPionInEvent->DrawCopy("hist,e,same"),"signal not identified");
  leg_dphiPFO12->Draw();

  TCanvas* can_dEPFO12=setUpperCanvas("can_dEPFO12");
  can_dEPFO12->cd();

  TLegend* leg_dEPFO12=new TLegend(0.45,0.65,0.89,0.89);
  leg_dEPFO12->SetBorderSize(0);
  leg_dEPFO12->SetFillStyle(0);
  leg_dEPFO12->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_dEPFO12->AddEntry(h_dEPFO12_pionInEvent->DrawCopy("hist,e"),"signal identified");
  leg_dEPFO12->AddEntry(h_dEPFO12_noPionInEvent->DrawCopy("hist,e,same"),"signal not identified");
  leg_dEPFO12->Draw();

  TCanvas* can_dthetaPFO12=setUpperCanvas("can_dthetaPFO12");
  can_dthetaPFO12->cd();

  TLegend* leg_dthetaPFO12=new TLegend(0.45,0.65,0.89,0.89);
  leg_dthetaPFO12->SetBorderSize(0);
  leg_dthetaPFO12->SetFillStyle(0);
  leg_dthetaPFO12->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_dthetaPFO12->AddEntry(h_dthetaPFO12_pionInEvent->DrawCopy("hist,e"),"signal identified");
  leg_dthetaPFO12->AddEntry(h_dthetaPFO12_noPionInEvent->DrawCopy("hist,e,same"),"all tracks, signal not IDd");
  leg_dthetaPFO12->Draw();

    */
  /*
  TCanvas* can_PiVeto_eff_vs_trueCosTheta_noPionInEvent=setUpperCanvas("can_PiVeto_eff_vs_trueCosTheta_noPionInEvent");
  can_PiVeto_eff_vs_trueCosTheta_noPionInEvent->cd();
  PionVetoEfficiencyVsTrueCosTheta_noPionInEvent->Draw();

  TCanvas* can_Pi_eff_vs_trueCosTheta_noPionInEvent=setUpperCanvas("can_Pi_eff_vs_trueCosTheta_noPionInEvent");
  can_Pi_eff_vs_trueCosTheta_noPionInEvent->cd();
  PionEfficiencyVsTrueCosTheta_noPionInEvent->Draw();

  TCanvas* can_El_eff_vs_trueCosTheta_noElonInEvent=setUpperCanvas("can_El_eff_vs_trueCosTheta_noElonInEvent");
  can_El_eff_vs_trueCosTheta_noElonInEvent->cd();
  muonEfficiencyVsTrueCosTheta_noPionInEvent->Draw();

  TCanvas* can_ElMu_eff_vs_trueCosTheta_noElMuonInEvent=setUpperCanvas("can_ElMu_eff_vs_trueCosTheta_noElMuonInEvent");
  can_ElMu_eff_vs_trueCosTheta_noElMuonInEvent->cd();
  muonMuonEfficiencyVsTrueCosTheta_noPionInEvent->Draw();
  */
    /*
  TCanvas* can_SignalID_eff_vs_trueCosTheta=setUpperCanvas("can_SignalID_eff_vs_trueCosTheta");
  can_SignalID_eff_vs_trueCosTheta->cd();
  SignalIDEfficiencyVsTrueCosTheta->Draw();
    
    std::cout<<"signal id efficiency"<< SignalIDEfficiencyVsTrueCosTheta->GetEfficiency(50)<<std::endl;
  
  

  TCanvas* can_reco0_E_MB=setUpperCanvas("can_reco0_E_MB");
  can_reco0_E_MB->cd();

  TLegend* leg_reco0_E_MB=new TLegend(0.45,0.65,0.89,0.89);
  leg_reco0_E_MB->SetBorderSize(0);
  leg_reco0_E_MB->SetFillStyle(0);
  leg_reco0_E_MB->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_reco0_E_MB->AddEntry(h_E_MB_pionInEvent->DrawCopy("hist,e"),"signal identified");
  leg_reco0_E_MB->AddEntry(h_E_MB_noPionInEvent->DrawCopy("hist,e,same"),"signal not identified");
  leg_reco0_E_MB->Draw();

  TCanvas* can_reco0_nHitsMB=setUpperCanvas("can_reco0_nHitsMB");
  can_reco0_nHitsMB->cd();

  TLegend* leg_reco0_nHitsMB=new TLegend(0.45,0.65,0.89,0.89);
  leg_reco0_nHitsMB->SetBorderSize(0);
  leg_reco0_nHitsMB->SetFillStyle(0);
  leg_reco0_nHitsMB->SetHeader("PFA e^{-},E=5 GeV","h");
  leg_reco0_nHitsMB->AddEntry(h_nHitsMB_pionInEvent->DrawCopy("hist,e"),"signal identified");
  leg_reco0_nHitsMB->AddEntry(h_nHitsMB_noPionInEvent->DrawCopy("hist,e,same"),"signal not identified");
  leg_reco0_nHitsMB->Draw();

    TTree* tree_SWC_CLIC_n_K0L= (TTree*)file_photons_fitECAL->Get("showerData");
    
    vector<float> *true_Energy_SWC_CLIC_n_K0L=0;
    vector<float> *true_Px_SWC_CLIC_n_K0L=0;
    vector<float> *true_Py_SWC_CLIC_n_K0L=0;
    vector<float> *true_Pz_SWC_CLIC_n_K0L=0;
    vector<float> *true_CosTheta_SWC_CLIC_n_K0L=0;
    vector<int> *true_PDGID_SWC_CLIC_n_K0L=0;
    vector<int> *true_GenStatus_SWC_CLIC_n_K0L=0;
    
    vector<float> *true_ConvVtx_x=0;
    
    vector<float> *reco_Energy_SWC_CLIC_n_K0L=0;
    vector<float> *reco_Px_SWC_CLIC_n_K0L=0;
    vector<float> *reco_Py_SWC_CLIC_n_K0L=0;
    vector<float> *reco_Pz_SWC_CLIC_n_K0L=0;
    vector<int> *reco_PDGID_SWC_CLIC_n_K0L=0;
    
    vector<float> *reco_E_EB=0;
    vector<float> *reco_E_EE=0;
    vector<float> *reco_E_HB=0;
    vector<float> *reco_E_HE=0;
    
    vector<float> *cluster_Energy_SWC_CLIC_n_K0L=0;
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_Energy", &true_Energy_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_Px", &true_Px_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_Py", &true_Py_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_Pz", &true_Pz_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_CosTheta", &true_CosTheta_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_PDGID", &true_PDGID_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("true_GenStatus", &true_GenStatus_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("trueConvVtxx", &true_ConvVtx_x);
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Energy", &reco_Energy_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Px", &reco_Px_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Py", &reco_Py_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Pz", &reco_Pz_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_PDGID", &reco_PDGID_SWC_CLIC_n_K0L);
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_EB", &reco_E_EB);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_EE", &reco_E_EE);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_HB", &reco_E_HB);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_HE", &reco_E_HE);
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("cluster_energy", &cluster_Energy_SWC_CLIC_n_K0L);
    
    float lim_energy_low_rel=1.0;
    float lim_energy_high_rel=1.075;
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_cosThetaTR =new TH1F("h_SWC_CLIC_n_K0L_cosThetaTR","", n_bins100, 0, 1.0);
    h_SWC_CLIC_n_K0L_cosThetaTR ->Sumw2();
    h_SWC_CLIC_n_K0L_cosThetaTR ->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_cosThetaTR ->GetXaxis()->SetTitle("|cos#theta|");
    
    
    float true_energy_fillplot=0;
    float true_CosTheta_min =  1.5;
    float true_CosTheta_max = -1.5;
    unsigned int count_hadron_candidate_events=0;
    for (unsigned int i_entry=0;i_entry<tree_SWC_CLIC_n_K0L->GetEntries();i_entry++){
        tree_SWC_CLIC_n_K0L->GetEntry(i_entry);
        bool pass_hadron_true=false;
        float true_CosTheta=-1.5;
        float true_energy_hadron=-1;
        unsigned int true_hadron_counter=0;
        TLorentzVector HadTrue(0,0,0,0);
        if(true_ConvVtx_x->size()>0){
            //std::cout<<"entry "<<i_entry<<" is conversion "<<(*true_ConvVtx_x)[0]<<"/"<<true_ConvVtx_x->size()<<std::endl;
            continue;
        }
        //std::cout<<"for entry "<<i_entry<<" after continue"<<std::endl;
        for(unsigned int i_true=0;i_true<true_Energy_SWC_CLIC_n_K0L->size();i_true++){
            if((*true_GenStatus_SWC_CLIC_n_K0L)[i_true]==1 && abs((*true_PDGID_SWC_CLIC_n_K0L)[i_true])==22){
                HadTrue.SetPxPyPzE((*true_Px_SWC_CLIC_n_K0L)[i_true],(*true_Py_SWC_CLIC_n_K0L)[i_true],(*true_Pz_SWC_CLIC_n_K0L)[i_true],(*true_Energy_SWC_CLIC_n_K0L)[i_true]);
                true_energy_hadron=(*true_Energy_SWC_CLIC_n_K0L)[i_true];
                true_hadron_counter+=1;
                true_CosTheta=(*true_CosTheta_SWC_CLIC_n_K0L)[i_true];
                true_energy_fillplot=true_energy_hadron;
            }
        }
        if(true_hadron_counter!=1){
            std::cout<<"WTF not enough hadrons "<<true_hadron_counter<<std::endl;
        }else{
            //std::cout<<"find true hadron "<<true_CosTheta<<"/"<<true_energy_hadron<<std::endl;
        }
        bool has_E_Endcap=false;
        bool has_E_Barrel=false;
        float E_barrel=0;
        float E_endcap=0;
        float E_reco_total=0;
        bool find_reco_candidate=false;
        float E_had_candidate_max=0;
        TLorentzVector HadReco(0,0,0,0);
        for(unsigned int i_reco=0;i_reco<reco_Energy_SWC_CLIC_n_K0L->size();i_reco++){
            E_reco_total+=(*reco_Energy_SWC_CLIC_n_K0L)[i_reco];
            if((*reco_PDGID_SWC_CLIC_n_K0L)[i_reco]==22){//check for neutral hadrons ONLY
                TLorentzVector temp(0,0,0,0);
                temp.SetPxPyPzE((*reco_Px_SWC_CLIC_n_K0L)[i_reco],(*reco_Py_SWC_CLIC_n_K0L)[i_reco],(*reco_Pz_SWC_CLIC_n_K0L)[i_reco],(*reco_Energy_SWC_CLIC_n_K0L)[i_reco]);
                double cosAngcheck=((HadTrue.Px()*temp.Px()+HadTrue.Py()*temp.Py()+HadTrue.Pz()*temp.Pz())/(HadTrue.P()*temp.P()));
                if(cosAngcheck>0.99999999994){//allow 8 % offset for neutral hadron candidate)
                    if(temp.Energy()>E_had_candidate_max){
                        E_had_candidate_max=temp.Energy();
                        HadReco=temp;
                        std::cout<<" i get to this loop "<<E_had_candidate_max<<"/"<<std::endl;
                    }
                }
            }
            if((*reco_E_EB)[i_reco]>0 || (*reco_E_HB)[i_reco]>0){
                has_E_Barrel=true;
                E_barrel+=(*reco_E_EB)[i_reco]+(*reco_E_HB)[i_reco];
            }
            if((*reco_E_EE)[i_reco]>0 || (*reco_E_HE)[i_reco]>0){
                has_E_Endcap=true;
                E_endcap+=(*reco_E_EE)[i_reco]+(*reco_E_HE)[i_reco];
            }
        }
        //end of reco loop check candidate, if existing for energy domination about all the rest
        if(E_had_candidate_max>0 && (E_had_candidate_max/E_reco_total)>0.99){
            find_reco_candidate=true;
            count_hadron_candidate_events+=1;
            //if((E_had_candidate_max/E_reco_total)<1){
            //}
            std::cout<<"candidate blabla "<<E_had_candidate_max/E_reco_total<<"/"<<fabs(true_CosTheta)<<std::endl;
            if(fabs(true_CosTheta)>=0.83 && fabs(true_CosTheta)<0.94){
                h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->Fill(E_had_candidate_max/true_energy_hadron);
            }else if (fabs(true_CosTheta)>=0.78 && fabs(true_CosTheta)<0.83){
                h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->Fill(E_had_candidate_max/true_energy_hadron);
            }else{
                std::cout<<"fill something in there "<<E_had_candidate_max/true_energy_hadron<<std::endl;
                if (fabs(true_CosTheta)<0.60){
                    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Fill(E_had_candidate_max/true_energy_hadron);
                }
                if (fabs(true_CosTheta)<0.78){
                    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->Fill(E_had_candidate_max/true_energy_hadron);
                }
            }
        }
        if(has_E_Endcap && has_E_Barrel){
            if(E_endcap>(0.01*(E_endcap+E_barrel)) && E_barrel>(0.01*(E_endcap+E_barrel))){
                h_SWC_CLIC_n_K0L_cosThetaTR->Fill(true_CosTheta);
                //std::cout<<"E_E/E_B/costheta/min/max "<<E_barrel<<"/"<<E_endcap<<"/"<<true_CosTheta<<"/"<<true_CosTheta_min<<"/"<<true_CosTheta_max<<std::endl;
                if(fabs(true_CosTheta)<true_CosTheta_min){
                    true_CosTheta_min=fabs(true_CosTheta);
                }
                if(fabs(true_CosTheta)>true_CosTheta_max){
                    true_CosTheta_max=fabs(true_CosTheta);
                }
            }
        }
    }
    
    //value 1,5,10,15,30,50,100,200,500,1000,1500
    TCanvas* can_res_rel_had_0_60=setUpperCanvas("can_h_res_rel_had_0_60");
    can_res_rel_had_0_60->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Draw("P,e,same");
    TF1* fit_RECO_0_60 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->GetFunction("gaus");
    
    std::cout<<"gre->SetPoint(10,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_60->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre->SetPointError(10,0,"<<100.*fit_RECO_0_60->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_78=setUpperCanvas("can_h_res_rel_had_0_78");
    can_res_rel_had_0_78->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->Draw("P,e,same");
    TF1* fit_RECO_0_78 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78->GetFunction("gaus");
    
    std::cout<<"gre->SetPoint(10,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_78->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre->SetPointError(10,0,"<<100.*fit_RECO_0_78->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_78_to_0_83=setUpperCanvas("can_h_res_rel_had_0_78_to_0_83");
    can_res_rel_had_0_78_to_0_83->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->Draw("P,e,same");
    TF1* fit_RECO_0_78_to_0_83 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_78_to_0_83->GetFunction("gaus");
    
    std::cout<<"gre->SetPoint(10,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_78_to_0_83->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre->SetPointError(10,0,"<<100.*fit_RECO_0_78_to_0_83->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_83_to_0_94=setUpperCanvas("can_h_res_rel_had_0_83_to_0_94");
    can_res_rel_had_0_83_to_0_94->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->Draw("P,e,same");
    TF1* fit_RECO_0_83_to_0_94 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_83_to_0_94->GetFunction("gaus");
    
    std::cout<<"gre->SetPoint(10,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_83_to_0_94->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre->SetPointError(10,0,"<<100.*fit_RECO_0_83_to_0_94->GetParError(2)<<");"<<std::endl;
    
    
    std::cout<<"gre_mean->SetPoint(10,"<<true_energy_fillplot<<","<<fit_RECO_0_60->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPointError(10,0,"<<fit_RECO_0_60->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPoint(10,"<<true_energy_fillplot<<","<<fit_RECO_0_78->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPointError(10,0,"<<fit_RECO_0_78->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPoint(10,"<<true_energy_fillplot<<","<<fit_RECO_0_78_to_0_83->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPointError(10,0,"<<fit_RECO_0_78_to_0_83->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPoint(10,"<<true_energy_fillplot<<","<<fit_RECO_0_83_to_0_94->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_mean->SetPointError(10,0,"<<fit_RECO_0_83_to_0_94->GetParError(1)<<");"<<std::endl;
    
    
    TCanvas* can_cosThetaTR=setUpperCanvas("can_cosThetaTR");
    can_cosThetaTR->cd();
    h_SWC_CLIC_n_K0L_cosThetaTR->Draw();
    
    std::cout<<"file end of this costheta_min/max "<<true_CosTheta_min<<"/"<<true_CosTheta_max<<" cand found/tot "<< count_hadron_candidate_events <<"/"<< tree_SWC_CLIC_n_K0L->GetEntries() <<std::endl;
    
    */
    
    //try barrel up to 0.6, endcap 0.82 to 0.92 TR 0.65 to 0.80
    

    TFile* file_histogram=new TFile(final_histo_name,"recreate");

    //now photon 0.8 to 0.82
    
    //=========Macro generated from canvas: resolutionGraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph = new TCanvas("resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary", "resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->Range(-186.894,-0.873515,1682.046,6.114605);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetFillColor(0);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetBorderMode(0);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetBorderSize(2);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetGridx();
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetGridy();
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetRightMargin(0.0172);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetTopMargin(0.055);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetBottomMargin(0.138);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetFrameBorderMode(0);
    resolutionGraphCanvas_CLIC_n_K0L_res_rel_ph->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre = new TGraphErrors(11);
    gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_60");
    gre->SetTitle("");
    gre->SetFillColor(1);
    gre->SetMarkerColor(kBlack);
    gre->SetMarkerStyle(kOpenCircle);
    gre->SetPoint(0,1,14.7341);
    gre->SetPointError(0,0,0.166047);
    gre->SetPoint(1,5,6.89549);
    gre->SetPointError(1,0,0.0710082);
    gre->SetPoint(2,10,4.97863);
    gre->SetPointError(2,0,0.0475384);
    gre->SetPoint(3,15,4.23218);
    gre->SetPointError(3,0,0.0397799);
    gre->SetPoint(4,30,2.97557);
    gre->SetPointError(4,0,0.0290383);
    gre->SetPoint(5,50,2.39);
    gre->SetPointError(5,0,0.0254124);
    gre->SetPoint(6,100,1.78017);
    gre->SetPointError(6,0,0.0180125);
    gre->SetPoint(7,200,1.34563);
    gre->SetPointError(7,0,0.014343);
    gre->SetPoint(8,500,1.07258);
    gre->SetPointError(8,0,0.0124886);
    gre->SetPoint(9,1000,0.940544);
    gre->SetPointError(9,0,0.0122502);
    gre->SetPoint(10,1500,0.906059);
    gre->SetPointError(10,0,0.0109578);
    
    TH1F *h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(0.5);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(17.5);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre->SetHistogram(h_ph_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    gre->Draw("ape");
    /*
    TF1* fit_sigma_ph_0_60 = new TF1("fit_sigma_ph_0_60", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",4,75);
    fit_sigma_ph_0_60->SetParameter(0,50.);
    fit_sigma_ph_0_60->SetParameter(1,0.10);
    fit_sigma_ph_0_60->SetParameter(2,3.0);
    fit_sigma_ph_0_60->SetParNames("a1_0_60","a2_0_60","const_0_60");
    gre->Fit("fit_sigma_ph_0_60","R");
    
    TF1* fit_sigma_ph_0_60_only_E = new TF1("fit_sigma_ph_0_60_only_E", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",4,75.);
    fit_sigma_ph_0_60_only_E->SetParameter(0,50.0);
    fit_sigma_ph_0_60_only_E->SetParameter(1,3.0);
    fit_sigma_ph_0_60_only_E->SetParNames("a1_0_60_only_E","const_0_60_only_E");
    gre->Fit("fit_sigma_ph_0_60_only_E","R");
    */
    gre->Draw("pesame");
    
    gre = new TGraphErrors(11);
    gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78");
    gre->SetTitle("");
    gre->SetFillColor(1);
    gre->SetMarkerColor(kBlack);
    gre->SetMarkerStyle(kOpenCircle);
    gre->SetPoint(0,1,14.9688);
    gre->SetPointError(0,0,0.151503);
    gre->SetPoint(1,5,7.00205);
    gre->SetPointError(1,0,0.063763);
    gre->SetPoint(2,10,5.02102);
    gre->SetPointError(2,0,0.0428576);
    gre->SetPoint(3,15,4.25148);
    gre->SetPointError(3,0,0.0352131);
    gre->SetPoint(4,30,2.99076);
    gre->SetPointError(4,0,0.0252701);
    gre->SetPoint(5,50,2.42809);
    gre->SetPointError(5,0,0.0221697);
    gre->SetPoint(6,100,1.78763);
    gre->SetPointError(6,0,0.01637);
    gre->SetPoint(7,200,1.33247);
    gre->SetPointError(7,0,0.0125256);
    gre->SetPoint(8,500,1.03128);
    gre->SetPointError(8,0,0.0105665);
    gre->SetPoint(9,1000,0.894994);
    gre->SetPointError(9,0,0.010169);
    gre->SetPoint(10,1500,0.860321);
    gre->SetPointError(10,0,0.00923972);
    
    TH1F *h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14","",100,0,415.);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetMinimum(0.5);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetMaximum(17.5);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetDirectory(0);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetStats(0);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetLineColor(kBlack);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre->SetHistogram(h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14);
    gre->Draw("ape");
    /* 
       TF1* fit_sigma_ph_0_78 = new TF1("fit_sigma_ph_0_78", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",4,75);
    fit_sigma_ph_0_78->SetParameter(0,50.);
    fit_sigma_ph_0_78->SetParameter(1,0.10);
    fit_sigma_ph_0_78->SetParameter(2,3.0);
    fit_sigma_ph_0_78->SetParNames("a1_0_78","a2_0_78","const_0_78");
    gre->Fit("fit_sigma_ph_0_78","R");
    
    TF1* fit_sigma_ph_0_78_only_E = new TF1("fit_sigma_ph_0_78_only_E", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",4,75);
    fit_sigma_ph_0_78_only_E->SetParameter(0,50.0);
    fit_sigma_ph_0_78_only_E->SetParameter(1,3.0);
    fit_sigma_ph_0_78_only_E->SetParNames("a1_0_78_only_E","const_0_78_only_E");
    gre->Fit("fit_sigma_ph_0_78_only_E","R");
    */
    gre->Draw("pesame");
    
    gre = new TGraphErrors(11);
    gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83");
    gre->SetTitle("");
    gre->SetFillColor(1);
    gre->SetMarkerColor(kRed-7);
    gre->SetMarkerStyle(kOpenSquare);
    gre->SetPoint(0,1,16.5011);
    gre->SetPointError(0,0,0.865444);
    gre->SetPoint(1,5,7.50338);
    gre->SetPointError(1,0,0.421356);
    gre->SetPoint(2,10,5.92026);
    gre->SetPointError(2,0,0.277318);
    gre->SetPoint(3,15,4.55544);
    gre->SetPointError(3,0,0.192466);
    gre->SetPoint(4,30,3.10068);
    gre->SetPointError(4,0,0.127407);
    gre->SetPoint(5,50,2.49905);
    gre->SetPointError(5,0,0.115118);
    gre->SetPoint(6,100,1.81596);
    gre->SetPointError(6,0,0.0904222);
    gre->SetPoint(7,200,1.27697);
    gre->SetPointError(7,0,0.0588278);
    gre->SetPoint(8,500,1.19124);
    gre->SetPointError(8,0,0.109471);
    gre->SetPoint(9,1000,0.797591);
    gre->SetPointError(9,0,0.0538643);
    gre->SetPoint(10,1500,0.761724);
    gre->SetPointError(10,0,0.0436831);
    
    TH1F *h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","",100,0,415);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetMinimum(0.5);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetMaximum(17.5);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetDirectory(0);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetStats(0);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetLineColor(kRed-7);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre->SetHistogram(h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14);
    gre->Draw("pe");
    /*
    TF1* fit_sigma_ph_0_78_to_0_83 = new TF1("fit_sigma_ph_0_78_to_0_83", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",4,75);
    fit_sigma_ph_0_78_to_0_83->SetParameter(0,50.);
    fit_sigma_ph_0_78_to_0_83->SetParameter(1,0.10);
    fit_sigma_ph_0_78_to_0_83->SetParameter(2,3.0);
    fit_sigma_ph_0_78_to_0_83->SetParNames("a1_0_78_to_0_83","a2_0_78_to_0_83","const_0_78_to_0_83");
    gre->Fit("fit_sigma_ph_0_78_to_0_83","R");
    
    TF1* fit_sigma_ph_0_78_to_0_83_only_E = new TF1("fit_sigma_ph_0_78_to_0_83_only_E", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",4,75);
    fit_sigma_ph_0_78_to_0_83_only_E->SetParameter(0,50.0);
    fit_sigma_ph_0_78_to_0_83_only_E->SetParameter(1,3.0);
    fit_sigma_ph_0_78_to_0_83_only_E->SetParNames("a1_0_78_to_0_83_only_E","const_0_78_to_0_83_only_E");
    gre->Fit("fit_sigma_ph_0_78_to_0_83_only_E","R");
    */
    gre->Draw("pesame");
    
    gre = new TGraphErrors(11);
    gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94");
    gre->SetTitle("");
    gre->SetFillColor(1);
    gre->SetMarkerColor(kBlue);
    gre->SetMarkerStyle(kOpenTriangleUp);
    gre->SetPoint(0,1,14.9116);
    gre->SetPointError(0,0,0.429003);
    gre->SetPoint(1,5,6.97093);
    gre->SetPointError(1,0,0.189273);
    gre->SetPoint(2,10,4.74473);
    gre->SetPointError(2,0,0.131706);
    gre->SetPoint(3,15,3.88792);
    gre->SetPointError(3,0,0.0957513);
    gre->SetPoint(4,30,2.78804);
    gre->SetPointError(4,0,0.0697395);
    gre->SetPoint(5,50,2.2429);
    gre->SetPointError(5,0,0.0687657);
    gre->SetPoint(6,100,1.59499);
    gre->SetPointError(6,0,0.040368);
    gre->SetPoint(7,200,1.15164);
    gre->SetPointError(7,0,0.0276968);
    gre->SetPoint(8,500,0.735529);
    gre->SetPointError(8,0,0.0180497);
    gre->SetPoint(9,1000,0.551053);
    gre->SetPointError(9,0,0.016265);
    gre->SetPoint(10,1500,0.451682);
    gre->SetPointError(10,0,0.0127654);
    
    TH1F *h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14 = new TH1F("h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14","",100,0,215.);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetMinimum(0.5);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetMaximum(17.5);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetStats(0);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre->SetHistogram(h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14);
    gre->Draw("pe");
    /*
    TF1* fit_sigma_ph_0_83_to_0_92 = new TF1("fit_sigma_ph_0_83_to_0_92", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",4,75);
    fit_sigma_ph_0_83_to_0_92->SetParameter(0,50.0);
    fit_sigma_ph_0_83_to_0_92->SetParameter(1,0.10);
    fit_sigma_ph_0_83_to_0_92->SetParameter(2,3.0);
    fit_sigma_ph_0_83_to_0_92->SetParNames("a1_0_83_to_0_92","a2_0_83_to_0_92","const_0_83_to_0_92");
    gre->Fit("fit_sigma_ph_0_83_to_0_92","R");
    
    TF1* fit_sigma_ph_0_83_to_0_92_only_E = new TF1("fit_sigma_ph_0_83_to_0_92_only_E", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",4,75);
    fit_sigma_ph_0_83_to_0_92_only_E->SetParameter(0,50.0);
    fit_sigma_ph_0_83_to_0_92_only_E->SetParameter(1,3.0);
    fit_sigma_ph_0_83_to_0_92_only_E->SetParNames("a1_0_83_to_0_92_only_E","const_0_83_to_0_92_only_E");
    gre->Fit("fit_sigma_ph_0_83_to_0_92_only_E","R");
    gre->Draw("pesame");
    */
    
    
    
    TLegend *leg_ph_E_reco_vs_E_true = new TLegend(0.60,0.60,0.94,0.91,"ph","brNDC");
    leg_ph_E_reco_vs_E_true->SetBorderSize(1);
    leg_ph_E_reco_vs_E_true->SetTextSize(0.027);
    leg_ph_E_reco_vs_E_true->SetLineColor(1);
    leg_ph_E_reco_vs_E_true->SetLineStyle(1);
    leg_ph_E_reco_vs_E_true->SetLineWidth(1);
    leg_ph_E_reco_vs_E_true->SetFillColor(0);
    leg_ph_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_ph_E_reco_vs_E_true=leg_ph_E_reco_vs_E_true->AddEntry("h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v14","|cos#theta_{#gamma}^{true}|<0.78","pe");
    leg_entry_ph_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_ph_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_ph_E_reco_vs_E_true=leg_ph_E_reco_vs_E_true->AddEntry("h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","0.78<|cos#theta_{#gamma}^{true}|<0.83","pe");
    leg_entry_ph_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_ph_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_ph_E_reco_vs_E_true=leg_ph_E_reco_vs_E_true->AddEntry("h_ph_E_reco_vs_E_true_0_83_to_094_CLIC_o3_v14","0.80<|cos#theta_{#gamma}^{true}|<0.94","pe");
    leg_entry_ph_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_ph_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_ph_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_E_reco_vs_E_true->SetTextFont(42);
    leg_ph_E_reco_vs_E_true->Draw();
    
    //now for mean
    
    //=========Macro generated from canvas: resolutionGraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph = new TCanvas("resolution_mean_GraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary", "resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetFillColor(0);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetBorderMode(0);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetBorderSize(2);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetGridx();
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetGridy();
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetRightMargin(0.0172);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetTopMargin(0.055);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetBottomMargin(0.138);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetFrameBorderMode(0);
    resolution_mean_GraphCanvas_CLIC_n_K0L_res_rel_ph->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre_mean = new TGraphErrors(11);
    gre_mean->SetName("CLIC_o3_v14_res_graph_ph_mean_rel_E_reco_vs_E_true_0_60");
    gre_mean->SetTitle("");
    gre_mean->SetFillColor(1);
    gre_mean->SetMarkerColor(kBlack);
    gre_mean->SetMarkerStyle(kOpenCircle);
    gre_mean->SetPoint(0,1,0.980207);
    gre_mean->SetPointError(0,0,0.00214853);
    gre_mean->SetPoint(1,5,1.00043);
    gre_mean->SetPointError(1,0,0.000933415);
    gre_mean->SetPoint(2,10,1.00595);
    gre_mean->SetPointError(2,0,0.000677113);
    gre_mean->SetPoint(3,15,1.00939);
    gre_mean->SetPointError(3,0,0.000580124);
    gre_mean->SetPoint(4,30,1.01406);
    gre_mean->SetPointError(4,0,0.000409721);
    gre_mean->SetPoint(5,50,1.01662);
    gre_mean->SetPointError(5,0,0.000327186);
    gre_mean->SetPoint(6,100,1.02189);
    gre_mean->SetPointError(6,0,0.000243184);
    gre_mean->SetPoint(7,200,1.02627);
    gre_mean->SetPointError(7,0,0.000185809);
    gre_mean->SetPoint(8,500,1.03181);
    gre_mean->SetPointError(8,0,0.000161659);
    gre_mean->SetPoint(9,1000,1.03607);
    gre_mean->SetPointError(9,0,0.000166876);
    gre_mean->SetPoint(10,1500,1.03865);
    gre_mean->SetPointError(10,0,0.000154155);
    
    TH1F *h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(0.9);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(1.1);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{#gamma}");
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_mean->SetHistogram(h_ph_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    gre_mean->Draw("ape");
    
    gre_mean = new TGraphErrors(11);
    gre_mean->SetName("CLIC_o3_v14_res_graph_ph_mean_rel_E_reco_vs_E_true_0_78");
    gre_mean->SetTitle("");
    gre_mean->SetFillColor(1);
    gre_mean->SetMarkerColor(kBlack);
    gre_mean->SetMarkerStyle(kOpenCircle);
    gre_mean->SetPoint(0,1,0.981614);
    gre_mean->SetPointError(0,0,0.001918);
    gre_mean->SetPoint(1,5,1.00154);
    gre_mean->SetPointError(1,0,0.000845255);
    gre_mean->SetPoint(2,10,1.0068);
    gre_mean->SetPointError(2,0,0.000602767);
    gre_mean->SetPoint(3,15,1.0105);
    gre_mean->SetPointError(3,0,0.000512398);
    gre_mean->SetPoint(4,30,1.01483);
    gre_mean->SetPointError(4,0,0.000363289);
    gre_mean->SetPoint(5,50,1.01778);
    gre_mean->SetPointError(5,0,0.000294353);
    gre_mean->SetPoint(6,100,1.02237);
    gre_mean->SetPointError(6,0,0.000214437);
    gre_mean->SetPoint(7,200,1.02651);
    gre_mean->SetPointError(7,0,0.000161854);
    gre_mean->SetPoint(8,500,1.03173);
    gre_mean->SetPointError(8,0,0.00013758);
    gre_mean->SetPoint(9,1000,1.03552);
    gre_mean->SetPointError(9,0,0.000137781);
    gre_mean->SetPoint(10,1500,1.03768);
    gre_mean->SetPointError(10,0,0.00012845);
    
    TH1F *h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14","",100,0,415.);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetMinimum(0.9);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetMaximum(1.1);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetDirectory(0);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetStats(0);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->SetLineColor(kBlack);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{#gamma}");
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_mean->SetHistogram(h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14);
    gre_mean->Draw("ape");
    
    gre_mean = new TGraphErrors(11);
    gre_mean->SetName("CLIC_o3_v14_res_graph_ph_mean_rel_E_reco_vs_E_true_0_78_to_0_83");
    gre_mean->SetTitle("");
    gre_mean->SetFillColor(1);
    gre_mean->SetMarkerColor(kBlue);
    gre_mean->SetMarkerStyle(kOpenSquare);
    gre_mean->SetPoint(0,1,0.983034);
    gre_mean->SetPointError(0,0,0.0096614);
    gre_mean->SetPoint(1,5,0.997193);
    gre_mean->SetPointError(1,0,0.00412138);
    gre_mean->SetPoint(2,10,1.00949);
    gre_mean->SetPointError(2,0,0.00325247);
    gre_mean->SetPoint(3,15,1.01304);
    gre_mean->SetPointError(3,0,0.00236323);
    gre_mean->SetPoint(4,30,1.01391);
    gre_mean->SetPointError(4,0,0.00169211);
    gre_mean->SetPoint(5,50,1.01995);
    gre_mean->SetPointError(5,0,0.00143809);
    gre_mean->SetPoint(6,100,1.02337);
    gre_mean->SetPointError(6,0,0.00107785);
    gre_mean->SetPoint(7,200,1.0267);
    gre_mean->SetPointError(7,0,0.00080888);
    gre_mean->SetPoint(8,500,1.03039);
    gre_mean->SetPointError(8,0,0.000800417);
    gre_mean->SetPoint(9,1000,1.03514);
    gre_mean->SetPointError(9,0,0.000580246);
    gre_mean->SetPoint(10,1500,1.0356);
    gre_mean->SetPointError(10,0,0.000515327);
    
    TH1F *h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","",100,0,415.);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetMinimum(0.9);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetMaximum(1.1);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetDirectory(0);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetStats(0);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->SetLineColor(kBlue);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{#gamma}");
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_mean->SetHistogram(h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14);
    gre_mean->Draw("pe");
    
    gre_mean = new TGraphErrors(11);
    gre_mean->SetName("CLIC_o3_v14_res_graph_ph_mean_rel_E_reco_vs_E_true_0_83_to_0_94");
    gre_mean->SetTitle("");
    gre_mean->SetFillColor(1);
    gre_mean->SetMarkerColor(kRed-7);
    gre_mean->SetMarkerStyle(kOpenTriangleUp);
    gre_mean->SetPoint(0,1,0.995968);
    gre_mean->SetPointError(0,0,0.00539088);
    gre_mean->SetPoint(1,5,1.00741);
    gre_mean->SetPointError(1,0,0.00246618);
    gre_mean->SetPoint(2,10,1.01483);
    gre_mean->SetPointError(2,0,0.00169471);
    gre_mean->SetPoint(3,15,1.01406);
    gre_mean->SetPointError(3,0,0.00134921);
    gre_mean->SetPoint(4,30,1.01696);
    gre_mean->SetPointError(4,0,0.000956952);
    gre_mean->SetPoint(5,50,1.01775);
    gre_mean->SetPointError(5,0,0.000782208);
    gre_mean->SetPoint(6,100,1.02219);
    gre_mean->SetPointError(6,0,0.000545847);
    gre_mean->SetPoint(7,200,1.02461);
    gre_mean->SetPointError(7,0,0.000388901);
    gre_mean->SetPoint(8,500,1.02826);
    gre_mean->SetPointError(8,0,0.000271027);
    gre_mean->SetPoint(9,1000,1.03057);
    gre_mean->SetPointError(9,0,0.000231621);
    gre_mean->SetPoint(10,1500,1.03145);
    gre_mean->SetPointError(10,0,0.000183809);
    
    TH1F *h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14 = new TH1F("h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14","",100,0,415.);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetMinimum(0.9);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetMaximum(1.1);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetStats(0);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{#gamma}");
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_mean->SetHistogram(h_ph_mean_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v14);
    gre_mean->Draw("pe");
    
    TLegend *leg_ph_mean_E_reco_vs_E_true = new TLegend(0.60,0.20,0.94,0.60,"K^{0}_{L}","brNDC");
    leg_ph_mean_E_reco_vs_E_true->SetBorderSize(1);
    leg_ph_mean_E_reco_vs_E_true->SetTextSize(0.027);
    leg_ph_mean_E_reco_vs_E_true->SetLineColor(1);
    leg_ph_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_ph_mean_E_reco_vs_E_true->SetLineWidth(1);
    leg_ph_mean_E_reco_vs_E_true->SetFillColor(0);
    leg_ph_mean_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_ph_mean_E_reco_vs_E_true=leg_ph_mean_E_reco_vs_E_true->AddEntry("h_ph_mean_E_reco_vs_E_true_0_78_CLIC_o3_v14","|cos#theta_{#gamma}^{true}|<0.78","pe");
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_ph_mean_E_reco_vs_E_true=leg_ph_mean_E_reco_vs_E_true->AddEntry("h_ph_mean_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v14","0.78<|cos#theta_{#gamma}^{true}|<0.83","pe");
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_ph_mean_E_reco_vs_E_true=leg_ph_mean_E_reco_vs_E_true->AddEntry("h_ph_mean_E_reco_vs_E_true_0_83_to_094_CLIC_o3_v14","0.83<|cos#theta_{#gamma}^{true}|<0.94","pe");
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_ph_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_ph_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_ph_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_ph_mean_E_reco_vs_E_true->Draw();
  
    

  std::cout<<"signal/no signal events "<<signal_events<<"/"<<no_signal_events<<std::endl;
  
  //tEff_SignalPion_CosTheta->Write();
  //h_energyTot_pionInEvent->Write();
  //h_energy_pions->Write();

  file_histogram->Write();
  file_histogram->Close();

    //500 new tune
    //  mean 500.076, sigma 4.84/rel 0.968 %
    // mean 9.999, sigma 0.0310, real 0.3125
    //mean 0.9986, sigma 0.003887, 0.389
    //mean 19.999, sigma 0.0600, 0.300
    // mean 49.99, sigma 0.165, 0.329
    //mean 75.001, sigma 0.263, 0.3510
    //mean 99.998, sigma 0.3714, 0.3715
    
  return 1;

}

