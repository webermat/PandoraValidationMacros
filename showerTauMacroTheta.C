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
#include "TLegend.h"
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

 
void setTDRStyle() {
    TStyle *tdrStyle = new TStyle("tdrStyle","StyleA for P-TDR");
    
    // For the canvas:
    tdrStyle->SetCanvasBorderMode(0); //tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasColor(0);
    tdrStyle->SetCanvasDefH(450); //Height of canvas
    tdrStyle->SetCanvasDefW(590); //Width of canvas
   //10,50,600,350
    tdrStyle->SetCanvasDefX(10);   //POsition on screen
    tdrStyle->SetCanvasDefY(50);
    tdrStyle->SetCanvasBorderMode(0);
    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
    // tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadBorderSize(0); 
    tdrStyle->SetPadColor(1);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);
    
    // For the frame: tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);
    
    // For the h_o:
    // tdrStyle->SetH_FillColor(1);
    // tdrStyle->SetH_FillStyle(0);
    tdrStyle->SetHistLineColor(1);
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(3);
    // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
    // tdrStyle->SetNumberContours(Int_t number = 20);
    
    tdrStyle->SetEndErrorSize(2);
    //tdrStyle->SetErrorMarker(20);
    tdrStyle->SetErrorX(0);
    
    //tdrStyle->SetMarkerStyle(20);
    
    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);
    
    //For the date:
    tdrStyle->SetOptDate(0);
    // tdrStyle->SetDateX(Float_t x = 0.01);
    // tdrStyle->SetDateY(Float_t y = 0.01);
    
    // For the statistics box:
    tdrStyle->SetOptFile(0);
    //tdrStyle->SetOptStat(0);
    //tdrStyle->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.025);
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat("6.4g");
    tdrStyle->SetStatBorderSize(1);
    tdrStyle->SetStatH(0.1);
    tdrStyle->SetStatW(0.15);
    tdrStyle->SetPaintTextFormat("5.2f");
    // tdrStyle->SetStatStyle(Style_t style = 1001);
    // tdrStyle->SetStatX(Float_t x = 0);
    // tdrStyle->SetStatY(Float_t y = 0);
    
    // Margins:
    tdrStyle->SetPadTopMargin(0.05);
    tdrStyle->SetPadBottomMargin(0.14);
    tdrStyle->SetPadLeftMargin(0.18);
    tdrStyle->SetPadRightMargin(0.02);
    
    // For the Global title:
    
    tdrStyle->SetOptTitle(1);
    tdrStyle->SetTitleFont(42);
    tdrStyle->SetTitleColor(1);
    tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.05);
    // tdrStyle->SetTitleH(0); // Set the height of the title box
    //tdrStyle->SetTitleW(0.8f); // Set the width of the title box
    //tdrStyle->SetTitleX(1.3); // Set the position of the title box
    //tdrStyle->SetTitleY(0.985); // Set the position of the title box
    // tdrStyle->SetTitleStyle(Style_t style = 1001);
    tdrStyle->SetTitleBorderSize(0);
    
    // For the axis titles:
    
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(42, "XYZ");
    tdrStyle->SetTitleBorderSize(0);
    tdrStyle->SetTitleSize(0.06, "XYZ");
    //tdrStyle->SetTitleXSize(0.055); // Another way to set the size?
    //tdrStyle->SetTitleYSize(0.06);
    //tdrStyle->SetTitleXOffset(1.95);
    //tdrStyle->SetTitleYOffset(1.35);
    //tdrStyle->SetTitleAlign(33);
    tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    
    tdrStyle->SetLabelColor(1, "XYZ");
    //font with ID4 (helvetica and precision 
    tdrStyle->SetLabelFont(42, "XYZ");
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.05, "XZ");
    tdrStyle->SetLabelSize(0.0425, "Y");
    
    
    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    //tdrStyle->SetTickLength(0.06, "Y");
    //tdrStyle->SetTickLength(0.05, "XZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);
    
    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);
    
    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
    tdrStyle->SetPalette(1,0);
    // tdrStyle->SetLineScalePS(Float_t scale = 3);
    // tdrStyle->SetLineStyleString(Int_t i, const char* text);
    // tdrStyle->SetHeaderPS(const char* header);
    // tdrStyle->SetTitlePS(const char* pstitle);
    
    // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
    // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
    // tdrStyle->SetPaintTextFormat(const char* format = "g");
    // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // tdrStyle->SetTimeOffset(Float_t toffset);
    // tdrStyle->SetH_MinimumZero(kTRUE);
    
    tdrStyle->cd();
    //SC5
    
    std::cout<<"tdrstyle set"<<std::endl;
    
    //return 1;
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
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}

TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,250);
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
//}


int plotTauShower(){

  gROOT->ProcessLine("#include <vector>");

    string label_legend= "#tau, E_{true}=200 GeV";
  //TFile* file_photons=TFile::Open("/home/weberma2/data/tauMacroFiles/tauStudy_whz_tautau_1000_ILC170823_TT_FitFW_CLIC_o3_v13_SWC_9Ebin_500K0LCalib.root");
  TFile* file_photons=TFile::Open("/eos/user/w/weberma2/data/tauMacroFiles/181011_gcc62/tauStudy_whz_tautau_2000_CLIC_o3_v14_CT.root");


  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/RootOutputFiles/tauFiles/pionStudy_pyt6_Zuds200_PhLhd_ClicPerf_Steer_TruthTrack_FitFW_ILC170202_CLIC_o3_v08_170217_status170404_correctTruthStatus1.root"); 


  setTDRStyle();

  unsigned int charged_pion_ID=211;

  float lim_energy_true_low = 10.0;
  float lim_energy_true_high = 1000;

  int test_taup = -15;
  int test_taum = 15;
  const char*  header_label_taup="tau^{+}, E=500 GeV";
  const char*  header_label_taum="tau^{-}, E=500 GeV";
  const char*  signal_label_taum="tau^{-} identified";
  const char*  signal_label_taup="tau^{-} identified";

  TTree* tree_photons = (TTree*)file_photons->Get("showerData");

  int eventNumber=0;

  vector<float> *tauJetE=0;
  vector<float> *tauJetPx=0;
  vector<float> *tauJetPy=0;
  vector<float> *tauJetPz=0;
  //vector<float> *tauJetPhi=0;
  vector<float> *tauJetCosTheta=0;
  vector<int> *tauJetCharge=0;

  vector<float> *tauJetPartE=0;
  vector<float> *tauJetPartPx=0;
  vector<float> *tauJetPartPy=0;
  vector<float> *tauJetPartPz=0;
  vector<int> *tauJetPartPDGID=0;
  vector<int> *tauJetPartJetIndex=0;
  vector<int> *tauJetPartCharge=0;


  vector<float> *trueE=0;
  vector<float> *truePx=0;
  vector<float> *truePy=0;
  vector<float> *truePz=0;
  vector<float> *truePhi=0;
  vector<float> *trueCosTheta=0;
  //vector<float> *trueTheta=0;
  vector<int> *truePDGID=0;
  vector<int> *trueIndex=0;
  vector<int> *trueGenStatus=0;
  vector<int> *trueNumDaughters=0;
  vector<int> *trueDecayTrackerCalo=0;//1 if decayed in tracker, 2 if decayed in tracker, 3 if coming from backscatter, 4 if left the detector
  vector<int> *trueMotherDecayTrackerCalo=0;//1 if decayed in tracker, 2 if decayed in calo

  vector<float> *trueTauDaughterE=0;
  vector<float> *trueTauDaughterPx=0;
  vector<float> *trueTauDaughterPy=0;
  vector<float> *trueTauDaughterPz=0;
  vector<int> *trueTauDaughterPDGID=0;
  vector<int> *trueTauDaughterCharge=0;
  vector<int> *trueTauDaughterTauIndex=0;
  vector<int> *trueTauDaughterMotherPDGID=0;
  vector<float> *trueTauDaughterMotherEnergy=0;

  vector<float> *recoE=0;
  vector<float> *recoPx=0;
  vector<float> *recoPy=0;
  vector<float> *recoPz=0;
  vector<float> *recoPhi=0;
  vector<float> *recoCosTheta=0;
  vector<float> *recoTheta=0;
  vector<float> *recoPhi_logERW=0;
  vector<float> *recoCosTheta_logERW=0;
  vector<float> *recoTheta_logERW=0;
  vector<int> *recoCharge=0;
  vector<int> *recoPDGID=0;
  vector<float> *recotrack0p=0;
  vector<float> *recotrack0pt=0;
  vector<float> *recotrack0Chi2OverNdof=0;
  vector<float> *recotrack0nHits=0;
  vector<float> *recoEEB=0;
  vector<float> *recoEEE=0;
  vector<float> *recoEEP=0;
  vector<float> *recoEHB=0;
  vector<float> *recoEHE=0;
  vector<float> *recoEHP=0;

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
  //tree_photons->SetBranchAddress("true_Theta",&trueTheta);
  tree_photons->SetBranchAddress("true_CosTheta",&trueCosTheta);
  tree_photons->SetBranchAddress("true_PDGID",&truePDGID);
  tree_photons->SetBranchAddress("true_GenStatus",&trueGenStatus);
  tree_photons->SetBranchAddress("true_decayTrackerCalo",&trueDecayTrackerCalo);
  tree_photons->SetBranchAddress("true_motherDecayTrackerCalo",&trueMotherDecayTrackerCalo);
  tree_photons->SetBranchAddress("true_CosTheta",&trueCosTheta);
  //tree_photons->SetBranchAddress("true_Theta",&trueTheta);
  tree_photons->SetBranchAddress("true_PDGID",&truePDGID);
  tree_photons->SetBranchAddress("true_index",&trueIndex);
  tree_photons->SetBranchAddress("true_numDaughters",&trueNumDaughters);
  tree_photons->SetBranchAddress("trueTauDaughterE",&trueTauDaughterE);
  tree_photons->SetBranchAddress("trueTauDaughterPx",&trueTauDaughterPx);
  tree_photons->SetBranchAddress("trueTauDaughterPy",&trueTauDaughterPy);
  tree_photons->SetBranchAddress("trueTauDaughterPz",&trueTauDaughterPz);
  tree_photons->SetBranchAddress("trueTauDaughterPDGID",&trueTauDaughterPDGID);
  tree_photons->SetBranchAddress("trueTauDaughterCharge",&trueTauDaughterCharge);
  tree_photons->SetBranchAddress("trueTauDaughterTauIndex",&trueTauDaughterTauIndex);
  tree_photons->SetBranchAddress("trueTauDaughterMotherPDGID",&trueTauDaughterMotherPDGID);
  tree_photons->SetBranchAddress("trueTauDaughterMotherEnergy",&trueTauDaughterMotherEnergy);

  tree_photons->SetBranchAddress("tauJetPx",&tauJetPx);
  tree_photons->SetBranchAddress("tauJetPy",&tauJetPy);
  tree_photons->SetBranchAddress("tauJetPz",&tauJetPz);
  tree_photons->SetBranchAddress("tauJetE",&tauJetE);
  //tree_photons->SetBranchAddress("tauJetPhi",&tauJetPhi);
  tree_photons->SetBranchAddress("tauJetcharge",&tauJetCharge);
  tree_photons->SetBranchAddress("tauJetCosTheta",&tauJetCosTheta);

  tree_photons->SetBranchAddress("tauJet_Part_Px",&tauJetPartPx);
  tree_photons->SetBranchAddress("tauJet_Part_Py",&tauJetPartPy);
  tree_photons->SetBranchAddress("tauJet_Part_Pz",&tauJetPartPz);
  tree_photons->SetBranchAddress("tauJet_Part_E",&tauJetPartE);
  tree_photons->SetBranchAddress("tauJet_Part_PDGID",&tauJetPartPDGID);
  tree_photons->SetBranchAddress("tauJet_Part_charge",&tauJetPartCharge);
  tree_photons->SetBranchAddress("tauJet_Part_JetIndex",&tauJetPartJetIndex);

  tree_photons->SetBranchAddress("reco_Energy",&recoE);
  tree_photons->SetBranchAddress("reco_Px",&recoPx);
  tree_photons->SetBranchAddress("reco_Py",&recoPy);
  tree_photons->SetBranchAddress("reco_Pz",&recoPz);
  tree_photons->SetBranchAddress("reco_Phi",&recoPhi);
  tree_photons->SetBranchAddress("reco_CosTheta",&recoCosTheta);
  tree_photons->SetBranchAddress("reco_Theta",&recoTheta);
  tree_photons->SetBranchAddress("reco_Phi_logERW",&recoPhi_logERW);
  tree_photons->SetBranchAddress("reco_CosTheta_logERW",&recoCosTheta_logERW);
  tree_photons->SetBranchAddress("reco_Theta_logERW",&recoTheta_logERW);
  tree_photons->SetBranchAddress("reco_PDGID",&recoPDGID);
  tree_photons->SetBranchAddress("reco_Charge",&recoCharge);
  tree_photons->SetBranchAddress("reco_track0_p",&recotrack0p);
  tree_photons->SetBranchAddress("reco_track0_pt",&recotrack0pt);
  tree_photons->SetBranchAddress("reco_track0_chi2OverNdof",&recotrack0Chi2OverNdof);
  tree_photons->SetBranchAddress("reco_track0_nHits",&recotrack0nHits);
  tree_photons->SetBranchAddress("reco_nClusters",&recoNclusters);
  tree_photons->SetBranchAddress("reco_nTracks",&recoNtracks);

  tree_photons->SetBranchAddress("reco_E_EB",&recoEEB);
  tree_photons->SetBranchAddress("reco_E_EE",&recoEEE);
  tree_photons->SetBranchAddress("reco_E_EP",&recoEEP);
  tree_photons->SetBranchAddress("reco_E_HB",&recoEHB);
  tree_photons->SetBranchAddress("reco_E_HE",&recoEHE);
  tree_photons->SetBranchAddress("reco_E_HP",&recoEHP);
  tree_photons->SetBranchAddress("reco_nhitsEB",&recoNhitsEB);
  tree_photons->SetBranchAddress("reco_nhitsEE",&recoNhitsEE);
  tree_photons->SetBranchAddress("reco_nhitsEP",&recoNhitsEP);
  tree_photons->SetBranchAddress("reco_nhitsHB",&recoNhitsHB);
  tree_photons->SetBranchAddress("reco_nhitsHE",&recoNhitsHE);
  tree_photons->SetBranchAddress("reco_nhitsHP",&recoNhitsHP);

  TH1F* h_taum_true_1p_dPhi_reco = new TH1F("taum_true_1prong_dPhi_reco","",1000,-0.51,TMath::Pi());
  h_taum_true_1p_dPhi_reco->SetLineWidth(2);
  h_taum_true_1p_dPhi_reco->GetXaxis()->SetTitle("#Delta#Phi(#tau_{true},#tau_{reco})");
  TH1F* h_taup_true_1p_dPhi_reco = new TH1F("taup_true_1prong_dPhi_reco","",1000,-0.51,TMath::Pi());
  h_taup_true_1p_dPhi_reco->SetLineWidth(2);
  h_taup_true_1p_dPhi_reco->SetLineColor(kRed);
  TH1F* h_taum_true_3p_dPhi_reco = new TH1F("taum_true_3prong_dPhi_reco","",1000,-0.51,TMath::Pi());
  h_taum_true_3p_dPhi_reco->SetLineWidth(2);
  h_taum_true_3p_dPhi_reco->SetLineColor(kBlue);
  TH1F* h_taup_true_3p_dPhi_reco = new TH1F("taup_true_3prong_dPhi_reco","",1000,-0.51,TMath::Pi());
  h_taup_true_3p_dPhi_reco->SetLineWidth(2);
  h_taup_true_3p_dPhi_reco->SetLineColor(kMagenta);

  int n_bins100 = 100;
  float lim_phi_low = -TMath::Pi();
  float lim_phi_high = TMath::Pi();

  TH1F* h_taum_nneut_true_1prong_and_reco = new TH1F("taum_nneut_true_1prong_and_reco","",10,-0.5,9.5);
  TH1F* h_taum_nchother_true_1prong_and_reco = new TH1F("taum_nchother_true_1prong_and_reco","",10,-0.5,9.5);

  TH1F* h_taum_nph_true_1prong_and_reco = new TH1F("taum_nph_true_1prong_and_reco","",10,-0.5,9.5);
  TH1F* h_taum_npi0_true_1prong_and_reco = new TH1F("taum_npi0_true_1prong_and_reco","",10,-0.5,9.5);

  TH1F* h_taum_dphi_min_pim_else_true_1prong_and_reco = new TH1F("h_taum_dphi_min_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);
  TH1F* h_taum_dtheta_min_pim_else_true_1prong_and_reco = new TH1F("h_taum_dtheta_min_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);
  TH1F* h_taum_dangle_min_pim_else_true_1prong_and_reco = new TH1F("h_taum_dangle_min_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);

  TH1F* h_taum_dphi_max_pim_else_true_1prong_and_reco = new TH1F("h_taum_dphi_max_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);
  TH1F* h_taum_dtheta_max_pim_else_true_1prong_and_reco = new TH1F("h_taum_dtheta_max_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);
  TH1F* h_taum_dangle_max_pim_else_true_1prong_and_reco = new TH1F("h_taum_dangle_max_pim_else_true_1prong_and_reco","",n_bins100,0,0.05);

  TH1F* h_taum_jetmass_true_1prong_and_reco = new TH1F("h_taum_jetmass_neut_true_1prong_and_reco","",n_bins100,0,10.0);
  TH1F* h_taum_jetenergy_true_1prong_and_reco = new TH1F("h_taum_jetenergy_neut_true_1prong_and_reco","",10*n_bins100,0,1100.);

  TH1F* h_taum_pim_true_1prong_no_reco_dangleMin_pim = new TH1F("h_taum_pim_true_1prong_no_reco_dangleMin_pim","",100,0,TMath::Pi()+0.5);
  h_taum_pim_true_1prong_no_reco_dangleMin_pim->SetLineColor(kRed);
  h_taum_pim_true_1prong_no_reco_dangleMin_pim->GetXaxis()->SetTitle("#Delta#alpha(#pi^{-}_{true}, #pi^{#pm}_{true})");

  TH1F* h_taum_pim_true_1prong_no_reco_dangleMin_track = new TH1F("h_taum_pim_true_1prong_no_reco_dangleMin_track","",100,0,TMath::Pi()+0.5);
  h_taum_pim_true_1prong_no_reco_dangleMin_track->SetLineColor(kRed);
  h_taum_pim_true_1prong_no_reco_dangleMin_track->GetXaxis()->SetTitle("#Delta#alpha(#pi^{-}_{true},any track)");

  TH1F* h_taum_nneut_true_1prong_no_reco = new TH1F("taum_nneut_true_1prong_no_reco","",10,-0.5,9.5);
  h_taum_nneut_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_nneut_true_1prong_no_reco->GetXaxis()->SetTitle("N_{neutrals}");
  TH1F* h_taum_nchother_true_1prong_no_reco = new TH1F("taum_nchother_true_1prong_no_reco","",10,-0.5,9.5);
  h_taum_nchother_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_nchother_true_1prong_no_reco->GetXaxis()->SetTitle("N_{charged} (no #pi's)");
  TH1F* h_taum_nph_true_1prong_no_reco = new TH1F("taum_nph_true_1prong_no_reco","",10,-0.5,9.5);
  h_taum_nph_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_nph_true_1prong_no_reco->GetXaxis()->SetTitle("N_{photons}");
  TH1F* h_taum_npi0_true_1prong_no_reco = new TH1F("taum_npi0_true_1prong_no_reco","",10,-0.5,9.5);
  h_taum_npi0_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_npi0_true_1prong_no_reco->GetXaxis()->SetTitle("N_{pi^{0}}");

  TH1F* h_taum_dphi_min_pim_else_true_1prong_no_reco = new TH1F("h_taum_dphi_min_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dphi_min_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dphi_min_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#phi_{min}{pi^{-},else}");
  TH1F* h_taum_dtheta_min_pim_else_true_1prong_no_reco = new TH1F("h_taum_dtheta_min_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dtheta_min_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dtheta_min_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#theta_{min}{pi^{-},else}");
  TH1F* h_taum_dangle_min_pim_else_true_1prong_no_reco = new TH1F("h_taum_dangle_min_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dangle_min_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dangle_min_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#alpha_{min}{pi^{-},else}");
  TH1F* h_taum_dphi_max_pim_else_true_1prong_no_reco = new TH1F("h_taum_dphi_max_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dphi_max_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dphi_max_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#phi_{max}{pi^{-},else}");
  TH1F* h_taum_dtheta_max_pim_else_true_1prong_no_reco = new TH1F("h_taum_dtheta_max_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dtheta_max_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dtheta_max_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#theta_{max}{pi^{-},else}");
  TH1F* h_taum_dangle_max_pim_else_true_1prong_no_reco = new TH1F("h_taum_dangle_max_pim_else_true_1prong_no_reco","",n_bins100,0,0.05);
  h_taum_dangle_max_pim_else_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_dangle_max_pim_else_true_1prong_no_reco->GetXaxis()->SetTitle("#Delta#alpha_{max}{pi^{-},else}");
  TH1F* h_taum_jetmass_true_1prong_no_reco = new TH1F("h_taum_jetmass_neut_true_1prong_no_reco","",n_bins100,0,10.0);
  h_taum_jetmass_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_jetmass_true_1prong_no_reco->GetXaxis()->SetTitle("mass (#tau^{-}-jet)");
  TH1F* h_taum_jetenergy_true_1prong_no_reco = new TH1F("h_taum_jetenergy_neut_true_1prong_no_reco","",10*n_bins100,0,1100.);
  h_taum_jetenergy_true_1prong_no_reco->SetLineColor(kRed);
  h_taum_jetenergy_true_1prong_no_reco->GetXaxis()->SetTitle("energy (#tau^{-}-jet)");

  int n_bins2Deff = 20;

  float lim_theta_low = 80./180.*TMath::Pi();
  float lim_theta_high = 100./180.*TMath::Pi();

  TH2F* h_taum_nch_true_vs_nch_reco = new TH2F("taum_nch_true_vs_nch_reco","",10,-0.5,9.5,10,-0.5,9.5);
  TH2F* h_taum_npi_true_vs_npi_reco = new TH2F("taum_npi_true_vs_npi_reco","",10,-0.5,9.5,10,-0.5,9.5);
  //TH2F* taupLep_truePhiVsTheta = new TH2F("taum_NoPi_true_vs_NoPi_reco","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);

  TEfficiency* TaumLepCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("TaumLepCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* TaupLepCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("TaumLepCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TH2F* taumLep_truePhiVsTheta = new TH2F("taumLep_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taupLep_truePhiVsTheta = new TH2F("taupLep_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);

  //0.99 for acos theta corresponds to DR of 0.141539, around 8.1 degrees
  //checked vs tau which has neutrino in decay --> take out possible high energy brem photon
  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTruePhiTheta_ATJ","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTruePhiTheta_ATJ","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  //checked vs tau which has neutrino in decay --> take out possible high energy brem photon
  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTruePhiTheta_CH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTruePhiTheta_CH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetTitle("#tau^{-} 1prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTrueTheta = new TEfficiency("MatchedTaum1prongCosTheta0_99EfficiencyVSTrueTheta","",50,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTrueTheta = new TEfficiency("MatchedTaup1prongCosTheta0_99EfficiencyVSTrueTheta","",50,lim_theta_low,lim_theta_high);


  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhi = new TEfficiency("MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhi","",50,lim_phi_low,lim_phi_high);
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhi = new TEfficiency("MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhi","",50,lim_phi_low,lim_phi_high);

  TH2F* taum1prong_truePhiVsTheta = new TH2F("taum1prong_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taup1prong_truePhiVsTheta = new TH2F("taup1prong_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taum1prong_recoPhiVsTheta = new TH2F("taum1prong_recoPhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taup1prong_recoPhiVsTheta = new TH2F("taup1prong_recoPhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);

  //checked vs tau which has neutrino in decay --> take out possible high energy brem photon
  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta = new TEfficiency("MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ = new TEfficiency("MatchedTau3prongpCosTheta0_99EfficiencyVSTruePhiTheta_ATJ","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  //checked vs tau which has neutrino in decay --> take out possible high energy brem photon
  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH = new TEfficiency("MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetTitle("#tau^{-} 3prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH = new TEfficiency("MatchedTau3prongpCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetTitle("#tau^{+} 3prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");


  TH1F* h_taum1p_pi1_energy = new TH1F("h_taum1p_pi1_energy","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taum1p_pi1_energy->GetXaxis()->SetTitle("#pi_{1}^{true} energy [GeV]");
  TH1F* h_taup1p_pi1_energy = new TH1F("h_taup1p_pi1_energy","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taup1p_pi1_energy->SetLineColor(kRed);
  h_taup1p_pi1_energy->GetXaxis()->SetTitle("#pi_{1}^{true} energy [GeV]");
  TH1F* h_taum3p_pi1_energy = new TH1F("h_taum3p_pi1_energy","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taum3p_pi1_energy->SetLineColor(kBlue);
  h_taum3p_pi1_energy->GetXaxis()->SetTitle("#pi_{1}^{true} energy [GeV]");
  TH1F* h_taup3p_pi1_energy = new TH1F("h_taup3p_pi1_energy","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taup3p_pi1_energy->SetLineColor(kGreen+2);
  h_taup3p_pi1_energy->GetXaxis()->SetTitle("#pi_{1}^{true} energy [GeV]");

  TH1F* h_taum1p_pi1_pt = new TH1F("h_taum1p_pi1_pt","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taum1p_pi1_pt->GetXaxis()->SetTitle("#pi_{1}^{true} p_{T} [GeV]");
  TH1F* h_taup1p_pi1_pt = new TH1F("h_taup1p_pi1_pt","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taup1p_pi1_pt->SetLineColor(kRed);
  h_taup1p_pi1_pt->GetXaxis()->SetTitle("#pi_{1}^{true} p_{T} [GeV]");
  TH1F* h_taum3p_pi1_pt = new TH1F("h_taum3p_pi1_pt","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taum3p_pi1_pt->SetLineColor(kBlue);
  h_taum3p_pi1_pt->GetXaxis()->SetTitle("#pi_{1}^{true} p_{T} [GeV]");
  TH1F* h_taup3p_pi1_pt = new TH1F("h_taup3p_pi1_pt","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  h_taup3p_pi1_pt->SetLineColor(kGreen+2);
  h_taup3p_pi1_pt->GetXaxis()->SetTitle("#pi_{1}^{true} p_{T} [GeV]");


  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->SetTitle("#tau^{-} 3prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#pi_{1}^{true}  p_{T} [GeV];");
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH = new TEfficiency("MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->SetTitle("#tau^{+} 3prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#pi_{1}^{true}  p_{T} [GeV];");

  TEfficiency* MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH = new TEfficiency("MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->SetTitle("#tau^{-} 1prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#pi_{1}^{true}  p_{T} [GeV];");
  TEfficiency* MatchedTaup1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH = new TEfficiency("MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH","",n_bins2Deff,lim_energy_true_low,lim_energy_true_high);
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->SetTitle("#tau^{+} 1prong Efficiency #tau^{-}#tau^{+}, E=500 GeV;#pi_{1}^{true}  p_{T} [GeV];");
 
  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTrueTheta = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTrueTheta","",50,lim_theta_low,lim_theta_high);
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTrueTheta = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTrueTheta","",50,lim_theta_low,lim_theta_high);


  TEfficiency* MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhi = new TEfficiency("MatchedTaumCosTheta0_99EfficiencyVSTruePhi","",50,lim_phi_low,lim_phi_high);
  TEfficiency* MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhi = new TEfficiency("MatchedTaupCosTheta0_99EfficiencyVSTruePhi","",50,lim_phi_low,lim_phi_high);

  TEfficiency* TrueTaum3prongEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaum3prongEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  //TrueTaum3prongEfficiencyVSTruePhiTheta->SetTitle("#tau^{-} 3prong Efficiency Z#rightarrow#tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* TrueTaum1prongEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaum1prongEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  //TrueTaum1prongEfficiencyVSTruePhiTheta->SetTitle("#tau^{-} 1prong Efficiency Z#rightarrow#tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* TrueTaumLepElEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaumLepElEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* TrueTaumLepMuEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaumLepMuEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  TEfficiency* TrueTaup3prongEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaup3prongEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  //TrueTaup3prongEfficiencyVSTruePhiTheta->SetTitle("#tau^{+} 3prong Efficiency Z#rightarrow#tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* TrueTaup1prongEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaup1prongEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  //TrueTaup1prongEfficiencyVSTruePhiTheta->SetTitle("#tau^{+} 1prong Efficiency Z#rightarrow#tau^{-}#tau^{+}, E=500 GeV;#phi_{true};#theta_{true}");
  TEfficiency* TrueTaupLepElEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaupLepElEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);
  TEfficiency* TrueTaupLepMuEfficiencyVSTruePhiTheta = new TEfficiency("TrueTaupLepMuEfficiencyVSTruePhiTheta","",n_bins2Deff,lim_phi_low,lim_phi_high,n_bins2Deff,lim_theta_low,lim_theta_high);

  /*
  TH2F* taum3prong_truePhiVsTheta = new TH2F("taum3prong_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taup3prong_truePhiVsTheta = new TH2F("taup3prong_truePhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taum3prong_recoPhiVsTheta = new TH2F("taum3prong_recohiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  TH2F* taup3prong_recoPhiVsTheta = new TH2F("taup3prong_recoPhiVsTheta ","",50,lim_phi_low,lim_phi_high,50,lim_theta_low,lim_theta_high);
  */

  unsigned int taum3prong_true=0;
  unsigned int taum3prong_reco=0;
  unsigned int taum3prong_true_and_reco=0;

  unsigned int taum1prong_true=0;

  unsigned int taumlep_true=0;

  unsigned int taum1prong_reco=0;
  unsigned int taum1prong_true_and_reco=0;

  unsigned int taum1prong_true_and_any_reco=0;
  unsigned int taum3prong_true_and_any_reco=0;
  unsigned int taumlep_true_and_any_reco=0;
  unsigned int taum_any_reco=0;


  unsigned int count_1prong_taum_reco_standalone=0;
  unsigned int count_3prong_taum_reco_standalone=0;

  unsigned int taup3prong_true=0;
  unsigned int taup3prong_reco=0;
  unsigned int taup3prong_true_and_reco=0;

  unsigned int taup1prong_true=0;

  unsigned int tauplep_true=0;

  unsigned int taup1prong_reco=0;
  unsigned int taup1prong_true_and_reco=0;

  unsigned int taup1prong_true_and_any_reco=0;
  unsigned int taup3prong_true_and_any_reco=0;
  unsigned int tauplep_true_and_any_reco=0;
  unsigned int taup_any_reco=0;

  unsigned int count_1prong_taup_reco_standalone=0;
  unsigned int count_3prong_taup_reco_standalone=0;

  for(unsigned int i_entry=0;i_entry<tree_photons->GetEntries();i_entry++){
  //for(unsigned int i_entry=0;i_entry<500;i_entry++){
    tree_photons->GetEntry(i_entry);
    if(i_entry%45000==0){
      std::cout<<"entry "<<i_entry<<" tot "<<tree_photons->GetEntries()<<std::endl;
    }
    bool found_true_taum=false;
    bool found_true_taup=false;
    bool found_true_taum1prong=false;
    bool found_true_taup1prong=false;

    bool found_true_taum3prong=false;
    bool found_true_taup3prong=false;

    bool found_true_taumLepEl=false;
    bool found_true_taupLepEl=false;

    bool found_true_taumLepMu=false;
    bool found_true_taupLepMu=false;

    bool any_taum_reco=false;
    bool any_taup_reco=false;

    TLorentzVector TTauM_true(0,0,0,0);
    TLorentzVector TTauP_true(0,0,0,0);

    TLorentzVector TTauM_JetTrue(0,0,0,0);
    TLorentzVector TTauP_JetTrue(0,0,0,0);

    TLorentzVector TTauM_NuTrue(0,0,0,0);
    TLorentzVector TTauP_NuTrue(0,0,0,0);

    bool found_reco_taum=false;
    bool found_reco_taup=false;

    bool found_reco_taum1prong=false;
    bool found_reco_taup1prong=false;

    bool found_reco_taum1prong_CH=false;
    bool found_reco_taup1prong_CH=false;

    bool found_reco_taum1prong_ATJ=false;
    bool found_reco_taup1prong_ATJ=false;

    bool found_reco_taum1prong_ATJCH=false;
    bool found_reco_taup1prong_ATJCH=false;

    bool found_reco_taum3prong=false;
    bool found_reco_taup3prong=false;

    bool found_reco_taum3prong_CH=false;
    bool found_reco_taup3prong_CH=false;

    bool found_reco_taum3prong_ATJ=false;
    bool found_reco_taup3prong_ATJ=false;

    bool found_reco_taum3prong_ATJCH=false;
    bool found_reco_taup3prong_ATJCH=false;


    unsigned int count_charged_pions_taum=0;
    unsigned int count_charged_electrons_taum=0;
    unsigned int count_charged_muons_taum=0;
    unsigned int count_particles_taum=0;

    unsigned int count_charged_pions_taup=0;
    unsigned int count_charged_electrons_taup=0;
    unsigned int count_charged_muons_taup=0;
    unsigned int count_particles_taup=0;

    float taum_pion1pt=0;
    float taup_pion1pt=0;

    float taum_pion1E=0;
    float taup_pion1E=0;

    for(unsigned int m=0;m<truePhi->size();m++){
    //for(unsigned int m=0;m<500;m++){
      unsigned int count_charged_pions=0;
      unsigned int count_charged_electrons=0;
      unsigned int count_charged_muons=0;
      unsigned int count_particles=0;
      int tau_charge_pion_check=0;
      int tau_charge_lepton_check=0;
      bool neutrino_check=false;
      if(!found_true_taum){
	TTauM_JetTrue.SetPxPyPzE(0,0,0,0);
	TTauM_NuTrue.SetPxPyPzE(0,0,0,0);
      }
      if(!found_true_taup){
	TTauP_JetTrue.SetPxPyPzE(0,0,0,0);
	TTauP_NuTrue.SetPxPyPzE(0,0,0,0);
      }
      if(!found_true_taum && (*truePDGID)[m]==15 && (*trueGenStatus)[m]==2){
	for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	  //std::cout<<"at start of taum d loop "<<d<<std::endl;
	  TLorentzVector TMCParttemp(0,0,0,0);
	  TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	  bool daughter_double_count=false;
	  //for(unsigned int d1=0;d1<d;d1++){
	  //TLorentzVector TMCD1Parttemp(0,0,0,0);
	  //TMCD1Parttemp.SetPxPyPzE((*trueTauDaughterPx)[d1],(*trueTauDaughterPy)[d1],(*trueTauDaughterPz)[d1],(*trueTauDaughterE)[d1]);
	  //}
	  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
	    count_particles+=1;
	    count_particles_taum+=1;
	    if((*trueTauDaughterPDGID)[d]==11){
	      count_charged_electrons+=1;
	      count_charged_electrons_taum+=1;
	    }
	    if((*trueTauDaughterPDGID)[d]==13){
	      count_charged_muons+=1;
	      count_charged_muons_taum+=1;
	    }
	    if(abs((*trueTauDaughterPDGID)[d])==charged_pion_ID){
	      if((*trueTauDaughterE)[d]>taum_pion1E){
		taum_pion1E=(*trueTauDaughterE)[d];
		taum_pion1pt=sqrt((*trueTauDaughterPx)[d]*(*trueTauDaughterPx)[d]+(*trueTauDaughterPy)[d]*(*trueTauDaughterPy)[d]);
	      }
	      count_charged_pions+=1;
	      count_charged_pions_taum+=1;
	      tau_charge_pion_check+=(*trueTauDaughterCharge)[d];
	    }//take neutrinos out of MCtau_trueJetVector
	    else{
	      tau_charge_lepton_check+=(*trueTauDaughterCharge)[d];
	    }
	    if(abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
	      TTauM_JetTrue+=TMCParttemp;
	    }else if((*trueTauDaughterPDGID)[d]==16){
	      neutrino_check=true;
	      TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	      TTauM_NuTrue+=TMCParttemp;
	      //std::cout<<(*trueIndex)[m]<<" is real taum"<<std::endl;
	    }else{
	      TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	      TTauM_NuTrue+=TMCParttemp;
	    }
	  }
	}
	if(neutrino_check){
	  found_true_taum=true;
	  TTauM_true.SetPxPyPzE((*truePx)[m],(*truePy)[m],(*truePz)[m],(*trueE)[m]);
	  if(fabs((TTauM_NuTrue+TTauM_JetTrue).E()-TTauM_true.E())>(0.01*TTauM_true.E())){
	    std::cout<<"tauM and tau jet calc off "<<fabs((TTauM_NuTrue+TTauM_JetTrue).E()-TTauM_true.E())/TTauM_true.E()<<"/"<<(TTauM_NuTrue+TTauM_JetTrue).E()<<"/"<<TTauM_true.E()<<"/"<<TTauM_NuTrue.E()<<"/"<<TTauM_JetTrue.E()<<std::endl;
	  }
	  //taum_truePhiVsTheta->Fill(TTauM_true.Phi(),TTauM_true.Theta());
	  if(count_charged_pions==1){
	    found_true_taum1prong=true;
	  }else if (count_charged_pions==3){
	    found_true_taum3prong=true;
	  }
	  if(count_charged_muons==1){
	    found_true_taumLepMu=true;
	  }
	  if(count_charged_electrons==1){
	    found_true_taumLepEl=true;
	  }
	  /*
	  if(found_true_taumLepEl && found_true_taumLepMu){
	    std::cout<<"would expect only el or mu taum "<<found_true_taumLepEl<<"/"<<found_true_taumLepMu<<std::endl;
	    for(unsigned int m=0;m<truePhi->size();m++){
	      if((*truePDGID)[m]==-15){
		for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		    std::cout<<"taum "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		  }
	  	}
	      }
	    }
	  }
	  if(found_true_taumLepEl && (found_true_taum1prong || found_true_taum3prong)){
	    //std::cout<<"taum EL and 1prong or 3 prong "<<found_true_taumLepEl<<"/"<<found_true_taum1prong<<"/"<<found_true_taum3prong<<std::endl;
	    for(unsigned int m=0;m<truePhi->size();m++){
	      if((*truePDGID)[m]==-15){
		for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		    std::cout<<"taum "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		  }
		}
	      }
	    }
	  }
	  */
	}
	/*else{
	  std::cout<<"neutrino check of taum failed"<<std::endl;
	  for(unsigned int m=0;m<truePhi->size();m++){
	    if((*truePDGID)[m]==-15){
	      for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		  std::cout<<"taum "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		}
	      }
	    }
	  }
	  }*/
      }//tau minus code 

      //now go to tau plus
      if(!found_true_taup && (*truePDGID)[m]==-15 && (*trueGenStatus)[m]==2){
	for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
	    count_particles+=1;
	    count_particles_taup+=1;
	    if((*trueTauDaughterPDGID)[d]==11){
	      count_charged_electrons+=1;
	      count_charged_electrons_taup+=1;
	    }
	    if((*trueTauDaughterPDGID)[d]==13){
	      count_charged_muons+=1;
	      count_charged_muons_taup+=1;
	    }
	    if(abs((*trueTauDaughterPDGID)[d])==charged_pion_ID){
	      if((*trueTauDaughterE)[d]>taup_pion1E){
		taup_pion1E=(*trueTauDaughterE)[d];
		taup_pion1pt=sqrt((*trueTauDaughterPx)[d]*(*trueTauDaughterPx)[d]+(*trueTauDaughterPy)[d]*(*trueTauDaughterPy)[d]);
	      }
	      count_charged_pions+=1;
	      count_charged_pions_taup+=1;
	      tau_charge_pion_check+=(*trueTauDaughterCharge)[d];
	    }//take neutrinos out of MCtau_trueJetVector
	    else{
	      tau_charge_lepton_check+=(*trueTauDaughterCharge)[d];
	    }
	    TLorentzVector TMCParttemp(0,0,0,0);
	    if(abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
	      TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	      TTauP_JetTrue+=TMCParttemp;
	    }else if((*trueTauDaughterPDGID)[d]==-16){
	      neutrino_check=true;
	      TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	      TTauP_NuTrue+=TMCParttemp;
	    }else{
	      TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
	      TTauP_NuTrue+=TMCParttemp;
	    }
	  }
	}
	if(neutrino_check){
	  found_true_taup=true;
	  TTauP_true.SetPxPyPzE((*truePx)[m],(*truePy)[m],(*truePz)[m],(*trueE)[m]);
	  if(fabs((TTauP_NuTrue+TTauP_JetTrue).E()-TTauP_true.E())>(0.01*TTauP_true.E())){
	    std::cout<<"tauP and tau jet calc off "<<fabs((TTauP_NuTrue+TTauP_JetTrue).E()-TTauP_true.E())/TTauP_true.E()<<"/"<<(TTauP_NuTrue+TTauP_JetTrue).E()<<"/"<<TTauP_true.E()<<"/"<<TTauP_NuTrue.E()<<"/"<<TTauP_JetTrue.E()<<std::endl;
	  }
	  //taup_truePhiVsTheta->Fill(TTaup_true.Phi(),TTaup_true.Theta());
	  if(count_charged_pions==1){
	    found_true_taup1prong=true;
	  }else if (count_charged_pions==3){
	    found_true_taup3prong=true;
	  }
	  /*
	  if(found_true_taupLepEl && found_true_taupLepMu){
	    //std::cout<<"would expect only el or mu taup "<<found_true_taupLepEl<<"/"<<found_true_taupLepMu<<std::endl;
	    for(unsigned int m=0;m<truePhi->size();m++){
	      if((*truePDGID)[m]==15){
		for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		    std::cout<<"taup "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		  }
		}
	      }
	    }
	  }
	  if(found_true_taupLepEl && (found_true_taup1prong || found_true_taup3prong)){
	    std::cout<<"taup EL and 1prong or 3 prong "<<found_true_taupLepEl<<"/"<<found_true_taup1prong<<"/"<<found_true_taup3prong<<std::endl;
	    for(unsigned int m=0;m<truePhi->size();m++){
	      if((*truePDGID)[m]==15){
		for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		  if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		    std::cout<<"taup "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		  }
		}
	      }
	    }
	    }*/
	}
	/*else{
	  std::cout<<"neutrino check of taup failed"<<std::endl;
	  for(unsigned int m=0;m<truePhi->size();m++){
	    if((*truePDGID)[m]==15){
	      for(unsigned int d=0;d<trueTauDaughterE->size();d++){
		if((*trueTauDaughterTauIndex)[d]==(*trueIndex)[m]){
		  std::cout<<"taup "<<m<<" daughter "<<d<<" PDG "<<(*trueTauDaughterPDGID)[d]<<std::endl;
		}
	      }
	    }
	  }
	}
	*/
      }//tau plus code 
    }
    if(found_true_taupLepEl || found_true_taupLepMu){
      tauplep_true+=1;
      if(any_taup_reco){
	tauplep_true_and_any_reco+=1;	
      }
    }
    

    if(found_reco_taup1prong_ATJ){
      taup1prong_reco+=1;
    }




    if(!found_true_taup || !found_true_taum){
      std::cout<<"at end of true loop no taum/taup? "<<found_true_taum<<"/"<<found_true_taup<<std::endl;
    }
    if(found_true_taum){
      TrueTaum3prongEfficiencyVSTruePhiTheta->Fill(found_true_taum3prong,TTauM_true.Phi(),TTauM_true.Theta());
      TrueTaum1prongEfficiencyVSTruePhiTheta->Fill(found_true_taum1prong,TTauM_true.Phi(),TTauM_true.Theta());
      TrueTaumLepElEfficiencyVSTruePhiTheta->Fill(found_true_taumLepEl,TTauM_true.Phi(),TTauM_true.Theta());
      TrueTaumLepMuEfficiencyVSTruePhiTheta->Fill(found_true_taumLepMu,TTauM_true.Phi(),TTauM_true.Theta());
    }
   if(found_true_taup){
      TrueTaup3prongEfficiencyVSTruePhiTheta->Fill(found_true_taup3prong,TTauP_true.Phi(),TTauP_true.Theta());
      TrueTaup1prongEfficiencyVSTruePhiTheta->Fill(found_true_taup1prong,TTauP_true.Phi(),TTauP_true.Theta());
      TrueTaupLepElEfficiencyVSTruePhiTheta->Fill(found_true_taupLepEl,TTauP_true.Phi(),TTauP_true.Theta());
      TrueTaupLepMuEfficiencyVSTruePhiTheta->Fill(found_true_taupLepMu,TTauP_true.Phi(),TTauP_true.Theta());
    }
    if(!found_true_taup || !found_true_taum){
      std::cout<<"taup or taum not found, what the hell "<<found_true_taup<<"/"<<found_true_taum<<std::endl;
    }
    unsigned int count_charged_pions_taum_reco=0;
    unsigned int count_charged_electrons_taum_reco=0;
    unsigned int count_charged_muons_taum_reco=0;
    unsigned int count_particles_taum_reco=0;
    unsigned int count_charged_particles_taum_reco=0;

    unsigned int count_charged_pions_taup_reco=0;
    unsigned int count_charged_electrons_taup_reco=0;
    unsigned int count_charged_muons_taup_reco=0;
    unsigned int count_particles_taup_reco=0;
    unsigned int count_charged_particles_taup_reco=0;

    if(tauJetPx->size()>2){
      std::cout<<"tau jet real big "<<tauJetPx->size()<<std::endl;
    }
    std::vector<float>dPhi_tauP_true;
    std::vector<float>dPhi_tauM_true;
    for(unsigned int m=0;m<tauJetPx->size();m++){
      TLorentzVector TTauTestReco(0,0,0,0);
      TTauTestReco.SetPxPyPzE((*tauJetPx)[m],(*tauJetPy)[m],(*tauJetPz)[m],(*tauJetE)[m]);
      if((*tauJetCharge)[m]==-1){
	taum_any_reco+=1;
	any_taum_reco=true;
	//reco says we should be in taum
	for(unsigned int d=0;d<tauJetPartPDGID->size();d++){
	  if((*tauJetPartJetIndex)[d]==m){
	    if(abs((*tauJetPartPDGID)[d])==charged_pion_ID){
	      count_charged_pions_taum_reco+=1;
	      count_charged_particles_taum_reco+=1;
	    }
	    if(abs((*tauJetPartPDGID)[d])==11){
	      count_charged_electrons_taum_reco+=1;
	      count_charged_particles_taum_reco+=1;
	    }
	    if(abs((*tauJetPartPDGID)[d])==15){
	      count_charged_muons_taum_reco+=1;
	      count_charged_particles_taum_reco+=1;
	    }
	  }
	}
      }
      if((*tauJetCharge)[m]==1){
	taup_any_reco+=1;
	any_taup_reco=true;
	//reco says we should be in taum
	for(unsigned int d=0;d<tauJetPartPDGID->size();d++){
	  if((*tauJetPartJetIndex)[d]==m){
	    if(abs((*tauJetPartPDGID)[d])==charged_pion_ID){
	      count_charged_pions_taup_reco+=1;
	      count_charged_particles_taup_reco+=1;
	    }
	    if(abs((*tauJetPartPDGID)[d])==11){
	      count_charged_electrons_taup_reco+=1;
	      count_charged_particles_taup_reco+=1;
	    }
	    if(abs((*tauJetPartPDGID)[d])==15){
	      count_charged_muons_taup_reco+=1;
	      count_charged_particles_taup_reco+=1;
	    }
	  }
	}
      }

      double cosAngcheckM=((TTauM_true.Px()*TTauTestReco.Px()+TTauM_true.Py()*TTauTestReco.Py()+TTauTestReco.Pz()*TTauM_true.Pz())/(TTauM_true.P()*TTauTestReco.P()));
      double cosAngcheckP=((TTauP_true.Px()*TTauTestReco.Px()+TTauP_true.Py()*TTauTestReco.Py()+TTauTestReco.Pz()*TTauP_true.Pz())/(TTauP_true.P()*TTauTestReco.P()));

      double cosAngcheckM_TrueJet=((TTauM_JetTrue.Px()*TTauTestReco.Px()+TTauM_JetTrue.Py()*TTauTestReco.Py()+TTauTestReco.Pz()*TTauM_JetTrue.Pz())/(TTauM_JetTrue.P()*TTauTestReco.P()));
      double cosAngcheckP_TrueJet=((TTauP_JetTrue.Px()*TTauTestReco.Px()+TTauP_JetTrue.Py()*TTauTestReco.Py()+TTauTestReco.Pz()*TTauP_JetTrue.Pz())/(TTauP_JetTrue.P()*TTauTestReco.P()));
      
      dPhi_tauP_true.push_back(acos(cosAngcheckP_TrueJet));
      dPhi_tauM_true.push_back(acos(cosAngcheckM_TrueJet));


      unsigned int count_charged_pions=0;
      int jet_charge_calc=0;//BUG in charge count
      for(unsigned int d=0;d<tauJetPartPDGID->size();d++){
	if((*tauJetPartJetIndex)[d]==m){
	  jet_charge_calc+=(*tauJetPartCharge)[d];
	  if(abs((*tauJetPartPDGID)[d])==charged_pion_ID){
	    //std::cout<<"jet "<<m<<" pion "<<(*tauJetPartPDGID)[d]<<" charge "<<(*tauJetPartCharge)[d]<< "jet charge "<<(*tauJetCharge)[m]<<std::endl;
	    count_charged_pions+=1;
	  }//else{
	  //std::cout<<"jet "<<m<<" no pion "<<(*tauJetPartPDGID)[d]<<" charge "<<(*tauJetPartCharge)[d]<< "jet charge "<<(*tauJetCharge)[m]<<std::endl;
	  //}
	}
      }
      if(count_charged_pions==3){//check 3 prongs
	//std::cout<<"we are in 3 prong "<<std::endl;
	if(cosAngcheckM>0.995){
	  //std::cout<<"we are in 3 prong check minus passed "<<std::endl;
	  found_reco_taum3prong=true;
	  //if((*tauJetCharge)[m]==-1){
	  if(jet_charge_calc==-1){
	    //std::cout<<"we are in 3 prong check minus passed and charge too "<<std::endl;
	    found_reco_taum3prong_CH=true;
	  }
	}
	if(cosAngcheckM_TrueJet>0.995){
	  //std::cout<<"we are in 3 prong check minus passed for true tau jet angle"<<std::endl;
	  found_reco_taum3prong_ATJ=true;
	  //if((*tauJetCharge)[m]==-1){
	  if(jet_charge_calc==-1){
	    //std::cout<<"we are in 3 prong check minus passed for true tau jet angle and charge too"<<std::endl;
	    found_reco_taum3prong_ATJCH=true;
	  }
	}
	if(cosAngcheckP>0.995){
	  found_reco_taup3prong=true;
	  //if((*tauJetCharge)[m]==1){
	  if(jet_charge_calc==1){
	    found_reco_taup3prong_CH=true;
	  }
	}
	if(cosAngcheckP_TrueJet>0.995){
	  found_reco_taup3prong_ATJ=true;
	  //if((*tauJetCharge)[m]==1){
	  if(jet_charge_calc==1){
	    found_reco_taup3prong_ATJCH=true;
	  }
	}
      }
      if(count_charged_pions==1){//check 1 prongs
	//now for tauminus
	if(cosAngcheckM>0.995){
	  found_reco_taum1prong=true;
	  if(jet_charge_calc==-1){
	    //if(count_charges_pions!=(*tauJetCharge)[m]){
	    //  std::cout<<"taum 1p charge from pi's positive "<<count_charges_pions<<std::endl;
	    //}
	    found_reco_taum1prong_CH=true;
	  }
	}
	if(cosAngcheckM_TrueJet>0.995){
	  found_reco_taum1prong_ATJ=true;
	  //if((*tauJetCharge)[m]==-1){
	  if(jet_charge_calc==-1){
	    found_reco_taum1prong_ATJCH=true;
	  }
	}
	//now for tauplus
	if(cosAngcheckP>0.995){
	  found_reco_taup1prong=true;
	  //if((*tauJetCharge)[m]==1){
	  if(jet_charge_calc==1){
	    //std::cout<<"taup 1p charge from pi's negative "<<count_charges_pions<<std::endl;
	    found_reco_taup1prong_CH=true;
	  }
	}
	if(cosAngcheckP_TrueJet>0.995){
	  found_reco_taup1prong_ATJ=true;
	  //if((*tauJetCharge)[m]==1){
	  if(jet_charge_calc==1){
	    found_reco_taup1prong_ATJCH=true;
	  }	
	}
	//std::cout<<"found 1 prong charge calc from pi's "<<jet_charge_calc<<" cos min/minJ "<<cosAngcheckM<<"/"<<cosAngcheckM_TrueJet<<" cos P/Pjet "<<cosAngcheckP<<"/"<<cosAngcheckP_TrueJet<<std::endl;
	//std::cout<<" true P/M 1 and 3 prongs "<<found_true_taup1prong<<"/"<<found_true_taum1prong<<"/"<<found_true_taup3prong<<"/"<<found_true_taum3prong<<std::endl;
	//std::cout<<"m/p "<<found_reco_taum1prong<<"/"<<found_reco_taup1prong<<" mCH/pCH "<<found_reco_taum1prong_CH<<"/"<<found_reco_taup1prong_CH<<std::endl;
      }
    }
    if(dPhi_tauP_true.size()!=dPhi_tauM_true.size()){
      std::cout<<"sizes of angle vectors different, strange "<<dPhi_tauP_true.size()<<"/"<<dPhi_tauM_true.size()<<std::endl;
    }
    if(dPhi_tauP_true.size()==0){
      if(found_true_taum1prong){
	h_taum_true_1p_dPhi_reco->Fill(-0.5);
      }
      if(found_true_taum3prong){
	h_taum_true_3p_dPhi_reco->Fill(-0.5);
      }
      if(found_true_taup1prong){
	h_taup_true_1p_dPhi_reco->Fill(-0.5);
      }
      if(found_true_taup3prong){
	h_taup_true_3p_dPhi_reco->Fill(-0.5);
      }
    }else if(dPhi_tauP_true.size()==1){
      if(dPhi_tauP_true[0]<dPhi_tauM_true[0]){
	if(found_true_taup1prong){
	  h_taup_true_1p_dPhi_reco->Fill(dPhi_tauP_true[0]);
	}
	if(found_true_taup3prong){
	  h_taup_true_3p_dPhi_reco->Fill(dPhi_tauP_true[0]);
	}
	if(found_true_taum1prong){
	  h_taum_true_1p_dPhi_reco->Fill(-0.25);
	}
	if(found_true_taum3prong){
	  h_taum_true_3p_dPhi_reco->Fill(-0.25);
	}
      }else{
	if(found_true_taum1prong){
	  h_taum_true_1p_dPhi_reco->Fill(dPhi_tauM_true[0]);
	}
	if(found_true_taum3prong){
	  h_taum_true_3p_dPhi_reco->Fill(dPhi_tauM_true[0]);
	}
	if(found_true_taup1prong){
	  h_taup_true_1p_dPhi_reco->Fill(-0.25);
	}
	if(found_true_taup3prong){
	  h_taup_true_3p_dPhi_reco->Fill(-0.25);
	}
      }
    }else if(dPhi_tauP_true.size()==2){
      if(dPhi_tauP_true[0]<dPhi_tauM_true[0]){
	if(found_true_taup1prong){
	  h_taup_true_1p_dPhi_reco->Fill(dPhi_tauP_true[0]);
	}
	if(found_true_taup3prong){
	  h_taup_true_3p_dPhi_reco->Fill(dPhi_tauP_true[0]);
	}
	if(found_true_taum1prong){
	  h_taum_true_1p_dPhi_reco->Fill(dPhi_tauM_true[1]);
	}
	if(found_true_taum3prong){
	  h_taum_true_3p_dPhi_reco->Fill(dPhi_tauM_true[1]);
	}
      }else{
	if(found_true_taup1prong){
	  h_taup_true_1p_dPhi_reco->Fill(dPhi_tauP_true[1]);
	}
	if(found_true_taup3prong){
	  h_taup_true_3p_dPhi_reco->Fill(dPhi_tauP_true[1]);
	}
	if(found_true_taum1prong){
	  h_taum_true_1p_dPhi_reco->Fill(dPhi_tauM_true[0]);
	}
	if(found_true_taum3prong){
	  h_taum_true_3p_dPhi_reco->Fill(dPhi_tauM_true[0]);
	}
      }
    }

    if(found_true_taumLepEl || found_true_taumLepMu){
      taumlep_true+=1;
      if(any_taum_reco){
	taumlep_true_and_any_reco+=1;	
      }
    }
    

    if(count_charged_pions_taum_reco==1){
      count_1prong_taum_reco_standalone+=1;
    }else if(count_charged_pions_taum_reco==3){
      count_3prong_taum_reco_standalone+=1;
    }

    if(count_charged_pions_taup_reco==1){
      count_1prong_taup_reco_standalone+=1;
    }else if(count_charged_pions_taup_reco==3){
      count_3prong_taup_reco_standalone+=1;
    }


    if(found_reco_taum1prong_ATJ){
      taum1prong_reco+=1;
    }

    if(found_true_taum1prong || found_true_taum3prong){
      if(any_taum_reco){
	h_taum_nch_true_vs_nch_reco->Fill(count_charged_pions_taum+count_charged_electrons_taum+count_charged_muons_taum,count_charged_particles_taum_reco); 
	h_taum_npi_true_vs_npi_reco->Fill(count_charged_pions_taum,count_charged_pions_taum_reco); 
      }
    }

    if(found_true_taum1prong){
      taum1prong_true+=1;
      if(any_taum_reco){
	taum1prong_true_and_any_reco+=1;
	TLorentzVector TTauM_GenJetTrue(0,0,0,0);
	TLorentzVector TTauM_PimTrue(0,0,0,0);
	bool found_true_taum_fill=false;
	for(unsigned int m=0;m<truePhi->size();m++){
	  bool neutrino_check_d=false;
	  if(!found_true_taum_fill && (*truePDGID)[m]==15 && (*trueGenStatus)[m]==2){
	    for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	      if((*trueIndex)[m]==(*trueTauDaughterTauIndex)[d]){
		if(abs((*trueTauDaughterPDGID)[d])==211){
		  TTauM_PimTrue.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		}
		if((*trueTauDaughterPDGID)[d]==16){
		  neutrino_check_d=true;
		}
	      }
	    }
	  }
	  if(neutrino_check_d){
	    found_true_taum_fill=true;
	  }
	}
	found_true_taum_fill=false;
	for(unsigned int m=0;m<truePhi->size();m++){
	  std::set<float> tau_daughters_Pi0_Energies;
	  unsigned int count_charged_pi0=0;
	  unsigned int count_charged_other=0;
	  unsigned int count_neutrals=0;//not counting neutrinos
	  unsigned int count_photons=0;
	  unsigned int count_particles=0;
	  float dangle_min=2*TMath::Pi();
	  float dphi_min=2*TMath::Pi();
	  float dtheta_min=2*TMath::Pi();
	  float dangle_max=-2*TMath::Pi();
	  float dphi_max=-2*TMath::Pi();
	  float dtheta_max=-2*TMath::Pi();
	  bool neutrino_check_d=false;
	  if(!found_true_taum_fill){
	    TTauM_GenJetTrue.SetPxPyPzE(0,0,0,0);
	  }
	  if(!found_true_taum_fill && (*truePDGID)[m]==15 && (*trueGenStatus)[m]==2){
	    for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	      if((*trueIndex)[m]==(*trueTauDaughterTauIndex)[d]){
		TLorentzVector TMCParttemp(0,0,0,0);
		TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		//exclude neutrinos and pion itself
		if(abs((*trueTauDaughterPDGID)[d])!=211 && abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  if(fabs(TMCParttemp.DeltaPhi(TTauM_PimTrue))<dphi_min){
		    dphi_min=TMCParttemp.DeltaPhi(TTauM_PimTrue);
		  }else if(fabs(TMCParttemp.DeltaPhi(TTauM_PimTrue))>dphi_max){
		    dphi_max=TMCParttemp.DeltaPhi(TTauM_PimTrue);
		  }
		  if(fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta())<dtheta_min){
		    dtheta_min=fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta());
		  }else if(fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta())>dtheta_max){
		    dtheta_max=fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta());
		  }
		  float dangle=fabs(TMCParttemp.Angle(TTauM_PimTrue.Vect()));
		  if(dangle<dangle_min){
		    dangle_min=dangle;
		  }else if(dangle>dangle_max){
		    dangle_max=dangle;
		  }
		}
		if((*trueTauDaughterMotherPDGID)[d]==111){
		  tau_daughters_Pi0_Energies.insert((*trueTauDaughterMotherEnergy)[d]);
		}
		if(abs((*trueTauDaughterPDGID)[d])==211){
		  TTauM_PimTrue.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		}
		if((*trueTauDaughterCharge)[d]!=0 && abs((*trueTauDaughterPDGID)[d])!=211){
		  count_charged_other+=1;
		}else if( (*trueTauDaughterCharge)[d]==0 && abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  count_neutrals+=1;
		  if((*trueTauDaughterPDGID)[d]==22){
		    count_photons+=1;
		  }
		}
		if(abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  TTauM_GenJetTrue+=TMCParttemp;
		}else if((*trueTauDaughterPDGID)[d]==16){
		  neutrino_check_d=true;
		}
	      }
	    }
	    if(neutrino_check_d){
	      found_true_taum_fill=true;
	      h_taum_nph_true_1prong_and_reco->Fill(count_photons);
	      h_taum_nneut_true_1prong_and_reco->Fill(count_neutrals);
	      h_taum_nchother_true_1prong_and_reco->Fill(count_charged_other);
	      h_taum_npi0_true_1prong_and_reco->Fill(tau_daughters_Pi0_Energies.size());
	      h_taum_jetmass_true_1prong_and_reco->Fill(TTauM_GenJetTrue.M());
	      h_taum_jetenergy_true_1prong_and_reco->Fill(TTauM_GenJetTrue.Energy());
	      if(dphi_min<(2*TMath::Pi())){
		h_taum_dphi_min_pim_else_true_1prong_and_reco->Fill(dphi_min);
		h_taum_dtheta_min_pim_else_true_1prong_and_reco->Fill(dtheta_min);
		h_taum_dangle_min_pim_else_true_1prong_and_reco->Fill(dangle_min);
		h_taum_dphi_max_pim_else_true_1prong_and_reco->Fill(dphi_max);
		h_taum_dtheta_max_pim_else_true_1prong_and_reco->Fill(dtheta_max);
		h_taum_dangle_max_pim_else_true_1prong_and_reco->Fill(dangle_max);
	      }
	    }else{
	      std::cout<<"fail filling neutrino check wth reco and true"<<std::endl;
	    }
	  }
	}//1 prong found and any type of reco
      }else{//not ANY taum found
	TLorentzVector TTauM_PimTrue(0,0,0,0);
	bool found_true_taum_fill=false;
	for(unsigned int m=0;m<truePhi->size();m++){
	  bool neutrino_check_d=false;
	  if(!found_true_taum_fill && (*truePDGID)[m]==15 && (*trueGenStatus)[m]==2){
	    for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	      if((*trueIndex)[m]==(*trueTauDaughterTauIndex)[d]){
		if(abs((*trueTauDaughterPDGID)[d])==211){
		  TTauM_PimTrue.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		}
		if((*trueTauDaughterPDGID)[d]==16){
		  neutrino_check_d=true;
		}
	      }
	    }
	  }
	  if(neutrino_check_d){
	    found_true_taum_fill=true;
	  }
	}
	float dangleMin_pi_reco=3.5;
	float dangleMin_track=3.5;

	TLorentzVector TTauM_PimReco(0,0,0,0);
	for(unsigned int m=0;m<recoPhi->size();m++){
	  TTauM_PimReco.SetPxPyPzE((*recoPx)[m],(*recoPy)[m],(*recoPz)[m],(*recoE)[m]);
	  if(abs((*recoPDGID)[m])==211){
	    if(TTauM_PimReco.Angle(TTauM_PimTrue.Vect())<dangleMin_pi_reco){
	      dangleMin_pi_reco=TTauM_PimReco.Angle(TTauM_PimTrue.Vect());
	    }
	  }
	}
	h_taum_pim_true_1prong_no_reco_dangleMin_pim->Fill(dangleMin_pi_reco);

	/*
    if(found_true_taum1prong && !found_reco_taum1prong_ATJ){
      TH1F* h_taum_pim_true_1prong_no_reco_dangleMin_pim = new TH1F("h_taum_pim_true_1prong_no_reco_dangleMin_pim","",100,0,,TMath::Pi()+0.5);
      h_taum_pim_true_1prong_no_reco_dangleMin_pim->SetLineColor(kRed);
      h_taum_pim_true_1prong_no_reco_dangleMin_pim->GetXaxis()->SetTitle("#Delta#alpha(#pi^{-}_{true}, #pi^{#pm}_{true})");
      
      TH1F* h_taum_pim_true_1prong_no_reco_dangleMin_track = new TH1F("h_taum_pim_true_1prong_no_reco_dangleMin_track","",100,0,,TMath::Pi()+0.5);
      h_taum_pim_true_1prong_no_reco_dangleMin_track->SetLineColor(kRed);
      h_taum_pim_true_1prong_no_reco_dangleMin_track->GetXaxis()->SetTitle("#Delta#alpha(#pi^{-}_{true},any track)");
      

      h_taum_pim_true_1prong_no_reco_dangleMin_pim
  }
	*/



	found_true_taum_fill=false;
	TLorentzVector TTauM_GenJetTrue(0,0,0,0);
	for(unsigned int m=0;m<truePhi->size();m++){
	  std::set<float> tau_daughters_Pi0_Energies;
	  //for(unsigned int m=0;m<500;m++){
	  unsigned int count_charged_pi0=0;
	  unsigned int count_charged_other=0;
	  unsigned int count_neutrals=0;//not counting neutrinosx
	  unsigned int count_photons=0;
	  unsigned int count_particles=0;
	  float dangle_min=2*TMath::Pi();
	  float dphi_min=2*TMath::Pi();
	  float dtheta_min=2*TMath::Pi();
	  float dangle_max=-2*TMath::Pi();
	  float dphi_max=-2*TMath::Pi();
	  float dtheta_max=-2*TMath::Pi();
	  bool neutrino_check_d=false;
	  if(!found_true_taum_fill){
	    TTauM_GenJetTrue.SetPxPyPzE(0,0,0,0);
	  }
	  if(!found_true_taum_fill && (*truePDGID)[m]==15 && (*trueGenStatus)[m]==2){
	    for(unsigned int d=0;d<trueTauDaughterE->size();d++){
	      if((*trueIndex)[m]==(*trueTauDaughterTauIndex)[d]){
		if(abs((*trueTauDaughterPDGID)[d])==211){
		  TTauM_PimTrue.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		}
		TLorentzVector TMCParttemp(0,0,0,0);
		TMCParttemp.SetPxPyPzE((*trueTauDaughterPx)[d],(*trueTauDaughterPy)[d],(*trueTauDaughterPz)[d],(*trueTauDaughterE)[d]);
		//exclude neutrinos and pion itself
		if(abs((*trueTauDaughterPDGID)[d])!=211 && abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  if(fabs(TMCParttemp.DeltaPhi(TTauM_PimTrue))<dphi_min){
		    dphi_min=TMCParttemp.DeltaPhi(TTauM_PimTrue);
		  }else if(fabs(TMCParttemp.DeltaPhi(TTauM_PimTrue))>dphi_max){
		    dphi_max=TMCParttemp.DeltaPhi(TTauM_PimTrue);
		  }
		  if(fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta())<dtheta_min){
		    dtheta_min=fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta());
		  }else if(fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta())>dtheta_max){
		    dtheta_max=fabs(TMCParttemp.Theta()-TTauM_PimTrue.Theta());
		  }
		  float dangle=fabs(TMCParttemp.Angle(TTauM_PimTrue.Vect()));
		  if(dangle<dangle_min){
		    dangle_min=dangle;
		  }else if(dangle>dangle_max){
		    dangle_max=dangle;
		  }
		}
		if((*trueTauDaughterMotherPDGID)[d]==111){
		  tau_daughters_Pi0_Energies.insert((*trueTauDaughterMotherEnergy)[d]);
		}
		if((*trueTauDaughterCharge)[d]!=0 && abs((*trueTauDaughterPDGID)[d])!=211){
		  count_charged_other+=1;
		}else if( (*trueTauDaughterCharge)[d]==0 && abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  count_neutrals+=1;
		  if((*trueTauDaughterPDGID)[d]==22){
		    count_photons+=1;
		  }
		}
		if(abs((*trueTauDaughterPDGID)[d])!=12 && abs((*trueTauDaughterPDGID)[d])!=14 && abs((*trueTauDaughterPDGID)[d])!=16){
		  TTauM_GenJetTrue+=TMCParttemp;
		}else if((*trueTauDaughterPDGID)[d]==16){
		  neutrino_check_d=true;
		}
	      }	 
	    }
	    if(neutrino_check_d){
	      found_true_taum_fill=true;
	      h_taum_nph_true_1prong_no_reco->Fill(count_photons);
	      h_taum_nneut_true_1prong_no_reco->Fill(count_neutrals);
	      h_taum_nchother_true_1prong_no_reco->Fill(count_charged_other);
	      h_taum_npi0_true_1prong_no_reco->Fill(tau_daughters_Pi0_Energies.size());
	      h_taum_jetmass_true_1prong_no_reco->Fill(TTauM_GenJetTrue.M());
	      h_taum_jetenergy_true_1prong_no_reco->Fill(TTauM_GenJetTrue.Energy());
	      if(dphi_min<(2*TMath::Pi())){
		h_taum_dphi_min_pim_else_true_1prong_no_reco->Fill(dphi_min);
		h_taum_dtheta_min_pim_else_true_1prong_no_reco->Fill(dtheta_min);
	      h_taum_dangle_min_pim_else_true_1prong_no_reco->Fill(dangle_min);
	      h_taum_dphi_max_pim_else_true_1prong_no_reco->Fill(dphi_max);
	      h_taum_dtheta_max_pim_else_true_1prong_no_reco->Fill(dtheta_max);
	      h_taum_dangle_max_pim_else_true_1prong_no_reco->Fill(dangle_max);
	      }
	    }else{
	      std::cout<<"fail filling neutrino check wth "<<std::endl;
	    }
	  }
	}
      }//1 prong but not any reco
      //1 prong found
      if(found_reco_taum1prong){
	taum1prong_true_and_reco+=1;
      }
    }
    if(found_reco_taum3prong_ATJ){
      taum3prong_reco+=1;
    }

    if(found_reco_taup3prong_ATJ){
      taup3prong_reco+=1;
    }

   if(found_true_taup3prong){
      taup3prong_true+=1;
      if(any_taup_reco){
	taup3prong_true_and_any_reco+=1;
      }
      if(found_reco_taup3prong){
	taup3prong_true_and_reco+=1;
      }
    }

    if(found_true_taum3prong){
      taum3prong_true+=1;
      if(any_taum_reco){
	taum3prong_true_and_any_reco+=1;
      }
      if(found_reco_taum3prong){
	taum3prong_true_and_reco+=1;
      }
    }
    if(found_true_taum1prong){
      h_taum1p_pi1_energy->Fill(taum_pion1E);
      h_taum1p_pi1_pt->Fill(taum_pion1pt);
      MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH ->Fill(found_reco_taum1prong,taum_pion1pt);
      //std::cout<<"true taum1, reco m1/1AJ, m3, 3AJ "<<found_reco_taum1prong<<"/"<<found_reco_taum1prong_ATJ<<"/"<<found_reco_taum3prong<<"/"<<found_reco_taum3prong_ATJ<<std::endl;
      MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta->Fill(found_reco_taum1prong,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Fill(found_reco_taum1prong_ATJ,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Fill(found_reco_taum1prong_CH,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Fill(found_reco_taum1prong_ATJCH,TTauM_true.Phi(),TTauM_true.Theta());
    }
    if(found_true_taum3prong){
      h_taum3p_pi1_energy->Fill(taum_pion1E);
      h_taum3p_pi1_pt->Fill(taum_pion1pt);
      MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH ->Fill(found_reco_taum3prong,taum_pion1pt);
      //std::cout<<"true taum3, reco m1/1AJ, m3, 3AJ "<<found_reco_taum1prong<<"/"<<found_reco_taum1prong_ATJ<<"/"<<found_reco_taum3prong<<"/"<<found_reco_taum3prong_ATJ<<std::endl;
      MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta->Fill(found_reco_taum3prong,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Fill(found_reco_taum3prong_ATJ,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Fill(found_reco_taum3prong_CH,TTauM_true.Phi(),TTauM_true.Theta());
      MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Fill(found_reco_taum3prong_ATJCH,TTauM_true.Phi(),TTauM_true.Theta());
    }
    if(found_true_taup1prong){
      h_taup1p_pi1_energy->Fill(taup_pion1E);
      h_taup1p_pi1_pt->Fill(taup_pion1pt);
      MatchedTaup1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH ->Fill(found_reco_taup1prong,taup_pion1pt);
      MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta->Fill(found_reco_taup1prong_CH,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Fill(found_reco_taup1prong_ATJ,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Fill(found_reco_taup1prong_CH,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Fill(found_reco_taup1prong_ATJCH,TTauP_true.Phi(),TTauP_true.Theta());
    }
    if(found_true_taup3prong){
      h_taup3p_pi1_energy->Fill(taup_pion1E);
      h_taup3p_pi1_pt->Fill(taup_pion1pt);
      MatchedTaup3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH ->Fill(found_reco_taup3prong,taup_pion1pt);
      MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta->Fill(found_reco_taup3prong_CH,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Fill(found_reco_taup3prong_ATJ,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Fill(found_reco_taup3prong_CH,TTauP_true.Phi(),TTauP_true.Theta());
      MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Fill(found_reco_taup3prong_ATJCH,TTauP_true.Phi(),TTauP_true.Theta());
    }
  }
  std::cout<<"taum1p_t/1p_r/1p_tr "<<taum1prong_true<<"/"<<taum1prong_reco<<"/"<<taum1prong_true_and_reco<<" RECO SA "<< count_1prong_taum_reco_standalone<<" true 1p and any RECO "<<taum1prong_true_and_any_reco<<std::endl;
  std::cout<<"taum3p_t/3p_r/3p_tr "<<taum3prong_true<<"/"<<taum3prong_reco<<"/"<<taum3prong_true_and_reco<<" RECO SA "<< count_3prong_taum_reco_standalone<<" true 3p and any RECO "<<taum3prong_true_and_any_reco<<std::endl;

  std::cout<<"taumlep t/r "<< taumlep_true<<"/"<< taumlep_true_and_any_reco<< "rat 3p/1p/l "<<(float)taum3prong_true/(float)tree_photons->GetEntries()<<"/"<< (float)taum1prong_true/(float)tree_photons->GetEntries()<<"/"<<(float)taumlep_true/(float)tree_photons->GetEntries() <<" any taum reco "<<taum_any_reco<<" total entries "<<tree_photons->GetEntries()<<std::endl;

  std::cout<<"taup1p_t/1p_r/1p_tr "<<taup1prong_true<<"/"<<taup1prong_reco<<"/"<<taup1prong_true_and_reco<<" RECO SA "<< count_1prong_taup_reco_standalone<<" true 1p and any RECO "<<taup1prong_true_and_any_reco<<std::endl;
  std::cout<<"taup3p_t/3p_r/3p_tr "<<taup3prong_true<<"/"<<taup3prong_reco<<"/"<<taup3prong_true_and_reco<<" RECO SA "<< count_3prong_taup_reco_standalone<<" true 3p and any RECO "<<taup3prong_true_and_any_reco<<std::endl;

  std::cout<<"tauplep t/r "<< tauplep_true<<"/"<< tauplep_true_and_any_reco<< "rat 3p/1p/l "<<(float)taup3prong_true/(float)tree_photons->GetEntries()<<"/"<< (float)taup1prong_true/(float)tree_photons->GetEntries()<<"/"<<(float)tauplep_true/(float)tree_photons->GetEntries() <<" any taup reco "<<taup_any_reco<<" total entries "<<tree_photons->GetEntries()<<std::endl;



  TCanvas* can_taupeff1prongVsPhiTheta=setUpperCanvas("can_taupeff1prongVsPhiTheta");
  can_taupeff1prongVsPhiTheta->cd();
  can_taupeff1prongVsPhiTheta->SetRightMargin(0.18);
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taupeff1prongVsPhiTheta_ATJ=setUpperCanvas("can_taupeff1prongVsPhiTheta_ATJ");
  can_taupeff1prongVsPhiTheta_ATJ->cd();
  can_taupeff1prongVsPhiTheta_ATJ->SetRightMargin(0.18);
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Draw("colz");

  TCanvas* can_taupeff1prongVsPhiTheta_CH=setUpperCanvas("can_taupeff1prongVsPhiTheta_CH");
  can_taupeff1prongVsPhiTheta_CH->cd();
  can_taupeff1prongVsPhiTheta_CH->SetRightMargin(0.18);
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Draw("colz");

  TCanvas* can_taupeff1prongVsPhiTheta_ATJCH=setUpperCanvas("can_taupeff1prongVsPhiTheta_ATJCH");
  can_taupeff1prongVsPhiTheta_ATJCH->cd();
  can_taupeff1prongVsPhiTheta_ATJCH->SetRightMargin(0.18);
  //MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMinimum(0.60);
  //MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMaximum(0.90);
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Draw("colz");

  TCanvas* can_taupeff3prongVsPhiTheta=setUpperCanvas("can_taupeff3prongVsPhiTheta");
  can_taupeff3prongVsPhiTheta->cd();
  can_taupeff3prongVsPhiTheta->SetRightMargin(0.18);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taupeff3prongVsPhiTheta_ATJ=setUpperCanvas("can_taupeff3prongVsPhiTheta_ATJ");
  can_taupeff3prongVsPhiTheta_ATJ->cd();
  can_taupeff3prongVsPhiTheta_ATJ->SetRightMargin(0.18);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Draw("colz");

  TCanvas* can_taupeff3prongVsPhiTheta_CH=setUpperCanvas("can_taupeff3prongVsPhiTheta_CH");
  can_taupeff3prongVsPhiTheta_CH->cd();
  can_taupeff3prongVsPhiTheta_CH->SetRightMargin(0.18);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Draw("colz");

  TCanvas* can_taupeff3prongVsPhiTheta_ATJCH=setUpperCanvas("can_taupeff3prongVsPhiTheta_ATJCH");
  can_taupeff3prongVsPhiTheta_ATJCH->cd();
  can_taupeff3prongVsPhiTheta_ATJCH->SetRightMargin(0.18);
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Draw("colz");
  //MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMinimum(0.20);
  //MatchedTaup3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMaximum(0.80);

  TCanvas* can_taumeff1prongVsPhiTheta=setUpperCanvas("can_taumeff1prongVsPhiTheta");
  can_taumeff1prongVsPhiTheta->cd();
  can_taumeff1prongVsPhiTheta->SetRightMargin(0.18);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumeff1prongVsPhiTheta_ATJ=setUpperCanvas("can_taumeff1prongVsPhiTheta_ATJ");
  can_taumeff1prongVsPhiTheta_ATJ->cd();
  can_taumeff1prongVsPhiTheta_ATJ->SetRightMargin(0.18);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Draw("colz");

  TCanvas* can_taumeff1prongVsPhiTheta_CH=setUpperCanvas("can_taumeff1prongVsPhiTheta_CH");
  can_taumeff1prongVsPhiTheta_CH->cd();
  can_taumeff1prongVsPhiTheta_CH->SetRightMargin(0.18);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Draw("colz");

  TCanvas* can_taumeff1prongVsPhiTheta_ATJCH=setUpperCanvas("can_taumeff1prongVsPhiTheta_ATJCH");
  can_taumeff1prongVsPhiTheta_ATJCH->cd();
  can_taumeff1prongVsPhiTheta_ATJCH->SetRightMargin(0.18);
  //MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMinimum(0.30);
  //MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMaximum(0.90);
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Draw("colz");

  TCanvas* can_taumeff3prongVsPhiTheta=setUpperCanvas("can_taumeff3prongVsPhiTheta");
  can_taumeff3prongVsPhiTheta->cd();
  can_taumeff3prongVsPhiTheta->SetRightMargin(0.18);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumeff3prongVsPhiTheta_ATJ=setUpperCanvas("can_taumeff3prongVsPhiTheta_ATJ");
  can_taumeff3prongVsPhiTheta_ATJ->cd();
  can_taumeff3prongVsPhiTheta_ATJ->SetRightMargin(0.18);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJ->Draw("colz");

  TCanvas* can_taumeff3prongVsPhiTheta_CH=setUpperCanvas("can_taumeff3prongVsPhiTheta_CH");
  can_taumeff3prongVsPhiTheta_CH->cd();
  can_taumeff3prongVsPhiTheta_CH->SetRightMargin(0.18);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_CH->Draw("colz");

  TCanvas* can_taumeff3prongVsPhiTheta_ATJCH=setUpperCanvas("can_taumeff3prongVsPhiTheta_ATJCH");
  can_taumeff3prongVsPhiTheta_ATJCH->cd();
  can_taumeff3prongVsPhiTheta_ATJCH->SetRightMargin(0.18);
  //MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMinimum(0.20);
  //MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->SetMaximum(0.80);
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePhiTheta_ATJCH->Draw("colz");

  TCanvas* can_taumtrueeff3prongVsPhiTheta=setUpperCanvas("can_truetaumeff3prongVsPhiTheta");
  can_taumtrueeff3prongVsPhiTheta->cd();
  can_taumtrueeff3prongVsPhiTheta->SetRightMargin(0.18);
  TrueTaum3prongEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumtrueeff1prongVsPhiTheta=setUpperCanvas("can_truetaumeff1prongVsPhiTheta");
  can_taumtrueeff1prongVsPhiTheta->cd();
  can_taumtrueeff1prongVsPhiTheta->SetRightMargin(0.18);
  TrueTaum1prongEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumtrueeffLepElVsPhiTheta=setUpperCanvas("can_truetaumeffLepElVsPhiTheta");
  can_taumtrueeffLepElVsPhiTheta->cd();
  can_taumtrueeffLepElVsPhiTheta->SetRightMargin(0.18);
  TrueTaumLepElEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumtrueeffLepMuVsPhiTheta=setUpperCanvas("can_truetaumeffLepMuVsPhiTheta");
  can_taumtrueeffLepMuVsPhiTheta->cd();
  can_taumtrueeffLepMuVsPhiTheta->SetRightMargin(0.18);
  TrueTaumLepMuEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_tauptrueeff3prongVsPhiTheta=setUpperCanvas("can_truetaupeff3prongVsPhiTheta");
  can_tauptrueeff3prongVsPhiTheta->cd();
  can_tauptrueeff3prongVsPhiTheta->SetRightMargin(0.18);
  TrueTaup3prongEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_tauptrueeff1prongVsPhiTheta=setUpperCanvas("can_truetaupeff1prongVsPhiTheta");
  can_tauptrueeff1prongVsPhiTheta->cd();
  can_tauptrueeff1prongVsPhiTheta->SetRightMargin(0.18);
  TrueTaup1prongEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_tauptrueeffLepElVsPhiTheta=setUpperCanvas("can_truetaupeffLepElVsPhiTheta");
  can_tauptrueeffLepElVsPhiTheta->cd();
  can_tauptrueeffLepElVsPhiTheta->SetRightMargin(0.18);
  TrueTaupLepElEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_tauptrueeffLepMuVsPhiTheta=setUpperCanvas("can_truetaupeffLepMuVsPhiTheta");
  can_tauptrueeffLepMuVsPhiTheta->cd();
  can_tauptrueeffLepMuVsPhiTheta->SetRightMargin(0.18);
  TrueTaupLepMuEfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taum_nch_true_vs_nch_reco=setUpperCanvas("can_taum_nch_true_vs_nch_reco");
  can_taum_nch_true_vs_nch_reco->cd();
  can_taum_nch_true_vs_nch_reco->SetRightMargin(0.18);
  h_taum_nch_true_vs_nch_reco->Draw("colz");

  TCanvas* can_taum_npi_true_vs_npi_reco=setUpperCanvas("can_taum_npi_true_vs_npi_reco");
  can_taum_npi_true_vs_npi_reco->cd();
  can_taum_npi_true_vs_npi_reco->SetRightMargin(0.18);
  h_taum_npi_true_vs_npi_reco->Draw("colz");

  /*
  TCanvas* can_taupeffVsPhiTheta=setUpperCanvas("can_taupeffVsPhiTheta");
  can_taupeffVsPhiTheta->cd();
  can_taupeffVsPhiTheta->SetRightMargin(0.18);
  MatchedTaupCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");
 
  
  TCanvas* can_taumeffVsPhiTheta=setUpperCanvas("can_taumeffVsPhiTheta");
  can_taumeffVsPhiTheta->cd();
  can_taumeffVsPhiTheta->SetRightMargin(0.18);
  MatchedTaumCosTheta0_99EfficiencyVSTruePhiTheta->Draw("colz");

  TCanvas* can_taumeffVsTheta=setUpperCanvas("can_taumeffVsTheta");
  can_taumeffVsTheta->cd();
  MatchedTaumCosTheta0_99EfficiencyVSTrueTheta->Draw();

  TCanvas* can_taupeffVsTheta=setUpperCanvas("can_taupeffVsTheta");
  can_taupeffVsTheta->cd();
  MatchedTaupCosTheta0_99EfficiencyVSTrueTheta->Draw();

  TCanvas* can_taumeffVsPhi=setUpperCanvas("can_taumeffVsPhi");
  can_taumeffVsPhi->cd();
  MatchedTaumCosTheta0_99EfficiencyVSTruePhi->Draw();

  TCanvas* can_taupeffVsPhi=setUpperCanvas("can_taupeffVsPhi");
  can_taupeffVsPhi->cd();
  MatchedTaupCosTheta0_99EfficiencyVSTruePhi->Draw();
  
  //TEfficiency* matched2DHistTauM_PhiVsTheta = new TEfficiency(taum_recoPhiVsTheta,taum_truePhiVsTheta);
  //TEfficiency* matched2DHistTauP_PhiVsTheta = new TEfficiency(taup_recoPhiVsTheta,taup_truePhiVsTheta);

  //std::cout<<"test thing "<<taum_recoPhiVsTheta->GetBinContent(10)/taum_truePhiVsTheta->GetBinContent(10);
  TCanvas* can_taumeffVsPhiThetaHist=setUpperCanvas("can_taumeffVsPhiThetaHist");
  can_taumeffVsPhiThetaHist->cd();
  can_taumeffVsPhiThetaHist->SetRightMargin(0.18);
  taum_recoPhiVsTheta->Divide(taum_recoPhiVsTheta,taum_truePhiVsTheta,1.,1.,"B");
  taum_recoPhiVsTheta->Draw("colz");
  
  TCanvas* can_taupeffVsPhiThetaHist=setUpperCanvas("can_taupeffVsPhiThetaHist");
  can_taupeffVsPhiThetaHist->cd();
  can_taupeffVsPhiThetaHist->SetRightMargin(0.18);
  taup_recoPhiVsTheta->Divide(taup_recoPhiVsTheta,taup_truePhiVsTheta,1.,1.,"B");
  taup_recoPhiVsTheta->Draw("colz");
  */
  /*
  TCanvas* can_taum_nneut_true_1prong=setUpperCanvas("can_taum_nneut_true_1prong");
  can_taum_nneut_true_1prong->cd();
  TLegend* leg_taum_nneut_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_nneut_true_1prong->SetBorderSize(0);
  leg_taum_nneut_true_1prong->SetFillStyle(0);
  leg_taum_nneut_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_nneut_true_1prong->AddEntry(h_taum_nneut_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_nneut_true_1prong->AddEntry(h_taum_nneut_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_nneut_true_1prong->Draw();

  TCanvas* can_taum_nchother_true_1prong=setUpperCanvas("can_taum_nchother_true_1prong");
  can_taum_nchother_true_1prong->cd();
  TLegend* leg_taum_nchother_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_nchother_true_1prong->SetBorderSize(0);
  leg_taum_nchother_true_1prong->SetFillStyle(0);
  leg_taum_nchother_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_nchother_true_1prong->AddEntry(h_taum_nchother_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_nchother_true_1prong->AddEntry(h_taum_nchother_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_nchother_true_1prong->Draw();

  TCanvas* can_taum_nph_true_1prong=setUpperCanvas("can_taum_nph_true_1prong");
  can_taum_nph_true_1prong->cd();
  TLegend* leg_taum_nph_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_nph_true_1prong->SetBorderSize(0);
  leg_taum_nph_true_1prong->SetFillStyle(0);
  leg_taum_nph_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_nph_true_1prong->AddEntry(h_taum_nph_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_nph_true_1prong->AddEntry(h_taum_nph_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_nph_true_1prong->Draw();

  TCanvas* can_taum_npi0_true_1prong=setUpperCanvas("can_taum_npi0_true_1prong");
  can_taum_npi0_true_1prong->cd();
  TLegend* leg_taum_npi0_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_npi0_true_1prong->SetBorderSize(0);
  leg_taum_npi0_true_1prong->SetFillStyle(0);
  leg_taum_npi0_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_npi0_true_1prong->AddEntry(h_taum_npi0_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_npi0_true_1prong->AddEntry(h_taum_npi0_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_npi0_true_1prong->Draw();

  TCanvas* can_taum_dphi_min_pim_else_true_1prong=setUpperCanvas("can_taum_dphi_min_pim_else_true_1prong");
  can_taum_dphi_min_pim_else_true_1prong->cd();
  TLegend* leg_taum_dphi_min_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dphi_min_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dphi_min_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dphi_min_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dphi_min_pim_else_true_1prong->AddEntry(h_taum_dphi_min_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dphi_min_pim_else_true_1prong->AddEntry(h_taum_dphi_min_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dphi_min_pim_else_true_1prong->Draw();

  TCanvas* can_taum_dtheta_min_pim_else_true_1prong=setUpperCanvas("can_taum_dtheta_min_pim_else_true_1prong");
  can_taum_dtheta_min_pim_else_true_1prong->cd();
  TLegend* leg_taum_dtheta_min_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dtheta_min_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dtheta_min_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dtheta_min_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dtheta_min_pim_else_true_1prong->AddEntry(h_taum_dtheta_min_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dtheta_min_pim_else_true_1prong->AddEntry(h_taum_dtheta_min_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dtheta_min_pim_else_true_1prong->Draw();

  TCanvas* can_taum_dangle_min_pim_else_true_1prong=setUpperCanvas("can_taum_dangle_min_pim_else_true_1prong");
  can_taum_dangle_min_pim_else_true_1prong->cd();
  TLegend* leg_taum_dangle_min_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dangle_min_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dangle_min_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dangle_min_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dangle_min_pim_else_true_1prong->AddEntry(h_taum_dangle_min_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dangle_min_pim_else_true_1prong->AddEntry(h_taum_dangle_min_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dangle_min_pim_else_true_1prong->Draw();

  TCanvas* can_taum_dphi_max_pim_else_true_1prong=setUpperCanvas("can_taum_dphi_max_pim_else_true_1prong");
  can_taum_dphi_max_pim_else_true_1prong->cd();
  TLegend* leg_taum_dphi_max_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dphi_max_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dphi_max_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dphi_max_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dphi_max_pim_else_true_1prong->AddEntry(h_taum_dphi_max_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dphi_max_pim_else_true_1prong->AddEntry(h_taum_dphi_max_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dphi_max_pim_else_true_1prong->Draw();

  TCanvas* can_taum_dtheta_max_pim_else_true_1prong=setUpperCanvas("can_taum_dtheta_max_pim_else_true_1prong");
  can_taum_dtheta_max_pim_else_true_1prong->cd();
  TLegend* leg_taum_dtheta_max_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dtheta_max_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dtheta_max_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dtheta_max_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dtheta_max_pim_else_true_1prong->AddEntry(h_taum_dtheta_max_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dtheta_max_pim_else_true_1prong->AddEntry(h_taum_dtheta_max_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dtheta_max_pim_else_true_1prong->Draw();

  TCanvas* can_taum_dangle_max_pim_else_true_1prong=setUpperCanvas("can_taum_dangle_max_pim_else_true_1prong");
  can_taum_dangle_max_pim_else_true_1prong->cd();
  TLegend* leg_taum_dangle_max_pim_else_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_dangle_max_pim_else_true_1prong->SetBorderSize(0);
  leg_taum_dangle_max_pim_else_true_1prong->SetFillStyle(0);
  leg_taum_dangle_max_pim_else_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_dangle_max_pim_else_true_1prong->AddEntry(h_taum_dangle_max_pim_else_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_dangle_max_pim_else_true_1prong->AddEntry(h_taum_dangle_max_pim_else_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_dangle_max_pim_else_true_1prong->Draw();

  TCanvas* can_taum_jetmass_true_1prong=setUpperCanvas("can_taum_jetmass_true_1prong");
  can_taum_jetmass_true_1prong->cd();
  TLegend* leg_taum_jetmass_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_jetmass_true_1prong->SetBorderSize(0);
  leg_taum_jetmass_true_1prong->SetFillStyle(0);
  leg_taum_jetmass_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_jetmass_true_1prong->AddEntry(h_taum_jetmass_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_jetmass_true_1prong->AddEntry(h_taum_jetmass_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_jetmass_true_1prong->Draw();

  TCanvas* can_taum_jetenergy_true_1prong=setUpperCanvas("can_taum_jetenergy_true_1prong");
  can_taum_jetenergy_true_1prong->cd();
  TLegend* leg_taum_jetenergy_true_1prong=new TLegend(0.45,0.65,0.89,0.89);
  leg_taum_jetenergy_true_1prong->SetBorderSize(0);
  leg_taum_jetenergy_true_1prong->SetFillStyle(0);
  leg_taum_jetenergy_true_1prong->SetHeader("true #tau^{-} 1 prong","h");
  leg_taum_jetenergy_true_1prong->AddEntry(h_taum_jetenergy_true_1prong_no_reco ->DrawCopy("hist,e"),"no reco tau-jet");
  leg_taum_jetenergy_true_1prong->AddEntry(h_taum_jetenergy_true_1prong_and_reco ->DrawCopy("hist,e,same"),"with reco tau-jet");
  leg_taum_jetenergy_true_1prong->Draw();
*/
  TCanvas* can_taum_dangle_pim_true_any_pich_reco_1prong=setUpperCanvas("can_taum_dangle_pim_true_any_pich_reco_1prong");
  can_taum_dangle_pim_true_any_pich_reco_1prong->cd();
  h_taum_pim_true_1prong_no_reco_dangleMin_pim->Draw();

  TCanvas* can_tau_true_vs_reco_angles=setUpperCanvas("can_tau_true_vs_reco_angles");
  can_tau_true_vs_reco_angles->cd();

  TLegend* leg_tau_true_vs_reco_angle=new TLegend(0.45,0.65,0.89,0.89);
  leg_tau_true_vs_reco_angle->SetBorderSize(0);
  leg_tau_true_vs_reco_angle->SetFillStyle(0);
  leg_tau_true_vs_reco_angle->SetHeader("whizard #tau#tau, #sqrt{s}=200 GeV");
  leg_tau_true_vs_reco_angle->AddEntry(h_taum_true_1p_dPhi_reco->DrawCopy("hist,e"),"true #tau^{-}, 1prong");
  leg_tau_true_vs_reco_angle->AddEntry(h_taup_true_1p_dPhi_reco->DrawCopy("hist,e,same"),"true #tau^{+}, 1prong");
  leg_tau_true_vs_reco_angle->AddEntry(h_taum_true_3p_dPhi_reco->DrawCopy("hist,e,same"),"true #tau^{-}, 3prong");
  leg_tau_true_vs_reco_angle->AddEntry(h_taup_true_3p_dPhi_reco->DrawCopy("hist,e,same"),"true #tau^{+}, 3prong");
  leg_tau_true_vs_reco_angle->Draw();


  TCanvas* can_tau_true_pi1_energy=setUpperCanvas("can_tau_true_pi1_energy");
  can_tau_true_pi1_energy->cd();

  TLegend* leg_tau_true_pi1_energy=new TLegend(0.45,0.65,0.89,0.89);
  leg_tau_true_pi1_energy->SetBorderSize(0);
  leg_tau_true_pi1_energy->SetFillStyle(0);
  leg_tau_true_pi1_energy->SetHeader("whizard #tau#tau, #sqrt{s}=200 GeV");
  leg_tau_true_pi1_energy->AddEntry(h_taum1p_pi1_energy->DrawCopy("hist,e"),"true #tau^{-}, 1prong");
  leg_tau_true_pi1_energy->AddEntry(h_taup1p_pi1_energy->DrawCopy("hist,e,same"),"true #tau^{+}, 1prong");
  leg_tau_true_pi1_energy->AddEntry(h_taum3p_pi1_energy->DrawCopy("hist,e,same"),"true #tau^{-}, 3prong");
  leg_tau_true_pi1_energy->AddEntry(h_taup3p_pi1_energy->DrawCopy("hist,e,same"),"true #tau^{+}, 3prong");
  leg_tau_true_pi1_energy->Draw();

  TCanvas* can_tau_true_pi1_pt=setUpperCanvas("can_tau_true_pi1_pt");
  can_tau_true_pi1_pt->cd();

  TLegend* leg_tau_true_pi1_pt=new TLegend(0.45,0.65,0.89,0.89);
  leg_tau_true_pi1_pt->SetBorderSize(0);
  leg_tau_true_pi1_pt->SetFillStyle(0);
  leg_tau_true_pi1_pt->SetHeader("whizard #tau#tau, #sqrt{s}=200 GeV");
  leg_tau_true_pi1_pt->AddEntry(h_taum1p_pi1_pt->DrawCopy("hist,e"),"true #tau^{-}, 1prong");
  leg_tau_true_pi1_pt->AddEntry(h_taup1p_pi1_pt->DrawCopy("hist,e,same"),"true #tau^{+}, 1prong");
  leg_tau_true_pi1_pt->AddEntry(h_taum3p_pi1_pt->DrawCopy("hist,e,same"),"true #tau^{-}, 3prong");
  leg_tau_true_pi1_pt->AddEntry(h_taup3p_pi1_pt->DrawCopy("hist,e,same"),"true #tau^{+}, 3prong");
  leg_tau_true_pi1_pt->Draw();

  
  TCanvas* can_taup1p_eff_vs_pi1pt=setUpperCanvas("can_taup1p_eff_vs_pi1pt");
  can_taup1p_eff_vs_pi1pt->cd();
  MatchedTaup1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->Draw();
    
    gPad->Update();
    auto graph_tau1prong_ATJCH = MatchedTaup1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->GetPaintedGraph();
    graph_tau1prong_ATJCH->SetMinimum(0.6);
    graph_tau1prong_ATJCH->SetMaximum(1.0);
    gPad->Update();

  TCanvas* can_taup3p_eff_vs_pi1pt=setUpperCanvas("can_taup3p_eff_vs_pi1pt");
  can_taup3p_eff_vs_pi1pt->cd();
  MatchedTaup3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->Draw();

  TCanvas* can_taum1p_eff_vs_pi1pt=setUpperCanvas("can_taum1p_eff_vs_pi1pt");
  can_taum1p_eff_vs_pi1pt->cd();
  MatchedTaum1prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->Draw();


  TCanvas* can_taum3p_eff_vs_pi1pt=setUpperCanvas("can_taum3p_eff_vs_pi1pt");
  can_taum3p_eff_vs_pi1pt->cd();
  MatchedTaum3prongCosTheta0_99EfficiencyVSTruePi1Pt_ATJCH->Draw();


  return 1;

}

