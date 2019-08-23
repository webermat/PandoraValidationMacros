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


int plotPhRecoShower(){

  gROOT->ProcessLine("#include <vector>");


  int n_bins_phi=200;
  float limit_phi_bins_low= -3.2;
  float limit_phi_bins_high = 3.2;

  int n_bins2D_E=30;
  float limit_E_bins_low= 400;
  float limit_E_bins_high = 510;

  float limit_phi_low= -0.0015;//0.01 for 1,5 0.004 for 10,30, 0.0025 for 50,100,200,500, 0.0015 for 1000,1500
  float limit_phi_high = 0.0015;
  float limit_theta_low= -0.0052;//0.04 for 1,5 0.002 fo 10,30 0.0010 for 50,100,200,500, 0.00075 for 1000,1500 
  float limit_theta_high = 0.0052;

  float limit_hitE_bins_low= 0;
  float limit_hitE_bins_high = 2;

  float limit_hitE_bins_raw_low= 0;
  float limit_hitE_bins_raw_high = 0.2;

  int nbins_hits=80;
  float limit_hit_bins_low= 0.5;
  float limit_hit_bins_high = 2000;

  int nbins_hits_2D= 200;
  float limit_hit_2D_low= -1720;
  float limit_hit_2D_high = 1720;

  float limit_E_bins_high_HCAL = 200;
  float limit_E_bins_high_CAL = 1575;//80 bins
  float limit_E_bins_rel_low = 0.75; //0.9 to 1.1 for 1500,1000,500,200,100, 0.85 to 1.1 for 50,30 , 0.8 to 1.2 for 10,15, 0.7 to 1.3 for 5, 0 to 1.6 for 1
  float limit_E_bins_rel_high = 1.25;

  float limit_E_bis_leakage_low = 0.0;
  float limit_E_bins_leakage_high = 0.4;

  //phi0, theta 90 weights
  float weight_ecal=1;//0.991577;//0.987952 
  float weight_hcal=1;//0.881699;//1.00866; //weight should be H_CAL = p1 ECAL+ p0 -> third one is maybe really off?

  //float weight_ecal_B=1;
  //float weight_hcal_B=1;
  float weight_ecal_B=0.987952;//0.990771;//0.987952 
  float weight_hcal_B=1.00866;//0.938532;//1.00866
  float weight_ecal_E=0.989012;//0.989012;
  float weight_hcal_E=1.12815;//1.12815;

  //FINAL 1109 numbers
  //definition for detector note -> 80 bins on each side for the 2D plots, 20 for bin 0.50-0.75
  //1500 0-0.25, ECAL  1200-1550, HCAL 0-350 p0=1598.24+/-10.1898 p1=-1.05566+/-0.00715996 ---> w_e=0.990771, w_h=0.938532
  //1500 0.25-0.50 ECAL 1300-1550, HCAL 0-300, p0=1487.12+/-12.3071 p1=-0.979469+/-0.00846033 --->w_e=0.987966, w_h=1.00866
  //1500, 0.50-0.75, ECAL 1375-1550, HCAL 0-150, p0=1028.52+/-68.6626 p1=-0.671517+/0.0457402 ---> w_e=0.979345, w_h --> seems pretty unstable, maybe take numbers from above or just below, nothing for the barrel anyway
  //1500 0.87 to 0.98, ECAL 1000-1550, H 0-450, p0=1329.61+/-9.42788  p1=-0.876667+/-0.00648104 -->w_e=, w_h=1.12815

  //CosTheta07 weights
  //float weight_ecal=0.991;
  //float weight_hcal=0.950;//weight should be H_CAL = p1 ECAL+ p0

  //for 1500 sample
  //reweight HCAL: HCAL'= 1500 / p0 
  //reweihgt ECAL: ECAL' = -p1/p0 *1500

  string label_legend= "#gamma, E_{true}=1 GeV";
  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_1024FIX_mcparticlesph500_uniformCosTheta_160927_CLIC_o3_V06_newPhLikelihood_wPLug_1024FIX.root");
  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_1024FIX__mcparticlesph1_uniformCosTheta_161109_CLIC_o3_V06_newPhLikelihood161109_wPLug_1024FIX_ddmarlinpandora_minus_FlipMinCosSin.root");
  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_1024FIX__mcparticlesph1_uniformCosTheta_161109_CLIC_o3_V06_newPhLikelihood161109_wPLug_1024FIX_newDDMarlinPandoraSVN_full.root");
  TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_mcparticlesph10_uniformCosTheta_161215_CLIC_o3_V07_newPhLikelihood161122.root");
  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_1024FIX__mcparticlesph1_uniformCosTheta_161109_CLIC_o3_V06_newPhLikelihood161109_wPLug_1024FIX_ddmarlinpandora_testModuleMinus.root");
  //TFile* file_photons=TFile::Open("/afs/cern.ch/work/w/weberma2/ClicPerformance/showerStudyNEW_1024FIX__mcparticlesph1_uniformCosTheta_161109_CLIC_o3_V06_newPhLikelihood161109_wPLug_1024FIX_ddmarlinpandora_orig.root");
  
  TTree* tree_photons = (TTree*)file_photons->Get("showerData");
  // TTree* 

  float truePhE;
  float truePhPx;
  float truePhPy;
  float truePhPz;
  float truePhPhi;
  float truePhCosTheta;
  float truePhTheta;


  float totalECALPhE;
  float totalHCALPhE;
  float totalEndcapECALPhE;
  float totalEndcapHCALPhE;
  float totalPlugECALPhE;
  float totalPlugHCALPhE;

  unsigned int nHitsECALPh;
  unsigned int nHitsHCALPh;

  vector<float> *recoPhE=0;
  vector<float> *recoPhEHitCorrection=0;
  vector<float> *recoPhPx=0;
  vector<float> *recoPhPy=0;
  vector<float> *recoPhPz=0;
  vector<float> *recoPhPhi=0;
  vector<float> *recoPhCosTheta=0;
  vector<float> *recoPhTheta=0;
  vector<float> *recoPhPhi_logERW=0;
  vector<float> *recoPhCosTheta_logERW=0;
  vector<float> *recoPhTheta_logERW=0;
  vector<float> *recoPhPDGID=0;

  vector<float> *recoPhEEB=0;
  vector<float> *recoPhEEE=0;
  vector<float> *recoPhEEP=0;
  vector<float> *recoPhEHB=0;
  vector<float> *recoPhEHE=0;
  vector<float> *recoPhEHP=0;

  vector<int> *recoPhNhitsEB=0;
  vector<int> *recoPhNhitsEE=0;
  vector<int> *recoPhNhitsEP=0;
  vector<int> *recoPhNhitsHB=0;
  vector<int> *recoPhNhitsHE=0;
  vector<int> *recoPhNhitsHP=0;

  vector<float> *trackNHits=0;
  vector<int> *trackNdf=0;
  vector<int> *trackNTrackerBarrelHits=0;
  vector<int> *trackNTrackerEndcapHits=0;
  vector<int> *trackNTrackerVertexBarrelHits=0;
  vector<int> *trackNTrackerVertexEndcapHits=0;
  vector<float> *trackNExpectedTrackerHits=0;

  vector<float> *track_z_innermostHit=0;
  vector<float> *track_x_innermostHit=0;
  vector<float> *track_y_innermostHit=0;
  vector<float> *track_r_innermostHit=0;

  vector<float> *track_z_outermostZHit=0;
  vector<float> *track_x_outermostZHit=0;
  vector<float> *track_y_outermostZHit=0;
  vector<float> *track_r_outermostZHit=0;

  vector<float> *track_z_outermostRHit=0;
  vector<float> *track_x_outermostRHit=0;
  vector<float> *track_y_outermostRHit=0;
  vector<float> *track_r_outermostRHit=0;

  vector<float> *track_x_atCalo=0;
  vector<float> *track_y_atCalo=0;
  vector<float> *track_z_atCalo=0;

  vector<float> *trackd0=0;
  vector<float> *trackz0=0;
  vector<float> *trackChi2=0;
  vector<float> *trackzMin=0;
  vector<float> *trackPhiAtIP=0;
  vector<float> *trackThetaAtIP=0;
  vector<float> *trackPtAtIP=0;
  vector<float> *trackPAtIP=0;
  vector<float> *trackPhiAtCalo=0;
  vector<float> *trackThetaAtCalo=0;
  vector<float> *trackPtAtCalo=0;
  vector<float> *trackPAtCalo=0;
  vector<float> *trackClusterPosAngleMin=0;
  vector<float> *trackClusterEOverP=0;
  vector<float> *tracksigmaPOverP=0;
  vector<float> *track_cluster_minPFAdistance_atCalo=0;
  vector<float> *track_cluster_min_parPFAdistance_atCalo=0;
  vector<float> *track_cluster_DirAngle_minPFAdistance_atCalo=0;

  vector<float> *track_cluster_minPFAdistance_atCalo_tinyHad=0;
  vector<float> *track_cluster_min_parPFAdistance_atCalo_tinyHad=0;
  vector<float> *track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad=0;


  tree_photons->SetBranchAddress("track_nHits",&trackNHits);
  tree_photons->SetBranchAddress("track_ndf",&trackNdf);

  tree_photons->SetBranchAddress("track_nTrackerBarrelHits",&trackNTrackerBarrelHits);
  tree_photons->SetBranchAddress("track_nTrackerEndcapHits",&trackNTrackerEndcapHits);
  tree_photons->SetBranchAddress("track_nVertexBarrelHits",&trackNTrackerVertexBarrelHits);
  tree_photons->SetBranchAddress("track_nVertexEndcapHits",&trackNTrackerVertexEndcapHits);
  tree_photons->SetBranchAddress("track_nExpectedTrackerHits",&trackNExpectedTrackerHits);
  tree_photons->SetBranchAddress("track_cluster_minPFAdistance_atCalo", &track_cluster_minPFAdistance_atCalo);
  tree_photons->SetBranchAddress("track_cluster_min_parPFAdistance_atCalo", &track_cluster_min_parPFAdistance_atCalo);
  tree_photons->SetBranchAddress("track_cluster_DirAngle_minPFAdistance_atCalo", &track_cluster_DirAngle_minPFAdistance_atCalo);
  tree_photons->SetBranchAddress("track_cluster_minPFAdistance_atCalo_tinyHad", &track_cluster_minPFAdistance_atCalo_tinyHad);
  tree_photons->SetBranchAddress("track_cluster_min_parPFAdistance_atCalo_tinyHad", &track_cluster_min_parPFAdistance_atCalo_tinyHad);
  tree_photons->SetBranchAddress("track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad", &track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad);

  tree_photons->SetBranchAddress("track_chi2",&trackChi2);
  tree_photons->SetBranchAddress("track_d0",&trackd0);
  tree_photons->SetBranchAddress("track_z0",&trackz0);
  tree_photons->SetBranchAddress("track_zMin",&trackzMin);
  tree_photons->SetBranchAddress("track_sigmaPOverP", &tracksigmaPOverP);

  tree_photons->SetBranchAddress("track_x_atCalo",&track_x_atCalo);
  tree_photons->SetBranchAddress("track_y_atCalo",&track_y_atCalo);
  tree_photons->SetBranchAddress("track_z_atCalo",&track_z_atCalo);

  tree_photons->SetBranchAddress("track_x_innermostHit",&track_x_innermostHit);
  tree_photons->SetBranchAddress("track_y_innermostHit",&track_y_innermostHit);
  tree_photons->SetBranchAddress("track_z_innermostHit",&track_z_innermostHit);
  tree_photons->SetBranchAddress("track_r_innermostHit",&track_r_innermostHit);
  tree_photons->SetBranchAddress("track_x_outermostRHit",&track_x_outermostRHit);
  tree_photons->SetBranchAddress("track_y_outermostRHit",&track_y_outermostRHit);
  tree_photons->SetBranchAddress("track_z_outermostRHit",&track_z_outermostRHit);
  tree_photons->SetBranchAddress("track_r_outermostRHit",&track_r_outermostRHit);
  tree_photons->SetBranchAddress("track_x_outermostZHit",&track_x_outermostZHit);
  tree_photons->SetBranchAddress("track_y_outermostZHit",&track_y_outermostZHit);
  tree_photons->SetBranchAddress("track_z_outermostZHit",&track_z_outermostZHit);
  tree_photons->SetBranchAddress("track_r_outermostZHit",&track_r_outermostZHit);


  tree_photons->SetBranchAddress("track_Phi_atIP",&trackPhiAtIP);
  tree_photons->SetBranchAddress("track_Phi_atIP",&trackPhiAtIP);
  tree_photons->SetBranchAddress("track_Phi_atIP",&trackPhiAtIP);
  tree_photons->SetBranchAddress("track_Theta_atIP",&trackThetaAtIP);
  tree_photons->SetBranchAddress("track_pt_atIP",&trackPtAtIP);
  tree_photons->SetBranchAddress("track_p_atIP",&trackPAtIP);
  tree_photons->SetBranchAddress("track_Phi_atCalo",&trackPhiAtCalo);
  tree_photons->SetBranchAddress("track_Theta_atCalo",&trackThetaAtCalo);
  tree_photons->SetBranchAddress("track_pt_atCalo",&trackPtAtCalo);
  tree_photons->SetBranchAddress("track_p_atCalo",&trackPAtCalo);
  tree_photons->SetBranchAddress("track_cluster_PosAngle_min_atCalo",&trackClusterPosAngleMin);
  tree_photons->SetBranchAddress("track_minDist_cluster_EOverP",&trackClusterEOverP);


  vector<float>   *hitECALE=0;
  vector<float>   *hitECALEraw=0;
  vector<float>   *hitECALx=0;
  vector<float>   *hitECALy=0;
  vector<float>   *hitECALz=0;
  vector<int>     *hitECALLayer=0;

  vector<float>   *hitHCALE=0;
  vector<float>   *hitHCALEraw=0;
  vector<float>   *hitHCALx=0;
  vector<float>   *hitHCALy=0;
  vector<float>   *hitHCALz=0;


  vector<float>   *hitECALEE=0;
  vector<float>   *hitECALEEraw=0;
  vector<float>   *hitECALEx=0;
  vector<float>   *hitECALEy=0;
  vector<float>   *hitECALEz=0;
  vector<int>     *hitECALELayer=0;

  vector<float>   *hitHCALEE=0;
  vector<float>   *hitHCALEEraw=0;
  vector<float>   *hitHCALEx=0;
  vector<float>   *hitHCALEy=0;
  vector<float>   *hitHCALEz=0;

  tree_photons->SetBranchAddress("trueEnergy",&truePhE);
  tree_photons->SetBranchAddress("truePx",&truePhPx);
  tree_photons->SetBranchAddress("truePy",&truePhPy);
  tree_photons->SetBranchAddress("truePz",&truePhPz);
  tree_photons->SetBranchAddress("truePhi",&truePhPhi);
  tree_photons->SetBranchAddress("trueCosTheta",&truePhCosTheta);
  tree_photons->SetBranchAddress("trueTheta",&truePhTheta);
  tree_photons->SetBranchAddress("recoPh_Energy",&recoPhE);
  tree_photons->SetBranchAddress("recoPh_hitEnergyCellCorrection",&recoPhEHitCorrection);
  tree_photons->SetBranchAddress("recoPh_Px",&recoPhPx);
  tree_photons->SetBranchAddress("recoPh_Py",&recoPhPy);
  tree_photons->SetBranchAddress("recoPh_Pz",&recoPhPz);
  tree_photons->SetBranchAddress("recoPh_Phi",&recoPhPhi);
  tree_photons->SetBranchAddress("recoPh_CosTheta",&recoPhCosTheta);
  tree_photons->SetBranchAddress("recoPh_Theta",&recoPhTheta);
  tree_photons->SetBranchAddress("recoPh_Phi_logERW",&recoPhPhi_logERW);
  tree_photons->SetBranchAddress("recoPh_CosTheta_logERW",&recoPhCosTheta_logERW);
  tree_photons->SetBranchAddress("recoPh_Theta_logERW",&recoPhTheta_logERW);
  tree_photons->SetBranchAddress("recoPh_PDGID",&recoPhPDGID);
  tree_photons->SetBranchAddress("recoPh_E_EB",&recoPhEEB);
  tree_photons->SetBranchAddress("recoPh_E_EE",&recoPhEEE);
  tree_photons->SetBranchAddress("recoPh_E_EP",&recoPhEEP);
  tree_photons->SetBranchAddress("recoPh_E_HB",&recoPhEHB);
  tree_photons->SetBranchAddress("recoPh_E_HE",&recoPhEHE);
  tree_photons->SetBranchAddress("recoPh_E_HP",&recoPhEHP);
  tree_photons->SetBranchAddress("recoPh_nhitsEB",&recoPhNhitsEB);
  tree_photons->SetBranchAddress("recoPh_nhitsEE",&recoPhNhitsEE);
  tree_photons->SetBranchAddress("recoPh_nhitsEP",&recoPhNhitsEP);
  tree_photons->SetBranchAddress("recoPh_nhitsHB",&recoPhNhitsHB);
  tree_photons->SetBranchAddress("recoPh_nhitsHE",&recoPhNhitsHE);
  tree_photons->SetBranchAddress("recoPh_nhitsHP",&recoPhNhitsHP);

  tree_photons->SetBranchAddress("totalEnergy",&totalECALPhE);
  tree_photons->SetBranchAddress("totalLeakEnergy",&totalHCALPhE);
  tree_photons->SetBranchAddress("totalEndcapEnergy",&totalEndcapECALPhE);
  tree_photons->SetBranchAddress("totalEndcapLeakEnergy",&totalEndcapHCALPhE);
  tree_photons->SetBranchAddress("totalEndcapEnergy",&totalEndcapECALPhE);
  tree_photons->SetBranchAddress("totalEndcapLeakEnergy",&totalEndcapHCALPhE);
  tree_photons->SetBranchAddress("totalPlugEnergy",&totalPlugECALPhE);
  tree_photons->SetBranchAddress("totalPlugLeakEnergy",&totalPlugHCALPhE);
  tree_photons->SetBranchAddress("hit_n",&nHitsECALPh);
  tree_photons->SetBranchAddress("hit_leak_n",&nHitsHCALPh);

  tree_photons->SetBranchAddress("hit_E",&hitECALE);
  tree_photons->SetBranchAddress("hit_rawE",&hitECALEraw);
  tree_photons->SetBranchAddress("hit_x",&hitECALx);
  tree_photons->SetBranchAddress("hit_y",&hitECALy);
  tree_photons->SetBranchAddress("hit_z",&hitECALz);
  tree_photons->SetBranchAddress("hit_layer",&hitECALLayer);

  tree_photons->SetBranchAddress("hit_leak_E",&hitHCALE);
  tree_photons->SetBranchAddress("hit_leak_rawE",&hitHCALEraw);
  tree_photons->SetBranchAddress("hit_leak_x",&hitHCALx);
  tree_photons->SetBranchAddress("hit_leak_y",&hitHCALy);
  tree_photons->SetBranchAddress("hit_leak_z",&hitHCALz);

  tree_photons->SetBranchAddress("Endcap_hit_E",&hitECALEE);
  tree_photons->SetBranchAddress("Endcap_hit_rawE",&hitECALEEraw);
  tree_photons->SetBranchAddress("Endcap_hit_x",&hitECALEx);
  tree_photons->SetBranchAddress("Endcap_hit_y",&hitECALEy);
  tree_photons->SetBranchAddress("Endcap_hit_z",&hitECALEz);
  tree_photons->SetBranchAddress("Endcap_hit_layer",&hitECALELayer);

  tree_photons->SetBranchAddress("Endcap_hit_leak_E",&hitHCALEE);
  tree_photons->SetBranchAddress("Endcap_hit_leak_rawE",&hitHCALEEraw);
  tree_photons->SetBranchAddress("Endcap_hit_leak_x",&hitHCALEx);
  tree_photons->SetBranchAddress("Endcap_hit_leak_y",&hitHCALEy);
  tree_photons->SetBranchAddress("Endcap_hit_leak_z",&hitHCALEz);

  TH1F* h_layerEB_E =  new TH1F("h_layerEB_E","", 40, -0.5, 39.5);
  h_layerEB_E  ->Sumw2();
  TH1F* h_layerEE_E =  new TH1F("h_layerEE_E","", 40, -0.5, 39.5);
  h_layerEE_E  ->Sumw2();

  TH1F* h_layerEB_E_ph =  new TH1F("h_layerEB_E_ph","", 40, -0.5, 39.5);
  h_layerEB_E_ph  ->Sumw2();
  h_layerEB_E_ph->SetLineColor(kBlue);
  TH1F* h_layerEE_E_ph =  new TH1F("h_layerEE_E_ph","", 40, -0.5, 39.5);
  h_layerEE_E_ph  ->Sumw2();
  h_layerEE_E_ph->SetLineColor(kBlue);

  TH1F* h_layerEB_E_n =  new TH1F("h_layerEB_E_n","", 40, -0.5, 39.5);
  h_layerEB_E_n  ->Sumw2();
  h_layerEB_E_n->SetLineColor(kRed);
  TH1F* h_layerEE_E_n =  new TH1F("h_layerEE_E_n","", 40, -0.5, 39.5);
  h_layerEE_E_n  ->Sumw2();
  h_layerEE_E_n->SetLineColor(kRed);

  TH1F* h_layerHB_E =  new TH1F("h_layerHB_E","", 60, -0.5, 59.5);
  h_layerHB_E  ->Sumw2();
  TH1F* h_layerHE_E =  new TH1F("h_layerHE_E","", 60, -0.5, 59.5);
  h_layerHE_E  ->Sumw2();

  TH1F* h_layerEB =  new TH1F("h_nHit_perEBLayer","", 40, -0.5, 39.5);
  h_layerEB->Sumw2();
  h_layerEB->SetLineWidth(2);
  h_layerEB->SetLineColor(2);
  TH1F* h_layerEE =  new TH1F("h_nHit_perEELayer","", 40, -0.5, 39.5);
  h_layerEE->Sumw2();
  h_layerEE->SetLineWidth(2);
  h_layerEE->SetLineColor(4);

  TH1F* h_phi_phReco =  new TH1F("h_phi_phReco","", n_bins_phi/6., limit_phi_bins_low, limit_phi_bins_high);
  //h_phi_phReco ->Sumw2();
  //h_phi_phReco ->SetLineWidth(2);
  //h_phi_phReco ->SetMarkerStyle(20);

  TH1F* h_phi_nReco =  new TH1F("h_phi_nReco","", n_bins_phi/6., limit_phi_bins_low, limit_phi_bins_high);
  //h_phi_nReco ->Sumw2();
  //h_phi_nReco ->SetLineWidth(2);
  //h_phi_nReco ->SetLineColor(2);
  //h_phi_nReco ->SetMarkerStyle(20);

  TH2F* h_theta_phReco_vs_theta_true =  new TH2F("h_theta_phReco_vs_theta_true","", n_bins_phi/6.,0, TMath::Pi() ,n_bins_phi/6., 0, TMath::Pi());
  h_theta_phReco_vs_theta_true ->Sumw2();
  h_theta_phReco_vs_theta_true->GetYaxis()->SetTitle("theta_{true}");
  h_theta_phReco_vs_theta_true->GetXaxis()->SetTitle("phi_{PFO}");

  TH2F* h_dtheta_vs_theta_true =  new TH2F("h_dtheta_phReco_vs_theta_true","", n_bins_phi/6., limit_theta_low, limit_theta_high,n_bins_phi/6., 0, TMath::Pi());
  h_dtheta_vs_theta_true ->Sumw2();
  h_dtheta_vs_theta_true->GetYaxis()->SetTitle("theta_{true}");
  h_dtheta_vs_theta_true->GetXaxis()->SetTitle("theta^{PFA}-theta_{true}");

  TH2F* h_dtheta_vs_costheta_true =  new TH2F("h_dtheta_phReco_vs_costheta_true","", n_bins_phi/6., limit_theta_low, limit_theta_high,n_bins_phi/6., -1.0, 1.0);
  h_dtheta_vs_costheta_true ->Sumw2();
  h_dtheta_vs_costheta_true->GetYaxis()->SetTitle("cos#theta_{true}");
  h_dtheta_vs_costheta_true->GetXaxis()->SetTitle("theta^{PFA}-theta_{true}");

  TH2F* h_dtheta_vs_phi_true =  new TH2F("h_dtheta_phReco_vs_phi_true","", n_bins_phi/6., limit_theta_low, limit_theta_high,n_bins_phi/6., -TMath::Pi(), TMath::Pi());
  h_dtheta_vs_phi_true ->Sumw2();
  h_dtheta_vs_phi_true->GetYaxis()->SetTitle("phi_{true}");
  h_dtheta_vs_phi_true->GetXaxis()->SetTitle("theta^{PFA}-theta_{true}");

  TH2F* h_dphi_vs_phi_true =  new TH2F("h_dphi_vs_phi_true","", n_bins_phi/6., limit_phi_low, limit_phi_high,n_bins_phi/6., -TMath::Pi(), TMath::Pi());
  h_dphi_vs_phi_true ->Sumw2();
  h_dphi_vs_phi_true->GetYaxis()->SetTitle("phi_{true}");
  h_dphi_vs_phi_true->GetXaxis()->SetTitle("phi_{PFO}-phi_{true}");

  TH2F* h_dphi_vs_theta_true =  new TH2F("h_dphi_vs_theta_true","", n_bins_phi/6., limit_phi_low, limit_phi_high,n_bins_phi/6., 0, TMath::Pi());
  h_dphi_vs_theta_true ->Sumw2();
  h_dphi_vs_theta_true->GetYaxis()->SetTitle("theta_{true}");
  h_dphi_vs_theta_true->GetXaxis()->SetTitle("phi_{PFO}-phi_{true}");

  TH2F* h_dphi_vs_costheta_true =  new TH2F("h_dphi_vs_costheta_true","", n_bins_phi/6., limit_phi_low, limit_phi_high,n_bins_phi/6., -1.0,1.0);
  h_dphi_vs_costheta_true ->Sumw2();
  h_dphi_vs_costheta_true->GetYaxis()->SetTitle("cos#theta_{true}");
  h_dphi_vs_costheta_true->GetXaxis()->SetTitle("phi_{PFO}-phi_{true}");

  TH2F* h_phi_phReco_vs_phi_true =  new TH2F("h_phi_phReco_vs_phi_true","", n_bins_phi/6.,-TMath::Pi(), TMath::Pi() ,n_bins_phi/6., -TMath::Pi(), TMath::Pi());
  h_phi_phReco_vs_phi_true ->Sumw2();
  h_phi_phReco_vs_phi_true->GetYaxis()->SetTitle("phi_{true}");
  h_phi_phReco_vs_phi_true->GetXaxis()->SetTitle("phi_{PFO}");

  TH1F* h_dtheta_phReco =  new TH1F("h_dtheta_phReco","", n_bins_phi/4., limit_theta_low, limit_theta_high);
  h_dtheta_phReco ->Sumw2();
  h_dtheta_phReco ->SetLineWidth(2);
  h_dtheta_phReco ->SetLineColor(2);
  h_dtheta_phReco ->SetMarkerStyle(20);

  TH1F* h_dphi_phReco =  new TH1F("h_dphi_phReco","", n_bins_phi/4., limit_phi_low, limit_phi_high);
  h_dphi_phReco ->Sumw2();
  h_dphi_phReco ->SetLineWidth(2);
  h_dphi_phReco ->SetLineColor(2);
  h_dphi_phReco ->SetMarkerStyle(20);

 TH1F* h_dtheta_logERW_phReco =  new TH1F("h_dtheta_logERW_phReco","", n_bins_phi/4., limit_theta_low, limit_theta_high);
  h_dtheta_logERW_phReco ->Sumw2();
  h_dtheta_logERW_phReco ->SetLineWidth(2);
  h_dtheta_logERW_phReco ->SetLineColor(2);
  h_dtheta_logERW_phReco ->SetMarkerStyle(20);

  TH1F* h_dphi_logERW_phReco =  new TH1F("h_dphi_logERW_phReco","", n_bins_phi/4., limit_phi_low, limit_phi_high);
  h_dphi_logERW_phReco ->Sumw2();
  h_dphi_logERW_phReco ->SetLineWidth(2);
  h_dphi_logERW_phReco ->SetLineColor(2);
  h_dphi_logERW_phReco ->SetMarkerStyle(20);

  TH1F* h_theta_phReco =  new TH1F("h_theta_phReco","", n_bins_phi/4., 0, limit_phi_bins_high);
  h_theta_phReco ->Sumw2();
  h_theta_phReco ->SetLineWidth(2);
  h_theta_phReco ->SetMarkerStyle(20);

  TH1F* h_theta_nReco =  new TH1F("h_theta_nReco","", n_bins_phi/4., 0, limit_phi_bins_high);
  h_theta_nReco ->Sumw2();
  h_theta_nReco ->SetLineWidth(2);
  h_theta_nReco ->SetLineColor(2);
  h_theta_nReco ->SetMarkerStyle(20);

  TH1F* h_CosTheta_phReco =  new TH1F("h_CosTheta_phReco","", n_bins_phi/6., -1,1);
  h_CosTheta_phReco ->Sumw2();
  h_CosTheta_phReco ->SetLineWidth(2);
  h_CosTheta_phReco ->SetMarkerStyle(20);

  TH1F* h_CosTheta_nReco =  new TH1F("h_CosTheta_nReco","", n_bins_phi/6., -1, 1);
  h_CosTheta_nReco ->Sumw2();
  h_CosTheta_nReco ->SetLineWidth(2);
  h_CosTheta_nReco ->SetLineColor(2);
  h_CosTheta_nReco ->SetMarkerStyle(20);

  int nbins_E=40;

  TH1F* h_E_phReco =  new TH1F("h_E_phReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_phReco ->Sumw2();
  h_E_phReco ->SetLineWidth(2);
  //h_E_phReco ->SetMarkerStyle(20);

  TH1F* h_E_hitsum =  new TH1F("h_E_hitsum","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_E_hitsum ->Sumw2();
  h_E_hitsum ->SetLineWidth(2);
  //h_E_phReco ->SetMarkerStyle(20);

  TH1F* h_EE_hitsum =  new TH1F("h_EE_hitsum","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_EE_hitsum ->Sumw2();
  h_EE_hitsum ->SetLineWidth(2);
  h_EE_hitsum ->SetLineColor(2);
  //h_E_phReco ->SetMarkerStyle(20);

  TH1F* h_EB_hitsum =  new TH1F("h_EB_hitsum","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_EB_hitsum ->Sumw2();
  h_EB_hitsum ->SetLineWidth(2);
  h_EB_hitsum ->SetLineColor(4);
  //h_E_phReco ->SetMarkerStyle(20);

  TH1F* h_E_hitsumH =  new TH1F("h_E_hitsumH","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_E_hitsumH ->Sumw2();
  h_E_hitsumH ->SetLineWidth(2);

  TH1F* h_HE_hitsum =  new TH1F("h_HE_hitsum","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_HE_hitsum ->Sumw2();
  h_HE_hitsum ->SetLineWidth(2);
  h_HE_hitsum ->SetLineColor(2);

  TH1F* h_HB_hitsum =  new TH1F("h_HB_hitsum","", nbins_E, limit_E_bins_low/(36.+5), limit_E_bins_high/36.);
  h_HB_hitsum ->Sumw2();
  h_HB_hitsum ->SetLineWidth(2);
  h_HB_hitsum ->SetLineColor(4);


  TH1F* h_E_HB_Reco =  new TH1F("h_E_HB_Reco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HB_Reco ->Sumw2();
  h_E_HB_Reco ->SetLineColor(2);
  h_E_HB_Reco ->SetLineWidth(2);

  TH1F* h_E_HE_Reco =  new TH1F("h_E_HE_Reco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HE_Reco ->Sumw2();
  h_E_HE_Reco ->SetLineColor(4);
  h_E_HE_Reco ->SetLineWidth(2);


  TH1F* h_E_nReco =  new TH1F("h_E_nReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_nReco ->Sumw2();
  h_E_nReco ->SetLineColor(2);
  h_E_nReco ->SetLineWidth(2);
  //h_E_nReco ->SetMarkerStyle(20);

  TH1F* h_E_EB_phReco =  new TH1F("h_E_EB_phReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_EB_phReco ->Sumw2();
  h_E_EB_phReco ->SetLineWidth(2);
  //h_E_EB_phReco ->SetMarkerStyle(20);

  TH1F* h_E_EB_nReco =  new TH1F("h_E_EB_nReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_EB_nReco ->Sumw2();
  h_E_EB_nReco ->SetLineColor(2);
  h_E_EB_nReco ->SetLineWidth(2);
  //h_E_EB_nReco ->SetMarkerStyle(20);

  TH1F* h_E_EE_phReco =  new TH1F("h_E_EE_phReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_EE_phReco ->Sumw2();
  h_E_EE_phReco ->SetLineWidth(2);
  //h_E_EE_phReco ->SetMarkerStyle(20);

  TH1F* h_E_EE_nReco =  new TH1F("h_E_EE_nReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_EE_nReco ->Sumw2();
  h_E_EE_nReco ->SetLineColor(2);
  h_E_EE_nReco ->SetLineWidth(2);
  //h_E_EE_nReco ->SetMarkerStyle(20);

 TH1F* h_E_HB_phReco =  new TH1F("h_E_HB_phReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HB_phReco ->Sumw2();
  h_E_HB_phReco ->SetLineWidth(2);
  //h_E_HB_phReco ->SetMarkerStyle(20);

  TH1F* h_E_HB_nReco =  new TH1F("h_E_HB_nReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HB_nReco ->Sumw2();
  h_E_HB_nReco ->SetLineColor(2);
  h_E_HB_nReco ->SetLineWidth(2);
  //h_E_HB_nReco ->SetMarkerStyle(20);

  TH1F* h_E_HE_phReco =  new TH1F("h_E_HE_phReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HE_phReco ->Sumw2();
  h_E_HE_phReco ->SetLineWidth(2);
  //h_E_HE_phReco ->SetMarkerStyle(20);

  TH1F* h_E_HE_nReco =  new TH1F("h_E_HE_nReco","", nbins_E, limit_E_bins_low, limit_E_bins_high);
  h_E_HE_nReco ->Sumw2();
  h_E_HE_nReco ->SetLineColor(2);
  h_E_HE_nReco ->SetLineWidth(2);
  //h_E_HE_nReco ->SetMarkerStyle(20);

  TH1F* h_nhits_phReco =  new TH1F("h_nhits_phReco","", nbins_hits/2., limit_hit_bins_low, limit_hit_bins_high);
  h_nhits_phReco ->Sumw2();
  h_nhits_phReco ->SetLineWidth(2);

  TH1F* h_nhits_nReco =  new TH1F("h_nhits_nReco","", nbins_hits/2., limit_hit_bins_low, limit_hit_bins_high);
  h_nhits_nReco ->Sumw2();
  h_nhits_nReco ->SetLineColor(2);
  h_nhits_nReco ->SetLineWidth(2);

  TH2F* h_2D_hit_x_y_phReco =  new TH2F("h_2D_hit_x_y_phReco","",nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high);
  h_2D_hit_x_y_phReco ->Sumw2();
  //h_2D_hit_x_y_phReco ->SetLineWidth(2);
  //h_2D_hit_x_y_phReco ->SetMarkerStyle(20);
  h_2D_hit_x_y_phReco ->GetXaxis()->SetTitle("hit_x");
  h_2D_hit_x_y_phReco ->GetYaxis()->SetTitle("hit_y");
  //h_2D_hit_x_y_phReco ->SetMaximum(200);

  TH2F* h_2D_hit_x_y_nReco =  new TH2F("h_2D_hit_x_y_nReco","",nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high);
  h_2D_hit_x_y_nReco ->Sumw2();
  //h_2D_hit_x_y_nReco ->SetLineWidth(2);
  //h_2D_hit_x_y_nReco ->SetMarkerStyle(20);
  h_2D_hit_x_y_nReco ->GetXaxis()->SetTitle("hit_x");
  h_2D_hit_x_y_nReco ->GetYaxis()->SetTitle("hit_y");
  //h_2D_hit_x_y_nReco ->SetMaximum(200);

  TH2F* h_2D_hit_x_y =  new TH2F("h_2D_hit_x_y","",nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high);
  h_2D_hit_x_y ->Sumw2();
  //h_2D_hit_x_y ->SetLineWidth(2);
  //h_2D_hit_x_y ->SetMarkerStyle(20);
  h_2D_hit_x_y ->GetXaxis()->SetTitle("hit_x");
  h_2D_hit_x_y ->GetYaxis()->SetTitle("hit_y");
  //h_2D_hit_x_y ->SetMaximum(200);

  //TH3F* h_3D_hit_x_y_z_barrel =  new TH3F("h_3D_hit_x_y_z_barrel","",nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D, limit_hit_2D_low, limit_hit_2D_high,350, 190,890);
  TH3F* h_3D_hit_x_y_z_barrel =  new TH3F("h_3D_hit_x_y_z_barrel","",nbins_hits_2D/2., -400, 400,nbins_hits_2D/2., 1500, 1700,350, 300,1000);
  h_3D_hit_x_y_z_barrel ->Sumw2();
  //h_2D_hit_x_y_barrel ->SetLineWidth(2);
  //h_2D_hit_x_y_barrel ->SetMarkerStyle(20);
  h_3D_hit_x_y_z_barrel ->GetXaxis()->SetTitle("hit_x");
  h_3D_hit_x_y_z_barrel ->GetYaxis()->SetTitle("hit_y");
  h_3D_hit_x_y_z_barrel ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_barrel ->SetMaximum(200);

  TH2F* h_2D_hit_z_r_barrel =  new TH2F("h_2D_hit_z_r_barrel","",200, 300,1000,200, 1500, 1702);
  h_2D_hit_z_r_barrel ->Sumw2();
  //h_2D_hit_x_y_barrel ->SetLineWidth(2);
  //h_2D_hit_x_y_barrel ->SetMarkerStyle(20);
  h_2D_hit_z_r_barrel ->GetXaxis()->SetTitle("hit_r");
  h_2D_hit_z_r_barrel ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_barrel ->SetMaximum(200);

 TH2F* h_2D_hitE_z_r_barrel =  new TH2F("h_2D_hitE_z_r_barrel","",200, 300,1000,200, 1500, 1702);
  h_2D_hitE_z_r_barrel ->Sumw2();
  //h_2D_hit_x_y_barrel ->SetLineWidth(2);
  //h_2D_hit_x_y_barrel ->SetMarkerStyle(20);
  h_2D_hitE_z_r_barrel ->GetXaxis()->SetTitle("hit_r");
  h_2D_hitE_z_r_barrel ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_barrel ->SetMaximum(200);

  TH2F* h_2D_hit_x_y_endcap =  new TH2F("h_2D_hit_x_y_endcap","",nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D , limit_hit_2D_low, limit_hit_2D_high);
  h_2D_hit_x_y_endcap ->Sumw2();
  //h_2D_hit_x_y_endcap ->SetLineWidth(2);
  //h_2D_hit_x_y_endcap ->SetMarkerStyle(20);
  h_2D_hit_x_y_endcap ->GetXaxis()->SetTitle("hit_x");
  h_2D_hit_x_y_endcap ->GetYaxis()->SetTitle("hit_y");
  //h_2D_hit_x_y_endcap ->SetMaximum(200);

  TH3F* h_3D_hit_x_y_z_endcap =  new TH3F("h_3D_hit_x_y_z_endcap","",nbins_hits_2D/4. , limit_hit_2D_low, limit_hit_2D_high,nbins_hits_2D/4., limit_hit_2D_low, limit_hit_2D_high, 230, 2290,2520);
  h_3D_hit_x_y_z_endcap ->Sumw2();
  //h_2D_hit_x_y_endcap ->SetLineWidth(2);
  //h_2D_hit_x_y_endcap ->SetMarkerStyle(20);
  h_3D_hit_x_y_z_endcap ->GetXaxis()->SetTitle("hit_x");
  h_3D_hit_x_y_z_endcap ->GetYaxis()->SetTitle("hit_y");

  h_3D_hit_x_y_z_endcap ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_endcap ->SetMaximum(200);

  TH2F* h_2D_hit_r_z_endcap =  new TH2F("h_2D_hit_r_z_endcap","",200 , 410, 1700,200, 2307,2509);
  h_2D_hit_r_z_endcap ->Sumw2();
  //h_2D_hit_r_y_endcap ->SetLineWidth(2);
  //h_2D_hit_r_y_endcap ->SetMarkerStyle(20);
  h_2D_hit_r_z_endcap ->GetXaxis()->SetTitle("hit_r");
  h_2D_hit_r_z_endcap ->GetYaxis()->SetTitle("hit_z");
  //h_3D_hit_x_z_endcap ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_endcap ->SetMaximum(200);

 TH2F* h_2D_hitE_r_z_endcap =  new TH2F("h_2D_hitE_r_z_endcap","",200 , 410, 1700,200, 2307,2509);
  h_2D_hitE_r_z_endcap ->Sumw2();
  //h_2D_hit_r_y_endcap ->SetLineWidth(2);
  //h_2D_hit_r_y_endcap ->SetMarkerStyle(20);
  h_2D_hitE_r_z_endcap ->GetXaxis()->SetTitle("hit_r");
  h_2D_hitE_r_z_endcap ->GetYaxis()->SetTitle("hit_z");
  //h_3D_hit_x_z_endcap ->GetYaxis()->SetTitle("hit_z");
  //h_2D_hit_x_y_endcap ->SetMaximum(200);


  TH1F* h_hit_E_phReco =  new TH1F("h_hit_E_phReco" ,"", n_bins_phi, limit_hitE_bins_low, limit_hitE_bins_high);
  h_hit_E_phReco  ->Sumw2();
  h_hit_E_phReco  ->SetLineWidth(2);
  //h_hit_E_phReco  ->SetMarkerStyle(20);

  TH1F* h_hit_E_nReco =  new TH1F("h_hit_E_nReco" ,"", n_bins_phi, limit_hitE_bins_low, limit_hitE_bins_high);
  h_hit_E_nReco  ->Sumw2();
  h_hit_E_nReco ->SetLineColor(2);
  h_hit_E_nReco  ->SetLineWidth(2);
  //h_hit_E_nReco  ->SetMarkerStyle(20);

  TH1F* h_hit_Eraw_truePhEE =  new TH1F("h_hit_Eraw_truePhEE" ,"", n_bins_phi, limit_hitE_bins_raw_low, limit_hitE_bins_raw_high);
  h_hit_Eraw_truePhEE ->Sumw2();
  h_hit_Eraw_truePhEE  ->SetLineWidth(2);
  h_hit_Eraw_truePhEE  ->SetLineColor(2);
  //h_hit_Eraw_phReco  ->SetMarkerStyle(20);

  TH1F* h_hit_Eraw_truePhEB =  new TH1F("h_hit_Eraw_truePhEB" ,"", n_bins_phi, limit_hitE_bins_raw_low, limit_hitE_bins_raw_high);
  h_hit_Eraw_truePhEB ->Sumw2();
  h_hit_Eraw_truePhEB  ->SetLineWidth(2);
  h_hit_Eraw_truePhEB  ->SetLineColor(2);
  //h_hit_Eraw_phReco  ->SetMarkerStyle(20);

  TH1F* h_hit_Eraw_phReco =  new TH1F("h_hit_Eraw_phReco" ,"", n_bins_phi, limit_hitE_bins_raw_low, limit_hitE_bins_raw_high);
  h_hit_Eraw_phReco  ->Sumw2();
  h_hit_Eraw_phReco  ->SetLineWidth(2);
  //h_hit_Eraw_phReco  ->SetMarkerStyle(20);

  TH1F* h_hit_Eraw_nReco =  new TH1F("h_hit_Eraw_nReco" ,"", n_bins_phi, limit_hitE_bins_raw_low, limit_hitE_bins_raw_high);
  h_hit_Eraw_nReco  ->Sumw2();
  h_hit_Eraw_nReco ->SetLineColor(2);
  h_hit_Eraw_nReco  ->SetLineWidth(2);
  //h_hit_Eraw_nReco  ->SetMarkerStyle(20);

  TH1F* h_hit_Eraw_totSum_over_trueE  =  new TH1F("h_hit_Eraw_totSum_over_trueE" ,"", 200, 0.020, 0.035);
  h_hit_Eraw_totSum_over_trueE   ->Sumw2();
  h_hit_Eraw_totSum_over_trueE  ->SetLineWidth(2);

  TH1F* h_phiReco_phiTrue0 =  new TH1F("h_phiReco_phiTrue0 ","", n_bins_phi, limit_phi_bins_low, limit_phi_bins_high);
  h_phiReco_phiTrue0 ->Sumw2();
  h_phiReco_phiTrue0 ->SetLineWidth(2);
  h_phiReco_phiTrue0 ->SetMarkerStyle(20);

  TH1F* h_phiTrueOnly =  new TH1F("h_phiTrueOnly ","", n_bins_phi, limit_phi_bins_low, limit_phi_bins_high);
  h_phiTrueOnly ->Sumw2();
  h_phiTrueOnly ->SetLineWidth(2);
  h_phiTrueOnly ->SetMarkerStyle(20);

  TH1F* h_ThetaTrueOnly =  new TH1F("h_ThetaTrueOnly ","", n_bins_phi, 0, limit_phi_bins_high);
  h_ThetaTrueOnly ->Sumw2();
  h_ThetaTrueOnly ->SetLineWidth(2);
  h_ThetaTrueOnly ->SetMarkerStyle(20);


  TH1F* h_Delta_phiReco_phiTrue =  new TH1F("h_Delta_phiReco_phiTrue ","", n_bins_phi, limit_phi_bins_low/10., limit_phi_bins_high/10.);
  h_Delta_phiReco_phiTrue ->Sumw2();
  h_Delta_phiReco_phiTrue ->SetLineWidth(2);
  h_Delta_phiReco_phiTrue ->SetMarkerStyle(20);

  TH2F* h_E_ECAL_vs_E_HCAL =  new TH2F("h_E_ECAL_vs_E_HCAL","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_vs_E_HCAL->Sumw2();
  TH2F* h_E_ECAL_vs_E_HCALwHE =  new TH2F("h_E_ECAL_vs_E_HCALwHE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_vs_E_HCALwHE->Sumw2();
  TH2F* h_E_ECAL_vs_E_HCAL_HE =  new TH2F("h_E_ECAL_vs_E_HCAL_HE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_vs_E_HCAL_HE->Sumw2();

  TH2F* h_E_ECAL_EE_vs_E_HCAL =  new TH2F("h_E_ECAL_EE_vs_E_HCAL","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_EE_vs_E_HCAL->Sumw2();
  TH2F* h_E_ECAL_EE_vs_E_HCAL_HE =  new TH2F("h_E_ECAL_EE_vs_E_HCAL_HE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_EE_vs_E_HCAL_HE->Sumw2();
  TH2F* h_E_ECAL_EE_vs_E_HCALwHE =  new TH2F("h_E_ECAL_EE_vs_E_HCALwHE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECAL_EE_vs_E_HCALwHE->Sumw2();

  TH2F* h_E_ECALwEE_vs_E_HCAL =  new TH2F("h_E_ECALwEE_vs_E_HCAL","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECALwEE_vs_E_HCAL->Sumw2();
  TH2F* h_E_ECALwEE_vs_E_HCAL_HE =  new TH2F("h_E_ECALwEE_vs_E_HCAL_HE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECALwEE_vs_E_HCAL_HE->Sumw2();
  TH2F* h_E_ECALwEE_vs_E_HCALwHE =  new TH2F("h_E_ECALwEE_vs_E_HCALwHE","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECALwEE_vs_E_HCALwHE->Sumw2();

  TH2F* h_E_ECALwEEP_vs_E_HCALwHEP =  new TH2F("h_E_ECALwEEP_vs_E_HCALwHEP","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high, n_bins2D_E, 0,limit_E_bins_high_HCAL);
  h_E_ECALwEEP_vs_E_HCALwHEP->Sumw2();

  TH1F* h_E_ECAL =  new TH1F("h_E_ECAL","", n_bins2D_E, limit_E_bins_low,limit_E_bins_high);
  h_E_ECAL->Sumw2();
  h_E_ECAL->SetLineWidth(2);
  h_E_ECAL->SetMarkerStyle(20);

  TH1F* h_E_HCAL =  new TH1F("h_H_ECAL","", n_bins2D_E, 0,(limit_E_bins_high/2.));
  h_E_HCAL->Sumw2();
  h_E_HCAL->SetLineWidth(2);
  h_E_HCAL->SetMarkerStyle(20);

  TH1F* h_E_CAL =  new TH1F("h_E_CAL","", n_bins2D_E, limit_E_bins_low,(limit_E_bins_high+1.));
  h_E_CAL->Sumw2();
  h_E_CAL->SetLineWidth(2);
  h_E_CAL->SetMarkerStyle(20);

  TH1F* h_E_CAL_over_true =  new TH1F("h_E_CAL_over_true","", n_bins2D_E, limit_E_bins_rel_low ,limit_E_bins_rel_high );
  h_E_CAL_over_true->Sumw2();
  h_E_CAL_over_true->SetLineWidth(2);
  h_E_CAL_over_true->SetMarkerStyle(20);

  TH1F* h_E_RECO_over_true =  new TH1F("h_E_RECO_over_true","", n_bins2D_E, limit_E_bins_rel_low ,limit_E_bins_rel_high );
  h_E_RECO_over_true->Sumw2();
  h_E_RECO_over_true->SetLineWidth(2);
  h_E_RECO_over_true->SetMarkerStyle(20);

  TH1F* h_E_phRECO_over_true =  new TH1F("h_E_phRECO_over_true","", n_bins2D_E, limit_E_bins_rel_low ,limit_E_bins_rel_high );
  h_E_phRECO_over_true->Sumw2();
  h_E_phRECO_over_true->SetLineWidth(2);
  h_E_phRECO_over_true->SetMarkerStyle(20);


  TH1F* h_E_leakage_over_tot =  new TH1F("h_E_leakage_over_tot","", (int)(n_bins2D_E/2), limit_E_bis_leakage_low,limit_E_bins_leakage_high );
  h_E_leakage_over_tot->Sumw2();
  h_E_leakage_over_tot->SetLineWidth(2);
  h_E_leakage_over_tot->SetMarkerStyle(20);


  std::vector<float>hitEE_per_layer(40,0);
  std::vector<float>hitEB_per_layer(40,0);

  int count_EB=0;
  int count_EE=0;

      int count_fill_EB=0;
      int count_fill_EE=0;

      int delta_count_EB=0;
      int delta_count_EE=0;

      double theta_min=0.5;

  for(unsigned int i_entry=0;i_entry<tree_photons->GetEntries();i_entry++){
    tree_photons->GetEntry(i_entry);
    if(i_entry%1000==0){
      std::cout<<"in entry "<<i_entry<<std::endl;
    }
    float correction=1;

    bool veto_transition=false;
    for(unsigned int i=0;i<trueTheta->size();i++){
      if((*trueGenStatus)[i]==1 &&((*trueTheta)[i]<0.597 || ((*trueTheta)[i]>0.6125 && (*trueTheta)[i]<2.529) || (*trueTheta)[i]>2.545)){
	veto_transition=true;
	break;
      }
    }
    if(veto_transition){
      continue;
    }



    //transition region stuff


    //correction=std::sqrt(truePhPx*truePhPx+truePhPy*truePhPy)/sqrt(truePhPx*truePhPx+truePhPy*truePhPy+truePhPz*truePhPz);
    //check endcap inclination and barrel inclination of theta

    //check if in barrel 
    /*
    bool inner_barrel_pass_theta=false;
    bool inner_barrel_pass_phi=false;

    if((truePhTheta>(0.5*TMath::Pi()+0.175) && truePhTheta<(0.5*TMath::Pi()+0.55)) || (truePhTheta>(0.5*TMath::Pi()-0.55) && truePhTheta<(0.5*TMath::Pi()-0.175))){
      inner_barrel_pass_theta=true;
    }
    for(unsigned int i=0;i<13;i++){//for 12 periodic
      if( (fabs(truePhPhi-(i*TMath::Pi()/6.-TMath::Pi()))<0.05) || (fabs((i*TMath::Pi()/6.-TMath::Pi())-truePhPhi)<0.05)){
	inner_barrel_pass_phi=true;
      }
    }
    //inner_barrel_pass_phi=true;

    bool is_barrel=false;
    bool is_endcap=false;
    if(truePhE>0 && recoPhE->size()==1 && ( (truePhTheta>0.175 && truePhTheta<0.55) || (truePhTheta>(TMath::Pi()-0.55) && truePhTheta<(TMath::Pi()-0.175)))){
      count_EE+=1;
      //std::cout<<i_entry<<"?"<<truePhTheta<<" endcap first "<<0.175<<" to "<<0.55<<" OR from "<<(TMath::Pi()-0.55)<<" to "<<(TMath::Pi()-0.175)<<"/"<<count_EE<<std::endl;
      is_endcap=true;
    }
    if(truePhE>0 && recoPhE->size()==1 && inner_barrel_pass_theta && inner_barrel_pass_phi){
      count_EB+=1;
      is_barrel=true;
      //std::cout<<i_entry<<"?"<<truePhTheta<<" barrel theta first "<<(0.5*TMath::Pi()+0.175)<<" to "<<(0.5*TMath::Pi()+0.55)<<" OR from "<<(0.5*TMath::Pi()-0.55)<<" to "<<(0.5*TMath::Pi()-0.175)<<"/"<<count_EB<<std::endl;
    }
    
    if(truePhE>0 && recoPhE->size()==1 &&
       //(is_endcap || is_barrel)){
       fabs(truePhCosTheta)>0.25 && fabs(truePhCosTheta)<0.50){
      float photonE=0;
      float recoE=0;
      for(unsigned int s=0; s < recoPhE->size();s++){
	if((*recoPhPDGID)[s]==22){
	  photonE+=(*recoPhE)[s];
	}
	recoE+=(*recoPhE)[s];
	
	//h_E_RECO_over_true->Fill(recoE/truePhE);
	//h_E_phRECO_over_true->Fill(photonE/truePhE);
	if(recoPhE->size()==1){
	  //correction=(*recoPhEHitCorrection)[0];
	  h_E_RECO_over_true->Fill((*recoPhE)[0]/truePhE);
	}
	if((*recoPhPDGID)[0]==22 && recoPhE->size()==1){
	  h_E_phRECO_over_true->Fill((*recoPhE)[0]/truePhE);
	}
      }
     float totalHitESum=0;
     float totalHitEBSum=0;
     float totalHitEESum=0;
     
     float totalHitESumH=0;
     float totalHitHBSum=0;
     float totalHitHESum=0;
     
     
     //if(is_barrel){
       for(unsigned int h=0;h<hitHCALx->size();h++){
	 //std::cout<<"hit raw for HCAL HB "<<(*hitHCALEraw)[h]<<"/"<<0.3*0.000485<<std::endl;
	 //if((*hitHCALEraw)[h]>(0.3*0.000485)){
	 totalHitESumH+=(*hitHCALEraw)[h];
	 totalHitHBSum+=(*hitHCALEraw)[h];
	 //}//MPI cut
       }
       h_HB_hitsum->Fill(totalHitHBSum);

       //}//barrel restriction
       //if(count_EB<610 && is_barrel){
       for(unsigned int h=0;h<hitECALx->size();h++){
	 hitEB_per_layer[(*hitECALLayer)[h]]=hitEB_per_layer[(*hitECALLayer)[h]]+(*hitECALE)[h];
	 if((*hitECALEraw)[h]>(0.5*0.00054)){
	   totalHitESum+=(*hitECALEraw)[h];
	   totalHitEBSum+=(*hitECALEraw)[h];
	 }
	 h_layerEB->Fill((*hitECALLayer)[h]);
	 h_2D_hit_x_y->Fill((*hitECALx)[h],(*hitECALy)[h]);
	 h_2D_hit_z_r_barrel->Fill((*hitECALz)[h],sqrt((*hitECALy)[h]*(*hitECALy)[h]+(*hitECALx)[h]*(*hitECALx)[h]));
	 h_2D_hitE_z_r_barrel->Fill((*hitECALz)[h],sqrt((*hitECALy)[h]*(*hitECALy)[h]+(*hitECALx)[h]*(*hitECALx)[h]),(*hitECALE)[h]);
	 h_3D_hit_x_y_z_barrel->Fill((*hitECALx)[h],(*hitECALy)[h],(*hitECALz)[h]);
	 //if(count_fill_EB<7500){
	 h_layerEB_E->Fill((*hitECALLayer)[h],(*hitECALEraw)[h]);
	 //}
	 if((*recoPhPDGID)[0]==22){
	   h_layerEB_E_ph->Fill((*hitECALLayer)[h],(*hitECALEraw)[h]);
	   h_hit_E_phReco->Fill((*hitECALE)[h]);
	   h_2D_hit_x_y_phReco->Fill((*hitECALx)[h],(*hitECALy)[h]);
	   h_hit_Eraw_phReco->Fill((*hitECALEraw)[h]);
	 }else if((*recoPhPDGID)[0]==2112){
	   h_layerEB_E_n->Fill((*hitECALLayer)[h],(*hitECALEraw)[h]);
	   h_2D_hit_x_y_nReco->Fill((*hitECALx)[h],(*hitECALy)[h]);
	   h_hit_E_nReco->Fill((*hitECALE)[h]);
	   h_hit_Eraw_nReco->Fill((*hitECALEraw)[h]);
	 }
       }
       h_EB_hitsum->Fill(totalHitEBSum);
       //}//barrel restriction
       //if(count_EE<610 && is_endcap){
       for(unsigned int h=0;h<hitECALEx->size();h++){
	 hitEE_per_layer[(*hitECALELayer)[h]]=hitEE_per_layer[(*hitECALELayer)[h]]+(*hitECALEE)[h];
	 //if((*hitECALEEraw)[h]>(0.5*0.00054)){
	   totalHitESum+=(*hitECALEEraw)[h];
	   totalHitEESum+=(*hitECALEEraw)[h];
	   //}
	 h_layerEE->Fill((*hitECALELayer)[h]);
	 //if(count_fill_EE<7500){
	   h_layerEE_E->Fill((*hitECALELayer)[h],(*hitECALEEraw)[h]);
	 if((*recoPhPDGID)[0]==22){
	   h_layerEE_E_ph->Fill((*hitECALELayer)[h],(*hitECALEEraw)[h]);
	 }else if((*recoPhPDGID)[0]==2112){
	   h_layerEE_E_n->Fill((*hitECALELayer)[h],(*hitECALEEraw)[h]);
	 }
	   //}
	 h_2D_hit_x_y_endcap->Fill((*hitECALEx)[h],(*hitECALEy)[h]);
	 h_3D_hit_x_y_z_endcap->Fill((*hitECALEx)[h],(*hitECALEy)[h],(*hitECALEz)[h]);
	 h_2D_hit_r_z_endcap->Fill(sqrt((*hitECALEy)[h]*(*hitECALEy)[h]+(*hitECALEx)[h]*(*hitECALEx)[h]),(*hitECALEz)[h]);
	 h_2D_hitE_r_z_endcap->Fill(sqrt((*hitECALEy)[h]*(*hitECALEy)[h]+(*hitECALEx)[h]*(*hitECALEx)[h]),(*hitECALEz)[h],(*hitECALEE)[h]);
       }     
       h_EE_hitsum->Fill(totalHitEESum);



       //}//endcap restriction
       //if(is_endcap){
       for(unsigned int h=0;h<hitHCALEx->size();h++){
	 // std::cout<<"hit raw for HCAL HE "<<(*hitHCALEEraw)[h]<<"/"<<0.3*0.00047000<<std::endl;
	 //if((*hitHCALEEraw)[h]>(0.3*0.00047000)){	    
	 totalHitESumH+=(*hitHCALEEraw)[h];
	 totalHitHESum+=(*hitHCALEEraw)[h];	  
	 //}
       }     
       h_HE_hitsum->Fill(totalHitHESum);
       //}//endcap restriction
     if(is_endcap || is_barrel){
       h_E_hitsumH->Fill(totalHitESumH);
       h_E_hitsum->Fill(totalHitESum);
     }
     
     //for(unsigned int s=0; s < recoPhE->size();s++){
     for(unsigned int s=0; s < 1;s++){
       //if( ((*recoPhEHE)[0]>0 ||(*recoPhEEE)[0]>0) && ((*recoPhEHB)[0]>0 ||(*recoPhEEB)[0]>0)){
       if((*recoPhPDGID)[s]==22){
	 
	 h_dtheta_vs_phi_true->Fill((*recoPhTheta)[s]-truePhTheta,truePhPhi);
	 h_dphi_vs_theta_true->Fill(DeltaPhiDir((*recoPhPhi)[s],truePhPhi),truePhTheta);
	 h_dphi_vs_costheta_true->Fill(DeltaPhiDir((*recoPhPhi)[s],truePhPhi),cos(truePhTheta));
	 h_phi_phReco_vs_phi_true->Fill((*recoPhPhi)[s],truePhPhi);
	 h_dphi_vs_phi_true->Fill(DeltaPhiDir((*recoPhPhi)[s],truePhPhi),truePhPhi);
	 h_theta_phReco_vs_theta_true->Fill((*recoPhTheta)[s],truePhTheta);
	 h_dtheta_vs_theta_true->Fill((*recoPhTheta)[s]-truePhTheta,truePhTheta);
	 h_dtheta_vs_costheta_true->Fill((*recoPhTheta)[s]-truePhTheta,cos(truePhTheta));
	 h_dphi_phReco->Fill(DeltaPhiDir((*recoPhPhi)[s],truePhPhi));
	 h_dphi_logERW_phReco->Fill(DeltaPhiDir((*recoPhPhi_logERW)[s],truePhPhi));
	 h_dtheta_phReco->Fill((*recoPhTheta)[s]-truePhTheta);
	 h_dtheta_logERW_phReco->Fill((*recoPhTheta_logERW)[s]-truePhTheta);
	 h_phi_phReco->Fill((*recoPhPhi)[s]);
	 h_theta_phReco->Fill((*recoPhTheta)[s]);
	 h_CosTheta_phReco->Fill((*recoPhCosTheta)[s]);
	 h_E_phReco->Fill((*recoPhE)[s]);

	 //if(count_EB<610 && is_barrel){
	   if((*recoPhNhitsHB)[s]>0){
	     h_E_HB_Reco->Fill((*recoPhEHB)[s]);
	   }
	   if((*recoPhNhitsEB)[s]>0){
	     count_fill_EB+=1;
	     //std::cout<<i_entry<<" count EB/count fill "<<count_EB<<"/"<<count_fill_EB<<"/"<<truePhTheta<<"/"<<(*recoPhNhitsEB)[s]<<std::endl;
	     h_E_EB_phReco->Fill((*recoPhEEB)[s]);
	   }else{
	     std::cout<<i_entry<<"no count EB/count fill "<<count_EB<<"/"<<count_fill_EB<<"/"<<truePhTheta<<"/"<<(*recoPhNhitsEB)[s]<<std::endl;
	     h_E_EB_phReco->Fill((*recoPhEEB)[s]);
	   }
	   if((*recoPhNhitsHB)[s]>0){
	     h_E_HB_phReco->Fill((*recoPhEHB)[s]);
	   }
	   //}//barrel cuts
	   //if(count_EE<610 && is_endcap){
	   if((*recoPhNhitsHE)[s]>0){
	     h_E_HE_phReco->Fill((*recoPhEHE)[s]);
	   }
	   if((*recoPhNhitsHE)[s]>0){
	     h_E_HE_Reco->Fill((*recoPhEHE)[s]);
	   }
	   if((*recoPhNhitsEE)[s]>0){
	     if(truePhTheta<theta_min){
	       theta_min=truePhTheta;
	     }else if ( (TMath::Pi()-truePhTheta)<theta_min){
	       theta_min=TMath::Pi()-truePhTheta;
	     }
	     count_fill_EE+=1;
	     //std::cout<<" count EE/count fill "<<count_EE<<"/"<<count_fill_EE<<"/"<<truePhTheta<<"/"<<(*recoPhNhitsEE)[s] <<std::endl;
	     h_E_EE_phReco->Fill((*recoPhEEE)[s]);
	     //std::cout<<i_entry<<" count EE/count fill "<<count_EE<<"/"<<count_fill_EE<<"/"<<truePhTheta<<"/"<<(*recoPhNhitsEE)[s]<<std::endl;
	     h_E_EB_phReco->Fill((*recoPhEEB)[s]);
	   }//else{
	   //std::cout<<i_entry<<" no count EE/count fill "<<count_EE<<"/"<<count_fill_EE<<"/"<<truePhTheta<<"/"<<(*recoPhNhitsEE)[s] <<std::endl;
	   //}
	   //}//endcap cuts
	}else if((*recoPhPDGID)[s]==2112 && recoPhE->size()==1){//if split photons and neutrons
	 h_phi_nReco->Fill((*recoPhPhi)[s]);
	 h_theta_nReco->Fill((*recoPhTheta)[s]);
	 h_CosTheta_nReco->Fill((*recoPhCosTheta)[s]);
	 h_E_nReco->Fill((*recoPhE)[s]);
	 if((*recoPhNhitsEB)[s]>0){
	   h_E_EB_nReco->Fill((*recoPhEEB)[s]);
	 }
	 if((*recoPhNhitsEE)[s]>0){
	   h_E_EE_nReco->Fill((*recoPhEEE)[s]);
	 }
	 if((*recoPhNhitsHB)[s]>0){
	   h_E_HB_nReco->Fill((*recoPhEHB)[s]);
	 }
	 if((*recoPhNhitsHE)[s]>0){
	   h_E_HE_nReco->Fill((*recoPhEHE)[s]);
	 }
       }//neutron if
       
       if((count_EB- count_fill_EB)>delta_count_EB && count_EB<610){
	 delta_count_EB+=1;
	 //for(unsigned int t=0;t<recoPhE->size();t++){
	 //std::cout<<i_entry<<"event EB not there "<<t<<"/"<<truePhTheta<<"/"<<truePhPhi<<"/"<<count_EB<<"/"<<count_fill_EB<<"/"<<delta_count_EB<<"/"<<count_EE - count_fill_EE<<" N hits EB/EE "<<(*recoPhNhitsEB)[t]<<"/"<<(*recoPhNhitsEE)[t]<<"/"<<recoPhE->size()<<is_barrel<<std::endl;
	 //}
	 //std::cout<<"event EB not there "<<truePhTheta<<"/"<<truePhPhi<<"/"<<count_EB<<"/"<<count_fill_EB<<"/"<<delta_count_EB<<std::endl;
       }
       if((count_EE - count_fill_EE)>delta_count_EE && count_EE<610){
	 delta_count_EE+=1;
	 //for(unsigned int t=0;t<recoPhE->size();t++){
	 //std::cout<<i_entry<<"event EE not there "<<t<<"/"<<truePhTheta<<"/"<<truePhPhi<<"/"<<count_EE<<"/"<<count_fill_EE<<"/"<<delta_count_EE<<"/"<<count_EE - count_fill_EE<<" N hits EB/EE "<<(*recoPhNhitsEB)[t]<<"/"<<(*recoPhNhitsEE)[t]<<"/"<<recoPhE->size()<<is_endcap<<std::endl;
	 //}
       }
       
       if(fabs(truePhPhi)<0.5){
	 h_phiReco_phiTrue0->Fill((*recoPhPhi)[0]);
       }
     }//loop over pfo particles
     h_Delta_phiReco_phiTrue->Fill(DeltaPhiDir((*recoPhPhi)[0],truePhPhi));
     //DeltaPhiDi
     if(truePhE>0  &&  //(is_endcap|| (inner_barrel_pass_theta && inner_barrel_pass_phi))){
       //(truePhCosTheta)>0.50 && (truePhCosTheta)<0.75){
       fabs(truePhCosTheta)>0.25 && fabs(truePhCosTheta)<0.50){
       float hit_Eraw_sum_EB=0;
       float hit_Eraw_sum_EE=0;
       float hit_Eraw_sum_HB=0;
       float hit_Eraw_sum_HE=0;
       for(unsigned int i=0;i<hitECALEraw->size();i++){
	 hit_Eraw_sum_EB+=(*hitECALEraw)[i];
       }
       for(unsigned int i=0;i<hitECALEEraw->size();i++){
	 hit_Eraw_sum_EE+=(*hitECALEEraw)[i];
       }
       for(unsigned int i=0;i<hitHCALEraw->size();i++){
	 hit_Eraw_sum_HB+=(*hitHCALEraw)[i];
       }
       for(unsigned int i=0;i<hitHCALEEraw->size();i++){
	 hit_Eraw_sum_HE+=(*hitHCALEEraw)[i];
       }
       //std::cout<<"E raw/normed "<<hit_Eraw_sum_EB<<"/"<<hit_Eraw_sum_EE<<"/"<<hit_Eraw_sum_HB<<"/"<<hit_Eraw_sum_HE<<"/"<<(hit_Eraw_sum_EB+hit_Eraw_sum_EE+hit_Eraw_sum_HB+hit_Eraw_sum_HE)/truePhE<<"/"<< truePhE/(hit_Eraw_sum_EB+hit_Eraw_sum_EE+hit_Eraw_sum_HB+hit_Eraw_sum_HE)<<std::endl;
       h_hit_Eraw_totSum_over_trueE->Fill((hit_Eraw_sum_EB+hit_Eraw_sum_EE)/truePhE);
       //std::cout<<"correction factor "<<correction<<"/"<<truePhCosTheta<<"/"<<totalECALPhE<<"/"<<totalHCALPhE<<"/"<<totalECALPhE<<"/"<<totalHCALPhE*correction<<"/"<<(totalECALPhE+totalHCALPhE)<<"/"<<(totalECALPhE+correction*totalHCALPhE)<<std::endl;
       h_E_ECAL_vs_E_HCAL->Fill(totalECALPhE*weight_ecal_B,totalHCALPhE*weight_hcal_B);
       h_E_ECAL_vs_E_HCAL_HE->Fill(totalECALPhE*weight_ecal_B,totalEndcapHCALPhE*weight_hcal_E);
       h_E_ECAL_vs_E_HCALwHE->Fill(totalECALPhE*weight_ecal_B,totalHCALPhE*weight_hcal_B+totalEndcapHCALPhE*weight_hcal_E);
       h_E_ECAL_EE_vs_E_HCAL->Fill(totalEndcapECALPhE*weight_ecal_E,totalHCALPhE*weight_hcal_B);
       h_E_ECAL_EE_vs_E_HCAL_HE->Fill(totalEndcapECALPhE*weight_ecal_E,totalEndcapHCALPhE*weight_hcal_E);
       h_E_ECAL_EE_vs_E_HCALwHE->Fill(totalEndcapECALPhE*weight_ecal_E,totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B);
       h_E_ECALwEE_vs_E_HCAL->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E,totalHCALPhE*weight_hcal_B);
       h_E_ECALwEE_vs_E_HCAL_HE->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E,totalEndcapHCALPhE*weight_hcal_E);
       h_E_ECALwEE_vs_E_HCALwHE->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E,totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B);
       h_E_ECALwEEP_vs_E_HCALwHEP->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E+totalPlugECALPhE,totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B+totalPlugHCALPhE);
       
       h_E_ECAL->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E);
       h_E_HCAL->Fill(totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B);
       h_E_CAL->Fill(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E);
       h_E_CAL_over_true->Fill((totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E+totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B)/truePhE);
       h_E_leakage_over_tot->Fill( (totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B)/(totalECALPhE*weight_ecal_B+totalEndcapECALPhE*weight_ecal_E+totalEndcapHCALPhE*weight_hcal_E+totalHCALPhE*weight_hcal_B));      
     }
    }//true photon in certain region
    */
  }//end of tree filling
  float EE_sum=0;
  float EB_sum=0;
  for(unsigned int h=0;h<hitEE_per_layer.size();h++){
    EE_sum+=hitEE_per_layer[h];
    EB_sum+=hitEB_per_layer[h];
    //std::cout<<"layer/hitEB-E_sum/hitEE-E_sum "<<h<<"/"<<hitEB_per_layer[h]<<"/"<<hitEE_per_layer[h]<<std::endl;
  }

  for(unsigned int h=0;h<hitEE_per_layer.size();h++){
    std::cout<<"layer/hitEB-E_sum/hitEE-E_sum /EE_raw"<<h<<"/"<<hitEB_per_layer[h]<<"/"<<hitEE_per_layer[h]*EB_sum/EE_sum<<"/"<<hitEE_per_layer[h]<<" "<< theta_min<<std::endl;
  }

  setTDRStyle();

  std::cout<<"fill EB-E vs HB-E"<<count_fill_EB<<"/"<<count_fill_EE<<std::endl;
  
  std::cout<<"count EB-E vs HB-E"<<count_EB<<"/"<<count_EE<<std::endl;

  //h_EE_hitsum->Scale(h_E_hitsum->Integral()/h_EE_hitsum->Integral());
  //h_EB_hitsum->Scale(h_E_hitsum->Integral()/h_EB_hitsum->Integral());
  /*
  TCanvas* can_E_hitsum_HCAL_HE=setUpperCanvas("can_E_hitsum_HCAL_HE");
  can_E_hitsum_HCAL_HE->cd();
  h_HE_hitsum->Draw("hist");
  h_HE_hitsum->Fit("gaus");

  TCanvas* can_E_hitsum_HCAL_HB=setUpperCanvas("can_E_hitsum_HCAL_HB");
  can_E_hitsum_HCAL_HB->cd();
  h_HB_hitsum->Draw("hist");
  h_HB_hitsum->Fit("gaus");
  
  
  TCanvas* can_E_hitsum_ECAL_EE=setUpperCanvas("can_E_hitsum_ECAL_EE");
  can_E_hitsum_ECAL_EE->cd();
  h_EE_hitsum->Draw("hist");
  h_EE_hitsum->Fit("gaus");

  TCanvas* can_E_hitsum_ECAL_EB=setUpperCanvas("can_E_hitsum_ECAL_EB");
  can_E_hitsum_ECAL_EB->cd();
  h_EB_hitsum->Draw("hist");
  h_EB_hitsum->Fit("gaus");
  
  
  TCanvas* can_E_hitsum_ECAL=setUpperCanvas("can_E_hitsum_ECAL");
  can_E_hitsum_ECAL->cd();

  //h_E_hitsum->Draw("hist");
  TLegend* leg_E_EB_vs_E_EE_hitsum=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_EB_vs_E_EE_hitsum->SetBorderSize(0);
  leg_E_EB_vs_E_EE_hitsum->SetFillStyle(0);
  leg_E_EB_vs_E_EE_hitsum->AddEntry(h_EE_hitsum->DrawCopy("hist"),"EE raw hit-E sum");
  leg_E_EB_vs_E_EE_hitsum->AddEntry(h_EB_hitsum->DrawCopy("hist,same"),"EB raw hit-E sum");
  leg_E_EB_vs_E_EE_hitsum->Draw();

  h_E_EE_phReco->SetLineColor(kRed);

  TCanvas* can_E_ECAL_split=setUpperCanvas("can_E_ECAL_split");
  can_E_ECAL_split->cd();
  //h_E_PFO->Draw("hist");
  TLegend* leg_E_EB_vs_E_EE_PFO=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_EB_vs_E_EE_PFO->SetBorderSize(0);
  leg_E_EB_vs_E_EE_PFO->SetFillStyle(0);
  leg_E_EB_vs_E_EE_PFO->AddEntry(h_E_EE_phReco->DrawCopy("hist"),"endcap photon");
  leg_E_EB_vs_E_EE_PFO->AddEntry(h_E_EB_phReco->DrawCopy("hist,same"),"barrel photon");
  leg_E_EB_vs_E_EE_PFO->Draw();

   TCanvas* can_E_EE=setUpperCanvas("can_E_EE");
    can_E_EE->cd();
    h_E_EE_phReco->Draw("h");

    TCanvas* can_E_EB=setUpperCanvas("can_E_EB");
    can_E_EB->cd();
    h_E_EB_phReco->Draw("h");
  

  
  TCanvas* can_E_hitsum_HCAL=setUpperCanvas("can_E_hitsum_HCAL");
  can_E_hitsum_HCAL->cd();

  h_E_hitsumH->Draw("hist");
  h_HB_hitsum->Draw("hist,same");
  h_HE_hitsum->Draw("hist,same");
  */  
/*
  //h_layerEE->Scale(h_layerEB->Integral()/h_layerEE->Integral());
  TCanvas* can_hitsum_Layer=setUpperCanvas("can_hitsum_Layer");
  can_hitsum_Layer->cd();
  h_layerEB->Draw("hist");
  h_layerEE->Draw("hist,same");
  */
  /*
  TCanvas* can_E_ECAL_vs_E_HCAL=setUpperCanvas("can_E_ECAL_vs_E_HCAL");
  can_E_ECAL_vs_E_HCAL->cd();

  h_E_ECAL_vs_E_HCAL->GetXaxis()->SetTitle("E_ECAL [GeV]");
  h_E_ECAL_vs_E_HCAL->GetYaxis()->SetTitle("E_HCAL [GeV]");
  h_E_ECAL_vs_E_HCAL->GetYaxis()->SetTitleOffset(1.05);

  TLegend* leg_E_ECAL_vs_E_HCAL=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_vs_E_HCAL->SetBorderSize(0);
  leg_E_ECAL_vs_E_HCAL->SetFillStyle(0);
  leg_E_ECAL_vs_E_HCAL->AddEntry(h_E_ECAL_vs_E_HCAL->DrawCopy("colz"),"ECAL-E vs HCAL-E");
  leg_E_ECAL_vs_E_HCAL->Draw();

  h_E_ECAL_vs_E_HCAL->Fit("pol1");

  std::cout<<"EB-E vs HE-E"<<std::endl;
  */
  /*
  TCanvas* can_E_ECAL_vs_E_HCAL_HE=setUpperCanvas("can_E_ECAL_vs_E_HCAL_HE");
  can_E_ECAL_vs_E_HCAL_HE->cd();

  TLegend* leg_E_ECAL_vs_E_HCAL_HE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_vs_E_HCAL_HE->SetBorderSize(0);
  leg_E_ECAL_vs_E_HCAL_HE->SetFillStyle(0);
  leg_E_ECAL_vs_E_HCAL_HE->AddEntry(h_E_ECAL_vs_E_HCAL_HE->DrawCopy("colz"),"ECAL-E vs HCAL-HE");
  leg_E_ECAL_vs_E_HCAL_HE->Draw();

  h_E_ECAL_vs_E_HCAL_HE->Fit("pol1");
  */
  std::cout<<"EB-E vs HB_HE-E"<<std::endl;
  /*
  TCanvas* can_E_ECAL_vs_E_HCALwHE=setUpperCanvas("can_E_ECAL_vs_E_HCALwHE");
  can_E_ECAL_vs_E_HCALwHE->cd();

  TLegend* leg_E_ECAL_vs_E_HCALwHE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_vs_E_HCALwHE->SetBorderSize(0);
  leg_E_ECAL_vs_E_HCALwHE->SetFillStyle(0);
  leg_E_ECAL_vs_E_HCALwHE->AddEntry(h_E_ECAL_vs_E_HCALwHE->DrawCopy("colz"),"ECAL-E vs HCALwHE");
  leg_E_ECAL_vs_E_HCALwHE->Draw();

  h_E_ECAL_vs_E_HCALwHE->Fit("pol1");

  std::cout<<"EE-E vs HB-E"<<std::endl;

  TCanvas* can_E_ECAL_EE_vs_E_HCAL=setUpperCanvas("can_E_ECAL_EE_vs_E_HCAL");
  can_E_ECAL_EE_vs_E_HCAL->cd();

  h_E_ECAL_EE_vs_E_HCAL->GetXaxis()->SetTitle("E_ECAL_EE [GeV]");
  h_E_ECAL_EE_vs_E_HCAL->GetYaxis()->SetTitle("E_HCAL [GeV]");
  h_E_ECAL_EE_vs_E_HCAL->GetYaxis()->SetTitleOffset(1.05);

  TLegend* leg_E_ECAL_EE_vs_E_HCAL=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_EE_vs_E_HCAL->SetBorderSize(0);
  leg_E_ECAL_EE_vs_E_HCAL->SetFillStyle(0);
  leg_E_ECAL_EE_vs_E_HCAL->AddEntry(h_E_ECAL_EE_vs_E_HCAL->DrawCopy("colz"),"ECAL-E vs HCAL-E");
  leg_E_ECAL_EE_vs_E_HCAL->Draw();

  h_E_ECAL_EE_vs_E_HCAL->Fit("pol1");
  */
  /*
  std::cout<<"EE-E vs HE-E"<<std::endl;
  
  TCanvas* can_E_ECAL_EE_vs_E_HCAL_HE=setUpperCanvas("can_E_ECAL_EE_vs_E_HCAL_HE");
  can_E_ECAL_EE_vs_E_HCAL_HE->cd();

  TLegend* leg_E_ECAL_EE_vs_E_HCAL_HE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_EE_vs_E_HCAL_HE->SetBorderSize(0);
  leg_E_ECAL_EE_vs_E_HCAL_HE->SetFillStyle(0);
  leg_E_ECAL_EE_vs_E_HCAL_HE->AddEntry(h_E_ECAL_EE_vs_E_HCAL_HE->DrawCopy("colz"),"ECAL-E vs HCAL-HE");
  leg_E_ECAL_EE_vs_E_HCAL_HE->Draw();

  h_E_ECAL_EE_vs_E_HCAL_HE->Fit("pol1");

  std::cout<<"EE-E vs HB_HE-E"<<std::endl;
  */
  /*
  TCanvas* can_E_ECAL_EE_vs_E_HCALwHE=setUpperCanvas("can_E_ECAL_EE_vs_E_HCALwHE");
  can_E_ECAL_EE_vs_E_HCALwHE->cd();

  TLegend* leg_E_ECAL_EE_vs_E_HCALwHE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_EE_vs_E_HCALwHE->SetBorderSize(0);
  leg_E_ECAL_EE_vs_E_HCALwHE->SetFillStyle(0);
  leg_E_ECAL_EE_vs_E_HCALwHE->AddEntry(h_E_ECAL_EE_vs_E_HCALwHE->DrawCopy("colz"),"ECAL-E vs HCALwHE");
  leg_E_ECAL_EE_vs_E_HCALwHE->Draw();

  h_E_ECAL_EE_vs_E_HCALwHE->Fit("pol1");





  std::cout<<"EB+EE-E vs HB-E"<<std::endl;

  TCanvas* can_E_ECALwEE_vs_E_HCAL=setUpperCanvas("can_E_ECALwEE_vs_E_HCAL");
  can_E_ECALwEE_vs_E_HCAL->cd();

  h_E_ECALwEE_vs_E_HCAL->GetXaxis()->SetTitle("E_ECALwEE [GeV]");
  h_E_ECALwEE_vs_E_HCAL->GetYaxis()->SetTitle("E_HCAL [GeV]");
  h_E_ECALwEE_vs_E_HCAL->GetYaxis()->SetTitleOffset(1.05);

  TLegend* leg_E_ECALwEE_vs_E_HCAL=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECALwEE_vs_E_HCAL->SetBorderSize(0);
  leg_E_ECALwEE_vs_E_HCAL->SetFillStyle(0);
  leg_E_ECALwEE_vs_E_HCAL->AddEntry(h_E_ECALwEE_vs_E_HCAL->DrawCopy("colz"),"ECAL-E vs HCAL-E");
  leg_E_ECALwEE_vs_E_HCAL->Draw();

  h_E_ECALwEE_vs_E_HCAL->Fit("pol1");

  std::cout<<"EB+EE-E vs HE-E"<<std::endl;

  TCanvas* can_E_ECALwEE_vs_E_HCAL_HE=setUpperCanvas("can_E_ECALwEE_vs_E_HCAL_HE");
  can_E_ECALwEE_vs_E_HCAL_HE->cd();

  TLegend* leg_E_ECALwEE_vs_E_HCAL_HE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECALwEE_vs_E_HCAL_HE->SetBorderSize(0);
  leg_E_ECALwEE_vs_E_HCAL_HE->SetFillStyle(0);
  leg_E_ECALwEE_vs_E_HCAL_HE->AddEntry(h_E_ECALwEE_vs_E_HCAL_HE->DrawCopy("colz"),"ECAL-E vs HCAL-HE");
  leg_E_ECALwEE_vs_E_HCAL_HE->Draw();

  h_E_ECALwEE_vs_E_HCAL_HE->Fit("pol1");
 
 
  std::cout<<"EB+EE-E vs HB_HE-E"<<std::endl;

  TCanvas* can_E_ECALwEE_vs_E_HCALwHE=setUpperCanvas("can_E_ECALwEE_vs_E_HCALwHE");
  can_E_ECALwEE_vs_E_HCALwHE->cd();

  TLegend* leg_E_ECALwEE_vs_E_HCALwHE=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECALwEE_vs_E_HCALwHE->SetBorderSize(0);
  leg_E_ECALwEE_vs_E_HCALwHE->SetFillStyle(0);
  leg_E_ECALwEE_vs_E_HCALwHE->AddEntry(h_E_ECALwEE_vs_E_HCALwHE->DrawCopy("colz"),"ECAL-E vs HCALwHE");
  leg_E_ECALwEE_vs_E_HCALwHE->Draw();

  h_E_ECALwEE_vs_E_HCALwHE->Fit("pol1");

  std::cout<<"EB+EE+EP-E vs HB_HE_HP-E"<<std::endl;
  */

  /*
  TCanvas* can_E_ECALwEEP_vs_E_HCALwHEP=setUpperCanvas("can_E_ECALwEEP_vs_E_HCALwHEP");
  can_E_ECALwEEP_vs_E_HCALwHEP->cd();

  TLegend* leg_E_ECALwEEP_vs_E_HCALwHEP=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECALwEEP_vs_E_HCALwHEP->SetBorderSize(0);
  leg_E_ECALwEEP_vs_E_HCALwHEP->SetFillStyle(0);
  leg_E_ECALwEEP_vs_E_HCALwHEP->AddEntry(h_E_ECALwEEP_vs_E_HCALwHEP->DrawCopy("colz"),"ECAL-EP vs HCALwHEP");
  leg_E_ECALwEEP_vs_E_HCALwHEP->Draw();

  h_E_ECALwEEP_vs_E_HCALwHEP->Fit("pol1");
  
  */

    /*
   TCanvas* can_E_HE=setUpperCanvas("can_E_HE");
    can_E_HE->cd();
    h_E_HE_Reco->Draw("h");

    TCanvas* can_E_HB=setUpperCanvas("can_E_HB");
    can_E_HB->cd();
    h_E_HB_Reco->Draw("h");
    */
  /*
  TCanvas* can_hits=setUpperCanvas("can_hits");
  can_hits->cd();
  h_2D_hit_x_y->Draw("colz");

  TCanvas* can_hits_endcap=setUpperCanvas("can_hits_endcap");
  can_hits_endcap->cd();
  h_2D_hit_x_y_endcap->Draw("colz");
  */
/*
  TCanvas* can_hits_endcap_3D=setUpperCanvas("can_hits_endcap_3D");
  can_hits_endcap_3D->cd();
  //h_3D_hit_x_y_z_endcap->Draw("ISO");

  TCanvas* can_hits_barrel_3D=setUpperCanvas("can_hits_barrel_3D");
  can_hits_barrel_3D->cd();
  //h_3D_hit_x_y_z_barrel->Draw("ISO");

  TCanvas* can_hits_nReco=setUpperCanvas("can_hits_nReco");
  can_hits_nReco->cd();
  h_2D_hit_x_y_nReco->Draw("colz");

  TCanvas* can_hits_phReco=setUpperCanvas("can_hits_phReco");
  can_hits_phReco->cd();
  h_2D_hit_x_y_phReco->Draw("colz");
    */
  /*
  THStack* thetastack= new THStack("theta-stack", "");
  //h_theta_phReco->SetFillColor(kBlue);
  //h_theta_nReco->SetFillColor(kRed);
  h_theta_phReco->GetXaxis()->SetTitle("#theta");
  thetastack->Add(h_theta_phReco);
  thetastack->Add(h_theta_nReco);

  TCanvas* can_theta_ph_and_n=setUpperCanvas("can_theta_ph_and_n");
  thetastack->Draw();

  TLegend* leg_theta_ph_and_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_theta_ph_and_n->SetBorderSize(0);
  leg_theta_ph_and_n->SetFillStyle(0);
  leg_theta_ph_and_n->AddEntry(h_theta_phReco,"PFA-photons");
  leg_theta_ph_and_n->AddEntry(h_theta_nReco,"PFA-neutrons");
  leg_theta_ph_and_n->Draw();

  std::cout<<"ph/n numbers "<<h_theta_phReco->Integral()<<"/"<<h_theta_nReco->Integral()<<" ratio "<<h_theta_nReco->Integral()/h_theta_phReco->Integral()<<std::endl;
  
  TCanvas* can_theta_ph_vs_n=setUpperCanvas("can_theta_ph_vs_n");
  can_theta_ph_vs_n->cd();



  TLegend* leg_theta_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_theta_ph_vs_n->SetBorderSize(0);
  leg_theta_ph_vs_n->SetFillStyle(0);
  leg_theta_ph_vs_n->AddEntry(h_theta_phReco->DrawCopy("hist"),"PFA-photons");
  leg_theta_ph_vs_n->AddEntry(h_theta_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_theta_ph_vs_n->Draw();


  TCanvas* can_CosTheta_ph_vs_n=setUpperCanvas("can_CosTheta_ph_vs_n");
  can_CosTheta_ph_vs_n->cd();

  TLegend* leg_CosTheta_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_CosTheta_ph_vs_n->SetBorderSize(0);
  leg_CosTheta_ph_vs_n->SetFillStyle(0);
  leg_CosTheta_ph_vs_n->AddEntry(h_CosTheta_phReco->DrawCopy("hist"),"PFA-photons");
  leg_CosTheta_ph_vs_n->AddEntry(h_CosTheta_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_CosTheta_ph_vs_n->Draw();

  TCanvas* can_ratio_CosTheta_ph_vs_n=setUpperCanvas("can_ratio_CosTheta_n_vs_ph");
  can_ratio_CosTheta_ph_vs_n->cd();
  TH1F* h_CosTheta_n_vs_phRECO=(TH1F*)h_CosTheta_nReco->Clone("ratio_costheta_n_vs_ph");
  h_CosTheta_n_vs_phRECO->Divide(h_CosTheta_phReco);
  h_CosTheta_n_vs_phRECO->DrawCopy("hist");

  THStack* phistack= new THStack("phi-stack", "");
  //h_phi_phReco->SetFillColor(kBlue);
  //h_phi_nReco->SetFillColor(kRed);
  h_phi_phReco->GetXaxis()->SetTitle("#phi");
  phistack->Add(h_phi_phReco);
  phistack->Add(h_phi_nReco);
  
  TCanvas* can_phi_ph_and_n=setUpperCanvas("can_phi_ph_and_n");
  phistack->Draw();

  TLegend* leg_phi_ph_and_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_phi_ph_and_n->SetBorderSize(0);
  leg_phi_ph_and_n->SetFillStyle(0);
  leg_phi_ph_and_n->AddEntry(h_phi_phReco,"PFA-photons");
  leg_phi_ph_and_n->AddEntry(h_phi_nReco,"PFA-neutrons");
  leg_phi_ph_and_n->Draw();
  
  TCanvas* can_phi_ph_vs_n=setUpperCanvas("can_phi_ph_vs_n");
  can_phi_ph_vs_n->cd();

  TLegend* leg_phi_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_phi_ph_vs_n->SetBorderSize(0);
  leg_phi_ph_vs_n->SetFillStyle(0);
  leg_phi_ph_vs_n->AddEntry(h_phi_phReco->DrawCopy("hist"),"PFA-photons");
  leg_phi_ph_vs_n->AddEntry(h_phi_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_phi_ph_vs_n->Draw();
  */  

/*
  TCanvas* can_E_ph_vs_n=setUpperCanvas("can_E_ph_vs_n");
  can_E_ph_vs_n->cd();

  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral());
  
  TLegend* leg_E_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ph_vs_n->SetBorderSize(0);
  leg_E_ph_vs_n->SetFillStyle(0);
  leg_E_ph_vs_n->AddEntry(h_E_phReco->DrawCopy("hist"),"PFA-photons");
  leg_E_ph_vs_n->AddEntry(h_E_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_E_ph_vs_n->Draw();

  TCanvas* can_E_EB_ph_vs_n=setUpperCanvas("can_E_EB_ph_vs_n");
  can_E_EB_ph_vs_n->cd();

  //h_E_EB_nReco->Scale(h_E_EB_phReco->Integral()/h_E_EB_nReco->Integral());
  
  TLegend* leg_E_EB_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_EB_ph_vs_n->SetBorderSize(0);
  leg_E_EB_ph_vs_n->SetFillStyle(0);
  leg_E_EB_ph_vs_n->AddEntry(h_E_EB_phReco->DrawCopy("hist"),"PFA-photons");
  leg_E_EB_ph_vs_n->AddEntry(h_E_EB_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_E_EB_ph_vs_n->Draw();

  TCanvas* can_E_HB_ph_vs_n=setUpperCanvas("can_E_HB_ph_vs_n");
  can_E_HB_ph_vs_n->cd();

  //h_E_HB_nReco->Scale(h_E_HB_phReco->Integral()/h_E_HB_nReco->Integral());
  
  TLegend* leg_E_HB_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_HB_ph_vs_n->SetBorderSize(0);
  leg_E_HB_ph_vs_n->SetFillStyle(0);
  leg_E_HB_ph_vs_n->AddEntry(h_E_HB_phReco->DrawCopy("hist"),"PFA-photons");
  leg_E_HB_ph_vs_n->AddEntry(h_E_HB_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_E_HB_ph_vs_n->Draw();

  TCanvas* can_hit_E_ph_vs_n=setUpperCanvas("can_hit_E_ph_vs_n");
  can_hit_E_ph_vs_n->cd();

  //h_hit_E_nReco->Scale(h_hit_E_phReco->Integral()/h_hit_E_nReco->Integral());
  
  TLegend* leg_hit_E_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_hit_E_ph_vs_n->SetBorderSize(0);
  leg_hit_E_ph_vs_n->SetFillStyle(0);
  leg_hit_E_ph_vs_n->AddEntry(h_hit_E_phReco->DrawCopy("hist"),"PFA-photons");
  leg_hit_E_ph_vs_n->AddEntry(h_hit_E_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_hit_E_ph_vs_n->Draw();

  TCanvas* can_hit_Eraw_ph_vs_n=setUpperCanvas("can_hit_Eraw_ph_vs_n");
  can_hit_Eraw_ph_vs_n->cd();

  //h_hit_Eraw_nReco->Scale(h_hit_Eraw_phReco->Integral()/h_hit_Eraw_nReco->Integral());
  
  TLegend* leg_hit_Eraw_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_hit_Eraw_ph_vs_n->SetBorderSize(0);
  leg_hit_Eraw_ph_vs_n->SetFillStyle(0);
  leg_hit_Eraw_ph_vs_n->AddEntry(h_hit_Eraw_phReco->DrawCopy("hist"),"PFA-photons");
  leg_hit_Eraw_ph_vs_n->AddEntry(h_hit_Eraw_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_hit_Eraw_ph_vs_n->Draw();

  TCanvas* can_hit_Eraw_totSum_over_trueE=setUpperCanvas("can_hit_Eraw_totSum_over_trueE");
  can_hit_Eraw_totSum_over_trueE->cd();
  h_hit_Eraw_totSum_over_trueE->Draw();


  TCanvas* can_nhits_ph_vs_n=setUpperCanvas("can_nhits_ph_vs_n");
  can_nhits_ph_vs_n->cd();

  //h_nhits_nReco->Scale(h_nhits_phReco->Integral()/h_nhits_nReco->Integral());
  
  TLegend* leg_nhits_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_nhits_ph_vs_n->SetBorderSize(0);
  leg_nhits_ph_vs_n->SetFillStyle(0);
  leg_nhits_ph_vs_n->AddEntry(h_nhits_phReco->DrawCopy("hist"),"PFA-photons");
  leg_nhits_ph_vs_n->AddEntry(h_nhits_nReco->DrawCopy("hist,same"),"PFA-neutrons");
  leg_nhits_ph_vs_n->Draw();


    */
  /*
  TCanvas* can_E_ECAL_reco=setUpperCanvas("can_E_ECAL");
  can_E_ECAL_reco->cd();

 h_E_ECAL->GetXaxis()->SetTitle("E_ECAL [GeV]");

  TLegend* leg_E_ECAL_reco=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_ECAL_reco->SetBorderSize(0);
  leg_E_ECAL_reco->SetFillStyle(0);
  leg_E_ECAL_reco->AddEntry(h_E_ECAL->DrawCopy("P,e,l"),"ECAL-E");
  leg_E_ECAL_reco->Draw();

  TCanvas* can_E_HCAL_reco=setUpperCanvas("can_E_HCAL");
  can_E_HCAL_reco->cd();

  h_E_HCAL->GetXaxis()->SetTitle("E_HCAL [GeV]");

  TLegend* leg_E_HCAL_reco=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_HCAL_reco->SetBorderSize(0);
  leg_E_HCAL_reco->SetFillStyle(0);
  leg_E_HCAL_reco->AddEntry(h_E_HCAL->DrawCopy("P,e,l"),"HCAL-E");
  leg_E_HCAL_reco->Draw();

  TCanvas* can_E_CAL_reco=setUpperCanvas("can_E_CAL");
  can_E_CAL_reco->cd();
  */
  /*
  TLegend* leg_E_CAL_reco=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_CAL_reco->SetBorderSize(0);
  leg_E_CAL_reco->SetFillStyle(0);
  leg_E_CAL_reco->AddEntry(h_E_CAL->DrawCopy("l"),"ECAL-E");
  leg_E_CAL_reco->Draw();

  h_E_CAL->Fit("gaus");
  h_E_CAL->Draw("P,E,same");
  h_E_CAL->GetXaxis()->SetTitle("E_ECAL+E_HCAL [GeV]");
  */
  /*
  TCanvas* can_dphi_vs_phi_true=setUpperCanvas("can_dphi_vs_phi_true");
  can_dphi_vs_phi_true->cd();
  h_dphi_vs_phi_true->Draw("colz");

  TCanvas* can_dtheta_vs_theta_true=setUpperCanvas("can_dtheta_vs_theta_true");
  can_dtheta_vs_theta_true->cd();
  h_dtheta_vs_theta_true->Draw("colz");

  TCanvas* can_dtheta_vs_costheta_true=setUpperCanvas("can_dtheta_vs_costheta_true");
  can_dtheta_vs_costheta_true->cd();
  h_dtheta_vs_costheta_true->Draw("colz");

  TCanvas* can_phi_phReco_vs_phi_true=setUpperCanvas("can_phi_phReco_vs_phi_true");
  can_phi_phReco_vs_phi_true->cd();
  h_phi_phReco_vs_phi_true->Draw("colz");

  TCanvas* can_theta_phReco_vs_theta_true=setUpperCanvas("can_theta_phReco_vs_theta_true");
  can_theta_phReco_vs_theta_true->cd();
  h_theta_phReco_vs_theta_true->Draw("colz");

  TCanvas* can_dphi_vs_costheta_true=setUpperCanvas("can_dphi_vs_costheta_true");
  can_dphi_vs_costheta_true->cd();
  h_dphi_vs_costheta_true->Draw("colz");


  TCanvas* can_dtheta_vs_phi_true=setUpperCanvas("can_dtheta_vs_phi_true");
  can_dtheta_vs_phi_true->cd();
  h_dtheta_vs_phi_true->Draw("colz");

  TCanvas* can_dphi_vs_theta_true=setUpperCanvas("can_dphi_vs_theta_true");
  can_dphi_vs_theta_true->cd();
  h_dphi_vs_theta_true->Draw("colz");
  
  TCanvas* can_E_RECO_vs_true=setUpperCanvas("can_E_RECO_vs_true");
  can_E_RECO_vs_true->cd();

  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral()); 
  TLegend* leg_E_RECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_RECO_vs_true->SetBorderSize(0);
  leg_E_RECO_vs_true->SetFillStyle(0);
  leg_E_RECO_vs_true->AddEntry(h_E_RECO_over_true->DrawCopy("P,e,l"),"E_{#gamma}^{PFA}/E_{true}");
  leg_E_RECO_vs_true->Draw();
  
  h_E_RECO_over_true->Fit("gaus");
  h_E_RECO_over_true->Draw("P,E,same");
  h_E_RECO_over_true->GetXaxis()->SetTitle("E_{#gamma}^{PFA}/E_true");
  
  TF1* fit_RECO_resolution = (TF1*)h_E_RECO_over_true->GetFunction("gaus");

  //std::cout<<"RECO rms/error/fit mean/fit sigma "<< h_E_RECO_over_true->GetRMS()<<"/"<<h_E_RECO_over_true->GetRMSError()<<"/"<<"/"<<fit_RECO_resolution->GetParameter(0)<<"/"<<fit_RECO_resolution->GetParameter(1)<<"/"<<fit_RECO_resolution->GetParError(1)<<"/"<<fit_RECO_resolution->GetParameter(2)<<"/"<<fit_RECO_resolution->GetParError(2)<<std::endl;
  */
  TCanvas* can_dphi_phRECO_vs_true=setUpperCanvas("can_dphi_phRECO_vs_true");
  can_dphi_phRECO_vs_true->cd();
  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral());                                                                                                           
  TLegend* leg_dphi_phRECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_dphi_phRECO_vs_true->SetBorderSize(0);
  leg_dphi_phRECO_vs_true->SetFillStyle(0);
  leg_dphi_phRECO_vs_true->AddEntry(h_dphi_phReco->DrawCopy("P,e,l"),"#Phi^{PFA}-#Phi_{true}");
  leg_dphi_phRECO_vs_true->Draw();

  h_dphi_phReco->Fit("gaus");
  h_dphi_phReco->Draw("P,E,same");
  h_dphi_phReco->GetXaxis()->SetTitle("#Phi(true)-#Phi(Ph-PFO)");

  TF1* fit_phRECO_dphi_resolution = (TF1*)h_dphi_phReco->GetFunction("gaus");

  std::cout<<"phRECO dphi rms/error/fit mean/fit sigma "<< h_dphi_phReco->GetRMS()<<"/"<<h_dphi_phReco->GetRMSError()<<"/"<<fit_phRECO_dphi_resolution->GetParameter(1)<<"/"<<fit_phRECO_dphi_resolution->GetParameter(2)<<"/"<<fit_phRECO_dphi_resolution->GetParError(2)<<std::endl;


  TCanvas* can_dphi_logERW_phRECO_vs_true=setUpperCanvas("can_dphi_logERW_phRECO_vs_true");
  can_dphi_logERW_phRECO_vs_true->cd();
  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral());                                                                                                           
  TLegend* leg_dphi_logERW_phRECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_dphi_logERW_phRECO_vs_true->SetBorderSize(0);
  leg_dphi_logERW_phRECO_vs_true->SetFillStyle(0);
  leg_dphi_logERW_phRECO_vs_true->AddEntry(h_dphi_logERW_phReco->DrawCopy("P,e,l"),"#Phi^{PFA}-#Phi_{true}");
  leg_dphi_logERW_phRECO_vs_true->Draw();

  h_dphi_logERW_phReco->Fit("gaus");
  h_dphi_logERW_phReco->Draw("P,E,same");
  h_dphi_logERW_phReco->GetXaxis()->SetTitle("#Phi(true)-#Phi(Ph-PFO)");

  TF1* fit_phRECO_dphi_logERW_resolution = (TF1*)h_dphi_logERW_phReco->GetFunction("gaus");

  std::cout<<"phRECO dphi_logERW rms/error/fit mean/fit sigma "<< h_dphi_logERW_phReco->GetRMS()<<"/"<<h_dphi_logERW_phReco->GetRMSError()<<"/"<<fit_phRECO_dphi_logERW_resolution->GetParameter(1)<<"/"<<fit_phRECO_dphi_logERW_resolution->GetParameter(2)<<"/"<<fit_phRECO_dphi_logERW_resolution->GetParError(2)<<std::endl;

  TCanvas* can_dtheta_phRECO_vs_true=setUpperCanvas("can_dtheta_phRECO_vs_true");
  can_dtheta_phRECO_vs_true->cd();
  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral());              
  TLegend* leg_dtheta_phRECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_dtheta_phRECO_vs_true->SetBorderSize(0);
  leg_dtheta_phRECO_vs_true->SetFillStyle(0);
  leg_dtheta_phRECO_vs_true->AddEntry(h_dtheta_phReco->DrawCopy("P,e,l"),"#Theta^{PFA}-#Theta_{true}");
  leg_dtheta_phRECO_vs_true->Draw();

  h_dtheta_phReco->Fit("gaus");
  h_dtheta_phReco->Draw("P,E,same");
  h_dtheta_phReco->GetXaxis()->SetTitle("#Theta(true)-#Theta(Ph-PFO)");

  TF1* fit_phRECO_dtheta_resolution = (TF1*)h_dtheta_phReco->GetFunction("gaus");

  std::cout<<"phRECO dtheta rms/error/fit mean/fit sigma "<< h_dtheta_phReco->GetRMS()<<"/"<<h_dtheta_phReco->GetRMSError()<<"/"<<fit_phRECO_dtheta_resolution->GetParameter(1)<<"/"<<fit_phRECO_dtheta_resolution->GetParameter(2)<<"/"<<fit_phRECO_dtheta_resolution->GetParError(2)<<std::endl;

  TCanvas* can_dtheta_logERW_phRECO_vs_true=setUpperCanvas("can_dtheta_logERW_phRECO_vs_true");
  can_dtheta_logERW_phRECO_vs_true->cd();
  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral());              
  TLegend* leg_dtheta_logERW_phRECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_dtheta_logERW_phRECO_vs_true->SetBorderSize(0);
  leg_dtheta_logERW_phRECO_vs_true->SetFillStyle(0);
  leg_dtheta_logERW_phRECO_vs_true->AddEntry(h_dtheta_logERW_phReco->DrawCopy("P,e,l"),"#Theta^{PFA}-#Theta_{true}");
  leg_dtheta_logERW_phRECO_vs_true->Draw();

  h_dtheta_logERW_phReco->Fit("gaus");
  h_dtheta_logERW_phReco->Draw("P,E,same");
  h_dtheta_logERW_phReco->GetXaxis()->SetTitle("#Theta(true)-#Theta(Ph-PFO)");

  TF1* fit_phRECO_dtheta_logERW_resolution = (TF1*)h_dtheta_logERW_phReco->GetFunction("gaus");

  std::cout<<"phRECO dtheta_logERW rms/error/fit mean/fit sigma "<< h_dtheta_logERW_phReco->GetRMS()<<"/"<<h_dtheta_logERW_phReco->GetRMSError()<<"/"<<fit_phRECO_dtheta_logERW_resolution->GetParameter(1)<<"/"<<fit_phRECO_dtheta_logERW_resolution->GetParameter(2)<<"/"<<fit_phRECO_dtheta_logERW_resolution->GetParError(2)<<std::endl;


  
  TCanvas* can_E_phRECO_vs_true=setUpperCanvas("can_E_phRECO_vs_true");
  can_E_phRECO_vs_true->cd();

  //h_E_nReco->Scale(h_E_phReco->Integral()/h_E_nReco->Integral()); 
  TLegend* leg_E_phRECO_vs_true=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_phRECO_vs_true->SetBorderSize(0);
  leg_E_phRECO_vs_true->SetFillStyle(0);
  leg_E_phRECO_vs_true->AddEntry(h_E_phRECO_over_true->DrawCopy("P,e,l"),"E_{#gamma}^{PFA}/E_{true}");
  leg_E_phRECO_vs_true->Draw();
  
  h_E_phRECO_over_true->Fit("gaus");
  h_E_phRECO_over_true->Draw("P,E,same");
  h_E_phRECO_over_true->GetXaxis()->SetTitle("E_{#gamma}^{PFA}/E_true");
  
  TF1* fit_phRECO_resolution = (TF1*)h_E_phRECO_over_true->GetFunction("gaus");

  std::cout<<"phRECO rms/error/fit mean/fit sigma "<< 100.*h_E_phRECO_over_true->GetRMS()<<"/"<<100.*h_E_phRECO_over_true->GetRMSError()<<"/"<<fit_phRECO_resolution->GetParameter(1)<<"/"<<100.*fit_phRECO_resolution->GetParameter(2)<<"/"<<100.*fit_phRECO_resolution->GetParError(2)<<std::endl;

  
  TCanvas* can_E_CAL_over_true_reco=setUpperCanvas("can_E_CAL_over_true");
  can_E_CAL_over_true_reco->cd();

  TLegend* leg_E_CAL_over_true_reco=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_CAL_over_true_reco->SetBorderSize(0);
  leg_E_CAL_over_true_reco->SetFillStyle(0);
  leg_E_CAL_over_true_reco->AddEntry(h_E_CAL_over_true->DrawCopy("P,e,l"),"CAL-E");
  leg_E_CAL_over_true_reco->Draw();

  h_E_CAL_over_true->Fit("gaus");
  h_E_CAL_over_true->Draw("P,E,same");
  h_E_CAL_over_true->GetXaxis()->SetTitle("(E_ECAL+E_HCAL)/E_true");

  TF1* fit_resolution = (TF1*)h_E_CAL_over_true->GetFunction("gaus");

  std::cout<<"rms/error/fit mean/fit sigma "<< 100.*h_E_CAL_over_true->GetRMS()<<"/"<<100.*h_E_CAL_over_true->GetRMSError()<<"/"<<100.*fit_resolution->GetParameter(1)<<"/"<<100.*fit_resolution->GetParameter(2)<<"/"<<100.*fit_resolution->GetParError(2)<<std::endl;




  

  /*
  TCanvas* can_E_leakage_over_tot_reco=setUpperCanvas("can_E_leakage_over_tot");
  can_E_leakage_over_tot_reco->cd();

  h_E_leakage_over_tot->GetXaxis()->SetTitle("E_HCAL/(E_ECAL+E_HCAL)");

  TLegend* leg_E_leakage_over_tot_reco=new TLegend(0.45,0.65,0.89,0.89,label_legend.c_str());
  leg_E_leakage_over_tot_reco->SetBorderSize(0);
  leg_E_leakage_over_tot_reco->SetFillStyle(0);
  leg_E_leakage_over_tot_reco->AddEntry(h_E_leakage_over_tot->DrawCopy("P,e,l"),"leakage");
  leg_E_leakage_over_tot_reco->Draw();

  TCanvas* can_sigma_overview_ECAL=setUpperCanvas("can_sigma_overview_ECAL");
  can_sigma_overview_ECAL->cd();


  Int_t n_sigma=11;
  Double_t x_sigma[n_sigma],y_sigma[n_sigma];
  
  x_sigma[0]=1;
  y_sigma[0]=15.14;
  x_sigma[1]=5;
  y_sigma[1]=6.95;
  x_sigma[2]=10;
  y_sigma[2]=4.873;
  x_sigma[3]=15;
  y_sigma[3]=4.00;
  x_sigma[4]=30;
  y_sigma[4]=2.829;
  x_sigma[5]=50;
  y_sigma[5]=2.281;
  x_sigma[6]=100;
  y_sigma[6]=1.69;
  x_sigma[7]=200;
  y_sigma[7]=1.27;
  x_sigma[8]=500;
  y_sigma[8]=0.955;
  x_sigma[9]=1000;
  y_sigma[9]=0.810;
  x_sigma[10]=1500;
  y_sigma[10]=0.775;

  TGraph* h_sigma_overview =  new TGraph(n_sigma,x_sigma,y_sigma);
  h_sigma_overview->SetLineWidth(2);
  h_sigma_overview->SetMarkerStyle(20);

   TGraphErrors* gre = new TGraphErrors(11);
   gre->SetName("CLICdet40_o3_V05X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineColor(kGreen);
   gre->SetMarkerColor(kGreen);
   gre->SetMarkerStyle(24);
  */
   /*
   gre->SetPoint(0,1,17.5813);
   gre->SetPointError(0,0,0.1464013);
   gre->SetPoint(1,10,5.89554);
   gre->SetPointError(1,0,0.04637848);
   gre->SetPoint(2,50,2.893315);
   gre->SetPointError(2,0,0.02121596);
   gre->SetPoint(3,100,2.254636);
   gre->SetPointError(3,0,0.01670716);
   gre->SetPoint(4,200,1.641406);
   gre->SetPointError(4,0,0.01258645);
   gre->SetPoint(5,500,1.144497);
   gre->SetPointError(5,0,0.008509384);
   gre->SetPoint(6,1000,0.8163775);
   gre->SetPointError(6,0,0.006066276);
   gre->SetPoint(7,1500,0.6745621);
   gre->SetPointError(7,0,0.00512168);
   */







  /*
   0.0-0.25 PFO sigma
   gre->SetPoint(0,1,16.595);
   gre->SetPointError(0,0,0.27116);
   gre->SetPoint(1,5,7.22371);
   gre->SetPointError(1,0,0.116276);
   gre->SetPoint(2,10,5/13621);
   gre->SetPointError(2,0,0.0962115);
   gre->SetPoint(3,15,4.30681);
   gre->SetPointError(3,0,0.0832983);
   gre->SetPoint(4,30,2.99232);
   gre->SetPointError(4,0,0.0503611);
   gre->SetPoint(5,50,2.42814);
   gre->SetPointError(5,0,0.0436082);
   gre->SetPoint(6,100,1.75307);
   gre->SetPointError(6,0,0.0309657);
   gre->SetPoint(7,200,1.28701);
   gre->SetPointError(7,0,0.021063);
   gre->SetPoint(8,500,1.03635);
   gre->SetPointError(8,0,0.0190124);
   gre->SetPoint(9,1000,0.916961);
   gre->SetPointError(9,0,0.0179033);
   gre->SetPoint(10,1500,0.885898);
   gre->SetPointError(10,0,0.016179);


   //0.0-0.25 hit level RMS
   gre->SetPoint(0,1,20.2433);
   gre->SetPointError(0,0,0.202433);
   gre->SetPoint(1,5,6.69045);
   gre->SetPointError(1,0,0.0932835);
   gre->SetPoint(2,10,4.72093);
   gre->SetPointError(2,0,0.0667908);
   gre->SetPoint(3,15,3.97856);
   gre->SetPointError(3,0,0.0570348);
   gre->SetPoint(4,30,2.92142);
   gre->SetPointError(4,0,0.0417091);
   gre->SetPoint(5,50,2.26957);
   gre->SetPointError52,0,0.0318366);
   gre->SetPoint(6,100,1.70603);
   gre->SetPointError(6,0,0.0240932);
   gre->SetPoint(7,200,1.2974);
   gre->SetPointError(7,0,0.0183333);
   gre->SetPoint(8,500,1.09089);
   gre->SetPointError(8,0,0.0153631);
   gre->SetPoint(9,1000,0.946449);
   gre->SetPointError(9,0,0.0135595);
   gre->SetPoint(10,1500,0.941685);
   gre->SetPointError(10,0,0.0133228);

   //0.0-0.25 hit level sigma
   gre->SetPoint(0,1,15.2607);
   gre->SetPointError(0,0,0.246774);
   gre->SetPoint(1,5,6.46757);
   gre->SetPointError(1,0,0.0941419);
   gre->SetPoint(2,10,4.54804);
   gre->SetPointError(2,0,0.0718447);
   gre->SetPoint(3,15,3.89062);
   gre->SetPointError(3,0,0.0624742);
   gre->SetPoint(4,30,2.89322);
   gre->SetPointError(4,0,0.0453367);
   gre->SetPoint(5,50,2.21127);
   gre->SetPointError(5,0,0.0329923);
   gre->SetPoint(6,100,1.66706);
   gre->SetPointError(6,0,0.0253022);
   gre->SetPoint(7,200,1.23588);
   gre->SetPointError(7,0,0.018947);
   gre->SetPoint(8,500,0.942072);
   gre->SetPointError(8,0,0.014366);
   gre->SetPoint(9,1000,0.788785);
   gre->SetPointError(9,0,0.010108437);
   gre->SetPoint(10,1500,0.772411);
   gre->SetPointError(10,0,0.0100857);
   */









 /*
   0.25-0.50 PFO sigma
   gre->SetPoint(0,1,16.4653);
   gre->SetPointError(0,0,0.283278);
   gre->SetPoint(1,5,7.2512);
   gre->SetPointError(1,0,0.115721);
   gre->SetPoint(2,10,5.14994);
   gre->SetPointError(2,0,0.0910081);
   gre->SetPoint(3,15,4.34142);
   gre->SetPointError(3,0,0.078066);
   gre->SetPoint(4,30,2.92558);
   gre->SetPointError(4,0,0.0492753);
   gre->SetPoint(5,50,2.34268);
   gre->SetPointError(5,0,0.0412328);
   gre->SetPoint(6,100,1.72443);
   gre->SetPointError(6,0,0.0312565);
   gre->SetPoint(7,200,1.26026);
   gre->SetPointError(7,0,0.0209845);
   gre->SetPoint(8,500,0.970558);
   gre->SetPointError(8,0,0.0152366);
   gre->SetPoint(9,1000,0.823838);
   gre->SetPointError(9,0,0.0137154);
   gre->SetPoint(10,1500,0.797728);
   gre->SetPointError(10,0,0.0138269);


   //0.25-0.50 hit level RMS
   gre->SetPoint(0,1,17.2938);
   gre->SetPointError(0,0,0.247206);
   gre->SetPoint(1,5,6.87257);
   gre->SetPointError(1,0,0.0967104);
   gre->SetPoint(2,10,4.91462);
   gre->SetPointError(2,0,0.0695589);
   gre->SetPoint(3,15,3.90733);
   gre->SetPointError(3,0,0.0551258);
   gre->SetPoint(4,30,2.88322);
   gre->SetPointError(4,0,0.0405146);
   gre->SetPoint(5,50,2.228695);
   gre->SetPointError52,0,0.0320614);
   gre->SetPoint(6,100,1.66163);
   gre->SetPointError(6,0,0.0236366);
   gre->SetPoint(7,200,1.37256);
   gre->SetPointError(7,0,0.0192423);
   gre->SetPoint(8,500,1.09778);
   gre->SetPointError(8,0,0.0154909);
   gre->SetPoint(9,1000,0.997316);
   gre->SetPointError(9,0,0.039844);
   gre->SetPoint(10,1500,1.0031);
   gre->SetPointError(10,0,0.014718);

   //0.25-0.50 hit level sigma
   gre->SetPoint(0,1,14.6227);
   gre->SetPointError(0,0,0.242287);
   gre->SetPoint(1,5,6.65832);
   gre->SetPointError(1,0,0.0969049);
   gre->SetPoint(2,10,4.83659);
   gre->SetPointError(2,0,0.0750541);
   gre->SetPoint(3,15,3.88679);
   gre->SetPointError(3,0,0.0595577);
   gre->SetPoint(4,30,2.78325);
   gre->SetPointError(4,0,0.042218);
   gre->SetPoint(5,50,2.2245);
   gre->SetPointError(5,0,0.0345251);
   gre->SetPoint(6,100,1.62065);
   gre->SetPointError(6,0,0.0247537);
   gre->SetPoint(7,200,1.22783);
   gre->SetPointError(7,0,0.017457);
   gre->SetPoint(8,500,0.936447);
   gre->SetPointError(8,0,0.0140159);
   gre->SetPoint(9,1000,0.807871);
   gre->SetPointError(9,0,0.0118402);
   gre->SetPoint(10,1500,0.777224);
   gre->SetPointError(10,0,0.0113027);
   */



   /*
   0.5-0.75 PFO sigma
   gre->SetPoint(0,1,17.1982);
   gre->SetPointError(0,0,0.286913);
   gre->SetPoint(1,5,7.11654);
   gre->SetPointError(1,0,0.126187);
   gre->SetPoint(2,10,5.41554);
   gre->SetPointError(2,0,0.0906201);
   gre->SetPoint(3,15,4.29473);
   gre->SetPointError(3,0,0.0723595);
   gre->SetPoint(4,30,3.08629);
   gre->SetPointError(4,0,0.0522377);
   gre->SetPoint(5,50,2.0239209);
   gre->SetPointError(5,0,0.0410864);
   gre->SetPoint(6,100,1.81468);
   gre->SetPointError(6,0,0.035585);
   gre->SetPoint(7,200,1.25641);
   gre->SetPointError(7,0,0.0225986);
   gre->SetPoint(8,500,0.923744);
   gre->SetPointError(8,0,0.0132948);
   gre->SetPoint(9,1000,0.732591);
   gre->SetPointError(9,0,0.010258);
   gre->SetPoint(10,1500,0.725751);
   gre->SetPointError(10,0,0.0115379);


   //0.5-0.75 hit level RMS
   gre->SetPoint(0,1,18.8278);
   gre->SetPointError(0,0,0.264891);
   gre->SetPoint(1,5,6.93684);
   gre->SetPointError(1,0,0.0979256);
   gre->SetPoint(2,10,5.13338);
   gre->SetPointError(2,0,0.0721795);
   gre->SetPoint(3,15,4.12133);
   gre->SetPointError(3,0,0.0579379);
   gre->SetPoint(4,30,2.86023);
   gre->SetPointError(4,0,0.0404659);
   gre->SetPoint(5,50,2.33018);
   gre->SetPointError52,0,0.0335703);
   gre->SetPoint(6,100,1.76228);
   gre->SetPointError(6,0,0.0246817);
   gre->SetPoint(7,200,1.3126);
   gre->SetPointError(7,0,0.0185481);
   gre->SetPoint(8,500,1.00049);
   gre->SetPointError(8,0,0.0142261);
   gre->SetPoint(9,1000,0.860736);
   gre->SetPointError(9,0,0.0123265);
   gre->SetPoint(10,1500,0.949787);
   gre->SetPointError(10,0,0.0134536);

   //0.5-0.75 hit level sigma
   gre->SetPoint(0,1,15.1148);
   gre->SetPointError(0,0,0.239315);
   gre->SetPoint(1,5,6.90736);
   gre->SetPointError(1,0,0.126187);
   gre->SetPoint(2,10,5.03977);
   gre->SetPointError(2,0,0.0769003);
   gre->SetPoint(3,15,4.10118);
   gre->SetPointError(3,0,0.0630573);
   gre->SetPoint(4,30,2.8244);
   gre->SetPointError(4,0,0.0430629);
   gre->SetPoint(5,50,2.26938);
   gre->SetPointError(5,0,0.0354298);
   gre->SetPoint(6,100,1.6927);
   gre->SetPointError(6,0,0.026556);
   gre->SetPoint(7,200,1.23372);
   gre->SetPointError(7,0,0.0199157);
   gre->SetPoint(8,500,0.899584);
   gre->SetPointError(8,0,0.0123605);
   gre->SetPoint(9,1000,0.744748);
   gre->SetPointError(9,0,0.0104553);
   gre->SetPoint(10,1500,0.715817);
   gre->SetPointError(10,0,0.0106226);
   */


   //none working for PFO at 1500,1000,500,200,100 -> double peak structure there, at 30, 50 maybe start to move together
   //seems to work at 15
   /*
   0.75-0.87 PFO sigma
   gre->SetPoint(0,1,14.728);
   gre->SetPointError(0,0,0.401521);
   gre->SetPoint(1,5,7.12988);
   gre->SetPointError(1,0,0.197205);
   gre->SetPoint(2,10,5.1438);
   gre->SetPointError(2,0,0.134202);
   gre->SetPoint(3,15,4.51798);
   gre->SetPointError(3,0,0.117806);
   gre->SetPoint(4,30,3.60117);
   gre->SetPointError(4,0,0.0100872);
   gre->SetPoint(5,50,3.21945);
   gre->SetPointError(5,0,0.09467);


   //0.75-0.87 hit level RMS
   gre->SetPoint(0,1,17.8683);
   gre->SetPointError(0,0,0.352603);
   gre->SetPoint(1,5,6.90276);
   gre->SetPointError(1,0,0.145847);
   gre->SetPoint(2,10,5.16694);
   gre->SetPointError(2,0,0.00106225);
   gre->SetPoint(3,15,4.32176);
   gre->SetPointError(3,0,0.0875992);
   gre->SetPoint(4,30,3.18725);
   gre->SetPointError(4,0,0.0643135);
   gre->SetPoint(5,50,2.6055);
   gre->SetPointError(5,0,0.0534075);
   gre->SetPoint(6,100,1.98919);
   gre->SetPointError(6,0,0.0419358);
   gre->SetPoint(7,200,1.79579);
   gre->SetPointError(7,0,0.037187);
   gre->SetPoint(8,500,1.68346);
   gre->SetPointError(8,0,0.0344786);
   gre->SetPoint(9,1000,1.61946);
   gre->SetPointError(9,0,0.0330021);
   gre->SetPoint(10,1500,1.8044);
   gre->SetPointError(10,0,0.0367861);

   //0.75-0.87 hit level sigma
   gre->SetPoint(0,1,14.8995);
   gre->SetPointError(0,0,0.349464);
   gre->SetPoint(1,5,6.83121);
   gre->SetPointError(1,0,0.185628);
   gre->SetPoint(2,10,5.18225);
   gre->SetPointError(2,0,0.12082);
   gre->SetPoint(3,15,4.23033);
   gre->SetPointError(3,0,0.0987889);
   gre->SetPoint(4,30,3.13727);
   gre->SetPointError(4,0,0.0781122);
   gre->SetPoint(5,50,2.5165);
   gre->SetPointError(5,0,0.0587726);
   gre->SetPoint(6,100,1.78123);
   gre->SetPointError(6,0,0.0447898);
   gre->SetPoint(7,200,1.3652);
   gre->SetPointError(7,0,0.041008);
   gre->SetPoint(8,500,0.915652);
   gre->SetPointError(8,0,0.0246656);
   gre->SetPoint(9,1000,0.751177);
   gre->SetPointError(9,0,0.0203032);
   gre->SetPoint(10,1500,0.745737);
   gre->SetPointError(10,0,0.0194221);
   */






  /*
   0.87-0.98 PFO sigma
   gre->SetPoint(0,1,14.3402);
   gre->SetPointError(0,0,0.430262);
   gre->SetPoint(1,5,6.58956);
   gre->SetPointError(1,0,0.181129);
   gre->SetPoint(2,10,4.69535);
   gre->SetPointError(2,0,0.12569);
   gre->SetPoint(3,15,3.71234);
   gre->SetPointError(3,0,0.0976756);
   gre->SetPoint(4,30,2.71595);
   gre->SetPointError(4,0,0.0661738);
   gre->SetPoint(5,50,2.15928);
   gre->SetPointError(5,0,0.0536925);
   gre->SetPoint(6,100,1.59802);
   gre->SetPointError(6,0,0.045252);
   gre->SetPoint(7,200,1.16731);
   gre->SetPointError(7,0,0.03692);
   gre->SetPoint(8,500,0.809776);
   gre->SetPointError(8,0,0.0274741);
   gre->SetPoint(9,1000,0.540026);
   gre->SetPointError(9,0,0.020106);
   gre->SetPoint(10,1500,0.512441);
   gre->SetPointError(10,0,0.0243851);


   //0.87-0.98 hit level RMS
   gre->SetPoint(0,1,22.331);
   gre->SetPointError(0,0,0.479158);
   gre->SetPoint(1,5,6.89553);
   gre->SetPointError(1,0,0.150258);
   gre->SetPoint(2,10,4.88916);
   gre->SetPointError(2,0,0.105738);
   gre->SetPoint(3,15,3.84296);
   gre->SetPointError(3,0,0.0818951);
   gre->SetPoint(4,30,2.80064);
   gre->SetPointError(4,0,0.0603439);
   gre->SetPoint(5,50,2.21815);
   gre->SetPointError52,0,0.0473557);
   gre->SetPoint(6,100,1.61195);
   gre->SetPointError(6,0,0.0340587);
   gre->SetPoint(7,200,1.17336);
   gre->SetPointError(7,0,0.0251421);
   gre->SetPoint(8,500,0.794015);
   gre->SetPointError(8,0,0.0170294);
   gre->SetPoint(9,1000,0.644118);
   gre->SetPointError(9,0,0.0133901);
   gre->SetPoint(10,1500,0.652979);
   gre->SetPointError(10,0,0.0140304);

   //0.87-0.98 hit level sigma
   gre->SetPoint(0,1,14.1549);
   gre->SetPointError(0,0,0.00412641);
   gre->SetPoint(1,5,6.7801);
   gre->SetPointError(1,0,0.182961);
   gre->SetPoint(2,10,4.80903);
   gre->SetPointError(2,0,0.12431);
   gre->SetPoint(3,15,3.73019);
   gre->SetPointError(3,0,0.086425);
   gre->SetPoint(4,30,2.72688);
   gre->SetPointError(4,0,0.063397);
   gre->SetPoint(5,50,2.15597);
   gre->SetPointError(5,0,0.0482383);
   gre->SetPoint(6,100,1.57722);
   gre->SetPointError(6,0,0.0384446);
   gre->SetPoint(7,200,1.13882);
   gre->SetPointError(7,0,0.0277035);
   gre->SetPoint(8,500,0.761246);
   gre->SetPointError(8,0,0.0180269);
   gre->SetPoint(9,1000,0.56314);
   gre->SetPointError(9,0,0.0125595);
   gre->SetPoint(10,1500,0.495647);
   gre->SetPointError(10,0,0.0103928);
   */

 
  /*
  TF1* fit_sigma= new TF1 ("fit_sigma","sqrt([0]*[0]+[1]*[1]/x)",70,1600);
  fit_sigma->SetParameter(0,0.7);
  fit_sigma->SetParameter(1,15.2);
  fit_sigma->SetLineColor(kRed);

  
  Int_t n=11;
  Double_t x[n],y[n];
  x[0]=1;
  y[0]=22.94;
  x[1]=5;
  y[1]=7.175;
  x[2]=10;
  y[2]=4.953;
  x[3]=15;
  y[3]=4.056;
  x[4]=30;
  y[4]=2.914;
  x[5]=50;
  y[5]=2.355;
  x[6]=100;
  y[6]=1.773;
  x[7]=200;
  y[7]=1.366;
  x[8]=500;
  y[8]=1.133;
  x[9]=1000;
  y[9]=1.065;
  x[10]=1500;
  y[10]=0.953;

  TGraph* h_RMS_overview =  new TGraph(n,x,y);
  //h_RMS_overview->Sumw2();
  //h_RMS_overview->SetLineWidth(2);
  h_RMS_overview->SetMarkerStyle(21);
  h_RMS_overview->SetMarkerColor(kBlue);

  TF1* fit_sigma2= new TF1 ("fit_sigma2","sqrt([0]*[0]+[1]*[1]/x)",5,1600);
  fit_sigma2->SetParameter(0,0.0055*0.0055);
  fit_sigma2->SetParameter(1,0);
 

  //h_RMS_overview->Draw("AP");
  //h_RMS_overview->Fit("fit_sigma");

 
  gre->Draw("AP,same");
  gre->Fit("fit_sigma2");

  //h_sigma_overview->Draw("AP,same");
  //h_sigma_overview->Fit("fit_sigma2");

  //for(int i=0;i<h_RMS_overview->GetNbinsX();i++){
  //h_RMS_overview->SetBinError(i,0.01);
  //}

  TLegend* leg_RMS=new TLegend(0.45,0.65,0.89,0.89);
  leg_RMS->SetBorderSize(0);
  leg_RMS->SetFillStyle(0);
  // leg_RMS->AddEntry(h_sigma_overview->DrawCopy("P,same"),"#sigma(E_{rel})");
  //leg_RMS->AddEntry(h_RMS_overview->DrawCopy("P,same"),"RMS (E_{rel})");
  //leg_RMS->Draw();

  TCanvas* can_hits_barrel_2D_z_r=setUpperCanvas("can_hits_barrel_2D_z_r");
  can_hits_barrel_2D_z_r->cd();
  h_2D_hit_z_r_barrel->Draw("colz");

  TCanvas* can_hits_endcap_2D_r_z=setUpperCanvas("can_hits_endcap_2D_r_z");
  can_hits_endcap_2D_r_z->cd();
  h_2D_hit_r_z_endcap->Draw("colz");

  TCanvas* can_hitsE_barrel_2D_z_r=setUpperCanvas("can_hitsE_barrel_2D_z_r");
  can_hitsE_barrel_2D_z_r->cd();
  h_2D_hitE_z_r_barrel->Draw("colz");

  TCanvas* can_hitsE_endcap_2D_r_z=setUpperCanvas("can_hitsE_endcap_2D_r_z");
  can_hitsE_endcap_2D_r_z->cd();
  h_2D_hitE_r_z_endcap->Draw("colz");
  */
/*
  TCanvas* can_layerEB_E=setUpperCanvas("can_layerEB_E");
  can_layerEB_E->cd();
  
  h_layerEE_E->SetLineColor(kRed);
  //h_layerEE_E->Scale(h_layerEB_E->Integral()/h_layerEE_E->Integral());
  h_layerEE_E->DrawCopy("h,same");

  TLegend* leg_E_layers=new TLegend(0.45,0.65,0.89,0.89,"E_{photon}=10 GeV");
  leg_E_layers->SetBorderSize(0);
  leg_E_layers->SetFillStyle(0);
  leg_E_layers->AddEntry(h_layerEB_E->DrawCopy("h"),"ECAL barrel");
  leg_E_layers->AddEntry(h_layerEE_E->DrawCopy("h,same"),"ECAL Endcap");
  leg_E_layers->Draw();

  TCanvas* can_ratio_layerEB_E=setRatioCanvas("can_layer_EB_EE_ratio");
  can_ratio_layerEB_E->cd();
  TH1F* h_layerEB_EE_ratio=(TH1F*)h_layerEE_E->Clone("h_ratio_Layer_EE_vs_EB");
  h_layerEB_EE_ratio->Divide(h_layerEB_E);
  h_layerEB_EE_ratio->SetLineWidth(2);
  h_layerEB_EE_ratio->SetMaximum(0.95);
  h_layerEB_EE_ratio->SetMinimum(0.80);
  h_layerEB_EE_ratio->GetYaxis()->SetTitle("ratio Endcap/Barrel");
  h_layerEB_EE_ratio->Draw("hist");

  TCanvas* can_layerEB_E_ph_vs_n=setUpperCanvas("can_layerEB_E_ph_vs_n");
  can_layerEB_E_ph_vs_n->cd();
  h_layerEB_E_n->Scale(h_layerEB_E_ph->Integral()/h_layerEB_E_n->Integral());

  TLegend* leg_E_layersEB_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,"E_{photon}=1 GeV");
  leg_E_layersEB_ph_vs_n->SetBorderSize(0);
  leg_E_layersEB_ph_vs_n->SetFillStyle(0);
  leg_E_layersEB_ph_vs_n->AddEntry(h_layerEB_E_ph->DrawCopy("h"),"ECAL barrel, photon");
  leg_E_layersEB_ph_vs_n->AddEntry(h_layerEB_E_n->DrawCopy("h,same"),"ECAL barrel, neutron");
  leg_E_layersEB_ph_vs_n->Draw();

  TCanvas* can_layerEE_E_ph_vs_n=setUpperCanvas("can_layerEE_E_ph_vs_n");
  can_layerEE_E_ph_vs_n->cd();
  h_layerEE_E_n->Scale(h_layerEE_E_ph->Integral()/h_layerEE_E_n->Integral());

  TLegend* leg_E_layersEE_ph_vs_n=new TLegend(0.45,0.65,0.89,0.89,"E_{photon}=1 GeV");
  leg_E_layersEE_ph_vs_n->SetBorderSize(0);
  leg_E_layersEE_ph_vs_n->SetFillStyle(0);
  leg_E_layersEE_ph_vs_n->AddEntry(h_layerEE_E_ph->DrawCopy("h"),"ECAL endcap, photon");
  leg_E_layersEE_ph_vs_n->AddEntry(h_layerEE_E_n->DrawCopy("h,same"),"ECAL endcap, neutron");
  leg_E_layersEE_ph_vs_n->Draw();
*/
  return 1;

}

