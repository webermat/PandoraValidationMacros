#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TList.h"

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

TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fillHistograms(TFile* file_training, std::vector<TH1F*> t_hist_vec){

  TTree* tree_train= (TTree*)file_training->Get("showerData");

  vector<float> *true_Energy=0;
  vector<float> *true_Px=0;
  vector<float> *true_Py=0;
  vector<float> *true_Pz=0;
  vector<int> *true_PDGID=0;
  vector<int> *true_GenStatus=0;

  vector<float> *reco_Energy=0;
  vector<float> *reco_Px=0;
  vector<float> *reco_Py=0;
  vector<float> *reco_Pz=0;
  vector<int> *reco_PDGID=0;

  vector<float> *cluster_Energy=0;

  tree_train->SetBranchAddress("true_Energy", &true_Energy);
  tree_train->SetBranchAddress("true_Px", &true_Px);
  tree_train->SetBranchAddress("true_Py", &true_Py);
  tree_train->SetBranchAddress("true_Pz", &true_Pz);
  tree_train->SetBranchAddress("true_PDGID", &true_PDGID);
  tree_train->SetBranchAddress("true_GenStatus", &true_GenStatus);

  tree_train->SetBranchAddress("reco_Energy", &reco_Energy);
  tree_train->SetBranchAddress("reco_Px", &reco_Px);
  tree_train->SetBranchAddress("reco_Py", &reco_Py);
  tree_train->SetBranchAddress("reco_Pz", &reco_Pz);
  tree_train->SetBranchAddress("reco_PDGID", &reco_PDGID);

  //tree_train->SetBranchAddress("cluster_energy", &cluster_Energy);

  for (unsigned int i_entry=0;i_entry<tree_train->GetEntries();i_entry++){
    tree_train->GetEntry(i_entry);
    bool pass_hadron=false;
    float true_energy=-1;
    for(unsigned int i_true=0;i_true<true_Energy->size();i_true++){
      if((*true_GenStatus)[i_true]==1){
	true_energy=((*true_Energy)[i_true]);
      }
    }
    float E_cluster_reco_sum=0;
    /*
    for(unsigned int i_clu=0;i_clu<cluster_Energy->size();i_clu++){
      E_cluster_reco_sum+=(*reco_Energy)[i_clu];
    }
    */

    float E_PFO_reco_sum=0;
    float E_PFO_reco_lead=-1;
    int PDF_lead_part=-1;
    for(unsigned int i_PFO=0;i_PFO<reco_Energy->size();i_PFO++){
      E_PFO_reco_sum+=(*reco_Energy)[i_PFO];
      if((*reco_Energy)[i_PFO]>E_PFO_reco_lead){
	E_PFO_reco_lead=(*reco_Energy)[i_PFO];
	PDF_lead_part=(*reco_PDGID)[i_PFO];
      }
    }
    if(PDF_lead_part==2112){
      if(reco_Energy->size()!=0){
	t_hist_vec[0]->Fill(E_PFO_reco_sum/true_energy);
	t_hist_vec[3]->Fill(E_PFO_reco_lead/true_energy);
	t_hist_vec[6]->Fill(E_cluster_reco_sum/true_energy);
      }
      if(reco_Energy->size()==1){
	t_hist_vec[1]->Fill(E_PFO_reco_sum/true_energy);
	t_hist_vec[4]->Fill(E_PFO_reco_lead/true_energy);
	t_hist_vec[7]->Fill(E_cluster_reco_sum/true_energy);
      }else if (reco_Energy->size()>1){
	//std::cout<<i_entry<<" "<<PDF_lead_part<<" 0/1 types "<<(*reco_PDGID)[0]<<"/"<<(*reco_PDGID)[1]<<" "<<(*reco_Energy)[0]<<" "<<(*reco_Energy)[1] <<" "<<E_PFO_reco_lead<< " "<<reco_Energy->size()<<std::endl;
	t_hist_vec[2]->Fill(E_PFO_reco_sum/true_energy);
	t_hist_vec[5]->Fill(E_PFO_reco_lead/true_energy);
	t_hist_vec[8]->Fill(E_cluster_reco_sum/true_energy);
      }
    }
  }
}

void plotSWCTuning_Summary(){

 
  CLICdpStyle();

    //2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
  /*

  //TFile* file_SWC_def=TFile::Open("/Users/matthiasweber/rootfiles180106/swcComparisonFiles/allK0L_CLIC_n_K0L_newTune/pionStudy_K0L10_ILC171109_CT_FitFW_wRefit_CLIC_o3_v13_SWC_CLIC_n_K0L_newTune.root");
  //TFile* file_SWC_def=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/171221_gcc7/pionStudy_K0L1000_ILC171221_gcc7_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC.root");
  TFile* file_SWC_def=TFile::Open("/eos/user/w/weberma2/data/swcComparisonFiles/171109/pionStudy_K0L1500_ILC171109_CT_FitFW_wRefit_CLIC_o3_v13_SWC_CLIC_n_K0L_newTune.root");
  TFile* file_SWC_CLIC_n=TFile::Open("/eos/user/w/weberma2/data/swcComparisonFiles/171109/pionStudy_K0L1500_ILC171109_CT_FitFW_wRefit_CLIC_o3_v13_SWC_CLIC_n_K0L_newTune.root");
  TFile* file_SWC_CLIC_K0L=TFile::Open("/eos/user/w/weberma2/data/swcComparisonFiles/171109/pionStudy_K0L1500_ILC171109_CT_FitFW_wRefit_CLIC_o3_v13_SWC_CLIC_n_K0L_newTune.root");
  TFile* file_SWC_CLIC_n_K0L=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/171221_gcc7/pionStudy_K0L1500_ILC171221_gcc7_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC.root");
  //TFile* file_SWC_CLIC_n_K0L=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/171221_gcc7/pionStudy_K0L1000_ILC171221_gcc7_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC.root");



  int n_bins100=100;
  double lim_energy_low=0.75;
  double lim_energy_high=1.25;

  TH1F* h_SWC_def_E_sum=new TH1F("h_SWC_def_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_sum->Sumw2();
  h_SWC_def_E_sum->SetLineWidth(2);
  h_SWC_def_E_sum->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_E_sum=new TH1F("h_SWC_CLIC_n_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_sum->Sumw2();
  h_SWC_CLIC_n_E_sum->SetLineWidth(2);
  h_SWC_CLIC_n_E_sum->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_sum->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_sum=new TH1F("h_SWC_CLIC_K0L_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_sum->Sumw2();
  h_SWC_CLIC_K0L_E_sum->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_sum->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_sum->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_sum=new TH1F("h_SWC_CLIC_n_K0L_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_sum->Sumw2();
  h_SWC_CLIC_n_K0L_E_sum->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_sum->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_sum->SetLineColor(kRed);

  TH1F* h_SWC_def_E_sum_one_Part=new TH1F("h_SWC_def_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_sum_one_Part->Sumw2();
  h_SWC_def_E_sum_one_Part->SetLineWidth(2);
  h_SWC_def_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_E_sum_one_Part=new TH1F("h_SWC_CLIC_n_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_n_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_sum_one_Part->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_sum_one_Part=new TH1F("h_SWC_CLIC_K0L_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_K0L_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_sum_one_Part->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_sum_one_Part=new TH1F("h_SWC_CLIC_n_K0L_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_n_K0L_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_sum_one_Part->SetLineColor(kRed);


  TH1F* h_SWC_def_E_sum_two_PartPlus=new TH1F("h_SWC_def_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_sum_two_PartPlus->Sumw2();
  h_SWC_def_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_def_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_n_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_sum_two_PartPlus->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_K0L_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_K0L_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_sum_two_PartPlus->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_n_K0L_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_K0L_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_sum_two_PartPlus->SetLineColor(kRed);

  TH1F* h_SWC_def_E_PFO_0=new TH1F("h_SWC_def_E_PFO_0","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_PFO_0->Sumw2();
  h_SWC_def_E_PFO_0->SetLineWidth(2);
  h_SWC_def_E_PFO_0->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_E_PFO_0=new TH1F("h_SWC_CLIC_n_E_PFO_0","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_PFO_0->Sumw2();
  h_SWC_CLIC_n_E_PFO_0->SetLineWidth(2);
  h_SWC_CLIC_n_E_PFO_0->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_PFO_0->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_PFO_0=new TH1F("h_SWC_CLIC_K0L_E_PFO_0","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_PFO_0->Sumw2();
  h_SWC_CLIC_K0L_E_PFO_0->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_PFO_0->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_PFO_0->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_PFO_0=new TH1F("h_SWC_CLIC_n_K0L_E_PFO_0","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_PFO_0->Sumw2();
  h_SWC_CLIC_n_K0L_E_PFO_0->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_PFO_0->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_PFO_0->SetLineColor(kRed);


  TH1F* h_SWC_def_E_PFO_0_one_Part=new TH1F("h_SWC_def_E_PFO_0_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_PFO_0_one_Part->Sumw2();
  h_SWC_def_E_PFO_0_one_Part->SetLineWidth(2);
  h_SWC_def_E_PFO_0_one_Part->GetXaxis()->SetTitle("E_{0}^{PFO/E_{true}^{K0L}}");

  TH1F* h_SWC_CLIC_n_E_PFO_0_one_Part=new TH1F("h_SWC_CLIC_n_E_PFO_0_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_PFO_0_one_Part->Sumw2();
  h_SWC_CLIC_n_E_PFO_0_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_E_PFO_0_one_Part->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_PFO_0_one_Part->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_PFO_0_one_Part=new TH1F("h_SWC_CLIC_K0L_E_PFO_0_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_PFO_0_one_Part->Sumw2();
  h_SWC_CLIC_K0L_E_PFO_0_one_Part->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_PFO_0_one_Part->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_PFO_0_one_Part->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_PFO_0_one_Part=new TH1F("h_SWC_CLIC_n_K0L_E_PFO_0_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_PFO_0_one_Part->Sumw2();
  h_SWC_CLIC_n_K0L_E_PFO_0_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_PFO_0_one_Part->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_PFO_0_one_Part->SetLineColor(kRed);

  TH1F* h_SWC_def_E_PFO_0_two_PartPlus=new TH1F("h_SWC_def_E_PFO_0_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_E_PFO_0_two_PartPlus->Sumw2();
  h_SWC_def_E_PFO_0_two_PartPlus->SetLineWidth(2);
  h_SWC_def_E_PFO_0_two_PartPlus->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_E_PFO_0_two_PartPlus=new TH1F("h_SWC_CLIC_n_E_PFO_0_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_E_PFO_0_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_E_PFO_0_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_E_PFO_0_two_PartPlus->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_E_PFO_0_two_PartPlus->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus=new TH1F("h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus->Sumw2();
  h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus=new TH1F("h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus->GetXaxis()->SetTitle("E_{0}^{PFO}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus->SetLineColor(kRed);

  TH1F* h_SWC_def_clust_E_sum=new TH1F("h_SWC_def_clust_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_clust_E_sum->Sumw2();
  h_SWC_def_clust_E_sum->SetLineWidth(2);
  h_SWC_def_clust_E_sum->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_clust_E_sum=new TH1F("h_SWC_CLIC_n_clust_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_clust_E_sum->Sumw2();
  h_SWC_CLIC_n_clust_E_sum->SetLineWidth(2);
  h_SWC_CLIC_n_clust_E_sum->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_clust_E_sum->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_clust_E_sum=new TH1F("h_SWC_CLIC_K0L_clust_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_clust_E_sum->Sumw2();
  h_SWC_CLIC_K0L_clust_E_sum->SetLineWidth(2);
  h_SWC_CLIC_K0L_clust_E_sum->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_clust_E_sum->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_clust_E_sum=new TH1F("h_SWC_CLIC_n_K0L_clust_E_sum","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_clust_E_sum->Sumw2();
  h_SWC_CLIC_n_K0L_clust_E_sum->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_clust_E_sum->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_clust_E_sum->SetLineColor(kRed);


  TH1F* h_SWC_def_clust_E_sum_one_Part=new TH1F("h_SWC_def_clust_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_clust_E_sum_one_Part->Sumw2();
  h_SWC_def_clust_E_sum_one_Part->SetLineWidth(2);
  h_SWC_def_clust_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_clust_E_sum_one_Part=new TH1F("h_SWC_CLIC_n_clust_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_clust_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_n_clust_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_clust_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_clust_E_sum_one_Part->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_clust_E_sum_one_Part=new TH1F("h_SWC_CLIC_K0L_clust_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_clust_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_K0L_clust_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_K0L_clust_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_clust_E_sum_one_Part->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_clust_E_sum_one_Part=new TH1F("h_SWC_CLIC_n_K0L_clust_E_sum_one_Part","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_clust_E_sum_one_Part->Sumw2();
  h_SWC_CLIC_n_K0L_clust_E_sum_one_Part->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_clust_E_sum_one_Part->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_clust_E_sum_one_Part->SetLineColor(kRed);

  TH1F* h_SWC_def_clust_E_sum_two_PartPlus=new TH1F("h_SWC_def_clust_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_def_clust_E_sum_two_PartPlus->Sumw2();
  h_SWC_def_clust_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_def_clust_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");

  TH1F* h_SWC_CLIC_n_clust_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_n_clust_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_clust_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_clust_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_clust_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_clust_E_sum_two_PartPlus->SetLineColor(kBlue);

  TH1F* h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus->SetLineColor(kMagenta);

  TH1F* h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus=new TH1F("h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus","", n_bins100, lim_energy_low,lim_energy_high);
  h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus->Sumw2();
  h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus->SetLineWidth(2);
  h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus->GetXaxis()->SetTitle("E_{sum}^{cluster}/E_{true}^{K0L}");
  h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus->SetLineColor(kRed);

  std::vector<TH1F*> hist_vector_SWC_def;
  hist_vector_SWC_def.push_back(h_SWC_def_E_sum);
  hist_vector_SWC_def.push_back(h_SWC_def_E_sum_one_Part);
  hist_vector_SWC_def.push_back(h_SWC_def_E_sum_two_PartPlus);
  hist_vector_SWC_def.push_back(h_SWC_def_E_PFO_0);
  hist_vector_SWC_def.push_back(h_SWC_def_E_PFO_0_one_Part);
  hist_vector_SWC_def.push_back(h_SWC_def_E_PFO_0_two_PartPlus);
  hist_vector_SWC_def.push_back(h_SWC_def_clust_E_sum);
  hist_vector_SWC_def.push_back(h_SWC_def_clust_E_sum_one_Part);
  hist_vector_SWC_def.push_back(h_SWC_def_clust_E_sum_two_PartPlus);
  fillHistograms(file_SWC_def,hist_vector_SWC_def);

   std::cout<<"n"<<std::endl;
  std::vector<TH1F*> hist_vector_SWC_CLIC_n;
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_sum);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_sum_one_Part);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_sum_two_PartPlus);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_PFO_0);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_PFO_0_one_Part);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_E_PFO_0_two_PartPlus);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_clust_E_sum);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_clust_E_sum_one_Part);
  hist_vector_SWC_CLIC_n.push_back(h_SWC_CLIC_n_clust_E_sum_two_PartPlus);
  fillHistograms(file_SWC_CLIC_n,hist_vector_SWC_CLIC_n);

  std::cout<<"K0L"<<std::endl;

  std::vector<TH1F*> hist_vector_SWC_CLIC_K0L;
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_sum);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_sum_one_Part);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_sum_two_PartPlus);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_PFO_0);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_PFO_0_one_Part);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_clust_E_sum);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_clust_E_sum_one_Part);
  hist_vector_SWC_CLIC_K0L.push_back(h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus);
  fillHistograms(file_SWC_CLIC_K0L,hist_vector_SWC_CLIC_K0L);

  std::cout<<"n_K0L"<<std::endl;

  std::vector<TH1F*> hist_vector_SWC_CLIC_n_K0L;
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_sum);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_sum_one_Part);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_sum_two_PartPlus);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_PFO_0);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_PFO_0_one_Part);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_clust_E_sum);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_clust_E_sum_one_Part);
  hist_vector_SWC_CLIC_n_K0L.push_back(h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus);
  fillHistograms(file_SWC_CLIC_n_K0L,hist_vector_SWC_CLIC_n_K0L);


  for (unsigned int i=0;i<hist_vector_SWC_def.size();i++){
    hist_vector_SWC_CLIC_n_K0L[i]->Scale(hist_vector_SWC_def[i]->Integral()/hist_vector_SWC_CLIC_n_K0L[i]->Integral());
    hist_vector_SWC_CLIC_n[i]->Scale(hist_vector_SWC_def[i]->Integral()/hist_vector_SWC_CLIC_n[i]->Integral());
    hist_vector_SWC_CLIC_K0L[i]->Scale(hist_vector_SWC_def[i]->Integral()/hist_vector_SWC_CLIC_K0L[i]->Integral());
    hist_vector_SWC_def[i]->SetMaximum(hist_vector_SWC_def[i]->GetMaximum()*1.45);
    std::cout<<"vector "<<i<<" K0L, 1500 GeV "<<std::endl;
    std::cout<<"ILD/n/K/n+K mean "<<hist_vector_SWC_def[i]->GetMean()<<"/"<<hist_vector_SWC_CLIC_n[i]->GetMean()<<"/"<<hist_vector_SWC_CLIC_K0L[i]->GetMean()<<"/"<<hist_vector_SWC_CLIC_n_K0L[i]->GetMean()<<std::endl;
    std::cout<<"ILD/n/K/n+K mean error "<<hist_vector_SWC_def[i]->GetMeanError()<<"/"<<hist_vector_SWC_CLIC_n[i]->GetMeanError()<<"/"<<hist_vector_SWC_CLIC_K0L[i]->GetMeanError()<<"/"<<hist_vector_SWC_CLIC_n_K0L[i]->GetMeanError()<<"/"<<hist_vector_SWC_CLIC_n[i]->GetMeanError()<<std::endl;
    std::cout<<"ILD/n/K/n+K rms "<<hist_vector_SWC_def[i]->GetRMS()<<"/"<<hist_vector_SWC_CLIC_n[i]->GetRMS()<<"/"<<hist_vector_SWC_CLIC_K0L[i]->GetRMS()<<"/"<<hist_vector_SWC_CLIC_n_K0L[i]->GetRMS()<<std::endl;
    std::cout<<"ILD/n/K/n+K rms error "<<hist_vector_SWC_def[i]->GetRMSError()<<"/"<<hist_vector_SWC_CLIC_n[i]->GetRMSError()<<"/"<<hist_vector_SWC_CLIC_n_K0L[i]->GetRMSError()<<"/"<<hist_vector_SWC_CLIC_K0L[i]->GetRMSError()<<std::endl;
  }

  TCanvas* can_E_PFO_sum=setUpperCanvas("can_E_PFO_sum");
  can_E_PFO_sum->cd();

  TLegend* leg_E_PFO_sum=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_sum->SetBorderSize(0);
  leg_E_PFO_sum->SetFillStyle(0);
  leg_E_PFO_sum->SetHeader("K0L, 1500 GeV, all events, K0L ID","h");
  leg_E_PFO_sum->AddEntry(h_SWC_def_E_sum->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_sum->AddEntry(h_SWC_CLIC_n_E_sum->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_sum->AddEntry(h_SWC_CLIC_K0L_E_sum->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_sum->AddEntry(h_SWC_CLIC_n_K0L_E_sum->DrawCopy("hist,e,same"),"CLIC n+K0L new weights");
  leg_E_PFO_sum->Draw();

  TCanvas* can_E_PFO_sum_one_Part=setUpperCanvas("can_E_PFO_sum_one_Part");
  can_E_PFO_sum_one_Part->cd();

  TLegend* leg_E_PFO_sum_one_Part=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_sum_one_Part->SetBorderSize(0);
  leg_E_PFO_sum_one_Part->SetFillStyle(0);
  leg_E_PFO_sum_one_Part->SetHeader("K0L, 1500 GeV, one reco PFO","h");
  leg_E_PFO_sum_one_Part->AddEntry(h_SWC_def_E_sum_one_Part->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_sum_one_Part->AddEntry(h_SWC_CLIC_n_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_sum_one_Part->AddEntry(h_SWC_CLIC_K0L_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_sum_one_Part->AddEntry(h_SWC_CLIC_n_K0L_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_E_PFO_sum_one_Part->Draw();

  TCanvas* can_E_PFO_sum_two_PartPlus=setUpperCanvas("can_E_PFO_sum_two_PartPlus");
  can_E_PFO_sum_two_PartPlus->cd();

  TLegend* leg_E_PFO_sum_two_PartPlus=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_sum_two_PartPlus->SetBorderSize(0);
  leg_E_PFO_sum_two_PartPlus->SetFillStyle(0);
  leg_E_PFO_sum_two_PartPlus->SetHeader("K0L, 1500 GeV, two or reco PFOs","h");
  leg_E_PFO_sum_two_PartPlus->AddEntry(h_SWC_def_E_sum_two_PartPlus->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_sum_two_PartPlus->AddEntry(h_SWC_CLIC_n_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_sum_two_PartPlus->AddEntry(h_SWC_CLIC_K0L_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_sum_two_PartPlus->AddEntry(h_SWC_CLIC_n_K0L_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_E_PFO_sum_two_PartPlus->Draw();


  TCanvas* can_E_PFO_0=setUpperCanvas("can_E_PFO_0");
  can_E_PFO_0->cd();

  TLegend* leg_E_PFO_0=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_0->SetBorderSize(0);
  leg_E_PFO_0->SetFillStyle(0);
  leg_E_PFO_0->SetHeader("K0L, 1500 GeV, all events","h");
  leg_E_PFO_0->AddEntry(h_SWC_def_E_PFO_0->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_0->AddEntry(h_SWC_CLIC_n_E_PFO_0->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_0->AddEntry(h_SWC_CLIC_K0L_E_PFO_0->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_0->AddEntry(h_SWC_CLIC_n_K0L_E_PFO_0->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_E_PFO_0->Draw();

  TCanvas* can_E_PFO_0_one_Part=setUpperCanvas("can_E_PFO_0_one_Part");
  can_E_PFO_0_one_Part->cd();

  TLegend* leg_E_PFO_0_one_Part=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_0_one_Part->SetBorderSize(0);
  leg_E_PFO_0_one_Part->SetFillStyle(0);
  leg_E_PFO_0_one_Part->SetHeader("K0L, 1500 GeV, one reco PFO","h");
  leg_E_PFO_0_one_Part->AddEntry(h_SWC_def_E_PFO_0_one_Part->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_0_one_Part->AddEntry(h_SWC_CLIC_n_E_PFO_0_one_Part->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_0_one_Part->AddEntry(h_SWC_CLIC_K0L_E_PFO_0_one_Part->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_0_one_Part->AddEntry(h_SWC_CLIC_n_K0L_E_PFO_0_one_Part->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_E_PFO_0_one_Part->Draw();

  TCanvas* can_E_PFO_0_two_PartPlus=setUpperCanvas("can_E_PFO_0_two_PartPlus");
  can_E_PFO_0_two_PartPlus->cd();

  TLegend* leg_E_PFO_0_two_PartPlus=new TLegend(0.18,0.70,0.60,0.89);
  leg_E_PFO_0_two_PartPlus->SetBorderSize(0);
  leg_E_PFO_0_two_PartPlus->SetFillStyle(0);
  leg_E_PFO_0_two_PartPlus->SetHeader("K0L, 1500 GeV, two or reco PFOs","h");
  leg_E_PFO_0_two_PartPlus->AddEntry(h_SWC_def_E_PFO_0_two_PartPlus->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_E_PFO_0_two_PartPlus->AddEntry(h_SWC_CLIC_n_E_PFO_0_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_E_PFO_0_two_PartPlus->AddEntry(h_SWC_CLIC_K0L_E_PFO_0_two_PartPlus->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_E_PFO_0_two_PartPlus->AddEntry(h_SWC_CLIC_n_K0L_E_PFO_0_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_E_PFO_0_two_PartPlus->Draw();

  TCanvas* can_clust_E_sum=setUpperCanvas("can_clust_E_sum");
  can_clust_E_sum->cd();

  TLegend* leg_clust_E_sum=new TLegend(0.18,0.70,0.60,0.89);
  leg_clust_E_sum->SetBorderSize(0);
  leg_clust_E_sum->SetFillStyle(0);
  leg_clust_E_sum->SetHeader("K0L, 1500 GeV, all events, K0L ID","h");
  leg_clust_E_sum->AddEntry(h_SWC_def_clust_E_sum->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_clust_E_sum->AddEntry(h_SWC_CLIC_n_clust_E_sum->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_clust_E_sum->AddEntry(h_SWC_CLIC_K0L_clust_E_sum->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_clust_E_sum->AddEntry(h_SWC_CLIC_n_K0L_clust_E_sum->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_clust_E_sum->Draw();

  TCanvas* can_clust_E_sum_one_Part=setUpperCanvas("can_clust_E_sum_one_Part");
  can_clust_E_sum_one_Part->cd();

  TLegend* leg_clust_E_sum_one_Part=new TLegend(0.18,0.70,0.60,0.89);
  leg_clust_E_sum_one_Part->SetBorderSize(0);
  leg_clust_E_sum_one_Part->SetFillStyle(0);
  leg_clust_E_sum_one_Part->SetHeader("K0L, 1500 GeV, one reco PFO","h");
  leg_clust_E_sum_one_Part->AddEntry(h_SWC_def_clust_E_sum_one_Part->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_clust_E_sum_one_Part->AddEntry(h_SWC_CLIC_n_clust_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_clust_E_sum_one_Part->AddEntry(h_SWC_CLIC_K0L_clust_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_clust_E_sum_one_Part->AddEntry(h_SWC_CLIC_n_K0L_clust_E_sum_one_Part->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_clust_E_sum_one_Part->Draw();

  TCanvas* can_clust_E_sum_two_PartPlus=setUpperCanvas("can_clust_E_sum_two_PartPlus");
  can_clust_E_sum_two_PartPlus->cd();

  TLegend* leg_clust_E_sum_two_PartPlus=new TLegend(0.18,0.70,0.60,0.89);
  leg_clust_E_sum_two_PartPlus->SetBorderSize(0);
  leg_clust_E_sum_two_PartPlus->SetFillStyle(0);
  leg_clust_E_sum_two_PartPlus->SetHeader("K0L, 1500 GeV, two or reco PFOs","h");
  leg_clust_E_sum_two_PartPlus->AddEntry(h_SWC_def_clust_E_sum_two_PartPlus->DrawCopy("hist,e"),"def weights, with 100 GeV cutoff");
  leg_clust_E_sum_two_PartPlus->AddEntry(h_SWC_CLIC_n_clust_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n weights");
  leg_clust_E_sum_two_PartPlus->AddEntry(h_SWC_CLIC_K0L_clust_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC K0L weights");
  leg_clust_E_sum_two_PartPlus->AddEntry(h_SWC_CLIC_n_K0L_clust_E_sum_two_PartPlus->DrawCopy("hist,e,same"),"CLIC n+K0L weights");
  leg_clust_E_sum_two_PartPlus->Draw();
    
    
    TTree* tree_SWC_CLIC_n_K0L= (TTree*)file_SWC_CLIC_n_K0L->Get("showerData");
    
    vector<float> *true_Energy_SWC_CLIC_n_K0L=0;
    vector<float> *true_Px_SWC_CLIC_n_K0L=0;
    vector<float> *true_Py_SWC_CLIC_n_K0L=0;
    vector<float> *true_Pz_SWC_CLIC_n_K0L=0;
    vector<float> *true_CosTheta_SWC_CLIC_n_K0L=0;
    vector<int> *true_PDGID_SWC_CLIC_n_K0L=0;
    vector<int> *true_GenStatus_SWC_CLIC_n_K0L=0;
    
    
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
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Energy", &reco_Energy_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Px", &reco_Px_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Py", &reco_Py_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_Pz", &reco_Pz_SWC_CLIC_n_K0L);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_PDGID", &reco_PDGID_SWC_CLIC_n_K0L);
    
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_EB", &reco_E_EB);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_EE", &reco_E_EE);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_HB", &reco_E_HB);
    tree_SWC_CLIC_n_K0L->SetBranchAddress("reco_E_HE", &reco_E_HE);
    
    //tree_SWC_CLIC_n_K0L->SetBranchAddress("cluster_energy", &cluster_Energy_SWC_CLIC_n_K0L);

    float lim_energy_low_rel=0.00;
    float lim_energy_high_rel=1.40;
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");

    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");

    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
    TH1F* h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94 =new TH1F("h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94","", n_bins100, lim_energy_low_rel,lim_energy_high_rel);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->Sumw2();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->SetLineWidth(2);
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->GetXaxis()->SetTitle("E^{reco}_{neut}/E_{true}^{MC}");
    
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
        for(unsigned int i_true=0;i_true<true_Energy_SWC_CLIC_n_K0L->size();i_true++){
            if((*true_GenStatus_SWC_CLIC_n_K0L)[i_true]==1 && abs((*true_PDGID_SWC_CLIC_n_K0L)[i_true])==130){
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
            if((*reco_PDGID_SWC_CLIC_n_K0L)[i_reco]==2112){//check for neutral hadrons ONLY
              TLorentzVector temp(0,0,0,0);
              temp.SetPxPyPzE((*reco_Px_SWC_CLIC_n_K0L)[i_reco],(*reco_Py_SWC_CLIC_n_K0L)[i_reco],(*reco_Pz_SWC_CLIC_n_K0L)[i_reco],(*reco_Energy_SWC_CLIC_n_K0L)[i_reco]);
              double cosAngcheck=((HadTrue.Px()*temp.Px()+HadTrue.Py()*temp.Py()+HadTrue.Pz()*temp.Pz())/(HadTrue.P()*temp.P()));
              if(cosAngcheck>0.99){//allow 8 % offset for neutral hadron candidate)
                  if(temp.Energy()>E_had_candidate_max){
                    E_had_candidate_max=temp.Energy();
                    HadReco=temp;
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
            if(fabs(true_CosTheta)>=0.80 && fabs(true_CosTheta)<0.94){
                h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->Fill(E_had_candidate_max/true_energy_hadron);
            }else if (fabs(true_CosTheta)>=0.65 && fabs(true_CosTheta)<0.80){
                h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->Fill(E_had_candidate_max/true_energy_hadron);
            }else{
	      if (fabs(true_CosTheta)>=0.30 && fabs(true_CosTheta)<0.65){
                h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->Fill(E_had_candidate_max/true_energy_hadron);
	      }else{
		h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->Fill(E_had_candidate_max/true_energy_hadron);
	      }
	      if (fabs(true_CosTheta)<0.60){
		h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Fill(E_had_candidate_max/true_energy_hadron);
	      }
	      if (fabs(true_CosTheta)<0.65){
		h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->Fill(E_had_candidate_max/true_energy_hadron);
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
    
    TCanvas* can_res_rel_had_0_30=setUpperCanvas("can_h_res_rel_had_0_30");
    can_res_rel_had_0_30->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->Draw("P,e,same");
    TF1* fit_RECO_0_30 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30->GetFunction("gaus");

    TCanvas* can_res_rel_had_0_30_to_0_65=setUpperCanvas("can_h_res_rel_had_0_30_to_0_65");
    can_res_rel_had_0_30_to_0_65->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->Draw("P,e,same");
    TF1* fit_RECO_0_30_to_0_65 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_30_to_0_65->GetFunction("gaus");

    std::cout<<"gre_gcc62->SetPoint(0,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_30->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc62->SetPointError(0,0,"<<100.*fit_RECO_0_30->GetParError(2)<<");"<<std::endl;
    std::cout<<"gre_gcc62->SetPoint(0,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_30_to_0_65->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc62->SetPointError(0,0,"<<100.*fit_RECO_0_30_to_0_65->GetParError(2)<<");"<<std::endl;

    std::cout<<"gre_gcc62_mean->SetPoint(0,"<<true_energy_fillplot<<","<<fit_RECO_0_30->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc62_mean->SetPointError(0,0,"<<fit_RECO_0_30->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_gcc62_mean->SetPoint(0,"<<true_energy_fillplot<<","<<fit_RECO_0_30_to_0_65->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc62_mean->SetPointError(0,0,"<<fit_RECO_0_30_to_0_65->GetParError(1)<<");"<<std::endl;
 


    
    TCanvas* can_res_rel_had_0_60=setUpperCanvas("can_h_res_rel_had_0_60");
    can_res_rel_had_0_60->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->Draw("P,e,same");
    TF1* fit_RECO_0_60 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_60->GetFunction("gaus");
    
    std::cout<<"gre_gcc7->SetPoint(17,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_60->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc7->SetPointError(17,0,"<<100.*fit_RECO_0_60->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_65=setUpperCanvas("can_h_res_rel_had_0_65");
    can_res_rel_had_0_65->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->Draw("P,e,same");
    TF1* fit_RECO_0_65 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65->GetFunction("gaus");
    
    std::cout<<"gre_gcc7->SetPoint(17,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_65->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc7->SetPointError(17,0,"<<100.*fit_RECO_0_65->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_65_to_0_80=setUpperCanvas("can_h_res_rel_had_0_65_to_0_80");
    can_res_rel_had_0_65_to_0_80->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->Draw("P,e,same");
    TF1* fit_RECO_0_65_to_0_80 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_65_to_0_80->GetFunction("gaus");
    
    std::cout<<"gre_gcc7->SetPoint(17,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_65_to_0_80->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc7->SetPointError(17,0,"<<100.*fit_RECO_0_65_to_0_80->GetParError(2)<<");"<<std::endl;
    
    TCanvas* can_res_rel_had_0_80_to_0_94=setUpperCanvas("can_h_res_rel_had_0_80_to_0_94");
    can_res_rel_had_0_80_to_0_94->cd();
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->Draw("hist,e");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->Fit("gaus");
    h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->Draw("P,e,same");
    TF1* fit_RECO_0_80_to_0_94 = (TF1*)h_SWC_CLIC_n_K0L_E_reco_over_E_gen_0_80_to_0_94->GetFunction("gaus");
    
    std::cout<<"gre_gcc7->SetPoint(17,"<<true_energy_fillplot<<","<<100.*fit_RECO_0_80_to_0_94->GetParameter(2)<<");"<<std::endl;
    std::cout<<"gre_gcc7->SetPointError(17,0,"<<100.*fit_RECO_0_80_to_0_94->GetParError(2)<<");"<<std::endl;
    
    
    std::cout<<"gre_gcc7_mean->SetPoint(17,"<<true_energy_fillplot<<","<<fit_RECO_0_60->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPointError(17,0,"<<fit_RECO_0_60->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPoint(17,"<<true_energy_fillplot<<","<<fit_RECO_0_65->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPointError(17,0,"<<fit_RECO_0_65->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPoint(17,"<<true_energy_fillplot<<","<<fit_RECO_0_65_to_0_80->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPointError(17,0,"<<fit_RECO_0_65_to_0_80->GetParError(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPoint(17,"<<true_energy_fillplot<<","<<fit_RECO_0_80_to_0_94->GetParameter(1)<<");"<<std::endl;
    std::cout<<"gre_gcc7_mean->SetPointError(17,0,"<<fit_RECO_0_80_to_0_94->GetParError(1)<<");"<<std::endl;
    
    
    TCanvas* can_cosThetaTR_gcc62=setUpperCanvas("can_cosThetaTR_gcc62");
    can_cosThetaTR_gcc62->cd();
    h_SWC_CLIC_n_K0L_cosThetaTR->Draw();
    
    std::cout<<"file end of this costheta_min/max "<<true_CosTheta_min<<"/"<<true_CosTheta_max<<" cand found/tot "<< count_hadron_candidate_events <<"/"<< tree_SWC_CLIC_n_K0L->GetEntries() <<std::endl;
    
  */
    //try barrel up to 0.6, endcap 0.82 to 0.92 TR 0.65 to 0.80
    
    //now photon 0.8 to 0.8
    //=========Macro generated from canvas: resolution_gcc62_GraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L = new TCanvas("resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFillColor(0);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderMode(0);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderSize(2);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridx();
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridy();
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetRightMargin(0.0172);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetTopMargin(0.055);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBottomMargin(0.138);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    resolution_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre_gcc62 = new TGraphErrors(18);
    gre_gcc62->SetName("CLIC_o3_v14_res_graph_K0L_rel_E_reco_vs_E_true_0_60");
    /*  gre_gcc62->SetTitle("");
    gre_gcc62->SetFillColor(1);
    gre_gcc62->SetMarkerColor(kBlack);
    gre_gcc62->SetMarkerStyle(kOpenCircle);
    gre_gcc62->SetPoint(0,2,36.4348);
    gre_gcc62->SetPointError(0,0,0.766461);
    gre_gcc62->SetPoint(1,5,25.66);
    gre_gcc62->SetPointError(1,0,0.247524);
    gre_gcc62->SetPoint(2,10,17.2658);
    gre_gcc62->SetPointError(2,0,0.0907592);
    gre_gcc62->SetPoint(3,20,11.1261);
    gre_gcc62->SetPointError(3,0,0.056156);
    gre_gcc62->SetPoint(4,30,8.77582);
    gre_gcc62->SetPointError(4,0,0.0414576);
    gre_gcc62->SetPoint(5,40,7.64806);
    gre_gcc62->SetPointError(5,0,0.0378284);
    gre_gcc62->SetPoint(6,50,7.02321);
    gre_gcc62->SetPointError(6,0,0.0299756);
    gre_gcc62->SetPoint(7,60,6.69508);
    gre_gcc62->SetPointError(7,0,0.0313796);
    gre_gcc62->SetPoint(8,75,6.36158);
    gre_gcc62->SetPointError(8,0,0.0265109);
    gre_gcc62->SetPoint(9,90,6.17088);
    gre_gcc62->SetPointError(9,0,0.0262288);
    gre_gcc62->SetPoint(10,100,5.9945);
    gre_gcc62->SetPointError(10,0,0.0253894);
    gre_gcc62->SetPoint(11,150,5.75223);
    gre_gcc62->SetPointError(11,0,0.0249522);
    gre_gcc62->SetPoint(12,200,5.73667);
    gre_gcc62->SetPointError(12,0,0.0245236);
    gre_gcc62->SetPoint(13,250,5.72887);
    gre_gcc62->SetPointError(13,0,0.0252957);
    gre_gcc62->SetPoint(14,400,5.59588);
    gre_gcc62->SetPointError(14,0,0.024268);
    gre_gcc62->SetPoint(15,500,5.52831);
    gre_gcc62->SetPointError(15,0,0.0235202);
    gre_gcc62->SetPoint(16,1000,5.3636);
    gre_gcc62->SetPointError(16,0,0.0265645);
    gre_gcc62->SetPoint(17,1500,5.47889);
    gre_gcc62->SetPointError(17,0,0.0955246);

    
    //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    //gre_gcc62->Draw("ape");
    
    TF1* fit_sigma_K0L_0_60 = new TF1("fit_sigma_K0L_0_60", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_K0L_0_60->SetParameter(0,50.);
    fit_sigma_K0L_0_60->SetParameter(1,0.10);
    fit_sigma_K0L_0_60->SetParameter(2,3.0);
    fit_sigma_K0L_0_60->SetParNames("a1_0_60","a2_0_60","const_0_60");
    //gre_gcc62->Fit("fit_sigma_K0L_0_60","R");
    
    TF1* fit_sigma_K0L_0_60_only_H = new TF1("fit_sigma_K0L_0_60_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415.);
    fit_sigma_K0L_0_60_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_60_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_60_only_H->SetParNames("a1_0_60_only_H","const_0_60_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_60_only_H","R");
    
    //gre_gcc62->Draw("pesame");
    
    //0-30 
    gre_gcc62 = new TGraphErrors(18);
    gre_gcc62->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_30");
    gre_gcc62->SetTitle("");
    gre_gcc62->SetFillColor(1);
    gre_gcc62->SetMarkerColor(kBlack);
    gre_gcc62->SetMarkerStyle(kOpenCircle);

    gre_gcc62->SetPoint(0,2,36.8365);
    gre_gcc62->SetPointError(0,0,1.01606);
    gre_gcc62->SetPoint(1,5,25.9393);
    gre_gcc62->SetPointError(1,0,0.327916);
    gre_gcc62->SetPoint(2,10,17.5172);
    gre_gcc62->SetPointError(2,0,0.121292);
    gre_gcc62->SetPoint(3,20,11.1017);
    gre_gcc62->SetPointError(3,0,0.0722639);
    gre_gcc62->SetPoint(4,30,8.87741);
    gre_gcc62->SetPointError(4,0,0.0538629);
    gre_gcc62->SetPoint(5,40,7.73324);
    gre_gcc62->SetPointError(5,0,0.050739);
    gre_gcc62->SetPoint(6,50,7.11303);
    gre_gcc62->SetPointError(6,0,0.0393035);
    gre_gcc62->SetPoint(7,60,6.7594);
    gre_gcc62->SetPointError(7,0,0.0404086);
    gre_gcc62->SetPoint(8,75,6.48149);
    gre_gcc62->SetPointError(8,0,0.0347541);
    gre_gcc62->SetPoint(9,90,6.26101);
    gre_gcc62->SetPointError(9,0,0.0338716);
    gre_gcc62->SetPoint(10,100,6.08244);
    gre_gcc62->SetPointError(10,0,0.0340026);
    gre_gcc62->SetPoint(11,150,5.86645);
    gre_gcc62->SetPointError(11,0,0.0334827);
    gre_gcc62->SetPoint(12,200,5.89237);
    gre_gcc62->SetPointError(12,0,0.0319589);
    gre_gcc62->SetPoint(13,250,5.93018);
    gre_gcc62->SetPointError(13,0,0.0338517);
    gre_gcc62->SetPoint(14,400,5.84802);
    gre_gcc62->SetPointError(14,0,0.0326581);
    gre_gcc62->SetPoint(15,500,5.89874);
    gre_gcc62->SetPointError(15,0,0.032282);
    gre_gcc62->SetPoint(16,1000,5.79678);
    gre_gcc62->SetPointError(16,0,0.0358639);
    gre_gcc62->SetPoint(17,1500,5.62048);
    gre_gcc62->SetPointError(17,0,0.127887);


    
     TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14);
    //gre_gcc62->Draw("apesame");
    TF1* fit_sigma_K0L_0_30 = new TF1("fit_sigma_K0L_0_30", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,425);
    fit_sigma_K0L_0_30->SetParameter(0,50.);
    fit_sigma_K0L_0_30->SetParameter(1,0.10);
    fit_sigma_K0L_0_30->SetParameter(2,3.0);
    fit_sigma_K0L_0_30->SetParNames("a1_0_30","a2_0_30","const_0_30");
    //gre_gcc62->Fit("fit_sigma_K0L_0_30","R");
    
    TF1* fit_sigma_K0L_0_30_only_H = new TF1("fit_sigma_K0L_0_30_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_K0L_0_30_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_30_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_30_only_H->SetParNames("a1_0_30_only_H","const_0_30_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_30_only_H","R");
    
    //gre_gcc62->Draw("pesame");


    gre_gcc62 = new TGraphErrors(18);
    gre_gcc62->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_30_to_0_65");
    gre_gcc62->SetTitle("");
    gre_gcc62->SetFillColor(1);
    gre_gcc62->SetMarkerColor(kRed-7);
    gre_gcc62->SetMarkerStyle(kOpenSquare);

    gre_gcc62->SetPoint(0,2,36.5934);
    gre_gcc62->SetPointError(0,0,1.04088);
    gre_gcc62->SetPoint(1,5,25.4568);
    gre_gcc62->SetPointError(1,0,0.322502);
    gre_gcc62->SetPoint(2,10,16.7441);
    gre_gcc62->SetPointError(2,0,0.112789);
    gre_gcc62->SetPoint(3,20,10.8584);
    gre_gcc62->SetPointError(3,0,0.0693925);
    gre_gcc62->SetPoint(4,30,8.64518);
    gre_gcc62->SetPointError(4,0,0.053494);
    gre_gcc62->SetPoint(5,40,7.60212);
    gre_gcc62->SetPointError(5,0,0.049676);
    gre_gcc62->SetPoint(6,50,7.00106);
    gre_gcc62->SetPointError(6,0,0.0395396);
    gre_gcc62->SetPoint(7,60,6.69625);
    gre_gcc62->SetPointError(7,0,0.0414138);
    gre_gcc62->SetPoint(8,75,6.3295);
    gre_gcc62->SetPointError(8,0,0.0345101);
    gre_gcc62->SetPoint(9,90,6.21146);
    gre_gcc62->SetPointError(9,0,0.0349074);
    gre_gcc62->SetPoint(10,100,5.99391);
    gre_gcc62->SetPointError(10,0,0.0334758);
    gre_gcc62->SetPoint(11,150,5.74977);
    gre_gcc62->SetPointError(11,0,0.0328801);
    gre_gcc62->SetPoint(12,200,5.69434);
    gre_gcc62->SetPointError(12,0,0.0317499);
    gre_gcc62->SetPoint(13,250,5.68273);
    gre_gcc62->SetPointError(13,0,0.0338116);
    gre_gcc62->SetPoint(14,400,5.50603);
    gre_gcc62->SetPointError(14,0,0.0308761);
    gre_gcc62->SetPoint(15,500,5.42815);
    gre_gcc62->SetPointError(15,0,0.0297664);
    gre_gcc62->SetPoint(16,1000,5.2198);
    gre_gcc62->SetPointError(16,0,0.0349558);
    gre_gcc62->SetPoint(17,1500,5.46587);
    gre_gcc62->SetPointError(17,0,0.116198);

    
    TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14","",100,0,415);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetLineColor(kRed-7);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14);
    //gre_gcc62->Draw("pesame");
    TF1* fit_sigma_K0L_0_30_to_0_65 = new TF1("fit_sigma_K0L_0_30_to_0_65", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,425);
    fit_sigma_K0L_0_30_to_0_65->SetParameter(0,50.);
    fit_sigma_K0L_0_30_to_0_65->SetParameter(1,0.10);
    fit_sigma_K0L_0_30_to_0_65->SetParameter(2,3.0);
    fit_sigma_K0L_0_30_to_0_65->SetParNames("a1_0_30_to_0_65","a2_0_30_to_0_65","const_0_30_to_0_65");
    //gre_gcc62->Fit("fit_sigma_K0L_0_30_to_0_65","R");
    
    TF1* fit_sigma_K0L_0_30_to_0_65_only_H = new TF1("fit_sigma_K0L_0_30_to_0_65_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_K0L_0_30_to_0_65_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_30_to_0_65_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_30_to_0_65_only_H->SetParNames("a1_0_30_to_0_65_only_H","const_0_30_to_0_65_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_30_to_0_65_only_H","R");
    
    //gre_gcc62->Draw("pesame");
    */

    TGraphErrors* gre_gcc62_0_65 = new TGraphErrors(18);
    gre_gcc62_0_65->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65");
    gre_gcc62_0_65->SetTitle("");
    gre_gcc62_0_65->SetFillColor(1);
    gre_gcc62_0_65->SetMarkerColor(kBlack);
    gre_gcc62_0_65->SetMarkerStyle(kOpenCircle);
    gre_gcc62_0_65->SetPoint(0,2,36.4344);
    gre_gcc62_0_65->SetPointError(0,0,0.739078);
    gre_gcc62_0_65->SetPoint(1,5,25.6253);
    gre_gcc62_0_65->SetPointError(1,0,0.237696);
    gre_gcc62_0_65->SetPoint(2,10,17.2049);
    gre_gcc62_0_65->SetPointError(2,0,0.0868783);
    gre_gcc62_0_65->SetPoint(3,20,11.1274);
    gre_gcc62_0_65->SetPointError(3,0,0.054239);
    gre_gcc62_0_65->SetPoint(4,30,8.79933);
    gre_gcc62_0_65->SetPointError(4,0,0.0401531);
    gre_gcc62_0_65->SetPoint(5,40,7.66443);
    gre_gcc62_0_65->SetPointError(5,0,0.0367902);
    gre_gcc62_0_65->SetPoint(6,50,7.03673);
    gre_gcc62_0_65->SetPointError(6,0,0.0291822);
    gre_gcc62_0_65->SetPoint(7,60,6.70435);
    gre_gcc62_0_65->SetPointError(7,0,0.0301525);
    gre_gcc62_0_65->SetPoint(8,75,6.35791);
    gre_gcc62_0_65->SetPointError(8,0,0.0255388);
    gre_gcc62_0_65->SetPoint(9,90,6.19407);
    gre_gcc62_0_65->SetPointError(9,0,0.0253676);
    gre_gcc62_0_65->SetPoint(10,100,6.00422);
    gre_gcc62_0_65->SetPointError(10,0,0.0245256);
    gre_gcc62_0_65->SetPoint(11,150,5.76138);
    gre_gcc62_0_65->SetPointError(11,0,0.0240602);
    gre_gcc62_0_65->SetPoint(12,200,5.75447);
    gre_gcc62_0_65->SetPointError(12,0,0.0237341);
    gre_gcc62_0_65->SetPoint(13,250,5.72028);
    gre_gcc62_0_65->SetPointError(13,0,0.0240146);
    gre_gcc62_0_65->SetPoint(14,400,5.59529);
    gre_gcc62_0_65->SetPointError(14,0,0.0232333);
    gre_gcc62_0_65->SetPoint(15,500,5.53228);
    gre_gcc62_0_65->SetPointError(15,0,0.0226888);
    gre_gcc62_0_65->SetPoint(16,1000,5.3759);
    gre_gcc62_0_65->SetPointError(16,0,0.0258426);
    gre_gcc62_0_65->SetPoint(17,1500,5.48218);
    gre_gcc62_0_65->SetPointError(17,0,0.0910431);
    
    TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_0_65->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14);
    gre_gcc62_0_65->Draw("ape");
    TF1* fit_sigma_K0L_0_65 = new TF1("fit_sigma_K0L_0_65", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_K0L_0_65->SetParameter(0,50.);
    fit_sigma_K0L_0_65->SetParameter(1,0.10);
    fit_sigma_K0L_0_65->SetParameter(2,3.0);
    fit_sigma_K0L_0_65->SetParNames("a1_0_65","a2_0_65","const_0_65");
    //gre_gcc62->Fit("fit_sigma_K0L_0_65","R");
    
    TF1* fit_sigma_K0L_0_65_only_H = new TF1("fit_sigma_K0L_0_65_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_K0L_0_65_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_65_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_65_only_H->SetParNames("a1_0_65_only_H","const_0_65_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_65_only_H","R");
    
    //gre_gcc62->Draw("pesame");
    
    gre_gcc62 = new TGraphErrors(18);
    gre_gcc62->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80");
    gre_gcc62->SetTitle("");
    gre_gcc62->SetFillColor(1);
    gre_gcc62->SetMarkerColor(kRed-7);
    gre_gcc62->SetMarkerStyle(kOpenSquare);
    gre_gcc62->SetPoint(0,2,36.4833);
    gre_gcc62->SetPointError(0,0,1.7139);
    gre_gcc62->SetPoint(1,5,25.8934);
    gre_gcc62->SetPointError(1,0,0.559949);
    gre_gcc62->SetPoint(2,10,16.4662);
    gre_gcc62->SetPointError(2,0,0.180735);
    gre_gcc62->SetPoint(3,20,11.089);
    gre_gcc62->SetPointError(3,0,0.131805);
    gre_gcc62->SetPoint(4,30,8.6318);
    gre_gcc62->SetPointError(4,0,0.100182);
    gre_gcc62->SetPoint(5,40,7.87732);
    gre_gcc62->SetPointError(5,0,0.105985);
    gre_gcc62->SetPoint(6,50,7.26023);
    gre_gcc62->SetPointError(6,0,0.075536);
    gre_gcc62->SetPoint(7,60,6.73712);
    gre_gcc62->SetPointError(7,0,0.0799938);
    gre_gcc62->SetPoint(8,75,6.56918);
    gre_gcc62->SetPointError(8,0,0.0706449);
    gre_gcc62->SetPoint(9,90,6.11101);
    gre_gcc62->SetPointError(9,0,0.0666521);
    gre_gcc62->SetPoint(10,100,6.19851);
    gre_gcc62->SetPointError(10,0,0.0745047);
    gre_gcc62->SetPoint(11,150,5.65581);
    gre_gcc62->SetPointError(11,0,0.0623192);
    gre_gcc62->SetPoint(12,200,5.56125);
    gre_gcc62->SetPointError(12,0,0.069142);
    gre_gcc62->SetPoint(13,250,5.6564);
    gre_gcc62->SetPointError(13,0,0.0697636);
    gre_gcc62->SetPoint(14,400,5.56311);
    gre_gcc62->SetPointError(14,0,0.0739029);
    gre_gcc62->SetPoint(15,500,5.4082);
    gre_gcc62->SetPointError(15,0,0.0683279);
    gre_gcc62->SetPoint(16,1000,5.15672);
    gre_gcc62->SetPointError(16,0,0.0787054);
    gre_gcc62->SetPoint(17,1500,5.13651);
    gre_gcc62->SetPointError(17,0,0.355858);
    
    TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,1649.9);
    //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,415);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetLineColor(kRed-7);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14);
    //gre_gcc62->Draw("pe");
    TF1* fit_sigma_K0L_0_65_to_0_80 = new TF1("fit_sigma_K0L_0_65_to_0_80", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_K0L_0_65_to_0_80->SetParameter(0,50.);
    fit_sigma_K0L_0_65_to_0_80->SetParameter(1,0.10);
    fit_sigma_K0L_0_65_to_0_80->SetParameter(2,3.0);
    fit_sigma_K0L_0_65_to_0_80->SetParNames("a1_0_65_to_0_80","a2_0_65_to_0_80","const_0_65_to_0_80");
    //gre_gcc62->Fit("fit_sigma_K0L_0_65_to_0_80","R");
    
    TF1* fit_sigma_K0L_0_65_to_0_80_only_H = new TF1("fit_sigma_K0L_0_65_to_0_80_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_K0L_0_65_to_0_80_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_65_to_0_80_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_65_to_0_80_only_H->SetParNames("a1_0_65_to_0_80_only_H","const_0_65_to_0_80_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_65_to_0_80_only_H","R");
    
    gre_gcc62->Draw("pesame");
    
    gre_gcc62 = new TGraphErrors(18);
    gre_gcc62->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94");
    gre_gcc62->SetTitle("");
    gre_gcc62->SetFillColor(1);
    gre_gcc62->SetMarkerColor(kBlue);
    gre_gcc62->SetMarkerStyle(kOpenTriangleUp);
    gre_gcc62->SetPoint(0,2,37.6859);
    gre_gcc62->SetPointError(0,0,1.79119);
    gre_gcc62->SetPoint(1,5,24.5551);
    gre_gcc62->SetPointError(1,0,0.453611);
    gre_gcc62->SetPoint(2,10,15.7351);
    gre_gcc62->SetPointError(2,0,0.167803);
    gre_gcc62->SetPoint(3,20,10.4049);
    gre_gcc62->SetPointError(3,0,0.114396);
    gre_gcc62->SetPoint(4,30,8.75853);
    gre_gcc62->SetPointError(4,0,0.0836679);
    gre_gcc62->SetPoint(5,40,7.51838);
    gre_gcc62->SetPointError(5,0,0.0756588);
    gre_gcc62->SetPoint(6,50,7.3408);
    gre_gcc62->SetPointError(6,0,0.0643899);
    gre_gcc62->SetPoint(7,60,6.8912);
    gre_gcc62->SetPointError(7,0,0.0660442);
    gre_gcc62->SetPoint(8,75,6.5579);
    gre_gcc62->SetPointError(8,0,0.0562063);
    gre_gcc62->SetPoint(9,90,6.37767);
    gre_gcc62->SetPointError(9,0,0.0569723);
    gre_gcc62->SetPoint(10,100,6.12034);
    gre_gcc62->SetPointError(10,0,0.0598237);
    gre_gcc62->SetPoint(11,150,5.6872);
    gre_gcc62->SetPointError(11,0,0.0537094);
    gre_gcc62->SetPoint(12,200,5.6244);
    gre_gcc62->SetPointError(12,0,0.053084);
    gre_gcc62->SetPoint(13,250,5.57162);
    gre_gcc62->SetPointError(13,0,0.0609708);
    gre_gcc62->SetPoint(14,400,5.56826);
    gre_gcc62->SetPointError(14,0,0.0558178);
    gre_gcc62->SetPoint(15,500,5.5777);
    gre_gcc62->SetPointError(15,0,0.0534595);
    gre_gcc62->SetPoint(16,1000,5.42587);
    gre_gcc62->SetPointError(16,0,0.0634342);
    gre_gcc62->SetPoint(17,1500,5.13964);
    gre_gcc62->SetPointError(17,0,0.196675);
    
    TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,1649.9);
   //TH1F *h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62->SetHistogram(h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14);
    //gre_gcc62->Draw("pe");
    TF1* fit_sigma_K0L_0_80_to_0_92 = new TF1("fit_sigma_K0L_0_80_to_0_92", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_K0L_0_80_to_0_92->SetParameter(0,50.0);
    fit_sigma_K0L_0_80_to_0_92->SetParameter(1,0.10);
    fit_sigma_K0L_0_80_to_0_92->SetParameter(2,3.0);
    fit_sigma_K0L_0_80_to_0_92->SetParNames("a1_0_80_to_0_92","a2_0_80_to_0_92","const_0_80_to_0_92");
    //gre_gcc62->Fit("fit_sigma_K0L_0_80_to_0_92","R");
    
    TF1* fit_sigma_K0L_0_80_to_0_92_only_H = new TF1("fit_sigma_K0L_0_80_to_0_92_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_K0L_0_80_to_0_92_only_H->SetParameter(0,50.0);
    fit_sigma_K0L_0_80_to_0_92_only_H->SetParameter(1,3.0);
    fit_sigma_K0L_0_80_to_0_92_only_H->SetParNames("a1_0_80_to_0_92_only_H","const_0_80_to_0_92_only_H");
    //gre_gcc62->Fit("fit_sigma_K0L_0_80_to_0_92_only_H","R");
    gre_gcc62->Draw("pesame");
    
    
    
    
    TLegend *leg_K0L_E_reco_vs_E_true = new TLegend(0.60,0.60,0.94,0.91,"K0L","brNDC");
    leg_K0L_E_reco_vs_E_true->SetBorderSize(1);
    leg_K0L_E_reco_vs_E_true->SetTextSize(0.027);
    leg_K0L_E_reco_vs_E_true->SetLineColor(1);
    leg_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_K0L_E_reco_vs_E_true->SetLineWidth(1);
    leg_K0L_E_reco_vs_E_true->SetFillColor(0);
    leg_K0L_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_K0L_E_reco_vs_E_true=leg_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_K0L_E_reco_vs_E_true=leg_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","0.65<|cos#theta_{K0L}^{true}|<0.80","pe");
    leg_entry_K0L_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_K0L_E_reco_vs_E_true=leg_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_80_to_094_CLIC_o3_v14","0.80<|cos#theta_{K0L}^{true}|<0.94","pe");
    leg_entry_K0L_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_K0L_E_reco_vs_E_true->Draw();
    
    /*

   TLegend *leg_K0L_E_reco_vs_E_true_V2 = new TLegend(0.60,0.60,0.94,0.91,"K0L","brNDC");
    leg_K0L_E_reco_vs_E_true_V2->SetBorderSize(1);
    leg_K0L_E_reco_vs_E_true_V2->SetTextSize(0.027);
    leg_K0L_E_reco_vs_E_true_V2->SetLineColor(1);
    leg_K0L_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_K0L_E_reco_vs_E_true_V2->SetLineWidth(1);
    leg_K0L_E_reco_vs_E_true_V2->SetFillColor(0);
    leg_K0L_E_reco_vs_E_true_V2->SetFillStyle(1001);
    TLegendEntry *leg_entry_K0L_E_reco_vs_E_true_V2=leg_K0L_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineColor(kBlack);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerColor(kBlack);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenCircle);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetTextFont(42);
    leg_entry_K0L_E_reco_vs_E_true_V2=leg_K0L_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_30_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.30","pe");
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineColor(kRed-7);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerColor(kRed-7);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenSquare);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetTextFont(42);
    leg_entry_K0L_E_reco_vs_E_true_V2=leg_K0L_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_30_to_065_CLIC_o3_v14","0.30<|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineColor(kBlue);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerColor(kBlue);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_E_reco_vs_E_true_V2->SetTextFont(42);
    //leg_K0L_E_reco_vs_E_true_V2->Draw();
    //now for mean
    
    //=========Macro generated from canvas: resolution_gcc62_GraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L = new TCanvas("resolution_gcc62_mean_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFillColor(0);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderMode(0);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderSize(2);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridx();
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridy();
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetRightMargin(0.0172);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetTopMargin(0.055);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBottomMargin(0.138);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    resolution_gcc62_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre_gcc62_mean_0_60 = new TGraphErrors(18);
    gre_gcc62_mean_0_60->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_60");
    gre_gcc62_mean_0_60->SetTitle("");
    gre_gcc62_mean_0_60->SetFillColor(1);
    gre_gcc62_mean_0_60->SetMarkerColor(kBlack);
    gre_gcc62_mean_0_60->SetMarkerStyle(kOpenCircle);
    gre_gcc62_mean_0_60->SetPoint(0,2,0.761003);
    gre_gcc62_mean_0_60->SetPointError(0,0,0.00932005);
    gre_gcc62_mean_0_60->SetPoint(1,5,0.928753);
    gre_gcc62_mean_0_60->SetPointError(1,0,0.00213603);
    gre_gcc62_mean_0_60->SetPoint(2,10,1.03667);
    gre_gcc62_mean_0_60->SetPointError(2,0,0.0010708);
    gre_gcc62_mean_0_60->SetPoint(3,20,1.04368);
    gre_gcc62_mean_0_60->SetPointError(3,0,0.000755486);
    gre_gcc62_mean_0_60->SetPoint(4,30,1.03477);
    gre_gcc62_mean_0_60->SetPointError(4,0,0.000561433);
    gre_gcc62_mean_0_60->SetPoint(5,40,1.02531);
    gre_gcc62_mean_0_60->SetPointError(5,0,0.0005218);
    gre_gcc62_mean_0_60->SetPoint(6,50,1.01999);
    gre_gcc62_mean_0_60->SetPointError(6,0,0.000401458);
    gre_gcc62_mean_0_60->SetPoint(7,60,1.01817);
    gre_gcc62_mean_0_60->SetPointError(7,0,0.000423026);
    gre_gcc62_mean_0_60->SetPoint(8,75,1.01802);
    gre_gcc62_mean_0_60->SetPointError(8,0,0.000360629);
    gre_gcc62_mean_0_60->SetPoint(9,90,1.01826);
    gre_gcc62_mean_0_60->SetPointError(9,0,0.000349172);
    gre_gcc62_mean_0_60->SetPoint(10,100,1.01807);
    gre_gcc62_mean_0_60->SetPointError(10,0,0.00034728);
    gre_gcc62_mean_0_60->SetPoint(11,150,1.01372);
    gre_gcc62_mean_0_60->SetPointError(11,0,0.000327611);
    gre_gcc62_mean_0_60->SetPoint(12,200,1.00885);
    gre_gcc62_mean_0_60->SetPointError(12,0,0.000317127);
    gre_gcc62_mean_0_60->SetPoint(13,250,1.00645);
    gre_gcc62_mean_0_60->SetPointError(13,0,0.000335543);
    gre_gcc62_mean_0_60->SetPoint(14,400,1.01153);
    gre_gcc62_mean_0_60->SetPointError(14,0,0.000324681);
    gre_gcc62_mean_0_60->SetPoint(15,500,1.01621);
    gre_gcc62_mean_0_60->SetPointError(15,0,0.00030627);
    gre_gcc62_mean_0_60->SetPoint(16,1000,0.995991);
    gre_gcc62_mean_0_60->SetPointError(16,0,0.000326217);
    gre_gcc62_mean_0_60->SetPoint(17,1500,0.957566);
    gre_gcc62_mean_0_60->SetPointError(17,0,0.0010725);
    
    //TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean_0_60->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    gre_gcc62_mean_0_60->Draw("ape");
    

   //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors* gre_gcc62_mean = new TGraphErrors(18);
    gre_gcc62_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_30");
    gre_gcc62_mean->SetTitle("");
    gre_gcc62_mean->SetFillColor(1);
    gre_gcc62_mean->SetMarkerColor(kBlack);
    gre_gcc62_mean->SetMarkerStyle(kOpenCircle);
    
    gre_gcc62_mean->SetPoint(0,2,0.776688);
    gre_gcc62_mean->SetPointError(0,0,0.0116593);
    gre_gcc62_mean->SetPoint(1,5,0.942629);
    gre_gcc62_mean->SetPointError(1,0,0.00276543);
    gre_gcc62_mean->SetPoint(2,10,1.04583);
    gre_gcc62_mean->SetPointError(2,0,0.00142457);
    gre_gcc62_mean->SetPoint(3,20,1.04824);
    gre_gcc62_mean->SetPointError(3,0,0.000985258);
    gre_gcc62_mean->SetPoint(4,30,1.03739);
    gre_gcc62_mean->SetPointError(4,0,0.000732433);
    gre_gcc62_mean->SetPoint(5,40,1.02689);
    gre_gcc62_mean->SetPointError(5,0,0.000685053);
    gre_gcc62_mean->SetPoint(6,50,1.01965);
    gre_gcc62_mean->SetPointError(6,0,0.000525965);
    gre_gcc62_mean->SetPoint(7,60,1.01771);
    gre_gcc62_mean->SetPointError(7,0,0.000555794);
    gre_gcc62_mean->SetPoint(8,75,1.01623);
    gre_gcc62_mean->SetPointError(8,0,0.000477194);
    gre_gcc62_mean->SetPoint(9,90,1.01638);
    gre_gcc62_mean->SetPointError(9,0,0.000460068);
    gre_gcc62_mean->SetPoint(10,100,1.01558);
    gre_gcc62_mean->SetPointError(10,0,0.000458638);
    gre_gcc62_mean->SetPoint(11,150,1.00984);
    gre_gcc62_mean->SetPointError(11,0,0.000434103);
    gre_gcc62_mean->SetPoint(12,200,1.00291);
    gre_gcc62_mean->SetPointError(12,0,0.000425818);
    gre_gcc62_mean->SetPoint(13,250,1.00067);
    gre_gcc62_mean->SetPointError(13,0,0.00045644);
    gre_gcc62_mean->SetPoint(14,400,1.0035);
    gre_gcc62_mean->SetPointError(14,0,0.000449702);
    gre_gcc62_mean->SetPoint(15,500,1.00784);
    gre_gcc62_mean->SetPointError(15,0,0.000433105);
    gre_gcc62_mean->SetPoint(16,1000,0.985106);
    gre_gcc62_mean->SetPointError(16,0,0.000469488);
    gre_gcc62_mean->SetPoint(17,1500,0.94578);
    gre_gcc62_mean->SetPointError(17,0,0.00145928);

    //TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14);
    gre_gcc62_mean->Draw("ape");



    gre_gcc62_mean = new TGraphErrors(18);
    gre_gcc62_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_30_to_0_65");
    gre_gcc62_mean->SetTitle("");
    gre_gcc62_mean->SetFillColor(1);
    gre_gcc62_mean->SetMarkerColor(kBlue);
    gre_gcc62_mean->SetMarkerStyle(kOpenSquare);

    gre_gcc62_mean->SetPoint(0,2,0.73919);
    gre_gcc62_mean->SetPointError(0,0,0.0134323);
    gre_gcc62_mean->SetPoint(1,5,0.917637);
    gre_gcc62_mean->SetPointError(1,0,0.00285225);
    gre_gcc62_mean->SetPoint(2,10,1.02535);
    gre_gcc62_mean->SetPointError(2,0,0.00133765);
    gre_gcc62_mean->SetPoint(3,20,1.03403);
    gre_gcc62_mean->SetPointError(3,0,0.000984791);
    gre_gcc62_mean->SetPoint(4,30,1.02702);
    gre_gcc62_mean->SetPointError(4,0,0.000742764);
    gre_gcc62_mean->SetPoint(5,40,1.01855);
    gre_gcc62_mean->SetPointError(5,0,0.00069332);
    gre_gcc62_mean->SetPoint(6,50,1.01555);
    gre_gcc62_mean->SetPointError(6,0,0.00053296);
    gre_gcc62_mean->SetPoint(7,60,1.0142);
    gre_gcc62_mean->SetPointError(7,0,0.000566336);
    gre_gcc62_mean->SetPoint(8,75,1.01474);
    gre_gcc62_mean->SetPointError(8,0,0.000479769);
    gre_gcc62_mean->SetPoint(9,90,1.01542);
    gre_gcc62_mean->SetPointError(9,0,0.000468723);
    gre_gcc62_mean->SetPoint(10,100,1.01616);
    gre_gcc62_mean->SetPointError(10,0,0.000460277);
    gre_gcc62_mean->SetPoint(11,150,1.01321);
    gre_gcc62_mean->SetPointError(11,0,0.000436144);
    gre_gcc62_mean->SetPoint(12,200,1.00926);
    gre_gcc62_mean->SetPointError(12,0,0.00041671);
    gre_gcc62_mean->SetPoint(13,250,1.00676);
    gre_gcc62_mean->SetPointError(13,0,0.000438275);
    gre_gcc62_mean->SetPoint(14,400,1.01329);
    gre_gcc62_mean->SetPointError(14,0,0.000419007);
    gre_gcc62_mean->SetPoint(15,500,1.01739);
    gre_gcc62_mean->SetPointError(15,0,0.000398619);
    gre_gcc62_mean->SetPoint(16,1000,0.99967);
    gre_gcc62_mean->SetPointError(16,0,0.000415391);
    gre_gcc62_mean->SetPoint(17,1500,0.96289);
    gre_gcc62_mean->SetPointError(17,0,0.00140091);

    
   // TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14);
    gre_gcc62_mean->Draw("pe");




    gre_gcc62_mean = new TGraphErrors(18);
    gre_gcc62_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_65");
    gre_gcc62_mean->SetTitle("");
    gre_gcc62_mean->SetFillColor(1);
    gre_gcc62_mean->SetMarkerColor(kBlack);
    gre_gcc62_mean->SetMarkerStyle(kOpenCircle);
    gre_gcc62_mean->SetPoint(0,2,0.757477);
    gre_gcc62_mean->SetPointError(0,0,0.00907269);
    gre_gcc62_mean->SetPoint(1,5,0.926023);
    gre_gcc62_mean->SetPointError(1,0,0.00206455);
    gre_gcc62_mean->SetPoint(2,10,1.03439);
    gre_gcc62_mean->SetPointError(2,0,0.0010222);
    gre_gcc62_mean->SetPoint(3,20,1.04188);
    gre_gcc62_mean->SetPointError(3,0,0.000729367);
    gre_gcc62_mean->SetPoint(4,30,1.03319);
    gre_gcc62_mean->SetPointError(4,0,0.000543671);
    gre_gcc62_mean->SetPoint(5,40,1.02397);
    gre_gcc62_mean->SetPointError(5,0,0.000505532);
    gre_gcc62_mean->SetPoint(6,50,1.01896);
    gre_gcc62_mean->SetPointError(6,0,0.000388274);
    gre_gcc62_mean->SetPoint(7,60,1.01743);
    gre_gcc62_mean->SetPointError(7,0,0.000409934);
    gre_gcc62_mean->SetPoint(8,75,1.01729);
    gre_gcc62_mean->SetPointError(8,0,0.000348357);
    gre_gcc62_mean->SetPoint(9,90,1.01759);
    gre_gcc62_mean->SetPointError(9,0,0.000338707);
    gre_gcc62_mean->SetPoint(10,100,1.01753);
    gre_gcc62_mean->SetPointError(10,0,0.000335537);
    gre_gcc62_mean->SetPoint(11,150,1.01367);
    gre_gcc62_mean->SetPointError(11,0,0.000317173);
    gre_gcc62_mean->SetPoint(12,200,1.00851);
    gre_gcc62_mean->SetPointError(12,0,0.000307128);
    gre_gcc62_mean->SetPoint(13,250,1.00632);
    gre_gcc62_mean->SetPointError(13,0,0.000324169);
    gre_gcc62_mean->SetPoint(14,400,1.01136);
    gre_gcc62_mean->SetPointError(14,0,0.000312844);
    gre_gcc62_mean->SetPoint(15,500,1.01628);
    gre_gcc62_mean->SetPointError(15,0,0.000295874);
    gre_gcc62_mean->SetPoint(16,1000,0.996249);
    gre_gcc62_mean->SetPointError(16,0,0.000314423);
    gre_gcc62_mean->SetPoint(17,1500,0.958239);
    gre_gcc62_mean->SetPointError(17,0,0.00103222);

    
    //TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_0_65_CLIC_o3_v14);
    gre_gcc62_mean->Draw("ape");
    
    gre_gcc62_mean = new TGraphErrors(18);
    gre_gcc62_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_65_to_0_80");
    gre_gcc62_mean->SetTitle("");
    gre_gcc62_mean->SetFillColor(1);
    gre_gcc62_mean->SetMarkerColor(kBlue);
    gre_gcc62_mean->SetMarkerStyle(kOpenSquare);
    gre_gcc62_mean->SetPoint(0,2,0.68446);
    gre_gcc62_mean->SetPointError(0,0,0.0264486);
    gre_gcc62_mean->SetPoint(1,5,0.867045);
    gre_gcc62_mean->SetPointError(1,0,0.00579954);
    gre_gcc62_mean->SetPoint(2,10,0.988972);
    gre_gcc62_mean->SetPointError(2,0,0.00221837);
    gre_gcc62_mean->SetPoint(3,20,1.01336);
    gre_gcc62_mean->SetPointError(3,0,0.00180493);
    gre_gcc62_mean->SetPoint(4,30,1.01244);
    gre_gcc62_mean->SetPointError(4,0,0.00134642);
    gre_gcc62_mean->SetPoint(5,40,1.00996);
    gre_gcc62_mean->SetPointError(5,0,0.00137738);
    gre_gcc62_mean->SetPoint(6,50,1.00906);
    gre_gcc62_mean->SetPointError(6,0,0.00105222);
    gre_gcc62_mean->SetPoint(7,60,1.01039);
    gre_gcc62_mean->SetPointError(7,0,0.00115058);
    gre_gcc62_mean->SetPoint(8,75,1.01389);
    gre_gcc62_mean->SetPointError(8,0,0.000981681);
    gre_gcc62_mean->SetPoint(9,90,1.01535);
    gre_gcc62_mean->SetPointError(9,0,0.0009243);
    gre_gcc62_mean->SetPoint(10,100,1.01289);
    gre_gcc62_mean->SetPointError(10,0,0.000946568);
    gre_gcc62_mean->SetPoint(11,150,1.01401);
    gre_gcc62_mean->SetPointError(11,0,0.000863373);
    gre_gcc62_mean->SetPoint(12,200,1.00975);
    gre_gcc62_mean->SetPointError(12,0,0.000848057);
    gre_gcc62_mean->SetPoint(13,250,1.00772);
    gre_gcc62_mean->SetPointError(13,0,0.000914872);
    gre_gcc62_mean->SetPoint(14,400,1.01161);
    gre_gcc62_mean->SetPointError(14,0,0.000882629);
    gre_gcc62_mean->SetPoint(15,500,1.01508);
    gre_gcc62_mean->SetPointError(15,0,0.00082358);
    gre_gcc62_mean->SetPoint(16,1000,0.995081);
    gre_gcc62_mean->SetPointError(16,0,0.000970177);
    gre_gcc62_mean->SetPoint(17,1500,0.953178);
    gre_gcc62_mean->SetPointError(17,0,0.00330987);
    
   // TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14);
    gre_gcc62_mean->Draw("pe");
    
    gre_gcc62_mean = new TGraphErrors(18);
    gre_gcc62_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_mean_rel_E_reco_vs_E_true_0_80_to_0_94");
    gre_gcc62_mean->SetTitle("");
    gre_gcc62_mean->SetFillColor(1);
    gre_gcc62_mean->SetMarkerColor(kRed-7);
    gre_gcc62_mean->SetMarkerStyle(kOpenTriangleUp);
    gre_gcc62_mean->SetPoint(0,2,0.768566);
    gre_gcc62_mean->SetPointError(0,0,0.021228);
    gre_gcc62_mean->SetPoint(1,5,0.958651);
    gre_gcc62_mean->SetPointError(1,0,0.00411343);
    gre_gcc62_mean->SetPoint(2,10,1.05522);
    gre_gcc62_mean->SetPointError(2,0,0.0020215);
    gre_gcc62_mean->SetPoint(3,20,1.04806);
    gre_gcc62_mean->SetPointError(3,0,0.0015267);
    gre_gcc62_mean->SetPoint(4,30,1.03442);
    gre_gcc62_mean->SetPointError(4,0,0.00119396);
    gre_gcc62_mean->SetPoint(5,40,1.02199);
    gre_gcc62_mean->SetPointError(5,0,0.00108963);
    gre_gcc62_mean->SetPoint(6,50,1.01385);
    gre_gcc62_mean->SetPointError(6,0,0.000891492);
    gre_gcc62_mean->SetPoint(7,60,1.01188);
    gre_gcc62_mean->SetPointError(7,0,0.000946572);
    gre_gcc62_mean->SetPoint(8,75,1.00919);
    gre_gcc62_mean->SetPointError(8,0,0.000802718);
    gre_gcc62_mean->SetPoint(9,90,1.00865);
    gre_gcc62_mean->SetPointError(9,0,0.000799345);
    gre_gcc62_mean->SetPoint(10,100,1.00627);
    gre_gcc62_mean->SetPointError(10,0,0.000769762);
    gre_gcc62_mean->SetPoint(11,150,0.998351);
    gre_gcc62_mean->SetPointError(11,0,0.00071463);
    gre_gcc62_mean->SetPoint(12,200,0.992859);
    gre_gcc62_mean->SetPointError(12,0,0.000690747);
    gre_gcc62_mean->SetPoint(13,250,0.987722);
    gre_gcc62_mean->SetPointError(13,0,0.000738937);
    gre_gcc62_mean->SetPoint(14,400,0.988587);
    gre_gcc62_mean->SetPointError(14,0,0.000709825);
    gre_gcc62_mean->SetPoint(15,500,0.989561);
    gre_gcc62_mean->SetPointError(15,0,0.000669009);
    gre_gcc62_mean->SetPoint(16,1000,0.965582);
    gre_gcc62_mean->SetPointError(16,0,0.000732338);
    gre_gcc62_mean->SetPoint(17,1500,0.930193);
    gre_gcc62_mean->SetPointError(17,0,0.00245952);
    
    //TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc62_mean->SetHistogram(h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14);
    gre_gcc62_mean->Draw("pe");
    
    
    TLegend *leg_K0L_gcc62_mean_E_reco_vs_E_true = new TLegend(0.60,0.20,0.94,0.60,"K^{0}_{L}","brNDC");
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetBorderSize(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetTextSize(0.027);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetLineColor(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetLineWidth(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetFillColor(0);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_K0L_gcc62_mean_E_reco_vs_E_true=leg_K0L_gcc62_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true=leg_K0L_gcc62_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","0.65<|cos#theta_{K0L}^{true}|<0.80","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true=leg_K0L_gcc62_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_80_to_094_CLIC_o3_v14","0.80<|cos#theta_{K0L}^{true}|<0.94","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_K0L_gcc62_mean_E_reco_vs_E_true->Draw();
                

    TLegend *leg_K0L_gcc62_mean_E_reco_vs_E_true_V2 = new TLegend(0.60,0.20,0.94,0.60,"K^{0}_{L}","brNDC");
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetBorderSize(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetTextSize(0.027);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineColor(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineWidth(1);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetFillColor(0);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetFillStyle(1001);
    TLegendEntry *leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2=leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineColor(kBlack);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerColor(kBlack);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenCircle);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetTextFont(42);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2=leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.30","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineColor(kRed-7);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerColor(kRed-7);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenSquare);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetTextFont(42);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2=leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->AddEntry("h_K0L_gcc62_mean_E_reco_vs_E_true_0_30_to_065_CLIC_o3_v14","0.30<|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineColor(kBlue);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineStyle(1);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetLineWidth(2);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerColor(kBlue);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetMarkerSize(1.5);
    leg_entry_K0L_gcc62_mean_E_reco_vs_E_true_V2->SetTextFont(42);
    leg_K0L_gcc62_mean_E_reco_vs_E_true_V2->Draw();                 
    //now the gcc7 versions

    
   //try barrel up to 0.6, endcap 0.82 to 0.92 TR 0.65 to 0.80
    
    //now photon 0.8 to 0.8
    //=========Macro generated from canvas: resolution_gcc7_GraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L = new TCanvas("resolution_gcc7_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc7_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFillColor(0);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderMode(0);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderSize(2);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridx();
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridy();
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetRightMargin(0.0172);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetTopMargin(0.055);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBottomMargin(0.138);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    resolution_gcc7_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre_gcc7 = new TGraphErrors(18);
    gre_gcc7->SetName("CLIC_o3_v14_res_graph_K0L_rel_E_reco_vs_E_true_0_60");
    gre_gcc7->SetTitle("");
    gre_gcc7->SetFillColor(1);
    gre_gcc7->SetMarkerColor(kBlack);
    gre_gcc7->SetMarkerStyle(kOpenCircle);
    gre_gcc7->SetPoint(0,2,30.6591);
    gre_gcc7->SetPointError(0,0,0.180819);
    gre_gcc7->SetPoint(1,5,23.7907);
    gre_gcc7->SetPointError(1,0,0.106618);
    gre_gcc7->SetPoint(2,10,17.5111);
    gre_gcc7->SetPointError(2,0,0.0848294);
    gre_gcc7->SetPoint(3,20,10.6126);
    gre_gcc7->SetPointError(3,0,0.0449012);
    gre_gcc7->SetPoint(4,30,8.3018);
    gre_gcc7->SetPointError(4,0,0.0365336);
    gre_gcc7->SetPoint(5,40,7.14691);
    gre_gcc7->SetPointError(5,0,0.0316106);
    gre_gcc7->SetPoint(6,50,6.4392);
    gre_gcc7->SetPointError(6,0,0.0287832);
    gre_gcc7->SetPoint(7,60,5.97085);
    gre_gcc7->SetPointError(7,0,0.0282135);
    gre_gcc7->SetPoint(8,75,5.47232);
    gre_gcc7->SetPointError(8,0,0.0252764);
    gre_gcc7->SetPoint(9,90,5.1965);
    gre_gcc7->SetPointError(9,0,0.0241039);
    gre_gcc7->SetPoint(10,100,4.99579);
    gre_gcc7->SetPointError(10,0,0.0219897);
    gre_gcc7->SetPoint(11,150,4.4526);
    gre_gcc7->SetPointError(11,0,0.0192098);
    gre_gcc7->SetPoint(12,200,4.24156);
    gre_gcc7->SetPointError(12,0,0.0186166);
    gre_gcc7->SetPoint(13,250,3.98141);
    gre_gcc7->SetPointError(13,0,0.0174288);
    gre_gcc7->SetPoint(14,400,3.52371);
    gre_gcc7->SetPointError(14,0,0.0155333);
    
    gre_gcc7->SetPoint(16,1000,3.51693);
    gre_gcc7->SetPointError(16,0,0.0150795);
    gre_gcc7->SetPoint(17,1500,3.40879);
    gre_gcc7->SetPointError(17,0,0.039613);


    
    //TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7->SetHistogram(h_K0L_gcc7_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    gre_gcc7->Draw("ape");
    
    TF1* fit_sigma_gcc7_K0L_0_60 = new TF1("fit_sigma_gcc7_K0L_0_60", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_gcc7_K0L_0_60->SetParameter(0,50.);
    fit_sigma_gcc7_K0L_0_60->SetParameter(1,0.10);
    fit_sigma_gcc7_K0L_0_60->SetParameter(2,3.0);
    fit_sigma_gcc7_K0L_0_60->SetParNames("a1_0_60","a2_0_60","const_0_60");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_60","R");
    
    TF1* fit_sigma_gcc7_K0L_0_60_only_H = new TF1("fit_sigma_gcc7_K0L_0_60_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415.);
    fit_sigma_gcc7_K0L_0_60_only_H->SetParameter(0,50.0);
    fit_sigma_gcc7_K0L_0_60_only_H->SetParameter(1,3.0);
    fit_sigma_gcc7_K0L_0_60_only_H->SetParNames("a1_0_60_only_H","const_0_60_only_H");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_60_only_H","R");
    
    gre_gcc7->Draw("pesame");
    
    TGraphErrors* gre_gcc7_0_65 = new TGraphErrors(18);
    gre_gcc7_0_65->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_rel_E_reco_vs_E_true_0_65");
    gre_gcc7_0_65->SetTitle("");
    gre_gcc7_0_65->SetFillColor(1);
    gre_gcc7_0_65->SetMarkerColor(kBlack);
    gre_gcc7_0_65->SetMarkerStyle(kOpenCircle);
    gre_gcc7_0_65->SetPoint(0,2,30.5452);
    gre_gcc7_0_65->SetPointError(0,0,0.172326);
    gre_gcc7_0_65->SetPoint(1,5,23.7958);
    gre_gcc7_0_65->SetPointError(1,0,0.102803);
    gre_gcc7_0_65->SetPoint(2,10,17.5079);
    gre_gcc7_0_65->SetPointError(2,0,0.0810838);
    gre_gcc7_0_65->SetPoint(3,20,10.6373);
    gre_gcc7_0_65->SetPointError(3,0,0.0435838);
    gre_gcc7_0_65->SetPoint(4,30,8.32413);
    gre_gcc7_0_65->SetPointError(4,0,0.0352956);
    gre_gcc7_0_65->SetPoint(5,40,7.16134);
    gre_gcc7_0_65->SetPointError(5,0,0.030495);
    gre_gcc7_0_65->SetPoint(6,50,6.45598);
    gre_gcc7_0_65->SetPointError(6,0,0.0278967);
    gre_gcc7_0_65->SetPoint(7,60,6.00043);
    gre_gcc7_0_65->SetPointError(7,0,0.0273518);
    gre_gcc7_0_65->SetPoint(8,75,5.49203);
    gre_gcc7_0_65->SetPointError(8,0,0.0245114);
    gre_gcc7_0_65->SetPoint(9,90,5.20002);
    gre_gcc7_0_65->SetPointError(9,0,0.0232637);
    gre_gcc7_0_65->SetPoint(10,100,5.01372);
    gre_gcc7_0_65->SetPointError(10,0,0.0213542);
    gre_gcc7_0_65->SetPoint(11,150,4.45648);
    gre_gcc7_0_65->SetPointError(11,0,0.018577);
    gre_gcc7_0_65->SetPoint(12,200,4.24537);
    gre_gcc7_0_65->SetPointError(12,0,0.0179416);
    gre_gcc7_0_65->SetPoint(13,250,3.98625);
    gre_gcc7_0_65->SetPointError(13,0,0.0164062);
    gre_gcc7_0_65->SetPoint(14,400,3.52353);
    gre_gcc7_0_65->SetPointError(14,0,0.0150266);
    
    gre_gcc7_0_65->SetPoint(16,1000,3.50932);
    gre_gcc7_0_65->SetPointError(16,0,0.0143775);
    gre_gcc7_0_65->SetPoint(17,1500,3.40218);
    gre_gcc7_0_65->SetPointError(17,0,0.0379265);

    
   // TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7_0_65->SetHistogram(h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14);
    gre_gcc7_0_65->Draw("ape");
    TF1* fit_sigma_gcc7_K0L_0_65 = new TF1("fit_sigma_gcc7_K0L_0_65", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_gcc7_K0L_0_65->SetParameter(0,50.);
    fit_sigma_gcc7_K0L_0_65->SetParameter(1,0.10);
    fit_sigma_gcc7_K0L_0_65->SetParameter(2,3.0);
    fit_sigma_gcc7_K0L_0_65->SetParNames("a1_0_65","a2_0_65","const_0_65");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_65","R");
    
    TF1* fit_sigma_gcc7_K0L_0_65_only_H = new TF1("fit_sigma_gcc7_K0L_0_65_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_gcc7_K0L_0_65_only_H->SetParameter(0,50.0);
    fit_sigma_gcc7_K0L_0_65_only_H->SetParameter(1,3.0);
    fit_sigma_gcc7_K0L_0_65_only_H->SetParNames("a1_0_65_only_H","const_0_65_only_H");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_65_only_H","R");
    
    gre_gcc7->Draw("pesame");
    
    gre_gcc7 = new TGraphErrors(18);
    gre_gcc7->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_rel_E_reco_vs_E_true_0_65_to_0_80");
    gre_gcc7->SetTitle("");
    gre_gcc7->SetFillColor(1);
    gre_gcc7->SetMarkerColor(kRed-7);
    gre_gcc7->SetMarkerStyle(kOpenSquare);
    gre_gcc7->SetPoint(0,2,28.0016);
    gre_gcc7->SetPointError(0,0,0.312513);
    gre_gcc7->SetPoint(1,5,23.1505);
    gre_gcc7->SetPointError(1,0,0.229874);
    gre_gcc7->SetPoint(2,10,16.8278);
    gre_gcc7->SetPointError(2,0,0.170017);
    gre_gcc7->SetPoint(3,20,11.0449);
    gre_gcc7->SetPointError(3,0,0.117717);
    gre_gcc7->SetPoint(4,30,8.71406);
    gre_gcc7->SetPointError(4,0,0.0968127);
    gre_gcc7->SetPoint(5,40,7.64218);
    gre_gcc7->SetPointError(5,0,0.081591);
    gre_gcc7->SetPoint(6,50,7.0374);
    gre_gcc7->SetPointError(6,0,0.0790884);
    gre_gcc7->SetPoint(7,60,6.498);
    gre_gcc7->SetPointError(7,0,0.0807282);
    gre_gcc7->SetPoint(8,75,5.92743);
    gre_gcc7->SetPointError(8,0,0.0700541);
    gre_gcc7->SetPoint(9,90,5.5929);
    gre_gcc7->SetPointError(9,0,0.0694742);
    gre_gcc7->SetPoint(10,100,5.49488);
    gre_gcc7->SetPointError(10,0,0.0678995);
    gre_gcc7->SetPoint(11,150,4.83972);
    gre_gcc7->SetPointError(11,0,0.0577685);
    gre_gcc7->SetPoint(12,200,4.39405);
    gre_gcc7->SetPointError(12,0,0.0536676);
    gre_gcc7->SetPoint(13,250,3.95761);
    gre_gcc7->SetPointError(13,0,0.0515883);
    gre_gcc7->SetPoint(14,400,3.4251);
    gre_gcc7->SetPointError(14,0,0.041055);
    
    gre_gcc7->SetPoint(16,1000,3.34558);
    gre_gcc7->SetPointError(16,0,0.0471084);
    gre_gcc7->SetPoint(17,1500,3.23992);
    gre_gcc7->SetPointError(17,0,0.118434);


    
    //TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,415);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetLineColor(kRed-7);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7->SetHistogram(h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14);
    gre_gcc7->Draw("pe");
    TF1* fit_sigma_gcc7_K0L_0_65_to_0_80 = new TF1("fit_sigma_gcc7_K0L_0_65_to_0_80", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_gcc7_K0L_0_65_to_0_80->SetParameter(0,50.);
    fit_sigma_gcc7_K0L_0_65_to_0_80->SetParameter(1,0.10);
    fit_sigma_gcc7_K0L_0_65_to_0_80->SetParameter(2,3.0);
    fit_sigma_gcc7_K0L_0_65_to_0_80->SetParNames("a1_0_65_to_0_80","a2_0_65_to_0_80","const_0_65_to_0_80");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_65_to_0_80","R");
    
    TF1* fit_sigma_gcc7_K0L_0_65_to_0_80_only_H = new TF1("fit_sigma_gcc7_K0L_0_65_to_0_80_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_gcc7_K0L_0_65_to_0_80_only_H->SetParameter(0,50.0);
    fit_sigma_gcc7_K0L_0_65_to_0_80_only_H->SetParameter(1,3.0);
    fit_sigma_gcc7_K0L_0_65_to_0_80_only_H->SetParNames("a1_0_65_to_0_80_only_H","const_0_65_to_0_80_only_H");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_65_to_0_80_only_H","R");
    
    gre_gcc7->Draw("pesame");
    
    gre_gcc7 = new TGraphErrors(18);
    gre_gcc7->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_rel_E_reco_vs_E_true_0_80_to_0_94");
    gre_gcc7->SetTitle("");
    gre_gcc7->SetFillColor(1);
    gre_gcc7->SetMarkerColor(kBlue);
    gre_gcc7->SetMarkerStyle(kOpenTriangleUp);
    gre_gcc7->SetPoint(0,2,32.8076);
    gre_gcc7->SetPointError(0,0,0.472874);
    gre_gcc7->SetPoint(1,5,24.7126);
    gre_gcc7->SetPointError(1,0,0.281666);
    gre_gcc7->SetPoint(2,10,16.3954);
    gre_gcc7->SetPointError(2,0,0.170378);
    gre_gcc7->SetPoint(3,20,10.7602);
    gre_gcc7->SetPointError(3,0,0.0982656);
    gre_gcc7->SetPoint(4,30,8.7965);
    gre_gcc7->SetPointError(4,0,0.086739);
    gre_gcc7->SetPoint(5,40,8.03994);
    gre_gcc7->SetPointError(5,0,0.079595);
    gre_gcc7->SetPoint(6,50,7.24163);
    gre_gcc7->SetPointError(6,0,0.068804);
    gre_gcc7->SetPoint(7,60,6.89813);
    gre_gcc7->SetPointError(7,0,0.0668595);
    gre_gcc7->SetPoint(8,75,6.28702);
    gre_gcc7->SetPointError(8,0,0.0607709);
    gre_gcc7->SetPoint(9,90,5.63505);
    gre_gcc7->SetPointError(9,0,0.051814);
    gre_gcc7->SetPoint(10,100,5.24042);
    gre_gcc7->SetPointError(10,0,0.0499932);
    gre_gcc7->SetPoint(11,150,4.47019);
    gre_gcc7->SetPointError(11,0,0.0425189);
    gre_gcc7->SetPoint(12,200,4.10166);
    gre_gcc7->SetPointError(12,0,0.0389849);
    gre_gcc7->SetPoint(13,250,3.79353);
    gre_gcc7->SetPointError(13,0,0.0367435);
    gre_gcc7->SetPoint(14,400,3.45924);
    gre_gcc7->SetPointError(14,0,0.0309687);
    
    gre_gcc7->SetPoint(16,1000,3.39242);
    gre_gcc7->SetPointError(16,0,0.0321176);
    gre_gcc7->SetPoint(17,1500,3.26134);
    gre_gcc7->SetPointError(17,0,0.0911823);

    
    //TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMinimum(4.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMaximum(17.5);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7->SetHistogram(h_K0L_gcc7_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14);
    gre_gcc7->Draw("pe");
    TF1* fit_sigma_gcc7_K0L_0_80_to_0_92 = new TF1("fit_sigma_gcc7_K0L_0_80_to_0_92", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+pow([2],2))",6,85);
    fit_sigma_gcc7_K0L_0_80_to_0_92->SetParameter(0,50.0);
    fit_sigma_gcc7_K0L_0_80_to_0_92->SetParameter(1,0.10);
    fit_sigma_gcc7_K0L_0_80_to_0_92->SetParameter(2,3.0);
    fit_sigma_gcc7_K0L_0_80_to_0_92->SetParNames("a1_0_80_to_0_92","a2_0_80_to_0_92","const_0_80_to_0_92");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_80_to_0_92","R");
    
    TF1* fit_sigma_gcc7_K0L_0_80_to_0_92_only_H = new TF1("fit_sigma_gcc7_K0L_0_80_to_0_92_only_H", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))",6,415);
    fit_sigma_gcc7_K0L_0_80_to_0_92_only_H->SetParameter(0,50.0);
    fit_sigma_gcc7_K0L_0_80_to_0_92_only_H->SetParameter(1,3.0);
    fit_sigma_gcc7_K0L_0_80_to_0_92_only_H->SetParNames("a1_0_80_to_0_92_only_H","const_0_80_to_0_92_only_H");
    gre_gcc7->Fit("fit_sigma_gcc7_K0L_0_80_to_0_92_only_H","R");
    gre_gcc7->Draw("pesame");
    

    
    
    TLegend *leg_gcc7_K0L_E_reco_vs_E_true = new TLegend(0.60,0.60,0.94,0.91,"K0L","brNDC");
    leg_gcc7_K0L_E_reco_vs_E_true->SetBorderSize(1);
    leg_gcc7_K0L_E_reco_vs_E_true->SetTextSize(0.027);
    leg_gcc7_K0L_E_reco_vs_E_true->SetLineColor(1);
    leg_gcc7_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_gcc7_K0L_E_reco_vs_E_true->SetLineWidth(1);
    leg_gcc7_K0L_E_reco_vs_E_true->SetFillColor(0);
    leg_gcc7_K0L_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_gcc7_K0L_E_reco_vs_E_true=leg_gcc7_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_gcc7_K0L_E_reco_vs_E_true=leg_gcc7_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","0.65<|cos#theta_{K0L}^{true}|<0.80","pe");
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_gcc7_K0L_E_reco_vs_E_true=leg_gcc7_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_E_reco_vs_E_true_0_80_to_094_CLIC_o3_v14","0.80<|cos#theta_{K0L}^{true}|<0.94","pe");
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_gcc7_K0L_E_reco_vs_E_true->Draw();
    
    //now for mean
    
    //=========Macro generated from canvas: resolution_gcc7_GraphCanvas/
    //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
    TCanvas *resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L = new TCanvas("resolution_gcc7_mean_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc7_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFillColor(0);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderMode(0);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderSize(2);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridx();
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridy();
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetRightMargin(0.0172);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetTopMargin(0.055);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBottomMargin(0.138);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    resolution_gcc7_mean_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    
    
    //points 2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
    TGraphErrors *gre_gcc7_mean = new TGraphErrors(18);
    gre_gcc7_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_mean_rel_E_reco_vs_E_true_0_60");
    gre_gcc7_mean->SetTitle("");
    gre_gcc7_mean->SetFillColor(1);
    gre_gcc7_mean->SetMarkerColor(kBlack);
    gre_gcc7_mean->SetMarkerStyle(kOpenCircle);
    gre_gcc7_mean->SetPoint(0,2,0.827197);
    gre_gcc7_mean->SetPointError(0,0,0.00224196);
    gre_gcc7_mean->SetPoint(1,5,0.950192);
    gre_gcc7_mean->SetPointError(1,0,0.00168022);
    gre_gcc7_mean->SetPoint(2,10,1.02976);
    gre_gcc7_mean->SetPointError(2,0,0.0010693);
    gre_gcc7_mean->SetPoint(3,20,1.03632);
    gre_gcc7_mean->SetPointError(3,0,0.000609081);
    gre_gcc7_mean->SetPoint(4,30,1.01891);
    gre_gcc7_mean->SetPointError(4,0,0.00048039);
    gre_gcc7_mean->SetPoint(5,40,1.00686);
    gre_gcc7_mean->SetPointError(5,0,0.000412386);
    gre_gcc7_mean->SetPoint(6,50,1.00152);
    gre_gcc7_mean->SetPointError(6,0,0.00037424);
    gre_gcc7_mean->SetPoint(7,60,1.00167);
    gre_gcc7_mean->SetPointError(7,0,0.000349307);
    gre_gcc7_mean->SetPoint(8,75,1.00411);
    gre_gcc7_mean->SetPointError(8,0,0.000315247);
    gre_gcc7_mean->SetPoint(9,90,1.00446);
    gre_gcc7_mean->SetPointError(9,0,0.000300152);
    gre_gcc7_mean->SetPoint(10,100,1.00475);
    gre_gcc7_mean->SetPointError(10,0,0.000288671);
    gre_gcc7_mean->SetPoint(11,150,1.00259);
    gre_gcc7_mean->SetPointError(11,0,0.000260141);
    gre_gcc7_mean->SetPoint(12,200,1.0025);
    gre_gcc7_mean->SetPointError(12,0,0.000253453);
    gre_gcc7_mean->SetPoint(13,250,1.00439);
    gre_gcc7_mean->SetPointError(13,0,0.000252411);
    gre_gcc7_mean->SetPoint(14,400,0.999447);
    gre_gcc7_mean->SetPointError(14,0,0.000214818);
    
    gre_gcc7_mean->SetPoint(16,1000,0.985535);
    gre_gcc7_mean->SetPointError(16,0,0.000224705);
    gre_gcc7_mean->SetPoint(17,1500,0.989918);
    gre_gcc7_mean->SetPointError(17,0,0.000576204);
    
    
    
    //TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7_mean->SetHistogram(h_K0L_gcc7_mean_E_reco_vs_E_true_0_60_CLIC_o3_v14);
    gre_gcc7_mean->Draw("ape");
    
    gre_gcc7_mean = new TGraphErrors(18);
    gre_gcc7_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_mean_rel_E_reco_vs_E_true_0_65");
    gre_gcc7_mean->SetTitle("");
    gre_gcc7_mean->SetFillColor(1);
    gre_gcc7_mean->SetMarkerColor(kBlack);
    gre_gcc7_mean->SetMarkerStyle(kOpenCircle);
    gre_gcc7_mean->SetPoint(0,2,0.823476);
    gre_gcc7_mean->SetPointError(0,0,0.00213594);
    gre_gcc7_mean->SetPoint(1,5,0.94779);
    gre_gcc7_mean->SetPointError(1,0,0.00160805);
    gre_gcc7_mean->SetPoint(3,20,1.0341);
    gre_gcc7_mean->SetPointError(3,0,0.000589192);
    gre_gcc7_mean->SetPoint(4,30,1.01745);
    gre_gcc7_mean->SetPointError(4,0,0.000464955);
    gre_gcc7_mean->SetPoint(5,40,1.0055);
    gre_gcc7_mean->SetPointError(5,0,0.000398643);
    gre_gcc7_mean->SetPoint(6,50,1.00031);
    gre_gcc7_mean->SetPointError(6,0,0.000362429);
    gre_gcc7_mean->SetPoint(7,60,1.00069);
    gre_gcc7_mean->SetPointError(7,0,0.000339381);
    gre_gcc7_mean->SetPoint(8,75,1.00319);
    gre_gcc7_mean->SetPointError(8,0,0.000306024);
    gre_gcc7_mean->SetPoint(9,90,1.00378);
    gre_gcc7_mean->SetPointError(9,0,0.000290035);
    gre_gcc7_mean->SetPoint(10,100,1.00413);
    gre_gcc7_mean->SetPointError(10,0,0.000279997);
    gre_gcc7_mean->SetPoint(11,150,1.00212);
    gre_gcc7_mean->SetPointError(11,0,0.000251579);
    gre_gcc7_mean->SetPoint(12,200,1.002);
    gre_gcc7_mean->SetPointError(12,0,0.000245791);
    gre_gcc7_mean->SetPoint(13,250,1.00405);
    gre_gcc7_mean->SetPointError(13,0,0.000238177);
    gre_gcc7_mean->SetPoint(14,400,0.999157);
    gre_gcc7_mean->SetPointError(14,0,0.000206509);
    
    gre_gcc7_mean->SetPoint(16,1000,0.985578);
    gre_gcc7_mean->SetPointError(16,0,0.000215609);
    gre_gcc7_mean->SetPoint(17,1500,0.990115);
    gre_gcc7_mean->SetPointError(17,0,0.000558356);

    
    //TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetLineColor(kBlack);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7_mean->SetHistogram(h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14);
    gre_gcc7_mean->Draw("ape");
    
    gre_gcc7_mean = new TGraphErrors(18);
    gre_gcc7_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_mean_rel_E_reco_vs_E_true_0_65_to_0_80");
    gre_gcc7_mean->SetTitle("");
    gre_gcc7_mean->SetFillColor(1);
    gre_gcc7_mean->SetMarkerColor(kBlue);
    gre_gcc7_mean->SetMarkerStyle(kOpenSquare);
    gre_gcc7_mean->SetPoint(0,2,0.759228);
    gre_gcc7_mean->SetPointError(0,0,0.00418471);
    gre_gcc7_mean->SetPoint(1,5,0.87814);
    gre_gcc7_mean->SetPointError(1,0,0.00326874);
    gre_gcc7_mean->SetPoint(2,10,0.967817);
    gre_gcc7_mean->SetPointError(2,0,0.00228689);
    gre_gcc7_mean->SetPoint(3,20,0.996529);
    gre_gcc7_mean->SetPointError(3,0,0.00151895);
    gre_gcc7_mean->SetPoint(4,30,0.991879);
    gre_gcc7_mean->SetPointError(4,0,0.001257);
    gre_gcc7_mean->SetPoint(5,40,0.988838);
    gre_gcc7_mean->SetPointError(5,0,0.0011573);
    gre_gcc7_mean->SetPoint(6,50,0.987606);
    gre_gcc7_mean->SetPointError(6,0,0.00105337);
    gre_gcc7_mean->SetPoint(7,60,0.990945);
    gre_gcc7_mean->SetPointError(7,0,0.00100961);
    gre_gcc7_mean->SetPoint(8,75,0.998912);
    gre_gcc7_mean->SetPointError(8,0,0.000919534);
    gre_gcc7_mean->SetPoint(9,90,1.00138);
    gre_gcc7_mean->SetPointError(9,0,0.000865912);
    gre_gcc7_mean->SetPoint(10,100,1.00365);
    gre_gcc7_mean->SetPointError(10,0,0.000845955);
    gre_gcc7_mean->SetPoint(11,150,1.00574);
    gre_gcc7_mean->SetPointError(11,0,0.000752807);
    gre_gcc7_mean->SetPoint(12,200,1.00648);
    gre_gcc7_mean->SetPointError(12,0,0.000714594);
    gre_gcc7_mean->SetPoint(13,250,1.00887);
    gre_gcc7_mean->SetPointError(13,0,0.000648246);
    gre_gcc7_mean->SetPoint(14,400,1.00325);
    gre_gcc7_mean->SetPointError(14,0,0.000566569);
    
    gre_gcc7_mean->SetPoint(16,1000,0.990667);
    gre_gcc7_mean->SetPointError(16,0,0.000565527);
    gre_gcc7_mean->SetPoint(17,1500,0.991918);
    gre_gcc7_mean->SetPointError(17,0,0.00154058);


    
   // TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7_mean->SetHistogram(h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14);
    gre_gcc7_mean->Draw("pe");
    
    gre_gcc7_mean = new TGraphErrors(18);
    gre_gcc7_mean->SetName("CLIC_o3_v14_res_graph_K0L_gcc7_mean_rel_E_reco_vs_E_true_0_80_to_0_94");
    gre_gcc7_mean->SetTitle("");
    gre_gcc7_mean->SetFillColor(1);
    gre_gcc7_mean->SetMarkerColor(kRed-7);
    gre_gcc7_mean->SetMarkerStyle(kOpenTriangleUp);
    gre_gcc7_mean->SetPoint(0,2,0.85404);
    gre_gcc7_mean->SetPointError(0,0,0.00543767);
    gre_gcc7_mean->SetPoint(1,5,0.984003);
    gre_gcc7_mean->SetPointError(1,0,0.00380124);
    gre_gcc7_mean->SetPoint(2,10,1.05078);
    gre_gcc7_mean->SetPointError(2,0,0.00213107);
    gre_gcc7_mean->SetPoint(3,20,1.04919);
    gre_gcc7_mean->SetPointError(3,0,0.00133534);
    gre_gcc7_mean->SetPoint(4,30,1.02996);
    gre_gcc7_mean->SetPointError(4,0,0.00109388);
    gre_gcc7_mean->SetPoint(5,40,1.0141);
    gre_gcc7_mean->SetPointError(5,0,0.000989474);
    gre_gcc7_mean->SetPoint(6,50,1.00934);
    gre_gcc7_mean->SetPointError(6,0,0.000891629);
    gre_gcc7_mean->SetPoint(7,60,1.00662);
    gre_gcc7_mean->SetPointError(7,0,0.000868459);
    gre_gcc7_mean->SetPoint(8,75,1.00923);
    gre_gcc7_mean->SetPointError(8,0,0.000779802);
    gre_gcc7_mean->SetPoint(9,90,1.00794);
    gre_gcc7_mean->SetPointError(9,0,0.000693572);
    gre_gcc7_mean->SetPoint(10,100,1.00861);
    gre_gcc7_mean->SetPointError(10,0,0.000649895);
    gre_gcc7_mean->SetPoint(11,150,1.00342);
    gre_gcc7_mean->SetPointError(11,0,0.0005485);
    gre_gcc7_mean->SetPoint(12,200,0.999522);
    gre_gcc7_mean->SetPointError(12,0,0.000501063);
    gre_gcc7_mean->SetPoint(13,250,0.999122);
    gre_gcc7_mean->SetPointError(13,0,0.000462173);
    gre_gcc7_mean->SetPoint(14,400,0.987724);
    gre_gcc7_mean->SetPointError(14,0,0.000416037);
    
    gre_gcc7_mean->SetPoint(16,1000,0.965478);
    gre_gcc7_mean->SetPointError(16,0,0.000408232);
    gre_gcc7_mean->SetPoint(17,1500,0.966004);
    gre_gcc7_mean->SetPointError(17,0,0.00105072);
    
    
    //TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,1649.9);
    TH1F *h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14 = new TH1F("h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14","",100,0,415.);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMinimum(0.9);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetMaximum(1.1);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetDirectory(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetStats(0);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->SetLineColor(kBlue);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetXaxis()->SetTitleOffset(1.18);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitle("E_{PFO}^{reco}/E_{true}^{K0L}");
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetLabelSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleSize(0.05);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleOffset(1.20);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetYaxis()->SetTitleFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelFont(42);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetLabelSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleSize(0.035);
    h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14->GetZaxis()->SetTitleFont(42);
    gre_gcc7_mean->SetHistogram(h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v14);
    gre_gcc7_mean->Draw("pe");
    
    TLegend *leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true = new TLegend(0.60,0.20,0.94,0.60,"K^{0}_{L}","brNDC");
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetBorderSize(1);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetTextSize(0.027);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineColor(1);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineWidth(1);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetFillColor(0);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true=leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65","pe");
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true=leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v14","0.65<|cos#theta_{K0L}^{true}|<0.80","pe");
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true=leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_mean_E_reco_vs_E_true_0_80_to_094_CLIC_o3_v14","0.80<|cos#theta_{K0L}^{true}|<0.94","pe");
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineColor(kBlue);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerColor(kBlue);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerStyle(kOpenTriangleUp);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->SetTextFont(42);
    leg_gcc7_K0L_gcc7_mean_E_reco_vs_E_true->Draw();
                                          
        
    TCanvas *resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L = new TCanvas("resolution_gcc7_vs_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc7_vs_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",0,0,800,700);
    gStyle->SetOptStat(0);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->Range(-186.894,-0.873515,1682.046,6.114605);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFillColor(0);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderMode(0);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBorderSize(2);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridx();
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetGridy();
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetRightMargin(0.0172);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetTopMargin(0.055);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetBottomMargin(0.138);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);
    resolution_gcc7_vs_gcc62_GraphCanvas_CLIC_n_K0L_res_rel_K0L->SetFrameBorderMode(0);

    //h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->Draw("p,e,a");
    //gre_gcc62_mean_0_60->Draw("ape");
    //gre_gcc62_0_65->Add()
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMarkerColor(kRed-7);
    h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14->SetMarkerStyle(kOpenSquare);
    gre_gcc62_0_65->SetMarkerColor(kRed-7);
    gre_gcc62_0_65->SetMarkerStyle(kOpenSquare);
    gre_gcc62_0_65->SetMinimum(3.0);
    gre_gcc62_0_65->SetMaximum(17.5);
    //gre_gcc62_0_65->Add(h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14);
    gre_gcc62_0_65->Draw("ape");
    gre_gcc7_0_65->Draw("pe,same");
    //h_K0L_gcc7_mean_E_reco_vs_E_true_0_65_CLIC_o3_v14->Draw("p,e");


    TLegend *leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true = new TLegend(0.60,0.60,0.94,0.91,"K0L","brNDC");
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetBorderSize(1);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetTextSize(0.027);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineColor(1);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineWidth(1);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetFillColor(0);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetFillStyle(1001);
    TLegendEntry *leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true=leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc7_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65,gcc7,10.03.p03","pe");
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineColor(kBlack);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerColor(kBlack);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenCircle);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true=leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->AddEntry("h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v14","|cos#theta_{K0L}^{true}|<0.65,gcc62,10.02.p02","pe");
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineColor(kRed-7);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineStyle(1);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetLineWidth(2);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerColor(kRed-7);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerStyle(kOpenSquare);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetMarkerSize(1.5);
    leg_entry_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->SetTextFont(42);
    leg_gcc7_vs_gcc62_K0L_E_reco_vs_E_true->Draw();   
    */
                                 
}
