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
#include "TEfficiency.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TColor.h"
#include <vector>
#include "TGraphAsymmErrors.h"


Color_t color1 = kBlack;
Color_t color2 = kRed-7;
Color_t color3 = kBlue;
Color_t color4 = kGreen+2;
Color_t color5 = kCyan+1;
Color_t color6 = kRed+2;
Color_t color7 = kOrange;
Color_t color8 = kViolet+2;


void CLICdpLabel(std::string status);
void CLICdpLabel(Double_t x=0.20, Double_t y=0.86, std::string status="work in progress", Color_t color=kBlack);

void CLICdpStyle()
{
  gROOT->SetStyle("Plain"); /*Default white background for all plots*/

  TStyle *clicdpStyle = gROOT->GetStyle("Plain");

  /* set bkg color of all to kWhite: white, but not 0*/
  clicdpStyle->SetCanvasColor(kWhite);
  clicdpStyle->SetFrameFillColor(kWhite);
  clicdpStyle->SetStatColor(kWhite);
  clicdpStyle->SetPadColor(kWhite);
  clicdpStyle->SetFillColor(10);
  clicdpStyle->SetTitleFillColor(kWhite);
  
  
   /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
   clicdpStyle->SetPaperSize(20, 26); 
   /* No yellow border around histogram*/
   clicdpStyle->SetDrawBorder(0);
   /* remove border of canvas*/
   clicdpStyle->SetCanvasBorderMode(0);
   /* remove border of pads*/
   clicdpStyle->SetPadBorderMode(0);
   clicdpStyle->SetFrameBorderMode(0);
   clicdpStyle->SetLegendBorderSize(0);
  
   /* default text size*/
   clicdpStyle->SetTextSize(0.05);
   clicdpStyle->SetTitleSize(0.07,"xyz");
   clicdpStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   clicdpStyle->SetLabelOffset(0.015,"xyz");
   clicdpStyle->SetTitleOffset(1.2,"yz"); //equivalent to: clicdpStyle->SetTitleYOffset(1.2);
   clicdpStyle->SetTitleOffset(1.0,"x");



   /* Use visible font for all text*/
   int font = 42; 
   clicdpStyle->SetTitleFont(font);
   clicdpStyle->SetTitleFontSize(0.06);
   clicdpStyle->SetStatFont(font);
   clicdpStyle->SetStatFontSize(0.07);
   clicdpStyle->SetTextFont(font);
   clicdpStyle->SetLabelFont(font,"xyz");
   clicdpStyle->SetTitleFont(font,"xyz");
   clicdpStyle->SetTitleBorderSize(0);
   clicdpStyle->SetStatBorderSize(1);
   //clicdpStyle->SetLegendFont(font);

   /* big marker points*/
   clicdpStyle->SetMarkerStyle(1);
   clicdpStyle->SetLineWidth(2);  
   clicdpStyle->SetMarkerSize(1.2);
   /*set palette in 2d histogram to nice and colorful one*/
   clicdpStyle->SetPalette(1,0); 

   /*No title for histograms*/
   clicdpStyle->SetOptTitle(0);
   /* show the errors on the stat box */
   clicdpStyle->SetOptStat(0); 
   /* show errors on fitted parameters*/
   clicdpStyle->SetOptFit(0); 
   /* number of decimals used for errors*/
   clicdpStyle->SetEndErrorSize(5);   

   /* set line width to 2 by default so that histograms are visible when printed small
      idea: emphasize the data, not the frame around*/
   clicdpStyle->SetHistLineWidth(2);
   clicdpStyle->SetFrameLineWidth(2);
   clicdpStyle->SetFuncWidth(2);
   clicdpStyle->SetHistLineColor(kBlack);
   clicdpStyle->SetFuncColor(kRed);
   clicdpStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   clicdpStyle->SetPadBottomMargin(0.18);
   clicdpStyle->SetPadTopMargin(0.08);
   clicdpStyle->SetPadRightMargin(0.08);
   clicdpStyle->SetPadLeftMargin(0.17);
   
   //set the number of divisions to show
   clicdpStyle->SetNdivisions(506, "xy");
   
   //turn off xy grids
   clicdpStyle->SetPadGridX(0);
   clicdpStyle->SetPadGridY(0);
   
   //set the tick mark style
   clicdpStyle->SetPadTickX(1);
   clicdpStyle->SetPadTickY(1);

   clicdpStyle->SetCanvasDefW(800);
   clicdpStyle->SetCanvasDefH(700);

   //gROOT->ForceStyle();
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 2 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw2Legend(TH1 *histo1, 
                 TH1 *histo2,
                 const Char_t *label1, 
                 const Char_t *label2,
                 const Char_t *header="")
{

  Float_t max1 = histo1->GetMaximum();
  Float_t max2 = histo2->GetMaximum();
  if (max1 >= max2) histo1->SetMaximum(max1 * 1.3);
  else histo1->SetMaximum(max2 * 1.3);

  TLegend *legend = new TLegend(0.1436782,0.8072034,0.6408046,0.9851695,NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.06);
  legend->SetTextFont(42);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LPF");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"LPF");
  entry2->SetTextColor(histo2->GetLineColor());

  if (std::string(header) != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw();

  gPad->Update();
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 7 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw7Legend(TGraph *histo1, 
                 TGraph *histo2,
                 TGraph *histo3,
                 TGraph *histo4,
                 TGraph *histo5,
                 TGraph *histo6,
                 TGraph *histo7,
		 const Char_t *label1, 
                 const Char_t *label2,
		 const Char_t *label3, 
                 const Char_t *label4,
		 const Char_t *label5, 
                 const Char_t *label6,
		 const Char_t *label7, 
                 const Char_t *header="")
{

  TLegend *legend = new TLegend(0.006231824,0.3648567,0.8683008,0.9190962,NULL,"brNDC");
  legend->SetTextAlign(12);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.1);
  legend->SetTextFont(42);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"L");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"L");
  entry2->SetTextColor(histo2->GetLineColor());

  TLegendEntry* entry3 = legend->AddEntry(histo3,label3,"L");
  entry3->SetTextColor(histo3->GetLineColor());

  TLegendEntry* entry4 = legend->AddEntry(histo4,label4,"L");
  entry4->SetTextColor(histo4->GetLineColor());

  TLegendEntry* entry5 = legend->AddEntry(histo5,label5,"L");
  entry5->SetTextColor(histo5->GetLineColor());

  TLegendEntry* entry6 = legend->AddEntry(histo6,label6,"L");
  entry6->SetTextColor(histo6->GetLineColor());

  TLegendEntry* entry7 = legend->AddEntry(histo7,label7,"L");
  entry7->SetTextColor(histo7->GetLineColor());

  if (std::string(header) != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw();

  gPad->Update();
}

/***********************************************************************/
/*                                                                     */
/*       draw statistics boxes for 2 histos                            */
/*       (needs h2->Draw("sames"))                                     */
/*                                                                     */
/***********************************************************************/
bool Draw2StatsBoxes(TH1 *histo1, TH1 *histo2)
{
  /*Get the stats box and set the right colors*/
  TPaveStats *statsh1 = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  if (statsh1 == NULL)
    {
      std::cout<<"\n Could not find statistics box for histo " << histo1->GetName()<<std::endl;
      return false;
    }
  statsh1->SetLineColor(histo1->GetLineColor());
  statsh1->SetTextColor(histo1->GetLineColor());
  statsh1->SetX1NDC(0.72);
  statsh1->SetY1NDC(0.68);
  statsh1->SetX2NDC(0.92);
  statsh1->SetY2NDC(0.92);

  TPaveStats *statsh2 = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  if (statsh2 == NULL)
    {
      std::cout<<"\n Could not find statistics box for histo " << histo2->GetName()<<std::endl;
      std::cout<<"You need to call: histo->Draw(\"sames\")"<<std::endl;
      return false;
    }
  statsh2->SetLineColor(histo2->GetLineColor());
  statsh2->SetTextColor(histo2->GetLineColor());

  /*need to move the stats box below the previous one, as by default they are stacked*/
  statsh2->SetX1NDC(0.72);
  statsh2->SetY1NDC(0.44);
  statsh2->SetX2NDC(0.92);
  statsh2->SetY2NDC(0.67);
  statsh2->Draw();
  
  return true;
}


void CLICdpLabel(std::string status){
  CLICdpLabel(0.20,0.86,status,kBlack);
}//end CLICdpLabel



// void CLICdpLabel(Double_t x=0.20, Double_t y=0.86, std::string status="work in progress", Color_t color=kBlack){

//   TLatex l; 
//   l.SetNDC();
//   l.SetTextFont(42);
//   l.SetTextColor(color);
//   l.SetTextSize(0.05);

//   std::string label = std::string("CLICdp ");
//   l.DrawLatex(x,y,label.c_str());

//   TLatex l2; 
//   l2.SetNDC();
//   l2.SetTextFont(42);
//   l2.SetTextColor(color);
//   l2.SetTextSize(0.035);

//   double dy = l2.GetTextSize()+0.005;

//   l2.DrawLatex(x,y-dy,status.c_str());

//   return;
// }//end CLICdpLabel


void CLICdpLabel(Double_t x, Double_t y, std::string status, Color_t color){

  TLatex l; 
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(kBlack);
  l.SetTextSize(0.05);

  std::string label = std::string("CLICdp ");
  l.DrawLatex(x,y,label.c_str());

  TLatex l2; 
  l2.SetNDC();
  l2.SetTextFont(42);
  l2.SetTextColor(color);
  l2.SetTextSize(0.035);

  //double dx = l2.GetTextSize()+0.005;
  double dx = 0.4*l.GetTextSize()*label.length();

  l2.DrawLatex(x+dx,y,status.c_str());

  return;
}//end CLICdpLabel


TCanvas* setUpperCanvas(const char* canvas_name) {
  TCanvas* c1= new TCanvas(canvas_name,canvas_name,0.,0.,800.,700.);
  c1->cd();
  gStyle->SetOptStat(0);
  //c1->Range(-186.894,-0.873515,1682.046,6.114605);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  //c1->SetGridx();
  //c1->SetGridy();
  //c1->SetRightMargin(0.0172);
  //c1->SetTopMargin(0.055);
  //c1->SetBottomMargin(0.138);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  return c1;
}

void drawFancyPlots(){

  CLICdpStyle();

  //TLegend *leg_label_CLICdp = new TLegend(0.20,0.90,0.40,0.94, "CLICdp");
  //lewwlabel_CLICdp->SetBorderSize(0);
 
  TLatex* l = new TLatex();
  l->SetNDC();
  l->SetTextFont(42);
  l->SetTextColor(kBlack);
  l->SetTextSize(0.045);
  std::string label("CLICdp");
  double x=0.20, y=0.93;
  
  TText* leg_label_CLICdp= new TText();
  leg_label_CLICdp->SetTextSize(0.045);
  leg_label_CLICdp->SetTextFont(42);

  //380 overlay sample defined below
  //TFile* file_histos_normalRange_jets=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_180810_and_180724_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER100_range10_WW_ZZ_mass_histos_pyt380_MeanReal.root");

  //TFile* file_histos_normalRange_jets=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER100_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC.root");
  //TFile* file_histos_normalRange_jets=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181101_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi300_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat_allSCpoints.root");
  TFile* file_histos_normalRange_jets=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi500_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi300_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat_allSCpoints_realSamples.root");
  //allSCpoints.root


  ///afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");
//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_180724_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBins75DPhiDTheta_range10_WW_ZZ_mass_histos_pyt380.root");
//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx.root");
//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER100_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC.root");



  //380 overlay sample defined below
  TFile* file_histos_normalRange_jets_SCStudies=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi500_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

  ///afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx.root");

//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER100_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC.root");

  //TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_180810_and_180724_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER100_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC.root");

  //do DrawManyBinsJER200ddphi500 for dphi plots and CB fits, for region of 250 and dphi vs E plot use DrawManyBinsJER200ddphi300
TFile* file_histos_normalRange_jets_drawbins=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi500_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi300_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

  //afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi300_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");
//afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx.root");

  //TFile* file_histos_normalRange_jets_drawbins=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx.root");
  //TFile* file_histos_normalRange_jets_drawbins=TFile::Open("afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_180724_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBins200DPhiDTheta_range10_WW_ZZ_mass_histos_pyt380.root");

  TH1F* h_Zuds100_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_100");
 //h_Zuds100_JER_DR07->SetLineColor(kRed);
  h_Zuds100_JER_DR07->SetLineStyle(1);
  h_Zuds100_JER_DR07->SetLineWidth(2);
  h_Zuds100_JER_DR07->SetMinimum(1.75);
  h_Zuds100_JER_DR07->SetMaximum(10.1);
  h_Zuds100_JER_DR07->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_200");
 //h_Zuds100_JER_DR07->SetLineColor(kRed);
  h_Zuds200_JER_DR07->SetLineStyle(1);
  h_Zuds200_JER_DR07->SetLineWidth(2);
  h_Zuds200_JER_DR07->SetMinimum(1.75);
  h_Zuds200_JER_DR07->SetMaximum(10.1);
  TH1F* h_Zuds380_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_380");
 //h_Zuds380_JER_DR07->SetLineColor(kRed);
  h_Zuds380_JER_DR07->SetLineStyle(1);
  h_Zuds380_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_500");
  h_Zuds500_JER_DR07->SetLineStyle(1);
  h_Zuds500_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_1500");
  h_Zuds1500_JER_DR07->SetLineStyle(1);
  h_Zuds1500_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR07_3000");
  h_Zuds3000_JER_DR07->SetLineStyle(1);
  h_Zuds3000_JER_DR07->SetLineWidth(2);

  TCanvas *resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy = setUpperCanvas("resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy");
  //resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy->cd();
  //TLegend *leg_JER_DR07_RMS90_CT_FullSummary = resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_DR07_RMS90_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_DR07_RMS90_CT_FullSummary->SetBorderSize(0);
  leg_JER_DR07_RMS90_CT_FullSummary->SetTextAlign(12);
  leg_JER_DR07_RMS90_CT_FullSummary->SetTextSize(0.050);
  leg_JER_DR07_RMS90_CT_FullSummary->SetTextFont(42);
  leg_JER_DR07_RMS90_CT_FullSummary->SetMargin(0.15);
  leg_JER_DR07_RMS90_CT_FullSummary->SetLineColor(1);
  leg_JER_DR07_RMS90_CT_FullSummary->SetLineStyle(1);
  leg_JER_DR07_RMS90_CT_FullSummary->SetLineWidth(1);
  leg_JER_DR07_RMS90_CT_FullSummary->SetFillColor(0);
  leg_JER_DR07_RMS90_CT_FullSummary->SetFillStyle(1001);
  leg_JER_DR07_RMS90_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds100_JER_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds200_JER_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds380_JER_DR07->DrawCopy("h,e,same"),"#approx 190 GeV Jets");
  leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds500_JER_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds1500_JER_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds3000_JER_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  h_Zuds100_JER_DR07->SetMaximum(15.0);

  TCanvas *resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy_zoomedOut = setUpperCanvas("resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy_zoomedOut");
  //resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy->cd();
  //TLegend *leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut = resolutionGraphCanvas_JER_DR07_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetBorderSize(0);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetTextAlign(12);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetTextSize(0.050);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetTextFont(42);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetMargin(0.15);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetLineColor(1);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetLineStyle(1);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetLineWidth(1);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetFillColor(0);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetFillStyle(1001);
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->SetHeader("VLC7 Jets");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds100_JER_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds200_JER_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds380_JER_DR07->DrawCopy("h,e,same"),"#approx 190 GeV Jets");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds500_JER_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds1500_JER_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->AddEntry(h_Zuds3000_JER_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_JER_DR07_RMS90_CT_FullSummary_zoomedOut->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds100_JER_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR05_100");
  h_Zuds100_JER_DR05->SetLineColor(kBlue);
  h_Zuds100_JER_DR05->SetLineStyle(1);
  h_Zuds100_JER_DR05->SetMinimum(1.750);
  h_Zuds100_JER_DR05->SetMaximum(10.1);
  h_Zuds100_JER_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_DR05->GetXaxis()->SetTitle("|cos#theta|");

  TH1F* h_Zuds100_JER_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR06_100");
  h_Zuds100_JER_DR06->SetLineColor(kRed);
  h_Zuds100_JER_DR06->SetLineStyle(1);
  h_Zuds100_JER_DR06->SetLineWidth(2);

  h_Zuds100_JER_DR07->SetLineColor(kAzure+7);
  h_Zuds100_JER_DR07->SetLineStyle(1);
  h_Zuds100_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds100_JER_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR08_100");
  h_Zuds100_JER_DR08->SetLineColor(kViolet);
  h_Zuds100_JER_DR08->SetLineStyle(1);
  h_Zuds100_JER_DR08->SetLineWidth(2);
  TH1F* h_Zuds100_JER_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR09_100");
  h_Zuds100_JER_DR09->SetLineColor(kOrange);
  h_Zuds100_JER_DR09->SetLineStyle(1);
  h_Zuds100_JER_DR09->SetLineWidth(2);
  TH1F* h_Zuds100_JER_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR10_100");
  h_Zuds100_JER_DR10->SetLineColor(kYellow+1);
  h_Zuds100_JER_DR10->SetLineStyle(1);
  h_Zuds100_JER_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds100_JetDR05_10 = setUpperCanvas("canvas_Zuds100_JetDR05_10");
  
  TLegend* leg_Zuds100_JetDR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds100_JetDR05_10->SetBorderSize(0);
  leg_Zuds100_JetDR05_10->SetTextAlign(12);
  leg_Zuds100_JetDR05_10->SetTextSize(0.050);
  leg_Zuds100_JetDR05_10->SetTextFont(42);
  leg_Zuds100_JetDR05_10->SetMargin(0.15);
  leg_Zuds100_JetDR05_10->SetLineColor(1);
  leg_Zuds100_JetDR05_10->SetLineStyle(1);
  leg_Zuds100_JetDR05_10->SetLineWidth(1);
  leg_Zuds100_JetDR05_10->SetFillColor(0);
  //leg_Zuds100_JetDR05_10->SetFillStyle(1001);
  leg_Zuds100_JetDR05_10->SetFillStyle(0);
  leg_Zuds100_JetDR05_10->SetHeader("Z#rightarrow uds, 91 GeV, #approx 45.5 GeV Jets");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds100_JetDR05_10->AddEntry(h_Zuds100_JER_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds100_JetDR05_10->Draw();


  TH1F* h_Zuds500_JER_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR05_500");
  h_Zuds500_JER_DR05->SetLineColor(kBlue);
  h_Zuds500_JER_DR05->SetLineStyle(1);
  h_Zuds500_JER_DR05->SetMinimum(1.75);
  h_Zuds500_JER_DR05->SetMaximum(10.1);
  h_Zuds500_JER_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds500_JER_DR05->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds500_JER_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR06_500");
  h_Zuds500_JER_DR06->SetLineColor(kRed);
  h_Zuds500_JER_DR06->SetLineStyle(1);
  h_Zuds500_JER_DR06->SetLineWidth(2);

  h_Zuds500_JER_DR07->SetLineColor(kAzure+7);
  h_Zuds500_JER_DR07->SetLineStyle(1);
  h_Zuds500_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR08_500");
  h_Zuds500_JER_DR08->SetLineColor(kViolet);
  h_Zuds500_JER_DR08->SetLineStyle(1);
  h_Zuds500_JER_DR08->SetLineWidth(2);
  TH1F* h_Zuds500_JER_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR09_500");
  h_Zuds500_JER_DR09->SetLineColor(kOrange);
  h_Zuds500_JER_DR09->SetLineStyle(1);
  h_Zuds500_JER_DR09->SetLineWidth(2);
  TH1F* h_Zuds500_JER_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_DR10_500");
  h_Zuds500_JER_DR10->SetLineColor(kYellow+1);
  h_Zuds500_JER_DR10->SetLineStyle(1);
  h_Zuds500_JER_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds500_JetDR05_10 = setUpperCanvas("canvas_Zuds500_JetDR05_10");
  
  TLegend* leg_Zuds500_JetDR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds500_JetDR05_10->SetBorderSize(0);
  leg_Zuds500_JetDR05_10->SetTextAlign(12);
  leg_Zuds500_JetDR05_10->SetTextSize(0.050);
  leg_Zuds500_JetDR05_10->SetTextFont(42);
  leg_Zuds500_JetDR05_10->SetMargin(0.15);
  leg_Zuds500_JetDR05_10->SetLineColor(1);
  leg_Zuds500_JetDR05_10->SetLineStyle(1);
  leg_Zuds500_JetDR05_10->SetLineWidth(1);
  leg_Zuds500_JetDR05_10->SetFillColor(0);
  //leg_Zuds500_JetDR05_10->SetFillStyle(1001);
  leg_Zuds500_JetDR05_10->SetFillStyle(0);
  leg_Zuds500_JetDR05_10->SetHeader("Z#rightarrow uds, 500 GeV, #approx 250 GeV Jets");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds500_JetDR05_10->AddEntry(h_Zuds500_JER_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds500_JetDR05_10->Draw();


  TH1F* h_Zuds3000_JER_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR05_3000");
  h_Zuds3000_JER_DR05->SetLineColor(kBlue);
  h_Zuds3000_JER_DR05->SetLineStyle(1);
  h_Zuds3000_JER_DR05->SetMinimum(1.750);
  h_Zuds3000_JER_DR05->SetMaximum(10.1);
  h_Zuds3000_JER_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds3000_JER_DR05->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds3000_JER_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR06_3000");
  h_Zuds3000_JER_DR06->SetLineColor(kRed);
  h_Zuds3000_JER_DR06->SetLineStyle(1);
  h_Zuds3000_JER_DR06->SetLineWidth(2);

  h_Zuds3000_JER_DR07->SetLineColor(kAzure+7);
  h_Zuds3000_JER_DR07->SetLineStyle(1);
  h_Zuds3000_JER_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR08_3000");
  h_Zuds3000_JER_DR08->SetLineColor(kViolet);
  h_Zuds3000_JER_DR08->SetLineStyle(1);
  h_Zuds3000_JER_DR08->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR09_3000");
  h_Zuds3000_JER_DR09->SetLineColor(kOrange);
  h_Zuds3000_JER_DR09->SetLineStyle(1);
  h_Zuds3000_JER_DR09->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR10_3000");
  h_Zuds3000_JER_DR10->SetLineColor(kYellow+1);
  h_Zuds3000_JER_DR10->SetLineStyle(1);
  h_Zuds3000_JER_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds3000_JetDR05_10 = setUpperCanvas("canvas_Zuds3000_JetDR05_10");
  
  TLegend* leg_Zuds3000_JetDR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds3000_JetDR05_10->SetBorderSize(0);
  leg_Zuds3000_JetDR05_10->SetTextAlign(12);
  leg_Zuds3000_JetDR05_10->SetTextSize(0.050);
  leg_Zuds3000_JetDR05_10->SetTextFont(42);
  leg_Zuds3000_JetDR05_10->SetMargin(0.15);
  leg_Zuds3000_JetDR05_10->SetLineColor(1);
  leg_Zuds3000_JetDR05_10->SetLineStyle(1);
  leg_Zuds3000_JetDR05_10->SetLineWidth(1);
  leg_Zuds3000_JetDR05_10->SetFillColor(0);
  //leg_Zuds3000_JetDR05_10->SetFillStyle(1001);
  leg_Zuds3000_JetDR05_10->SetFillStyle(0);
  leg_Zuds3000_JetDR05_10->SetHeader("Z#rightarrow uds, 3 TeV BG, #approx 1.5 TeV Jets");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds3000_JetDR05_10->AddEntry(h_Zuds3000_JER_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds3000_JetDR05_10->Draw();



  //NOW compare things WITH overlay

  TH1F* h_Zuds100_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_JER_wO_DR07->SetLineColor(kRed);
  h_Zuds100_JER_wO_DR07->SetLineColor(kCyan+1);
  h_Zuds100_JER_wO_DR07->SetLineStyle(1);
  h_Zuds100_JER_wO_DR07->SetLineWidth(2);
  h_Zuds100_JER_wO_DR07->SetMinimum(1.75);
  h_Zuds100_JER_wO_DR07->SetMaximum(15.);
  h_Zuds100_JER_wO_DR07->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_wO_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_200");
  h_Zuds200_JER_wO_DR07->SetLineColor(kYellow+1);
  h_Zuds200_JER_wO_DR07->SetLineStyle(1);
  h_Zuds200_JER_wO_DR07->SetLineWidth(2);
  h_Zuds200_JER_wO_DR07->SetMinimum(1.75);
  h_Zuds200_JER_wO_DR07->SetMaximum(15.0);
  TH1F* h_Zuds380_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_380");
  h_Zuds380_JER_wO_DR07->SetLineColor(kOrange);
  h_Zuds380_JER_wO_DR07->SetLineStyle(1);
  h_Zuds380_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_500");
  h_Zuds500_JER_wO_DR07->SetLineStyle(1);
  h_Zuds500_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_JER_wO_DR07->SetLineColor(kBlue);
  h_Zuds1500_JER_wO_DR07->SetLineStyle(1);
  h_Zuds1500_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_wO_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_JER_wO_DR07->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR07->SetLineWidth(2);
  h_Zuds3000_JER_wO_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_JER_wO_DR07_RMS90_CT_fancy = setUpperCanvas("resolutionGraphCanvas_JER_wO_DR07_RMS90_CT_fancy");
  //resolutionGraphCanvas_JER_wO_DR07_RMS90_CT_fancy->cd();
  //TLegend *leg_JER_wO_DR07_RMS90_CT_FullSummary = resolutionGraphCanvas_JER_wO_DR07_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_wO_DR07_RMS90_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetBorderSize(0);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetTextAlign(12);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetTextSize(0.050);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetTextFont(42);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetMargin(0.15);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetLineColor(1);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetLineStyle(1);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetLineWidth(1);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetFillColor(0);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetFillStyle(1001);
  leg_JER_wO_DR07_RMS90_CT_FullSummary->SetHeader("VLC7 Jets, with 3TeV BG");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds100_JER_wO_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds200_JER_wO_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds380_JER_wO_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds500_JER_wO_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds1500_JER_wO_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds3000_JER_wO_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_JER_wO_DR07_RMS90_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  h_Zuds100_JER_DR07->SetMaximum(15.0);



  TH1F* h_Zuds100_JER_wO_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR05_100");
  h_Zuds100_JER_wO_DR05->SetLineColor(kBlue);
  h_Zuds100_JER_wO_DR05->SetLineStyle(1);
  h_Zuds100_JER_wO_DR05->SetMinimum(1.750);
  h_Zuds100_JER_wO_DR05->SetMaximum(15.0);
  h_Zuds100_JER_wO_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_wO_DR05->GetXaxis()->SetTitle("|cos#theta|");

  TH1F* h_Zuds100_JER_wO_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR06_100");
  h_Zuds100_JER_wO_DR06->SetLineColor(kRed);
  h_Zuds100_JER_wO_DR06->SetLineStyle(1);
  h_Zuds100_JER_wO_DR06->SetLineWidth(2);

  h_Zuds100_JER_wO_DR07->SetLineColor(kAzure+7);
  h_Zuds100_JER_wO_DR07->SetLineStyle(1);
  h_Zuds100_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds100_JER_wO_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR08_100");
  h_Zuds100_JER_wO_DR08->SetLineColor(kViolet);
  h_Zuds100_JER_wO_DR08->SetLineStyle(1);
  h_Zuds100_JER_wO_DR08->SetLineWidth(2);
  TH1F* h_Zuds100_JER_wO_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR09_100");
  h_Zuds100_JER_wO_DR09->SetLineColor(kOrange);
  h_Zuds100_JER_wO_DR09->SetLineStyle(1);
  h_Zuds100_JER_wO_DR09->SetLineWidth(2);
  TH1F* h_Zuds100_JER_wO_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR10_100");
  h_Zuds100_JER_wO_DR10->SetLineColor(kYellow+1);
  h_Zuds100_JER_wO_DR10->SetLineStyle(1);
  h_Zuds100_JER_wO_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds100_JetwODR05_10 = setUpperCanvas("canvas_Zuds100_JetwODR05_10");
  
  TLegend* leg_Zuds100_JetwODR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds100_JetwODR05_10->SetBorderSize(0);
  leg_Zuds100_JetwODR05_10->SetTextAlign(12);
  leg_Zuds100_JetwODR05_10->SetTextSize(0.050);
  leg_Zuds100_JetwODR05_10->SetTextFont(42);
  leg_Zuds100_JetwODR05_10->SetMargin(0.15);
  leg_Zuds100_JetwODR05_10->SetLineColor(1);
  leg_Zuds100_JetwODR05_10->SetLineStyle(1);
  leg_Zuds100_JetwODR05_10->SetLineWidth(1);
  leg_Zuds100_JetwODR05_10->SetFillColor(0);
  //leg_Zuds100_JetwODR05_10->SetFillStyle(1001);
  leg_Zuds100_JetwODR05_10->SetFillStyle(0);
  leg_Zuds100_JetwODR05_10->SetHeader("Z#rightarrow uds, 91 GeV, #approx 45.5 GeV Jets");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds100_JetwODR05_10->AddEntry(h_Zuds100_JER_wO_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds100_JetwODR05_10->Draw();


  TH1F* h_Zuds500_JER_wO_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR05_500");
  h_Zuds500_JER_wO_DR05->SetLineColor(kBlue);
  h_Zuds500_JER_wO_DR05->SetLineStyle(1);
  h_Zuds500_JER_wO_DR05->SetMinimum(1.75);
  h_Zuds500_JER_wO_DR05->SetMaximum(15.0);
  h_Zuds500_JER_wO_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds500_JER_wO_DR05->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds500_JER_wO_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR06_500");
  h_Zuds500_JER_wO_DR06->SetLineColor(kRed);
  h_Zuds500_JER_wO_DR06->SetLineStyle(1);
  h_Zuds500_JER_wO_DR06->SetLineWidth(2);

  h_Zuds500_JER_wO_DR07->SetLineColor(kAzure+7);
  h_Zuds500_JER_wO_DR07->SetLineStyle(1);
  h_Zuds500_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_wO_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR08_500");
  h_Zuds500_JER_wO_DR08->SetLineColor(kViolet);
  h_Zuds500_JER_wO_DR08->SetLineStyle(1);
  h_Zuds500_JER_wO_DR08->SetLineWidth(2);
  TH1F* h_Zuds500_JER_wO_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR09_500");
  h_Zuds500_JER_wO_DR09->SetLineColor(kOrange);
  h_Zuds500_JER_wO_DR09->SetLineStyle(1);
  h_Zuds500_JER_wO_DR09->SetLineWidth(2);
  TH1F* h_Zuds500_JER_wO_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR10_500");
  h_Zuds500_JER_wO_DR10->SetLineColor(kYellow+1);
  h_Zuds500_JER_wO_DR10->SetLineStyle(1);
  h_Zuds500_JER_wO_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds500_JetwODR05_10 = setUpperCanvas("canvas_Zuds500_JetwODR05_10");
  
  TLegend* leg_Zuds500_JetwODR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds500_JetwODR05_10->SetBorderSize(0);
  leg_Zuds500_JetwODR05_10->SetTextAlign(12);
  leg_Zuds500_JetwODR05_10->SetTextSize(0.050);
  leg_Zuds500_JetwODR05_10->SetTextFont(42);
  leg_Zuds500_JetwODR05_10->SetMargin(0.15);
  leg_Zuds500_JetwODR05_10->SetLineColor(1);
  leg_Zuds500_JetwODR05_10->SetLineStyle(1);
  leg_Zuds500_JetwODR05_10->SetLineWidth(1);
  leg_Zuds500_JetwODR05_10->SetFillColor(0);
  //leg_Zuds500_JetwODR05_10->SetFillStyle(1001);
  leg_Zuds500_JetwODR05_10->SetFillStyle(0);
  leg_Zuds500_JetwODR05_10->SetHeader("Z#rightarrow uds, 500 GeV, #approx 250 GeV Jets");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds500_JetwODR05_10->AddEntry(h_Zuds500_JER_wO_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds500_JetwODR05_10->Draw();

  TH1F* h_Zuds3000_JER_wO_DR05 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR05_3000");
  h_Zuds3000_JER_wO_DR05->SetLineColor(kBlue);
  h_Zuds3000_JER_wO_DR05->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR05->SetMinimum(1.750);
  h_Zuds3000_JER_wO_DR05->SetMaximum(10.1);
  h_Zuds3000_JER_wO_DR05->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds3000_JER_wO_DR05->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds3000_JER_wO_DR06 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR06_3000");
  h_Zuds3000_JER_wO_DR06->SetLineColor(kRed);
  h_Zuds3000_JER_wO_DR06->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR06->SetLineWidth(2);

  h_Zuds3000_JER_wO_DR07->SetLineColor(kAzure+7);
  h_Zuds3000_JER_wO_DR07->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_wO_DR08 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR08_3000");
  h_Zuds3000_JER_wO_DR08->SetLineColor(kViolet);
  h_Zuds3000_JER_wO_DR08->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR08->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_wO_DR09 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR09_3000");
  h_Zuds3000_JER_wO_DR09->SetLineColor(kOrange);
  h_Zuds3000_JER_wO_DR09->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR09->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_wO_DR10 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR10_3000");
  h_Zuds3000_JER_wO_DR10->SetLineColor(kYellow+1);
  h_Zuds3000_JER_wO_DR10->SetLineStyle(1);
  h_Zuds3000_JER_wO_DR10->SetLineWidth(2);

  TCanvas *canvas_Zuds3000_JetwODR05_10 = setUpperCanvas("canvas_Zuds3000_JetwODR05_10");
  
  TLegend* leg_Zuds3000_JetwODR05_10=new TLegend(0.345,0.546,0.545,0.87);
  leg_Zuds3000_JetwODR05_10->SetBorderSize(0);
  leg_Zuds3000_JetwODR05_10->SetTextAlign(12);
  leg_Zuds3000_JetwODR05_10->SetTextSize(0.050);
  leg_Zuds3000_JetwODR05_10->SetTextFont(42);
  leg_Zuds3000_JetwODR05_10->SetMargin(0.15);
  leg_Zuds3000_JetwODR05_10->SetLineColor(1);
  leg_Zuds3000_JetwODR05_10->SetLineStyle(1);
  leg_Zuds3000_JetwODR05_10->SetLineWidth(1);
  leg_Zuds3000_JetwODR05_10->SetFillColor(0);
  //leg_Zuds3000_JetwODR05_10->SetFillStyle(1001);
  leg_Zuds3000_JetwODR05_10->SetFillStyle(0);
  leg_Zuds3000_JetwODR05_10->SetHeader("Z#rightarrow uds, 3 TeV BG, #approx 1.5 TeV Jets");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR05->DrawCopy("hist,e"),"#Delta R=0.5");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR06->DrawCopy("hist,e,same"),"#Delta R= 0.6");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR07->DrawCopy("hist,e,same"),"#Delta R= 0.7");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR08->DrawCopy("hist,e,same"),"#Delta R= 0.8");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR09->DrawCopy("hist,e,same"),"#Delta R= 0.9");
  leg_Zuds3000_JetwODR05_10->AddEntry(h_Zuds3000_JER_wO_DR10->DrawCopy("hist,e,same"),"#Delta R= 1.0");
  leg_Zuds3000_JetwODR05_10->Draw();



  //phi resolution
  TH1F* h_JER_DR07_wO_Zuds100 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_100");
  h_JER_DR07_wO_Zuds100->SetLineColor(kCyan+1);
  h_JER_DR07_wO_Zuds100->SetLineStyle(1);
  h_JER_DR07_wO_Zuds100->SetLineWidth(2);
  h_JER_DR07_wO_Zuds100->SetMinimum(1.75);
  h_JER_DR07_wO_Zuds100->SetMaximum(15.0);
  h_JER_DR07_wO_Zuds100->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_DR07_wO_Zuds100->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_JER_DR07_wO_Zuds200 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_200");
  h_JER_DR07_wO_Zuds200->SetLineColor(kYellow+1);
  h_JER_DR07_wO_Zuds200->SetLineStyle(1);
  h_JER_DR07_wO_Zuds200->SetLineWidth(2);
  h_JER_DR07_wO_Zuds200->SetMinimum(1.75);
  h_JER_DR07_wO_Zuds200->SetMaximum(15.0);
  TH1F* h_JER_DR07_wO_Zuds380 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_380");
  h_JER_DR07_wO_Zuds380->SetLineColor(kOrange);
  h_JER_DR07_wO_Zuds380->SetLineStyle(1);
  h_JER_DR07_wO_Zuds380->SetLineWidth(2);
  TH1F* h_JER_DR07_wO_Zuds500 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_500");
  h_JER_DR07_wO_Zuds500->SetLineStyle(1);
  h_JER_DR07_wO_Zuds500->SetLineWidth(2);
  h_JER_DR07_wO_Zuds500->SetMinimum(1.75);
  h_JER_DR07_wO_Zuds500->SetMaximum(10.1);
  h_JER_DR07_wO_Zuds500->SetLineColor(kGreen-2);
  TH1F* h_JER_DR07_wO_Zuds1500 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_1500");
  h_JER_DR07_wO_Zuds1500->SetLineStyle(1);
  h_JER_DR07_wO_Zuds1500->SetLineWidth(2);
  h_JER_DR07_wO_Zuds1500->SetLineColor(kBlue);
  TH1F* h_JER_DR07_wO_Zuds3000 = (TH1F*)file_histos_normalRange_jets->Get("JER_RMS90VsCosTheta_CT_wO_DR07_3000");
  h_JER_DR07_wO_Zuds3000->SetLineStyle(1);
  h_JER_DR07_wO_Zuds3000->SetLineWidth(2);
  h_JER_DR07_wO_Zuds3000->SetLineColor(kRed);
  /*
  TCanvas *canvas_JER_DR07_Zuds_wO = setUpperCanvas("canvas_Zuds_JER_DR07_wO_Energy");
  
  TLegend* leg_JER_DR07_Zuds_wO=new TLegend(0.345,0.546,0.545,0.87);
  leg_JER_DR07_Zuds_wO->SetBorderSize(0);
  leg_JER_DR07_Zuds_wO->SetTextAlign(12);
  leg_JER_DR07_Zuds_wO->SetTextSize(0.050);
  leg_JER_DR07_Zuds_wO->SetTextFont(42);
  leg_JER_DR07_Zuds_wO->SetMargin(0.15);
  leg_JER_DR07_Zuds_wO->SetLineColor(1);
  leg_JER_DR07_Zuds_wO->SetLineStyle(1);
  leg_JER_DR07_Zuds_wO->SetLineWidth(1);
  leg_JER_DR07_Zuds_wO->SetFillColor(0);
  //leg_JER_DR07_Zuds_wO->SetFillStyle(1001);
  leg_JER_DR07_Zuds_wO->SetFillStyle(0);
  leg_JER_DR07_Zuds_wO->SetHeader("VLC, #Delta R=0.7, #gamma=#beta=1, 3 TeV BKG");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds100->DrawCopy("hist,e"),"#approx 50 GeV Jets");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds200->DrawCopy("hist,e,same"),"#approx 100 GeV Jets");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds380->DrawCopy("hist,e,same"),"#approx 190 GeV Jets");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds500->DrawCopy("hist,e,same"),"#approx 250 GeV Jets");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds1500->DrawCopy("hist,e,same"),"#approx 750 GeV Jets");
  leg_JER_DR07_Zuds_wO->AddEntry(h_JER_DR07_wO_Zuds3000->DrawCopy("hist,e,same"),"#approx 1500 GeV Jets");
  leg_JER_DR07_Zuds_wO->Draw();
  */
  //phi resolution
  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetLineColor(kCyan+1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetLineWidth(2);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetMinimum(1.75);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetMaximum(10.1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->GetYaxis()->SetTitle("A.U.");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) / E_{tot}^{true}");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetMaximum(2995);

  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200 = (TH1F*)file_histos_normalRange_jets->Get("h_200_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->SetLineColor(kYellow+1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->SetLineWidth(2);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->SetMinimum(1.75);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->SetMaximum(10.1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->Scale(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->Integral()/h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->Integral());
  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380 = (TH1F*)file_histos_normalRange_jets->Get("h_380_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->SetLineColor(kOrange);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->SetLineWidth(2);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->Scale(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->Integral()/h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->Integral());
  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->SetLineWidth(2);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->SetMinimum(1.75);
  //h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->SetMaximum(10.1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->SetLineColor(kGreen-2);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->Scale(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->Integral()/h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->Integral());
  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500 = (TH1F*)file_histos_normalRange_jets->Get("h_1500_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->SetLineWidth(2);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->SetLineColor(kBlue);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->Scale(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->Integral()/h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->Integral());
  TH1F* h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_genjet07_E1_E2_over_E_True");
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->SetLineStyle(1);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->SetLineWidth(2);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->SetLineColor(kRed);
  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->Scale(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->Integral()/h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->Integral());

  TCanvas *canvas_genjets_E1_E2_over_E_True_DR07_Zuds_wO = setUpperCanvas("canvas_Zuds_genjets_E1_E2_over_E_True_DR07_wO_Energy");
  
  TLegend* leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO=new TLegend(0.345,0.546,0.545,0.87);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetBorderSize(0);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextAlign(12);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextSize(0.050);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextFont(42);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetMargin(0.15);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineColor(1);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineStyle(1);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineWidth(1);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillColor(0);
  //leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillStyle(1001);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillStyle(0);
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetHeader("VLC7 GenJets");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->DrawCopy("hist,e,same"),"#approx 100 GeV");
  //leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->DrawCopy("hist,e,same"),"#approx 250 GeV");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->DrawCopy("hist,e,same"),"#approx 750 GeV");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->DrawCopy("hist,e,same"),"#approx 1500 GeV");
  leg_genjets_E1_E2_over_E_True_DR07_Zuds_wO->Draw();

  l->DrawLatex(x,y,label.c_str());


  h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->SetMaximum(11995);
  TCanvas *canvas_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO = setUpperCanvas("canvas_Zuds_zoomedout_genjets_E1_E2_over_E_True_DR07_wO_Energy");
  
  TLegend* leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO=new TLegend(0.345,0.546,0.545,0.87);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetBorderSize(0);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextAlign(12);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextSize(0.050);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetTextFont(42);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetMargin(0.15);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineColor(1);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineStyle(1);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetLineWidth(1);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillColor(0);
  //leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillStyle(1001);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetFillStyle(0);
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->SetHeader("VLC7 GenJets");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds100->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds200->DrawCopy("hist,e,same"),"#approx 100 GeV");
  //leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds380->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds500->DrawCopy("hist,e,same"),"#approx 250 GeV");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds1500->DrawCopy("hist,e,same"),"#approx 750 GeV");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->AddEntry(h_genjets_E1_E2_over_E_True_DR07_wO_Zuds3000->DrawCopy("hist,e,same"),"#approx 1500 GeV");
  leg_zoomedout_genjets_E1_E2_over_E_True_DR07_Zuds_wO->Draw();

  l->DrawLatex(x,y,label.c_str());

  //phi resolution
  TH1F* h_Jet_reco_over_gen_DR07_Zuds100 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_100");
  h_Jet_reco_over_gen_DR07_Zuds100->SetLineColor(kCyan+1);
  h_Jet_reco_over_gen_DR07_Zuds100->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds100->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_Zuds100->SetMinimum(0.9);
  h_Jet_reco_over_gen_DR07_Zuds100->SetMaximum(1.1);
  h_Jet_reco_over_gen_DR07_Zuds100->GetYaxis()->SetTitle("Mean(E_{j}^{R}/E_{j}^{G})");
  h_Jet_reco_over_gen_DR07_Zuds100->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Jet_reco_over_gen_DR07_Zuds200 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_200");
  h_Jet_reco_over_gen_DR07_Zuds200->SetLineColor(kYellow+1);
  h_Jet_reco_over_gen_DR07_Zuds200->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds200->SetLineWidth(2);
  //h_Jet_reco_over_gen_DR07_Zuds200->SetMinimum(1.75);
  //h_Jet_reco_over_gen_DR07_Zuds200->SetMaximum(10.1);
  TH1F* h_Jet_reco_over_gen_DR07_Zuds380 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_380");
  h_Jet_reco_over_gen_DR07_Zuds380->SetLineColor(kOrange);
  h_Jet_reco_over_gen_DR07_Zuds380->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds380->SetLineWidth(2);
  TH1F* h_Jet_reco_over_gen_DR07_Zuds500 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_500");
  h_Jet_reco_over_gen_DR07_Zuds500->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds500->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_Zuds500->SetMinimum(0.9);
  h_Jet_reco_over_gen_DR07_Zuds500->SetMaximum(1.2);
  h_Jet_reco_over_gen_DR07_Zuds500->GetYaxis()->SetTitle("Mean(E^{R}_{j}/E^{G}_{j})");
  TH1F* h_Jet_reco_over_gen_DR07_Zuds1500 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_1500");
  h_Jet_reco_over_gen_DR07_Zuds1500->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds1500->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_Zuds1500->SetLineColor(kBlue);
  TH1F* h_Jet_reco_over_gen_DR07_Zuds3000 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_DR07_3000");
  h_Jet_reco_over_gen_DR07_Zuds3000->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_Zuds3000->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_Zuds3000->SetLineColor(kRed);


  TCanvas *canvas_Jet_reco_over_gen_DR07_Zuds = setUpperCanvas("canvas_Zuds_Jet_reco_over_gen_DR07_Energy");
  
  h_Jet_reco_over_gen_DR07_Zuds100->GetYaxis()->SetTitleOffset(1.35);
  h_Jet_reco_over_gen_DR07_Zuds100->GetYaxis()->SetLabelSize(0.05);

  TLegend* leg_Jet_reco_over_gen_DR07_Zuds=new TLegend(0.345,0.57,0.545,0.90);
  leg_Jet_reco_over_gen_DR07_Zuds->SetBorderSize(0);
  leg_Jet_reco_over_gen_DR07_Zuds->SetTextAlign(12);
  leg_Jet_reco_over_gen_DR07_Zuds->SetTextSize(0.050);
  leg_Jet_reco_over_gen_DR07_Zuds->SetTextFont(42);
  leg_Jet_reco_over_gen_DR07_Zuds->SetMargin(0.15);
  leg_Jet_reco_over_gen_DR07_Zuds->SetLineColor(1);
  leg_Jet_reco_over_gen_DR07_Zuds->SetLineStyle(1);
  leg_Jet_reco_over_gen_DR07_Zuds->SetLineWidth(1);
  leg_Jet_reco_over_gen_DR07_Zuds->SetFillColor(0);
  //leg_Jet_reco_over_gen_DR07_Zuds->SetFillStyle(1001);
  leg_Jet_reco_over_gen_DR07_Zuds->SetFillStyle(0);
  leg_Jet_reco_over_gen_DR07_Zuds->SetHeader("VLC7 Jets");
  leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds100->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds200->DrawCopy("hist,e,same"),"#approx 100 GeV");
  //leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds380->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds500->DrawCopy("hist,e,same"),"#approx 250 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds1500->DrawCopy("hist,e,same"),"#approx 750 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds->AddEntry(h_Jet_reco_over_gen_DR07_Zuds3000->DrawCopy("hist,e,same"),"#approx 1500 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds->Draw();

  l->DrawLatex(x,y,label.c_str());

  //phi resolution
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds100 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_100");
  h_Jet_reco_over_gen_DR07_wO_Zuds100->SetLineColor(kCyan+1);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->SetMinimum(0.9);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->SetMaximum(1.1);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->GetYaxis()->SetTitle("Mean(E_{j}^{R}/E_{j}^{G})");
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds200 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_200");
  h_Jet_reco_over_gen_DR07_wO_Zuds200->SetLineColor(kYellow+1);
  h_Jet_reco_over_gen_DR07_wO_Zuds200->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds200->SetLineWidth(2);
  //h_Jet_reco_over_gen_DR07_wO_Zuds200->SetMinimum(1.75);
  //h_Jet_reco_over_gen_DR07_wO_Zuds200->SetMaximum(10.1);
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds380 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_380");
  h_Jet_reco_over_gen_DR07_wO_Zuds380->SetLineColor(kOrange);
  h_Jet_reco_over_gen_DR07_wO_Zuds380->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds380->SetLineWidth(2);
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds500 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_500");
  h_Jet_reco_over_gen_DR07_wO_Zuds500->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds500->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_wO_Zuds500->SetMinimum(0.9);
  h_Jet_reco_over_gen_DR07_wO_Zuds500->SetMaximum(1.2);
  h_Jet_reco_over_gen_DR07_wO_Zuds500->GetYaxis()->SetTitle("Mean(E_{j}^{R}/E_{j}^{G})");
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds1500 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_1500");
  h_Jet_reco_over_gen_DR07_wO_Zuds1500->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds1500->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_wO_Zuds1500->SetLineColor(kBlue);
  TH1F* h_Jet_reco_over_gen_DR07_wO_Zuds3000 = (TH1F*)file_histos_normalRange_jets->Get("JER_MeanVsCosTheta_CT_wO_DR07_3000");
  h_Jet_reco_over_gen_DR07_wO_Zuds3000->SetLineStyle(1);
  h_Jet_reco_over_gen_DR07_wO_Zuds3000->SetLineWidth(2);
  h_Jet_reco_over_gen_DR07_wO_Zuds3000->SetLineColor(kRed);


  TCanvas *canvas_Jet_reco_over_gen_DR07_Zuds_wO = setUpperCanvas("canvas_Zuds_Jet_reco_over_gen_DR07_wO_Energy");
  
  h_Jet_reco_over_gen_DR07_wO_Zuds100->GetYaxis()->SetTitleOffset(1.35);
  h_Jet_reco_over_gen_DR07_wO_Zuds100->GetYaxis()->SetLabelSize(0.05);

  TLegend* leg_Jet_reco_over_gen_DR07_Zuds_wO=new TLegend(0.345,0.54,0.545,0.87);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetBorderSize(0);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetTextAlign(12);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetTextSize(0.050);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetTextFont(42);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetMargin(0.15);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetLineColor(1);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetLineStyle(1);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetLineWidth(1);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetFillColor(0);
  //leg_Jet_reco_over_gen_DR07_Zuds_wO->SetFillStyle(1001);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetFillStyle(0);
  leg_Jet_reco_over_gen_DR07_Zuds_wO->SetHeader("VLC7 Jets, with 3TeV BG");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds100->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds200->DrawCopy("hist,e,same"),"#approx 100 GeV");
  //leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds380->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds500->DrawCopy("hist,e,same"),"#approx 250 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds1500->DrawCopy("hist,e,same"),"#approx 750 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->AddEntry(h_Jet_reco_over_gen_DR07_wO_Zuds3000->DrawCopy("hist,e,same"),"#approx 1500 GeV");
  leg_Jet_reco_over_gen_DR07_Zuds_wO->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_100");
 //h_Zuds100_JER_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_JER_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_JER_FitCB_DR07->SetMinimum(1.75);
  h_Zuds100_JER_FitCB_DR07->SetMaximum(10.1);
  h_Zuds100_JER_FitCB_DR07->GetYaxis()->SetTitle("#sigma(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_200");
 //h_Zuds100_JER_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_JER_FitCB_DR07->SetLineWidth(2);
  h_Zuds200_JER_FitCB_DR07->SetMinimum(1.75);
  h_Zuds200_JER_FitCB_DR07->SetMaximum(10.1);
  TH1F* h_Zuds380_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_380");
 //h_Zuds380_JER_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_JER_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_500");
  h_Zuds500_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_JER_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_1500");
  h_Zuds1500_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_JER_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_DR07_3000");
  h_Zuds3000_JER_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_JER_FitCB_DR07->SetLineWidth(2);

  TCanvas *resolutionGraphCanvas_JER_FitCB_DR07_SigmaCB_CT_fancy = setUpperCanvas("resolutionGraphCanvas_JER_FitCB_DR07_SigmaCB_CT_fancy");
  //resolutionGraphCanvas_JER_FitCB_DR07_SigmaCB_CT_fancy->cd();
  //TLegend *leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary = resolutionGraphCanvas_JER_FitCB_DR07_SigmaCB_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetBorderSize(0);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetTextAlign(12);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetTextSize(0.050);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetTextFont(42);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetMargin(0.15);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetLineColor(1);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetLineStyle(1);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetLineWidth(1);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetFillColor(0);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetFillStyle(1001);
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds100_JER_FitCB_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds200_JER_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds380_JER_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds500_JER_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds1500_JER_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds3000_JER_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_JER_wO_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_JER_wO_FitCB_DR07->SetMinimum(1.75);
  h_Zuds100_JER_wO_FitCB_DR07->SetMaximum(15.0);
  h_Zuds100_JER_wO_FitCB_DR07->GetYaxis()->SetTitle("#sigma(E_{j}^{R}/E_{j}^{G})[%]");
  h_Zuds100_JER_wO_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_200");
 //h_Zuds100_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_JER_wO_FitCB_DR07->SetLineWidth(2);
  h_Zuds200_JER_wO_FitCB_DR07->SetMinimum(1.75);
  h_Zuds200_JER_wO_FitCB_DR07->SetMaximum(15.0);
  TH1F* h_Zuds380_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_380");
 //h_Zuds380_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_JER_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_500");
  h_Zuds500_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_JER_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_JER_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_JER_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets->Get("JER_SigmaCBVsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_JER_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_JER_wO_FitCB_DR07->SetLineWidth(2);

  TCanvas *resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy = setUpperCanvas("resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy");
  //resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->cd();
  //TLegend *leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary = resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetBorderSize(0);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetTextAlign(12);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetTextSize(0.050);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetTextFont(42);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetMargin(0.15);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetLineColor(1);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetLineStyle(1);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetLineWidth(1);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetFillColor(0);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetFillStyle(1001);
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->SetHeader("VLC7 Jets, with 3TeV BG");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds100_JER_wO_FitCB_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds200_JER_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds380_JER_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds500_JER_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds1500_JER_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds3000_JER_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_JER_wO_FitCB_DR07_SigmaCB_CT_wO_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dphi_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_dphi_0_to_065->SetLineStyle(1);
  h_Zuds500_dphi_0_to_065->SetLineWidth(2);
  h_Zuds500_dphi_0_to_065->SetLineColor(kBlack);
  h_Zuds500_dphi_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_dphi_0_to_065->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");
 
  TH1F* h_Zuds500_dphi_065_to_080 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dphi_reco_gen_cosTheta_0_65_to_0_80");
  h_Zuds500_dphi_065_to_080->SetLineStyle(1);
  h_Zuds500_dphi_065_to_080->SetLineWidth(2);
  h_Zuds500_dphi_065_to_080->SetLineColor(kRed);
  h_Zuds500_dphi_065_to_080->Scale(h_Zuds500_dphi_0_to_065->Integral()/h_Zuds500_dphi_065_to_080->Integral());


  TH1F* h_Zuds500_dphi_080_to_0925 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dphi_reco_gen_cosTheta_0_80_to_0_925");
  h_Zuds500_dphi_080_to_0925->SetLineStyle(1);
  h_Zuds500_dphi_080_to_0925->SetLineWidth(2);
  h_Zuds500_dphi_080_to_0925->SetLineColor(kBlue);
  h_Zuds500_dphi_080_to_0925->Scale(h_Zuds500_dphi_0_to_065->Integral()/h_Zuds500_dphi_080_to_0925->Integral());

  TH1F* h_Zuds500_dphi_0925_to_0975 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dphi_reco_gen_cosTheta_0_925_to_0_975");
  h_Zuds500_dphi_0925_to_0975->SetLineStyle(1);
  h_Zuds500_dphi_0925_to_0975->SetLineWidth(2);
  h_Zuds500_dphi_0925_to_0975->SetLineColor(kGreen-2);
  h_Zuds500_dphi_0925_to_0975->Scale(h_Zuds500_dphi_0_to_065->Integral()/h_Zuds500_dphi_0925_to_0975->Integral());



  TCanvas *resolutionGraphCanvas_dphi_regions_Zuds500 = setUpperCanvas("resolutionGraphCanvas_dphi_regions_Zuds500");
  //resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->cd();
  //TLegend *leg_dphi_regions_Zuds500 = resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphi_regions_Zuds500 = new TLegend(0.20,0.560,0.50,0.87);
  leg_dphi_regions_Zuds500->SetBorderSize(0);
  leg_dphi_regions_Zuds500->SetTextAlign(12);
  leg_dphi_regions_Zuds500->SetTextSize(0.050);
  leg_dphi_regions_Zuds500->SetTextFont(42);
  leg_dphi_regions_Zuds500->SetMargin(0.15);
  leg_dphi_regions_Zuds500->SetLineColor(1);
  leg_dphi_regions_Zuds500->SetLineStyle(1);
  leg_dphi_regions_Zuds500->SetLineWidth(1);
  leg_dphi_regions_Zuds500->SetFillColor(0);
  leg_dphi_regions_Zuds500->SetFillStyle(1001);
  leg_dphi_regions_Zuds500->SetHeader("VLC7 Jets,#approx 250GeV");
  leg_dphi_regions_Zuds500->AddEntry(h_Zuds500_dphi_0_to_065->DrawCopy("h,e"),"|cos#theta|<0.65");
  leg_dphi_regions_Zuds500->AddEntry(h_Zuds500_dphi_065_to_080->DrawCopy("h,e,same"),"0.65<|cos#theta|<0.80");
  leg_dphi_regions_Zuds500->AddEntry(h_Zuds500_dphi_080_to_0925->DrawCopy("h,e,same"),"0.80<|cos#theta|<0.925");
  leg_dphi_regions_Zuds500->AddEntry(h_Zuds500_dphi_0925_to_0975->DrawCopy("h,e,same"),"0.925<|cos#theta|<0.975");
  leg_dphi_regions_Zuds500->Draw();

  l->DrawLatex(x,y,label.c_str());



  TH1F* h_Zuds500_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dtheta_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds500_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds500_dtheta_0_to_065->SetLineColor(kBlack);
  h_Zuds500_dtheta_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_dtheta_0_to_065->GetXaxis()->SetTitle("#Delta#theta(j_{R},j_{G}) [#circ]");
 
  TH1F* h_Zuds500_dtheta_065_to_080 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dtheta_reco_gen_cosTheta_0_65_to_0_80");
  h_Zuds500_dtheta_065_to_080->SetLineStyle(1);
  h_Zuds500_dtheta_065_to_080->SetLineWidth(2);
  h_Zuds500_dtheta_065_to_080->SetLineColor(kRed);
  h_Zuds500_dtheta_065_to_080->Scale(h_Zuds500_dtheta_0_to_065->Integral()/h_Zuds500_dtheta_065_to_080->Integral());


  TH1F* h_Zuds500_dtheta_080_to_0925 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dtheta_reco_gen_cosTheta_0_80_to_0_925");
  h_Zuds500_dtheta_080_to_0925->SetLineStyle(1);
  h_Zuds500_dtheta_080_to_0925->SetLineWidth(2);
  h_Zuds500_dtheta_080_to_0925->SetLineColor(kBlue);
  h_Zuds500_dtheta_080_to_0925->Scale(h_Zuds500_dtheta_0_to_065->Integral()/h_Zuds500_dtheta_080_to_0925->Integral());

  TH1F* h_Zuds500_dtheta_0925_to_0975 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_DR07_dtheta_reco_gen_cosTheta_0_925_to_0_975");
  h_Zuds500_dtheta_0925_to_0975->SetLineStyle(1);
  h_Zuds500_dtheta_0925_to_0975->SetLineWidth(2);
  h_Zuds500_dtheta_0925_to_0975->SetLineColor(kGreen-2);
  h_Zuds500_dtheta_0925_to_0975->Scale(h_Zuds500_dtheta_0_to_065->Integral()/h_Zuds500_dtheta_0925_to_0975->Integral());

  TCanvas *resolutionGraphCanvas_dtheta_regions_Zuds500 = setUpperCanvas("resolutionGraphCanvas_dtheta_regions_Zuds500");
  //resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->cd();
  //TLegend *leg_dtheta_regions_Zuds500 = resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dtheta_regions_Zuds500 = new TLegend(0.20,0.560,0.50,0.87);
  leg_dtheta_regions_Zuds500->SetBorderSize(0);
  leg_dtheta_regions_Zuds500->SetTextAlign(12);
  leg_dtheta_regions_Zuds500->SetTextSize(0.050);
  leg_dtheta_regions_Zuds500->SetTextFont(42);
  leg_dtheta_regions_Zuds500->SetMargin(0.15);
  leg_dtheta_regions_Zuds500->SetLineColor(1);
  leg_dtheta_regions_Zuds500->SetLineStyle(1);
  leg_dtheta_regions_Zuds500->SetLineWidth(1);
  leg_dtheta_regions_Zuds500->SetFillColor(0);
  leg_dtheta_regions_Zuds500->SetFillStyle(1001);
  leg_dtheta_regions_Zuds500->SetHeader("VLC7 Jets,#approx 250GeV");
  leg_dtheta_regions_Zuds500->AddEntry(h_Zuds500_dtheta_0_to_065->DrawCopy("h,e"),"|cos#theta|<0.65");
  leg_dtheta_regions_Zuds500->AddEntry(h_Zuds500_dtheta_065_to_080->DrawCopy("h,e,same"),"0.65<|cos#theta|<0.80");
  leg_dtheta_regions_Zuds500->AddEntry(h_Zuds500_dtheta_080_to_0925->DrawCopy("h,e,same"),"0.80<|cos#theta|<0.925");
  leg_dtheta_regions_Zuds500->AddEntry(h_Zuds500_dtheta_0925_to_0975->DrawCopy("h,e,same"),"0.925<|cos#theta|<0.975");
  leg_dtheta_regions_Zuds500->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds500_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds500_wO_dphi_0_to_065->SetLineColor(kBlack);
  h_Zuds500_wO_dphi_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_dphi_0_to_065->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");

  TH1F* h_Zuds500_wO_dphi_065_to_080 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65_to_0_80");
  h_Zuds500_wO_dphi_065_to_080->SetLineStyle(1);
  h_Zuds500_wO_dphi_065_to_080->SetLineWidth(2);
  h_Zuds500_wO_dphi_065_to_080->SetLineColor(kRed);
  h_Zuds500_wO_dphi_065_to_080->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds500_wO_dphi_065_to_080->Integral());


  TH1F* h_Zuds500_wO_dphi_080_to_0925 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_80_to_0_925");
  h_Zuds500_wO_dphi_080_to_0925->SetLineStyle(1);
  h_Zuds500_wO_dphi_080_to_0925->SetLineWidth(2);
  h_Zuds500_wO_dphi_080_to_0925->SetLineColor(kBlue);
  h_Zuds500_wO_dphi_080_to_0925->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds500_wO_dphi_080_to_0925->Integral());

  TH1F* h_Zuds500_wO_dphi_0925_to_0975 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_925_to_0_975");
  h_Zuds500_wO_dphi_0925_to_0975->SetLineStyle(1);
  h_Zuds500_wO_dphi_0925_to_0975->SetLineWidth(2);
  h_Zuds500_wO_dphi_0925_to_0975->SetLineColor(kGreen-2);
  h_Zuds500_wO_dphi_0925_to_0975->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds500_wO_dphi_0925_to_0975->Integral());

  h_Zuds500_wO_dphi_0925_to_0975->GetXaxis()->SetRangeUser(-6.1,6.1);
  h_Zuds500_wO_dphi_0925_to_0975->SetMaximum(27500);
  h_Zuds500_wO_dphi_0925_to_0975->SetMinimum(0);
  h_Zuds500_wO_dphi_0925_to_0975->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_dphi_0925_to_0975->GetYaxis()->SetTitleOffset(1.55);
  h_Zuds500_wO_dphi_0925_to_0975->GetYaxis()->SetLabelSize(0.05);
  h_Zuds500_wO_dphi_0925_to_0975->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");

  TCanvas *resolutionGraphCanvas_wO_dphi_regions_Zuds500 = setUpperCanvas("resolutionGraphCanvas_wO_dphi_regions_Zuds500");
  //resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_wO_fancy->cd();
  //TLegend *leg_wO_dphi_regions_Zuds500 = resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_wO_dphi_regions_Zuds500 = new TLegend(0.20,0.590,0.50,0.90);
  leg_wO_dphi_regions_Zuds500->SetBorderSize(0);
  leg_wO_dphi_regions_Zuds500->SetTextAlign(12);
  leg_wO_dphi_regions_Zuds500->SetTextSize(0.050);
  leg_wO_dphi_regions_Zuds500->SetTextFont(42);
  leg_wO_dphi_regions_Zuds500->SetMargin(0.15);
  leg_wO_dphi_regions_Zuds500->SetLineColor(1);
  leg_wO_dphi_regions_Zuds500->SetLineStyle(1);
  leg_wO_dphi_regions_Zuds500->SetLineWidth(1);
  leg_wO_dphi_regions_Zuds500->SetFillColor(0);
  leg_wO_dphi_regions_Zuds500->SetFillStyle(1001);
  leg_wO_dphi_regions_Zuds500->SetHeader("VLC7 Jets,#approx 250GeV, with 3TeV BG");
  leg_wO_dphi_regions_Zuds500->AddEntry(h_Zuds500_wO_dphi_0925_to_0975->DrawCopy("h,e"),"0.925<|cos#theta|<0.975");
  leg_wO_dphi_regions_Zuds500->AddEntry(h_Zuds500_wO_dphi_080_to_0925->DrawCopy("h,e,same"),"0.80<|cos#theta|<0.925");
  leg_wO_dphi_regions_Zuds500->AddEntry(h_Zuds500_wO_dphi_065_to_080->DrawCopy("h,e,same"),"0.65<|cos#theta|<0.80");
  leg_wO_dphi_regions_Zuds500->AddEntry(h_Zuds500_wO_dphi_0_to_065->DrawCopy("h,e,same"),"|cos#theta|<0.65");
  leg_wO_dphi_regions_Zuds500->Draw();


  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds500_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds500_wO_dtheta_0_to_065->SetLineColor(kBlack);
  h_Zuds500_wO_dtheta_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_dtheta_0_to_065->GetXaxis()->SetTitle("#Delta#theta(j_{R},j_{G}) [#circ]");
 
  TH1F* h_Zuds500_wO_dtheta_065_to_080 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65_to_0_80");
  h_Zuds500_wO_dtheta_065_to_080->SetLineStyle(1);
  h_Zuds500_wO_dtheta_065_to_080->SetLineWidth(2);
  h_Zuds500_wO_dtheta_065_to_080->SetLineColor(kRed);
  h_Zuds500_wO_dtheta_065_to_080->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds500_wO_dtheta_065_to_080->Integral());


  TH1F* h_Zuds500_wO_dtheta_080_to_0925 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_80_to_0_925");
  h_Zuds500_wO_dtheta_080_to_0925->SetLineStyle(1);
  h_Zuds500_wO_dtheta_080_to_0925->SetLineWidth(2);
  h_Zuds500_wO_dtheta_080_to_0925->SetLineColor(kBlue);
  h_Zuds500_wO_dtheta_080_to_0925->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds500_wO_dtheta_080_to_0925->Integral());

  TH1F* h_Zuds500_wO_dtheta_0925_to_0975 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_925_to_0_975");
  h_Zuds500_wO_dtheta_0925_to_0975->SetLineStyle(1);
  h_Zuds500_wO_dtheta_0925_to_0975->SetLineWidth(2);
  h_Zuds500_wO_dtheta_0925_to_0975->SetLineColor(kGreen-2);
  h_Zuds500_wO_dtheta_0925_to_0975->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds500_wO_dtheta_0925_to_0975->Integral());

  TCanvas *resolutionGraphCanvas_wO_dtheta_regions_Zuds500 = setUpperCanvas("resolutionGraphCanvas_wO_dtheta_regions_Zuds500");
  //resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_wO_fancy->cd();
  //TLegend *leg_wO_dtheta_regions_Zuds500 = resolutionGraphCanvas_JER_wO_FitCB_DR07_SigmaCB_CT_wO_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_wO_dtheta_regions_Zuds500 = new TLegend(0.20,0.560,0.50,0.87);
  leg_wO_dtheta_regions_Zuds500->SetBorderSize(0);
  leg_wO_dtheta_regions_Zuds500->SetTextAlign(12);
  leg_wO_dtheta_regions_Zuds500->SetTextSize(0.050);
  leg_wO_dtheta_regions_Zuds500->SetTextFont(42);
  leg_wO_dtheta_regions_Zuds500->SetMargin(0.15);
  leg_wO_dtheta_regions_Zuds500->SetLineColor(1);
  leg_wO_dtheta_regions_Zuds500->SetLineStyle(1);
  leg_wO_dtheta_regions_Zuds500->SetLineWidth(1);
  leg_wO_dtheta_regions_Zuds500->SetFillColor(0);
  leg_wO_dtheta_regions_Zuds500->SetFillStyle(1001);
  leg_wO_dtheta_regions_Zuds500->SetHeader("VLC7 Jets,#approx 250GeV, with 3TeV BG");
  leg_wO_dtheta_regions_Zuds500->AddEntry(h_Zuds500_wO_dtheta_0_to_065->DrawCopy("h,e"),"|cos#theta|<0.65");
  leg_wO_dtheta_regions_Zuds500->AddEntry(h_Zuds500_wO_dtheta_065_to_080->DrawCopy("h,e,same"),"0.65<|cos#theta|<0.80");
  leg_wO_dtheta_regions_Zuds500->AddEntry(h_Zuds500_wO_dtheta_080_to_0925->DrawCopy("h,e,same"),"0.80<|cos#theta|<0.925");
  leg_wO_dtheta_regions_Zuds500->AddEntry(h_Zuds500_wO_dtheta_0925_to_0975->DrawCopy("h,e,same"),"0.925<|cos#theta|<0.975");
  leg_wO_dtheta_regions_Zuds500->Draw();

  l->DrawLatex(x,y,label.c_str());

  h_Zuds500_wO_dphi_0_to_065->SetLineColor(kGreen-2);
  TH1F* h_Zuds100_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_100_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
  h_Zuds100_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds100_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds100_wO_dphi_0_to_065->SetLineColor(kCyan+1);
  h_Zuds100_wO_dphi_0_to_065->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds100_wO_dphi_0_to_065->Integral());
  h_Zuds100_wO_dphi_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_wO_dphi_0_to_065->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");
  h_Zuds100_wO_dphi_0_to_065->SetMaximum (27000);
  TH1F* h_Zuds200_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_200_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
  h_Zuds200_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds200_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds200_wO_dphi_0_to_065->SetLineColor(kYellow+1);
  h_Zuds200_wO_dphi_0_to_065->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds200_wO_dphi_0_to_065->Integral());
  TH1F* h_Zuds380_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_380_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
  h_Zuds380_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds380_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds380_wO_dphi_0_to_065->SetLineColor(kOrange);
  h_Zuds380_wO_dphi_0_to_065->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds380_wO_dphi_0_to_065->Integral());
  TH1F* h_Zuds1500_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_1500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
  h_Zuds1500_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds1500_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds1500_wO_dphi_0_to_065->SetLineColor(kBlue);
  h_Zuds1500_wO_dphi_0_to_065->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds1500_wO_dphi_0_to_065->Integral());
  TH1F* h_Zuds3000_wO_dphi_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_3000_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
  h_Zuds3000_wO_dphi_0_to_065->SetLineStyle(1);
  h_Zuds3000_wO_dphi_0_to_065->SetLineWidth(2);
  h_Zuds3000_wO_dphi_0_to_065->SetLineColor(kRed);
  h_Zuds3000_wO_dphi_0_to_065->Scale(h_Zuds500_wO_dphi_0_to_065->Integral()/h_Zuds3000_wO_dphi_0_to_065->Integral());

  h_Zuds100_wO_dphi_0_to_065->GetXaxis()->SetRangeUser(-4.1,4.1);
  h_Zuds100_wO_dphi_0_to_065->GetYaxis()->SetRangeUser(0,24000);
  h_Zuds100_wO_dphi_0_to_065->GetYaxis()->SetLabelSize(0.05);

  TCanvas *resolutionGraphCanvas_wO_dphi_vs_jetE = setUpperCanvas("resolutionGraphCanvas_wO_dphi_vs_jetE");
  //resolutionGraphCanvas_wO_dphi_vs_jetE->cd();
  //TLegend *leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary = resolutionGraphCanvas_wO_dphi_vs_jetE->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_wO_dphi_vs_jetE = new TLegend(0.20,0.546,0.50,0.87);
  leg_wO_dphi_vs_jetE->SetBorderSize(0);
  leg_wO_dphi_vs_jetE->SetTextAlign(12);
  leg_wO_dphi_vs_jetE->SetTextSize(0.050);
  leg_wO_dphi_vs_jetE->SetTextFont(42);
  leg_wO_dphi_vs_jetE->SetMargin(0.15);
  leg_wO_dphi_vs_jetE->SetLineColor(1);
  leg_wO_dphi_vs_jetE->SetLineStyle(1);
  leg_wO_dphi_vs_jetE->SetLineWidth(1);
  leg_wO_dphi_vs_jetE->SetFillColor(0);
  leg_wO_dphi_vs_jetE->SetFillStyle(1001);
  leg_wO_dphi_vs_jetE->SetHeader("VLC7 Jets,|cos#theta_{jet}|<0.65, 3 TeV BG");
  leg_wO_dphi_vs_jetE->AddEntry(h_Zuds100_wO_dphi_0_to_065->DrawCopy("h,e"),"#approx 50 GeV");
  leg_wO_dphi_vs_jetE->AddEntry(h_Zuds200_wO_dphi_0_to_065->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_wO_dphi_vs_jetE->AddEntry(h_Zuds380_wO_dphi_0_to_065->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_wO_dphi_vs_jetE->AddEntry(h_Zuds500_wO_dphi_0_to_065->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_wO_dphi_vs_jetE->AddEntry(h_Zuds1500_wO_dphi_0_to_065->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_wO_dphi_vs_jetE->AddEntry(h_Zuds3000_wO_dphi_0_to_065->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_wO_dphi_vs_jetE->Draw();

  l->DrawLatex(x,y,label.c_str());

  h_Zuds500_wO_dtheta_0_to_065->SetLineColor(kGreen-2);
  TH1F* h_Zuds100_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_100_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
  h_Zuds100_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds100_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds100_wO_dtheta_0_to_065->SetLineColor(kCyan+1);
  h_Zuds100_wO_dtheta_0_to_065->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds100_wO_dtheta_0_to_065->Integral());
  h_Zuds100_wO_dtheta_0_to_065->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_wO_dtheta_0_to_065->GetXaxis()->SetTitle("#Delta#theta(j_{R},j_{G})");
  TH1F* h_Zuds200_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_200_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
  h_Zuds200_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds200_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds200_wO_dtheta_0_to_065->SetLineColor(kYellow+1);
  h_Zuds200_wO_dtheta_0_to_065->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds200_wO_dtheta_0_to_065->Integral());
  TH1F* h_Zuds380_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_380_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
  h_Zuds380_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds380_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds380_wO_dtheta_0_to_065->SetLineColor(kOrange);
  h_Zuds380_wO_dtheta_0_to_065->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds380_wO_dtheta_0_to_065->Integral());
  TH1F* h_Zuds1500_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_1500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
  h_Zuds1500_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds1500_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds1500_wO_dtheta_0_to_065->SetLineColor(kBlue);
  h_Zuds1500_wO_dtheta_0_to_065->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds1500_wO_dtheta_0_to_065->Integral());
  TH1F* h_Zuds3000_wO_dtheta_0_to_065 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("h_3000_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
  h_Zuds3000_wO_dtheta_0_to_065->SetLineStyle(1);
  h_Zuds3000_wO_dtheta_0_to_065->SetLineWidth(2);
  h_Zuds3000_wO_dtheta_0_to_065->SetLineColor(kRed);
  h_Zuds3000_wO_dtheta_0_to_065->Scale(h_Zuds500_wO_dtheta_0_to_065->Integral()/h_Zuds3000_wO_dtheta_0_to_065->Integral());

  TCanvas *resolutionGraphCanvas_wO_dtheta_vs_jetE = setUpperCanvas("resolutionGraphCanvas_wO_dtheta_vs_jetE");
  //resolutionGraphCanvas_wO_dtheta_vs_jetE->cd();
  //TLegend *leg_JER_FitCB_DR07_SigmaCB_CT_FullSummary = resolutionGraphCanvas_wO_dtheta_vs_jetE->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_wO_dtheta_vs_jetE = new TLegend(0.20,0.546,0.50,0.87);
  leg_wO_dtheta_vs_jetE->SetBorderSize(0);
  leg_wO_dtheta_vs_jetE->SetTextAlign(12);
  leg_wO_dtheta_vs_jetE->SetTextSize(0.050);
  leg_wO_dtheta_vs_jetE->SetTextFont(42);
  leg_wO_dtheta_vs_jetE->SetMargin(0.15);
  leg_wO_dtheta_vs_jetE->SetLineColor(1);
  leg_wO_dtheta_vs_jetE->SetLineStyle(1);
  leg_wO_dtheta_vs_jetE->SetLineWidth(1);
  leg_wO_dtheta_vs_jetE->SetFillColor(0);
  leg_wO_dtheta_vs_jetE->SetFillStyle(1001);
  leg_wO_dtheta_vs_jetE->SetHeader("VLC7 Jets, |cos#theta_{jet}|<0.65, with 3 TeV BG");
  leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds100_wO_dtheta_0_to_065->DrawCopy("h,e"),"#approx 50 GeV");
  leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds200_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds380_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds500_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds1500_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_wO_dtheta_vs_jetE->AddEntry(h_Zuds3000_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_wO_dtheta_vs_jetE->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_wO_dtheta_0_to_065_defBins = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_wO_DR07_dtheta_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_wO_dtheta_0_to_065_defBins->SetLineStyle(1);
  h_Zuds500_wO_dtheta_0_to_065_defBins->SetLineWidth(2);
  h_Zuds500_wO_dtheta_0_to_065_defBins->SetLineColor(kGreen-2);
  h_Zuds500_wO_dtheta_0_to_065_defBins->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_dtheta_0_to_065_defBins->GetXaxis()->SetTitle("#Delta#theta(j_{R},j_{G}) [#circ]");

  TCanvas *resolutionGraphCanvas_wO_dtheta_Zuds500_0_65 = setUpperCanvas("resolutionGraphCanvas_wO_dtheta_Zuds_0_65");
  //resolutionGraphCanvas_wO_dtheta_vs_jetE->cd();
  resolutionGraphCanvas_wO_dtheta_Zuds500_0_65->SetLogy();
  h_Zuds500_wO_dtheta_0_to_065_defBins->SetMinimum(0.5);
  h_Zuds500_wO_dtheta_0_to_065_defBins->SetMaximum(500000);
  h_Zuds500_wO_dtheta_0_to_065_defBins->DrawCopy("h,e");
  TLegend *leg_wO_dtheta_Zuds500_0_65 = new TLegend(0.60,0.61,0.90,0.87);
  leg_wO_dtheta_Zuds500_0_65->SetBorderSize(0);
  leg_wO_dtheta_Zuds500_0_65->SetTextAlign(12);
  leg_wO_dtheta_Zuds500_0_65->SetTextSize(0.050);
  leg_wO_dtheta_Zuds500_0_65->SetTextFont(42);
  leg_wO_dtheta_Zuds500_0_65->SetMargin(0.15);
  leg_wO_dtheta_Zuds500_0_65->SetLineColor(1);
  leg_wO_dtheta_Zuds500_0_65->SetLineStyle(1);
  leg_wO_dtheta_Zuds500_0_65->SetLineWidth(1);
  leg_wO_dtheta_Zuds500_0_65->SetFillColor(0);
  leg_wO_dtheta_Zuds500_0_65->SetFillStyle(1001);
  leg_wO_dtheta_Zuds500_0_65->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_wO_dtheta_Zuds500_0_65->AddEntry((TObject*)NULL,"#approx 250 GeV","h");
  leg_wO_dtheta_Zuds500_0_65->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.65","h");
  leg_wO_dtheta_Zuds500_0_65->AddEntry((TObject*)NULL,"3 TeV BG","h");
  //leg_wO_dtheta_Zuds500_0_65->SetHeader("VLC7,#approx250GeV,|cos#theta_{jet}|<0.65,3TeV BG");
  leg_wO_dtheta_Zuds500_0_65->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_wO_dphi_0_to_065_defBins = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_wO_DR07_dphi_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_wO_dphi_0_to_065_defBins->SetLineStyle(1);
  h_Zuds500_wO_dphi_0_to_065_defBins->SetLineWidth(2);
  h_Zuds500_wO_dphi_0_to_065_defBins->SetLineColor(kGreen-2);
  h_Zuds500_wO_dphi_0_to_065_defBins->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_dphi_0_to_065_defBins->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");

  TCanvas *resolutionGraphCanvas_wO_dphi_Zuds500_0_65 = setUpperCanvas("resolutionGraphCanvas_wO_dphi_Zuds_0_65");
  //resolutionGraphCanvas_wO_dphi_vs_jetE->cd();
  resolutionGraphCanvas_wO_dphi_Zuds500_0_65->SetLogy();
  h_Zuds500_wO_dphi_0_to_065_defBins->SetMinimum(0.5);
  h_Zuds500_wO_dphi_0_to_065_defBins->SetMaximum(500000);
  h_Zuds500_wO_dphi_0_to_065_defBins->DrawCopy("h,e");
  TLegend *leg_wO_dphi_Zuds500_0_65 = new TLegend(0.60,0.61,0.90,0.87);
  leg_wO_dphi_Zuds500_0_65->SetBorderSize(0);
  leg_wO_dphi_Zuds500_0_65->SetTextAlign(12);
  leg_wO_dphi_Zuds500_0_65->SetTextSize(0.050);
  leg_wO_dphi_Zuds500_0_65->SetTextFont(42);
  leg_wO_dphi_Zuds500_0_65->SetMargin(0.15);
  leg_wO_dphi_Zuds500_0_65->SetLineColor(1);
  leg_wO_dphi_Zuds500_0_65->SetLineStyle(1);
  leg_wO_dphi_Zuds500_0_65->SetLineWidth(1);
  leg_wO_dphi_Zuds500_0_65->SetFillColor(0);
  leg_wO_dphi_Zuds500_0_65->SetFillStyle(1001);
  leg_wO_dphi_Zuds500_0_65->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_wO_dphi_Zuds500_0_65->AddEntry((TObject*)NULL,"#approx 250 GeV","h");
  leg_wO_dphi_Zuds500_0_65->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.65","h");
  leg_wO_dphi_Zuds500_0_65->AddEntry((TObject*)NULL,"3 TeV BG","h");
  //leg_wO_dphi_Zuds500_0_65->SetHeader("VLC7,#approx250GeV,|cos#theta_{jet}|<0.65,3TeV BG");
  leg_wO_dphi_Zuds500_0_65->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds500_dphi_0_to_065_defBins = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_DR07_dphi_reco_gen_cosTheta_0_65");
 //h_Zuds500_JER_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds500_dphi_0_to_065_defBins->SetLineStyle(1);
  h_Zuds500_dphi_0_to_065_defBins->SetLineWidth(2);
  h_Zuds500_dphi_0_to_065_defBins->SetLineColor(kGreen-2);
  h_Zuds500_dphi_0_to_065_defBins->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_dphi_0_to_065_defBins->GetXaxis()->SetTitle("#Delta#phi(j_{R},j_{G}) [#circ]");
  TCanvas *resolutionGraphCanvas_dphi_Zuds500_0_65 = setUpperCanvas("resolutionGraphCanvas_dphi_Zuds_0_65");
  //resolutionGraphCanvas_dphi_vs_jetE->cd();
  h_Zuds500_dphi_0_to_065_defBins->SetMaximum(350000);
  h_Zuds500_dphi_0_to_065_defBins->SetMinimum(0.5);
  h_Zuds500_dphi_0_to_065_defBins->SetLineColor(kGreen-2);
  h_Zuds500_dphi_0_to_065_defBins->DrawCopy("h,e");
  TLegend *leg_dphi_Zuds500_0_65 = new TLegend(0.20,0.815,0.50,0.87);
  leg_dphi_Zuds500_0_65->SetBorderSize(0);
  leg_dphi_Zuds500_0_65->SetTextAlign(12);
  leg_dphi_Zuds500_0_65->SetTextSize(0.050);
  leg_dphi_Zuds500_0_65->SetTextFont(42);
  leg_dphi_Zuds500_0_65->SetMargin(0.15);
  leg_dphi_Zuds500_0_65->SetLineColor(1);
  leg_dphi_Zuds500_0_65->SetLineStyle(1);
  leg_dphi_Zuds500_0_65->SetLineWidth(1);
  leg_dphi_Zuds500_0_65->SetFillColor(0);
  leg_dphi_Zuds500_0_65->SetFillStyle(1001);
  leg_dphi_Zuds500_0_65->SetHeader("VLC7 Jets,#approx250GeV,|cos#theta_{jet}|<0.65");
  leg_dphi_Zuds500_0_65->Draw();

  l->DrawLatex(x,y,label.c_str());


 TH1F* h_Zuds100_wO_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds100_wO_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds100_wO_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds100_wO_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds100_wO_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_wO_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds100_wO_E_reco_over_E_gen_0_10->SetMaximum(700);
  h_Zuds100_wO_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_wO_0_10 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds100_0_10");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  resolutionGraphCanvas_deltaEjet_Zuds100_wO_0_10->SetLogy();
  h_Zuds100_wO_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds100_wO_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds100_wO_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds100_wO_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds100_wO_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_wO_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds100_wO_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds100_wO_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds100_wO_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds100_wO_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds100_wO_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds100_wO_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_wO_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds100_wO_0_10->AddEntry((TObject*)NULL,"#approx 50 GeV","h");
  leg_deltaEjet_Zuds100_wO_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds100_wO_0_10->AddEntry((TObject*)NULL,"3 TeV BG","h");
  //leg_deltaEjet_Zuds100_wO_0_10->SetHeader("VLC7 Jets,#approx250GeV,|cos#theta_{jet}|<0.10,3TeV BG");
  leg_deltaEjet_Zuds100_wO_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());

 TH1F* h_Zuds3000_wO_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->SetMaximum(700);
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_wO_0_10 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds3000_0_10");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  resolutionGraphCanvas_deltaEjet_Zuds3000_wO_0_10->SetLogy();
  h_Zuds3000_wO_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds3000_wO_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds3000_wO_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_wO_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_wO_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_wO_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds3000_wO_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_wO_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds3000_wO_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_wO_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_wO_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds3000_wO_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds3000_wO_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds3000_wO_0_10->AddEntry((TObject*)NULL,"#approx 1500 GeV","h");
  leg_deltaEjet_Zuds3000_wO_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds3000_wO_0_10->AddEntry((TObject*)NULL,"3 TeV BG","h");
  //leg_deltaEjet_Zuds3000_wO_0_10->SetHeader("VLC7 Jets,#approx250GeV,|cos#theta_{jet}|<0.10,3TeV BG");
  leg_deltaEjet_Zuds3000_wO_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_wO_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds500_wO_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds500_wO_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds500_wO_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds500_wO_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds500_wO_E_reco_over_E_gen_0_10->SetMaximum(3750);
  h_Zuds500_wO_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_wO_0_10 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds500_0_10");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  resolutionGraphCanvas_deltaEjet_Zuds500_wO_0_10->SetLogy();
  h_Zuds500_wO_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds500_wO_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds500_wO_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds500_wO_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds500_wO_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_wO_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds500_wO_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds500_wO_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds500_wO_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds500_wO_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds500_wO_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds500_wO_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_wO_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds500_wO_0_10->AddEntry((TObject*)NULL,"#approx 250 GeV","h");
  leg_deltaEjet_Zuds500_wO_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds500_wO_0_10->AddEntry((TObject*)NULL,"3 TeV BG","h");
  //leg_deltaEjet_Zuds500_wO_0_10->SetHeader("VLC7 Jets,#approx250GeV,|cos#theta_{jet}|<0.10,3TeV BG");
  leg_deltaEjet_Zuds500_wO_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds500_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds500_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds500_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds500_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds500_E_reco_over_E_gen_0_10->SetMaximum(3750);
  h_Zuds500_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_0_10 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds500_0_10");
  resolutionGraphCanvas_deltaEjet_Zuds500_0_10->SetLogy();
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds500_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds500_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds500_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds500_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds500_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds500_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds500_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds500_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds500_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds500_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds500_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds500_0_10->AddEntry((TObject*)NULL,"#approx 250 GeV","h");
  leg_deltaEjet_Zuds500_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds500_0_10->AddEntry((TObject*)NULL,"","h");
  leg_deltaEjet_Zuds500_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds100_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds100_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds100_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds100_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds100_E_reco_over_E_gen_0_10->SetMaximum(700);
  h_Zuds100_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_0_10 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds100_0_10");
  resolutionGraphCanvas_deltaEjet_Zuds100_0_10->SetLogy();
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds100_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds100_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds100_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds100_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds100_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds100_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds100_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds100_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds100_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds100_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds100_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds100_0_10->AddEntry((TObject*)NULL,"#approx 50 GeV","h");
  leg_deltaEjet_Zuds100_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds100_0_10->AddEntry((TObject*)NULL,"","h");
  leg_deltaEjet_Zuds100_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds3000_E_reco_over_E_gen_0_10 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_to_0_10");
  h_Zuds3000_E_reco_over_E_gen_0_10->SetLineStyle(1);
  h_Zuds3000_E_reco_over_E_gen_0_10->SetLineWidth(2);
  h_Zuds3000_E_reco_over_E_gen_0_10->SetLineColor(kGreen+2);
  h_Zuds3000_E_reco_over_E_gen_0_10->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_E_reco_over_E_gen_0_10->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");
  h_Zuds3000_E_reco_over_E_gen_0_10->SetMaximum(700);
  h_Zuds3000_E_reco_over_E_gen_0_10->SetMinimum(0.5);

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_0_10 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds3000_0_10");
  resolutionGraphCanvas_deltaEjet_Zuds3000_0_10->SetLogy();
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds3000_E_reco_over_E_gen_0_10->DrawCopy("h,e");
  TLegend *leg_deltaEjet_Zuds3000_0_10 = new TLegend(0.60,0.61,0.90,0.87);
  leg_deltaEjet_Zuds3000_0_10->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_0_10->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_0_10->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_0_10->SetTextFont(42);
  leg_deltaEjet_Zuds3000_0_10->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_0_10->SetLineColor(1);
  leg_deltaEjet_Zuds3000_0_10->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_0_10->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_0_10->SetFillColor(0);
  leg_deltaEjet_Zuds3000_0_10->SetFillStyle(1001);
  leg_deltaEjet_Zuds3000_0_10->AddEntry((TObject*)NULL,"VLC7 Jets","h");
  leg_deltaEjet_Zuds3000_0_10->AddEntry((TObject*)NULL,"#approx 1500 GeV","h");
  leg_deltaEjet_Zuds3000_0_10->AddEntry((TObject*)NULL,"|cos#theta_{jet}|<0.10","h");
  leg_deltaEjet_Zuds3000_0_10->AddEntry((TObject*)NULL,"","h");
  leg_deltaEjet_Zuds3000_0_10->Draw();

  l->DrawLatex(x,y,label.c_str());
  //now forward bu8ut not the acceptance 0_925_to_0_95

  TH1F* h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_wO_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds500_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds500_wO_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds500_wO_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->SetHeader("VLC07, wBG, 0.925<|cos#theta_{jet}|<0.95,Zuds500");
  leg_deltaEjet_Zuds500_wO_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds500_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds500_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds500_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_0_925_to_0_95->SetHeader("VLC07, no BG, 0.925<|cos#theta_{jet}|<0.95,Zuds500");
  leg_deltaEjet_Zuds500_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());


 TH1F* h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_wO_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds3000_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds3000_wO_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds3000_wO_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetFillStyle(30001);
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->SetHeader("VLC07, wBG, 0.925<|cos#theta_{jet}|<0.95,Zuds3000");
  leg_deltaEjet_Zuds3000_wO_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds3000_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds3000_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds3000_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetFillStyle(30001);
  leg_deltaEjet_Zuds3000_0_925_to_0_95->SetHeader("VLC07, no BG, 0.925<|cos#theta_{jet}|<0.95,Zuds3000");
  leg_deltaEjet_Zuds3000_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_wO_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds100_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds100_wO_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds100_wO_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->SetHeader("VLC07, wBG, 0.925<|cos#theta_{jet}|<0.95,Zuds100");
  leg_deltaEjet_Zuds100_wO_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_E_reco_over_E_gen_0_925_to_0_95 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_925_to_0_950");
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->SetLineStyle(1);
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->SetLineWidth(2);
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->SetLineColor(kGreen+2);
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_0_925_to_0_95 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds100_0_925_to_0_950");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds100_E_reco_over_E_gen_0_925_to_0_95->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds100_0_925_to_0_95 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetBorderSize(0);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetTextAlign(12);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetTextFont(42);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetMargin(0.15);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetLineColor(1);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetLineStyle(1);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetLineWidth(1);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetFillColor(0);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_0_925_to_0_95->SetHeader("VLC07, no BG, 0.925<|cos#theta_{jet}|<0.95,Zuds100");
  leg_deltaEjet_Zuds100_0_925_to_0_95->Draw();

  l->DrawLatex(x,y,label.c_str());

  //now the 0_975_to_1_00 histograms



  TH1F* h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_wO_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds500_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds500_wO_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds500_wO_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->SetHeader("VLC07, wBG, 0.975<|cos#theta_{jet}|<1.00,Zuds500");
  leg_deltaEjet_Zuds500_wO_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds500_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_500_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds500_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds500_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds500_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds500_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetFillStyle(1001);
  leg_deltaEjet_Zuds500_0_975_to_1_00->SetHeader("VLC07, no BG, 0.975<|cos#theta_{jet}|<1.00,Zuds500");
  leg_deltaEjet_Zuds500_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());


 TH1F* h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_wO_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds3000_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds3000_wO_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds3000_wO_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetFillStyle(30001);
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->SetHeader("VLC07, wBG, 0.975<|cos#theta_{jet}|<1.00,Zuds3000");
  leg_deltaEjet_Zuds3000_wO_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_3000_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds3000_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds3000_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds3000_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds3000_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetFillStyle(30001);
  leg_deltaEjet_Zuds3000_0_975_to_1_00->SetHeader("VLC07, no BG, 0.975<|cos#theta_{jet}|<1.00,Zuds3000");
  leg_deltaEjet_Zuds3000_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_wO_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_wO_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_wO_deltaEjet_Zuds100_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds100_wO_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");

  l->DrawLatex(x,y,label.c_str());

  TLegend *leg_deltaEjet_Zuds100_wO_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->SetHeader("VLC07, wBG, 0.975<|cos#theta_{jet}|<1.00,Zuds100");
  leg_deltaEjet_Zuds100_wO_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_E_reco_over_E_gen_0_975_to_1_00 = (TH1F*)file_histos_normalRange_jets->Get("h_100_CT_DR07_E_rel_totVis_RMS_d1_cosT_0_975_to_1_00");
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->SetLineStyle(1);
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->SetLineWidth(2);
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->SetLineColor(kGreen+2);
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->GetYaxis()->SetTitle("A.U.");
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->GetXaxis()->SetTitle("E_{recojet}/E_{genjet}");

  TCanvas *resolutionGraphCanvas_deltaEjet_Zuds100_0_975_to_1_00 = setUpperCanvas("resolutionGraphCanvas_deltaEjet_Zuds100_0_975_to_1_00");
  //resolutionGraphCanvas_deltaEjet_vs_jetE->cd();
  h_Zuds100_E_reco_over_E_gen_0_975_to_1_00->DrawCopy("h,e");



  TLegend *leg_deltaEjet_Zuds100_0_975_to_1_00 = new TLegend(0.20,0.815,0.50,0.87);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetBorderSize(0);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetTextAlign(12);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetTextSize(0.050);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetTextFont(42);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetMargin(0.15);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetLineColor(1);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetLineStyle(1);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetLineWidth(1);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetFillColor(0);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetFillStyle(1001);
  leg_deltaEjet_Zuds100_0_975_to_1_00->SetHeader("VLC07, no BG, 0.975<|cos#theta_{jet}|<1.00,Zuds100");
  leg_deltaEjet_Zuds100_0_975_to_1_00->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds100_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_100");
 //h_Zuds100_dphiRes_DR07->SetLineColor(kRed);
  h_Zuds100_dphiRes_DR07->SetLineStyle(1);
  h_Zuds100_dphiRes_DR07->SetLineWidth(2);
  h_Zuds100_dphiRes_DR07->SetMinimum(0);
  h_Zuds100_dphiRes_DR07->SetMaximum(1.75);
  h_Zuds100_dphiRes_DR07->GetYaxis()->SetTitle("RMS_{90}#Delta#phi(j_{R},j_{G}) [#circ]");
  h_Zuds100_dphiRes_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_200");
 //h_Zuds100_dphiRes_DR07->SetLineColor(kRed);
  h_Zuds200_dphiRes_DR07->SetLineStyle(1);
  h_Zuds200_dphiRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_380");
 //h_Zuds380_dphiRes_DR07->SetLineColor(kRed);
  h_Zuds380_dphiRes_DR07->SetLineStyle(1);
  h_Zuds380_dphiRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_500");
  h_Zuds500_dphiRes_DR07->SetLineStyle(1);
  h_Zuds500_dphiRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_1500");
  h_Zuds1500_dphiRes_DR07->SetLineStyle(1);
  h_Zuds1500_dphiRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dphiRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dphiRes_RMS90VsCosTheta_CT_DR07_3000");
  h_Zuds3000_dphiRes_DR07->SetLineStyle(1);
  h_Zuds3000_dphiRes_DR07->SetLineWidth(2);

  h_Zuds100_dphiRes_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dphiRes_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dphiRes_DR07->SetLineColor(kOrange);
  h_Zuds500_dphiRes_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dphiRes_DR07->SetLineColor(kBlue);
  h_Zuds3000_dphiRes_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dphiRes_DR07_RMS90_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dphiRes_DR07_RMS90_CT_fancy");
  //resolutionGraphCanvas_dphiRes_DR07_RMS90_CT_fancy->cd();
  //TLegend *leg_dphiRes_DR07_RMS90_CT_FullSummary = resolutionGraphCanvas_dphiRes_DR07_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphiRes_DR07_RMS90_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetBorderSize(0);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetTextAlign(12);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetTextSize(0.050);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetTextFont(42);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetMargin(0.15);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetLineColor(1);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetLineStyle(1);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetLineWidth(1);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetFillColor(0);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetFillStyle(1001);
  leg_dphiRes_DR07_RMS90_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds100_dphiRes_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds200_dphiRes_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds380_dphiRes_DR07->DrawCopy("h,e,same"),"#approx 190 GeV Jets");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds500_dphiRes_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds1500_dphiRes_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds3000_dphiRes_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dphiRes_DR07_RMS90_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds100_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_100");
 //h_Zuds100_dphiRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_dphiRes_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_dphiRes_FitCB_DR07->SetMinimum(0);
  h_Zuds100_dphiRes_FitCB_DR07->SetMaximum(1.75);
  h_Zuds100_dphiRes_FitCB_DR07->GetYaxis()->SetTitle("#sigma(#Delta#phi(j_{R},j_{G})) [#circ]");
  h_Zuds100_dphiRes_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_200");
 //h_Zuds100_dphiRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_dphiRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_380");
 //h_Zuds380_dphiRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_dphiRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_500");
  h_Zuds500_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_dphiRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_1500");
  h_Zuds1500_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_dphiRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dphiRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_DR07_3000");
  h_Zuds3000_dphiRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_dphiRes_FitCB_DR07->SetLineWidth(2);

  h_Zuds100_dphiRes_FitCB_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dphiRes_FitCB_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dphiRes_FitCB_DR07->SetLineColor(kOrange);
  h_Zuds500_dphiRes_FitCB_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dphiRes_FitCB_DR07->SetLineColor(kBlue);
  h_Zuds3000_dphiRes_FitCB_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dphiRes_DR07_SigmaCB_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dphiRes_DR07_SigmaCB_CT_fancy");
  //resolutionGraphCanvas_dphiRes_DR07_SigmaCB_CT_fancy->cd();
  //TLegend *leg_dphiRes_DR07_SigmaCB_CT_FullSummary = resolutionGraphCanvas_dphiRes_DR07_SigmaCB_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphiRes_DR07_SigmaCB_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetBorderSize(0);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetTextAlign(12);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetTextSize(0.050);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetTextFont(42);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetMargin(0.15);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetLineColor(1);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetLineStyle(1);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetLineWidth(1);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetFillColor(0);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetFillStyle(1001);
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds100_dphiRes_FitCB_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds200_dphiRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds380_dphiRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds500_dphiRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds1500_dphiRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds3000_dphiRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dphiRes_DR07_SigmaCB_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds100_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_100");
 //h_Zuds100_dthetaRes_DR07->SetLineColor(kRed);
  h_Zuds100_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds100_dthetaRes_DR07->SetLineWidth(2);
  h_Zuds100_dthetaRes_DR07->SetMinimum(0);
  h_Zuds100_dthetaRes_DR07->SetMaximum(1.75);
  h_Zuds100_dthetaRes_DR07->GetYaxis()->SetTitle("RMS_{90}#Delta#theta(j_{R},j_{G}) [#circ]");
  h_Zuds100_dthetaRes_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_200");
 //h_Zuds100_dthetaRes_DR07->SetLineColor(kRed);
  h_Zuds200_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds200_dthetaRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_380");
 //h_Zuds380_dthetaRes_DR07->SetLineColor(kRed);
  h_Zuds380_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds380_dthetaRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_500");
  h_Zuds500_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds500_dthetaRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_1500");
  h_Zuds1500_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds1500_dthetaRes_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dthetaRes_DR07 = (TH1F*)file_histos_normalRange_jets->Get("dthetaRes_RMS90VsCosTheta_CT_DR07_3000");
  h_Zuds3000_dthetaRes_DR07->SetLineStyle(1);
  h_Zuds3000_dthetaRes_DR07->SetLineWidth(2);


  h_Zuds100_dthetaRes_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dthetaRes_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dthetaRes_DR07->SetLineColor(kOrange);
  h_Zuds500_dthetaRes_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dthetaRes_DR07->SetLineColor(kBlue);
  h_Zuds3000_dthetaRes_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dthetaRes_DR07_RMS90_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dthetaRes_DR07_RMS90_CT_fancy");
  //resolutionGraphCanvas_dthetaRes_DR07_RMS90_CT_fancy->cd();
  //TLegend *leg_dthetaRes_DR07_RMS90_CT_FullSummary = resolutionGraphCanvas_dthetaRes_DR07_RMS90_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dthetaRes_DR07_RMS90_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetBorderSize(0);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetTextAlign(12);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetTextSize(0.050);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetTextFont(42);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetMargin(0.15);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetLineColor(1);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetLineStyle(1);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetLineWidth(1);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetFillColor(0);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetFillStyle(1001);
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds100_dthetaRes_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds200_dthetaRes_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds380_dthetaRes_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds500_dthetaRes_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds1500_dthetaRes_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->AddEntry(h_Zuds3000_dthetaRes_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dthetaRes_DR07_RMS90_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_100");
 //h_Zuds100_dthetaRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_dthetaRes_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_dthetaRes_FitCB_DR07->SetMinimum(0);
  h_Zuds100_dthetaRes_FitCB_DR07->SetMaximum(1.75);
  h_Zuds100_dthetaRes_FitCB_DR07->GetYaxis()->SetTitle("#sigma(#Delta#theta(j_{R},j_{G})) [#circ]");
  h_Zuds100_dthetaRes_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_200");
 //h_Zuds100_dthetaRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_dthetaRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_380");
 //h_Zuds380_dthetaRes_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_dthetaRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_500");
  h_Zuds500_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_dthetaRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_1500");
  h_Zuds1500_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_dthetaRes_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dthetaRes_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_DR07_3000");
  h_Zuds3000_dthetaRes_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_dthetaRes_FitCB_DR07->SetLineWidth(2);

  h_Zuds100_dthetaRes_FitCB_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dthetaRes_FitCB_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dthetaRes_FitCB_DR07->SetLineColor(kOrange);
  h_Zuds500_dthetaRes_FitCB_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dthetaRes_FitCB_DR07->SetLineColor(kBlue);
  h_Zuds3000_dthetaRes_FitCB_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dthetaRes_DR07_SigmaCB_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dthetaRes_DR07_SigmaCB_CT_fancy");
  //resolutionGraphCanvas_dthetaRes_DR07_SigmaCB_CT_fancy->cd();
  //TLegend *leg_dthetaRes_DR07_SigmaCB_CT_FullSummary = resolutionGraphCanvas_dthetaRes_DR07_SigmaCB_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dthetaRes_DR07_SigmaCB_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetBorderSize(0);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetTextAlign(12);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetTextSize(0.050);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetTextFont(42);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetMargin(0.15);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetLineColor(1);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetLineStyle(1);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetLineWidth(1);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetFillColor(0);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetFillStyle(1001);
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->SetHeader("VLC7 Jets");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds100_dthetaRes_FitCB_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds200_dthetaRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds380_dthetaRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds500_dthetaRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds1500_dthetaRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->AddEntry(h_Zuds3000_dthetaRes_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dthetaRes_DR07_SigmaCB_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  //NOW position resolution with overlay
  TH1F* h_Zuds100_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_dphiRes_wO_DR07->SetLineColor(kRed);
  h_Zuds100_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds100_dphiRes_wO_DR07->SetLineWidth(2);
  h_Zuds100_dphiRes_wO_DR07->SetMinimum(0);
  h_Zuds100_dphiRes_wO_DR07->SetMaximum(3.5);
  h_Zuds100_dphiRes_wO_DR07->GetYaxis()->SetTitle("RMS_{90}#Delta#phi(j_{R},j_{G}) [#circ]");
  h_Zuds100_dphiRes_wO_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_200");
 //h_Zuds100_dphiRes_wO_DR07->SetLineColor(kRed);
  h_Zuds200_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds200_dphiRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_380");
 //h_Zuds380_dphiRes_wO_DR07->SetLineColor(kRed);
  h_Zuds380_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds380_dphiRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_500");
  h_Zuds500_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds500_dphiRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds1500_dphiRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dphiRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_RMS90VsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_dphiRes_wO_DR07->SetLineStyle(1);
  h_Zuds3000_dphiRes_wO_DR07->SetLineWidth(2);

  h_Zuds100_dphiRes_wO_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dphiRes_wO_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dphiRes_wO_DR07->SetLineColor(kOrange);
  h_Zuds500_dphiRes_wO_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dphiRes_wO_DR07->SetLineColor(kBlue);
  h_Zuds3000_dphiRes_wO_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dphiRes_wO_DR07_RMS90_CT_wO_fancy = setUpperCanvas("resolutionGraphCanvas_dphiRes_wO_DR07_RMS90_CT_wO_fancy");
  //resolutionGraphCanvas_dphiRes_wO_DR07_RMS90_CT_wO_fancy->cd();
  //TLegend *leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary = resolutionGraphCanvas_dphiRes_wO_DR07_RMS90_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetBorderSize(0);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextAlign(12);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextSize(0.050);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextFont(42);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetMargin(0.15);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineColor(1);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineStyle(1);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineWidth(1);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetFillColor(0);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetFillStyle(1001);
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->SetHeader("VLC7 Jets, with 3TeV BG");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds100_dphiRes_wO_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds200_dphiRes_wO_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds380_dphiRes_wO_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds500_dphiRes_wO_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds1500_dphiRes_wO_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds3000_dphiRes_wO_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dphiRes_wO_DR07_RMS90_CT_wO_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());



  TH1F* h_Zuds100_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_dphiRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_dphiRes_wO_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_dphiRes_wO_FitCB_DR07->SetMinimum(0);
  h_Zuds100_dphiRes_wO_FitCB_DR07->SetMaximum(3.5);
  h_Zuds100_dphiRes_wO_FitCB_DR07->GetYaxis()->SetTitle("#sigma(#Delta#phi(j_{R},j_{G})) [#circ]");
  h_Zuds100_dphiRes_wO_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_200");
 //h_Zuds100_dphiRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_dphiRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_380");
 //h_Zuds380_dphiRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_dphiRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_500");
  h_Zuds500_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_dphiRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_dphiRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dphiRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dphiRes_SigmaCBVsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_dphiRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_dphiRes_wO_FitCB_DR07->SetLineWidth(2);

  h_Zuds100_dphiRes_wO_FitCB_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dphiRes_wO_FitCB_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dphiRes_wO_FitCB_DR07->SetLineColor(kOrange);
  h_Zuds500_dphiRes_wO_FitCB_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dphiRes_wO_FitCB_DR07->SetLineColor(kBlue);
  h_Zuds3000_dphiRes_wO_FitCB_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dphiRes_wO_DR07_SigmaCB_CT_wO_fancy = setUpperCanvas("resolutionGraphCanvas_dphiRes_wO_DR07_SigmaCB_CT_wO_fancy");
  //resolutionGraphCanvas_dphiRes_wO_DR07_SigmaCB_CT_wO_fancy->cd();
  //TLegend *leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary = resolutionGraphCanvas_dphiRes_wO_DR07_SigmaCB_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetBorderSize(0);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextAlign(12);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextSize(0.050);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextFont(42);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetMargin(0.15);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineColor(1);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineStyle(1);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineWidth(1);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetFillColor(0);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetFillStyle(1001);
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetHeader("VLC7 Jets, with 3TeV BG");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds100_dphiRes_wO_FitCB_DR07->DrawCopy("h,e"),"#approx 50 GeV");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds200_dphiRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds380_dphiRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds500_dphiRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds1500_dphiRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds3000_dphiRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dphiRes_wO_DR07_SigmaCB_CT_wO_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_Zuds100_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_dthetaRes_wO_DR07->SetLineColor(kRed);
  h_Zuds100_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds100_dthetaRes_wO_DR07->SetLineWidth(2);
  h_Zuds100_dthetaRes_wO_DR07->SetMinimum(0);
  h_Zuds100_dthetaRes_wO_DR07->SetMaximum(1.75);
  h_Zuds100_dthetaRes_wO_DR07->GetYaxis()->SetTitle("RMS_{90}#Delta#theta(j_{R},j_{G}) [#circ]");
  h_Zuds100_dthetaRes_wO_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_200");
 //h_Zuds100_dthetaRes_wO_DR07->SetLineColor(kRed);
  h_Zuds200_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds200_dthetaRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_380");
 //h_Zuds380_dthetaRes_wO_DR07->SetLineColor(kRed);
  h_Zuds380_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds380_dthetaRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_500");
  h_Zuds500_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds500_dthetaRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds1500_dthetaRes_wO_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dthetaRes_wO_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_RMS90VsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_dthetaRes_wO_DR07->SetLineStyle(1);
  h_Zuds3000_dthetaRes_wO_DR07->SetLineWidth(2);


  h_Zuds100_dthetaRes_wO_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dthetaRes_wO_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dthetaRes_wO_DR07->SetLineColor(kOrange);
  h_Zuds500_dthetaRes_wO_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dthetaRes_wO_DR07->SetLineColor(kBlue);
  h_Zuds3000_dthetaRes_wO_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dthetaRes_wO_DR07_RMS90_CT_wO_fancy = setUpperCanvas("resolutionGraphCanvas_dthetaRes_wO_DR07_RMS90_CT_wO_fancy");
  //resolutionGraphCanvas_dthetaRes_wO_DR07_RMS90_CT_wO_fancy->cd();
  //TLegend *leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary = resolutionGraphCanvas_dthetaRes_wO_DR07_RMS90_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetBorderSize(0);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextAlign(12);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextSize(0.050);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetTextFont(42);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetMargin(0.15);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineColor(1);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineStyle(1);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetLineWidth(1);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetFillColor(0);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetFillStyle(1001);
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->SetHeader("VLC7 Jets, with 3 TeV BG");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds100_dthetaRes_wO_DR07->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds200_dthetaRes_wO_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds380_dthetaRes_wO_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds500_dthetaRes_wO_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds1500_dthetaRes_wO_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->AddEntry(h_Zuds3000_dthetaRes_wO_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dthetaRes_wO_DR07_RMS90_CT_wO_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TH1F* h_Zuds100_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_100");
 //h_Zuds100_dthetaRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds100_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds100_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);
  h_Zuds100_dthetaRes_wO_FitCB_DR07->SetMinimum(0);
  h_Zuds100_dthetaRes_wO_FitCB_DR07->SetMaximum(1.75);
  h_Zuds100_dthetaRes_wO_FitCB_DR07->GetYaxis()->SetTitle("#sigma(#Delta#theta(j_{R},j_{G})) [#circ]");
  h_Zuds100_dthetaRes_wO_FitCB_DR07->GetXaxis()->SetTitle("|cos#theta|");
  TH1F* h_Zuds200_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_200");
 //h_Zuds100_dthetaRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds200_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds200_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds380_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_380");
 //h_Zuds380_dthetaRes_wO_FitCB_DR07->SetLineColor(kRed);
  h_Zuds380_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds380_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds500_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_500");
  h_Zuds500_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds500_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds1500_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_1500");
  h_Zuds1500_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds1500_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);
  TH1F* h_Zuds3000_dthetaRes_wO_FitCB_DR07 = (TH1F*)file_histos_normalRange_jets_drawbins->Get("dthetaRes_SigmaCBVsCosTheta_CT_wO_DR07_3000");
  h_Zuds3000_dthetaRes_wO_FitCB_DR07->SetLineStyle(1);
  h_Zuds3000_dthetaRes_wO_FitCB_DR07->SetLineWidth(2);

  h_Zuds100_dthetaRes_wO_FitCB_DR07->SetLineColor(kCyan+1);
  h_Zuds200_dthetaRes_wO_FitCB_DR07->SetLineColor(kYellow+1);
  h_Zuds380_dthetaRes_wO_FitCB_DR07->SetLineColor(kOrange);
  h_Zuds500_dthetaRes_wO_FitCB_DR07->SetLineColor(kGreen-2);
  h_Zuds1500_dthetaRes_wO_FitCB_DR07->SetLineColor(kBlue);
  h_Zuds3000_dthetaRes_wO_FitCB_DR07->SetLineColor(kRed);

  TCanvas *resolutionGraphCanvas_dthetaRes_wO_DR07_SigmaCB_CT_wO_fancy = setUpperCanvas("resolutionGraphCanvas_dthetaRes_wO_DR07_SigmaCB_CT_wO_fancy");
  //resolutionGraphCanvas_dthetaRes_wO_DR07_SigmaCB_CT_wO_fancy->cd();
  //TLegend *leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary = resolutionGraphCanvas_dthetaRes_wO_DR07_SigmaCB_CT_wO_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetBorderSize(0);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextAlign(12);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextSize(0.050);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetTextFont(42);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetMargin(0.15);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineColor(1);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineStyle(1);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetLineWidth(1);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetFillColor(0);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetFillStyle(1001);
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->SetHeader("VLC7 Jets, with 3 TeV BG");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds100_dthetaRes_wO_FitCB_DR07->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds200_dthetaRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 100 GeV");
  //leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds380_dthetaRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 190 GeV");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds500_dthetaRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 250 GeV");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds1500_dthetaRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 750 GeV");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->AddEntry(h_Zuds3000_dthetaRes_wO_FitCB_DR07->DrawCopy("h,e,same"),"#approx 1500 GeV");
  leg_dthetaRes_wO_DR07_SigmaCB_CT_wO_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  h_Zuds500_dphi_0_to_065->SetLineColor(kBlack);
  h_Zuds500_wO_dphi_0_to_065->Scale(h_Zuds500_dphi_0_to_065->Integral()/h_Zuds500_wO_dphi_0_to_065->Integral());

  TCanvas *resolutionGraphCanvas_dphiRes_Zuds500_DR07_FitCB_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dphiRes_Zuds500_DR07_FitCB_CT_fancy");
  //resolutionGraphCanvas_dphiRes_Zuds500_DR07_FitCB_CT_fancy->cd();
  //TLegend *leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary = resolutionGraphCanvas_dphiRes_Zuds500_DR07_FitCB_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetBorderSize(0);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextAlign(12);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextSize(0.050);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextFont(42);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetMargin(0.15);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineColor(1);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineStyle(1);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineWidth(1);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetFillColor(0);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetFillStyle(1001);
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->SetHeader("VLC7 Jets, #approx250GeV, |cos#theta|<0.65");
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->AddEntry(h_Zuds500_dphi_0_to_065->DrawCopy("h,e"),"no BG");
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->AddEntry(h_Zuds500_wO_dphi_0_to_065->DrawCopy("h,e,same"),"with 3 TeV BG");
  leg_dphiRes_Zuds500_DR07_FitCB_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  h_Zuds500_dtheta_0_to_065->SetLineColor(kBlack);
  h_Zuds500_wO_dtheta_0_to_065->Scale(h_Zuds500_dtheta_0_to_065->Integral()/h_Zuds500_wO_dtheta_0_to_065->Integral());

  TCanvas *resolutionGraphCanvas_dthetaRes_Zuds500_DR07_FitCB_CT_fancy = setUpperCanvas("resolutionGraphCanvas_dthetaRes_Zuds500_DR07_FitCB_CT_fancy");
  //resolutionGraphCanvas_dthetaRes_Zuds500_DR07_FitCB_CT_fancy->cd();
  //TLegend *leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary = resolutionGraphCanvas_dthetaRes_Zuds500_DR07_FitCB_CT_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary = new TLegend(0.20,0.546,0.50,0.87);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetBorderSize(0);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextAlign(12);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextSize(0.050);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetTextFont(42);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetMargin(0.15);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineColor(1);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineStyle(1);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetLineWidth(1);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetFillColor(0);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetFillStyle(1001);
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->SetHeader("VLC7 Jets, #approx250GeV, |cos#theta|<0.65");
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->AddEntry(h_Zuds500_dtheta_0_to_065->DrawCopy("h,e"),"no BG");
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->AddEntry(h_Zuds500_wO_dtheta_0_to_065->DrawCopy("h,e,same"),"with 3 TeV BG");
  leg_dthetaRes_Zuds500_DR07_FitCB_CT_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());



  TFile* file_histo_Jets380=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181101_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi300_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");

  //afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_wO_181011_histoFitRange_0_00_to_3_00_JetAlgorithms_angMatch_10Deg_for_dphi_dtheta_gj1_gj2_dphi_2_80_DrawManyBinsJER200dphi200_range10_WW_ZZ_mass_histos_pyt380_MeanReal_addNoSC_SetNpx_largeStat.root");



  TH1F* h_JER_Zuds100_whz = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_DR07_100");
  h_JER_Zuds100_whz->SetLineStyle(1);
  h_JER_Zuds100_whz->SetLineColor(kCyan+1);
  h_JER_Zuds100_whz->SetLineWidth(2);
  h_JER_Zuds100_whz->SetMinimum(1.75);
  h_JER_Zuds100_whz->SetMaximum(10.1);
  h_JER_Zuds100_whz->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_Zuds100_whz->GetXaxis()->SetTitle("|cos#theta|");          
  TH1F* h_JER_Zuds100_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_wO_380_DR07_100");
  h_JER_Zuds100_whz_wO_380->SetLineStyle(2);
  h_JER_Zuds100_whz_wO_380->SetLineColor(kCyan+1);
  h_JER_Zuds100_whz_wO_380->SetLineWidth(2);
  h_JER_Zuds100_whz_wO_380->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_Zuds100_whz_wO_380->GetXaxis()->SetTitle("|cos#theta|");         
                                                   
  TH1F* h_JER_Zuds200_whz = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_DR07_200");
  h_JER_Zuds200_whz->SetLineStyle(1);
  h_JER_Zuds200_whz->SetLineColor(kYellow+1);
  h_JER_Zuds200_whz->SetLineWidth(2);                              
  TH1F* h_JER_Zuds200_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_wO_380_DR07_200");
  h_JER_Zuds200_whz_wO_380->SetLineStyle(2);
  h_JER_Zuds200_whz_wO_380->SetLineColor(kYellow+1);
  h_JER_Zuds200_whz_wO_380->SetLineWidth(2);

  TH1F* h_JER_Zuds380_whz = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_DR07_380");
  h_JER_Zuds380_whz->SetLineStyle(1);
  h_JER_Zuds380_whz->SetLineColor(kBlack);
  h_JER_Zuds380_whz->SetLineWidth(2);
  TH1F* h_JER_Zuds380_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_RMS90VsCosTheta_CT_wO_380_DR07_380");
  h_JER_Zuds380_whz_wO_380->SetLineStyle(2);
  h_JER_Zuds380_whz_wO_380->SetLineWidth(2);
  h_JER_Zuds380_whz_wO_380->SetLineColor(kBlack);

  TCanvas *canvas_Jet_reco_over_gen_DR07_whz_DY_def_and_wO = setUpperCanvas("canvas_Jet_reco_over_gen_DR07_whz_DY_def_and_wO");
  
  TLegend* leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO=new TLegend(0.205,0.54,0.405,0.87);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetBorderSize(0);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetTextAlign(12);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetTextSize(0.050);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetTextFont(42);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetMargin(0.15);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetLineColor(1);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetLineStyle(1);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetLineWidth(1);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetFillColor(0);
  //leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetFillStyle(1001);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetFillStyle(0);
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->SetHeader("VLC7 Jets");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds100_whz->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds100_whz_wO_380->DrawCopy("hist,e,same"),"#approx 50 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds200_whz->DrawCopy("hist,e,same"),"#approx 100 GeV");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds200_whz_wO_380->DrawCopy("hist,e,same"),"#approx 100 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds380_whz->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->AddEntry(h_JER_Zuds380_whz_wO_380->DrawCopy("hist,e,same"),"#approx 190 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_DR07_whz_DY_def_and_wO->Draw();

  l->DrawLatex(x,y,label.c_str());


  TH1F* h_JER_FitCB_Zuds100_whz = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_DR07_100");
  h_JER_FitCB_Zuds100_whz->SetLineStyle(1);
  h_JER_FitCB_Zuds100_whz->SetLineColor(kCyan+1);
  h_JER_FitCB_Zuds100_whz->SetLineWidth(2);
  h_JER_FitCB_Zuds100_whz->SetMinimum(1.75);
  h_JER_FitCB_Zuds100_whz->SetMaximum(10.1);
  h_JER_FitCB_Zuds100_whz->GetYaxis()->SetTitle("#sigma(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_FitCB_Zuds100_whz->GetXaxis()->SetTitle("|cos#theta|");          
  TH1F* h_JER_FitCB_Zuds100_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_wO_380_DR07_100");
  h_JER_FitCB_Zuds100_whz_wO_380->SetLineStyle(2);
  h_JER_FitCB_Zuds100_whz_wO_380->SetLineColor(kCyan+1);
  h_JER_FitCB_Zuds100_whz_wO_380->SetLineWidth(2);
  h_JER_FitCB_Zuds100_whz_wO_380->GetYaxis()->SetTitle("#sigma(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_FitCB_Zuds100_whz_wO_380->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  h_JER_FitCB_Zuds100_whz_wO_380->GetXaxis()->SetTitle("|cos#theta|");         
                                                   
  TH1F* h_JER_FitCB_Zuds200_whz = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_DR07_200");
  h_JER_FitCB_Zuds200_whz->SetLineStyle(1);
  h_JER_FitCB_Zuds200_whz->SetLineColor(kYellow+1);
  h_JER_FitCB_Zuds200_whz->SetLineWidth(2);                              
  TH1F* h_JER_FitCB_Zuds200_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_wO_380_DR07_200");
  h_JER_FitCB_Zuds200_whz_wO_380->SetLineStyle(2);
  h_JER_FitCB_Zuds200_whz_wO_380->SetLineColor(kYellow+1);
  h_JER_FitCB_Zuds200_whz_wO_380->SetLineWidth(2);

  TH1F* h_JER_FitCB_Zuds380_whz = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_DR07_380");
  h_JER_FitCB_Zuds380_whz->SetLineStyle(1);
  h_JER_FitCB_Zuds380_whz->SetLineColor(kBlack);
  h_JER_FitCB_Zuds380_whz->SetLineWidth(2);
  TH1F* h_JER_FitCB_Zuds380_whz_wO_380 = (TH1F*)file_histo_Jets380->Get("JER_SigmaCBVsCosTheta_CT_wO_380_DR07_380");
  h_JER_FitCB_Zuds380_whz_wO_380->SetLineStyle(2);
  h_JER_FitCB_Zuds380_whz_wO_380->SetLineWidth(2);
  h_JER_FitCB_Zuds380_whz_wO_380->SetLineColor(kBlack);

  TCanvas *canvas_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO = setUpperCanvas("canvas_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO");
  
  TLegend* leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO=new TLegend(0.205,0.54,0.405,0.87);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetBorderSize(0);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetTextAlign(12);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetTextSize(0.050);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetTextFont(42);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetMargin(0.15);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetLineColor(1);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetLineStyle(1);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetLineWidth(1);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetFillColor(0);
  //leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetFillStyle(1001);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetFillStyle(0);
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->SetHeader("VLC7 Jets");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds100_whz->DrawCopy("hist,e"),"#approx 50 GeV");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds100_whz_wO_380->DrawCopy("hist,e,same"),"#approx 50 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds200_whz->DrawCopy("hist,e,same"),"#approx 100 GeV");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds200_whz_wO_380->DrawCopy("hist,e,same"),"#approx 100 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds380_whz->DrawCopy("hist,e,same"),"#approx 190 GeV");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->AddEntry(h_JER_FitCB_Zuds380_whz_wO_380->DrawCopy("hist,e,same"),"#approx 190 GeV, 380 GeV BG");
  leg_Jet_reco_over_gen_FitCB_DR07_whz_DY_def_and_wO->Draw();

  l->DrawLatex(x,y,label.c_str());

  TCanvas *resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_65 = setUpperCanvas("resolutionGraphCanvas_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65");
  resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_65->cd();
  TGraphErrors* TG_JER_0_65_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_wSC");
  if(TG_JER_0_65_wSC!=NULL){
    std::cout<<"at least graph should exist"<<std::endl;
  }
  TG_JER_0_65_wSC->GetXaxis()->SetTitle("E_{jet} [GeV]");
  TG_JER_0_65_wSC->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  TG_JER_0_65_wSC->GetYaxis()->SetTitleSize(0.06);
  //TG_JER_0_65_wSC->GetYaxis()->SetRangeUser(1.75,10.1);
  TG_JER_0_65_wSC->SetMaximum(10.5);
  TG_JER_0_65_wSC->SetMinimum(2.4);

  TGraphErrors* TG_JER_0_65_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_noSC");
  TGraphErrors* TG_JER_0_65_wO_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_wO_wSC");
  TGraphErrors* TG_JER_0_65_wO_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_wO_noSC");

  TG_JER_0_65_wSC->Draw("PA");
  TG_JER_0_65_wO_wSC->Draw("P,same");
  TG_JER_0_65_wO_noSC->Draw("P,same");
  TG_JER_0_65_noSC->Draw("P,same");

  
  //resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_fancy->cd();
  //TLegend *leg_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_FullSummary = resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary = new TLegend(0.30,0.546,0.65,0.87);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetBorderSize(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetTextAlign(12);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetTextSize(0.050);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetTextFont(42);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetMargin(0.15);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetLineColor(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetLineStyle(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetLineWidth(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetFillColor(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetFillStyle(1001);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->SetHeader("VLC7 Jets, |cos#theta|<0.65");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->AddEntry(TG_JER_0_65_wSC,"with software comp.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->AddEntry(TG_JER_0_65_noSC,"no energy corr.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->AddEntry(TG_JER_0_65_wO_wSC,"with software comp., 3TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->AddEntry(TG_JER_0_65_wO_noSC,"no energy corr., 3 TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


  TCanvas *resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_to_0_80 = setUpperCanvas("resolutionGraphCanvas_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80");
  resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_to_0_80->cd();
  TGraphErrors* TG_JER_0_65_to_0_80_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_to_0_80_wSC");
  if(TG_JER_0_65_to_0_80_wSC!=NULL){
    std::cout<<"at least graph should exist"<<std::endl;
  }
  TG_JER_0_65_to_0_80_wSC->GetXaxis()->SetTitle("E_{jet} [GeV]");
  TG_JER_0_65_to_0_80_wSC->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  TG_JER_0_65_to_0_80_wSC->GetYaxis()->SetTitleSize(0.06);
  //TG_JER_0_65_to_0_80_wSC->GetYaxis()->SetRangeUser(1.75,10.1);
  TG_JER_0_65_to_0_80_wSC->SetMaximum(10.5);
  TG_JER_0_65_to_0_80_wSC->SetMinimum(2.4);

  TGraphErrors* TG_JER_0_65_to_0_80_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_to_0_80_noSC");
  TGraphErrors* TG_JER_0_65_to_0_80_wO_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_to_0_80_wO_wSC");
  TGraphErrors* TG_JER_0_65_to_0_80_wO_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_65_to_0_80_wO_noSC");

  TG_JER_0_65_to_0_80_wSC->Draw("PA");
  TG_JER_0_65_to_0_80_wO_wSC->Draw("P,same");
  TG_JER_0_65_to_0_80_wO_noSC->Draw("P,same");
  TG_JER_0_65_to_0_80_noSC->Draw("P,same");

  
  //resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_fancy->cd();
  //TLegend *leg_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary = resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary = new TLegend(0.30,0.546,0.65,0.87);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetBorderSize(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetTextAlign(12);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetTextSize(0.050);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetTextFont(42);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetMargin(0.15);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetLineColor(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetLineStyle(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetLineWidth(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetFillColor(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetFillStyle(1001);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->SetHeader("VLC7 Jets, 0.65<|cos#theta|<0.80");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->AddEntry(TG_JER_0_65_to_0_80_wSC,"with software comp.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->AddEntry(TG_JER_0_65_to_0_80_noSC,"no energy corr.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->AddEntry(TG_JER_0_65_to_0_80_wO_wSC,"with software comp., 3TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->AddEntry(TG_JER_0_65_to_0_80_wO_noSC,"no energy corr., 3 TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_65_to_0_80_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TCanvas *resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925 = setUpperCanvas("resolutionGraphCanvas_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925");
  resolutionGraphCanvas_JER_100_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925->cd();
  TGraphErrors* TG_JER_0_80_to_0_925_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_80_to_0_925_wSC");
  if(TG_JER_0_80_to_0_925_wSC!=NULL){
    std::cout<<"at least graph should exist"<<std::endl;
  }
  TG_JER_0_80_to_0_925_wSC->GetXaxis()->SetTitle("E_{jet} [GeV]");
  TG_JER_0_80_to_0_925_wSC->GetYaxis()->SetTitle("RMS_{90}(E_{j}^{R}/E_{j}^{G})/Mean_{90}(E_{j}^{R}/E_{j}^{G})[%]");
  TG_JER_0_80_to_0_925_wSC->GetYaxis()->SetTitleSize(0.06);
  //TG_JER_0_80_to_0_925_wSC->GetYaxis()->SetRangeUser(1.75,10.1);
  TG_JER_0_80_to_0_925_wSC->SetMaximum(10.50);
  TG_JER_0_80_to_0_925_wSC->SetMinimum(2.4);

  TGraphErrors* TG_JER_0_80_to_0_925_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_80_to_0_925_noSC");
  TGraphErrors* TG_JER_0_80_to_0_925_wO_wSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_80_to_0_925_wO_wSC");
  TGraphErrors* TG_JER_0_80_to_0_925_wO_noSC=(TGraphErrors*) file_histos_normalRange_jets->Get("gre_DR07_RMS90_JER_0_80_to_0_925_wO_noSC");

  TG_JER_0_80_to_0_925_wSC->Draw("PA");
  TG_JER_0_80_to_0_925_wO_wSC->Draw("P,same");
  TG_JER_0_80_to_0_925_wO_noSC->Draw("P,same");
  TG_JER_0_80_to_0_925_noSC->Draw("P,same");

  
  //resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_fancy->cd();
  //TLegend *leg_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary = resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary = new TLegend(0.30,0.546,0.65,0.87);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetBorderSize(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetTextAlign(12);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetTextSize(0.050);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetTextFont(42);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetMargin(0.15);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetLineColor(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetLineStyle(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetLineWidth(1);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetFillColor(0);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetFillStyle(1001);
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->SetHeader("VLC7 Jets, 0.80<|cos#theta|<0.925");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->AddEntry(TG_JER_0_80_to_0_925_wSC,"with software comp.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->AddEntry(TG_JER_0_80_to_0_925_noSC,"no energy corr.","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->AddEntry(TG_JER_0_80_to_0_925_wO_wSC,"with software comp., 3TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->AddEntry(TG_JER_0_80_to_0_925_wO_noSC,"no energy corr., 3 TeV BG","PE");
  leg_JER_100_3000_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());

  TFile* file_totE_JER=TFile::Open("/afs/cern.ch/user/w/weberma2/performanceHistoFiles/Zuds_CT_181101_correctMC_histoFitRange_0_00_to_2_00_RedXRange_centralQuark_cosTheta_0_70_RealMeanHisto_RMS90_and_RMS_redBinsForGaussianFit.root");



  TCanvas *resolutionGraphCanvas_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70 = setUpperCanvas("resolutionGraphCanvas_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70");
  resolutionGraphCanvas_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70->cd();
  TGraphErrors* TG_JER_totE_0_70_wSC=(TGraphErrors*)file_totE_JER->Get("gre_DR07_RMS90_JER_0_70_wSC");
  if(TG_JER_totE_0_70_wSC!=NULL){
    std::cout<<"at least graph should exist"<<std::endl;
  }
  TG_JER_totE_0_70_wSC->GetXaxis()->SetTitle("E_{jet} [GeV]");
  TG_JER_totE_0_70_wSC->GetYaxis()->SetTitle("RMS_{90}(E_{j})/Mean_{90}(E_{j})[%]");
  TG_JER_totE_0_70_wSC->GetYaxis()->SetTitleSize(0.06);
  //TG_JER_totE_0_70_wSC->GetYaxis()->SetRangeUser(1.75,10.1);
  TG_JER_totE_0_70_wSC->SetMaximum(5.50);
  TG_JER_totE_0_70_wSC->SetMinimum(2.4);
  TG_JER_totE_0_70_wSC->Draw("PA");

  TGraphErrors* TG_JER_totE_0_70_noSC=(TGraphErrors*) file_totE_JER->Get("gre_DR07_RMS90_JER_0_70_noSC");
  TG_JER_totE_0_70_noSC->Draw("P,same");

  //resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_fancy->cd();
  //TLegend *leg_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_FullSummary = resolutionGraphCanvas_JER_3000_wO_DR07_RMS90_CT_SC_comparison_0_80_to_0_925_fancy->BuildLegend(0.20,0.65,0.50,0.85);
  TLegend *leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary = new TLegend(0.30,0.546,0.65,0.87);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetBorderSize(0);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetTextAlign(12);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetTextSize(0.050);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetTextFont(42);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetMargin(0.15);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetLineColor(1);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetLineStyle(1);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetLineWidth(1);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetFillColor(0);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetFillStyle(1001);
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->SetHeader("VLC7 Jets, 0.80<|cos#theta|<0.925");
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->AddEntry(TG_JER_totE_0_70_wSC,"with software comp.","PE");
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->AddEntry(TG_JER_totE_0_70_noSC,"no energy corr.","PE");
  leg_JER_100_3000_totE_RMS90_CT_SC_comparison_0_70_FullSummary->Draw();

  l->DrawLatex(x,y,label.c_str());


}




