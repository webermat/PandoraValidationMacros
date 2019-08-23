#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TLegendEntry.h"
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


void plotECALvsRECOSummary(){


 CLICdpStyle();


//=========Macro generated from canvas: resolutionGraphCanvas/
//=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
   TCanvas *resolutionGraphCanvas = new TCanvas("resolutionGraphCanvas", "resolutionGraphCanvas",497,349,698,650);
   gStyle->SetOptStat(0);
   resolutionGraphCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
   resolutionGraphCanvas->SetFillColor(0);
   resolutionGraphCanvas->SetBorderMode(0);
   resolutionGraphCanvas->SetBorderSize(2);
   resolutionGraphCanvas->SetGridx();
   resolutionGraphCanvas->SetGridy();
    resolutionGraphCanvas->SetRightMargin(0.0172);
   resolutionGraphCanvas->SetTopMargin(0.0164);
   resolutionGraphCanvas->SetBottomMargin(0.125);
   resolutionGraphCanvas->SetFrameBorderMode(0);
   resolutionGraphCanvas->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(8);
   gre->SetName("CLIC_o2_v04_CUR_22.85X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetPoint(0,1,14.5097);
   gre->SetPointError(0,0,0.115643);
   gre->SetPoint(1,10,4.66929);
   gre->SetPointError(1,0,0.0357393);
   gre->SetPoint(2,50,2.13389);
   gre->SetPointError(2,0,0.0162179);
   gre->SetPoint(3,100,1.52347);
   gre->SetPointError(3,0,0.0116793);
   gre->SetPoint(4,200,1.09412);
   gre->SetPointError(4,0,0.0085235);
   gre->SetPoint(5,500,0.735121);
   gre->SetPointError(5,0,0.00583176);
   gre->SetPoint(6,1000,0.568574);
   gre->SetPointError(6,0,0.00477162);
   gre->SetPoint(7,1500,0.507686);
   gre->SetPointError(7,0,0.00469561);
   
   TH1F *Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1 = new TH1F("Graph_CLIC_o2_v04_CUR_22.85X0_resolutionGraph1","",100,0,1649.9);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetMinimum(0);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetMaximum(6);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetDirectory(0);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetLineColor(ci);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetLabelFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetLabelSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleOffset(1.18);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetTitle("#sigma(E_{RECO})/E_{true}_{#gamma} [%]");
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetLabelFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetLabelSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetTitleSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetTitleOffset(1.20);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetTitleFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetZaxis()->SetLabelFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetZaxis()->SetLabelSize(0.035);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetZaxis()->SetTitleSize(0.035);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1);
   
   
   TF1 *resFunc = new TF1("resFunc","sqrt(sq([0])+sq([1])/x)",0.5,1600);
   resFunc->SetFillColor(19);
   resFunc->SetFillStyle(0);

   ci = TColor::GetColor("#0000ff");
   resFunc->SetLineColor(ci);
   resFunc->SetLineWidth(1);
   resFunc->SetLineStyle(2);
   resFunc->SetChisquare(77.68741);
   resFunc->SetNDF(3);
   resFunc->GetXaxis()->SetLabelFont(42);
   resFunc->GetXaxis()->SetLabelSize(0.035);
   resFunc->GetXaxis()->SetTitleSize(0.035);
   resFunc->GetXaxis()->SetTitleFont(42);
   resFunc->GetYaxis()->SetLabelFont(42);
   resFunc->GetYaxis()->SetLabelSize(0.035);
   resFunc->GetYaxis()->SetTitleSize(0.035);
   resFunc->GetYaxis()->SetTitleFont(42);
   resFunc->SetParameter(0,0.51234);
   resFunc->SetParError(0,0.00640);
   resFunc->SetParLimits(0,0,0);
   resFunc->SetParameter(1,19.903);
   resFunc->SetParError(1,0.07461);
   resFunc->SetParLimits(1,0,0);
   //gre->GetListOfFunctions()->Add(resFunc);
   gre->Draw("ape");
   

   gre = new TGraphErrors(8);
   gre->SetName("CLICdet30_23.216X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ffff");
   gre->SetLineColor(kOrange-3);

   //ci = TColor::GetColor("#00ffff");
   //ci = TColor::GetColor(kOrange-3);
   gre->SetMarkerColor(kOrange-3);
   gre->SetMarkerStyle(23);
   gre->SetPoint(0,1,17.88717);
   gre->SetPointError(0,0,0.1518231);
   gre->SetPoint(1,10,5.430001);
   gre->SetPointError(1,0,0.04173555);
   gre->SetPoint(2,50,2.439253);
   gre->SetPointError(2,0,0.01851144);
   gre->SetPoint(3,100,1.74868);
   gre->SetPointError(3,0,0.01342834);
   gre->SetPoint(4,200,1.231169);
   gre->SetPointError(4,0,0.00938356);
   gre->SetPoint(5,500,0.789896);
   gre->SetPointError(5,0,0.005813784);
   gre->SetPoint(6,1000,0.5719569);
   gre->SetPointError(6,0,0.004399268);
   gre->SetPoint(7,1500,0.485081);
   gre->SetPointError(7,0,0.003683288);
   
   TH1F *Graph_CLICdet30_23_216X0_resolutionGraph4 = new TH1F("Graph_CLICdet30_23_216X0_resolutionGraph4","",100,0,1649.9);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetMinimum(0);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetMaximum(19.79475);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetDirectory(0);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetStats(0);

   //ci = TColor::GetColor("#000099");
   ci = TColor::GetColor(kOrange-3);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetLineColor(kOrange-3);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetXaxis()->SetLabelFont(42);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetXaxis()->SetLabelSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetXaxis()->SetTitleSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetXaxis()->SetTitleFont(42);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetYaxis()->SetLabelFont(42);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetYaxis()->SetLabelSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetYaxis()->SetTitleSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetYaxis()->SetTitleFont(42);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetZaxis()->SetLabelFont(42);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetZaxis()->SetLabelSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetZaxis()->SetTitleSize(0.035);
   Graph_CLICdet30_23_216X0_resolutionGraph4->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_CLICdet30_23_216X0_resolutionGraph4);
   gre->Draw("pe");

   gre = new TGraphErrors(8);
   gre->SetName("CLICdet40_o3_V05X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);

   ci = TColor::GetColor(kBlack);
   gre->SetLineColor(ci);
   gre->SetMarkerColor(kBlack);
   gre->SetMarkerStyle(22);
   gre->SetPoint(0,1,15.08893);
   gre->SetPointError(0,0,0.1254097);
   gre->SetPoint(1,10,4.635);
   gre->SetPointError(1,0,0.0337);
   gre->SetPoint(2,50,2.084);
   gre->SetPointError(2,0,0.01492);
   gre->SetPoint(3,100,1.4734);
   gre->SetPointError(3,0,0.010492);
   gre->SetPoint(4,200,1.059);
   gre->SetPointError(4,0,0.00769);
   gre->SetPoint(5,500,0.6834);
   gre->SetPointError(5,0,0.00506);
   gre->SetPoint(6,1000,0.4912);
   gre->SetPointError(6,0,0.00347);
   gre->SetPoint(7,1500,0.4078);
   gre->SetPointError(7,0,0.003);
   
   TH1F *Graph_CLICdet40_o3_V05_resolutionGraph6 = new TH1F("Graph_CLICdet40_o3_V05X0_resolutionGraph6","",100,0,1649.9);
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetMinimum(0);
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetMaximum(16.71195);
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetDirectory(0);
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetStats(0);

   ci = TColor::GetColor("#000107");
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetLineColor(ci);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetXaxis()->SetLabelFont(42);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetXaxis()->SetLabelSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetXaxis()->SetTitleSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetXaxis()->SetTitleFont(42);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetYaxis()->SetLabelFont(42);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetYaxis()->SetLabelSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetYaxis()->SetTitleSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetYaxis()->SetTitleFont(42);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetZaxis()->SetLabelFont(42);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetZaxis()->SetLabelSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetZaxis()->SetTitleSize(0.035);
   Graph_CLICdet40_o3_V05_resolutionGraph6->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_CLICdet40_o3_V05_resolutionGraph6);
   
   
   TF1 *resFunc7 = new TF1("resFunc7","sqrt(sq([0])+sq([1])/x)",0.5,1600);
   resFunc7->SetFillColor(19);
   resFunc7->SetFillStyle(0);
   resFunc7->SetLineWidth(1);
   resFunc7->SetLineStyle(2);
   resFunc7->SetChisquare(2.8686);
   resFunc7->SetNDF(3);
   resFunc7->GetXaxis()->SetLabelFont(42);
   resFunc7->GetXaxis()->SetLabelSize(0.035);
   resFunc7->GetXaxis()->SetTitleSize(0.035);
   resFunc7->GetXaxis()->SetTitleFont(42);
   resFunc7->GetYaxis()->SetLabelFont(42);
   resFunc7->GetYaxis()->SetLabelSize(0.035);
   resFunc7->GetYaxis()->SetTitleSize(0.035);
   resFunc7->GetYaxis()->SetTitleFont(42);
   resFunc7->SetParameter(0,0.3505505);
   resFunc7->SetParError(0,0.01042596);
   resFunc7->SetParLimits(0,0,0);
   resFunc7->SetParameter(1,14.72119);
   resFunc7->SetParError(1,0.08035938);
   resFunc7->SetParLimits(1,0,0);
   //gre->GetListOfFunctions()->Add(resFunc7);
   gre->Draw("pe");


   
   TLegend *leg = new TLegend(0.30,0.60,0.94,0.94,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextSize(0.027);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("CLIC_o2_v04_CUR_22.85X0_resolutionGraph","CLICdet_40 PFAPhoton Energies","pe");
   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(kOrange-3);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor(kGreen);
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("CLICdet40_o3_V05X0_resolutionGraph","CLICdet_40 CaloHit Energies","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(14);
   entry->SetMarkerStyle(27);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);

   entry=leg->AddEntry("CLICdet30_23.216X0_resolutionGraph","CLICdet_30","pe");

   leg->Draw();
   /*
   TLatex *   tex = new TLatex(75,5.23,"#splitline{#frac{#sigma(E)}{E}=#sqrt{#alpha^{2}+#(){#frac{#beta}{#sqrt{E}}}^{2}}}{#splitline{Fit between}{80 and 1600 GeV}}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.027);
   tex->SetLineWidth(2);
   tex->Draw();
   resolutionGraphCanvas->Modified();
   resolutionGraphCanvas->cd();
   resolutionGraphCanvas->SetSelected(resolutionGraphCanvas);
   */


   //=========Macro generated from canvas: resolutionGraphCanvas/
   //=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34                                       
   TCanvas *PFOEnergyResolutionGraphCanvas = new TCanvas("PFOEnergyResolutionGraphCanvas", "PFOEnergyResolutionGraphCanvas",0,0,800,700);
   gStyle->SetOptStat(0);
   PFOEnergyResolutionGraphCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
   PFOEnergyResolutionGraphCanvas->SetFillColor(0);
   PFOEnergyResolutionGraphCanvas->SetBorderMode(0);
   PFOEnergyResolutionGraphCanvas->SetBorderSize(2);
   PFOEnergyResolutionGraphCanvas->SetGridx();
   PFOEnergyResolutionGraphCanvas->SetGridy();
   PFOEnergyResolutionGraphCanvas->SetRightMargin(0.0172);
   PFOEnergyResolutionGraphCanvas->SetTopMargin(0.0164);
   PFOEnergyResolutionGraphCanvas->SetBottomMargin(0.125);
   PFOEnergyResolutionGraphCanvas->SetFrameBorderMode(0);
   PFOEnergyResolutionGraphCanvas->SetFrameBorderMode(0);



   TGraphErrors *gre_E_PFO = new TGraphErrors(11);
   gre_E_PFO->SetName("E_ph_PFO_innerBarrel_025_Graph");
   gre_E_PFO->SetTitle("");
   gre_E_PFO->SetFillColor(1);
   gre_E_PFO->SetLineColor(kBlack);
   gre_E_PFO->SetMarkerColor(kBlack);
   gre_E_PFO->SetMarkerStyle(kOpenCircle);
   gre_E_PFO->SetPoint(0,1,14.1481);
   gre_E_PFO->SetPointError(0,0,0.283079);
   gre_E_PFO->SetPoint(1,5,6.59359);
   gre_E_PFO->SetPointError(1,0,0.102231);
   gre_E_PFO->SetPoint(2,10,4.91568);
   gre_E_PFO->SetPointError(2,0,0.0374898);
   gre_E_PFO->SetPoint(3,15,3.90517);
   gre_E_PFO->SetPointError(3,0,0.0668372);
   gre_E_PFO->SetPoint(4,30,2.88762);
   gre_E_PFO->SetPointError(4,0,0.04636);
   gre_E_PFO->SetPoint(5,50,2.2835);
   gre_E_PFO->SetPointError(5,0,0.0345492);
   gre_E_PFO->SetPoint(6,100,1.69843);
   gre_E_PFO->SetPointError(6,0,0.0285563);
   gre_E_PFO->SetPoint(7,200,1.27415);
   gre_E_PFO->SetPointError(7,0,0.0201239);
   gre_E_PFO->SetPoint(8,500,1.02881);
   gre_E_PFO->SetPointError(8,0,0.0167799);
   gre_E_PFO->SetPoint(9,1000,0.888444);
   gre_E_PFO->SetPointError(9,0,0.0146771);
   gre_E_PFO->SetPoint(10,1500,0.802634);
   gre_E_PFO->SetPointError(10,0,0.0138609);
   
   TH1F *hist_E_ph_PFO_innerBarrel_025 = new TH1F("hist_E_ph_PFO_innerBarrel_025","",100,0,1649.9);
   hist_E_ph_PFO_innerBarrel_025->SetMinimum(0);
   hist_E_ph_PFO_innerBarrel_025->SetMaximum(6);
   hist_E_ph_PFO_innerBarrel_025->SetDirectory(0);
   hist_E_ph_PFO_innerBarrel_025->SetStats(0);

   hist_E_ph_PFO_innerBarrel_025->SetLineColor(kBlack);
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelFont(42);
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelSize(0.05);
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleSize(0.05);
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleFont(42);
   hist_E_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleOffset(1.18);
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetTitle("#sigma(E_{RECO}^{PFO})/E_{true}_{#gamma} [%]");
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelFont(42);
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelSize(0.05);
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleSize(0.05);
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleOffset(1.20);
   hist_E_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleFont(42);
   hist_E_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelFont(42);
   hist_E_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleFont(42);
   gre_E_PFO->SetHistogram(hist_E_ph_PFO_innerBarrel_025);
   
   gre_E_PFO->Draw("ape");
   
   gre_E_PFO = new TGraphErrors(11);
   gre_E_PFO->SetName("Graph_E_ph_PFO_mediumBarrel_025_050");
   gre_E_PFO->SetTitle("");
   gre_E_PFO->SetLineColor(kRed-7);

   gre_E_PFO->SetMarkerColor(kRed-7);
   gre_E_PFO->SetMarkerStyle(kOpenSquare);

   gre_E_PFO->SetPoint(0,1,15.2042);
   gre_E_PFO->SetPointError(0,0,0.279179);
   gre_E_PFO->SetPoint(1,5,7.05405);
   gre_E_PFO->SetPointError(1,0,0.11127);
   gre_E_PFO->SetPoint(2,10,4.97435);
   gre_E_PFO->SetPointError(2,0,0.036666);
   gre_E_PFO->SetPoint(3,15,4.07722);
   gre_E_PFO->SetPointError(3,0,0.0669813);
   gre_E_PFO->SetPoint(4,30,2.99291);
   gre_E_PFO->SetPointError(4,0,0.048517);
   gre_E_PFO->SetPoint(5,50,2.33379);
   gre_E_PFO->SetPointError(5,0,0.0407161);
   gre_E_PFO->SetPoint(6,100,1.66169);
   gre_E_PFO->SetPointError(6,0,0.0270955);
   gre_E_PFO->SetPoint(7,200,1.29463);
   gre_E_PFO->SetPointError(7,0,0.0202295);
   gre_E_PFO->SetPoint(8,500,0.968115);
   gre_E_PFO->SetPointError(8,0,0.0154271);
   gre_E_PFO->SetPoint(9,1000,0.816857);
   gre_E_PFO->SetPointError(9,0,0.0127141);
   gre_E_PFO->SetPoint(10,1500,0.74401);
   gre_E_PFO->SetPointError(10,0,0.0105428);

   TH1F *hist_E_ph_PFO_mediumBarrel_025_050 = new TH1F("hist_E_ph_PFO_mediumBarrel_025_050","",100,0,1649.9);
   hist_E_ph_PFO_mediumBarrel_025_050->SetMinimum(0);
   hist_E_ph_PFO_mediumBarrel_025_050->SetMaximum(19.79475);
   hist_E_ph_PFO_mediumBarrel_025_050->SetDirectory(0);
   hist_E_ph_PFO_mediumBarrel_025_050->SetStats(0);

   hist_E_ph_PFO_mediumBarrel_025_050->SetLineColor(kRed-7);
   hist_E_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelFont(42);
   hist_E_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleFont(42);
   hist_E_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelFont(42);
   hist_E_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleFont(42);
   hist_E_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelFont(42);
   hist_E_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleFont(42);
   gre_E_PFO->SetHistogram(hist_E_ph_PFO_mediumBarrel_025_050);
   gre_E_PFO->Draw("pe");

  
   gre_E_PFO = new TGraphErrors(11);
   gre_E_PFO->SetName("Graph_E_ph_PFO_outerBarrel_050_075");
   gre_E_PFO->SetTitle("");
   gre_E_PFO->SetLineColor(kBlue);

   gre_E_PFO->SetMarkerColor(kBlue);
   gre_E_PFO->SetMarkerStyle(kOpenTriangleUp);
   gre_E_PFO->SetPoint(0,1,15.409);
   gre_E_PFO->SetPointError(0,0,0.289162);
   gre_E_PFO->SetPoint(1,5,7.38227);
   gre_E_PFO->SetPointError(1,0,0.122999);
   gre_E_PFO->SetPoint(2,10,5.00032);
   gre_E_PFO->SetPointError(2,0,0.0333707);
   gre_E_PFO->SetPoint(3,15,4.30284);
   gre_E_PFO->SetPointError(3,0,0.0656688);
   gre_E_PFO->SetPoint(4,30,3.03287);
   gre_E_PFO->SetPointError(4,0,0.0456417);
   gre_E_PFO->SetPoint(5,50,2.343);
   gre_E_PFO->SetPointError(5,0,0.0381346);
   gre_E_PFO->SetPoint(6,100,1.77049);
   gre_E_PFO->SetPointError(6,0,0.0289073);
   gre_E_PFO->SetPoint(7,200,1.32186);
   gre_E_PFO->SetPointError(7,0,0.0209187);
   gre_E_PFO->SetPoint(8,500,0.933383);
   gre_E_PFO->SetPointError(8,0,0.0144638);
   gre_E_PFO->SetPoint(9,1000,0.770742);
   gre_E_PFO->SetPointError(9,0,0.0114544);
   gre_E_PFO->SetPoint(10,1500,0.715001);
   gre_E_PFO->SetPointError(10,0,0.0107809);
 
   TH1F *hist_E_ph_PFO_outerBarrel_050_075 = new TH1F("hist_E_ph_PFO_outerBarrel_050_075","",100,0,1649.9);
   hist_E_ph_PFO_outerBarrel_050_075->SetMinimum(0);
   hist_E_ph_PFO_outerBarrel_050_075->SetMaximum(19.79475);
   hist_E_ph_PFO_outerBarrel_050_075->SetDirectory(0);
   hist_E_ph_PFO_outerBarrel_050_075->SetStats(0);

   hist_E_ph_PFO_outerBarrel_050_075->SetLineColor(kBlue);
   hist_E_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelFont(42);
   hist_E_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleFont(42);
   hist_E_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelFont(42);
   hist_E_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleFont(42);
   hist_E_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelFont(42);
   hist_E_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleFont(42);
   gre_E_PFO->SetHistogram(hist_E_ph_PFO_outerBarrel_050_075);
   gre_E_PFO->Draw("pe");

 
   gre_E_PFO = new TGraphErrors(11);
   gre_E_PFO->SetName("Graph_E_ph_PFO_Transition_075_087");
   gre_E_PFO->SetTitle("");
   gre_E_PFO->SetLineColor(kGreen+2);

   gre_E_PFO->SetMarkerColor(kGreen+2);
   gre_E_PFO->SetMarkerStyle(kOpenDiamond);
   gre_E_PFO->SetPoint(0,1,16.3793);
   gre_E_PFO->SetPointError(0,0,0.488282);
   gre_E_PFO->SetPoint(1,5,7.61453);
   gre_E_PFO->SetPointError(1,0,0.202691);
   gre_E_PFO->SetPoint(2,10,5.23531);
   gre_E_PFO->SetPointError(2,0,0.0579319);
   gre_E_PFO->SetPoint(3,15,4.45309);
   gre_E_PFO->SetPointError(3,0,0.117852);
   gre_E_PFO->SetPoint(4,30,3.1883);
   gre_E_PFO->SetPointError(4,0,0.0883245);
   gre_E_PFO->SetPoint(5,50,2.44644);
   gre_E_PFO->SetPointError(5,0,0.0632229);
   gre_E_PFO->SetPoint(6,100,1.72244);
   gre_E_PFO->SetPointError(6,0,0.0441727);
   gre_E_PFO->SetPoint(7,200,1.29086);
   gre_E_PFO->SetPointError(7,0,0.0387336);
   gre_E_PFO->SetPoint(8,500,0.891414);
   gre_E_PFO->SetPointError(8,0,0.0234621);
   gre_E_PFO->SetPoint(9,1000,0.676626);
   gre_E_PFO->SetPointError(9,0,0.0185027);
   gre_E_PFO->SetPoint(10,1500,0.622652);
   gre_E_PFO->SetPointError(10,0,0.0155014);
 
   TH1F *hist_E_ph_PFO_Transition_075_087 = new TH1F("hist_E_ph_PFO_Transition_075_087","",100,0,1649.9);
   hist_E_ph_PFO_Transition_075_087->SetMinimum(0);
   hist_E_ph_PFO_Transition_075_087->SetMaximum(19.79475);
   hist_E_ph_PFO_Transition_075_087->SetDirectory(0);
   hist_E_ph_PFO_Transition_075_087->SetStats(0);

   hist_E_ph_PFO_Transition_075_087->SetLineColor(kGreen+2);
   hist_E_ph_PFO_Transition_075_087->GetXaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Transition_075_087->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetXaxis()->SetTitleFont(42);
   hist_E_ph_PFO_Transition_075_087->GetYaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Transition_075_087->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetYaxis()->SetTitleFont(42);
   hist_E_ph_PFO_Transition_075_087->GetZaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Transition_075_087->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Transition_075_087->GetZaxis()->SetTitleFont(42);
   gre_E_PFO->SetHistogram(hist_E_ph_PFO_Transition_075_087);
   gre_E_PFO->Draw("pe");

   gre_E_PFO = new TGraphErrors(11);
   gre_E_PFO->SetName("Graph_E_ph_PFO_Endcap_087_098");
   gre_E_PFO->SetTitle("");
   gre_E_PFO->SetLineColor(kCyan+1);

   gre_E_PFO->SetMarkerColor(kCyan+1);
   gre_E_PFO->SetMarkerStyle(kOpenCross);
   gre_E_PFO->SetPoint(0,1,15.1341);
   gre_E_PFO->SetPointError(0,0,0.383339);
   gre_E_PFO->SetPoint(1,5,6.723);
   gre_E_PFO->SetPointError(1,0,0.183003);
   gre_E_PFO->SetPoint(2,10,4.79806);
   gre_E_PFO->SetPointError(2,0,0.0512466);
   gre_E_PFO->SetPoint(3,15,3.79981);
   gre_E_PFO->SetPointError(3,0,0.0846845);
   gre_E_PFO->SetPoint(4,30,2.86315);
   gre_E_PFO->SetPointError(4,0,0.0786549);
   gre_E_PFO->SetPoint(5,50,2.2102);
   gre_E_PFO->SetPointError(5,0,0.0594837);
   gre_E_PFO->SetPoint(6,100,1.54802);
   gre_E_PFO->SetPointError(6,0,0.0354124);
   gre_E_PFO->SetPoint(7,200,1.13692);
   gre_E_PFO->SetPointError(7,0,0.0301534);
   gre_E_PFO->SetPoint(8,500,0.772743);
   gre_E_PFO->SetPointError(8,0,0.017799);
   gre_E_PFO->SetPoint(9,1000,0.602535);
   gre_E_PFO->SetPointError(9,0,0.0130302);
   gre_E_PFO->SetPoint(10,1500,0.545755);
   gre_E_PFO->SetPointError(10,0,0.013064);
 
   TH1F *hist_E_ph_PFO_Endcap_087_098 = new TH1F("hist_E_ph_PFO_Endcap_087_098","",100,0,1649.9);
   hist_E_ph_PFO_Endcap_087_098->SetMinimum(0);
   hist_E_ph_PFO_Endcap_087_098->SetMaximum(19.79475);
   hist_E_ph_PFO_Endcap_087_098->SetDirectory(0);
   hist_E_ph_PFO_Endcap_087_098->SetStats(0);

   hist_E_ph_PFO_Endcap_087_098->SetLineColor(kCyan+1);
   hist_E_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleFont(42);
   hist_E_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleFont(42);
   hist_E_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelFont(42);
   hist_E_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleFont(42);
   gre_E_PFO->SetHistogram(hist_E_ph_PFO_Endcap_087_098);
   gre_E_PFO->Draw("pe");

   
   TLegend *leg_E_PFO = new TLegend(0.30,0.60,0.94,0.94,NULL,"brNDC");
   leg_E_PFO->SetBorderSize(1);
   leg_E_PFO->SetTextSize(0.027);
   leg_E_PFO->SetLineColor(1);
   leg_E_PFO->SetLineStyle(1);
   leg_E_PFO->SetLineWidth(1);
   leg_E_PFO->SetFillColor(0);
   leg_E_PFO->SetFillStyle(1001);
   TLegendEntry *entry_E_PFO=leg_E_PFO->AddEntry("hist_E_ph_PFO_innerBarrel_025","Photon |Cos#theta|<0.25","pe");
   entry_E_PFO->SetMarkerStyle(kOpenCircle);
   entry_E_PFO->SetMarkerColor(kBlack);
   entry_E_PFO->SetLineColor(kBlack);
   entry_E_PFO->SetLineStyle(1);
   entry_E_PFO->SetLineWidth(2);
   entry_E_PFO=leg_E_PFO->AddEntry("hist_E_ph_PFO_mediumBarrel_025_050","Photon 0.25<|Cos#theta|<0.50","pe");
   entry_E_PFO->SetMarkerStyle(kOpenSquare);
   entry_E_PFO->SetMarkerColor(kRed-7);
   entry_E_PFO->SetLineColor(kRed-7);
   entry_E_PFO->SetLineStyle(1);
   entry_E_PFO->SetLineWidth(2);

   entry_E_PFO=leg_E_PFO->AddEntry("hist_E_ph_PFO_outerBarrel_050_075","Photon 0.50<|Cos#theta|<0.75","pe");
   entry_E_PFO->SetMarkerStyle(kOpenTriangleUp);
   entry_E_PFO->SetMarkerColor(kBlue);
   entry_E_PFO->SetLineColor(kBlue);
   entry_E_PFO->SetLineStyle(1);
   entry_E_PFO->SetLineWidth(2);

   entry_E_PFO=leg_E_PFO->AddEntry("hist_E_ph_PFO_Transition_075_087","Photon 0.75<|Cos#theta|<0.87","pe");
   entry_E_PFO->SetMarkerStyle(kOpenDiamond);
   entry_E_PFO->SetMarkerColor(kGreen+2);
   entry_E_PFO->SetLineColor(kGreen+2);
   entry_E_PFO->SetLineStyle(1);
   entry_E_PFO->SetLineWidth(2);
   entry_E_PFO=leg_E_PFO->AddEntry("hist_E_ph_PFO_Endcap_087_098","Photon 0.87<|Cos#theta|<0.98","pe");
   entry_E_PFO->SetMarkerStyle(kOpenCross);
   entry_E_PFO->SetMarkerColor(kCyan+1);
   entry_E_PFO->SetLineColor(kCyan+1);
   entry_E_PFO->SetLineStyle(1);
   entry_E_PFO->SetLineWidth(2);
   leg_E_PFO->Draw();

                    
   TCanvas *HITEnergyResolutionGraphCanvas = new TCanvas("HITEnergyResolutionGraphCanvas", "HITEnergyResolutionGraphCanvas",0,0,800,700);
   gStyle->SetOptStat(0);
   HITEnergyResolutionGraphCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
   HITEnergyResolutionGraphCanvas->SetFillColor(0);
   HITEnergyResolutionGraphCanvas->SetBorderMode(0);
   HITEnergyResolutionGraphCanvas->SetBorderSize(2);
   HITEnergyResolutionGraphCanvas->SetGridx();
   HITEnergyResolutionGraphCanvas->SetGridy();
   HITEnergyResolutionGraphCanvas->SetRightMargin(0.0172);
   HITEnergyResolutionGraphCanvas->SetTopMargin(0.0164);
   HITEnergyResolutionGraphCanvas->SetBottomMargin(0.125);
   HITEnergyResolutionGraphCanvas->SetFrameBorderMode(0);
   HITEnergyResolutionGraphCanvas->SetFrameBorderMode(0);



   TGraphErrors *gre_E_HIT = new TGraphErrors(11);
   gre_E_HIT->SetName("E_ph_HIT_innerBarrel_025_Graph");
   gre_E_HIT->SetTitle("");
   gre_E_HIT->SetFillColor(1);
   gre_E_HIT->SetLineColor(kBlack);
   gre_E_HIT->SetMarkerColor(kBlack);
   gre_E_HIT->SetMarkerStyle(kOpenCircle);

   gre_E_HIT->SetPoint(0,1,14.6221);
   gre_E_HIT->SetPointError(0,0,0.246466);
   gre_E_HIT->SetPoint(1,5,6.35512);
   gre_E_HIT->SetPointError(1,0,0.088714);
   gre_E_HIT->SetPoint(2,10,4.75765);
   gre_E_HIT->SetPointError(2,0,0.0333619);
   gre_E_HIT->SetPoint(3,15,3.85014);
   gre_E_HIT->SetPointError(3,0,0.0609838);
   gre_E_HIT->SetPoint(4,30,2.79257);
   gre_E_HIT->SetPointError(4,0,0.0409409);
   gre_E_HIT->SetPoint(5,50,2.2439);
   gre_E_HIT->SetPointError(5,0,0.0311863);
   gre_E_HIT->SetPoint(6,100,1.64473);
   gre_E_HIT->SetPointError(6,0,0.0242098);
   gre_E_HIT->SetPoint(7,200,1.25174);
   gre_E_HIT->SetPointError(7,0,0.0189752);
   gre_E_HIT->SetPoint(8,500,0.968162);
   gre_E_HIT->SetPointError(8,0,0.0146818);
   gre_E_HIT->SetPoint(9,1000,0.814196);
   gre_E_HIT->SetPointError(9,0,0.0112405);
   gre_E_HIT->SetPoint(10,1500,0.774489);
   gre_E_HIT->SetPointError(10,0,0.010919);
   
   TH1F *hist_E_ph_HIT_innerBarrel_025 = new TH1F("hist_E_ph_HIT_innerBarrel_025","",100,0,1649.9);
   hist_E_ph_HIT_innerBarrel_025->SetMinimum(0);
   hist_E_ph_HIT_innerBarrel_025->SetMaximum(6);
   hist_E_ph_HIT_innerBarrel_025->SetDirectory(0);
   hist_E_ph_HIT_innerBarrel_025->SetStats(0);

   hist_E_ph_HIT_innerBarrel_025->SetLineColor(kBlack);
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetLabelFont(42);
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetLabelSize(0.05);
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetTitleSize(0.05);
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetTitleFont(42);
   hist_E_ph_HIT_innerBarrel_025->GetXaxis()->SetTitleOffset(1.18);
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetTitle("#sigma(E_{ECAL}^{hits}+E_{HCAL}^{hits})/E_{true}_{#gamma} [%]");
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetLabelFont(42);
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetLabelSize(0.05);
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetTitleSize(0.05);
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetTitleOffset(1.20);
   hist_E_ph_HIT_innerBarrel_025->GetYaxis()->SetTitleFont(42);
   hist_E_ph_HIT_innerBarrel_025->GetZaxis()->SetLabelFont(42);
   hist_E_ph_HIT_innerBarrel_025->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_innerBarrel_025->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_innerBarrel_025->GetZaxis()->SetTitleFont(42);
   gre_E_HIT->SetHistogram(hist_E_ph_HIT_innerBarrel_025);
   
   gre_E_HIT->Draw("ape");
   
   gre_E_HIT = new TGraphErrors(11);
   gre_E_HIT->SetName("Graph_E_ph_HIT_mediumBarrel_025_050");
   gre_E_HIT->SetTitle("");
   gre_E_HIT->SetLineColor(kRed-7);

   gre_E_HIT->SetMarkerColor(kRed-7);
   gre_E_HIT->SetMarkerStyle(kOpenSquare);

   gre_E_HIT->SetPoint(0,1,14.7807);
   gre_E_HIT->SetPointError(0,0,0.243956);
   gre_E_HIT->SetPoint(1,5,6.91781);
   gre_E_HIT->SetPointError(1,0,0.105245);
   gre_E_HIT->SetPoint(2,10,4.8771);
   gre_E_HIT->SetPointError(2,0,0.0343484);
   gre_E_HIT->SetPoint(3,15,3.98902);
   gre_E_HIT->SetPointError(3,0,0.0606569);
   gre_E_HIT->SetPoint(4,30,2.83459);
   gre_E_HIT->SetPointError(4,0,0.0406483);
   gre_E_HIT->SetPoint(5,50,2.23702);
   gre_E_HIT->SetPointError(5,0,0.0332875);
   gre_E_HIT->SetPoint(6,100,1.65662);
   gre_E_HIT->SetPointError(6,0,0.0252521);
   gre_E_HIT->SetPoint(7,200,1.26838);
   gre_E_HIT->SetPointError(7,0,0.0191601);
   gre_E_HIT->SetPoint(8,500,0.957154);
   gre_E_HIT->SetPointError(8,0,0.0143418);
   gre_E_HIT->SetPoint(9,1000,0.819906);
   gre_E_HIT->SetPointError(9,0,0.0121039);
   gre_E_HIT->SetPoint(10,1500,0.772114);
   gre_E_HIT->SetPointError(10,0,0.0104655);

   TH1F *hist_E_ph_HIT_mediumBarrel_025_050 = new TH1F("hist_E_ph_HIT_mediumBarrel_025_050","",100,0,1649.9);
   hist_E_ph_HIT_mediumBarrel_025_050->SetMinimum(0);
   hist_E_ph_HIT_mediumBarrel_025_050->SetMaximum(19.79475);
   hist_E_ph_HIT_mediumBarrel_025_050->SetDirectory(0);
   hist_E_ph_HIT_mediumBarrel_025_050->SetStats(0);

   hist_E_ph_HIT_mediumBarrel_025_050->SetLineColor(kRed-7);
   hist_E_ph_HIT_mediumBarrel_025_050->GetXaxis()->SetLabelFont(42);
   hist_E_ph_HIT_mediumBarrel_025_050->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetXaxis()->SetTitleFont(42);
   hist_E_ph_HIT_mediumBarrel_025_050->GetYaxis()->SetLabelFont(42);
   hist_E_ph_HIT_mediumBarrel_025_050->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetYaxis()->SetTitleFont(42);
   hist_E_ph_HIT_mediumBarrel_025_050->GetZaxis()->SetLabelFont(42);
   hist_E_ph_HIT_mediumBarrel_025_050->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_mediumBarrel_025_050->GetZaxis()->SetTitleFont(42);
   gre_E_HIT->SetHistogram(hist_E_ph_HIT_mediumBarrel_025_050);
   gre_E_HIT->Draw("pe");

  
   gre_E_HIT = new TGraphErrors(11);
   gre_E_HIT->SetName("Graph_E_ph_HIT_outerBarrel_050_075");
   gre_E_HIT->SetTitle("");
   gre_E_HIT->SetLineColor(kBlue);

   gre_E_HIT->SetMarkerColor(kBlue);
   gre_E_HIT->SetMarkerStyle(kOpenTriangleUp);

   gre_E_HIT->SetPoint(0,1,15.3154);
   gre_E_HIT->SetPointError(0,0,0.229863);
   gre_E_HIT->SetPoint(1,5,7.21758);
   gre_E_HIT->SetPointError(1,0,0.112546);
   gre_E_HIT->SetPoint(2,10,4.8771);
   gre_E_HIT->SetPointError(2,0,0.0343484);
   gre_E_HIT->SetPoint(3,15,4.15576);
   gre_E_HIT->SetPointError(3,0,0.0582111);
   gre_E_HIT->SetPoint(4,30,2.8801);
   gre_E_HIT->SetPointError(4,0,0.0400708);
   gre_E_HIT->SetPoint(5,50,2.31738);
   gre_E_HIT->SetPointError(5,0,0.0348981);
   gre_E_HIT->SetPoint(6,100,1.69441);
   gre_E_HIT->SetPointError(6,0,0.0244134);
   gre_E_HIT->SetPoint(7,200,1.27814);
   gre_E_HIT->SetPointError(7,0,0.0185205);
   gre_E_HIT->SetPoint(8,500,0.898104);
   gre_E_HIT->SetPointError(8,0,0.0140954);
   gre_E_HIT->SetPoint(9,1000,0.748763);
   gre_E_HIT->SetPointError(9,0,0.0114007);
   gre_E_HIT->SetPoint(10,1500,0.687051);
   gre_E_HIT->SetPointError(10,0,0.00925405);
 
   TH1F *hist_E_ph_HIT_outerBarrel_050_075 = new TH1F("hist_E_ph_HIT_outerBarrel_050_075","",100,0,1649.9);
   hist_E_ph_HIT_outerBarrel_050_075->SetMinimum(0);
   hist_E_ph_HIT_outerBarrel_050_075->SetMaximum(19.79475);
   hist_E_ph_HIT_outerBarrel_050_075->SetDirectory(0);
   hist_E_ph_HIT_outerBarrel_050_075->SetStats(0);

   hist_E_ph_HIT_outerBarrel_050_075->SetLineColor(kBlue);
   hist_E_ph_HIT_outerBarrel_050_075->GetXaxis()->SetLabelFont(42);
   hist_E_ph_HIT_outerBarrel_050_075->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetXaxis()->SetTitleFont(42);
   hist_E_ph_HIT_outerBarrel_050_075->GetYaxis()->SetLabelFont(42);
   hist_E_ph_HIT_outerBarrel_050_075->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetYaxis()->SetTitleFont(42);
   hist_E_ph_HIT_outerBarrel_050_075->GetZaxis()->SetLabelFont(42);
   hist_E_ph_HIT_outerBarrel_050_075->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_outerBarrel_050_075->GetZaxis()->SetTitleFont(42);
   gre_E_HIT->SetHistogram(hist_E_ph_HIT_outerBarrel_050_075);
   gre_E_HIT->Draw("pe");

 
   gre_E_HIT = new TGraphErrors(11);
   gre_E_HIT->SetName("Graph_E_ph_HIT_Transition_075_087");
   gre_E_HIT->SetTitle("");
   gre_E_HIT->SetLineColor(kGreen+2);

   gre_E_HIT->SetMarkerColor(kGreen+2);
   gre_E_HIT->SetMarkerStyle(kOpenDiamond);
   gre_E_HIT->SetPoint(0,1,15.6261);
   gre_E_HIT->SetPointError(0,0,0.369859);
   gre_E_HIT->SetPoint(1,5,7.4008);
   gre_E_HIT->SetPointError(1,0,0.178709);
   gre_E_HIT->SetPoint(2,10,5.07078);
   gre_E_HIT->SetPointError(2,0,0.0505584);
   gre_E_HIT->SetPoint(3,15,4.18517);
   gre_E_HIT->SetPointError(3,0,0.100256);
   gre_E_HIT->SetPoint(4,30,3.13539);
   gre_E_HIT->SetPointError(4,0,0.0709226);
   gre_E_HIT->SetPoint(5,50,2.3726);
   gre_E_HIT->SetPointError(5,0,0.0495531);
   gre_E_HIT->SetPoint(6,100,1.79833);
   gre_E_HIT->SetPointError(6,0,0.0428804);
   gre_E_HIT->SetPoint(7,200,1.42661);
   gre_E_HIT->SetPointError(7,0,0.0440726);
   gre_E_HIT->SetPoint(8,500,0.70655);
   gre_E_HIT->SetPointError(8,0,0.0210691);
   gre_E_HIT->SetPoint(9,1000,0.70655);
   gre_E_HIT->SetPointError(9,0,0.0210691);
   gre_E_HIT->SetPoint(10,1500,0.602679);
   gre_E_HIT->SetPointError(10,0,0.019165);
 
   TH1F *hist_E_ph_HIT_Transition_075_087 = new TH1F("hist_E_ph_HIT_Transition_075_087","",100,0,1649.9);
   hist_E_ph_HIT_Transition_075_087->SetMinimum(0);
   hist_E_ph_HIT_Transition_075_087->SetMaximum(19.79475);
   hist_E_ph_HIT_Transition_075_087->SetDirectory(0);
   hist_E_ph_HIT_Transition_075_087->SetStats(0);

   hist_E_ph_HIT_Transition_075_087->SetLineColor(kGreen+2);
   hist_E_ph_HIT_Transition_075_087->GetXaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Transition_075_087->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetXaxis()->SetTitleFont(42);
   hist_E_ph_HIT_Transition_075_087->GetYaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Transition_075_087->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetYaxis()->SetTitleFont(42);
   hist_E_ph_HIT_Transition_075_087->GetZaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Transition_075_087->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Transition_075_087->GetZaxis()->SetTitleFont(42);
   gre_E_HIT->SetHistogram(hist_E_ph_HIT_Transition_075_087);
   gre_E_HIT->Draw("pe");

   gre_E_HIT = new TGraphErrors(11);
   gre_E_HIT->SetName("Graph_E_ph_HIT_Endcap_087_098");
   gre_E_HIT->SetTitle("");
   gre_E_HIT->SetLineColor(kCyan+1);
   gre_E_HIT->SetMarkerColor(kCyan+1);
   gre_E_HIT->SetMarkerStyle(kOpenCross);


   gre_E_HIT->SetPoint(0,1,15.1341);
   gre_E_HIT->SetPointError(0,0,0.383339);
   gre_E_HIT->SetPoint(1,5,6.723);
   gre_E_HIT->SetPointError(1,0,0.183003);
   gre_E_HIT->SetPoint(2,10,4.79806);
   gre_E_HIT->SetPointError(2,0,0.0512466);
   gre_E_HIT->SetPoint(3,15,3.79981);
   gre_E_HIT->SetPointError(3,0,0.0846845);
   gre_E_HIT->SetPoint(4,30,2.86315);
   gre_E_HIT->SetPointError(4,0,0.0786549);
   gre_E_HIT->SetPoint(5,50,2.2102);
   gre_E_HIT->SetPointError(5,0,0.0594837);
   gre_E_HIT->SetPoint(6,100,1.54802);
   gre_E_HIT->SetPointError(6,0,0.0354124);
   gre_E_HIT->SetPoint(7,200,1.13692);
   gre_E_HIT->SetPointError(7,0,0.0301534);
   gre_E_HIT->SetPoint(8,500,0.772743);
   gre_E_HIT->SetPointError(8,0,0.017799);
   gre_E_HIT->SetPoint(9,1000,0.602535);
   gre_E_HIT->SetPointError(9,0,0.0130302);
   gre_E_HIT->SetPoint(10,1500,0.545755);
   gre_E_HIT->SetPointError(10,0,0.013064);
 
   TH1F *hist_E_ph_HIT_Endcap_087_098 = new TH1F("hist_E_ph_HIT_Endcap_087_098","",100,0,1649.9);
   hist_E_ph_HIT_Endcap_087_098->SetMinimum(0);
   hist_E_ph_HIT_Endcap_087_098->SetMaximum(19.79475);
   hist_E_ph_HIT_Endcap_087_098->SetDirectory(0);
   hist_E_ph_HIT_Endcap_087_098->SetStats(0);

   hist_E_ph_HIT_Endcap_087_098->SetLineColor(kCyan+1);
   hist_E_ph_HIT_Endcap_087_098->GetXaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Endcap_087_098->GetXaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetXaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetXaxis()->SetTitleFont(42);
   hist_E_ph_HIT_Endcap_087_098->GetYaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Endcap_087_098->GetYaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetYaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetYaxis()->SetTitleFont(42);
   hist_E_ph_HIT_Endcap_087_098->GetZaxis()->SetLabelFont(42);
   hist_E_ph_HIT_Endcap_087_098->GetZaxis()->SetLabelSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetZaxis()->SetTitleSize(0.035);
   hist_E_ph_HIT_Endcap_087_098->GetZaxis()->SetTitleFont(42);
   gre_E_HIT->SetHistogram(hist_E_ph_HIT_Endcap_087_098);
   gre_E_HIT->Draw("pe");

   
   TLegend *leg_E_HIT = new TLegend(0.30,0.60,0.94,0.94,NULL,"brNDC");
   leg_E_HIT->SetBorderSize(1);
   leg_E_HIT->SetTextSize(0.027);
   leg_E_HIT->SetLineColor(1);
   leg_E_HIT->SetLineStyle(1);
   leg_E_HIT->SetLineWidth(1);
   leg_E_HIT->SetFillColor(0);
   leg_E_HIT->SetFillStyle(1001);
   TLegendEntry *entry_E_HIT=leg_E_HIT->AddEntry("hist_E_ph_HIT_innerBarrel_025","Photon |Cos#theta|<0.25","pe");
   entry_E_HIT->SetMarkerStyle(kOpenCircle);
   entry_E_HIT->SetMarkerColor(kBlack);
   entry_E_HIT->SetLineColor(kBlack);
   entry_E_HIT->SetLineStyle(1);
   entry_E_HIT->SetLineWidth(2);
   entry_E_HIT=leg_E_HIT->AddEntry("hist_E_ph_HIT_mediumBarrel_025_050","Photon 0.25<|Cos#theta|<0.50","pe");
   entry_E_HIT->SetMarkerStyle(kOpenSquare);
   entry_E_HIT->SetMarkerColor(kRed-7);
   entry_E_HIT->SetLineColor(kRed-7);
   entry_E_HIT->SetLineStyle(1);
   entry_E_HIT->SetLineWidth(2);

   entry_E_HIT=leg_E_HIT->AddEntry("hist_E_ph_HIT_outerBarrel_050_075","Photon 0.50<|Cos#theta|<0.75","pe");
   entry_E_HIT->SetMarkerStyle(kOpenTriangleUp);
   entry_E_HIT->SetMarkerColor(kBlue);
   entry_E_HIT->SetLineColor(kBlue);
   entry_E_HIT->SetLineStyle(1);
   entry_E_HIT->SetLineWidth(2);

   entry_E_HIT=leg_E_HIT->AddEntry("hist_E_ph_HIT_Transition_075_087","Photon 0.75<|Cos#theta|<0.87","pe");
   entry_E_HIT->SetMarkerStyle(kOpenDiamond);
   entry_E_HIT->SetMarkerColor(kGreen+2);
   entry_E_HIT->SetLineColor(kGreen+2);
   entry_E_HIT->SetLineStyle(1);
   entry_E_HIT->SetLineWidth(2);
   entry_E_HIT=leg_E_HIT->AddEntry("hist_E_ph_HIT_Endcap_087_098","Photon 0.87<|Cos#theta|<0.98","pe");
   entry_E_HIT->SetMarkerStyle(kOpenCross);
   entry_E_HIT->SetMarkerColor(kCyan+1);
   entry_E_HIT->SetLineColor(kCyan+1);
   entry_E_HIT->SetLineStyle(1);
   entry_E_HIT->SetLineWidth(2);
   leg_E_HIT->Draw();


   TCanvas *PFOPhiResolutionGraphCanvas = new TCanvas("PFOPhiResolutionGraphCanvas", "",0,0,800,700);
   gStyle->SetOptStat(0);
   PFOPhiResolutionGraphCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
   PFOPhiResolutionGraphCanvas->SetFillColor(0);
   PFOPhiResolutionGraphCanvas->SetBorderMode(0);
   PFOPhiResolutionGraphCanvas->SetBorderSize(2);
   PFOPhiResolutionGraphCanvas->SetGridx();
   PFOPhiResolutionGraphCanvas->SetGridy();
   PFOPhiResolutionGraphCanvas->SetRightMargin(0.0172);
   PFOPhiResolutionGraphCanvas->SetLeftMargin(0.20);
   PFOPhiResolutionGraphCanvas->SetTopMargin(0.018);
   PFOPhiResolutionGraphCanvas->SetBottomMargin(0.125);
   PFOPhiResolutionGraphCanvas->SetFrameBorderMode(0);
   PFOPhiResolutionGraphCanvas->SetFrameBorderMode(0);



   TGraphErrors *gre_phi_PFO = new TGraphErrors(11);
   gre_phi_PFO->SetName("phi_ph_PFO_innerBarrel_025_Graph");
   gre_phi_PFO->SetTitle("");
   gre_phi_PFO->SetFillColor(1);
   gre_phi_PFO->SetLineColor(kBlack);
   gre_phi_PFO->SetMarkerColor(kBlack);
   gre_phi_PFO->SetMarkerStyle(kOpenCircle);
   gre_phi_PFO->SetPoint(0,1,0.00124702);
   gre_phi_PFO->SetPointError(0,0,2.12122e-05);
   gre_phi_PFO->SetPoint(1,5,0.000683088);
   gre_phi_PFO->SetPointError(1,0,1.16664e-05);
   gre_phi_PFO->SetPoint(2,10,0.000534335);
   gre_phi_PFO->SetPointError(2,0,4.66041e-06);
   gre_phi_PFO->SetPoint(3,15,0.000469351);
   gre_phi_PFO->SetPointError(3,0,7.9935e-06);
   gre_phi_PFO->SetPoint(4,30,0.000342378);
   gre_phi_PFO->SetPointError(4,0,6.34721e-06);
   gre_phi_PFO->SetPoint(5,50,0.0002980);
   gre_phi_PFO->SetPointError(5,0,5.09314e-06);
   gre_phi_PFO->SetPoint(6,100,0.000213762);
   gre_phi_PFO->SetPointError(6,0,3.95409e-06);
   gre_phi_PFO->SetPoint(7,200,0.000163548);
   gre_phi_PFO->SetPointError(7,0,2.88893e-06);
   gre_phi_PFO->SetPoint(8,500,0.000121504);
   gre_phi_PFO->SetPointError(8,0,1.91637e-06);
   //seems off for all
   gre_phi_PFO->SetPoint(9,1000,9.6102e-05);
   gre_phi_PFO->SetPointError(9,0,1.97216e-06);
   //hift 0.0002658
   gre_phi_PFO->SetPoint(10,1500,9.06809e-05);
   gre_phi_PFO->SetPointError(10,0,1.9182e-06);
   
   TH1F *hist_phi_ph_PFO_innerBarrel_025 = new TH1F("hist_phi_ph_PFO_innerBarrel_025","",100,0,1649.9);
   hist_phi_ph_PFO_innerBarrel_025->SetMinimum(0);
   hist_phi_ph_PFO_innerBarrel_025->SetMaximum(0.0010);
   hist_phi_ph_PFO_innerBarrel_025->SetDirectory(0);
   hist_phi_ph_PFO_innerBarrel_025->SetStats(0);

   hist_phi_ph_PFO_innerBarrel_025->SetLineColor(kBlack);
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelSize(0.05);
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleSize(0.05);
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleOffset(1.18);
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetTitle("#sigma(phi_{RECO}^{PFO}) [rad]");
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelSize(0.045);
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleSize(0.05);
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleOffset(2.00);
   hist_phi_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleFont(42);
   gre_phi_PFO->SetHistogram(hist_phi_ph_PFO_innerBarrel_025);
   
   gre_phi_PFO->Draw("ape");
   
   gre_phi_PFO = new TGraphErrors(11);
   gre_phi_PFO->SetName("Graph_phi_ph_PFO_mediumBarrel_025_050");
   gre_phi_PFO->SetTitle("");
   gre_phi_PFO->SetLineColor(kRed-7);
   gre_phi_PFO->SetMarkerColor(kRed-7);
   gre_phi_PFO->SetMarkerStyle(kOpenSquare);
   gre_phi_PFO->SetPoint(0,1,0.00120014);
   gre_phi_PFO->SetPointError(0,0,2.5325e-05);
   gre_phi_PFO->SetPoint(1,5,0.000707273);
   gre_phi_PFO->SetPointError(1,0,1.14843e-05);
   gre_phi_PFO->SetPoint(2,10,0.000522809);
   gre_phi_PFO->SetPointError(2,0,4.03839e-06);
   gre_phi_PFO->SetPoint(3,15,0.000451546);
   gre_phi_PFO->SetPointError(3,0,7.51032e-06);
   gre_phi_PFO->SetPoint(4,30,0.000339495);
   gre_phi_PFO->SetPointError(4,0,6.08925e-06);
   gre_phi_PFO->SetPoint(5,50,0.000265168);
   gre_phi_PFO->SetPointError(5,0,4.89253e-06);
   gre_phi_PFO->SetPoint(6,100,0.000209597);
   gre_phi_PFO->SetPointError(6,0,3.53651e-06);
   gre_phi_PFO->SetPoint(7,200,0.000164445);
   gre_phi_PFO->SetPointError(7,0,2.60944e-06);
   gre_phi_PFO->SetPoint(8,500,0.000119863);
   gre_phi_PFO->SetPointError(8,0,2.0811e-06);
   //seems off for all
   gre_phi_PFO->SetPoint(9,1000,9.3253e-05);
   gre_phi_PFO->SetPointError(9,0,1.77511e-06);
   //hift 0.0002658
   gre_phi_PFO->SetPoint(10,1500,8.38946e-05);
   gre_phi_PFO->SetPointError(10,0,1.64855e-06);


 

   TH1F *hist_phi_ph_PFO_mediumBarrel_025_050 = new TH1F("hist_phi_ph_PFO_mediumBarrel_025_050","hist_phi_ph_PFO_mediumBarrel_025_050",100,0,1649.9);
   hist_phi_ph_PFO_mediumBarrel_025_050->SetMinimum(0);
   hist_phi_ph_PFO_mediumBarrel_025_050->SetMaximum(19.79475);
   hist_phi_ph_PFO_mediumBarrel_025_050->SetDirectory(0);
   hist_phi_ph_PFO_mediumBarrel_025_050->SetStats(0);

   hist_phi_ph_PFO_mediumBarrel_025_050->SetLineColor(kRed-7);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleFont(42);
   gre_phi_PFO->SetHistogram(hist_phi_ph_PFO_mediumBarrel_025_050);
   gre_phi_PFO->Draw("pe");

  
   gre_phi_PFO = new TGraphErrors(11);
   gre_phi_PFO->SetName("Graph_phi_ph_PFO_outerBarrel_050_075");
   gre_phi_PFO->SetTitle("");
   gre_phi_PFO->SetLineColor(kBlue);

   gre_phi_PFO->SetMarkerColor(kBlue);
   gre_phi_PFO->SetMarkerStyle(kOpenTriangleUp);
 
   gre_phi_PFO->SetPoint(0,1,0.00117924);
   gre_phi_PFO->SetPointError(0,0,2.30306e-05);
   gre_phi_PFO->SetPoint(1,5,0.000652212);
   gre_phi_PFO->SetPointError(1,0,1.0828e-05);
   gre_phi_PFO->SetPoint(2,10,0.000488857);
   gre_phi_PFO->SetPointError(2,0,5.60632e-06);
   gre_phi_PFO->SetPoint(3,15,0.000438967);
   gre_phi_PFO->SetPointError(3,0,6.59672e-06);
   gre_phi_PFO->SetPoint(4,30,0.000323471);
   gre_phi_PFO->SetPointError(4,0,5.04381e-06);
   gre_phi_PFO->SetPoint(5,50,0.00026323);
   gre_phi_PFO->SetPointError(5,0,4.29467e-06);
   gre_phi_PFO->SetPoint(6,100,0.000197111);
   gre_phi_PFO->SetPointError(6,0,3.00498e-06);
   gre_phi_PFO->SetPoint(7,200,0.000155413);
   gre_phi_PFO->SetPointError(7,0,2.57391e-06);
   gre_phi_PFO->SetPoint(8,500,0.000112065);
   gre_phi_PFO->SetPointError(8,0,1.76583e-06);
   //seems off for all
   gre_phi_PFO->SetPoint(9,1000,9.32723e-05);
   gre_phi_PFO->SetPointError(9,0,1.38331e-06);
   //hift 0.0002658
   gre_phi_PFO->SetPoint(10,1500,8.70863e-05);
   gre_phi_PFO->SetPointError(10,0,1.37668e-06);

 
   TH1F *hist_phi_ph_PFO_outerBarrel_050_075 = new TH1F("hist_phi_ph_PFO_outerBarrel_050_075","",100,0,1649.9);
   hist_phi_ph_PFO_outerBarrel_050_075->SetMinimum(0);
   hist_phi_ph_PFO_outerBarrel_050_075->SetMaximum(19.79475);
   hist_phi_ph_PFO_outerBarrel_050_075->SetDirectory(0);
   hist_phi_ph_PFO_outerBarrel_050_075->SetStats(0);

   hist_phi_ph_PFO_outerBarrel_050_075->SetLineColor(kBlue);
   hist_phi_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleFont(42);
   gre_phi_PFO->SetHistogram(hist_phi_ph_PFO_outerBarrel_050_075);
   gre_phi_PFO->Draw("pe");

 
   gre_phi_PFO = new TGraphErrors(11);
   gre_phi_PFO->SetName("Graph_phi_ph_PFO_Transition_075_087");
   gre_phi_PFO->SetTitle("");
   gre_phi_PFO->SetLineColor(kGreen+2);
   gre_phi_PFO->SetMarkerColor(kGreen+2);
   gre_phi_PFO->SetMarkerStyle(kOpenDiamond);
   gre_phi_PFO->SetPoint(0,1,0.00128288);
   gre_phi_PFO->SetPointError(0,0,3.75724e-05);
   gre_phi_PFO->SetPoint(1,5,0.000753169);
   gre_phi_PFO->SetPointError(1,0,2.11694e-05);
   gre_phi_PFO->SetPoint(2,10,0.000536608);
   gre_phi_PFO->SetPointError(2,0,2.83893e-05);
   gre_phi_PFO->SetPoint(3,15,0.000433308);
   gre_phi_PFO->SetPointError(3,0,1.17656e-05);
   gre_phi_PFO->SetPoint(4,30,0.00035683);
   gre_phi_PFO->SetPointError(4,0,9.70757e-06);
   gre_phi_PFO->SetPoint(5,50,0.000275847);
   gre_phi_PFO->SetPointError(5,0,6.80619e-06);
   gre_phi_PFO->SetPoint(6,100,0.000228979);
   gre_phi_PFO->SetPointError(6,0,6.31546e-06);
   gre_phi_PFO->SetPoint(7,200,0.000183998);
   gre_phi_PFO->SetPointError(7,0,4.93067e-06);
   gre_phi_PFO->SetPoint(8,500,0.000124664);
   gre_phi_PFO->SetPointError(8,0,2.75498e-06);
   gre_phi_PFO->SetPoint(9,1000,0.000124664);
   gre_phi_PFO->SetPointError(9,0,2.75498e-06);
   gre_phi_PFO->SetPoint(10,1500,0.00011826);
   gre_phi_PFO->SetPointError(10,0,6.21883e-06);
 
   TH1F *hist_phi_ph_PFO_Transition_075_087 = new TH1F("hist_phi_ph_PFO_Transition_075_087","",100,0,1649.9);
   hist_phi_ph_PFO_Transition_075_087->SetMinimum(0);
   hist_phi_ph_PFO_Transition_075_087->SetMaximum(19.79475);
   hist_phi_ph_PFO_Transition_075_087->SetDirectory(0);
   hist_phi_ph_PFO_Transition_075_087->SetStats(0);
   hist_phi_ph_PFO_Transition_075_087->SetLineColor(kGreen+2);
   hist_phi_ph_PFO_Transition_075_087->GetXaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Transition_075_087->GetXaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetXaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetXaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_Transition_075_087->GetYaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Transition_075_087->GetYaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetYaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetYaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_Transition_075_087->GetZaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Transition_075_087->GetZaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetZaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Transition_075_087->GetZaxis()->SetTitleFont(42);
   gre_phi_PFO->SetHistogram(hist_phi_ph_PFO_Transition_075_087);
   gre_phi_PFO->Draw("pe");

   gre_phi_PFO = new TGraphErrors(11);
   gre_phi_PFO->SetName("Graph_phi_ph_PFO_Endcap_087_098");
   gre_phi_PFO->SetTitle("");
   gre_phi_PFO->SetLineColor(kCyan+1);

   gre_phi_PFO->SetMarkerColor(kCyan+1);
   gre_phi_PFO->SetMarkerStyle(kOpenCross);
   gre_phi_PFO->SetPoint(0,1,0.000480557);
   gre_phi_PFO->SetPointError(0,0,4.84346e-05);
   gre_phi_PFO->SetPoint(1,5,0.000346778);
   gre_phi_PFO->SetPointError(1,0,2.00065e-05);
   gre_phi_PFO->SetPoint(2,10,0.000266806);
   gre_phi_PFO->SetPointError(2,0,4.72632e-06);
   gre_phi_PFO->SetPoint(3,15,0.00024861);
   gre_phi_PFO->SetPointError(3,0,1.12344e-05);
   gre_phi_PFO->SetPoint(4,30,0.000157067);
   gre_phi_PFO->SetPointError(4,0,5.91438e-06);
   gre_phi_PFO->SetPoint(5,50,0.000137119);
   gre_phi_PFO->SetPointError(5,0,4.55531e-06);
   gre_phi_PFO->SetPoint(6,100,0.000112185);
   gre_phi_PFO->SetPointError(6,0,3.78071e-06);
   gre_phi_PFO->SetPoint(7,200,8.86089e-05);
   gre_phi_PFO->SetPointError(7,0,3.38888e-06);
   gre_phi_PFO->SetPoint(8,500,6.06863e-05);//maybe 8.06863e-05 ??
   gre_phi_PFO->SetPointError(8,0,2.38686e-06);
   gre_phi_PFO->SetPoint(9,1000,5.09302e-05);
   gre_phi_PFO->SetPointError(9,0,1.88043e-06);
   gre_phi_PFO->SetPoint(10,1500,4.22612e-05);
   gre_phi_PFO->SetPointError(10,0,1.55717e-06);

 
   TH1F *hist_phi_ph_PFO_Endcap_087_098 = new TH1F("hist_phi_ph_PFO_Endcap_087_098","",100,0,1649.9);
   hist_phi_ph_PFO_Endcap_087_098->SetMinimum(0);
   hist_phi_ph_PFO_Endcap_087_098->SetMaximum(19.79475);
   hist_phi_ph_PFO_Endcap_087_098->SetDirectory(0);
   hist_phi_ph_PFO_Endcap_087_098->SetStats(0);

   hist_phi_ph_PFO_Endcap_087_098->SetLineColor(kCyan+1);
   hist_phi_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleFont(42);
   hist_phi_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelFont(42);
   hist_phi_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleSize(0.035);
   hist_phi_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleFont(42);
   gre_phi_PFO->SetHistogram(hist_phi_ph_PFO_Endcap_087_098);
   gre_phi_PFO->Draw("pe");

   
   TLegend *leg_phi_PFO = new TLegend(0.30,0.60,0.94,0.94,NULL,"brNDC");
   leg_phi_PFO->SetBorderSize(1);
   leg_phi_PFO->SetTextSize(0.027);
   leg_phi_PFO->SetLineColor(1);
   leg_phi_PFO->SetLineStyle(1);
   leg_phi_PFO->SetLineWidth(1);
   leg_phi_PFO->SetFillColor(0);
   leg_phi_PFO->SetFillStyle(1001);
   TLegendEntry *entry_phi_PFO=leg_phi_PFO->AddEntry("hist_phi_ph_PFO_innerBarrel_025","Photon |Cos#theta|<0.25","pe");
   entry_phi_PFO->SetMarkerStyle(kOpenCircle);
   entry_phi_PFO->SetMarkerColor(kBlack);
   entry_phi_PFO->SetLineColor(kBlack);
   entry_phi_PFO->SetLineStyle(1);
   entry_phi_PFO->SetLineWidth(2);
   entry_phi_PFO=leg_phi_PFO->AddEntry("hist_phi_ph_PFO_mediumBarrel_025_050","Photon 0.25<|Cos#theta|<0.50","pe");
   entry_phi_PFO->SetMarkerStyle(kOpenSquare);
   entry_phi_PFO->SetMarkerColor(kRed-7);
   entry_phi_PFO->SetLineColor(kRed-7);
   entry_phi_PFO->SetLineStyle(1);
   entry_phi_PFO->SetLineWidth(2);

   entry_phi_PFO=leg_phi_PFO->AddEntry("hist_phi_ph_PFO_outerBarrel_050_075","Photon 0.50<|Cos#theta|<0.75","pe");
   entry_phi_PFO->SetMarkerStyle(kOpenTriangleUp);
   entry_phi_PFO->SetMarkerColor(kBlue);
   entry_phi_PFO->SetLineColor(kBlue);
   entry_phi_PFO->SetLineStyle(1);
   entry_phi_PFO->SetLineWidth(2);

   entry_phi_PFO=leg_phi_PFO->AddEntry("hist_phi_ph_PFO_Transition_075_087","Photon 0.75<|Cos#theta|<0.87","pe");
   entry_phi_PFO->SetMarkerStyle(kOpenDiamond);
   entry_phi_PFO->SetMarkerColor(kGreen+2);
   entry_phi_PFO->SetLineColor(kGreen+2);
   entry_phi_PFO->SetLineStyle(1);
   entry_phi_PFO->SetLineWidth(2);
   entry_phi_PFO=leg_phi_PFO->AddEntry("hist_phi_ph_PFO_Endcap_087_098","Photon 0.87<|Cos#theta|<0.98","pe");
   entry_phi_PFO->SetMarkerStyle(kOpenCross);
   entry_phi_PFO->SetMarkerColor(kCyan+1);
   entry_phi_PFO->SetLineColor(kCyan+1);
   entry_phi_PFO->SetLineStyle(1);
   entry_phi_PFO->SetLineWidth(2);
   leg_phi_PFO->Draw();



  TCanvas *PFOThetaResolutionGraphCanvas = new TCanvas("PFOThetaResolutionGraphCanvas", "",0,0,800,700);
   gStyle->SetOptStat(0);
   PFOThetaResolutionGraphCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
   PFOThetaResolutionGraphCanvas->SetFillColor(0);
   PFOThetaResolutionGraphCanvas->SetBorderMode(0);
   PFOThetaResolutionGraphCanvas->SetBorderSize(2);
   PFOThetaResolutionGraphCanvas->SetGridx();
   PFOThetaResolutionGraphCanvas->SetGridy();
   PFOThetaResolutionGraphCanvas->SetRightMargin(0.0172);
   PFOThetaResolutionGraphCanvas->SetLeftMargin(0.20);
   PFOThetaResolutionGraphCanvas->SetTopMargin(0.018);
   PFOThetaResolutionGraphCanvas->SetBottomMargin(0.125);
   PFOThetaResolutionGraphCanvas->SetFrameBorderMode(0);
   PFOThetaResolutionGraphCanvas->SetFrameBorderMode(0);



   TGraphErrors *gre_theta_PFO = new TGraphErrors(11);
   gre_theta_PFO->SetName("theta_ph_PFO_innerBarrel_025_Graph");
   gre_theta_PFO->SetTitle("");
   gre_theta_PFO->SetFillColor(1);
   gre_theta_PFO->SetLineColor(kBlack);
   gre_theta_PFO->SetMarkerColor(kBlack);
   gre_theta_PFO->SetMarkerStyle(kOpenCircle);

   gre_theta_PFO->SetPoint(0,1,0.0010880981);
   gre_theta_PFO->SetPointError(0,0,4.80693e-05);
   gre_theta_PFO->SetPoint(1,5,0.000577756);
   gre_theta_PFO->SetPointError(1,0,1.49089e-05);
   gre_theta_PFO->SetPoint(2,10,0.000462226);
   gre_theta_PFO->SetPointError(2,0,5.1543e-06);
   gre_theta_PFO->SetPoint(3,15,0.000408671);
   gre_theta_PFO->SetPointError(3,0,9.5075e-06);
   gre_theta_PFO->SetPoint(4,30,0.000305285);
   gre_theta_PFO->SetPointError(4,0,7.38444e-06);
   gre_theta_PFO->SetPoint(5,50,0.00024874);
   gre_theta_PFO->SetPointError(5,0,5.59047e-06);
   gre_theta_PFO->SetPoint(6,100,0.000199758);
   gre_theta_PFO->SetPointError(6,0,4.55222e-06);
   gre_theta_PFO->SetPoint(7,200,0.000153064);
   gre_theta_PFO->SetPointError(7,0,3.65218e-06);
   gre_theta_PFO->SetPoint(8,500,0.000119975);
   gre_theta_PFO->SetPointError(8,0,3.50617e-06);
   gre_theta_PFO->SetPoint(9,1000,0.882152);
   gre_theta_PFO->SetPointError(9,0,0.022244);
   gre_theta_PFO->SetPoint(10,1500,9.13229e-05);
   gre_theta_PFO->SetPointError(10,0,2.57551e-06);

   
   TH1F *hist_theta_ph_PFO_innerBarrel_025 = new TH1F("hist_theta_ph_PFO_innerBarrel_025","",100,0,1649.9);
   hist_theta_ph_PFO_innerBarrel_025->SetMinimum(0);
   hist_theta_ph_PFO_innerBarrel_025->SetMaximum(0.0010);
   hist_theta_ph_PFO_innerBarrel_025->SetDirectory(0);
   hist_theta_ph_PFO_innerBarrel_025->SetStats(0);

   hist_theta_ph_PFO_innerBarrel_025->SetLineColor(kBlack);
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetLabelSize(0.05);
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleSize(0.05);
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_innerBarrel_025->GetXaxis()->SetTitleOffset(1.18);
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetTitle("#sigma(theta_{RECO}^{PFO}) [rad]");
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetLabelSize(0.045);
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleSize(0.05);
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleOffset(2.00);
   hist_theta_ph_PFO_innerBarrel_025->GetYaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_innerBarrel_025->GetZaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_innerBarrel_025->GetZaxis()->SetTitleFont(42);
   gre_theta_PFO->SetHistogram(hist_theta_ph_PFO_innerBarrel_025);
   
   gre_theta_PFO->Draw("ape");
   
   gre_theta_PFO = new TGraphErrors(11);
   gre_theta_PFO->SetName("Graph_theta_ph_PFO_mediumBarrel_025_050");
   gre_theta_PFO->SetTitle("");
   gre_theta_PFO->SetLineColor(kRed-7);
   gre_theta_PFO->SetMarkerColor(kRed-7);
   gre_theta_PFO->SetMarkerStyle(kOpenSquare);
   gre_theta_PFO->SetPoint(0,1,0.00092057);
   gre_theta_PFO->SetPointError(0,0,2.4713e-05);
   gre_theta_PFO->SetPoint(1,5,0.000595767);
   gre_theta_PFO->SetPointError(1,0,1.58724e-05);
   gre_theta_PFO->SetPoint(2,10,0.000424997);
   gre_theta_PFO->SetPointError(2,0,3.98357e-06);
   gre_theta_PFO->SetPoint(3,15,0.000372653);
   gre_theta_PFO->SetPointError(3,0,8.79647e-06);
   gre_theta_PFO->SetPoint(4,30,0.000274631);
   gre_theta_PFO->SetPointError(4,0,6.65575e-06);
   gre_theta_PFO->SetPoint(5,50,0.000230315);
   gre_theta_PFO->SetPointError(5,0,5.46118e-06);
   gre_theta_PFO->SetPoint(6,100,0.000168057);
   gre_theta_PFO->SetPointError(6,0,4.10172e-06);
   gre_theta_PFO->SetPoint(7,200,0.000123637);
   gre_theta_PFO->SetPointError(7,0,3.15432e-06);
   gre_theta_PFO->SetPoint(8,500,9.85406e-05);
   gre_theta_PFO->SetPointError(8,0,2.33344e-06);
   gre_theta_PFO->SetPoint(9,1000, 7.86159e-05);
   gre_theta_PFO->SetPointError(9,0,1.87519e-06);
   gre_theta_PFO->SetPoint(10,1500,6.69543e-05);
   gre_theta_PFO->SetPointError(10,0,1.80124e-06);


 

   TH1F *hist_theta_ph_PFO_mediumBarrel_025_050 = new TH1F("hist_theta_ph_PFO_mediumBarrel_025_050","hist_theta_ph_PFO_mediumBarrel_025_050",100,0,1649.9);
   hist_theta_ph_PFO_mediumBarrel_025_050->SetMinimum(0);
   hist_theta_ph_PFO_mediumBarrel_025_050->SetMaximum(19.79475);
   hist_theta_ph_PFO_mediumBarrel_025_050->SetDirectory(0);
   hist_theta_ph_PFO_mediumBarrel_025_050->SetStats(0);

   hist_theta_ph_PFO_mediumBarrel_025_050->SetLineColor(kRed-7);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetXaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetYaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_mediumBarrel_025_050->GetZaxis()->SetTitleFont(42);
   gre_theta_PFO->SetHistogram(hist_theta_ph_PFO_mediumBarrel_025_050);
   gre_theta_PFO->Draw("pe");

  
   gre_theta_PFO = new TGraphErrors(11);
   gre_theta_PFO->SetName("Graph_theta_ph_PFO_outerBarrel_050_075");
   gre_theta_PFO->SetTitle("");
   gre_theta_PFO->SetLineColor(kBlue);

   gre_theta_PFO->SetMarkerColor(kBlue);
   gre_theta_PFO->SetMarkerStyle(kOpenTriangleUp);
   gre_theta_PFO->SetPoint(0,1,0.000861644);
   gre_theta_PFO->SetPointError(0,0,2.37433e-05);
   gre_theta_PFO->SetPoint(1,5,0.000497132);
   gre_theta_PFO->SetPointError(1,0,1.28861e-05);
   gre_theta_PFO->SetPoint(2,10,0.000348073);
   gre_theta_PFO->SetPointError(2,0,4.29101e-06);
   gre_theta_PFO->SetPoint(3,15,0.000300222);
   gre_theta_PFO->SetPointError(3,0,7.88762e-06);
   gre_theta_PFO->SetPoint(4,30,0.00022784);
   gre_theta_PFO->SetPointError(4,0,4.965e-06);
   gre_theta_PFO->SetPoint(5,50,0.000189182);
   gre_theta_PFO->SetPointError(5,0,4.64235e-06);
   gre_theta_PFO->SetPoint(6,100,0.000131703);
   gre_theta_PFO->SetPointError(6,0,2.97601e-06);
   gre_theta_PFO->SetPoint(7,200,0.000109247);
   gre_theta_PFO->SetPointError(7,0,2.46978e-06);
   gre_theta_PFO->SetPoint(8,500,7.29701e-05);
   gre_theta_PFO->SetPointError(8,0,1.70915e-06);
   gre_theta_PFO->SetPoint(9,1000,5.77957e-05);
   gre_theta_PFO->SetPointError(9,0,1.54542e-06);
   gre_theta_PFO->SetPoint(10,1500,4.99557e-05);
   gre_theta_PFO->SetPointError(10,0,1.22063e-06);


 
   TH1F *hist_theta_ph_PFO_outerBarrel_050_075 = new TH1F("hist_theta_ph_PFO_outerBarrel_050_075","",100,0,1649.9);
   hist_theta_ph_PFO_outerBarrel_050_075->SetMinimum(0);
   hist_theta_ph_PFO_outerBarrel_050_075->SetMaximum(19.79475);
   hist_theta_ph_PFO_outerBarrel_050_075->SetDirectory(0);
   hist_theta_ph_PFO_outerBarrel_050_075->SetStats(0);

   hist_theta_ph_PFO_outerBarrel_050_075->SetLineColor(kBlue);
   hist_theta_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_outerBarrel_050_075->GetXaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetXaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_outerBarrel_050_075->GetYaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetYaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_outerBarrel_050_075->GetZaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_outerBarrel_050_075->GetZaxis()->SetTitleFont(42);
   gre_theta_PFO->SetHistogram(hist_theta_ph_PFO_outerBarrel_050_075);
   gre_theta_PFO->Draw("pe");


   gre_theta_PFO = new TGraphErrors(11);
   gre_theta_PFO->SetName("Graph_theta_ph_PFO_Endcap_087_098");
   gre_theta_PFO->SetTitle("");
   gre_theta_PFO->SetLineColor(kCyan+1);

   gre_theta_PFO->SetMarkerColor(kCyan+1);
   gre_theta_PFO->SetMarkerStyle(kOpenCross);
   
   gre_theta_PFO->SetPoint(0,1,0.000480557);
   gre_theta_PFO->SetPointError(0,0,4.84346e-05);
   gre_theta_PFO->SetPoint(1,5,0.000346778);
   gre_theta_PFO->SetPointError(1,0,2.00065e-05);
   gre_theta_PFO->SetPoint(2,10,0.000266806);
   gre_theta_PFO->SetPointError(2,0,4.72632e-06);
   gre_theta_PFO->SetPoint(3,15,0.00024861);
   gre_theta_PFO->SetPointError(3,0,1.12344e-05);
   gre_theta_PFO->SetPoint(4,30,0.000157067);
   gre_theta_PFO->SetPointError(4,0,5.91438e-06);
   gre_theta_PFO->SetPoint(5,50,0.000137119);
   gre_theta_PFO->SetPointError(5,0,4.55531e-06);
   gre_theta_PFO->SetPoint(6,100,0.000112185);
   gre_theta_PFO->SetPointError(6,0,3.78071e-06);
   gre_theta_PFO->SetPoint(7,200,8.86089e-05);
   gre_theta_PFO->SetPointError(7,0,3.38888e-06);
   gre_theta_PFO->SetPoint(8,500,6.06863e-05);
   gre_theta_PFO->SetPointError(8,0,2.38686e-06);
   gre_theta_PFO->SetPoint(9,1000,5.09302e-05);
   gre_theta_PFO->SetPointError(9,0,1.88043e-06);
   gre_theta_PFO->SetPoint(10,1500,4.22612e-05);
   gre_theta_PFO->SetPointError(10,0,1.55717e-06);

 
   TH1F *hist_theta_ph_PFO_Endcap_087_098 = new TH1F("hist_theta_ph_PFO_Endcap_087_098","",100,0,1649.9);
   hist_theta_ph_PFO_Endcap_087_098->SetMinimum(0);
   hist_theta_ph_PFO_Endcap_087_098->SetMaximum(19.79475);
   hist_theta_ph_PFO_Endcap_087_098->SetDirectory(0);
   hist_theta_ph_PFO_Endcap_087_098->SetStats(0);

   hist_theta_ph_PFO_Endcap_087_098->SetLineColor(kCyan+1);
   hist_theta_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_Endcap_087_098->GetXaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetXaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_Endcap_087_098->GetYaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetYaxis()->SetTitleFont(42);
   hist_theta_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelFont(42);
   hist_theta_ph_PFO_Endcap_087_098->GetZaxis()->SetLabelSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleSize(0.035);
   hist_theta_ph_PFO_Endcap_087_098->GetZaxis()->SetTitleFont(42);
   gre_theta_PFO->SetHistogram(hist_theta_ph_PFO_Endcap_087_098);
   gre_theta_PFO->Draw("pe");

   
   TLegend *leg_theta_PFO = new TLegend(0.30,0.60,0.94,0.94,NULL,"brNDC");
   leg_theta_PFO->SetBorderSize(1);
   leg_theta_PFO->SetTextSize(0.027);
   leg_theta_PFO->SetLineColor(1);
   leg_theta_PFO->SetLineStyle(1);
   leg_theta_PFO->SetLineWidth(1);
   leg_theta_PFO->SetFillColor(0);
   leg_theta_PFO->SetFillStyle(1001);
   TLegendEntry *entry_theta_PFO=leg_theta_PFO->AddEntry("hist_theta_ph_PFO_innerBarrel_025","Photon |Cos#theta|<0.25","pe");
   entry_theta_PFO->SetMarkerStyle(kOpenCircle);
   entry_theta_PFO->SetMarkerColor(kBlack);
   entry_theta_PFO->SetLineColor(kBlack);
   entry_theta_PFO->SetLineStyle(1);
   entry_theta_PFO->SetLineWidth(2);
   entry_theta_PFO=leg_theta_PFO->AddEntry("hist_theta_ph_PFO_mediumBarrel_025_050","Photon 0.25<|Cos#theta|<0.50","pe");
   entry_theta_PFO->SetMarkerStyle(kOpenSquare);
   entry_theta_PFO->SetMarkerColor(kRed-7);
   entry_theta_PFO->SetLineColor(kRed-7);
   entry_theta_PFO->SetLineStyle(1);
   entry_theta_PFO->SetLineWidth(2);

   entry_theta_PFO=leg_theta_PFO->AddEntry("hist_theta_ph_PFO_outerBarrel_050_075","Photon 0.50<|Cos#theta|<0.75","pe");
   entry_theta_PFO->SetMarkerStyle(kOpenTriangleUp);
   entry_theta_PFO->SetMarkerColor(kBlue);
   entry_theta_PFO->SetLineColor(kBlue);
   entry_theta_PFO->SetLineStyle(1);
   entry_theta_PFO->SetLineWidth(2);
   entry_theta_PFO=leg_theta_PFO->AddEntry("hist_theta_ph_PFO_Endcap_087_098","Photon 0.87<|Cos#theta|<0.98","pe");
   entry_theta_PFO->SetMarkerStyle(kOpenCross);
   entry_theta_PFO->SetMarkerColor(kCyan+1);
   entry_theta_PFO->SetLineColor(kCyan+1);
   entry_theta_PFO->SetLineStyle(1);
   entry_theta_PFO->SetLineWidth(2);
   leg_theta_PFO->Draw();


   //pi 1 GeV TT, res values

   //rticle E signal0.998476/0.000236487/0.0224027/0.000167221
   //particle E tot1.02504/0.00110684/0.0980166/0.000782657

   //gamma 1 GeV
   //particle E signal0.973544/0.00189963/0.180445/0.00134324
   //particle E tot0.984541/0.0017681/0.168629/0.00125023
   //gamma 10 GeV
   //particle E signal1.00031/0.000283235/0.0582047/0.000200278
   //particle E tot1.00619/0.000256762/0.0541244/0.000181558
   //pi 10 GeV
   //particle E signal1.00004/3.52232e-05/0.0151956/2.49066e-05
   //particle E tot1.01704/0.000183417/0.0792035/0.000129695
   //n 10 GeV
   //particle E signal0.843485/0.000427909/0.167323/0.000302578
   //particle E tot0.859796/0.000396934/0.160954/0.000280675

}
