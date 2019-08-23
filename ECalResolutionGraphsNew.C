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


void plotECALSummary(){


 CLICdpStyle();


//=========Macro generated from canvas: resolutionGraphCanvas/
//=========  (Tue Aug  2 10:131:49 2016) by ROOT version5.34/34
   TCanvas *resolutionGraphCanvas = new TCanvas("resolutionGraphCanvas", "",0,0,800,700);
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
   ci = TColor::GetColor(kBlack);
   gre->SetLineColor(ci);

   ci = TColor::GetColor(kBlack);
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(kFullCircle);
   //gre->SetMarkerStyle(kOpenCircle);
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
   
   TH1F *Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1 = new TH1F("Graph_CLIC_o2_v04_CUR_22.85X0_resolutionGraph1","",100,0,1649.9);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetMinimum(0);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetMaximum(6);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetDirectory(0);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetStats(0);

   ci = TColor::GetColor(kBlack);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->SetLineColor(ci);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitle("E^{true}_{#gamma} [GeV]");
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetLabelFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetLabelSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleSize(0.05);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleFont(42);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetXaxis()->SetTitleOffset(1.18);
   Graph_CLIC_o2_v04_CUR_22_85X0_resolutionGraph1->GetYaxis()->SetTitle("#sigma(E_{ECal}+E_{HCal})/E [%]");
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
   gre->Draw("ape");
   
   gre = new TGraphErrors(8);
   gre->SetName("CLICdet20_10_CDR_23.35X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineColor(kRed-7);
   gre->SetMarkerColor(kRed-7);
   gre->SetMarkerStyle(kFullSquare);
   //gre->SetMarkerStyle(kOpenSquare);
   gre->SetPoint(0,1,16.34374);
   gre->SetPointError(0,0,0.1385529);
   gre->SetPoint(1,10,5.362804);
   gre->SetPointError(1,0,0.04106499);
   gre->SetPoint(2,50,2.660309);
   gre->SetPointError(2,0,0.01969936);
   gre->SetPoint(3,100,2.01346);
   gre->SetPointError(3,0,0.01553595);
   gre->SetPoint(4,200,1.461181);
   gre->SetPointError(4,0,0.01068301);
   gre->SetPoint(5,500,0.9808201);
   gre->SetPointError(5,0,0.007458539);
   gre->SetPoint(6,1000,0.7111181);
   gre->SetPointError(6,0,0.005086895);
   gre->SetPoint(7,1500,0.6003635);
   gre->SetPointError(7,0,0.004395324);
   
   TH1F *Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2 = new TH1F("Graph_CLICdet20_10_CDR_23.35X0_resolutionGraph2","",100,0,1649.9);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->SetMinimum(0);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->SetMaximum(18.07093);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->SetDirectory(0);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->SetStats(0);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->SetLineColor(kRed-7);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetXaxis()->SetLabelFont(42);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetXaxis()->SetLabelSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetXaxis()->SetTitleSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetXaxis()->SetTitleFont(42);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetYaxis()->SetLabelFont(42);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetYaxis()->SetLabelSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetYaxis()->SetTitleSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetYaxis()->SetTitleFont(42);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetZaxis()->SetLabelFont(42);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetZaxis()->SetLabelSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetZaxis()->SetTitleSize(0.035);
   Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_CLICdet20_10_CDR_23_35X0_resolutionGraph2);
   gre->Draw("pe");
   
   
   gre = new TGraphErrors(8);
   gre->SetName("CLICdet30_23.216X0_resolutionGraph");
   gre->SetTitle("");
   gre->SetFillColor(1);

   ci = TColor::GetColor(kBlue);
   gre->SetLineColor(kBlue);

   gre->SetMarkerColor(kBlue);
   //gre->SetMarkerStyle(kOpenTriangleUp);
   gre->SetMarkerStyle(kFullTriangleUp);
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

   ci = TColor::GetColor(kBlue);
   Graph_CLICdet30_23_216X0_resolutionGraph4->SetLineColor(kBlue);
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
   gre->SetFillColor(kGreen+2);

   ci = TColor::GetColor(kGreen+2);
   gre->SetLineColor(ci);
   gre->SetMarkerColor(kGreen+2);
   gre->SetMarkerStyle(kFullDiamond);
   //gre->SetMarkerStyle(kOpenDiamond);
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

   ci = TColor::GetColor(kGreen+2);
   Graph_CLICdet40_o3_V05_resolutionGraph6->SetLineColor(kGreen+2);
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
   gre->Draw("pe");


   
   TLegend *leg = new TLegend(0.60,0.60,0.94,0.94,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextSize(0.027);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   //TLegendEntry *entry=leg->AddEntry("Graph_CLIC_o2_v04_CUR_22.85X0_resolutionGraph1","CLICdet_17_8","pe");
   TLegendEntry *entry=leg->AddEntry("Graph_CLIC_o2_v04_CUR_22.85X0_resolutionGraph1","17+8 layers","pe");
   entry->SetLineColor(kBlack);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(kBlack);
   //entry->SetMarkerStyle(kOpenCircle);
   entry->SetMarkerStyle(kFullCircle);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   //entry=leg->AddEntry("Graph_CLICdet20_10_CDR_23.35X0_resolutionGraph2","CLICdet_20_10","pe");
   entry=leg->AddEntry("Graph_CLICdet20_10_CDR_23.35X0_resolutionGraph2","20+10 layers","pe");
   entry->SetLineColor(kRed-7);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(kRed-7);
   //entry->SetMarkerStyle(kOpenSquare);
   entry->SetMarkerStyle(kFullSquare);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);

   //entry=leg->AddEntry("Graph_CLICdet30_23_216X0_resolutionGraph4","CLICdet_30","pe");
   entry=leg->AddEntry("Graph_CLICdet30_23_216X0_resolutionGraph4","30 layers","pe");
   entry->SetLineColor(kBlue);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(kBlue);
   //entry->SetMarkerStyle(kOpenTriangleUp + 1000*10);
   entry->SetMarkerStyle(kFullTriangleUp + 1000*10);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   //entry=leg->AddEntry("Graph_CLICdet40_o3_V05X0_resolutionGraph6","CLICdet_40_b","pe");
   entry=leg->AddEntry("Graph_CLICdet40_o3_V05X0_resolutionGraph6","40 layers","pe");
   entry->SetLineColor(kGreen+2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(kGreen+2);
   //entry->SetMarkerStyle(kOpenDiamond);
   entry->SetMarkerStyle(kFullDiamond);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);

   leg->Draw();
}
