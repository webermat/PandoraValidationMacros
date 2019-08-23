void resolutionSummaryKaons_K0L_2_1500()
{
//=========Macro generated from canvas: resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary/resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary
//=========  (Thu Mar 22 19:47:34 2018) by ROOT version6.08/00
   TCanvas *resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary = new TCanvas("resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary", "resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary",2,51,800,700);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->Range(-345.0824,2.276952,1684.814,18.386);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetFillColor(0);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetBorderMode(0);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetBorderSize(2);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetGridx();
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetGridy();
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetTickx(1);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetTicky(1);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetLeftMargin(0.17);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetRightMargin(0.0172);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetTopMargin(0.055);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetBottomMargin(0.138);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetFrameLineWidth(2);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetFrameBorderMode(0);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetFrameLineWidth(2);
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetFrameBorderMode(0);
   
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fx1001[18] = {
   2,
   5,
   10,
   20,
   30,
   40,
   50,
   60,
   75,
   90,
   100,
   150,
   200,
   250,
   400,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fy1001[18] = {
   36.4344,
   25.6253,
   17.2049,
   11.1274,
   8.79933,
   7.66443,
   7.03673,
   6.70435,
   6.35791,
   6.19407,
   6.00422,
   5.76138,
   5.75447,
   5.72028,
   5.59529,
   5.53228,
   5.3759,
   5.48218};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fex1001[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fey1001[18] = {
   0.739078,
   0.237696,
   0.0868783,
   0.054239,
   0.0401531,
   0.0367902,
   0.0291822,
   0.0301525,
   0.0255388,
   0.0253676,
   0.0245256,
   0.0240602,
   0.0237341,
   0.0240146,
   0.0232333,
   0.0226888,
   0.0258426,
   0.0910431};
   TGraphErrors *gre = new TGraphErrors(18,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fx1001,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fy1001,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fex1001,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_fey1001);
   gre->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001 = new TH1F("Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001","",100,0,1649.9);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetMinimum(4.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetMaximum(17.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetDirectory(0);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetStats(0);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetLineWidth(2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->SetMarkerSize(1.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetXaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetYaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetLabelSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetTitleSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_CLIC_o3_v141001);
   
   gre->Draw("ape");
   
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fx1002[18] = {
   2,
   5,
   10,
   20,
   30,
   40,
   50,
   60,
   75,
   90,
   100,
   150,
   200,
   250,
   400,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fy1002[18] = {
   36.4833,
   25.8934,
   16.4662,
   11.089,
   8.6318,
   7.87732,
   7.26023,
   6.73712,
   6.56918,
   6.11101,
   6.19851,
   5.65581,
   5.56125,
   5.6564,
   5.56311,
   5.4082,
   5.15672,
   5.13651};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fex1002[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fey1002[18] = {
   1.7139,
   0.559949,
   0.180735,
   0.131805,
   0.100182,
   0.105985,
   0.075536,
   0.0799938,
   0.0706449,
   0.0666521,
   0.0745047,
   0.0623192,
   0.069142,
   0.0697636,
   0.0739029,
   0.0683279,
   0.0787054,
   0.355858};
   gre = new TGraphErrors(18,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fx1002,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fy1002,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fex1002,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80_fey1002);
   gre->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_65_to_0_80");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff6666");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002 = new TH1F("Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002","",100,0,1649.9);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetMinimum(4.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetMaximum(17.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetDirectory(0);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetStats(0);

   ci = TColor::GetColor("#ff6666");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetLineColor(ci);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetLineWidth(2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->SetMarkerSize(1.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetXaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetYaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetLabelSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetTitleSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_K0L_gcc62_E_reco_vs_E_true_0_65_to_0_80_CLIC_o3_v141002);
   
   gre->Draw("pe");
   
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fx1003[18] = {
   2,
   5,
   10,
   20,
   30,
   40,
   50,
   60,
   75,
   90,
   100,
   150,
   200,
   250,
   400,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fy1003[18] = {
   37.6859,
   24.5551,
   15.7351,
   10.4049,
   8.75853,
   7.51838,
   7.3408,
   6.8912,
   6.5579,
   6.37767,
   6.12034,
   5.6872,
   5.6244,
   5.57162,
   5.56826,
   5.5777,
   5.42587,
   5.13964};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fex1003[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fey1003[18] = {
   1.79119,
   0.453611,
   0.167803,
   0.114396,
   0.0836679,
   0.0756588,
   0.0643899,
   0.0660442,
   0.0562063,
   0.0569723,
   0.0598237,
   0.0537094,
   0.053084,
   0.0609708,
   0.0558178,
   0.0534595,
   0.0634342,
   0.196675};
   gre = new TGraphErrors(18,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fx1003,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fy1003,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fex1003,CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94_fey1003);
   gre->SetName("CLIC_o3_v14_res_graph_K0L_gcc62_rel_E_reco_vs_E_true_0_80_to_0_94");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003 = new TH1F("Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003","",100,0,1649.9);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetMinimum(4.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetMaximum(17.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetDirectory(0);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetLineColor(ci);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetLineWidth(2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->SetMarkerSize(1.5);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetTitle("E_{true}^{K0L} [GeV]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetXaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{K0L}) [%]");
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetNdivisions(506);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetLabelSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetTitleSize(0.05);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetYaxis()->SetTitleFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetLabelFont(42);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetLabelSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetTitleSize(0.035);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_K0L_gcc62_E_reco_vs_E_true_0_80_to_0_94_CLIC_o3_v141003);
   
   gre->Draw("pe");
   
   TLegend *leg = new TLegend(0.6,0.6,0.94,0.91,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetTextSize(0.027);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","K0L","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("NULL","|cos#theta_{K0L}^{true}|<0.65","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","0.65<|cos#theta_{K0L}^{true}|<0.80","pe");

   ci = TColor::GetColor("#ff6666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff6666");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","0.80<|cos#theta_{K0L}^{true}|<0.94","pe");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->Modified();
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->cd();
   resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary->SetSelected(resolution_gcc62_GraphCanvas_res_rel_K0L_CLIC_n_K0L_FullSummary);
}
