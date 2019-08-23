void resolutionValuesPhoton_ph1_1500()
{
//=========Macro generated from canvas: resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary/resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary
//=========  (Thu Mar 22 19:32:16 2018) by ROOT version6.08/00
   TCanvas *resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary = new TCanvas("resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary", "resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary",2,51,800,700);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->Range(-3.450824,-2.407063,16.84814,18.65861);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetFillColor(0);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetBorderMode(0);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetBorderSize(2);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetGridx();
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetGridy();
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetTickx(1);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetTicky(1);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetLeftMargin(0.17);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetRightMargin(0.0172);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetTopMargin(0.055);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetBottomMargin(0.138);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetFrameLineWidth(2);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetFrameBorderMode(0);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetFrameLineWidth(2);
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetFrameBorderMode(0);
   
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fx1001[11] = {
   1,
   5,
   10,
   15,
   30,
   50,
   100,
   200,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fy1001[11] = {
   14.9688,
   7.00205,
   5.02102,
   4.25148,
   2.99076,
   2.42809,
   1.78763,
   1.33247,
   1.03128,
   0.894994,
   0.860321};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fex1001[11] = {
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
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fey1001[11] = {
   0.151503,
   0.063763,
   0.0428576,
   0.0352131,
   0.0252701,
   0.0221697,
   0.01637,
   0.0125256,
   0.0105665,
   0.010169,
   0.00923972};
   TGraphErrors *gre = new TGraphErrors(11,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fx1001,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fy1001,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fex1001,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fey1001);
   gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001 = new TH1F("Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001","",100,0,1649.9);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetMinimum(0.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetMaximum(17.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetDirectory(0);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetStats(0);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetLineWidth(2);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->SetMarkerSize(1.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetRange(1,1);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetXaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetYaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetLabelSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetTitleSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v141001);
   
   gre->Draw("ape");
   
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fx1002[11] = {
   1,
   5,
   10,
   15,
   30,
   50,
   100,
   200,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fy1002[11] = {
   14.9688,
   7.00205,
   5.02102,
   4.25148,
   2.99076,
   2.42809,
   1.78763,
   1.33247,
   1.03128,
   0.894994,
   0.860321};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fex1002[11] = {
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
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fey1002[11] = {
   0.151503,
   0.063763,
   0.0428576,
   0.0352131,
   0.0252701,
   0.0221697,
   0.01637,
   0.0125256,
   0.0105665,
   0.010169,
   0.00923972};
   gre = new TGraphErrors(11,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fx1002,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fy1002,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fex1002,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_fey1002);
   gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002 = new TH1F("Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002","",100,0,1649.9);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetMinimum(0.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetMaximum(17.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetDirectory(0);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetStats(0);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetLineWidth(2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->SetMarkerSize(1.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetRange(1,1);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetNdivisions(506);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetTitleOffset(1.18);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetXaxis()->SetTitleFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetNdivisions(506);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetYaxis()->SetTitleFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_h_ph_E_reco_vs_E_true_0_78_CLIC_o3_v1410011002);
   
   gre->Draw("pe");
   
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fx1003[11] = {
   1,
   5,
   10,
   15,
   30,
   50,
   100,
   200,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fy1003[11] = {
   16.5011,
   7.50338,
   5.92026,
   4.55544,
   3.10068,
   2.49905,
   1.81596,
   1.27697,
   1.19124,
   0.797591,
   0.761724};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fex1003[11] = {
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
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fey1003[11] = {
   0.865444,
   0.421356,
   0.277318,
   0.192466,
   0.127407,
   0.115118,
   0.0904222,
   0.0588278,
   0.109471,
   0.0538643,
   0.0436831};
   gre = new TGraphErrors(11,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fx1003,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fy1003,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fex1003,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fey1003);
   gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff6666");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003 = new TH1F("Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003","",100,0,1649.9);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetMinimum(0.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetMaximum(17.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetDirectory(0);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetStats(0);

   ci = TColor::GetColor("#ff6666");
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetLineColor(ci);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetLineWidth(2);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->SetMarkerSize(1.5);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetXaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetYaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetLabelSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetTitleSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v141003);
   
   gre->Draw("pe");
   
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fx1004[11] = {
   1,
   5,
   10,
   15,
   30,
   50,
   100,
   200,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fy1004[11] = {
   16.5011,
   7.50338,
   5.92026,
   4.55544,
   3.10068,
   2.49905,
   1.81596,
   1.27697,
   1.19124,
   0.797591,
   0.761724};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fex1004[11] = {
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
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fey1004[11] = {
   0.865444,
   0.421356,
   0.277318,
   0.192466,
   0.127407,
   0.115118,
   0.0904222,
   0.0588278,
   0.109471,
   0.0538643,
   0.0436831};
   gre = new TGraphErrors(11,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fx1004,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fy1004,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fex1004,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83_fey1004);
   gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_78_to_0_83");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#ff6666");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004 = new TH1F("Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004","",100,0,1649.9);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetMinimum(0.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetMaximum(17.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetDirectory(0);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetStats(0);

   ci = TColor::GetColor("#ff6666");
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetLineColor(ci);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetLineWidth(2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->SetMarkerSize(1.5);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetNdivisions(506);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetTitleOffset(1.18);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetXaxis()->SetTitleFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetNdivisions(506);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetYaxis()->SetTitleFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetLabelFont(42);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_h_ph_E_reco_vs_E_true_0_78_to_0_83_CLIC_o3_v1410031004);
   
   gre->Draw("pe");
   
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fx1005[11] = {
   1,
   5,
   10,
   15,
   30,
   50,
   100,
   200,
   500,
   1000,
   1500};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fy1005[11] = {
   14.9116,
   6.97093,
   4.74473,
   3.88792,
   2.78804,
   2.2429,
   1.59499,
   1.15164,
   0.735529,
   0.551053,
   0.451682};
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fex1005[11] = {
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
   Double_t CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fey1005[11] = {
   0.429003,
   0.189273,
   0.131706,
   0.0957513,
   0.0697395,
   0.0687657,
   0.040368,
   0.0276968,
   0.0180497,
   0.016265,
   0.0127654};
   gre = new TGraphErrors(11,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fx1005,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fy1005,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fex1005,CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94_fey1005);
   gre->SetName("CLIC_o3_v14_res_graph_ph_rel_E_reco_vs_E_true_0_83_to_0_94");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005 = new TH1F("Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005","",100,0,1649.9);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetMinimum(0.5);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetMaximum(17.5);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetDirectory(0);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetLineColor(ci);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetLineWidth(2);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->SetMarkerSize(1.5);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetTitle("E_{true}^{#gamma} [GeV]");
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetTitleOffset(1.18);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetXaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetTitle("#sigma(E_{PFO}^{reco}/E_{true}^{#gamma}) [%]");
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetNdivisions(506);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetLabelSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetTitleSize(0.05);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetYaxis()->SetTitleFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetLabelFont(42);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetLabelOffset(0.015);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetLabelSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetTitleSize(0.035);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetTitleOffset(1.2);
   Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_h_ph_E_reco_vs_E_true_0_83_to_0_94_CLIC_o3_v141005);
   
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
   TLegendEntry *entry=leg->AddEntry("NULL","ph","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("NULL","|cos#theta_{#gamma}^{true}|<0.78","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","0.78<|cos#theta_{#gamma}^{true}|<0.83","pe");

   ci = TColor::GetColor("#ff6666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff6666");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","0.80<|cos#theta_{#gamma}^{true}|<0.94","pe");

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
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->Modified();
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->cd();
   resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary->SetSelected(resolutionGraphCanvas_res_rel_ph_CLIC_n_K0L_FullSummary);
}
