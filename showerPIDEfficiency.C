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
#include "TVector3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"

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
//}


float sigmaE_rel(float E_true,float cosTheta){
  float a_stoch=1.6;
  float a_const=0.99;
  float sigma=-1;
  if(fabs(cosTheta)<0.78){
    a_stoch=16.;
    a_const=0.99;
  }else if (fabs(cosTheta)<0.83){
    a_stoch=18.;
    a_const=2.2e-5;
  }else{
    a_stoch=15.1;
    a_const=0.57;
  }
   sigma=sqrt(pow(a_stoch/sqrt(E_true),2)+pow(a_const,2))*0.01;
   return sigma;
}

void fillRESTreeHistograms(TFile* file_tree, std::vector<TH1F*> hist, int particle_type_true, int particle_type_reco){
  TTree* tree = (TTree*)file_tree->Get("showerData");

  vector<float> *trueE=0;
  vector<float> *truePx=0;
  vector<float> *truePy=0;
  vector<float> *truePz=0;
  vector<float> *truePhi=0;
  vector<float> *trueCosTheta=0;
  vector<int> *truePDGID=0;
  vector<int> *trueGenStatus=0;

  vector<float> *recoE=0;
  vector<float> *recoPx=0;
  vector<float> *recoPy=0;
  vector<float> *recoPz=0;
  vector<float> *recoPhi=0;
  vector<float> *recoCosTheta=0;
  vector<int> *recoPDGID=0;

  //tree->SetBranchAddress("eventNumber",&eventNumber);
  tree->SetBranchAddress("true_Energy",&trueE);
  tree->SetBranchAddress("true_Px",&truePx);
  tree->SetBranchAddress("true_Py",&truePy);
  tree->SetBranchAddress("true_Pz",&truePz);
  tree->SetBranchAddress("true_Phi",&truePhi);
  tree->SetBranchAddress("true_GenStatus",&trueGenStatus);
  tree->SetBranchAddress("true_CosTheta",&trueCosTheta);
  tree->SetBranchAddress("true_PDGID",&truePDGID);

  tree->SetBranchAddress("reco_Energy",&recoE);
  tree->SetBranchAddress("reco_Px",&recoPx);
  tree->SetBranchAddress("reco_Py",&recoPy);
  tree->SetBranchAddress("reco_Pz",&recoPz);
  tree->SetBranchAddress("reco_Phi",&recoPhi);
  tree->SetBranchAddress("reco_CosTheta",&recoCosTheta);
  tree->SetBranchAddress("reco_PDGID",&recoPDGID);

  for(unsigned int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(i%3000==0){
      std::cout<<"i_entry "<<i<<std::endl;
    }

    //resolution histograms for detector note
    TLorentzVector TTrueRES(0,0,0,0);
    for(unsigned int i=0;i<trueE->size();i++){
      if((*trueGenStatus)[i]==1 && abs((*truePDGID)[i])==particle_type_true){
	TTrueRES.SetPxPyPzE((*truePx)[i],(*truePy)[i],(*truePz)[i],(*trueE)[i]);
      }
    }
    int index_recoRES=-1;
    double AngleMinRES=2*TMath::Pi();
    for(unsigned int i=0;i<recoE->size();i++){
      if(abs((*recoPDGID)[i])==particle_type_reco){
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE((*recoPx)[i],(*recoPy)[i],(*recoPz)[i],(*recoE)[i]);
	if(TTrueRES.Angle(temp.Vect())<AngleMinRES){
	  AngleMinRES=TTrueRES.Angle(temp.Vect());
	  index_recoRES=i;
	}
      }
    }
    //index now known
    TLorentzVector TRecoRES(0,0,0,0);
    if(index_recoRES!=-1){
      TRecoRES.SetPxPyPzE((*recoPx)[index_recoRES],(*recoPy)[index_recoRES],(*recoPz)[index_recoRES],(*recoE)[index_recoRES]);
    }
    if( (TTrueRES.Angle(TRecoRES.Vect())*TMath::RadToDeg())<2.0){
      if(particle_type_true==22){//photon resolution
	if(fabs(cos(TTrueRES.Theta()))<0.78){
	  hist[0]->Fill(TRecoRES.E()/TTrueRES.E());
	}else if(fabs(cos(TTrueRES.Theta()))<0.83){
	  hist[1]->Fill(TRecoRES.E()/TTrueRES.E());
	}else if(fabs(cos(TTrueRES.Theta()))<0.94){
	  hist[2]->Fill(TRecoRES.E()/TTrueRES.E());
	}
      }else if(particle_type_true==130){//Kaon resolution
	if(fabs(cos(TTrueRES.Theta()))<0.65){
	  hist[0]->Fill(TRecoRES.E()/TTrueRES.E());
	}else if(fabs(cos(TTrueRES.Theta()))<0.80){
	  hist[1]->Fill(TRecoRES.E()/TTrueRES.E());
	}else if(fabs(cos(TTrueRES.Theta()))<0.94){
	  hist[2]->Fill(TRecoRES.E()/TTrueRES.E());
	}
      }
    }
  }

  delete tree;
}

void fillTreeHistograms(TFile* file_tree, std::vector<TH1F*> histo_vector_oneE, int particle_type, std::vector<TEfficiency*> tEff_vector_allE,std::vector<TEfficiency*> tEff_vector_oneE, std::vector<TH1F*> histo_vector_allE){
  TTree* tree = (TTree*)file_tree->Get("showerData");
  
  
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
  vector<float> *true_conv_e1_E=0;
  vector<float> *true_conv_e1_Px=0;
  vector<float> *true_conv_e1_Py=0;
  vector<float> *true_conv_e1_Pz=0;
  vector<float> *true_conv_e2_E=0;
  vector<float> *true_conv_e2_Px=0;
  vector<float> *true_conv_e2_Py=0;
  vector<float> *true_conv_e2_Pz=0;

  vector<float> *true_conv_Vtx_x=0;
  vector<float> *true_conv_Vtx_y=0;
  vector<float> *true_conv_Vtx_z=0;
  
  vector<float> *recoE=0;
  vector<float> *recoPx=0;
  vector<float> *recoPy=0;
  vector<float> *recoPz=0;
  vector<float> *recoPhi=0;
  vector<float> *recoCosTheta=0;
  vector<int> *recoCharge=0;
  vector<int> *recoPDGID=0;
  
  //tree->SetBranchAddress("eventNumber",&eventNumber);
  tree->SetBranchAddress("true_Energy",&trueE);
  tree->SetBranchAddress("true_Px",&truePx);
  tree->SetBranchAddress("true_Py",&truePy);
  tree->SetBranchAddress("true_Pz",&truePz);
  tree->SetBranchAddress("true_Phi",&truePhi);
  tree->SetBranchAddress("true_GenStatus",&trueGenStatus);
  tree->SetBranchAddress("true_decayTrackerCalo",&trueDecayTrackerCalo);
  tree->SetBranchAddress("true_motherDecayTrackerCalo",&trueMotherDecayTrackerCalo);
  tree->SetBranchAddress("true_CosTheta",&trueCosTheta);
  //tree->SetBranchAddress("true_Theta",&trueTheta);
  tree->SetBranchAddress("true_PDGID",&truePDGID);
  tree->SetBranchAddress("true_numDaughters",&trueNumDaughters);
  if(particle_type==22){
    tree->SetBranchAddress("trueConvE1E",&true_conv_e1_E);
    tree->SetBranchAddress("trueConvE1Px",&true_conv_e1_Px);
    tree->SetBranchAddress("trueConvE1Py",&true_conv_e1_Py);
    tree->SetBranchAddress("trueConvE1Pz",&true_conv_e1_Pz);
    tree->SetBranchAddress("trueConvE2E",&true_conv_e2_E);
    tree->SetBranchAddress("trueConvE2Px",&true_conv_e2_Px);
    tree->SetBranchAddress("trueConvE2Py",&true_conv_e2_Py);
    tree->SetBranchAddress("trueConvE2Pz",&true_conv_e2_Pz);
    tree->SetBranchAddress("trueConvVtxx",&true_conv_Vtx_x);
    tree->SetBranchAddress("trueConvVtxy",&true_conv_Vtx_y);
    tree->SetBranchAddress("trueConvVtxz",&true_conv_Vtx_z);
  }
  
  tree->SetBranchAddress("reco_Energy",&recoE);
  tree->SetBranchAddress("reco_Px",&recoPx);
  tree->SetBranchAddress("reco_Py",&recoPy);
  tree->SetBranchAddress("reco_Pz",&recoPz);
  tree->SetBranchAddress("reco_Phi",&recoPhi);
  tree->SetBranchAddress("reco_CosTheta",&recoCosTheta);
  tree->SetBranchAddress("reco_PDGID",&recoPDGID);
  tree->SetBranchAddress("reco_Charge",&recoCharge);
  
  float angle_match_max=1.;//in degrees
  float photon_energy_term=0.16;//2 parameter fit for photons from 10-50 or 4-200 is between 15-15.6 %, for transition region 17.5 --> take 16 %

  float delta_energy_correction=0;
  float delta_energy_correction_barrel=0;
  float delta_energy_correction_TR=0;
  float delta_energy_correction_endcap=0;

  std::vector<float>count_all_events_perEBin;
  std::vector<float>count_all_ConvEvents_perEBin;

  std::vector<float>count_all_events_perThetaBin;
  std::vector<float>count_all_ConvEvents_perThetaBin;

  if(particle_type==22){
    
    for(int i=0;i<histo_vector_allE[0]->GetNbinsX();i++){
      count_all_events_perEBin.push_back(0);
      count_all_ConvEvents_perEBin.push_back(0);
    }
    for(int i=0;i<histo_vector_oneE[18]->GetNbinsX();i++){
      count_all_events_perThetaBin.push_back(0);
      count_all_ConvEvents_perThetaBin.push_back(0);
    }
  
    float true_energy=-10;
    TH1F* h_hist_energy_temp_barrel = new TH1F("hist_temp_energy_barrel","",1000,0.5,1.5);
    h_hist_energy_temp_barrel->Sumw2();
    TH1F* h_hist_energy_temp_endcap = new TH1F("hist_temp_energy_endcap","",1000,0.5,1.5);
    h_hist_energy_temp_endcap->Sumw2();
    TH1F* h_hist_energy_temp_TR = new TH1F("hist_temp_energy_TR","",1000,0.5,1.5);
    h_hist_energy_temp_TR->Sumw2();
    
    
    //if we check photons, the energy requirement is more tricky. We want that the photon energy is not far away regarding tails of 
    for(unsigned int i=0;i<tree->GetEntries();i++){
      tree->GetEntry(i);
      if(i%3000==0){
	std::cout<<"i_entry "<<i<<std::endl;
      }
      bool is_conversion=false;
      if(true_conv_e1_E->size()>0){
	is_conversion=true;
      }
      if(!is_conversion){
	TLorentzVector TTrue(0,0,0,0);
	for(unsigned int i=0;i<trueE->size();i++){
	  if((*trueGenStatus)[i]==1 && abs((*truePDGID)[i])==particle_type){
	    TTrue.SetPxPyPzE((*truePx)[i],(*truePy)[i],(*truePz)[i],(*trueE)[i]);
	    true_energy=(*trueE)[i];
	  }
	}
	int index_reco=-1;
	double AngleMin=2*TMath::Pi();
	for(unsigned int i=0;i<recoE->size();i++){
	  if(abs((*recoPDGID)[i])==particle_type){
	    TLorentzVector temp(0,0,0,0);
	    temp.SetPxPyPzE((*recoPx)[i],(*recoPy)[i],(*recoPz)[i],(*recoE)[i]);
	    if(TTrue.Angle(temp.Vect())<AngleMin){
	      AngleMin=TTrue.Angle(temp.Vect());
	      index_reco=i;
	    }
	  }
	}
	//index now known
	TLorentzVector TReco(0,0,0,0);
	if(index_reco!=-1){
	  TReco.SetPxPyPzE((*recoPx)[index_reco],(*recoPy)[index_reco],(*recoPz)[index_reco],(*recoE)[index_reco]);
	  if( (TTrue.Angle(TReco.Vect())*TMath::RadToDeg())<1.0){
	    if(fabs(TTrue.CosTheta())<0.78){
	      h_hist_energy_temp_barrel->Fill(TReco.Energy()/TTrue.Energy());
	    }else if (fabs(TTrue.CosTheta())<0.83){
	      h_hist_energy_temp_TR->Fill(TReco.Energy()/TTrue.Energy());
	    }else{
	      h_hist_energy_temp_endcap->Fill(TReco.Energy()/TTrue.Energy());
	    }
	  }
	}
      }
    }
    //now the histogram is filled with all entries, only for non converted photons, and for angular matched photons
    delta_energy_correction_barrel=(1.-h_hist_energy_temp_barrel->GetMean())*true_energy;
    delta_energy_correction_TR=(1.-h_hist_energy_temp_TR->GetMean())*true_energy;
    delta_energy_correction_endcap=(1.-h_hist_energy_temp_endcap->GetMean())*true_energy;
    
    std::cout<<"mean correction for energy response barrel/TR/endcap "<<delta_energy_correction_barrel<<"/"<<h_hist_energy_temp_barrel->GetMean()<<"/"<< delta_energy_correction_TR<<"/"<<h_hist_energy_temp_TR->GetMean()<<"/"<<delta_energy_correction_endcap<<"/"<<h_hist_energy_temp_endcap->GetMean()<<std::endl;
    delete h_hist_energy_temp_barrel;
    delete h_hist_energy_temp_TR;
    delete h_hist_energy_temp_endcap;
  }
  

  for(unsigned int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(i%3000==0){
      std::cout<<"i_entry "<<i<<std::endl;
    }

    bool is_conversion=false;
    if(particle_type==22){
      if(true_conv_e1_E->size()>0){
	//std::cout<<"photon is converted"<<std::endl;
	is_conversion=true;
      }else{    
	//std::cout<<"photon is not converted"<<std::endl;
      }
    }
    TLorentzVector TTrue(0,0,0,0);
    unsigned int ntrue=0;
    for(unsigned int i=0;i<trueE->size();i++){
      if((*trueGenStatus)[i]==1 && abs((*truePDGID)[i])==particle_type){
	TTrue.SetPxPyPzE((*truePx)[i],(*truePy)[i],(*truePz)[i],(*trueE)[i]);
	ntrue+=1;
      }
    }
    if(particle_type==22){
      for(int i=1;i<=histo_vector_allE[0]->GetNbinsX();i++){
	if(TTrue.E()>=histo_vector_allE[0]->GetBinLowEdge(i) && TTrue.E()<histo_vector_allE[0]->GetBinLowEdge(i+1)){
	  count_all_events_perEBin[i-1]+=1;
	  if(is_conversion){
	    count_all_ConvEvents_perEBin[i-1]+=1;
	  }
	} 
      }
      for(int i=1;i<=histo_vector_oneE[18]->GetNbinsX();i++){
	if((TTrue.Theta()*TMath::RadToDeg())>=histo_vector_oneE[18]->GetBinLowEdge(i) && (TTrue.Theta()*TMath::RadToDeg())<histo_vector_oneE[18]->GetBinLowEdge(i+1)){
	  count_all_events_perThetaBin[i-1]+=1;
	  if(is_conversion){
	    count_all_ConvEvents_perThetaBin[i-1]+=1;
	  } 
	}
      }
    }

    if(is_conversion){
      //we consider so far only particle gones, i.e. 1 particle present in the envet, don't need to loop over conversion vector
      TLorentzVector TTrueConv1(0,0,0,0);
      TTrueConv1.SetPxPyPzE((*true_conv_e1_Px)[0],(*true_conv_e1_Py)[0],(*true_conv_e1_Pz)[0],(*true_conv_e1_E)[0]);
      TLorentzVector TTrueConv2(0,0,0,0);
      TTrueConv2.SetPxPyPzE((*true_conv_e2_Px)[0],(*true_conv_e2_Py)[0],(*true_conv_e2_Pz)[0],(*true_conv_e2_E)[0]);
      if(fabs(TTrue.CosTheta())<0.78){
	histo_vector_oneE[2]->Fill((*true_conv_Vtx_z)[0]);
	histo_vector_oneE[4]->Fill(sqrt((*true_conv_Vtx_x)[0]*(*true_conv_Vtx_x)[0]+(*true_conv_Vtx_y)[0]*(*true_conv_Vtx_y)[0]));
	histo_vector_oneE[6]->Fill(TTrueConv1.Angle(TTrueConv2.Vect())*TMath::RadToDeg());
	if(TTrueConv1.E()>TTrueConv2.E()){
	  histo_vector_oneE[8]->Fill(TTrueConv1.E()/TTrue.E());
	  histo_vector_oneE[10]->Fill(TTrueConv2.E()/TTrue.E());
	}else{
	  histo_vector_oneE[8]->Fill(TTrueConv2.E()/TTrue.E());
	  histo_vector_oneE[10]->Fill(TTrueConv1.E()/TTrue.E());
	}
	if(TTrueConv1.Angle(TTrue.Vect())>TTrueConv2.Angle(TTrue.Vect())){
	  histo_vector_oneE[12]->Fill(TTrueConv1.Angle(TTrue.Vect())*TMath::RadToDeg());
	}else{
	  histo_vector_oneE[12]->Fill(TTrueConv2.Angle(TTrue.Vect())*TMath::RadToDeg());
	}
      }else if (fabs(TTrue.CosTheta())>0.82){
	histo_vector_oneE[3]->Fill((*true_conv_Vtx_z)[0]);
	histo_vector_oneE[5]->Fill(sqrt((*true_conv_Vtx_x)[0]*(*true_conv_Vtx_x)[0]+(*true_conv_Vtx_y)[0]*(*true_conv_Vtx_y)[0]));
	histo_vector_oneE[7]->Fill(TTrueConv1.Angle(TTrueConv2.Vect())*TMath::RadToDeg());
	if(TTrueConv1.E()>TTrueConv2.E()){
	  histo_vector_oneE[9]->Fill(TTrueConv1.E()/TTrue.E());
	  histo_vector_oneE[11]->Fill(TTrueConv2.E()/TTrue.E());
	}else{
	  histo_vector_oneE[9]->Fill(TTrueConv2.E()/TTrue.E());
	  histo_vector_oneE[11]->Fill(TTrueConv1.E()/TTrue.E());
	}
	if(TTrueConv1.Angle(TTrue.Vect())>TTrueConv2.Angle(TTrue.Vect())){
	  histo_vector_oneE[13]->Fill(TTrueConv1.Angle(TTrue.Vect())*TMath::RadToDeg());
	}else{
	  histo_vector_oneE[13]->Fill(TTrueConv2.Angle(TTrue.Vect())*TMath::RadToDeg());
	}
      }
    }



    delta_energy_correction=delta_energy_correction_endcap;
    if(fabs(TTrue.CosTheta())<0.78){
      delta_energy_correction=delta_energy_correction_barrel;
    }else if (fabs(TTrue.CosTheta())<0.83){
      delta_energy_correction=delta_energy_correction_TR;
    }

    if(ntrue>1){
      std::cout<<"claims to have more than one true particle "<<ntrue<<" "<<!is_conversion<<std::endl;
    }
    int index_reco=-1;
    int index_recoE=-1;
    double AngleMin=2*TMath::Pi();
    double EnergyOffsetMin=10.;
    unsigned int nreco=0;
    int index_reco_lead1=-1;
    double E_lead1=-1;
    double E_lead2=-1;
    int index_reco_lead2=-1;
    //for conversion recovery --> look at the two leading photons, in case those are less than two degrees appart from each other, sum them up 
    //then proceed to check the efficiency of the combined candidate --> in reality need exactly two photons and NOTHING else really around them
    //skip this check for now 
    for(unsigned int i=0;i<recoE->size();i++){
      if(abs((*recoPDGID)[i])==particle_type){
	if((*recoE)[i]>E_lead1){
	  index_reco_lead2=index_reco_lead1;
	  E_lead2=E_lead1;
	  index_reco_lead1=i;
	  E_lead1=(*recoE)[i];
	}else if( (*recoE)[i]>E_lead2 ){
	  index_reco_lead2=i;
	  E_lead2=(*recoE)[i];
	}
	/*	TLorentzVector temp(0,0,0,0);
		temp.SetPxPyPzE((*recoPx)[i],(*recoPy)[i],(*recoPz)[i],(*recoE)[i]);
		nreco+=1;
		if( (TTrue.Angle(temp.Vect())*TMath::RadToDeg())<1.0){
		if( (fabs(TTrue.E()-temp.E())/TTrue.E())<EnergyOffsetMin){
		EnergyOffsetMin=fabs(TTrue.E()-temp.E())/TTrue.E();
		index_recoE=i;
		}
		}
		if(TTrue.Angle(temp.Vect())<AngleMin){
		AngleMin=TTrue.Angle(temp.Vect());
		index_reco=i;
		}
	}*/
      }
    }
    //leading candidate particle type is defined as the default reco candidate
    index_reco=index_reco_lead1;
    TLorentzVector TRecoComb(0,0,0,0);
    TLorentzVector TReco1;
    TLorentzVector TReco2;
    if(index_reco_lead2!=-1 && particle_type == 22){
      TReco1.SetPxPyPzE((*recoPx)[index_reco_lead1],(*recoPy)[index_reco_lead1],(*recoPz)[index_reco_lead1],(*recoE)[index_reco_lead1]);
      TReco2.SetPxPyPzE((*recoPx)[index_reco_lead2],(*recoPy)[index_reco_lead2],(*recoPz)[index_reco_lead2],(*recoE)[index_reco_lead2]);
      //means two candidate particles have been found, now combine them
      //--> used for photon conversion recovery
      TRecoComb.SetPxPyPzE((*recoPx)[index_reco_lead1]+(*recoPx)[index_reco_lead2],(*recoPy)[index_reco_lead1]+(*recoPy)[index_reco_lead2],(*recoPz)[index_reco_lead1]+(*recoPz)[index_reco_lead2],(*recoE)[index_reco_lead1]+(*recoE)[index_reco_lead2]);
      if (fabs(TTrue.CosTheta())<0.78){
	histo_vector_oneE[14]->Fill(TReco1.Angle(TReco2.Vect())*TMath::RadToDeg());
	if(is_conversion){
	  histo_vector_oneE[16]->Fill(TReco1.Angle(TReco2.Vect())*TMath::RadToDeg());
	}
      }else if (fabs(TTrue.CosTheta())>0.82){
	histo_vector_oneE[15]->Fill(TReco1.Angle(TReco2.Vect())*TMath::RadToDeg());
	if(is_conversion){
	  histo_vector_oneE[17]->Fill(TReco1.Angle(TReco2.Vect())*TMath::RadToDeg());
	}
      }
    }

    //if(nreco>1){
    //std::cout<<"reco claims to have more than one true particle/conv "<<nreco<<" "<<is_conversion<<std::endl;
    //}
    if(index_reco!=index_recoE){
      TLorentzVector TRecoA(0,0,0,0);
      TLorentzVector TRecoE(0,0,0,0);
      if(index_reco!=-1){
	TRecoA.SetPxPyPzE((*recoPx)[index_reco],(*recoPy)[index_reco],(*recoPz)[index_reco],(*recoE)[index_reco]);
      }
      if(index_recoE!=-1){
	TRecoE.SetPxPyPzE((*recoPx)[index_recoE],(*recoPy)[index_recoE],(*recoPz)[index_recoE],(*recoE)[index_recoE]);
      }
    }
    bool pass_efficiency_1_degree=false;
    bool pass_efficiency_1_degree_Energy=false;
    
    bool pass_efficiency_1_degree_w_convMerged=false;
    bool pass_efficiency_1_degree_Energy_w_convMerged=false;
    
    bool pass_efficiency_1_degree_allPh1deg_Merged=false;
    bool pass_efficiency_1_degree_Energy_allPh1deg_Merged=false;

    //else would need to loop over two closeby particles
    if(!is_conversion){
      TLorentzVector TReco(0,0,0,0);
      if(index_reco!=-1){
	TReco.SetPxPyPzE((*recoPx)[index_reco],(*recoPy)[index_reco],(*recoPz)[index_reco],(*recoE)[index_reco]);
	if( (TTrue.Angle(TReco.Vect())*TMath::RadToDeg())<1.0){
	  pass_efficiency_1_degree=true;
	  histo_vector_oneE[0]->Fill((TReco.Energy()+delta_energy_correction)/TTrue.Energy());
	  histo_vector_oneE[1]->Fill(TReco.Pt()/TTrue.Pt());
	  if(pass_efficiency_1_degree){
	    histo_vector_oneE[0]->Fill((TReco.Energy()+delta_energy_correction)/TTrue.Energy());
	  }
	  if(particle_type==22){//photon energy requirement
	    if( fabs(TReco.Energy()+delta_energy_correction-TTrue.Energy())<(5.0*sigmaE_rel(TTrue.Energy(),TTrue.CosTheta())*TTrue.Energy()) ){//Delta E/E_true<(5*sigma(E))=5*0.18/sqrt{E}
	      pass_efficiency_1_degree_Energy=true;
	    }
	  }else{//charged particle energy requirement
	    if( fabs(TReco.Pt()-TTrue.Pt())<(0.05*TTrue.Pt())){//Delta E/E_true<(5*sigma(E))=5*0.18/sqrt{E}
	      pass_efficiency_1_degree_Energy=true;
	    }
	  }
	}
      }
      if(particle_type==22){
	//if reconstructed candidates far appart, don't do merging
	if(index_reco_lead2!=-1 && TReco1.Angle(TReco2.Vect())*TMath::RadToDeg()<2.0){
	  if( (TTrue.Angle(TRecoComb.Vect())*TMath::RadToDeg())<1.0){
	    pass_efficiency_1_degree_allPh1deg_Merged=true;
	    if( fabs(TRecoComb.Energy()+delta_energy_correction-TTrue.Energy())<(5.0*sigmaE_rel(TTrue.Energy(),TTrue.CosTheta())*TTrue.Energy()) ){//Delta E/E_true<(5*sigma(E))=5*0.18/sqrt{E}
	     pass_efficiency_1_degree_Energy_allPh1deg_Merged =true;
	    }
	  }
	}else{
	  pass_efficiency_1_degree_allPh1deg_Merged=pass_efficiency_1_degree;
	  pass_efficiency_1_degree_Energy_allPh1deg_Merged=pass_efficiency_1_degree_Energy;
	}
      }
      //energy dependent histograms, cut on cosTheta<0.95 to avoid artificial forward bias
      if(fabs(TTrue.CosTheta())<0.95){
	tEff_vector_allE[0]->Fill(pass_efficiency_1_degree,TTrue.Energy());
	tEff_vector_allE[1]->Fill(pass_efficiency_1_degree_Energy,TTrue.Energy());
      }
      tEff_vector_oneE[0]->Fill(pass_efficiency_1_degree,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[1]->Fill(pass_efficiency_1_degree_Energy,TTrue.Theta()*TMath::RadToDeg());
      //for true events without conversion, keep previous comparison to leading particle
      pass_efficiency_1_degree_w_convMerged=pass_efficiency_1_degree;
      pass_efficiency_1_degree_Energy_w_convMerged=pass_efficiency_1_degree_Energy;

    }else{
      //conversion are ONLY of interest (and flag filled) for photons
      TLorentzVector TReco(0,0,0,0);
      if(index_reco!=-1){
	TReco.SetPxPyPzE((*recoPx)[index_reco],(*recoPy)[index_reco],(*recoPz)[index_reco],(*recoE)[index_reco]);
	if( (TTrue.Angle(TReco.Vect())*TMath::RadToDeg())<1.0){
	  pass_efficiency_1_degree=true;
	  if( fabs(TReco.Energy()+delta_energy_correction-TTrue.Energy())<(5.0*sigmaE_rel(TTrue.Energy(),TTrue.CosTheta())*TTrue.Energy()) ){//Delta E/E_true<(5*sigma(E))=5*0.18/sqrt{E}
	    pass_efficiency_1_degree_Energy=true;
	  }
	}
      }
      //found two photons to recover, else maybe photons so close that they
      //merged in one candidate, then use results from above
      //if reconstructed candidates far appart, don't do merging
      if(index_reco_lead2!=-1 && TReco1.Angle(TReco2.Vect())*TMath::RadToDeg()<2.0){
	if( (TTrue.Angle(TRecoComb.Vect())*TMath::RadToDeg())<1.0){
	  pass_efficiency_1_degree_w_convMerged=true;
	  pass_efficiency_1_degree_allPh1deg_Merged=true;
	  if( fabs(TRecoComb.Energy()+delta_energy_correction-TTrue.Energy())<(5.0*sigmaE_rel(TTrue.Energy(),TTrue.CosTheta())*TTrue.Energy()) ){//Delta E/E_true<(5*sigma(E))=5*0.18/sqrt{E}
	    pass_efficiency_1_degree_Energy_w_convMerged=true;
	    pass_efficiency_1_degree_Energy_allPh1deg_Merged =true;
	  }
	}
      }else{
	pass_efficiency_1_degree_w_convMerged=pass_efficiency_1_degree;
	pass_efficiency_1_degree_Energy_w_convMerged=pass_efficiency_1_degree_Energy;
	pass_efficiency_1_degree_allPh1deg_Merged=pass_efficiency_1_degree;
	pass_efficiency_1_degree_Energy_allPh1deg_Merged=pass_efficiency_1_degree_Energy;
      }
      //TEfficiencies only for conversions
      if(fabs(TTrue.CosTheta())<0.95){
	tEff_vector_allE[8]->Fill(pass_efficiency_1_degree_allPh1deg_Merged,TTrue.Energy());
	tEff_vector_allE[9]->Fill(pass_efficiency_1_degree_Energy_allPh1deg_Merged,TTrue.Energy());
	tEff_vector_allE[10]->Fill(pass_efficiency_1_degree,TTrue.Energy());
	tEff_vector_allE[11]->Fill(pass_efficiency_1_degree_Energy,TTrue.Energy());
      }
      tEff_vector_oneE[8]->Fill(pass_efficiency_1_degree_allPh1deg_Merged,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[9]->Fill(pass_efficiency_1_degree_Energy_allPh1deg_Merged,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[10]->Fill(pass_efficiency_1_degree,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[11]->Fill(pass_efficiency_1_degree_Energy,TTrue.Theta()*TMath::RadToDeg());
    }
    //here both conversions and no conversions are filled, 
    //conversions and no conversions are filled by here
    if(particle_type==22 && fabs(TTrue.CosTheta())<0.95){
      tEff_vector_allE[2]->Fill(pass_efficiency_1_degree,TTrue.Energy());
      tEff_vector_allE[3]->Fill(pass_efficiency_1_degree_Energy,TTrue.Energy());
      //here compare leading reco vs gen for non conversions, but merged for conversions
      tEff_vector_allE[4]->Fill(pass_efficiency_1_degree_w_convMerged,TTrue.Energy());
      tEff_vector_allE[5]->Fill(pass_efficiency_1_degree_Energy_w_convMerged,TTrue.Energy());
      //merge leading candidates within 1 degree, even if even identified as not converted
      tEff_vector_allE[6]->Fill(pass_efficiency_1_degree_allPh1deg_Merged,TTrue.Energy());
      tEff_vector_allE[7]->Fill(pass_efficiency_1_degree_Energy_allPh1deg_Merged,TTrue.Energy());
    }
    if(particle_type==22){
      tEff_vector_oneE[2]->Fill(pass_efficiency_1_degree,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[3]->Fill(pass_efficiency_1_degree_Energy,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[4]->Fill(pass_efficiency_1_degree_w_convMerged,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[5]->Fill(pass_efficiency_1_degree_Energy_w_convMerged,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[6]->Fill(pass_efficiency_1_degree_allPh1deg_Merged,TTrue.Theta()*TMath::RadToDeg());
      tEff_vector_oneE[7]->Fill(pass_efficiency_1_degree_Energy_allPh1deg_Merged,TTrue.Theta()*TMath::RadToDeg());
    }
    
  }
  std::cout<<"at end of file, energy vector response "<<histo_vector_oneE[0]->GetMean()<<std::endl;

  if(particle_type==22){
    //here we should per file fill one bin, so nothing should be overwritten
    for(int i=1;i<=histo_vector_allE[0]->GetNbinsX();i++){
      if(count_all_events_perEBin[i-1]!=0){
	histo_vector_allE[0]->SetBinContent(i,count_all_ConvEvents_perEBin[i-1]/count_all_events_perEBin[i-1]);
      }
    }
    for(int i=1;i<=histo_vector_oneE[18]->GetNbinsX();i++){
      if(count_all_events_perThetaBin[i-1]!=0){
	histo_vector_oneE[18]->SetBinContent(i,count_all_ConvEvents_perThetaBin[i-1]/count_all_events_perThetaBin[i-1]);
      }
    }
  }
  delete tree;
}

std::vector<float> sigmaGaussianFits2Sigma(TH1F* hist, float limit_high, float limit_low){
  float sigma = 1.0;
  float sigmaError = 1.0;

  TF1* fit_raw_range_0_8_to_1_2 = new TF1("fit_sigma_range_0_8_to_1_2", "gaus",limit_high,limit_low);

  fit_raw_range_0_8_to_1_2->SetParameter(0,hist->GetBinContent(hist->GetMaximumBin()));
  std::cout<<"start with maximum "<<hist->GetMaximumBin()<<" value "<<hist->GetBinContent(hist->GetMaximumBin())<<std::endl;
  
  fit_raw_range_0_8_to_1_2->SetParLimits(0,0.50*hist->GetBinContent(hist->GetMaximumBin()),1.50*hist->GetBinContent(hist->GetMaximumBin()));


  hist->Fit("fit_sigma_range_0_8_to_1_2","R");
  //sigma is the parameter 2, mean is parameter 1
  float sigma_temp=fit_raw_range_0_8_to_1_2->GetParameter(2);
  float mean_temp=fit_raw_range_0_8_to_1_2->GetParameter(1);
  TF1* fit_raw_range_two_sigma = new TF1("fit_sigma_range_two_sigma", "gaus",mean_temp-2.*sigma_temp,mean_temp+2.*sigma_temp);
  fit_raw_range_0_8_to_1_2->SetParLimits(0,0.50*hist->GetBinContent(hist->GetMaximumBin()),1.50*hist->GetBinContent(hist->GetMaximumBin()));
  fit_raw_range_0_8_to_1_2->SetParLimits (1,mean_temp-2.*sigma_temp,mean_temp+2.*sigma_temp);
  hist->Fit("fit_sigma_range_two_sigma","R");
  sigma=fit_raw_range_two_sigma->GetParameter(2);
  sigmaError=fit_raw_range_two_sigma->GetParError(2);
  //gaussian mean
  //fit_raw_range_two_sigma->SetParLimits (1,0.75,1.25);

  std::vector<float> SigmaWError;
  SigmaWError.push_back(sigma);
  SigmaWError.push_back(sigmaError);
  float sigma_before=fit_raw_range_0_8_to_1_2->GetParameter(2);
  float sigma_beforeError=fit_raw_range_0_8_to_1_2->GetParError(2);
  float mean_before=fit_raw_range_0_8_to_1_2->GetParameter(1);
  for(unsigned int i=0;i<1000;i++){
    TF1* fit_raw_range_temp = new TF1("fit_sigma_range_temp", "gaus",mean_before-2.*sigma_before,mean_before+2.*sigma_before);
    hist->Fit("fit_sigma_range_temp","R");
    float sigma_after=fit_raw_range_temp->GetParameter(2);
    float sigma_afterError=fit_raw_range_temp->GetParError(2);
    if((sigma_after/sigma_before)<1.05 && (sigma_after/sigma_before)>0.95){
      std::cout<<"break point RES/mean/rescaledRES/err "<<sigma_after<<"/"<<fit_raw_range_temp->GetParameter(1)<<"/"<<sigma_after/fit_raw_range_temp->GetParameter(1)<<"/"<<sigma_afterError/fit_raw_range_temp->GetParameter(1)<<std::endl;
      sigma_before=sigma_after/fit_raw_range_temp->GetParameter(1);
      sigma_beforeError=sigma_afterError/fit_raw_range_temp->GetParameter(1);
      break;
    }
    sigma_before=sigma_after;
    sigma_beforeError=sigma_afterError;
    mean_before=fit_raw_range_temp->GetParameter(1);
 }
  SigmaWError.push_back(sigma_before);
  SigmaWError.push_back(sigma_beforeError);

  return SigmaWError;

}

int plotReducedShower(){
  CLICdpStyle();

  
    string label_legend= "#pi^{-}, E_{true}=10 GeV";
    
    TFile* file_photon_ph1=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph1_ILC181011_CT.root");
    TFile* file_photon_ph5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph5_ILC181011_CT.root");
    TFile* file_photon_ph10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph10_ILC181011_CT.root");
    TFile* file_photon_ph15=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph15_ILC181011_CT.root");
    TFile* file_photon_ph30=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph30_ILC181011_CT.root");
    TFile* file_photon_ph50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph50_ILC181011_CT.root");
    TFile* file_photon_ph100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph100_ILC181011_CT.root");
    TFile* file_photon_ph200=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph200_ILC181011_CT.root");
    TFile* file_photon_ph500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph500_ILC181011_CT.root");
    TFile* file_photon_ph1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph1000_ILC181011_CT.root");
    TFile* file_photon_ph1500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_ph1500_ILC181011_CT.root");
 
  //for photon events 1/5/10/15/30/50/100/200/500/1000/1500
  const char* final_histo_name_photons="/eos/user/w/weberma2/data/validation181011/pionStudy_ph1_to_ph1500_ILC181011_finalhistos_deg1_match.root";

  TFile* file_histogram_photons=new TFile(final_histo_name_photons,"recreate");
    
  int n_bins_highRES=100;
  int n_bins_highTRRES=75;
  float lim_rel_energy_lowRES=0.75;
  float lim_rel_energy_highRES=1.25;
  int n_bins_high=50;
  float lim_rel_energy_low=0.5;
  float lim_rel_energy_high=1.5;

  int test_signal_particle = 22;


  int test_signal_particle_true = 22;
  int test_signal_particle_reco = 22;

  TGraphErrors* gre_ph_E_reco_over_E_true_0_78= new TGraphErrors(11);
  gre_ph_E_reco_over_E_true_0_78->SetName("gre_ph_E_reco_over_E_true_0_78");
  gre_ph_E_reco_over_E_true_0_78->SetLineColor(kBlack); 
  gre_ph_E_reco_over_E_true_0_78->SetMarkerColor(kBlack); 
  gre_ph_E_reco_over_E_true_0_78->SetMarkerStyle(kFullCircle); 
  gre_ph_E_reco_over_E_true_0_78->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_ph_E_reco_over_E_true_0_78->GetXaxis()->SetTitle("E(#gamma_{true}) [GeV]");

  TGraphErrors* gre_ph_E_reco_over_E_true_0_78_to_0_83 = new TGraphErrors(11);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetName("gre_ph_E_reco_over_E_true_0_78_to_0_83");
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetLineColor(kRed); 
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetMarkerColor(kRed); 
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetMarkerStyle(kFullSquare); 
  gre_ph_E_reco_over_E_true_0_78_to_0_83->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_ph_E_reco_over_E_true_0_78_to_0_83->GetXaxis()->SetTitle("E(#gamma_{true}) [GeV]");

  TGraphErrors* gre_ph_E_reco_over_E_true_0_83_to_0_94 = new TGraphErrors(11);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetName("gre_ph_E_reco_over_E_true_0_83_to_0_94");
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetLineColor(kBlue); 
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetMarkerColor(kBlue); 
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetMarkerStyle(kFullTriangleUp); 
  gre_ph_E_reco_over_E_true_0_83_to_0_94->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_ph_E_reco_over_E_true_0_83_to_0_94->GetXaxis()->SetTitle("E(#gamma_{true}) [GeV]");


  std::cout<<"RES ph1"<<std::endl;
  TH1F* h_ph1_reco_E_over_true_E_match_0_78 = new TH1F("h_ph5_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph5_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph5_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph1;
  TH_vec_ERES_ph1.push_back(h_ph1_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph1.push_back(h_ph1_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph1.push_back(h_ph1_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph1.size();i++){
    TH_vec_ERES_ph1[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph1, TH_vec_ERES_ph1,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph1_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(0,1,100.*sigma_ph1_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(0,0,100.*sigma_ph1_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(0,1,100.*sigma_ph1_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(0,0,100.*sigma_ph1_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(0,1,100.*sigma_ph1_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(0,0,100.*sigma_ph1_0_83_to_0_94[1]);

  std::cout<<"RES ph5"<<std::endl;
  TH1F* h_ph5_reco_E_over_true_E_match_0_78 = new TH1F("h_ph5_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph5_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph5_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph5_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph5_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph5;
  TH_vec_ERES_ph5.push_back(h_ph5_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph5.push_back(h_ph5_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph5.push_back(h_ph5_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph5.size();i++){
    TH_vec_ERES_ph5[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph5, TH_vec_ERES_ph5,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph5_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph5[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph5_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph5[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph5_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph5[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(1,5,100.*sigma_ph5_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(1,0,100.*sigma_ph5_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(1,5,100.*sigma_ph5_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(1,0,100.*sigma_ph5_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(1,5,100.*sigma_ph5_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(1,0,100.*sigma_ph5_0_83_to_0_94[1]);

  std::cout<<"RES ph10"<<std::endl;
  TH1F* h_ph10_reco_E_over_true_E_match_0_78 = new TH1F("h_ph10_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph10_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph10_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph10_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph10_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph10;
  TH_vec_ERES_ph10.push_back(h_ph10_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph10.push_back(h_ph10_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph10.push_back(h_ph10_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph10.size();i++){
    TH_vec_ERES_ph10[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph10, TH_vec_ERES_ph10,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph10_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph10[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph10_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph10[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph10_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph10[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(2,10,100.*sigma_ph10_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(2,0,100.*sigma_ph10_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(2,10,100.*sigma_ph10_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(2,0,100.*sigma_ph10_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(2,10,100.*sigma_ph10_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(2,0,100.*sigma_ph10_0_83_to_0_94[1]);

  std::cout<<"RES ph15"<<std::endl;
  TH1F* h_ph15_reco_E_over_true_E_match_0_78 = new TH1F("h_ph15_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph15_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph15_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph15_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph15_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph15;
  TH_vec_ERES_ph15.push_back(h_ph15_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph15.push_back(h_ph15_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph15.push_back(h_ph15_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph15.size();i++){
    TH_vec_ERES_ph15[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph15, TH_vec_ERES_ph15,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph15_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph15[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph15_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph15[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph15_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph15[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(3,15,100.*sigma_ph15_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(3,0,100.*sigma_ph15_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(3,15,100.*sigma_ph15_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(3,0,100.*sigma_ph15_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(3,15,100.*sigma_ph15_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(3,0,100.*sigma_ph15_0_83_to_0_94[1]);

  std::cout<<"RES ph30"<<std::endl;
  TH1F* h_ph30_reco_E_over_true_E_match_0_78 = new TH1F("h_ph30_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph30_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph30_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph30_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph30_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph30;
  TH_vec_ERES_ph30.push_back(h_ph30_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph30.push_back(h_ph30_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph30.push_back(h_ph30_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph30.size();i++){
    TH_vec_ERES_ph30[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph30, TH_vec_ERES_ph30,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph30_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph30[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph30_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph30[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph30_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph30[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(4,30,100.*sigma_ph30_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(4,0,100.*sigma_ph30_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(4,30,100.*sigma_ph30_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(4,0,100.*sigma_ph30_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(4,30,100.*sigma_ph30_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(4,0,100.*sigma_ph30_0_83_to_0_94[1]);

  std::cout<<"RES ph50"<<std::endl;
  TH1F* h_ph50_reco_E_over_true_E_match_0_78 = new TH1F("h_ph50_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph50_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph50_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph50_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph50_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph50;
  TH_vec_ERES_ph50.push_back(h_ph50_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph50.push_back(h_ph50_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph50.push_back(h_ph50_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph50.size();i++){
    TH_vec_ERES_ph50[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph50, TH_vec_ERES_ph50,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph50_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph50[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph50_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph50[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph50_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph50[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(5,50,100.*sigma_ph50_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(5,0,100.*sigma_ph50_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(5,50,100.*sigma_ph50_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(5,0,100.*sigma_ph50_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(5,50,100.*sigma_ph50_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(5,0,100.*sigma_ph50_0_83_to_0_94[1]);


  std::cout<<"RES ph100"<<std::endl;
  TH1F* h_ph100_reco_E_over_true_E_match_0_78 = new TH1F("h_ph100_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph100_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph100_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph100_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph100_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph100;
  TH_vec_ERES_ph100.push_back(h_ph100_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph100.push_back(h_ph100_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph100.push_back(h_ph100_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph100.size();i++){
    TH_vec_ERES_ph100[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph100, TH_vec_ERES_ph100,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph100_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph100[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph100_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph100[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph100_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph100[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(6,100,100.*sigma_ph100_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(6,0,100.*sigma_ph100_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(6,100,100.*sigma_ph100_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(6,0,100.*sigma_ph100_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(6,100,100.*sigma_ph100_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(6,0,100.*sigma_ph100_0_83_to_0_94[1]);


  std::cout<<"RES ph200"<<std::endl;
  TH1F* h_ph200_reco_E_over_true_E_match_0_78 = new TH1F("h_ph200_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph200_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph200_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph200_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph200_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph200;
  TH_vec_ERES_ph200.push_back(h_ph200_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph200.push_back(h_ph200_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph200.push_back(h_ph200_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph200.size();i++){
    TH_vec_ERES_ph200[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph200, TH_vec_ERES_ph200,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph200_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph200[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph200_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph200[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph200_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph200[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(7,200,100.*sigma_ph200_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(7,0,100.*sigma_ph200_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(7,200,100.*sigma_ph200_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(7,0,100.*sigma_ph200_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(7,200,100.*sigma_ph200_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(7,0,100.*sigma_ph200_0_83_to_0_94[1]);

  //curves gets tighter, increase bins for more statistics in fit range
  n_bins_highRES=200;
  n_bins_highTRRES=150;
  std::cout<<"RES ph500"<<std::endl;
  TH1F* h_ph500_reco_E_over_true_E_match_0_78 = new TH1F("h_ph500_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph500_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph500_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph500_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph500_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph500;
  TH_vec_ERES_ph500.push_back(h_ph500_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph500.push_back(h_ph500_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph500.push_back(h_ph500_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph500.size();i++){
    TH_vec_ERES_ph500[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph500, TH_vec_ERES_ph500,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph500_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph500[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph500_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph500[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph500_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph500[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(8,500,100.*sigma_ph500_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(8,0,100.*sigma_ph500_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(8,500,100.*sigma_ph500_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(8,0,100.*sigma_ph500_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(8,500,100.*sigma_ph500_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(8,0,100.*sigma_ph500_0_83_to_0_94[1]);

  std::cout<<"RES ph1000"<<std::endl;
  TH1F* h_ph1000_reco_E_over_true_E_match_0_78 = new TH1F("h_ph1000_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1000_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph1000_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1000_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph1000_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph1000;
  TH_vec_ERES_ph1000.push_back(h_ph1000_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph1000.push_back(h_ph1000_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph1000.push_back(h_ph1000_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph1000.size();i++){
    TH_vec_ERES_ph1000[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph1000, TH_vec_ERES_ph1000,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph1000_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1000[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1000_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1000[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1000_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1000[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(9,1000,100.*sigma_ph1000_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(9,0,100.*sigma_ph1000_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(9,1000,100.*sigma_ph1000_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(9,0,100.*sigma_ph1000_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(9,1000,100.*sigma_ph1000_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(9,0,100.*sigma_ph1000_0_83_to_0_94[1]);

  std::cout<<"RES ph1500"<<std::endl;
  TH1F* h_ph1500_reco_E_over_true_E_match_0_78 = new TH1F("h_ph1500_reco_E_over_true_E_match_0_78","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1500_reco_E_over_true_E_match_0_78_to_0_83 = new TH1F("h_ph1500_reco_E_over_true_E_match_0_78_to_0_83","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_ph1500_reco_E_over_true_E_match_0_83_to_0_94 = new TH1F("h_ph1500_reco_E_over_true_E_match_0_83_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_ph1500;
  TH_vec_ERES_ph1500.push_back(h_ph1500_reco_E_over_true_E_match_0_78);
  TH_vec_ERES_ph1500.push_back(h_ph1500_reco_E_over_true_E_match_0_78_to_0_83);
  TH_vec_ERES_ph1500.push_back(h_ph1500_reco_E_over_true_E_match_0_83_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_ph1500.size();i++){
    TH_vec_ERES_ph1500[i]->Sumw2();
  }
  fillRESTreeHistograms(file_photon_ph1500, TH_vec_ERES_ph1500,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_ph1500_0_78=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1500[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1500_0_78_to_0_83=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1500[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_ph1500_0_83_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_ph1500[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_ph_E_reco_over_E_true_0_78->SetPoint(10,1500,100.*sigma_ph1500_0_78[0]);
  gre_ph_E_reco_over_E_true_0_78->SetPointError(10,0,100.*sigma_ph1500_0_78[1]);

  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPoint(10,1500,100.*sigma_ph1500_0_78_to_0_83[0]);
  gre_ph_E_reco_over_E_true_0_78_to_0_83->SetPointError(10,0,100.*sigma_ph1500_0_78_to_0_83[1]);

  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPoint(10,1500,100.*sigma_ph1500_0_83_to_0_94[0]);
  gre_ph_E_reco_over_E_true_0_83_to_0_94->SetPointError(10,0,100.*sigma_ph1500_0_83_to_0_94[1]);


    
  int n_bins=25;
  float lim_theta_low = 0;
  float lim_theta_high=180;
  

  const unsigned int nbins_teff_ph = 11; 
  double xbins_teff_ph[nbins_teff_ph+1]={0.75,2.5,7.5,12.5,20.,40.,75.,150.,300.,750.,1250.,1550.};
  
  TH1F* h_photonConvRateVsTrueE = new TH1F("h_PhotonConvRateVsTrueE_ph1_1500","", nbins_teff_ph,xbins_teff_ph);
  std::vector<TH1F*>hist_vec_all_ph;
  hist_vec_all_ph.push_back(h_photonConvRateVsTrueE);
  for(unsigned int i=0;i<hist_vec_all_ph.size();i++){
    hist_vec_all_ph[i]->Sumw2();
  }

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_wConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_wConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_wConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_wConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_wConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_wConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedConv_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedConv_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedConv_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");

  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay = new TEfficiency("tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay","", nbins_teff_ph,xbins_teff_ph);
  tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay->SetTitle(";true photon energy [GeV];photon ID efficiency");
  
  std::vector<TEfficiency*> TEff_vector_all_ph_noOverlay;
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_wConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_wConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedConv_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_w_mergedPh_1_deg_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvMergedPh_1_deg_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay);
  TEff_vector_all_ph_noOverlay.push_back(tEff_PhotonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_OnlyConvUnmergedPh_1_deg_noOverlay);

  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");


    
  std::vector<TEfficiency*> TEff_vector_ph1_noOverlay;
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1_noOverlay);
  TEff_vector_ph1_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1_noOverlay);

  TH1F* h_ph1_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph1_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph1_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph1_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);

  float lim_convVtx_Z_min=-2310;
  float lim_convVtx_Z_max=-2310;

  TH1F* h_ph1_convVtx_Z_barrel = new TH1F("h_ph1_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph1_convVtx_Z_endcap = new TH1F("h_ph1_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);

  float lim_convVtx_R_min=0;
  float lim_convVtx_R_max=1500;

  TH1F* h_ph1_convVtx_R_barrel = new TH1F("h_ph1_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph1_convVtx_R_endcap = new TH1F("h_ph1_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);

  //limits in degrees
  float lim_convAngle_min=0;
  float lim_convAngle_max=5.0;
  TH1F* h_ph1_convAngle_barrel = new TH1F("h_ph1_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph1_convAngle_endcap = new TH1F("h_ph1_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);

  //limits in degrees
  float lim_convE1rel_min=0.4;
  float lim_convE1rel_max=1.0;
  TH1F* h_ph1_convE1_over_ETrue_barrel = new TH1F("h_ph1_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph1_convE1_over_ETrue_endcap = new TH1F("h_ph1_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);

  float lim_convE2rel_min=0.0;
  float lim_convE2rel_max=0.6;
  TH1F* h_ph1_convE2_over_ETrue_barrel = new TH1F("h_ph1_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph1_convE2_over_ETrue_endcap = new TH1F("h_ph1_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);

  float lim_convDelta_min=0.0;
  float lim_convDelta_max=5.0;
  TH1F* h_ph1_conv_DeltaAng_max_barrel = new TH1F("h_ph1_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph1_conv_DeltaAng_max_endcap = new TH1F("h_ph1_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);

  float lim_recoAngle_min=0;
  float lim_recoAngle_max=15.0;

  TH1F* h_ph1_DeltaPh12_barrel = new TH1F("h_ph1_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1_DeltaPh12_endcap = new TH1F("h_ph1_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1_DeltaPh12_conv_barrel = new TH1F("h_ph1_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1_DeltaPh12_conv_endcap = new TH1F("h_ph1_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1_photonConvRateVsTrueTheta = new TH1F("h_ph1_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);


  std::vector<TH1F*>hist_vec_ph1;
  hist_vec_ph1.push_back(h_ph1_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph1.push_back(h_ph1_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph1.push_back(h_ph1_convVtx_Z_barrel);
  hist_vec_ph1.push_back(h_ph1_convVtx_Z_endcap);
  hist_vec_ph1.push_back(h_ph1_convVtx_R_barrel);
  hist_vec_ph1.push_back(h_ph1_convVtx_R_endcap);
  hist_vec_ph1.push_back(h_ph1_convAngle_barrel);
  hist_vec_ph1.push_back(h_ph1_convAngle_endcap);
  hist_vec_ph1.push_back(h_ph1_convE1_over_ETrue_barrel);
  hist_vec_ph1.push_back(h_ph1_convE1_over_ETrue_endcap);
  hist_vec_ph1.push_back(h_ph1_convE2_over_ETrue_barrel);
  hist_vec_ph1.push_back(h_ph1_convE2_over_ETrue_endcap);
  hist_vec_ph1.push_back(h_ph1_conv_DeltaAng_max_barrel);
  hist_vec_ph1.push_back(h_ph1_conv_DeltaAng_max_endcap);
  hist_vec_ph1.push_back(h_ph1_DeltaPh12_barrel);
  hist_vec_ph1.push_back(h_ph1_DeltaPh12_endcap);
  hist_vec_ph1.push_back(h_ph1_DeltaPh12_conv_barrel);
  hist_vec_ph1.push_back(h_ph1_DeltaPh12_conv_endcap);
  hist_vec_ph1.push_back(h_ph1_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph1.size();i++){
    hist_vec_ph1[i]->Sumw2();
  }
  std::cout<<"fill ph1"<<std::endl;
  fillTreeHistograms(file_photon_ph1,hist_vec_ph1,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph1_noOverlay,hist_vec_all_ph);
    
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph5_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph5_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph5_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph5_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph5_noOverlay;
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph5_noOverlay);
  TEff_vector_ph5_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph5_noOverlay);
  
  TH1F* h_ph5_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph5_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph5_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph5_convVtx_Z_barrel = new TH1F("h_ph5_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph5_convVtx_Z_endcap = new TH1F("h_ph5_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph5_convVtx_R_barrel = new TH1F("h_ph5_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph5_convVtx_R_endcap = new TH1F("h_ph5_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph5_convAngle_barrel = new TH1F("h_ph5_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph5_convAngle_endcap = new TH1F("h_ph5_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph5_convE1_over_ETrue_barrel = new TH1F("h_ph5_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph5_convE1_over_ETrue_endcap = new TH1F("h_ph5_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph5_convE2_over_ETrue_barrel = new TH1F("h_ph5_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph5_convE2_over_ETrue_endcap = new TH1F("h_ph5_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph5_conv_DeltaAng_max_barrel = new TH1F("h_ph5_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph5_conv_DeltaAng_max_endcap = new TH1F("h_ph5_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph5_DeltaPh12_barrel = new TH1F("h_ph5_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph5_DeltaPh12_endcap = new TH1F("h_ph5_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph5_DeltaPh12_conv_barrel = new TH1F("h_ph5_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph5_DeltaPh12_conv_endcap = new TH1F("h_ph5_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph5_photonConvRateVsTrueTheta = new TH1F("h_ph5_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph5;
  hist_vec_ph5.push_back(h_ph5_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph5.push_back(h_ph5_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph5.push_back(h_ph5_convVtx_Z_barrel);
  hist_vec_ph5.push_back(h_ph5_convVtx_Z_endcap);
  hist_vec_ph5.push_back(h_ph5_convVtx_R_barrel);
  hist_vec_ph5.push_back(h_ph5_convVtx_R_endcap);
  hist_vec_ph5.push_back(h_ph5_convAngle_barrel);
  hist_vec_ph5.push_back(h_ph5_convAngle_endcap);
  hist_vec_ph5.push_back(h_ph5_convE1_over_ETrue_barrel);
  hist_vec_ph5.push_back(h_ph5_convE1_over_ETrue_endcap);
  hist_vec_ph5.push_back(h_ph5_convE2_over_ETrue_barrel);
  hist_vec_ph5.push_back(h_ph5_convE2_over_ETrue_endcap);
  hist_vec_ph5.push_back(h_ph5_conv_DeltaAng_max_barrel);
  hist_vec_ph5.push_back(h_ph5_conv_DeltaAng_max_endcap);
  hist_vec_ph5.push_back(h_ph5_DeltaPh12_barrel);
  hist_vec_ph5.push_back(h_ph5_DeltaPh12_endcap);
  hist_vec_ph5.push_back(h_ph5_DeltaPh12_conv_barrel);
  hist_vec_ph5.push_back(h_ph5_DeltaPh12_conv_endcap);
  hist_vec_ph5.push_back(h_ph5_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph5.size();i++){
    hist_vec_ph5[i]->Sumw2();
  }
  std::cout<<"fill ph5"<<std::endl;
  fillTreeHistograms(file_photon_ph5,hist_vec_ph5,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph5_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph10_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph10_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph10_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph10_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph10_noOverlay;
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph10_noOverlay);
  TEff_vector_ph10_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph10_noOverlay);
    
  TH1F* h_ph10_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph10_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph10_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph10_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph10_convVtx_Z_barrel = new TH1F("h_ph10_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph10_convVtx_Z_endcap = new TH1F("h_ph10_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph10_convVtx_R_barrel = new TH1F("h_ph10_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph10_convVtx_R_endcap = new TH1F("h_ph10_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph10_convAngle_barrel = new TH1F("h_ph10_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph10_convAngle_endcap = new TH1F("h_ph10_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph10_convE1_over_ETrue_barrel = new TH1F("h_ph10_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph10_convE1_over_ETrue_endcap = new TH1F("h_ph10_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph10_convE2_over_ETrue_barrel = new TH1F("h_ph10_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph10_convE2_over_ETrue_endcap = new TH1F("h_ph10_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph10_conv_DeltaAng_max_barrel = new TH1F("h_ph10_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph10_conv_DeltaAng_max_endcap = new TH1F("h_ph10_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph10_DeltaPh12_barrel = new TH1F("h_ph10_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph10_DeltaPh12_endcap = new TH1F("h_ph10_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph10_DeltaPh12_conv_barrel = new TH1F("h_ph10_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph10_DeltaPh12_conv_endcap = new TH1F("h_ph10_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph10_photonConvRateVsTrueTheta = new TH1F("h_ph10_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph10;
  hist_vec_ph10.push_back(h_ph10_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph10.push_back(h_ph10_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph10.push_back(h_ph10_convVtx_Z_barrel);
  hist_vec_ph10.push_back(h_ph10_convVtx_Z_endcap);
  hist_vec_ph10.push_back(h_ph10_convVtx_R_barrel);
  hist_vec_ph10.push_back(h_ph10_convVtx_R_endcap);
  hist_vec_ph10.push_back(h_ph10_convAngle_barrel);
  hist_vec_ph10.push_back(h_ph10_convAngle_endcap);
  hist_vec_ph10.push_back(h_ph10_convE1_over_ETrue_barrel);
  hist_vec_ph10.push_back(h_ph10_convE1_over_ETrue_endcap);
  hist_vec_ph10.push_back(h_ph10_convE2_over_ETrue_barrel);
  hist_vec_ph10.push_back(h_ph10_convE2_over_ETrue_endcap);
  hist_vec_ph10.push_back(h_ph10_conv_DeltaAng_max_barrel);
  hist_vec_ph10.push_back(h_ph10_conv_DeltaAng_max_endcap);
  hist_vec_ph10.push_back(h_ph10_DeltaPh12_barrel);
  hist_vec_ph10.push_back(h_ph10_DeltaPh12_endcap);
  hist_vec_ph10.push_back(h_ph10_DeltaPh12_conv_barrel);
  hist_vec_ph10.push_back(h_ph10_DeltaPh12_conv_endcap);
  hist_vec_ph10.push_back(h_ph10_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph10.size();i++){
    hist_vec_ph10[i]->Sumw2();
  }
  std::cout<<"fill ph10"<<std::endl;
  fillTreeHistograms(file_photon_ph10,hist_vec_ph10,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph10_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph15_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph15_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph15_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph15_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph15_noOverlay;
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph15_noOverlay);
  TEff_vector_ph15_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph15_noOverlay);
  
  TH1F* h_ph15_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph15_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph15_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph15_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph15_convVtx_Z_barrel = new TH1F("h_ph15_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph15_convVtx_Z_endcap = new TH1F("h_ph15_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph15_convVtx_R_barrel = new TH1F("h_ph15_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph15_convVtx_R_endcap = new TH1F("h_ph15_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph15_convAngle_barrel = new TH1F("h_ph15_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph15_convAngle_endcap = new TH1F("h_ph15_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph15_convE1_over_ETrue_barrel = new TH1F("h_ph15_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph15_convE1_over_ETrue_endcap = new TH1F("h_ph15_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph15_convE2_over_ETrue_barrel = new TH1F("h_ph15_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph15_convE2_over_ETrue_endcap = new TH1F("h_ph15_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph15_conv_DeltaAng_max_barrel = new TH1F("h_ph15_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph15_conv_DeltaAng_max_endcap = new TH1F("h_ph15_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph15_DeltaPh12_barrel = new TH1F("h_ph15_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph15_DeltaPh12_endcap = new TH1F("h_ph15_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph15_DeltaPh12_conv_barrel = new TH1F("h_ph15_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph15_DeltaPh12_conv_endcap = new TH1F("h_ph15_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph15_photonConvRateVsTrueTheta = new TH1F("h_ph15_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph15;
  hist_vec_ph15.push_back(h_ph15_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph15.push_back(h_ph15_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph15.push_back(h_ph15_convVtx_Z_barrel);
  hist_vec_ph15.push_back(h_ph15_convVtx_Z_endcap);
  hist_vec_ph15.push_back(h_ph15_convVtx_R_barrel);
  hist_vec_ph15.push_back(h_ph15_convVtx_R_endcap);
  hist_vec_ph15.push_back(h_ph15_convAngle_barrel);
  hist_vec_ph15.push_back(h_ph15_convAngle_endcap);
  hist_vec_ph15.push_back(h_ph15_convE1_over_ETrue_barrel);
  hist_vec_ph15.push_back(h_ph15_convE1_over_ETrue_endcap);
  hist_vec_ph15.push_back(h_ph15_convE2_over_ETrue_barrel);
  hist_vec_ph15.push_back(h_ph15_convE2_over_ETrue_endcap);
  hist_vec_ph15.push_back(h_ph15_conv_DeltaAng_max_barrel);
  hist_vec_ph15.push_back(h_ph15_conv_DeltaAng_max_endcap);
  hist_vec_ph15.push_back(h_ph15_DeltaPh12_barrel);
  hist_vec_ph15.push_back(h_ph15_DeltaPh12_endcap);
  hist_vec_ph15.push_back(h_ph15_DeltaPh12_conv_barrel);
  hist_vec_ph15.push_back(h_ph15_DeltaPh12_conv_endcap);
  hist_vec_ph15.push_back(h_ph15_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph15.size();i++){
    hist_vec_ph15[i]->Sumw2();
  }
  std::cout<<"fill ph15"<<std::endl;
  fillTreeHistograms(file_photon_ph15,hist_vec_ph15,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph15_noOverlay,hist_vec_all_ph);
    
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph30_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph30_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph30_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph30_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph30_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph30_noOverlay;
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph30_noOverlay);
  TEff_vector_ph30_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph30_noOverlay);
  
  TH1F* h_ph30_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph30_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph30_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph30_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph30_convVtx_Z_barrel = new TH1F("h_ph30_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph30_convVtx_Z_endcap = new TH1F("h_ph30_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph30_convVtx_R_barrel = new TH1F("h_ph30_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph30_convVtx_R_endcap = new TH1F("h_ph30_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph30_convAngle_barrel = new TH1F("h_ph30_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph30_convAngle_endcap = new TH1F("h_ph30_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph30_convE1_over_ETrue_barrel = new TH1F("h_ph30_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph30_convE1_over_ETrue_endcap = new TH1F("h_ph30_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph30_convE2_over_ETrue_barrel = new TH1F("h_ph30_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph30_convE2_over_ETrue_endcap = new TH1F("h_ph30_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph30_conv_DeltaAng_max_barrel = new TH1F("h_ph30_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph30_conv_DeltaAng_max_endcap = new TH1F("h_ph30_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph30_DeltaPh12_barrel = new TH1F("h_ph30_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph30_DeltaPh12_endcap = new TH1F("h_ph30_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph30_DeltaPh12_conv_barrel = new TH1F("h_ph30_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph30_DeltaPh12_conv_endcap = new TH1F("h_ph30_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph30_photonConvRateVsTrueTheta = new TH1F("h_ph30_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph30;
  hist_vec_ph30.push_back(h_ph30_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph30.push_back(h_ph30_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph30.push_back(h_ph30_convVtx_Z_barrel);
  hist_vec_ph30.push_back(h_ph30_convVtx_Z_endcap);
  hist_vec_ph30.push_back(h_ph30_convVtx_R_barrel);
  hist_vec_ph30.push_back(h_ph30_convVtx_R_endcap);
  hist_vec_ph30.push_back(h_ph30_convAngle_barrel);
  hist_vec_ph30.push_back(h_ph30_convAngle_endcap);
  hist_vec_ph30.push_back(h_ph30_convE1_over_ETrue_barrel);
  hist_vec_ph30.push_back(h_ph30_convE1_over_ETrue_endcap);
  hist_vec_ph30.push_back(h_ph30_convE2_over_ETrue_barrel);
  hist_vec_ph30.push_back(h_ph30_convE2_over_ETrue_endcap);
  hist_vec_ph30.push_back(h_ph30_conv_DeltaAng_max_barrel);
  hist_vec_ph30.push_back(h_ph30_conv_DeltaAng_max_endcap);
  hist_vec_ph30.push_back(h_ph30_DeltaPh12_barrel);
  hist_vec_ph30.push_back(h_ph30_DeltaPh12_endcap);
  hist_vec_ph30.push_back(h_ph30_DeltaPh12_conv_barrel);
  hist_vec_ph30.push_back(h_ph30_DeltaPh12_conv_endcap);
  hist_vec_ph30.push_back(h_ph30_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph30.size();i++){
    hist_vec_ph30[i]->Sumw2();
  }

  std::cout<<"fill ph30"<<std::endl;
  fillTreeHistograms(file_photon_ph30,hist_vec_ph30,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph30_noOverlay,hist_vec_all_ph);

  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph50_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph50_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph50_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph50_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph50_noOverlay;
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph50_noOverlay);
  TEff_vector_ph50_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph50_noOverlay);
  
  TH1F* h_ph50_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph50_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph50_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph50_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph50_convVtx_Z_barrel = new TH1F("h_ph50_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph50_convVtx_Z_endcap = new TH1F("h_ph50_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph50_convVtx_R_barrel = new TH1F("h_ph50_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph50_convVtx_R_endcap = new TH1F("h_ph50_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph50_convAngle_barrel = new TH1F("h_ph50_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph50_convAngle_endcap = new TH1F("h_ph50_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph50_convE1_over_ETrue_barrel = new TH1F("h_ph50_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph50_convE1_over_ETrue_endcap = new TH1F("h_ph50_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph50_convE2_over_ETrue_barrel = new TH1F("h_ph50_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph50_convE2_over_ETrue_endcap = new TH1F("h_ph50_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph50_conv_DeltaAng_max_barrel = new TH1F("h_ph50_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph50_conv_DeltaAng_max_endcap = new TH1F("h_ph50_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph50_DeltaPh12_barrel = new TH1F("h_ph50_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph50_DeltaPh12_endcap = new TH1F("h_ph50_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph50_DeltaPh12_conv_barrel = new TH1F("h_ph50_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph50_DeltaPh12_conv_endcap = new TH1F("h_ph50_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph50_photonConvRateVsTrueTheta = new TH1F("h_ph50_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph50;
  hist_vec_ph50.push_back(h_ph50_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph50.push_back(h_ph50_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph50.push_back(h_ph50_convVtx_Z_barrel);
  hist_vec_ph50.push_back(h_ph50_convVtx_Z_endcap);
  hist_vec_ph50.push_back(h_ph50_convVtx_R_barrel);
  hist_vec_ph50.push_back(h_ph50_convVtx_R_endcap);
  hist_vec_ph50.push_back(h_ph50_convAngle_barrel);
  hist_vec_ph50.push_back(h_ph50_convAngle_endcap);
  hist_vec_ph50.push_back(h_ph50_convE1_over_ETrue_barrel);
  hist_vec_ph50.push_back(h_ph50_convE1_over_ETrue_endcap);
  hist_vec_ph50.push_back(h_ph50_convE2_over_ETrue_barrel);
  hist_vec_ph50.push_back(h_ph50_convE2_over_ETrue_endcap);
  hist_vec_ph50.push_back(h_ph50_conv_DeltaAng_max_barrel);
  hist_vec_ph50.push_back(h_ph50_conv_DeltaAng_max_endcap);
  hist_vec_ph50.push_back(h_ph50_DeltaPh12_barrel);
  hist_vec_ph50.push_back(h_ph50_DeltaPh12_endcap);
  hist_vec_ph50.push_back(h_ph50_DeltaPh12_conv_barrel);
  hist_vec_ph50.push_back(h_ph50_DeltaPh12_conv_endcap);
  hist_vec_ph50.push_back(h_ph50_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph50.size();i++){
    hist_vec_ph50[i]->Sumw2();
  }
  std::cout<<"fill ph50"<<std::endl;
  fillTreeHistograms(file_photon_ph50,hist_vec_ph50,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph50_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph100_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph100_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph100_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph100_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph100_noOverlay;
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph100_noOverlay);
  TEff_vector_ph100_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph100_noOverlay);
  
  TH1F* h_ph100_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph100_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph100_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph100_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph100_convVtx_Z_barrel = new TH1F("h_ph100_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph100_convVtx_Z_endcap = new TH1F("h_ph100_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph100_convVtx_R_barrel = new TH1F("h_ph100_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph100_convVtx_R_endcap = new TH1F("h_ph100_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph100_convAngle_barrel = new TH1F("h_ph100_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph100_convAngle_endcap = new TH1F("h_ph100_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph100_convE1_over_ETrue_barrel = new TH1F("h_ph100_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph100_convE1_over_ETrue_endcap = new TH1F("h_ph100_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph100_convE2_over_ETrue_barrel = new TH1F("h_ph100_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph100_convE2_over_ETrue_endcap = new TH1F("h_ph100_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph100_conv_DeltaAng_max_barrel = new TH1F("h_ph100_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph100_conv_DeltaAng_max_endcap = new TH1F("h_ph100_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph100_DeltaPh12_barrel = new TH1F("h_ph100_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph100_DeltaPh12_endcap = new TH1F("h_ph100_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph100_DeltaPh12_conv_barrel = new TH1F("h_ph100_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph100_DeltaPh12_conv_endcap = new TH1F("h_ph100_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph100_photonConvRateVsTrueTheta = new TH1F("h_ph100_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph100;
  hist_vec_ph100.push_back(h_ph100_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph100.push_back(h_ph100_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph100.push_back(h_ph100_convVtx_Z_barrel);
  hist_vec_ph100.push_back(h_ph100_convVtx_Z_endcap);
  hist_vec_ph100.push_back(h_ph100_convVtx_R_barrel);
  hist_vec_ph100.push_back(h_ph100_convVtx_R_endcap);
  hist_vec_ph100.push_back(h_ph100_convAngle_barrel);
  hist_vec_ph100.push_back(h_ph100_convAngle_endcap);
  hist_vec_ph100.push_back(h_ph100_convE1_over_ETrue_barrel);
  hist_vec_ph100.push_back(h_ph100_convE1_over_ETrue_endcap);
  hist_vec_ph100.push_back(h_ph100_convE2_over_ETrue_barrel);
  hist_vec_ph100.push_back(h_ph100_convE2_over_ETrue_endcap);
  hist_vec_ph100.push_back(h_ph100_conv_DeltaAng_max_barrel);
  hist_vec_ph100.push_back(h_ph100_conv_DeltaAng_max_endcap);
  hist_vec_ph100.push_back(h_ph100_DeltaPh12_barrel);
  hist_vec_ph100.push_back(h_ph100_DeltaPh12_endcap);
  hist_vec_ph100.push_back(h_ph100_DeltaPh12_conv_barrel);
  hist_vec_ph100.push_back(h_ph100_DeltaPh12_conv_endcap);
  hist_vec_ph100.push_back(h_ph100_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph100.size();i++){
    hist_vec_ph100[i]->Sumw2();
  }
  std::cout<<"fill ph100"<<std::endl;
  fillTreeHistograms(file_photon_ph100,hist_vec_ph100,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph100_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph200_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph200_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph200_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph200_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph200_noOverlay;
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph200_noOverlay);
  TEff_vector_ph200_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph200_noOverlay);

  TH1F* h_ph200_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph200_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph200_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph200_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph200_convVtx_Z_barrel = new TH1F("h_ph200_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph200_convVtx_Z_endcap = new TH1F("h_ph200_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph200_convVtx_R_barrel = new TH1F("h_ph200_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph200_convVtx_R_endcap = new TH1F("h_ph200_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph200_convAngle_barrel = new TH1F("h_ph200_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph200_convAngle_endcap = new TH1F("h_ph200_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph200_convE1_over_ETrue_barrel = new TH1F("h_ph200_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph200_convE1_over_ETrue_endcap = new TH1F("h_ph200_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph200_convE2_over_ETrue_barrel = new TH1F("h_ph200_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph200_convE2_over_ETrue_endcap = new TH1F("h_ph200_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph200_conv_DeltaAng_max_barrel = new TH1F("h_ph200_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph200_conv_DeltaAng_max_endcap = new TH1F("h_ph200_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph200_DeltaPh12_barrel = new TH1F("h_ph200_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph200_DeltaPh12_endcap = new TH1F("h_ph200_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph200_DeltaPh12_conv_barrel = new TH1F("h_ph200_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph200_DeltaPh12_conv_endcap = new TH1F("h_ph200_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph200_photonConvRateVsTrueTheta = new TH1F("h_ph200_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph200;
  hist_vec_ph200.push_back(h_ph200_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph200.push_back(h_ph200_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph200.push_back(h_ph200_convVtx_Z_barrel);
  hist_vec_ph200.push_back(h_ph200_convVtx_Z_endcap);
  hist_vec_ph200.push_back(h_ph200_convVtx_R_barrel);
  hist_vec_ph200.push_back(h_ph200_convVtx_R_endcap);
  hist_vec_ph200.push_back(h_ph200_convAngle_barrel);
  hist_vec_ph200.push_back(h_ph200_convAngle_endcap);
  hist_vec_ph200.push_back(h_ph200_convE1_over_ETrue_barrel);
  hist_vec_ph200.push_back(h_ph200_convE1_over_ETrue_endcap);
  hist_vec_ph200.push_back(h_ph200_convE2_over_ETrue_barrel);
  hist_vec_ph200.push_back(h_ph200_convE2_over_ETrue_endcap);
  hist_vec_ph200.push_back(h_ph200_conv_DeltaAng_max_barrel);
  hist_vec_ph200.push_back(h_ph200_conv_DeltaAng_max_endcap);
  hist_vec_ph200.push_back(h_ph200_DeltaPh12_barrel);
  hist_vec_ph200.push_back(h_ph200_DeltaPh12_endcap);
  hist_vec_ph200.push_back(h_ph200_DeltaPh12_conv_barrel);
  hist_vec_ph200.push_back(h_ph200_DeltaPh12_conv_endcap);
  hist_vec_ph200.push_back(h_ph200_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph200.size();i++){
    hist_vec_ph200[i]->Sumw2();
  }
  std::cout<<"fill ph200"<<std::endl;
  fillTreeHistograms(file_photon_ph200,hist_vec_ph200,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph200_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph500_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph500_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph500_noOverlay;
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph500_noOverlay);
  TEff_vector_ph500_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph500_noOverlay);
  
  TH1F* h_ph500_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph500_convVtx_Z_barrel = new TH1F("h_ph500_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph500_convVtx_Z_endcap = new TH1F("h_ph500_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph500_convVtx_R_barrel = new TH1F("h_ph500_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph500_convVtx_R_endcap = new TH1F("h_ph500_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph500_convAngle_barrel = new TH1F("h_ph500_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph500_convAngle_endcap = new TH1F("h_ph500_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph500_convE1_over_ETrue_barrel = new TH1F("h_ph500_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph500_convE1_over_ETrue_endcap = new TH1F("h_ph500_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph500_convE2_over_ETrue_barrel = new TH1F("h_ph500_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph500_convE2_over_ETrue_endcap = new TH1F("h_ph500_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph500_conv_DeltaAng_max_barrel = new TH1F("h_ph500_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph500_conv_DeltaAng_max_endcap = new TH1F("h_ph500_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph500_DeltaPh12_barrel = new TH1F("h_ph500_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph500_DeltaPh12_endcap = new TH1F("h_ph500_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph500_DeltaPh12_conv_barrel = new TH1F("h_ph500_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph500_DeltaPh12_conv_endcap = new TH1F("h_ph500_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph500_photonConvRateVsTrueTheta = new TH1F("h_ph500_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph500;
  hist_vec_ph500.push_back(h_ph500_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph500.push_back(h_ph500_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph500.push_back(h_ph500_convVtx_Z_barrel);
  hist_vec_ph500.push_back(h_ph500_convVtx_Z_endcap);
  hist_vec_ph500.push_back(h_ph500_convVtx_R_barrel);
  hist_vec_ph500.push_back(h_ph500_convVtx_R_endcap);
  hist_vec_ph500.push_back(h_ph500_convAngle_barrel);
  hist_vec_ph500.push_back(h_ph500_convAngle_endcap);
  hist_vec_ph500.push_back(h_ph500_convE1_over_ETrue_barrel);
  hist_vec_ph500.push_back(h_ph500_convE1_over_ETrue_endcap);
  hist_vec_ph500.push_back(h_ph500_convE2_over_ETrue_barrel);
  hist_vec_ph500.push_back(h_ph500_convE2_over_ETrue_endcap);
  hist_vec_ph500.push_back(h_ph500_conv_DeltaAng_max_barrel);
  hist_vec_ph500.push_back(h_ph500_conv_DeltaAng_max_endcap);
  hist_vec_ph500.push_back(h_ph500_DeltaPh12_barrel);
  hist_vec_ph500.push_back(h_ph500_DeltaPh12_endcap);
  hist_vec_ph500.push_back(h_ph500_DeltaPh12_conv_barrel);
  hist_vec_ph500.push_back(h_ph500_DeltaPh12_conv_endcap);
  hist_vec_ph500.push_back(h_ph500_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph500.size();i++){
    hist_vec_ph500[i]->Sumw2();
  }
  std::cout<<"fill ph500"<<std::endl;
  fillTreeHistograms(file_photon_ph500,hist_vec_ph500,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph500_noOverlay,hist_vec_all_ph);
  
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1000_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1000_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph1000_noOverlay;
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay);
  TEff_vector_ph1000_noOverlay.push_back( tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1000_noOverlay);
    
  TH1F* h_ph1000_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph1000_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph1000_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph1000_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph1000_convVtx_Z_barrel = new TH1F("h_ph1000_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph1000_convVtx_Z_endcap = new TH1F("h_ph1000_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph1000_convVtx_R_barrel = new TH1F("h_ph1000_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph1000_convVtx_R_endcap = new TH1F("h_ph1000_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph1000_convAngle_barrel = new TH1F("h_ph1000_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph1000_convAngle_endcap = new TH1F("h_ph1000_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph1000_convE1_over_ETrue_barrel = new TH1F("h_ph1000_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph1000_convE1_over_ETrue_endcap = new TH1F("h_ph1000_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph1000_convE2_over_ETrue_barrel = new TH1F("h_ph1000_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph1000_convE2_over_ETrue_endcap = new TH1F("h_ph1000_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph1000_conv_DeltaAng_max_barrel = new TH1F("h_ph1000_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph1000_conv_DeltaAng_max_endcap = new TH1F("h_ph1000_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph1000_DeltaPh12_barrel = new TH1F("h_ph1000_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1000_DeltaPh12_endcap = new TH1F("h_ph1000_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1000_DeltaPh12_conv_barrel = new TH1F("h_ph1000_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1000_DeltaPh12_conv_endcap = new TH1F("h_ph1000_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1000_photonConvRateVsTrueTheta = new TH1F("h_ph1000_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph1000;
  hist_vec_ph1000.push_back(h_ph1000_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph1000.push_back(h_ph1000_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph1000.push_back(h_ph1000_convVtx_Z_barrel);
  hist_vec_ph1000.push_back(h_ph1000_convVtx_Z_endcap);
  hist_vec_ph1000.push_back(h_ph1000_convVtx_R_barrel);
  hist_vec_ph1000.push_back(h_ph1000_convVtx_R_endcap);
  hist_vec_ph1000.push_back(h_ph1000_convAngle_barrel);
  hist_vec_ph1000.push_back(h_ph1000_convAngle_endcap);
  hist_vec_ph1000.push_back(h_ph1000_convE1_over_ETrue_barrel);
  hist_vec_ph1000.push_back(h_ph1000_convE1_over_ETrue_endcap);
  hist_vec_ph1000.push_back(h_ph1000_convE2_over_ETrue_barrel);
  hist_vec_ph1000.push_back(h_ph1000_convE2_over_ETrue_endcap);
  hist_vec_ph1000.push_back(h_ph1000_conv_DeltaAng_max_barrel);
  hist_vec_ph1000.push_back(h_ph1000_conv_DeltaAng_max_endcap);
  hist_vec_ph1000.push_back(h_ph1000_DeltaPh12_barrel);
  hist_vec_ph1000.push_back(h_ph1000_DeltaPh12_endcap);
  hist_vec_ph1000.push_back(h_ph1000_DeltaPh12_conv_barrel);
  hist_vec_ph1000.push_back(h_ph1000_DeltaPh12_conv_endcap);
  hist_vec_ph1000.push_back(h_ph1000_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph1000.size();i++){
    hist_vec_ph1000[i]->Sumw2();
  }
  std::cout<<"fill ph1000"<<std::endl;
  fillTreeHistograms(file_photon_ph1000,hist_vec_ph1000,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph1000_noOverlay,hist_vec_all_ph);

  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1500_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1500_noOverlay->SetTitle(";true photon #theta;photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  TEfficiency* tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay = new TEfficiency("tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay->SetTitle(";true photon theta [#circ];photon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_ph1500_noOverlay;
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_noConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_wConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_wConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedConv_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_w_mergedPh_1_deg_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_w_mergedPh_1_deg_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvMergedPh_1_deg_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvMergedPh_1_deg_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay);
  TEff_vector_ph1500_noOverlay.push_back(tEff_PhotonVsTrueTheta_AngMatch_1_deg_EMatch_OnlyConvUnmergedPh_1_deg_ph1500_noOverlay);


  TH1F* h_ph1500_reco_E_over_true_E_match_noOverlay = new TH1F("h_ph1500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph1500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_ph1500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_ph1500_convVtx_Z_barrel = new TH1F("h_ph1500_convVtx_Z_barrel","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph1500_convVtx_Z_endcap = new TH1F("h_ph1500_convVtx_Z_endcap","", n_bins_high,lim_convVtx_Z_min,lim_convVtx_Z_max);
  TH1F* h_ph1500_convVtx_R_barrel = new TH1F("h_ph1500_convVtx_R_barrel","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph1500_convVtx_R_endcap = new TH1F("h_ph1500_convVtx_R_endcap","", n_bins_high,lim_convVtx_R_min,lim_convVtx_R_max);
  TH1F* h_ph1500_convAngle_barrel = new TH1F("h_ph1500_convAngle_barrel","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph1500_convAngle_endcap = new TH1F("h_ph1500_convAngle_endcap","", n_bins_high,lim_convAngle_min,lim_convAngle_max);
  TH1F* h_ph1500_convE1_over_ETrue_barrel = new TH1F("h_ph1500_convE1_over_ETrue_barrel","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph1500_convE1_over_ETrue_endcap = new TH1F("h_ph1500_convE1_over_ETrue_endcap","", n_bins_high,lim_convE1rel_min,lim_convE1rel_max);
  TH1F* h_ph1500_convE2_over_ETrue_barrel = new TH1F("h_ph1500_convE2_over_ETrue_barrel","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph1500_convE2_over_ETrue_endcap = new TH1F("h_ph1500_convE2_over_ETrue_endcap","", n_bins_high,lim_convE2rel_min,lim_convE2rel_max);
  TH1F* h_ph1500_conv_DeltaAng_max_barrel = new TH1F("h_ph1500_conv_DeltaAng_barrel","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph1500_conv_DeltaAng_max_endcap = new TH1F("h_ph1500_conv_DeltaAng_endcap","", n_bins_high,lim_convDelta_min,lim_convDelta_max);
  TH1F* h_ph1500_DeltaPh12_barrel = new TH1F("h_ph1500_DeltaPh12_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1500_DeltaPh12_endcap = new TH1F("h_ph1500_DeltaPh12_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1500_DeltaPh12_conv_barrel = new TH1F("h_ph1500_DeltaPh12_conv_barrel","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1500_DeltaPh12_conv_endcap = new TH1F("h_ph1500_DeltaPh12_conv_endcap","", n_bins_high,lim_recoAngle_min,lim_recoAngle_max);
  TH1F* h_ph1500_photonConvRateVsTrueTheta = new TH1F("h_ph1500_PhotonConvRateVsTrueTheta","", n_bins,lim_theta_low,lim_theta_high);
  std::vector<TH1F*>hist_vec_ph1500;
  hist_vec_ph1500.push_back(h_ph1500_reco_E_over_true_E_match_noOverlay);
  hist_vec_ph1500.push_back(h_ph1500_reco_pt_over_true_pt_match_noOverlay);
  hist_vec_ph1500.push_back(h_ph1500_convVtx_Z_barrel);
  hist_vec_ph1500.push_back(h_ph1500_convVtx_Z_endcap);
  hist_vec_ph1500.push_back(h_ph1500_convVtx_R_barrel);
  hist_vec_ph1500.push_back(h_ph1500_convVtx_R_endcap);
  hist_vec_ph1500.push_back(h_ph1500_convAngle_barrel);
  hist_vec_ph1500.push_back(h_ph1500_convAngle_endcap);
  hist_vec_ph1500.push_back(h_ph1500_convE1_over_ETrue_barrel);
  hist_vec_ph1500.push_back(h_ph1500_convE1_over_ETrue_endcap);
  hist_vec_ph1500.push_back(h_ph1500_convE2_over_ETrue_barrel);
  hist_vec_ph1500.push_back(h_ph1500_convE2_over_ETrue_endcap);
  hist_vec_ph1500.push_back(h_ph1500_conv_DeltaAng_max_barrel);
  hist_vec_ph1500.push_back(h_ph1500_conv_DeltaAng_max_endcap);
  hist_vec_ph1500.push_back(h_ph1500_DeltaPh12_barrel);
  hist_vec_ph1500.push_back(h_ph1500_DeltaPh12_endcap);
  hist_vec_ph1500.push_back(h_ph1500_DeltaPh12_conv_barrel);
  hist_vec_ph1500.push_back(h_ph1500_DeltaPh12_conv_endcap);
  hist_vec_ph1500.push_back(h_ph1500_photonConvRateVsTrueTheta);
  for(unsigned int i=0;i<hist_vec_ph1500.size();i++){
    hist_vec_ph1500[i]->Sumw2();
  }
  std::cout<<"fill ph1500"<<std::endl;
  fillTreeHistograms(file_photon_ph1500,hist_vec_ph1500,test_signal_particle,TEff_vector_all_ph_noOverlay,TEff_vector_ph1500_noOverlay,hist_vec_all_ph);



  file_histogram_photons->cd();
  gre_ph_E_reco_over_E_true_0_78->Write();
  gre_ph_E_reco_over_E_true_0_78_to_0_83->Write();
  gre_ph_E_reco_over_E_true_0_83_to_0_94->Write();
  file_histogram_photons->Write();
  file_histogram_photons->Close();



  TFile* file_K0L5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L5_ILC181011_CT.root");
  TFile* file_K0L10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L10_ILC181011_CT.root");
  TFile* file_K0L20=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L20_ILC181011_CT.root");
  TFile* file_K0L30=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L30_ILC181011_CT.root");
  TFile* file_K0L40=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L40_ILC181011_CT.root");
  TFile* file_K0L50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L50_ILC181011_CT.root");
  TFile* file_K0L60=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L60_ILC181011_CT.root");
  TFile* file_K0L75=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L75_ILC181011_CT.root");
  TFile* file_K0L90=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L90_ILC181011_CT.root");
  TFile* file_K0L100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L100_ILC181011_CT.root");
  TFile* file_K0L150=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L150_ILC181011_CT.root");
  TFile* file_K0L200=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L200_ILC181011_CT.root");
  TFile* file_K0L250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L250_ILC181011_CT.root");
  TFile* file_K0L400=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181011_gcc62/pionStudy_K0L400_ILC181011_CT.root");
 
  //for K0L events 5/10/20/30/40/50/60/75/90/100/150/200/250/400
  const char* final_histo_name_K0L="/eos/user/w/weberma2/data/validation181011/pionStudy_K0L5_to_K0L400_ILC181011_finalhistos_deg2_match.root";

  TFile* file_histogram_K0L=new TFile(final_histo_name_K0L,"recreate");
    
  n_bins_highRES=100;
  n_bins_highTRRES=75;
  lim_rel_energy_lowRES=0.00;
  lim_rel_energy_highRES=2.00;

  test_signal_particle_true = 130;
  test_signal_particle_reco = 2112;

  TGraphErrors* gre_K0L_E_reco_over_E_true_0_65= new TGraphErrors(14);
  gre_K0L_E_reco_over_E_true_0_65->SetName("gre_K0L_E_reco_over_E_true_0_65");
  gre_K0L_E_reco_over_E_true_0_65->SetLineColor(kBlack); 
  gre_K0L_E_reco_over_E_true_0_65->SetMarkerColor(kBlack); 
  gre_K0L_E_reco_over_E_true_0_65->SetMarkerStyle(kFullCircle); 
  gre_K0L_E_reco_over_E_true_0_65->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_K0L_E_reco_over_E_true_0_65->GetXaxis()->SetTitle("E(K^{0}_{L,true}) [GeV]");

  TGraphErrors* gre_K0L_E_reco_over_E_true_0_65_to_0_80 = new TGraphErrors(14);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetName("gre_K0L_E_reco_over_E_true_0_65_to_0_80");
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetLineColor(kRed); 
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetMarkerColor(kRed); 
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetMarkerStyle(kFullSquare); 
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->GetXaxis()->SetTitle("E(K^{0}_{L,true}) [GeV]");

  TGraphErrors* gre_K0L_E_reco_over_E_true_0_80_to_0_94 = new TGraphErrors(14);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetName("gre_K0L_E_reco_over_E_true_0_80_to_0_94");
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetLineColor(kBlue); 
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetMarkerColor(kBlue); 
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetMarkerStyle(kFullTriangleUp); 
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->GetYaxis()->SetTitle("#sigma(E_{reco}/E_{true}) [%]");
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->GetXaxis()->SetTitle("E(K^{0}_{L,true}) [GeV]");


  std::cout<<"RES K0L5"<<std::endl;
  TH1F* h_K0L5_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L5_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L5_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L5_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L5_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L5_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L5;
  TH_vec_ERES_K0L5.push_back(h_K0L5_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L5.push_back(h_K0L5_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L5.push_back(h_K0L5_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L5.size();i++){
    TH_vec_ERES_K0L5[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L5, TH_vec_ERES_K0L5,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L5_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L5[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L5_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L5[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L5_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L5[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(0,5,100.*sigma_K0L5_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(0,0,100.*sigma_K0L5_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(0,5,100.*sigma_K0L5_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(0,0,100.*sigma_K0L5_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(0,5,100.*sigma_K0L5_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(0,0,100.*sigma_K0L5_0_80_to_0_94[1]);

  lim_rel_energy_lowRES=0.50;
  lim_rel_energy_highRES=1.50;

  std::cout<<"RES K0L10"<<std::endl;
  TH1F* h_K0L10_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L10_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L10_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L10_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L10_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L10_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L10;
  TH_vec_ERES_K0L10.push_back(h_K0L10_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L10.push_back(h_K0L10_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L10.push_back(h_K0L10_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L10.size();i++){
    TH_vec_ERES_K0L10[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L10, TH_vec_ERES_K0L10,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L10_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L10[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L10_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L10[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L10_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L10[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(1,10,100.*sigma_K0L10_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(1,0,100.*sigma_K0L10_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(1,10,100.*sigma_K0L10_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(1,0,100.*sigma_K0L10_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(1,10,100.*sigma_K0L10_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(1,0,100.*sigma_K0L10_0_80_to_0_94[1]);

  std::cout<<"RES K0L20"<<std::endl;
  TH1F* h_K0L20_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L20_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L20_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L20_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L20_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L20_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L20;
  TH_vec_ERES_K0L20.push_back(h_K0L20_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L20.push_back(h_K0L20_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L20.push_back(h_K0L20_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L20.size();i++){
    TH_vec_ERES_K0L20[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L20, TH_vec_ERES_K0L20,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L20_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L20[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L20_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L20[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L20_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L20[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(2,20,100.*sigma_K0L20_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(2,0,100.*sigma_K0L20_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(2,20,100.*sigma_K0L20_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(2,0,100.*sigma_K0L20_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(2,20,100.*sigma_K0L20_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(2,0,100.*sigma_K0L20_0_80_to_0_94[1]);

  std::cout<<"RES K0L30"<<std::endl;
  TH1F* h_K0L30_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L30_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L30_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L30_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L30_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L30_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L30;
  TH_vec_ERES_K0L30.push_back(h_K0L30_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L30.push_back(h_K0L30_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L30.push_back(h_K0L30_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L30.size();i++){
    TH_vec_ERES_K0L30[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L30, TH_vec_ERES_K0L30,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L30_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L30[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L30_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L30[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L30_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L30[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(3,30,100.*sigma_K0L30_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(3,0,100.*sigma_K0L30_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(3,30,100.*sigma_K0L30_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(3,0,100.*sigma_K0L30_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(3,30,100.*sigma_K0L30_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(3,0,100.*sigma_K0L30_0_80_to_0_94[1]);

  std::cout<<"RES K0L40"<<std::endl;
  TH1F* h_K0L40_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L40_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L40_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L40_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L40_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L40_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L40;
  TH_vec_ERES_K0L40.push_back(h_K0L40_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L40.push_back(h_K0L40_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L40.push_back(h_K0L40_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L40.size();i++){
    TH_vec_ERES_K0L40[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L40, TH_vec_ERES_K0L40,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L40_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L40[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L40_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L40[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L40_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L40[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(4,40,100.*sigma_K0L40_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(4,0,100.*sigma_K0L40_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(4,40,100.*sigma_K0L40_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(4,0,100.*sigma_K0L40_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(4,40,100.*sigma_K0L40_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(4,0,100.*sigma_K0L40_0_80_to_0_94[1]);

  std::cout<<"RES K0L50"<<std::endl;
  TH1F* h_K0L50_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L50_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L50_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L50_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L50_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L50_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L50;
  TH_vec_ERES_K0L50.push_back(h_K0L50_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L50.push_back(h_K0L50_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L50.push_back(h_K0L50_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L50.size();i++){
    TH_vec_ERES_K0L50[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L50, TH_vec_ERES_K0L50,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L50_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L50[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L50_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L50[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L50_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L50[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(5,50,100.*sigma_K0L50_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(5,0,100.*sigma_K0L50_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(5,50,100.*sigma_K0L50_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(5,0,100.*sigma_K0L50_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(5,50,100.*sigma_K0L50_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(5,0,100.*sigma_K0L50_0_80_to_0_94[1]);

  lim_rel_energy_lowRES=0.75;
  lim_rel_energy_highRES=1.25;

 std::cout<<"RES K0L60"<<std::endl;
  TH1F* h_K0L60_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L60_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L60_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L60_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L60_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L60_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L60;
  TH_vec_ERES_K0L60.push_back(h_K0L60_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L60.push_back(h_K0L60_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L60.push_back(h_K0L60_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L60.size();i++){
    TH_vec_ERES_K0L60[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L60, TH_vec_ERES_K0L60,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L60_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L60[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L60_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L60[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L60_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L60[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(6,60,100.*sigma_K0L60_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(6,0,100.*sigma_K0L60_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(6,60,100.*sigma_K0L60_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(6,0,100.*sigma_K0L60_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(6,60,100.*sigma_K0L60_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(6,0,100.*sigma_K0L60_0_80_to_0_94[1]);

  std::cout<<"RES K0L75"<<std::endl;
  TH1F* h_K0L75_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L75_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L75_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L75_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L75_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L75_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L75;
  TH_vec_ERES_K0L75.push_back(h_K0L75_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L75.push_back(h_K0L75_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L75.push_back(h_K0L75_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L75.size();i++){
    TH_vec_ERES_K0L75[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L75, TH_vec_ERES_K0L75,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L75_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L75[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L75_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L75[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L75_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L75[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(7,75,100.*sigma_K0L75_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(7,0,100.*sigma_K0L75_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(7,75,100.*sigma_K0L75_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(7,0,100.*sigma_K0L75_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(7,75,100.*sigma_K0L75_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(7,0,100.*sigma_K0L75_0_80_to_0_94[1]);

  std::cout<<"RES K0L90"<<std::endl;
  TH1F* h_K0L90_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L90_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L90_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L90_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L90_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L90_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L90;
  TH_vec_ERES_K0L90.push_back(h_K0L90_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L90.push_back(h_K0L90_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L90.push_back(h_K0L90_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L90.size();i++){
    TH_vec_ERES_K0L90[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L90, TH_vec_ERES_K0L90,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L90_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L90[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L90_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L90[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L90_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L90[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(8,90,100.*sigma_K0L90_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(8,0,100.*sigma_K0L90_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(8,90,100.*sigma_K0L90_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(8,0,100.*sigma_K0L90_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(8,90,100.*sigma_K0L90_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(8,0,100.*sigma_K0L90_0_80_to_0_94[1]);

  std::cout<<"RES K0L100"<<std::endl;
  TH1F* h_K0L100_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L100_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L100_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L100_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L100_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L100_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L100;
  TH_vec_ERES_K0L100.push_back(h_K0L100_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L100.push_back(h_K0L100_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L100.push_back(h_K0L100_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L100.size();i++){
    TH_vec_ERES_K0L100[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L100, TH_vec_ERES_K0L100,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L100_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L100[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L100_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L100[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L100_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L100[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(9,100,100.*sigma_K0L100_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(9,0,100.*sigma_K0L100_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(9,100,100.*sigma_K0L100_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(9,0,100.*sigma_K0L100_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(9,100,100.*sigma_K0L100_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(9,0,100.*sigma_K0L100_0_80_to_0_94[1]);

  std::cout<<"RES K0L150"<<std::endl;
  TH1F* h_K0L150_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L150_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L150_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L150_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L150_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L150_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L150;
  TH_vec_ERES_K0L150.push_back(h_K0L150_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L150.push_back(h_K0L150_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L150.push_back(h_K0L150_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L150.size();i++){
    TH_vec_ERES_K0L150[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L150, TH_vec_ERES_K0L150,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L150_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L150[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L150_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L150[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L150_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L150[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(10,150,100.*sigma_K0L150_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(10,0,100.*sigma_K0L150_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(10,150,100.*sigma_K0L150_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(10,0,100.*sigma_K0L150_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(10,150,100.*sigma_K0L150_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(10,0,100.*sigma_K0L150_0_80_to_0_94[1]);


  std::cout<<"RES K0L200"<<std::endl;
  TH1F* h_K0L200_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L200_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L200_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L200_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L200_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L200_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L200;
  TH_vec_ERES_K0L200.push_back(h_K0L200_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L200.push_back(h_K0L200_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L200.push_back(h_K0L200_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L200.size();i++){
    TH_vec_ERES_K0L200[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L200, TH_vec_ERES_K0L200,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L200_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L200[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L200_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L200[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L200_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L200[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(11,200,100.*sigma_K0L200_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(11,0,100.*sigma_K0L200_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(11,200,100.*sigma_K0L200_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(11,0,100.*sigma_K0L200_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(11,200,100.*sigma_K0L200_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(11,0,100.*sigma_K0L200_0_80_to_0_94[1]);

  std::cout<<"RES K0L250"<<std::endl;
  TH1F* h_K0L250_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L250_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L250_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L250_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L250_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L250_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L250;
  TH_vec_ERES_K0L250.push_back(h_K0L250_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L250.push_back(h_K0L250_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L250.push_back(h_K0L250_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L250.size();i++){
    TH_vec_ERES_K0L250[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L250, TH_vec_ERES_K0L250,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L250_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L250[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L250_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L250[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L250_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L250[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(12,250,100.*sigma_K0L250_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(12,0,100.*sigma_K0L250_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(12,250,100.*sigma_K0L250_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(12,0,100.*sigma_K0L250_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(12,250,100.*sigma_K0L250_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(12,0,100.*sigma_K0L250_0_80_to_0_94[1]);

  std::cout<<"RES K0L400"<<std::endl;
  TH1F* h_K0L400_reco_E_over_true_E_match_0_65 = new TH1F("h_K0L400_reco_E_over_true_E_match_0_65","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L400_reco_E_over_true_E_match_0_65_to_0_80 = new TH1F("h_K0L400_reco_E_over_true_E_match_0_65_to_0_80","", n_bins_highTRRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  TH1F* h_K0L400_reco_E_over_true_E_match_0_80_to_0_94 = new TH1F("h_K0L400_reco_E_over_true_E_match_0_80_to_0_94","", n_bins_highRES,lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<TH1F*> TH_vec_ERES_K0L400;
  TH_vec_ERES_K0L400.push_back(h_K0L400_reco_E_over_true_E_match_0_65);
  TH_vec_ERES_K0L400.push_back(h_K0L400_reco_E_over_true_E_match_0_65_to_0_80);
  TH_vec_ERES_K0L400.push_back(h_K0L400_reco_E_over_true_E_match_0_80_to_0_94);
  for(unsigned int i=0;i<TH_vec_ERES_K0L400.size();i++){
    TH_vec_ERES_K0L400[i]->Sumw2();
  }
  fillRESTreeHistograms(file_K0L400, TH_vec_ERES_K0L400,test_signal_particle_true,test_signal_particle_reco);
  std::vector<float>sigma_K0L400_0_65=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L400[0],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L400_0_65_to_0_80=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L400[1],lim_rel_energy_lowRES,lim_rel_energy_highRES);
  std::vector<float>sigma_K0L400_0_80_to_0_94=sigmaGaussianFits2Sigma(TH_vec_ERES_K0L400[2],lim_rel_energy_lowRES,lim_rel_energy_highRES);

  gre_K0L_E_reco_over_E_true_0_65->SetPoint(13,400,100.*sigma_K0L400_0_65[0]);
  gre_K0L_E_reco_over_E_true_0_65->SetPointError(13,0,100.*sigma_K0L400_0_65[1]);

  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPoint(13,400,100.*sigma_K0L400_0_65_to_0_80[0]);
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->SetPointError(13,0,100.*sigma_K0L400_0_65_to_0_80[1]);

  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPoint(13,400,100.*sigma_K0L400_0_80_to_0_94[0]);
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->SetPointError(13,0,100.*sigma_K0L400_0_80_to_0_94[1]);

  file_histogram_K0L->cd();
  gre_K0L_E_reco_over_E_true_0_65->Write();
  gre_K0L_E_reco_over_E_true_0_65_to_0_80->Write();
  gre_K0L_E_reco_over_E_true_0_80_to_0_94->Write();
  file_histogram_K0L->Write();
  file_histogram_K0L->Close();

  TFile* file_electron_em5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em5_ILC181101_CT.root");
  TFile* file_electron_em10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em10_ILC181101_CT.root");
  TFile* file_electron_em15=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em15_ILC181101_CT.root");
  TFile* file_electron_em20=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em20_ILC181101_CT.root");
  TFile* file_electron_em25=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em25_ILC181101_CT.root");
  TFile* file_electron_em40=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em40_ILC181101_CT.root");
  TFile* file_electron_em50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em50_ILC181101_CT.root");
  TFile* file_electron_em100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em100_ILC181101_CT.root");
  TFile* file_electron_em250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em250_ILC181101_CT.root");
  TFile* file_electron_em500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em500_ILC181101_CT.root");
  TFile* file_electron_em750=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em750_ILC181101_CT.root");
  TFile* file_electron_em1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em1000_ILC181101_CT.root");
  TFile* file_electron_em1500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_em1500_ILC181101_CT.root");

  test_signal_particle = 11;
  const unsigned int nbins_teff_em = 13; 
  double xbins_teff_em[nbins_teff_em+1]={2.5,7.5,12.5,17.5,22.5,30.,45.,75.,150.,300.,600.,800.,1250.,1550.};

  const char* final_histo_name_electrons="/eos/user/w/weberma2/data/validation181101/pionStudy_em5_to_em1500_ILC181101_finalhistos_deg1_match_EMatch_5.root";
  TFile* file_histogram_electrons=new TFile(final_histo_name_electrons,"recreate");

  TEfficiency* tEff_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_em,xbins_teff_em);
  tEff_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");
  TEfficiency* tEff_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_em,xbins_teff_em);
  tEff_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");

  std::vector<TEfficiency*> TEff_vector_all_em_noOverlay;
  TEff_vector_all_em_noOverlay.push_back(tEff_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_em_noOverlay.push_back(tEff_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em1_noOverlay;
  TEff_vector_em1_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay);
  TEff_vector_em1_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay);
    
  TH1F* h_em1_reco_E_over_true_E_match_noOverlay = new TH1F("h_em1_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em1_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em1_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em1;
  hist_vec_em1.push_back(h_em1_reco_E_over_true_E_match_noOverlay);
  hist_vec_em1.push_back(h_em1_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em1.size();i++){
    hist_vec_em1[i]->Sumw2();
  }

  std::cout<<"em5"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em5_noOverlay;
  TEff_vector_em5_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay);
  TEff_vector_em5_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay);
    
  TH1F* h_em5_reco_E_over_true_E_match_noOverlay = new TH1F("h_em5_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em5_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em5;
  hist_vec_em5.push_back(h_em5_reco_E_over_true_E_match_noOverlay);
  hist_vec_em5.push_back(h_em5_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em5.size();i++){
    hist_vec_em5[i]->Sumw2();
  }
  std::vector<TH1F*>hist_vec_all_em;


  fillTreeHistograms(file_electron_em5,hist_vec_em5,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em5_noOverlay,hist_vec_all_em);
  std::cout<<"em10"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em10_noOverlay;
  TEff_vector_em10_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay);
  TEff_vector_em10_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay);
  
  TH1F* h_em10_reco_E_over_true_E_match_noOverlay = new TH1F("h_em10_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em10_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em10_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em10;
  hist_vec_em10.push_back(h_em10_reco_E_over_true_E_match_noOverlay);
  hist_vec_em10.push_back(h_em10_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em10.size();i++){
    hist_vec_em10[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em10,hist_vec_em10,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em10_noOverlay,hist_vec_all_em);
  std::cout<<"em15"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em15_noOverlay;
  TEff_vector_em15_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay);
  TEff_vector_em15_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay);
    
  TH1F* h_em15_reco_E_over_true_E_match_noOverlay = new TH1F("h_em15_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em15_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em15_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em15;
  hist_vec_em15.push_back(h_em15_reco_E_over_true_E_match_noOverlay);
  hist_vec_em15.push_back(h_em15_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em15.size();i++){
    hist_vec_em15[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em15,hist_vec_em15,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em15_noOverlay,hist_vec_all_em);
  std::cout<<"em20"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em20_noOverlay;
  TEff_vector_em20_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay);
  TEff_vector_em20_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay);
  
  TH1F* h_em20_reco_E_over_true_E_match_noOverlay = new TH1F("h_em20_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em20_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em20_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em20;
  hist_vec_em20.push_back(h_em20_reco_E_over_true_E_match_noOverlay);
  hist_vec_em20.push_back(h_em20_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em20.size();i++){
    hist_vec_em20[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em20,hist_vec_em20,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em20_noOverlay,hist_vec_all_em);
  std::cout<<"em25"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em25_noOverlay;
  TEff_vector_em25_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay);
  TEff_vector_em25_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay);
  
  TH1F* h_em25_reco_E_over_true_E_match_noOverlay = new TH1F("h_em25_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em25_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em25_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em25;
  hist_vec_em25.push_back(h_em25_reco_E_over_true_E_match_noOverlay);
  hist_vec_em25.push_back(h_em25_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em25.size();i++){
    hist_vec_em25[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em25,hist_vec_em25,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em25_noOverlay,hist_vec_all_em);
  std::cout<<"em40"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em40_noOverlay;
  TEff_vector_em40_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay);
  TEff_vector_em40_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay);
  
  TH1F* h_em40_reco_E_over_true_E_match_noOverlay = new TH1F("h_em40_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em40_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em40_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em40;
  hist_vec_em40.push_back(h_em40_reco_E_over_true_E_match_noOverlay);
  hist_vec_em40.push_back(h_em40_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em40.size();i++){
    hist_vec_em40[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em40,hist_vec_em40,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em40_noOverlay,hist_vec_all_em);
  std::cout<<"em50"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em50_noOverlay;
  TEff_vector_em50_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay);
  TEff_vector_em50_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay);
    
  TH1F* h_em50_reco_E_over_true_E_match_noOverlay = new TH1F("h_em50_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em50_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em50_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em50;
  hist_vec_em50.push_back(h_em50_reco_E_over_true_E_match_noOverlay);
  hist_vec_em50.push_back(h_em50_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em50.size();i++){
    hist_vec_em50[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em50,hist_vec_em50,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em50_noOverlay,hist_vec_all_em);
  std::cout<<"em100"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em100_noOverlay;
  TEff_vector_em100_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay);
  TEff_vector_em100_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay);
    
  TH1F* h_em100_reco_E_over_true_E_match_noOverlay = new TH1F("h_em100_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em100_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em100_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em100;
  hist_vec_em100.push_back(h_em100_reco_E_over_true_E_match_noOverlay);
  hist_vec_em100.push_back(h_em100_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em100.size();i++){
    hist_vec_em100[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em100,hist_vec_em100,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em100_noOverlay,hist_vec_all_em);
  std::cout<<"em250"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em250_noOverlay;
  TEff_vector_em250_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay);
  TEff_vector_em250_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay);
    
  TH1F* h_em250_reco_E_over_true_E_match_noOverlay = new TH1F("h_em250_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em250_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em250_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em250;
  hist_vec_em250.push_back(h_em250_reco_E_over_true_E_match_noOverlay);
  hist_vec_em250.push_back(h_em250_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em250.size();i++){
    hist_vec_em250[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em250,hist_vec_em250,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em250_noOverlay,hist_vec_all_em);
  std::cout<<"em500"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em500_noOverlay;
  TEff_vector_em500_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay);
  TEff_vector_em500_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay);
  
  TH1F* h_em500_reco_E_over_true_E_match_noOverlay = new TH1F("h_em500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em500;
  hist_vec_em500.push_back(h_em500_reco_E_over_true_E_match_noOverlay);
  hist_vec_em500.push_back(h_em500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em500.size();i++){
    hist_vec_em500[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em500,hist_vec_em500,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em500_noOverlay,hist_vec_all_em);
  std::cout<<"em750"<<std::endl;
  
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em750_noOverlay;
  TEff_vector_em750_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay);
  TEff_vector_em750_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay);
  
  TH1F* h_em750_reco_E_over_true_E_match_noOverlay = new TH1F("h_em750_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em750_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em750_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em750;
  hist_vec_em750.push_back(h_em750_reco_E_over_true_E_match_noOverlay);
  hist_vec_em750.push_back(h_em750_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em750.size();i++){
    hist_vec_em750[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em750,hist_vec_em750,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em750_noOverlay,hist_vec_all_em);
  std::cout<<"em1000"<<std::endl;
    
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em1000_noOverlay;
  TEff_vector_em1000_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay);
  TEff_vector_em1000_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay);
  
  TH1F* h_em1000_reco_E_over_true_E_match_noOverlay = new TH1F("h_em1000_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em1000_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em1000_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em1000;
  hist_vec_em1000.push_back(h_em1000_reco_E_over_true_E_match_noOverlay);
  hist_vec_em1000.push_back(h_em1000_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em1000.size();i++){
    hist_vec_em1000[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em1000,hist_vec_em1000,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em1000_noOverlay,hist_vec_all_em);
  std::cout<<"em1500"<<std::endl;
  
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay = new TEfficiency("tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_em1500_noOverlay;
  TEff_vector_em1500_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay);
  TEff_vector_em1500_noOverlay.push_back(tEff_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay);
  
  TH1F* h_em1500_reco_E_over_true_E_match_noOverlay = new TH1F("h_em1500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_em1500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_em1500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_em1500;
  hist_vec_em1500.push_back(h_em1500_reco_E_over_true_E_match_noOverlay);
  hist_vec_em1500.push_back(h_em1500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_em1500.size();i++){
    hist_vec_em1500[i]->Sumw2();
  }

  fillTreeHistograms(file_electron_em1500,hist_vec_em1500,test_signal_particle,TEff_vector_all_em_noOverlay,TEff_vector_em1500_noOverlay,hist_vec_all_em);


  file_histogram_electrons->Write();
  file_histogram_electrons->Close();
  

  TFile* file_dressed_electron_em5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em5_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em10_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em15=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em15_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em20=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em20_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em25=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em25_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em40=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em40_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em50_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em100_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em250_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em500_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em750=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em750_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em1000_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");
  TFile* file_dressed_electron_em1500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180208/pionStudy_em1500_ILC180208_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_DressedLepton.root");

  const char* final_histo_name_dressed_electrons="/eos/user/w/weberma2/data/validation171212/pionStudy_em5_to_em1500_ILC180208_CT_FitFW_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune_default_finalhistos_deg1_match_DressedTest_EMatch_5.root";
  TFile* file_histogram_dressed_electrons=new TFile(final_histo_name_dressed_electrons,"recreate");

  TEfficiency* tEff_dressed_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_Dressed_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_em,xbins_teff_em);
  tEff_dressed_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_Dressed_ElectronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_em,xbins_teff_em);
  tEff_dressed_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");

  std::vector<TEfficiency*> TEff_vector_all_dressed_em_noOverlay;
  TEff_vector_all_dressed_em_noOverlay.push_back(tEff_dressed_electronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_dressed_em_noOverlay.push_back(tEff_dressed_electronVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em1_noOverlay;
  TEff_vector_dressed_em1_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1_noOverlay);
  TEff_vector_dressed_em1_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1_noOverlay);
    
  TH1F* h_dressed_em1_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em1_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em1_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em1_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em1;
  hist_vec_dressed_em1.push_back(h_dressed_em1_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em1.push_back(h_dressed_em1_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em1.size();i++){
    hist_vec_dressed_em1[i]->Sumw2();
  }
  std::vector<TH1F*>hist_dressed_vec_all_em;

  std::cout<<"em5"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em5_noOverlay;
  TEff_vector_dressed_em5_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em5_noOverlay);
  TEff_vector_dressed_em5_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em5_noOverlay);
    
  TH1F* h_dressed_em5_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em5_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em5_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em5;
  hist_vec_dressed_em5.push_back(h_dressed_em5_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em5.push_back(h_dressed_em5_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em5.size();i++){
    hist_vec_dressed_em5[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em5,hist_vec_dressed_em5,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em5_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em10"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em10_noOverlay;
  TEff_vector_dressed_em10_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em10_noOverlay);
  TEff_vector_dressed_em10_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em10_noOverlay);
  
  TH1F* h_dressed_em10_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em10_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em10_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em10_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em10;
  hist_vec_dressed_em10.push_back(h_dressed_em10_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em10.push_back(h_dressed_em10_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em10.size();i++){
    hist_vec_dressed_em10[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em10,hist_vec_dressed_em10,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em10_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em15"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em15_noOverlay;
  TEff_vector_dressed_em15_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em15_noOverlay);
  TEff_vector_dressed_em15_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em15_noOverlay);
    
  TH1F* h_dressed_em15_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em15_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em15_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em15_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em15;
  hist_vec_dressed_em15.push_back(h_dressed_em15_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em15.push_back(h_dressed_em15_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em15.size();i++){
    hist_vec_dressed_em15[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em15,hist_vec_dressed_em15,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em15_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em20"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em20_noOverlay;
  TEff_vector_dressed_em20_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em20_noOverlay);
  TEff_vector_dressed_em20_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em20_noOverlay);
  
  TH1F* h_dressed_em20_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em20_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em20_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em20_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em20;
  hist_vec_dressed_em20.push_back(h_dressed_em20_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em20.push_back(h_dressed_em20_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em20.size();i++){
    hist_vec_dressed_em20[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em20,hist_vec_dressed_em20,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em20_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em25"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em25_noOverlay;
  TEff_vector_dressed_em25_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em25_noOverlay);
  TEff_vector_dressed_em25_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em25_noOverlay);
  
  TH1F* h_dressed_em25_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em25_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em25_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em25_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em25;
  hist_vec_dressed_em25.push_back(h_dressed_em25_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em25.push_back(h_dressed_em25_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em25.size();i++){
    hist_vec_dressed_em25[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em25,hist_vec_dressed_em25,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em25_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em40"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em40_noOverlay;
  TEff_vector_dressed_em40_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em40_noOverlay);
  TEff_vector_dressed_em40_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em40_noOverlay);
  
  TH1F* h_dressed_em40_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em40_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em40_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em40_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em40;
  hist_vec_dressed_em40.push_back(h_dressed_em40_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em40.push_back(h_dressed_em40_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em40.size();i++){
    hist_vec_dressed_em40[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em40,hist_vec_dressed_em40,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em40_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em50"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em50_noOverlay;
  TEff_vector_dressed_em50_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em50_noOverlay);
  TEff_vector_dressed_em50_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em50_noOverlay);
    
  TH1F* h_dressed_em50_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em50_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em50_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em50_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em50;
  hist_vec_dressed_em50.push_back(h_dressed_em50_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em50.push_back(h_dressed_em50_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em50.size();i++){
    hist_vec_dressed_em50[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em50,hist_vec_dressed_em50,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em50_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em100"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em100_noOverlay;
  TEff_vector_dressed_em100_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em100_noOverlay);
  TEff_vector_dressed_em100_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em100_noOverlay);
    
  TH1F* h_dressed_em100_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em100_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em100_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em100_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em100;
  hist_vec_dressed_em100.push_back(h_dressed_em100_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em100.push_back(h_dressed_em100_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em100.size();i++){
    hist_vec_dressed_em100[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em100,hist_vec_dressed_em100,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em100_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em250"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em250_noOverlay;
  TEff_vector_dressed_em250_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em250_noOverlay);
  TEff_vector_dressed_em250_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em250_noOverlay);
    
  TH1F* h_dressed_em250_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em250_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em250_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em250_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em250;
  hist_vec_dressed_em250.push_back(h_dressed_em250_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em250.push_back(h_dressed_em250_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em250.size();i++){
    hist_vec_dressed_em250[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em250,hist_vec_dressed_em250,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em250_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em500"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em500_noOverlay;
  TEff_vector_dressed_em500_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em500_noOverlay);
  TEff_vector_dressed_em500_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em500_noOverlay);
  
  TH1F* h_dressed_em500_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em500;
  hist_vec_dressed_em500.push_back(h_dressed_em500_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em500.push_back(h_dressed_em500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em500.size();i++){
    hist_vec_dressed_em500[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em500,hist_vec_dressed_em500,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em500_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em750"<<std::endl;
  
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em750_noOverlay;
  TEff_vector_dressed_em750_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em750_noOverlay);
  TEff_vector_dressed_em750_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em750_noOverlay);
  
  TH1F* h_dressed_em750_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em750_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em750_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em750_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em750;
  hist_vec_dressed_em750.push_back(h_dressed_em750_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em750.push_back(h_dressed_em750_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em750.size();i++){
    hist_vec_dressed_em750[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em750,hist_vec_dressed_em750,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em750_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em1000"<<std::endl;
    
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em1000_noOverlay;
  TEff_vector_dressed_em1000_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1000_noOverlay);
  TEff_vector_dressed_em1000_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1000_noOverlay);
  
  TH1F* h_dressed_em1000_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em1000_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em1000_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em1000_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em1000;
  hist_vec_dressed_em1000.push_back(h_dressed_em1000_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em1000.push_back(h_dressed_em1000_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em1000.size();i++){
    hist_vec_dressed_em1000[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em1000,hist_vec_dressed_em1000,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em1000_noOverlay,hist_dressed_vec_all_em);
  std::cout<<"em1500"<<std::endl;
  
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  TEfficiency* tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay = new TEfficiency("tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay->SetTitle(";true electron #theta;electron ID efficiency");
  std::vector<TEfficiency*> TEff_vector_dressed_em1500_noOverlay;
  TEff_vector_dressed_em1500_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_noConv_em1500_noOverlay);
  TEff_vector_dressed_em1500_noOverlay.push_back(tEff_dressed_electronVsTrueTheta_AngMatch_1_deg_EMatch_noConv_em1500_noOverlay);
  
  TH1F* h_dressed_em1500_reco_E_over_true_E_match_noOverlay = new TH1F("h_dressed_em1500_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_dressed_em1500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_dressed_em1500_reco_pt_over_true_pt_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  std::vector<TH1F*>hist_vec_dressed_em1500;
  hist_vec_dressed_em1500.push_back(h_dressed_em1500_reco_E_over_true_E_match_noOverlay);
  hist_vec_dressed_em1500.push_back(h_dressed_em1500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_dressed_em1500.size();i++){
    hist_vec_dressed_em1500[i]->Sumw2();
  }

  fillTreeHistograms(file_dressed_electron_em1500,hist_vec_dressed_em1500,test_signal_particle,TEff_vector_all_dressed_em_noOverlay,TEff_vector_dressed_em1500_noOverlay,hist_dressed_vec_all_em);


  file_histogram_dressed_electrons->Write();
  file_histogram_dressed_electrons->Close();
  
  

  TFile* file_muon_mum1=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum1_ILC181101_CT.root");
  TFile* file_muon_mum2=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum2_ILC181101_CT.root");
  TFile* file_muon_mum3=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum3_ILC181101_CT.root");
  TFile* file_muon_mum4=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum4_ILC181101_CT.root");
  TFile* file_muon_mum5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum5_ILC181101_CT.root");
  TFile* file_muon_mum7_5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum7_5_ILC181101_CT.root");
  TFile* file_muon_mum10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum10_ILC181101_CT.root");
  TFile* file_muon_mum20=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum20_ILC181101_CT.root");
  TFile* file_muon_mum50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum50_ILC181101_CT.root");
  TFile* file_muon_mum100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum100_ILC181101_CT.root");
  TFile* file_muon_mum150=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum150_ILC181101_CT.root");
  TFile* file_muon_mum200=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum200_ILC181101_CT.root");
  TFile* file_muon_mum250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum250_ILC181101_CT.root");
  TFile* file_muon_mum400=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum400_ILC181101_CT.root");
  TFile* file_muon_mum500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum500_ILC181101_CT.root");
  TFile* file_muon_mum750=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum750_ILC181101_CT.root");
  TFile* file_muon_mum1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum1000_ILC181101_CT.root");
  TFile* file_muon_mum1250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum1250_ILC181101_CT.root");
  TFile* file_muon_mum1500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_mum1500_ILC181101_CT.root");

  float n_bins_muon=30;

  float lim_rel_energy_low_muon = 0.75;
  float lim_rel_energy_high_muon = 1.25;

  test_signal_particle = 13;
  const unsigned int nbins_teff_mum = 19; 
  double xbins_teff_mum[nbins_teff_mum+1]={0,1.5,2.5,3.5,4.5,6.,9.,12.5,25.,75.,125.,175.,225.,275.,450.,600.,800.,1100.,1300.,1550.};

  const char* final_histo_name_muons="/eos/user/w/weberma2/data/validation181101/pionStudy_mum1_to_mum1500_ILC181101_finalhistos_deg1_match_EMatch_5.root";
  TFile* file_histogram_muons=new TFile(final_histo_name_muons,"recreate");

  TEfficiency* tEff_muonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_mum,xbins_teff_mum);
  tEff_muonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true muon energy [GeV];muon ID efficiency");
  TEfficiency* tEff_muonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_mum,xbins_teff_mum);
  tEff_muonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true muon energy [GeV];muon ID efficiency");

  std::vector<TEfficiency*> TEff_vector_all_mum_noOverlay;
  TEff_vector_all_mum_noOverlay.push_back(tEff_muonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_mum_noOverlay.push_back(tEff_muonVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum1_noOverlay;
  TEff_vector_mum1_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1_noOverlay);
  TEff_vector_mum1_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1_noOverlay);
 
  TH1F* h_mum1_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum1_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum1_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum1_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum1;
  hist_vec_mum1.push_back(h_mum1_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum1.push_back(h_mum1_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum1.size();i++){
    hist_vec_mum1[i]->Sumw2();
  }
  std::vector<TH1F*>hist_vec_all_mu;

  std::cout<<"mum1"<<std::endl;
  fillTreeHistograms(file_muon_mum1,hist_vec_mum1,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum1_noOverlay,hist_vec_all_mu);
  std::cout<<"mum2"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum2_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum2_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum2_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum2_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum2_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum2_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum2_noOverlay;
  TEff_vector_mum2_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum2_noOverlay);
  TEff_vector_mum2_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum2_noOverlay);
    
  TH1F* h_mum2_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum2_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum2_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum2_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum2;
  hist_vec_mum2.push_back(h_mum2_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum2.push_back(h_mum2_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum2.size();i++){
    hist_vec_mum2[i]->Sumw2();
  }
  fillTreeHistograms(file_muon_mum2,hist_vec_mum2,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum2_noOverlay,hist_vec_all_mu);
  std::cout<<"mum3"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum3_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum3_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum3_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum3_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum3_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum3_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum3_noOverlay;
  TEff_vector_mum3_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum3_noOverlay);
  TEff_vector_mum3_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum3_noOverlay);
  
  TH1F* h_mum3_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum3_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum3_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum3_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum3;
  hist_vec_mum3.push_back(h_mum3_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum3.push_back(h_mum3_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum3.size();i++){
    hist_vec_mum3[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum3,hist_vec_mum3,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum3_noOverlay,hist_vec_all_mu);
  std::cout<<"mum4"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum4_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum4_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum4_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum4_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum4_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum4_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum4_noOverlay;
  TEff_vector_mum4_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum4_noOverlay);
  TEff_vector_mum4_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum4_noOverlay);
    
  TH1F* h_mum4_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum4_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum4_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum4_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum4;
  hist_vec_mum4.push_back(h_mum4_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum4.push_back(h_mum4_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum4.size();i++){
    hist_vec_mum4[i]->Sumw2();
  }
  fillTreeHistograms(file_muon_mum4,hist_vec_mum4,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum4_noOverlay,hist_vec_all_mu);
  std::cout<<"mum5"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum5_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum5_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum5_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum5_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum5_noOverlay;
  TEff_vector_mum5_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum5_noOverlay);
  TEff_vector_mum5_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum5_noOverlay); 
  
  TH1F* h_mum5_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum5_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum5_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum5;
  hist_vec_mum5.push_back(h_mum5_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum5.push_back(h_mum5_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum5.size();i++){
    hist_vec_mum5[i]->Sumw2();
  }
  fillTreeHistograms(file_muon_mum5,hist_vec_mum5,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum5_noOverlay,hist_vec_all_mu);
  std::cout<<"mum7.5"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum7_5_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum7_5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum7_5_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum7_5_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum7_5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum7_5_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum7_5_noOverlay;
  TEff_vector_mum7_5_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum7_5_noOverlay);
  TEff_vector_mum7_5_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum7_5_noOverlay);
  
  TH1F* h_mum7_5_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum7_5_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum7_5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum7_5_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum7_5;
  hist_vec_mum7_5.push_back(h_mum7_5_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum7_5.push_back(h_mum7_5_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum7_5.size();i++){
    hist_vec_mum7_5[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum7_5,hist_vec_mum7_5,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum7_5_noOverlay,hist_vec_all_mu);
  std::cout<<"mum10"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum10_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum10_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum10_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum10_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum10_noOverlay;
  TEff_vector_mum10_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum10_noOverlay);
  TEff_vector_mum10_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum10_noOverlay);
    
  TH1F* h_mum10_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum10_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum10_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum10_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum10;
  hist_vec_mum10.push_back(h_mum10_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum10.push_back(h_mum10_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum10.size();i++){
    hist_vec_mum10[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum10,hist_vec_mum10,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum10_noOverlay,hist_vec_all_mu);
  std::cout<<"mum30"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum20_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum20_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum20_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum20_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum20_noOverlay;
  TEff_vector_mum20_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum20_noOverlay);
  TEff_vector_mum20_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum20_noOverlay);
  
  TH1F* h_mum20_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum20_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum20_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum20_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum20;
  hist_vec_mum20.push_back(h_mum20_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum20.push_back(h_mum20_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum20.size();i++){
    hist_vec_mum20[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum20,hist_vec_mum20,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum20_noOverlay,hist_vec_all_mu);
  std::cout<<"mum50"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum50_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum50_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum50_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum50_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum50_noOverlay;
  TEff_vector_mum50_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum50_noOverlay);
  TEff_vector_mum50_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum50_noOverlay);
  
  TH1F* h_mum50_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum50_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum50_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum50_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum50;
  hist_vec_mum50.push_back(h_mum50_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum50.push_back(h_mum50_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum50.size();i++){
    hist_vec_mum50[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum50,hist_vec_mum50,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum50_noOverlay,hist_vec_all_mu);
  std::cout<<"mum100"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum100_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum100_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum100_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum100_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum100_noOverlay;
  TEff_vector_mum100_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum100_noOverlay);
  TEff_vector_mum100_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum100_noOverlay);
    
  TH1F* h_mum100_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum100_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum100_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum100_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum100;
  hist_vec_mum100.push_back(h_mum100_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum100.push_back(h_mum100_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum100.size();i++){
    hist_vec_mum100[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum100,hist_vec_mum100,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum100_noOverlay,hist_vec_all_mu);
  std::cout<<"mum150"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum150_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum150_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum150_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum150_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum150_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum150_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum150_noOverlay;
  TEff_vector_mum150_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum150_noOverlay);
  TEff_vector_mum150_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum150_noOverlay);
    
  TH1F* h_mum150_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum150_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum150_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum150_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum150;
  hist_vec_mum150.push_back(h_mum150_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum150.push_back(h_mum150_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum150.size();i++){
    hist_vec_mum150[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum150,hist_vec_mum150,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum150_noOverlay,hist_vec_all_mu);
  std::cout<<"mum200"<<std::endl;
    
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum200_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum200_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum200_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum200_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum200_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum200_noOverlay;
  TEff_vector_mum200_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum200_noOverlay);
  TEff_vector_mum200_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum200_noOverlay);
  
  TH1F* h_mum200_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum200_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum200_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum200_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum200;
  hist_vec_mum200.push_back(h_mum200_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum200.push_back(h_mum200_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum200.size();i++){
    hist_vec_mum200[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum200,hist_vec_mum200,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum200_noOverlay,hist_vec_all_mu);
  std::cout<<"mum250"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum250_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum250_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum250_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum250_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum250_noOverlay;
  TEff_vector_mum250_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum250_noOverlay);
  TEff_vector_mum250_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum250_noOverlay);
    
  TH1F* h_mum250_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum250_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum250_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum250_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum250;
  hist_vec_mum250.push_back(h_mum250_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum250.push_back(h_mum250_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum250.size();i++){
    hist_vec_mum250[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum250,hist_vec_mum250,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum250_noOverlay,hist_vec_all_mu);
  std::cout<<"mum400"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum400_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum400_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum400_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum400_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum400_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum400_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum400_noOverlay;
  TEff_vector_mum400_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum400_noOverlay);
  TEff_vector_mum400_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum400_noOverlay);
    
  TH1F* h_mum400_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum400_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum400_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum400_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum400;
  hist_vec_mum400.push_back(h_mum400_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum400.push_back(h_mum400_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum400.size();i++){
    hist_vec_mum400[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum400,hist_vec_mum400,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum400_noOverlay,hist_vec_all_mu);
  std::cout<<"mum500"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum500_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum500_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum500_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum500_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum500_noOverlay;
  TEff_vector_mum500_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum500_noOverlay);
  TEff_vector_mum500_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum500_noOverlay);
  
  TH1F* h_mum500_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum500_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum500_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum500;
  hist_vec_mum500.push_back(h_mum500_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum500.push_back(h_mum500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum500.size();i++){
    hist_vec_mum500[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum500,hist_vec_mum500,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum500_noOverlay,hist_vec_all_mu);
  std::cout<<"mum750"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum750_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum750_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum750_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum750_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum750_noOverlay;
  TEff_vector_mum750_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum750_noOverlay);
  TEff_vector_mum750_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum750_noOverlay);
    
  TH1F* h_mum750_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum750_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum750_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum750_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum750;
  hist_vec_mum750.push_back(h_mum750_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum750.push_back(h_mum750_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum750.size();i++){
    hist_vec_mum750[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum750,hist_vec_mum750,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum750_noOverlay,hist_vec_all_mu);
  std::cout<<"mum1000"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1000_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1000_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1000_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1000_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum1000_noOverlay;
  TEff_vector_mum1000_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1000_noOverlay);
  TEff_vector_mum1000_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1000_noOverlay);

  TH1F* h_mum1000_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum1000_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum1000_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum1000_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum1000;
  hist_vec_mum1000.push_back(h_mum1000_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum1000.push_back(h_mum1000_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum1000.size();i++){
    hist_vec_mum1000[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum1000,hist_vec_mum1000,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum1000_noOverlay,hist_vec_all_mu);
  std::cout<<"mum1250"<<std::endl;
    
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1250_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1250_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1250_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1250_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum1250_noOverlay;
  TEff_vector_mum1250_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1250_noOverlay);
  TEff_vector_mum1250_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1250_noOverlay);
    
  TH1F* h_mum1250_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum1250_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum1250_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum1250_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum1250;
  hist_vec_mum1250.push_back(h_mum1250_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum1250.push_back(h_mum1250_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum1250.size();i++){
    hist_vec_mum1250[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum1250,hist_vec_mum1250,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum1250_noOverlay,hist_vec_all_mu);
  std::cout<<"mum1500"<<std::endl; 
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1500_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1500_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  TEfficiency* tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1500_noOverlay = new TEfficiency("tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1500_noOverlay->SetTitle(";true muon #theta;muon ID efficiency");
  std::vector<TEfficiency*> TEff_vector_mum1500_noOverlay;
  TEff_vector_mum1500_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_noConv_mum1500_noOverlay);
  TEff_vector_mum1500_noOverlay.push_back(tEff_muonVsTrueTheta_AngMatch_1_deg_EMatch_noConv_mum1500_noOverlay);
    
  TH1F* h_mum1500_reco_E_over_true_E_match_noOverlay = new TH1F("h_mum1500_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_mum1500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_mum1500_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_mum1500;
  hist_vec_mum1500.push_back(h_mum1500_reco_E_over_true_E_match_noOverlay);
  hist_vec_mum1500.push_back(h_mum1500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_mum1500.size();i++){
    hist_vec_mum1500[i]->Sumw2();
  }

  fillTreeHistograms(file_muon_mum1500,hist_vec_mum1500,test_signal_particle,TEff_vector_all_mum_noOverlay,TEff_vector_mum1500_noOverlay,hist_vec_all_mu);


  file_histogram_muons->Write();
  file_histogram_muons->Close();



  //1,10,20,50,75,100,150,250,500,1000,1500

  TFile* file_pion_pim1=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim1_ILC181101_CT.root");
  TFile* file_pion_pim2=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim2_ILC181101_CT.root");
  TFile* file_pion_pim5=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim5_ILC181101_CT.root");
  TFile* file_pion_pim10=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim10_ILC181101_CT.root");
  TFile* file_pion_pim20=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim20_ILC181101_CT.root");
  TFile* file_pion_pim50=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim50_ILC181101_CT.root");
  TFile* file_pion_pim75=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim75_ILC181101_CT.root");
  TFile* file_pion_pim100=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim100_ILC181101_CT.root");
  TFile* file_pion_pim150=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim150_ILC181101_CT.root");
  TFile* file_pion_pim250=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim250_ILC181101_CT.root");
  TFile* file_pion_pim500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim500_ILC181101_CT.root");
  TFile* file_pion_pim750=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim750_ILC181101_CT.root");
  TFile* file_pion_pim1000=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim1000_ILC181101_CT.root");
  TFile* file_pion_pim1500=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/181101_gcc62/pionStudy_pim1500_ILC181101_CT.root");

  test_signal_particle = 211;
  const unsigned int nbins_teff_pim = 14; 
  double xbins_teff_pim[nbins_teff_pim+1]={0,1.5,4.0,4.5,7.5,15.,30.,85.,125.,200.,275.,600.,850.,1100.,1550.};

  const char* final_histo_name_pions="/eos/user/w/weberma2/data/validation181101/pionStudy_pim1_to_pim1500_ILC181101_finalhistos_deg1_match_EMatch_5.root";
  TFile* file_histogram_pions=new TFile(final_histo_name_pions,"recreate");

  TEfficiency* tEff_pionVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PionVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_pim,xbins_teff_pim);
  tEff_pionVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true pion energy [GeV];pion ID efficiency");
  TEfficiency* tEff_pionVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PionVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_pim,xbins_teff_pim);
  tEff_pionVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true pion energy [GeV];pion ID efficiency");

  std::vector<TEfficiency*> TEff_vector_all_pim_noOverlay;
  TEff_vector_all_pim_noOverlay.push_back(tEff_pionVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_pim_noOverlay.push_back(tEff_pionVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);

  std::cout<<"pim1"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim1_noOverlay;
  TEff_vector_pim1_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay);
  TEff_vector_pim1_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay);
  
  TH1F* h_pim1_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim1_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim1_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim1_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim1;
  hist_vec_pim1.push_back(h_pim1_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim1.push_back(h_pim1_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim1.size();i++){
    hist_vec_pim1[i]->Sumw2();
  }
  std::vector<TH1F*>hist_vec_all_pi;

  fillTreeHistograms(file_pion_pim1,hist_vec_pim1,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim1_noOverlay,hist_vec_all_pi);
  std::cout<<"pim2"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim2_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim2_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim2_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim2_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim2_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim2_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim2_noOverlay;
  TEff_vector_pim2_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim2_noOverlay);
  TEff_vector_pim2_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim2_noOverlay);
  
  TH1F* h_pim2_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim2_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim2_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim2_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim2;
  hist_vec_pim2.push_back(h_pim2_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim2.push_back(h_pim2_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim2.size();i++){
    hist_vec_pim2[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim2,hist_vec_pim2,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim2_noOverlay,hist_vec_all_pi);
  std::cout<<"pim5"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim5_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim5_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim5_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim5_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim5_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim5_noOverlay;
  TEff_vector_pim5_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim5_noOverlay);
  TEff_vector_pim5_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim5_noOverlay);
  
  TH1F* h_pim5_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim5_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim5_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim5_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim5;
  hist_vec_pim5.push_back(h_pim5_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim5.push_back(h_pim5_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim5.size();i++){
    hist_vec_pim5[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim5,hist_vec_pim5,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim5_noOverlay,hist_vec_all_pi);
  std::cout<<"pim10"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim10_noOverlay;
  TEff_vector_pim10_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay);
  TEff_vector_pim10_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay);
  
  TH1F* h_pim10_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim10_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim10_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim10_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim10;
  hist_vec_pim10.push_back(h_pim10_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim10.push_back(h_pim10_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim10.size();i++){
    hist_vec_pim10[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim10,hist_vec_pim10,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim10_noOverlay,hist_vec_all_pi);
  std::cout<<"pim20"<<std::endl;
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim20_noOverlay;
  TEff_vector_pim20_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay);
  TEff_vector_pim20_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay);
  
  TH1F* h_pim20_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim20_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim20_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim20_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim20;
  hist_vec_pim20.push_back(h_pim20_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim20.push_back(h_pim20_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim20.size();i++){
    hist_vec_pim20[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim20,hist_vec_pim20,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim20_noOverlay,hist_vec_all_pi);
  std::cout<<"pim50"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim50_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim50_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim50_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim50_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim50_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim50_noOverlay;
  TEff_vector_pim50_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim50_noOverlay);
  TEff_vector_pim50_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim50_noOverlay);
  
  TH1F* h_pim50_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim50_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim50_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim50_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim50;
  hist_vec_pim50.push_back(h_pim50_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim50.push_back(h_pim50_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim50.size();i++){
    hist_vec_pim50[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim50,hist_vec_pim50,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim50_noOverlay,hist_vec_all_pi);
  std::cout<<"pim75"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim75_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim75_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim75_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim75_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim75_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim75_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim75_noOverlay;
  TEff_vector_pim75_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim75_noOverlay);
  TEff_vector_pim75_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim75_noOverlay);
  
  TH1F* h_pim75_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim75_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim75_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim75_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim75;
  hist_vec_pim75.push_back(h_pim75_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim75.push_back(h_pim75_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim75.size();i++){
    hist_vec_pim75[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim75,hist_vec_pim75,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim75_noOverlay,hist_vec_all_pi);
  std::cout<<"pim100"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim100_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim100_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim100_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim100_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim100_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim100_noOverlay;
  TEff_vector_pim100_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim100_noOverlay);
  TEff_vector_pim100_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim100_noOverlay);
  
  TH1F* h_pim100_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim100_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim100_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim100_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim100;
  hist_vec_pim100.push_back(h_pim100_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim100.push_back(h_pim100_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim100.size();i++){
    hist_vec_pim100[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim100,hist_vec_pim100,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim100_noOverlay,hist_vec_all_pi);
  std::cout<<"pim150"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim150_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim150_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim150_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim150_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim150_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim150_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim150_noOverlay;
  TEff_vector_pim150_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim150_noOverlay);
  TEff_vector_pim150_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim150_noOverlay);
  
  fillTreeHistograms(file_pion_pim150,hist_vec_pim1,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim150_noOverlay,hist_vec_all_pi);
  std::cout<<"pim250"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim250_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim250_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim250_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim250_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim250_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim250_noOverlay;
  TEff_vector_pim250_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim250_noOverlay);
  TEff_vector_pim250_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim250_noOverlay);
  
  TH1F* h_pim250_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim250_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim250_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim250_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim250;
  hist_vec_pim250.push_back(h_pim250_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim250.push_back(h_pim250_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim250.size();i++){
    hist_vec_pim250[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim250,hist_vec_pim250,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim250_noOverlay,hist_vec_all_pi);
  std::cout<<"pim500"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim500_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim500_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim500_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim500_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim500_noOverlay;
  TEff_vector_pim500_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim500_noOverlay);
  TEff_vector_pim500_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim500_noOverlay);
  
  TH1F* h_pim500_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim500_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim500_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim500;
  hist_vec_pim500.push_back(h_pim500_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim500.push_back(h_pim500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim500.size();i++){
    hist_vec_pim500[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim500,hist_vec_pim500,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim500_noOverlay,hist_vec_all_pi);
  std::cout<<"pim750"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim750_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim750_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim750_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim750_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim750_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim750_noOverlay;
  TEff_vector_pim750_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim750_noOverlay);
  TEff_vector_pim750_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim750_noOverlay);
  
  TH1F* h_pim750_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim750_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim750_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim750_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim750;
  hist_vec_pim750.push_back(h_pim750_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim750.push_back(h_pim750_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim750.size();i++){
    hist_vec_pim750[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim750,hist_vec_pim750,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim750_noOverlay,hist_vec_all_pi);
  std::cout<<"pim1000"<<std::endl;
  
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1000_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1000_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1000_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1000_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1000_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim1000_noOverlay;
  TEff_vector_pim1000_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1000_noOverlay);
  TEff_vector_pim1000_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1000_noOverlay);
  
  TH1F* h_pim1000_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim1000_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim1000_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim1000_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim1000;
  hist_vec_pim1000.push_back(h_pim1000_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim1000.push_back(h_pim1000_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim1000.size();i++){
    hist_vec_pim1000[i]->Sumw2();
  }

  fillTreeHistograms(file_pion_pim1000,hist_vec_pim1000,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim1000_noOverlay,hist_vec_all_pi);
  std::cout<<"pim1500"<<std::endl;
    
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1500_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1500_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1500_noOverlay = new TEfficiency("tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1500_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1500_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim1500_noOverlay;
  TEff_vector_pim1500_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_noConv_pim1500_noOverlay);
  TEff_vector_pim1500_noOverlay.push_back(tEff_pionVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1500_noOverlay);

  TH1F* h_pim1500_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim1500_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim1500_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim1500_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim1500;
  hist_vec_pim1500.push_back(h_pim1500_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim1500.push_back(h_pim1500_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim1500.size();i++){
    hist_vec_pim1500[i]->Sumw2();
  }
  
  fillTreeHistograms(file_pion_pim1500,hist_vec_pim1500,test_signal_particle,TEff_vector_all_pim_noOverlay,TEff_vector_pim1500_noOverlay,hist_vec_all_pi);
  

  file_histogram_pions->Write();
  file_histogram_pions->Close();






  /*
  TFile* file_pion_pim1_FW=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180517/pionStudy_pim1_FW_180517_CT_CLIC_o3_v14_DR200_redInfo.root");
  TFile* file_pion_pim10_FW=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180517/pionStudy_pim10_FW_180517_CT_CLIC_o3_v14_DR200_redInfo.root");
  TFile* file_pion_pim20_FW=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180517/pionStudy_pim20_FW_180517_CT_CLIC_o3_v14_DR200_redInfo.root");
  TFile* file_pion_pim10_180517=TFile::Open("/eos/user/w/weberma2/data/pionMacroFiles/180517/pionStudy_pim10_180517_CT_CLIC_o3_v14_DR200_redInfo.root");

  const char* final_histo_name_pions_FW="/eos/user/w/weberma2/data/validation171212/pionStudy_pim1_to_pim10_forward_ILC180517_CT_CLIC_o3_v14_SWC_CLIC_default_finalhistos_deg1_match_EMatch_5.root";
  TFile* file_histogram_pions_FW=new TFile(final_histo_name_pions,"recreate");

  TEfficiency* tEff_pionFWVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PionFWVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_pim,xbins_teff_pim);
  tEff_pionFWVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true pion energy [GeV];pion ID efficiency");
  TEfficiency* tEff_pionFWVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay = new TEfficiency("tEff_PionFWVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay","", nbins_teff_pim,xbins_teff_pim);
  tEff_pionFWVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay->SetTitle(";true pion energy [GeV];pion ID efficiency");

  std::vector<TEfficiency*> TEff_vector_all_pim_FW_noOverlay;
  TEff_vector_all_pim_FW_noOverlay.push_back(tEff_pionFWVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noConv_noOverlay);
  TEff_vector_all_pim_FW_noOverlay.push_back(tEff_pionFWVsTrueE_AngMatch_1_deg_EMatch_CosTheta_True_0_95_noConv_noOverlay);

  std::cout<<"pim1 FW"<<std::endl;
  
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim1_FW_noOverlay;
  TEff_vector_pim1_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim1_noOverlay);
  TEff_vector_pim1_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim1_noOverlay);
  
  TH1F* h_pim1_FW_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim1_FW_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim1_FW_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim1_FW_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim1_FW;
  hist_vec_pim1_FW.push_back(h_pim1_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim1_FW.push_back(h_pim1_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim1_FW.size();i++){
    hist_vec_pim1_FW[i]->Sumw2();
  }
  std::vector<TH1F*>hist_vec_all_pi_FW;

  fillTreeHistograms(file_pion_pim1_FW,hist_vec_pim1_FW,test_signal_particle,TEff_vector_all_pim_FW_noOverlay,TEff_vector_pim1_FW_noOverlay,hist_vec_all_pi_FW);

  std::cout<<"pim10 FW"<<std::endl;
  
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim10_FW_noOverlay;
  TEff_vector_pim10_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim10_noOverlay);
  TEff_vector_pim10_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim10_noOverlay);
  
  TH1F* h_pim10_FW_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim10_FW_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim10_FW_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim10_FW_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim10_FW;
  hist_vec_pim10_FW.push_back(h_pim10_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim10_FW.push_back(h_pim10_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim10_FW.size();i++){
    hist_vec_pim10_FW[i]->Sumw2();
  }

 fillTreeHistograms(file_pion_pim20_FW,hist_vec_pim10_FW,test_signal_particle,TEff_vector_all_pim_FW_noOverlay,TEff_vector_pim10_FW_noOverlay,hist_vec_all_pi_FW);


  std::cout<<"pim20 FW"<<std::endl;
  
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  TEfficiency* tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay = new TEfficiency("tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay","", n_bins,lim_theta_low,lim_theta_high);
  tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay->SetTitle(";true pion #theta;pion ID efficiency");
  std::vector<TEfficiency*> TEff_vector_pim20_FW_noOverlay;
  TEff_vector_pim20_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_noConv_pim20_noOverlay);
  TEff_vector_pim20_FW_noOverlay.push_back(tEff_pionFWVsTrueTheta_AngMatch_1_deg_EMatch_noConv_pim20_noOverlay);
  
  TH1F* h_pim20_FW_reco_E_over_true_E_match_noOverlay = new TH1F("h_pim20_FW_reco_E_over_true_E_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  TH1F* h_pim20_FW_reco_pt_over_true_pt_match_noOverlay = new TH1F("h_pim20_FW_reco_pt_over_true_pt_match_noOverlay","", n_bins_muon,lim_rel_energy_low_muon,lim_rel_energy_high_muon);
  std::vector<TH1F*>hist_vec_pim20_FW;
  hist_vec_pim20_FW.push_back(h_pim20_reco_E_over_true_E_match_noOverlay);
  hist_vec_pim20_FW.push_back(h_pim20_reco_pt_over_true_pt_match_noOverlay);
  for(unsigned int i=0;i<hist_vec_pim20_FW.size();i++){
    hist_vec_pim20_FW[i]->Sumw2();
  }

 fillTreeHistograms(file_pion_pim20_FW,hist_vec_pim20_FW,test_signal_particle,TEff_vector_all_pim_FW_noOverlay,TEff_vector_pim20_FW_noOverlay,hist_vec_all_pi_FW);

  file_histogram_pions_FW->Write();
  file_histogram_pions_FW->Close();
  */
    //500 new tune    //  mean 500.076, sigma 4.84/rel 0.968 %
    // mean 9.999, sigma 0.0310, real 0.3125
    //mean 0.9986, sigma 0.003887, 0.389
    //mean 19.999, sigma 0.0600, 0.300
    // mean 49.99, sigma 0.165, 0.329
    //mean 75.001, sigma 0.263, 0.3510
    //mean 99.998, sigma 0.3714, 0.3715
    
  return 1;

}

