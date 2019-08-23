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
#include <TVirtualFitter.h>
#include "TProfile.h"
#include "TColor.h"
#include <vector>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
//#include "RooNovosibirsk.h"

using namespace std;

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

double fnc_dscb(double*xx,double*pp)
{
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];
  
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

std::vector<float> sigmaGaussianFits(TH1F* hist){
  float sigma = 1.0;
  float sigmaError = 1.0;

  TF1* fit_raw_range_0_8_to_1_2 = new TF1("fit_sigma_range_0_8_to_1_2", "gaus",0.0,3.0);
  hist->Fit("fit_sigma_range_0_8_to_1_2","R");
  //sigma is the parameter 2, mean is parameter 1
  float sigma_temp=fit_raw_range_0_8_to_1_2->GetParameter(2);
  TF1* fit_raw_range_two_sigma = new TF1("fit_sigma_range_two_sigma", "gaus",1.0-3.*sigma_temp,1.0+3.*sigma_temp);
  hist->Fit("fit_sigma_range_two_sigma","R");
  sigma=fit_raw_range_two_sigma->GetParameter(2);
  sigmaError=fit_raw_range_two_sigma->GetParError(2);

  std::vector<float> SigmaWError;
  SigmaWError.push_back(sigma);
  SigmaWError.push_back(sigmaError);
  float sigma_before=fit_raw_range_0_8_to_1_2->GetParameter(2);
  float sigma_beforeError=fit_raw_range_0_8_to_1_2->GetParError(2);

  for(unsigned int i=0;i<1000;i++){
    TF1* fit_raw_range_temp = new TF1("fit_sigma_range_temp", "gaus",1.0-3.*sigma_before,1.0+3.*sigma_before);
    hist->Fit("fit_sigma_range_temp","R");
    float sigma_after=fit_raw_range_temp->GetParameter(2);
    float sigma_afterError=fit_raw_range_temp->GetParError(2);
    if((sigma_after/sigma_before)<1.02 && (sigma_after/sigma_before)>0.98){
      sigma_before=sigma_after;
      sigma_beforeError=sigma_afterError;
      break;
    }
    sigma_before=sigma_after;
    sigma_beforeError=sigma_afterError;
 }
  SigmaWError.push_back(sigma_before);
  SigmaWError.push_back(sigma_beforeError);

  return SigmaWError;

}


std::vector<float> sigmaGaussianFits2Sigma(TH1F* hist){
  float sigma = 1.0;
  float sigmaError = 1.0;

  TF1* fit_raw_range_0_8_to_1_2 = new TF1("fit_sigma_range_0_8_to_1_2", "gaus",0.5,1.5);
  hist->Fit("fit_sigma_range_0_8_to_1_2","R");
  //sigma is the parameter 2, mean is parameter 1
  float sigma_temp=fit_raw_range_0_8_to_1_2->GetParameter(2);
  TF1* fit_raw_range_two_sigma = new TF1("fit_sigma_range_two_sigma", "gaus",1.0-2.*sigma_temp,1.0+2.*sigma_temp);
  hist->Fit("fit_sigma_range_two_sigma","R");
  sigma=fit_raw_range_two_sigma->GetParameter(2);
  sigmaError=fit_raw_range_two_sigma->GetParError(2);

  std::vector<float> SigmaWError;
  SigmaWError.push_back(sigma);
  SigmaWError.push_back(sigmaError);
  float sigma_before=fit_raw_range_0_8_to_1_2->GetParameter(2);
  float sigma_beforeError=fit_raw_range_0_8_to_1_2->GetParError(2);

  for(unsigned int i=0;i<1000;i++){
    TF1* fit_raw_range_temp = new TF1("fit_sigma_range_temp", "gaus",1.0-2.*sigma_before,1.0+2.*sigma_before);
    hist->Fit("fit_sigma_range_temp","R");
    float sigma_after=fit_raw_range_temp->GetParameter(2);
    float sigma_afterError=fit_raw_range_temp->GetParError(2);
    if((sigma_after/sigma_before)<1.02 && (sigma_after/sigma_before)>0.98){
      sigma_before=sigma_after;
      sigma_beforeError=sigma_afterError;
      break;
    }
    sigma_before=sigma_after;
    sigma_beforeError=sigma_afterError;
 }
  SigmaWError.push_back(sigma_before);
  SigmaWError.push_back(sigma_beforeError);

  return SigmaWError;

}

float DeltaPhi(float Phi1,float Phi2){
  float deltaphi=fabs(Phi1-Phi2);
  if(deltaphi>M_PI){
    deltaphi=2*M_PI-deltaphi;
  }
  return deltaphi;
}

float calculate_over90(TH1F* hist){

  if(hist->GetEntries()==0){
    return 0;
  }

  float sum_events=0;
  for (int i=0;i<=(hist->GetNbinsX()+1);i++){
    if(hist->GetBinLowEdge(i)>0.90){
      sum_events+=hist->GetBinContent(i);
    }
  }
  return sum_events/(hist->GetBinContent(0)+hist->Integral()+hist->GetBinContent(hist->GetNbinsX()+1));
}

void CalculatePerformance(const TH1F *const pTH1F, float &resolution, float &resolutionError, float &mean, float &meanError)
{
  //expects an average of two quantities --> originally from total event energy, thus multiply by a factor of sqrt(2) later
  static const float FLOAT_MAX(std::numeric_limits<float>::max());
  
  if (NULL == pTH1F)
    return;
  
  if (5 > pTH1F->GetEntries())
    {
      std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
      return;
    }
  
  // Calculate raw properties of distribution
  float sum = 0., total = 0.;
  double sx = 0., sxx = 0.;
  const unsigned int nbins(pTH1F->GetNbinsX());
  
  for (unsigned int i = 0; i <= nbins; ++i)
    {
      const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
      const float yi(pTH1F->GetBinContent(i));
      sx += yi * binx;
      sxx += yi * binx * binx;
      total += yi;
    }
  
  const float rawMean(sx / total);
  const float rawMeanSqua(sxx / total);
  const float rawRms(std::sqrt(rawMeanSqua - rawMean * rawMean));
  
  sum = 0.;
  unsigned int is0 = 0;
  
  for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
      sum += pTH1F->GetBinContent(i);
      is0 = i;
    }
  
  // Calculate truncated properties
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX)/*, mean(FLOAT_MAX)*/, low(FLOAT_MAX), rms(FLOAT_MAX);
  float high(0.f);
  
  for (unsigned int istart = 0; istart <= is0; ++istart)
    {
      double sumn = 0.;
      double csum = 0.;
      double sumx = 0.;
      double sumxx = 0.;
      unsigned int iend = 0;
      
      for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
	  const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
	  const float yi(pTH1F->GetBinContent(i));
	  csum += yi;
	  
	  if (sumn < 0.9 * total)
            {
	      sumn += yi;
	      sumx += yi * binx;
	      sumxx+= yi * binx * binx;
	      iend = i;
            }
        }
      
        const float localMean(sumx / sumn);
        const float localMeanSqua(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSqua - localMean * localMean));
	
        if (localRms < rmsmin)
	  {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;
	    frac = rms / mean * std::sqrt(2) * 100.;
	    efrac = frac / std::sqrt(total);
	    meanError = mean / std::sqrt(total); 
	  }
	resolution = frac;
	resolutionError = efrac;
    }
  std::cout<<"resolution/error "<<resolution<<"/"<<resolutionError<<" mean "<<mean<<std::endl;
}


void CalculatePerformanceNoMean(const TH1F *const pTH1F, float &resolution, float &resolutionError, float &mean, float &meanError)
{
  static const float FLOAT_MAX(std::numeric_limits<float>::max());
  
  if (NULL == pTH1F)
    return;
  
  if (5 > pTH1F->GetEntries())
    {
      std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
      return;
    }
  
  // Calculate raw properties of distribution
  float sum = 0., total = 0.;
  double sx = 0., sxx = 0.;
  const unsigned int nbins(pTH1F->GetNbinsX());
  
  for (unsigned int i = 0; i <= nbins; ++i)
    {
      const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
      const float yi(pTH1F->GetBinContent(i));
      sx += yi * binx;
      sxx += yi * binx * binx;
      total += yi;
    }
  
  const float rawMean(sx / total);
  const float rawMeanSqua(sxx / total);
  const float rawRms(std::sqrt(rawMeanSqua - rawMean * rawMean));

  sum = 0.;
  unsigned int is0 = 0;
  
  for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
      sum += pTH1F->GetBinContent(i);
      is0 = i;
    }
  
  // Calculate truncated properties
  float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX)/*, mean(FLOAT_MAX)*/, low(FLOAT_MAX), rms(FLOAT_MAX);
  float high(0.f);
  
  for (unsigned int istart = 0; istart <= is0; ++istart)
    {
      double sumn = 0.;
      double csum = 0.;
      double sumx = 0.;
      double sumxx = 0.;
      unsigned int iend = 0;
      
      for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
	  const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
	  const float yi(pTH1F->GetBinContent(i));
	  csum += yi;
	  
	  if (sumn < 0.9 * total)
            {
	      sumn += yi;
	      sumx += yi * binx;
	      sumxx+= yi * binx * binx;
	      iend = i;
            }
        }
      
        const float localMean(sumx / sumn);
        const float localMeanSqua(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSqua - localMean * localMean));
	
        if (localRms < rmsmin)
	  {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;
	    frac = rms;
	    efrac = frac / std::sqrt(total);
	    meanError = mean / std::sqrt(total);
	  }
	resolution = frac;
	resolutionError = efrac;
    }

}



TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


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
   gStyle->SetTitleSize(0.06,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   gStyle->SetLabelOffset(0.015,"xyz");
   gStyle->SetTitleOffset(1.2,"yz"); //equivalent to: gStyle->SetTitleYOffset(1.2);
   gStyle->SetTitleOffset(1.17,"x");



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
   gStyle->SetFuncColor(kBlack);
   gStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadTopMargin(0.11);
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


TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,0,50,800,300);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fill_mass_histograms(TFile* file, std::vector<TH1F*> h_hist_vec){

  TTree* tree = (TTree*)file->Get("showerData");

  float Px_trueAll=0;
  float Py_trueAll=0;
  float Pz_trueAll=0;
  float E_trueAll=0;

  float Px_totPFO=0;
  float Py_totPFO=0;
  float Pz_totPFO=0;
  float E_totPFO=0;

  vector<float> *trueME_Px=0;
  vector<float> *trueME_Py=0;
  vector<float> *trueME_Pz=0;
  vector<float> *trueME_E=0;
  vector<int> *trueME_PDGID=0;

  vector<float> *genJetR07Px=0;
  vector<float> *genJetR07Py=0;
  vector<float> *genJetR07Pz=0;
  vector<float> *genJetR07E=0;


  vector<float> *recoJetR07Px=0;
  vector<float> *recoJetR07Py=0;
  vector<float> *recoJetR07Pz=0;
  vector<float> *recoJetR07E=0;


  tree->SetBranchAddress("E_trueAll", &E_trueAll);
  tree->SetBranchAddress("Px_trueAll", &Px_trueAll);
  tree->SetBranchAddress("Py_trueAll", &Py_trueAll);
  tree->SetBranchAddress("Pz_trueAll", &Pz_trueAll);

  tree->SetBranchAddress("E_totPFO", &E_totPFO);
  tree->SetBranchAddress("Px_totPFO", &Px_totPFO);
  tree->SetBranchAddress("Py_totPFO", &Py_totPFO);
  tree->SetBranchAddress("Pz_totPFO", &Pz_totPFO);

  /*
  tree->SetBranchAddress("E_totPFO_cos09", &E_totPFO);
  tree->SetBranchAddress("Px_totPFO_cos09", &Px_totPFO);
  tree->SetBranchAddress("Py_totPFO_cos09", &Py_totPFO);
  tree->SetBranchAddress("Pz_totPFO_cos09", &Pz_totPFO);
  */
  
  tree->SetBranchAddress("trueME_E", &trueME_E);
  tree->SetBranchAddress("trueME_Px", &trueME_Px);
  tree->SetBranchAddress("trueME_Py", &trueME_Py);
  tree->SetBranchAddress("trueME_Pz", &trueME_Pz);
  tree->SetBranchAddress("trueME_PDGID", &trueME_PDGID);

  tree->SetBranchAddress("genJetR07E", &genJetR07E);
  tree->SetBranchAddress("genJetR07Px", &genJetR07Px);
  tree->SetBranchAddress("genJetR07Py", &genJetR07Py);
  tree->SetBranchAddress("genJetR07Pz", &genJetR07Pz);

  tree->SetBranchAddress("recoJetR07E", &recoJetR07E);
  tree->SetBranchAddress("recoJetR07Px", &recoJetR07Px);
  tree->SetBranchAddress("recoJetR07Py", &recoJetR07Py);
  tree->SetBranchAddress("recoJetR07Pz", &recoJetR07Pz);

  double m_costheta_max = 0.9;

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);
 
   //atan2(y,x) --> so for MET atan2(-y_vis,-x_vis)
    float phi_MET_gen=atan2(-Py_trueAll,-Px_trueAll);
    float phi_MET_reco=atan2(-Py_totPFO,-Px_totPFO);
    
    //reevaluate boson ME's using the diquarks
    TLorentzVector tempME_boson(0,0,0,0);
    unsigned int ind_quarks=0;
    for(int m=0;m<trueME_E->size();m++){
      if(abs((*trueME_PDGID)[m])<7){
	ind_quarks+=1;
	//consider the sum of two quarks
	TLorentzVector tempME_quarks(0,0,0,0);
	tempME_quarks.SetPxPyPzE((*trueME_Px)[0],(*trueME_Py)[0],(*trueME_Pz)[0],(*trueME_E)[0]);
	tempME_boson+=tempME_quarks;
      }
    }
    
    //if(ind_quarks!=2){
    //std::cout<<"should have found exactly two quarks, but found instead "<<ind_quarks<<std::endl;
    //}

    TLorentzVector tempGenJet1(0,0,0,0);
    TLorentzVector tempGenJet2(0,0,0,0);
    if(genJetR07E->size()>1){
      tempGenJet1.SetPxPyPzE((*genJetR07Px)[0],(*genJetR07Py)[0],(*genJetR07Pz)[0],(*genJetR07E)[0]);
      tempGenJet2.SetPxPyPzE((*genJetR07Px)[1],(*genJetR07Py)[1],(*genJetR07Pz)[1],(*genJetR07E)[1]);
    }

    
    // if((tempGenJet2.Angle(tempGenJet1.Vect())*TMath::RadToDeg())<-0.01){
    //continue;
    //}

    TLorentzVector tempMEVecQuarks(0,0,0,0);
    bool boson_is_Z=false;
    bool boson_is_W=false;
    float theta_min=180;
    float theta_max=90;
    float theta1=0;
    float theta2=0;
    float phi1=0;
    float phi2=0;
    float lep1E=0;
    unsigned int count_quark=0;
    for(unsigned int m=0;m<trueME_E->size();m++){
      if((*trueME_PDGID)[m]==23){
	boson_is_Z=true;
      }
      if(abs((*trueME_PDGID)[m])==24){
	boson_is_W=true;
      }
      if(abs((*trueME_PDGID)[m])==11 || abs((*trueME_PDGID)[m])==13 ){
	lep1E=(*trueME_E)[m];
      }

      if(abs((*trueME_PDGID)[m])<7){
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE((*trueME_Px)[m],(*trueME_Py)[m],(*trueME_Pz)[m],(*trueME_E)[m]);
	tempMEVecQuarks+=temp;
	if(count_quark==0){
	  theta1=temp.Theta();
	  phi1=temp.Phi();
	  count_quark+=1;
	}else{
	  theta2=temp.Theta();
	  phi2=temp.Phi();
	}
	if(fabs(temp.CosTheta())<fabs(cos(theta_min*TMath::DegToRad()))){
	  theta_min=temp.Theta()*TMath::RadToDeg();
	}
	if(fabs(temp.CosTheta())>fabs(cos(theta_max*TMath::DegToRad()))){
	  theta_max=temp.Theta()*TMath::RadToDeg();
	}
	//std::cout<<i_entry<<" q "<<m<<"/"<<temp.CosTheta()<<"/"<<theta_min<<"/"<<theta_max<<"/"<<cos(theta_min*TMath::DegToRad())<<"/"<<cos(theta_max*TMath::DegToRad())<<"/"<<std::endl;
      }
      if(abs((*trueME_PDGID)[m])==24 || (*trueME_PDGID)[m]==23){
	h_hist_vec[4]->Fill(sqrt((*trueME_E)[m]*(*trueME_E)[m]-(*trueME_Px)[m]*(*trueME_Px)[m]-(*trueME_Py)[m]*(*trueME_Py)[m]-(*trueME_Pz)[m]*(*trueME_Pz)[m]));
      }
    }

    if(boson_is_Z && ((tempGenJet1.E()+tempGenJet2.E())/E_trueAll)<0.9){
      //continue;
    }
    if(boson_is_W && ((tempGenJet1.E()+tempGenJet2.E())/(E_trueAll-lep1E))<0.9){
      //continue;
    }

    if(genJetR07E->size()>1){
      if(fabs(tempGenJet1.CosTheta())<m_costheta_max && fabs(tempGenJet2.CosTheta())<m_costheta_max){
	h_hist_vec[0]->Fill(tempGenJet1.E()+tempGenJet2.E());
	h_hist_vec[2]->Fill((tempGenJet1+tempGenJet2).M());
	if((tempGenJet1+tempGenJet2).M()<60){
	  h_hist_vec[6]->Fill(DeltaPhi(phi_MET_gen,(tempGenJet1+tempGenJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[6]->Fill(fabs(theta_max));
	  //h_hist_vec[6]->Fill(fabs(theta1-theta2)*TMath::RadToDeg());
	  //if(boson_is_Z){
	    //h_hist_vec[6]->Fill((tempGenJet1.E()+tempGenJet2.E())/E_trueAll);
	    //}else{
	    //h_hist_vec[6]->Fill((tempGenJet1.E()+tempGenJet2.E())/(E_trueAll-lep1E));
	    //}
	}
	if(boson_is_Z && (tempGenJet1+tempGenJet2).M()>85.){
	  h_hist_vec[7]->Fill(DeltaPhi(phi_MET_gen,(tempGenJet1+tempGenJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[7]->Fill(fabs(theta_max));
	  //h_hist_vec[7]->Fill(fabs(theta1-theta2)*TMath::RadToDeg());
	  //h_hist_vec[7]->Fill((tempGenJet1.E()+tempGenJet2.E())/E_trueAll);
	}
	if(boson_is_W && (tempGenJet1+tempGenJet2).M()>75.){
	  h_hist_vec[7]->Fill(DeltaPhi(phi_MET_gen,(tempGenJet1+tempGenJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[7]->Fill(theta_max);
	  //h_hist_vec[7]->Fill(fabs(theta1-theta2)*TMath::RadToDeg());
	  //h_hist_vec[7]->Fill((tempGenJet1.E()+tempGenJet2.E())/(E_trueAll-lep1E));
	}
      }
    }
    h_hist_vec[5]->Fill(tempMEVecQuarks.M());
    h_hist_vec[10]->Fill(tempMEVecQuarks.Theta()*TMath::RadToDeg());

    TLorentzVector tempRecoJet1(0,0,0,0);
    TLorentzVector tempRecoJet2(0,0,0,0);

    if(recoJetR07E->size()>1){
      tempRecoJet1.SetPxPyPzE((*recoJetR07Px)[0],(*recoJetR07Py)[0],(*recoJetR07Pz)[0],(*recoJetR07E)[0]);
      tempRecoJet2.SetPxPyPzE((*recoJetR07Px)[1],(*recoJetR07Py)[1],(*recoJetR07Pz)[1],(*recoJetR07E)[1]);
      if(fabs(tempRecoJet1.CosTheta())<m_costheta_max && fabs(tempRecoJet2.CosTheta())<m_costheta_max){
	h_hist_vec[1]->Fill(tempRecoJet1.E()+tempRecoJet2.E());
	h_hist_vec[3]->Fill((tempRecoJet1+tempRecoJet2).M());
	if((tempRecoJet1+tempRecoJet2).M()<60){	  
	  h_hist_vec[8]->Fill(DeltaPhi(phi_MET_reco,(tempRecoJet1+tempRecoJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[8]->Fill(fabs(theta_min));
	  //h_hist_vec[8]->Fill(DeltaPhi(phi1,phi2)*TMath::RadToDeg());
	  //if(boson_is_Z){
	    //h_hist_vec[8]->Fill((tempRecoJet1.E()+tempRecoJet2.E())/E_totPFO);
	    //}else{
	    //h_hist_vec[8]->Fill((tempRecoJet1.E()+tempRecoJet2.E())/(E_totPFO-lep1E));
	    //}
	}
	if(boson_is_Z && (tempRecoJet1+tempRecoJet2).M()>85.){
	  h_hist_vec[9]->Fill(DeltaPhi(phi_MET_reco,(tempRecoJet1+tempRecoJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[9]->Fill(fabs(theta_min));
	  //h_hist_vec[9]->Fill(DeltaPhi(phi1,phi2)*TMath::RadToDeg());
	  //h_hist_vec[9]->Fill((tempRecoJet1.E()+tempRecoJet2.E())/E_totPFO);
	}
	if(boson_is_W && (tempRecoJet1+tempRecoJet2).M()>75.){
	  h_hist_vec[9]->Fill(DeltaPhi(phi_MET_reco,(tempRecoJet1+tempRecoJet2).Phi())*TMath::RadToDeg());
	  //h_hist_vec[9]->Fill(fabs(theta_min));
	  //h_hist_vec[9]->Fill(DeltaPhi(phi1,phi2)*TMath::RadToDeg());
	  //h_hist_vec[9]->Fill((tempRecoJet1.E()+tempRecoJet2.E())/(E_totPFO-lep1E));
	}
      }//else{
	//std::cout<<"have at least one forward jet "<<fabs(tempRecoJet1.CosTheta())<<"/"<<fabs(tempRecoJet2.CosTheta())<<std::endl;
      //}      
    }
  }
}

  


void JERMassPlotsFull(){

  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

  //const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles/ZZ_WW_wO_181101_mass_histos_noDAngleMin_VLC7_cosjet09_no_cut_E12_over_E_tot_True_09_50bins_mass50_120_wO_380.root";
  const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles181123/ZZ_WW_wO_181011_newVLC_mass_histos_noDAngleMin_VLC7_cosjet09_no_cut_E12_over_E_tot_True_09_50bins_mass50_120_wO_380.root";

  //ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  //pMinimizer->SetMaxFunctionCalls(10000000);
  //pMinimizer->SetMaxIterations(1000000);

  ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  pMinimizer->SetMaxFunctionCalls(10000000);
  pMinimizer->SetMaxIterations(1000000);
  pMinimizer->SetTolerance(0.001);

  TVirtualFitter::SetDefaultFitter("Minuit2");
  TVirtualFitter::SetMaxIterations(20000);
  //TVirtualFitter::SetPrecision(0.01);

  string label_legend= "#gamma, E_{true}=1 GeV";

  const unsigned int nRegionBinsAngles(4);
  float pRegionBinEdgesAngles[nRegionBinsAngles + 1] = {0, 0.65, 0.8, 0.925, 0.975};

  const unsigned int nRegionBins(13);
  float pRegionBinEdges[nRegionBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};

  float resolutionTH1(0.f), resolutionErrorTH1(0.f),meanTH1(0.f),meanErrorTH1(0.f);
  resolutionTH1 = 0.0;
  resolutionErrorTH1 = 0.0;

  meanTH1 = 0.0;
  meanErrorTH1 = 0.0;
 


  TFile* file_CLIC_WW_250_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW250_11560_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11560_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_WW_250_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW250_11561_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  ///eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11561_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  //TFile* file_CLIC_WW_250_CT_wO_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11561_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_250_CT_wO_Selected = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11561_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_250_CT_wO_allPFOs = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11561_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");

  TFile* file_CLIC_WW_500_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW500_11564_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //"/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW500_11564_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_WW_500_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW500_11565_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW500_11565_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  //TFile* file_CLIC_WW_500_CT_wO_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW500_11565_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_500_CT_wO_Selected = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW500_11565_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_500_CT_wO_allPFOs = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW500_11565_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");

  TFile* file_CLIC_WW_1000_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW1000_11568_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW1000_11568_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_WW_1000_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW1000_11569_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW1000_11569_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  //TFile* file_CLIC_WW_1000_CT_wO_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW1000_11569_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_1000_CT_wO_Selected = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW1000_11569_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_1000_CT_wO_allPFOs = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW1000_11569_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");

  TFile* file_CLIC_WW_2000_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW2000_11572_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW2000_11572_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_WW_2000_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW2000_11573_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW2000_11573_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  //TFile* file_CLIC_WW_2000_CT_wO_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW2000_11573_3TeVOverlay_CLIC_o3_v14_CT_SelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_2000_CT_wO_Selected = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW2000_11573_3TeVOverlay_CLIC_o3_v14_CT_LooseSelectedPandoraPFOs.root");
  //TFile* file_CLIC_WW_2000_CT_wO_allPFOs = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW2000_11573_3TeVOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  
  TFile* file_CLIC_ZZ_250_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ250_11576_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ250_11576_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_250_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ250_11577_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ250_11577_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_500_CT    = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ500_11580_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ500_11580_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_500_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ500_11581_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ500_11581_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  
  TFile* file_CLIC_ZZ_1000_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ1000_11584_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ1000_11584_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_1000_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ1000_11585_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ1000_11585_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_2000_CT = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ2000_11588_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ2000_11588_noOverlay_CLIC_o3_v14_CT_PandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_2000_CT_wO = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ2000_11589_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ2000_11589_3TeVOverlay_CLIC_o3_v14_CT_TightSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  

  TFile* file_CLIC_WW_250_CT_wO_380_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_WW250_11608_380GeVOverlay_CLIC_o3_v14_CT_LE_LooseSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_WW250_11608_380GeVOverlay_CLIC_o3_v14_CT_LE_LooseSelectedPandoraPFOs_largeCones_DR12_DR15.root");
  TFile* file_CLIC_ZZ_250_CT_wO_380_Loose = TFile::Open("/eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62_newVLC/JetStudy_ZZ250_11610_380GeVOverlay_CLIC_o3_v14_CT_LE_LooseSelectedPandoraPFOs.root");
  //eos/user/w/weberma2/data/JetAnalyzerFiles/181011_gcc62/JetStudy_ZZ250_11610_380GeVOverlay_CLIC_o3_v14_CT_LE_LooseSelectedPandoraPFOs_largeCones_DR12_DR15.root");

  
  TFile* file_histogram=new TFile(final_histo_name,"recreate");

  unsigned int n_bins100=75;
  int n_binsJER = 50;
  float lim_E_rel_low=0;
  float lim_E_rel_high=3.00;

  float lim_angleDeg=10.00;

  int firstbin=59; 


  TH1F* h_ZZ_250_CT_genjet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_genjet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_recojet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_recojet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_genjet_M_jj =  new TH1F("h_ZZ_250_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_recojet_M_jj =  new TH1F("h_ZZ_250_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_M_ZBoson =  new TH1F("h_ZZ_250_CT_parton_M_ZBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_M_qq =  new TH1F("h_ZZ_250_CT_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_Z_qq_theta =  new TH1F("h_ZZ_250_CT_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_250_CT;
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_genjet_E1_plus_E2);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_recojet_E1_plus_E2);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_genjet_M_jj);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_recojet_M_jj);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_M_ZBoson);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_M_qq);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT.push_back(h_ZZ_250_CT_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_250_CT.size();i++){
    hist_vec_ZZ_250_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_250_CT,hist_vec_ZZ_250_CT);

  TH1F* h_ZZ_250_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_wO_genjet_M_jj =  new TH1F("h_ZZ_250_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_wO_recojet_M_jj =  new TH1F("h_ZZ_250_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_wO_M_ZBoson =  new TH1F("h_ZZ_250_CT_wO_parton_M_ZBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_wO_M_qq =  new TH1F("h_ZZ_250_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_wO_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_wO_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_Z_qq_theta =  new TH1F("h_ZZ_250_CT_wO_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_250_CT_wO;
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_genjet_E1_plus_E2);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_recojet_E1_plus_E2);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_genjet_M_jj);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_recojet_M_jj);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_M_ZBoson);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_M_qq);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT_wO.push_back(h_ZZ_250_CT_wO_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_250_CT_wO.size();i++){
    hist_vec_ZZ_250_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_250_CT_wO,hist_vec_ZZ_250_CT_wO);

  TH1F* h_ZZ_250_CT_wO_380_genjet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_wO_380_genjet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_wO_380_recojet_E1_plus_E2 =  new TH1F("h_ZZ_250_CT_wO_380_recojet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_ZZ_250_CT_wO_380_genjet_M_jj =  new TH1F("h_ZZ_250_CT_wO_380_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_wO_380_recojet_M_jj =  new TH1F("h_ZZ_250_CT_wO_380_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_250_CT_wO_380_M_ZBoson =  new TH1F("h_ZZ_250_CT_wO_380_parton_M_ZBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_wO_380_M_qq =  new TH1F("h_ZZ_250_CT_wO_380_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_ZZ_250_CT_wO_380_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_wO_380_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_380_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_wO_380_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_380_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_250_CT_wO_380_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_380_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_250_CT_wO_380_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_250_CT_wO_380_Z_qq_theta =  new TH1F("h_ZZ_250_CT_wO_380_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_250_CT_wO_380;
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_genjet_E1_plus_E2);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_recojet_E1_plus_E2);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_genjet_M_jj);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_recojet_M_jj);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_M_ZBoson);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_M_qq);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_250_CT_wO_380.push_back(h_ZZ_250_CT_wO_380_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_250_CT_wO_380.size();i++){
    hist_vec_ZZ_250_CT_wO_380[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_250_CT_wO_380_Loose,hist_vec_ZZ_250_CT_wO_380);

  std::cout<<"ZZ 250 done"<<std::endl;
 
  TH1F* h_ZZ_500_CT_genjet_E1_plus_E2 =  new TH1F("h_ZZ_500_CT_genjet_E1_plus_E2","", n_binsJER, 50,500);
  TH1F* h_ZZ_500_CT_recojet_E1_plus_E2 =  new TH1F("h_ZZ_500_CT_recojet_E1_plus_E2","", n_binsJER, 50,500);
  TH1F* h_ZZ_500_CT_genjet_M_jj =  new TH1F("h_ZZ_500_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_500_CT_recojet_M_jj =  new TH1F("h_ZZ_500_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_500_CT_M_ZBoson =  new TH1F("h_ZZ_500_CT_parton_M_ZBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_ZZ_500_CT_M_qq =  new TH1F("h_ZZ_500_CT_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_ZZ_500_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_500_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_500_CT_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_500_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_500_CT_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_Z_qq_theta =  new TH1F("h_ZZ_500_CT_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_500_CT;
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_genjet_E1_plus_E2);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_recojet_E1_plus_E2);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_genjet_M_jj);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_recojet_M_jj);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_M_ZBoson);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_M_qq);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_500_CT.push_back(h_ZZ_500_CT_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_500_CT.size();i++){
    hist_vec_ZZ_500_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_500_CT,hist_vec_ZZ_500_CT);

  TH1F* h_ZZ_500_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_ZZ_500_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,500);
  TH1F* h_ZZ_500_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_ZZ_500_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,500);
  TH1F* h_ZZ_500_CT_wO_genjet_M_jj =  new TH1F("h_ZZ_500_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_500_CT_wO_recojet_M_jj =  new TH1F("h_ZZ_500_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_500_CT_wO_M_ZBoson =  new TH1F("h_ZZ_500_CT_wO_parton_M_ZBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_ZZ_500_CT_wO_M_qq =  new TH1F("h_ZZ_500_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_ZZ_500_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_500_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_wO_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_500_CT_wO_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_500_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_wO_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_500_CT_wO_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_500_CT_wO_Z_qq_theta =  new TH1F("h_ZZ_500_CT_wO_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_500_CT_wO;
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_genjet_E1_plus_E2);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_recojet_E1_plus_E2);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_genjet_M_jj);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_recojet_M_jj);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_M_ZBoson);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_M_qq);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_500_CT_wO.push_back(h_ZZ_500_CT_wO_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_500_CT_wO.size();i++){
    hist_vec_ZZ_500_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_500_CT_wO,hist_vec_ZZ_500_CT_wO);
  std::cout<<"ZZ 500 done"<<std::endl;
 
  TH1F* h_ZZ_1000_CT_genjet_E1_plus_E2 =  new TH1F("h_ZZ_1000_CT_genjet_E1_plus_E2","", n_binsJER, 50,1000);
  TH1F* h_ZZ_1000_CT_recojet_E1_plus_E2 =  new TH1F("h_ZZ_1000_CT_recojet_E1_plus_E2","", n_binsJER, 50,1000);
  TH1F* h_ZZ_1000_CT_genjet_M_jj =  new TH1F("h_ZZ_1000_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_1000_CT_recojet_M_jj =  new TH1F("h_ZZ_1000_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_1000_CT_M_ZBoson =  new TH1F("h_ZZ_1000_CT_parton_M_ZBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_ZZ_1000_CT_M_qq =  new TH1F("h_ZZ_1000_CT_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_ZZ_1000_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_1000_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_1000_CT_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_1000_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_1000_CT_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_Z_qq_theta =  new TH1F("h_ZZ_1000_CT_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_1000_CT;
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_genjet_E1_plus_E2);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_recojet_E1_plus_E2);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_genjet_M_jj);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_recojet_M_jj);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_M_ZBoson);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_M_qq);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_1000_CT.push_back(h_ZZ_1000_CT_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_1000_CT.size();i++){
    hist_vec_ZZ_1000_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_1000_CT,hist_vec_ZZ_1000_CT);

  TH1F* h_ZZ_1000_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_ZZ_1000_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,1000);
  TH1F* h_ZZ_1000_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_ZZ_1000_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,1000);
  TH1F* h_ZZ_1000_CT_wO_genjet_M_jj =  new TH1F("h_ZZ_1000_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_1000_CT_wO_recojet_M_jj =  new TH1F("h_ZZ_1000_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_1000_CT_wO_M_ZBoson =  new TH1F("h_ZZ_1000_CT_wO_parton_M_ZBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_ZZ_1000_CT_wO_M_qq =  new TH1F("h_ZZ_1000_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_ZZ_1000_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_1000_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_wO_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_1000_CT_wO_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_1000_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_wO_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_1000_CT_wO_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_1000_CT_wO_Z_qq_theta =  new TH1F("h_ZZ_1000_CT_wO_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_1000_CT_wO;
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_genjet_E1_plus_E2);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_recojet_E1_plus_E2);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_genjet_M_jj);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_recojet_M_jj);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_M_ZBoson);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_M_qq);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_1000_CT_wO.push_back(h_ZZ_1000_CT_wO_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_1000_CT_wO.size();i++){
    hist_vec_ZZ_1000_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_1000_CT_wO,hist_vec_ZZ_1000_CT_wO);
  std::cout<<"ZZ 1000 done"<<std::endl;
 
  TH1F* h_ZZ_2000_CT_genjet_E1_plus_E2 =  new TH1F("h_ZZ_2000_CT_genjet_E1_plus_E2","", n_binsJER, 50,2000);
  TH1F* h_ZZ_2000_CT_recojet_E1_plus_E2 =  new TH1F("h_ZZ_2000_CT_recojet_E1_plus_E2","", n_binsJER, 50,2000);
  TH1F* h_ZZ_2000_CT_genjet_M_jj =  new TH1F("h_ZZ_2000_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_2000_CT_recojet_M_jj =  new TH1F("h_ZZ_2000_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_2000_CT_M_ZBoson =  new TH1F("h_ZZ_2000_CT_parton_M_ZBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_ZZ_2000_CT_M_qq =  new TH1F("h_ZZ_2000_CT_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_ZZ_2000_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_2000_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_2000_CT_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_2000_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_2000_CT_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_Z_qq_theta =  new TH1F("h_ZZ_2000_CT_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_2000_CT;
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_genjet_E1_plus_E2);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_recojet_E1_plus_E2);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_genjet_M_jj);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_recojet_M_jj);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_M_ZBoson);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_M_qq);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_2000_CT.push_back(h_ZZ_2000_CT_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_2000_CT.size();i++){
    hist_vec_ZZ_2000_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_2000_CT,hist_vec_ZZ_2000_CT);

  TH1F* h_ZZ_2000_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_ZZ_2000_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,2000);
  TH1F* h_ZZ_2000_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_ZZ_2000_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,2000);
  TH1F* h_ZZ_2000_CT_wO_genjet_M_jj =  new TH1F("h_ZZ_2000_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_2000_CT_wO_recojet_M_jj =  new TH1F("h_ZZ_2000_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_ZZ_2000_CT_wO_M_ZBoson =  new TH1F("h_ZZ_2000_CT_wO_parton_M_ZBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_ZZ_2000_CT_wO_M_qq =  new TH1F("h_ZZ_2000_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_ZZ_2000_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_2000_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_wO_genjet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_2000_CT_wO_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_ZZ_2000_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_wO_recojet_dPhijjMET_Mjj85  =  new TH1F("h_ZZ_2000_CT_wO_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_ZZ_2000_CT_wO_Z_qq_theta =  new TH1F("h_ZZ_2000_CT_wO_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_ZZ_2000_CT_wO;
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_genjet_E1_plus_E2);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_recojet_E1_plus_E2);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_genjet_M_jj);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_recojet_M_jj);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_M_ZBoson);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_M_qq);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_genjet_dPhijjMET_Mjj85);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_recojet_dPhijjMET_Mjj85);
  hist_vec_ZZ_2000_CT_wO.push_back(h_ZZ_2000_CT_wO_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_ZZ_2000_CT_wO.size();i++){
    hist_vec_ZZ_2000_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_ZZ_2000_CT_wO,hist_vec_ZZ_2000_CT_wO);
  std::cout<<"ZZ 2000 done"<<std::endl;
  TH1F* h_WW_250_CT_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_genjet_M_jj =  new TH1F("h_WW_250_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_recojet_M_jj =  new TH1F("h_WW_250_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_M_WBoson =  new TH1F("h_WW_250_CT_parton_M_WBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_M_qq =  new TH1F("h_WW_250_CT_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_W_qq_theta =  new TH1F("h_WW_250_CT_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT;
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_genjet_E1_plus_E2);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_recojet_E1_plus_E2);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_genjet_M_jj);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_recojet_M_jj);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_M_WBoson);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_M_qq);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT.push_back(h_WW_250_CT_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT.size();i++){
    hist_vec_WW_250_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT,hist_vec_WW_250_CT);

  TH1F* h_WW_250_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_genjet_M_jj =  new TH1F("h_WW_250_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_recojet_M_jj =  new TH1F("h_WW_250_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_M_WBoson =  new TH1F("h_WW_250_CT_wO_parton_M_WBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_M_qq =  new TH1F("h_WW_250_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_W_qq_theta =  new TH1F("h_WW_250_CT_wO_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT_wO;
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_genjet_E1_plus_E2);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_recojet_E1_plus_E2);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_genjet_M_jj);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_recojet_M_jj);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_M_WBoson);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_M_qq);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO.push_back(h_WW_250_CT_wO_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT_wO.size();i++){
    hist_vec_WW_250_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT_wO,hist_vec_WW_250_CT_wO);
  /*
  TH1F* h_WW_250_CT_wO_Loose_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_Loose_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_Loose_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_Loose_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_Loose_genjet_M_jj =  new TH1F("h_WW_250_CT_wO_Loose_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_Loose_recojet_M_jj =  new TH1F("h_WW_250_CT_wO_Loose_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_Loose_M_WBoson =  new TH1F("h_WW_250_CT_wO_Loose_parton_M_WBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_Loose_M_qq =  new TH1F("h_WW_250_CT_wO_Loose_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_Loose_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_Loose_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Loose_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_Loose_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Loose_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_Loose_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Loose_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_Loose_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Loose_W_qq_theta =  new TH1F("h_WW_250_CT_wO_Loose_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT_wO_Loose;
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_genjet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_recojet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_genjet_M_jj);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_recojet_M_jj);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_M_WBoson);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_M_qq);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_Loose.push_back(h_WW_250_CT_wO_Loose_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT_wO_Loose.size();i++){
    hist_vec_WW_250_CT_wO_Loose[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT_wO_Loose,hist_vec_WW_250_CT_wO_Loose);

  TH1F* h_WW_250_CT_wO_Selected_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_Selected_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_Selected_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_Selected_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_Selected_genjet_M_jj =  new TH1F("h_WW_250_CT_wO_Selected_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_Selected_recojet_M_jj =  new TH1F("h_WW_250_CT_wO_Selected_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_Selected_M_WBoson =  new TH1F("h_WW_250_CT_wO_Selected_parton_M_WBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_Selected_M_qq =  new TH1F("h_WW_250_CT_wO_Selected_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_Selected_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_Selected_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Selected_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_Selected_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Selected_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_Selected_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Selected_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_Selected_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_Selected_W_qq_theta =  new TH1F("h_WW_250_CT_wO_Selected_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT_wO_Selected;
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_genjet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_recojet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_genjet_M_jj);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_recojet_M_jj);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_M_WBoson);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_M_qq);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_Selected.push_back(h_WW_250_CT_wO_Selected_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT_wO_Selected.size();i++){
    hist_vec_WW_250_CT_wO_Selected[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT_wO_Selected,hist_vec_WW_250_CT_wO_Selected);

  TH1F* h_WW_250_CT_wO_allPFOs_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_allPFOs_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_allPFOs_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_allPFOs_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_250_CT_wO_allPFOs_genjet_M_jj =  new TH1F("h_WW_250_CT_wO_allPFOs_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_allPFOs_recojet_M_jj =  new TH1F("h_WW_250_CT_wO_allPFOs_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_allPFOs_M_WBoson =  new TH1F("h_WW_250_CT_wO_allPFOs_parton_M_WBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_allPFOs_M_qq =  new TH1F("h_WW_250_CT_wO_allPFOs_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_allPFOs_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_allPFOs_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_allPFOs_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_250_CT_wO_allPFOs_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_allPFOs_W_qq_theta =  new TH1F("h_WW_250_CT_wO_allPFOs_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT_wO_allPFOs;
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_genjet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_recojet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_genjet_M_jj);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_recojet_M_jj);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_M_WBoson);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_M_qq);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_250_CT_wO_allPFOs.push_back(h_WW_250_CT_wO_allPFOs_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT_wO_allPFOs.size();i++){
    hist_vec_WW_250_CT_wO_allPFOs[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT_wO_allPFOs,hist_vec_WW_250_CT_wO_allPFOs);
  */
  TH1F* h_WW_250_CT_wO_380_genjet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_380_genjet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_WW_250_CT_wO_380_recojet_E1_plus_E2 =  new TH1F("h_WW_250_CT_wO_380_recojet_E1_plus_E2","", n_binsJER, 50,250);
  TH1F* h_WW_250_CT_wO_380_genjet_M_jj =  new TH1F("h_WW_250_CT_wO_380_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_380_recojet_M_jj =  new TH1F("h_WW_250_CT_wO_380_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_250_CT_wO_380_M_ZBoson =  new TH1F("h_WW_250_CT_wO_380_parton_M_ZBoson","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_380_M_qq =  new TH1F("h_WW_250_CT_wO_380_parton_M_qq","", 10*n_binsJER, 20.,250.);
  TH1F* h_WW_250_CT_wO_380_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_380_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_380_genjet_dPhijjMET_Mjj85  =  new TH1F("h_WW_250_CT_wO_380_genjet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_380_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_250_CT_wO_380_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_380_recojet_dPhijjMET_Mjj85  =  new TH1F("h_WW_250_CT_wO_380_recojet_dphijjMET_Mjj85","", n_binsJER, 50.,180.);
  TH1F* h_WW_250_CT_wO_380_Z_qq_theta =  new TH1F("h_WW_250_CT_wO_380_Z_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_250_CT_wO_380;
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_genjet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_recojet_E1_plus_E2);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_genjet_M_jj);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_recojet_M_jj);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_M_ZBoson);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_M_qq);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_genjet_dPhijjMET_Mjj85);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_recojet_dPhijjMET_Mjj85);
  hist_vec_WW_250_CT_wO_380.push_back(h_WW_250_CT_wO_380_Z_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_250_CT_wO_380.size();i++){
    hist_vec_WW_250_CT_wO_380[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_250_CT_wO_380_Loose,hist_vec_WW_250_CT_wO_380);

  std::cout<<"WW_250 done"<<std::endl;
  TH1F* h_WW_500_CT_genjet_E1_plus_E2 =  new TH1F("h_WW_500_CT_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_recojet_E1_plus_E2 =  new TH1F("h_WW_500_CT_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_genjet_M_jj =  new TH1F("h_WW_500_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_recojet_M_jj =  new TH1F("h_WW_500_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_M_WBoson =  new TH1F("h_WW_500_CT_parton_M_WBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_M_qq =  new TH1F("h_WW_500_CT_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_W_qq_theta =  new TH1F("h_WW_500_CT_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_500_CT;
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_genjet_E1_plus_E2);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_recojet_E1_plus_E2);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_genjet_M_jj);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_recojet_M_jj);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_M_WBoson);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_M_qq);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT.push_back(h_WW_500_CT_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_500_CT.size();i++){
    hist_vec_WW_500_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_500_CT,hist_vec_WW_500_CT);

  TH1F* h_WW_500_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_genjet_M_jj =  new TH1F("h_WW_500_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_recojet_M_jj =  new TH1F("h_WW_500_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_M_WBoson =  new TH1F("h_WW_500_CT_wO_parton_M_WBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_M_qq =  new TH1F("h_WW_500_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_W_qq_theta =  new TH1F("h_WW_500_CT_wO_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_500_CT_wO;
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_genjet_E1_plus_E2);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_recojet_E1_plus_E2);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_genjet_M_jj);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_recojet_M_jj);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_M_WBoson);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_M_qq);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO.push_back(h_WW_500_CT_wO_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_500_CT_wO.size();i++){
    hist_vec_WW_500_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_500_CT_wO,hist_vec_WW_500_CT_wO);
  /*
  TH1F* h_WW_500_CT_wO_Loose_genjet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_Loose_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_Loose_recojet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_Loose_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_Loose_genjet_M_jj =  new TH1F("h_WW_500_CT_wO_Loose_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_Loose_recojet_M_jj =  new TH1F("h_WW_500_CT_wO_Loose_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_Loose_M_WBoson =  new TH1F("h_WW_500_CT_wO_Loose_parton_M_WBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_Loose_M_qq =  new TH1F("h_WW_500_CT_wO_Loose_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_Loose_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_Loose_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Loose_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_Loose_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Loose_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_Loose_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Loose_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_Loose_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Loose_W_qq_theta =  new TH1F("h_WW_500_CT_wO_Loose_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_500_CT_wO_Loose;
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_genjet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_recojet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_genjet_M_jj);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_recojet_M_jj);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_M_WBoson);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_M_qq);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_Loose.push_back(h_WW_500_CT_wO_Loose_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_500_CT_wO_Loose.size();i++){
    hist_vec_WW_500_CT_wO_Loose[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_500_CT_wO_Loose,hist_vec_WW_500_CT_wO_Loose);

  TH1F* h_WW_500_CT_wO_Selected_genjet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_Selected_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_Selected_recojet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_Selected_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_Selected_genjet_M_jj =  new TH1F("h_WW_500_CT_wO_Selected_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_Selected_recojet_M_jj =  new TH1F("h_WW_500_CT_wO_Selected_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_Selected_M_WBoson =  new TH1F("h_WW_500_CT_wO_Selected_parton_M_WBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_Selected_M_qq =  new TH1F("h_WW_500_CT_wO_Selected_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_Selected_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_Selected_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Selected_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_Selected_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Selected_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_Selected_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Selected_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_Selected_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_Selected_W_qq_theta =  new TH1F("h_WW_500_CT_wO_Selected_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_500_CT_wO_Selected;
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_genjet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_recojet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_genjet_M_jj);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_recojet_M_jj);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_M_WBoson);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_M_qq);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_Selected.push_back(h_WW_500_CT_wO_Selected_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_500_CT_wO_Selected.size();i++){
    hist_vec_WW_500_CT_wO_Selected[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_500_CT_wO_Selected,hist_vec_WW_500_CT_wO_Selected);

  TH1F* h_WW_500_CT_wO_allPFOs_genjet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_allPFOs_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_allPFOs_recojet_E1_plus_E2 =  new TH1F("h_WW_500_CT_wO_allPFOs_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_500_CT_wO_allPFOs_genjet_M_jj =  new TH1F("h_WW_500_CT_wO_allPFOs_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_allPFOs_recojet_M_jj =  new TH1F("h_WW_500_CT_wO_allPFOs_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_500_CT_wO_allPFOs_M_WBoson =  new TH1F("h_WW_500_CT_wO_allPFOs_parton_M_WBoson","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_allPFOs_M_qq =  new TH1F("h_WW_500_CT_wO_allPFOs_parton_M_qq","", 10*n_binsJER, 20.,500.);
  TH1F* h_WW_500_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_allPFOs_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_allPFOs_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_500_CT_wO_allPFOs_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_500_CT_wO_allPFOs_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_500_CT_wO_allPFOs_W_qq_theta =  new TH1F("h_WW_500_CT_wO_allPFOs_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_500_CT_wO_allPFOs;
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_genjet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_recojet_E1_plus_E2);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_genjet_M_jj);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_recojet_M_jj);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_M_WBoson);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_M_qq);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_500_CT_wO_allPFOs.push_back(h_WW_500_CT_wO_allPFOs_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_500_CT_wO_allPFOs.size();i++){
    hist_vec_WW_500_CT_wO_allPFOs[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_500_CT_wO_allPFOs,hist_vec_WW_500_CT_wO_allPFOs);
  */
  std::cout<<"WW_500 done"<<std::endl;
  TH1F* h_WW_1000_CT_genjet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_recojet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_genjet_M_jj =  new TH1F("h_WW_1000_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_recojet_M_jj =  new TH1F("h_WW_1000_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_M_WBoson =  new TH1F("h_WW_1000_CT_parton_M_WBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_M_qq =  new TH1F("h_WW_1000_CT_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_W_qq_theta =  new TH1F("h_WW_1000_CT_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_1000_CT;
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_genjet_E1_plus_E2);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_recojet_E1_plus_E2);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_genjet_M_jj);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_recojet_M_jj);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_M_WBoson);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_M_qq);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT.push_back(h_WW_1000_CT_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_1000_CT.size();i++){
    hist_vec_WW_1000_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_1000_CT,hist_vec_WW_1000_CT);

  TH1F* h_WW_1000_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_genjet_M_jj =  new TH1F("h_WW_1000_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_recojet_M_jj =  new TH1F("h_WW_1000_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_M_WBoson =  new TH1F("h_WW_1000_CT_wO_parton_M_WBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_M_qq =  new TH1F("h_WW_1000_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_W_qq_theta =  new TH1F("h_WW_1000_CT_wO_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_1000_CT_wO;
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_genjet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_recojet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_genjet_M_jj);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_recojet_M_jj);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_M_WBoson);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_M_qq);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO.push_back(h_WW_1000_CT_wO_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_1000_CT_wO.size();i++){
    hist_vec_WW_1000_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_1000_CT_wO,hist_vec_WW_1000_CT_wO);
  /*
  TH1F* h_WW_1000_CT_wO_Loose_genjet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_Loose_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_Loose_recojet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_Loose_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_Loose_genjet_M_jj =  new TH1F("h_WW_1000_CT_wO_Loose_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_Loose_recojet_M_jj =  new TH1F("h_WW_1000_CT_wO_Loose_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_Loose_M_WBoson =  new TH1F("h_WW_1000_CT_wO_Loose_parton_M_WBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_Loose_M_qq =  new TH1F("h_WW_1000_CT_wO_Loose_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_Loose_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_Loose_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Loose_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_Loose_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Loose_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_Loose_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Loose_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_Loose_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Loose_W_qq_theta =  new TH1F("h_WW_1000_CT_wO_Loose_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_1000_CT_wO_Loose;
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_genjet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_recojet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_genjet_M_jj);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_recojet_M_jj);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_M_WBoson);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_M_qq);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_Loose.push_back(h_WW_1000_CT_wO_Loose_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_1000_CT_wO_Loose.size();i++){
    hist_vec_WW_1000_CT_wO_Loose[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_1000_CT_wO_Loose,hist_vec_WW_1000_CT_wO_Loose);

  TH1F* h_WW_1000_CT_wO_Selected_genjet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_Selected_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_Selected_recojet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_Selected_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_Selected_genjet_M_jj =  new TH1F("h_WW_1000_CT_wO_Selected_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_Selected_recojet_M_jj =  new TH1F("h_WW_1000_CT_wO_Selected_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_Selected_M_WBoson =  new TH1F("h_WW_1000_CT_wO_Selected_parton_M_WBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_Selected_M_qq =  new TH1F("h_WW_1000_CT_wO_Selected_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_Selected_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_Selected_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Selected_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_Selected_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Selected_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_Selected_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Selected_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_Selected_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_Selected_W_qq_theta =  new TH1F("h_WW_1000_CT_wO_Selected_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_1000_CT_wO_Selected;
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_genjet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_recojet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_genjet_M_jj);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_recojet_M_jj);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_M_WBoson);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_M_qq);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_Selected.push_back(h_WW_1000_CT_wO_Selected_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_1000_CT_wO_Selected.size();i++){
    hist_vec_WW_1000_CT_wO_Selected[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_1000_CT_wO_Selected,hist_vec_WW_1000_CT_wO_Selected);

  TH1F* h_WW_1000_CT_wO_allPFOs_genjet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_allPFOs_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_allPFOs_recojet_E1_plus_E2 =  new TH1F("h_WW_1000_CT_wO_allPFOs_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_1000_CT_wO_allPFOs_genjet_M_jj =  new TH1F("h_WW_1000_CT_wO_allPFOs_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_allPFOs_recojet_M_jj =  new TH1F("h_WW_1000_CT_wO_allPFOs_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_1000_CT_wO_allPFOs_M_WBoson =  new TH1F("h_WW_1000_CT_wO_allPFOs_parton_M_WBoson","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_allPFOs_M_qq =  new TH1F("h_WW_1000_CT_wO_allPFOs_parton_M_qq","", 10*n_binsJER, 20.,1000.);
  TH1F* h_WW_1000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_allPFOs_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_allPFOs_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_1000_CT_wO_allPFOs_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_1000_CT_wO_allPFOs_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_1000_CT_wO_allPFOs_W_qq_theta =  new TH1F("h_WW_1000_CT_wO_allPFOs_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_1000_CT_wO_allPFOs;
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_genjet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_recojet_E1_plus_E2);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_genjet_M_jj);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_recojet_M_jj);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_M_WBoson);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_M_qq);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_1000_CT_wO_allPFOs.push_back(h_WW_1000_CT_wO_allPFOs_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_1000_CT_wO_allPFOs.size();i++){
    hist_vec_WW_1000_CT_wO_allPFOs[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_1000_CT_wO_allPFOs,hist_vec_WW_1000_CT_wO_allPFOs);
  */
  std::cout<<"WW_1000 done"<<std::endl;
  TH1F* h_WW_2000_CT_genjet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_recojet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_genjet_M_jj =  new TH1F("h_WW_2000_CT_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_recojet_M_jj =  new TH1F("h_WW_2000_CT_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_M_WBoson =  new TH1F("h_WW_2000_CT_parton_M_WBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_M_qq =  new TH1F("h_WW_2000_CT_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_W_qq_theta =  new TH1F("h_WW_2000_CT_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_2000_CT;
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_genjet_E1_plus_E2);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_recojet_E1_plus_E2);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_genjet_M_jj);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_recojet_M_jj);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_M_WBoson);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_M_qq);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT.push_back(h_WW_2000_CT_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_2000_CT.size();i++){
    hist_vec_WW_2000_CT[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_2000_CT,hist_vec_WW_2000_CT);

  TH1F* h_WW_2000_CT_wO_genjet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_recojet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_genjet_M_jj =  new TH1F("h_WW_2000_CT_wO_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_recojet_M_jj =  new TH1F("h_WW_2000_CT_wO_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_M_WBoson =  new TH1F("h_WW_2000_CT_wO_parton_M_WBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_M_qq =  new TH1F("h_WW_2000_CT_wO_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_W_qq_theta =  new TH1F("h_WW_2000_CT_wO_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_2000_CT_wO;
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_genjet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_recojet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_genjet_M_jj);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_recojet_M_jj);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_M_WBoson);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_M_qq);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO.push_back(h_WW_2000_CT_wO_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_2000_CT_wO.size();i++){
    hist_vec_WW_2000_CT_wO[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_2000_CT_wO,hist_vec_WW_2000_CT_wO);
  /*
  TH1F* h_WW_2000_CT_wO_Loose_genjet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_Loose_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_Loose_recojet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_Loose_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_Loose_genjet_M_jj =  new TH1F("h_WW_2000_CT_wO_Loose_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_Loose_recojet_M_jj =  new TH1F("h_WW_2000_CT_wO_Loose_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_Loose_M_WBoson =  new TH1F("h_WW_2000_CT_wO_Loose_parton_M_WBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_Loose_M_qq =  new TH1F("h_WW_2000_CT_wO_Loose_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_Loose_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_Loose_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Loose_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_Loose_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Loose_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_Loose_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Loose_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_Loose_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Loose_W_qq_theta =  new TH1F("h_WW_2000_CT_wO_Loose_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_2000_CT_wO_Loose;
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_genjet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_recojet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_genjet_M_jj);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_recojet_M_jj);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_M_WBoson);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_M_qq);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_Loose.push_back(h_WW_2000_CT_wO_Loose_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_2000_CT_wO_Loose.size();i++){
    hist_vec_WW_2000_CT_wO_Loose[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_2000_CT_wO_Loose,hist_vec_WW_2000_CT_wO_Loose);

  TH1F* h_WW_2000_CT_wO_Selected_genjet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_Selected_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_Selected_recojet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_Selected_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_Selected_genjet_M_jj =  new TH1F("h_WW_2000_CT_wO_Selected_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_Selected_recojet_M_jj =  new TH1F("h_WW_2000_CT_wO_Selected_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_Selected_M_WBoson =  new TH1F("h_WW_2000_CT_wO_Selected_parton_M_WBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_Selected_M_qq =  new TH1F("h_WW_2000_CT_wO_Selected_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_Selected_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_Selected_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Selected_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_Selected_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Selected_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_Selected_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Selected_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_Selected_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_Selected_W_qq_theta =  new TH1F("h_WW_2000_CT_wO_Selected_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_2000_CT_wO_Selected;
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_genjet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_recojet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_genjet_M_jj);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_recojet_M_jj);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_M_WBoson);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_M_qq);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_Selected.push_back(h_WW_2000_CT_wO_Selected_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_2000_CT_wO_Selected.size();i++){
    hist_vec_WW_2000_CT_wO_Selected[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_2000_CT_wO_Selected,hist_vec_WW_2000_CT_wO_Selected);

  TH1F* h_WW_2000_CT_wO_allPFOs_genjet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_allPFOs_genjet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_allPFOs_recojet_E1_plus_E2 =  new TH1F("h_WW_2000_CT_wO_allPFOs_recojet_E1_plus_E2","", n_binsJER, 50,2250);
  TH1F* h_WW_2000_CT_wO_allPFOs_genjet_M_jj =  new TH1F("h_WW_2000_CT_wO_allPFOs_genjet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_allPFOs_recojet_M_jj =  new TH1F("h_WW_2000_CT_wO_allPFOs_recojet_M_jj","", n_binsJER, 50.,120.);
  TH1F* h_WW_2000_CT_wO_allPFOs_M_WBoson =  new TH1F("h_WW_2000_CT_wO_allPFOs_parton_M_WBoson","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_allPFOs_M_qq =  new TH1F("h_WW_2000_CT_wO_allPFOs_parton_M_qq","", 10*n_binsJER, 20.,2000.);
  TH1F* h_WW_2000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_allPFOs_genjet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_allPFOs_genjet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60  =  new TH1F("h_WW_2000_CT_wO_allPFOs_recojet_dphijjMET_Mjj60","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75  =  new TH1F("h_WW_2000_CT_wO_allPFOs_recojet_dphijjMET_Mjj75","", n_binsJER, 50.,180.);
  TH1F* h_WW_2000_CT_wO_allPFOs_W_qq_theta =  new TH1F("h_WW_2000_CT_wO_allPFOs_W_qq_theta","", n_binsJER, 0.,180.);

  std::vector<TH1F*> hist_vec_WW_2000_CT_wO_allPFOs;
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_genjet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_recojet_E1_plus_E2);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_genjet_M_jj);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_recojet_M_jj);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_M_WBoson);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_M_qq);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_genjet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj60);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_recojet_dPhijjMET_Mjj75);
  hist_vec_WW_2000_CT_wO_allPFOs.push_back(h_WW_2000_CT_wO_allPFOs_W_qq_theta);

  for(unsigned int i=0;i<hist_vec_WW_2000_CT_wO_allPFOs.size();i++){
    hist_vec_WW_2000_CT_wO_allPFOs[i]->Sumw2();
  }
  fill_mass_histograms(file_CLIC_WW_2000_CT_wO_allPFOs,hist_vec_WW_2000_CT_wO_allPFOs);
  */
  std::cout<<"WW_2000 done"<<std::endl;


  file_histogram->Write();
  file_histogram->Close();

}
