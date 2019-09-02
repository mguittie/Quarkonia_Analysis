#include "TMath.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

#include <iostream>
#include <string>

using namespace std;

const int nBins          = 9;

Double_t MultClassesMaxRef[20] = {1.0, 8.0, 9.0, 14.0, 15.0, 20.0, 21.0, 25.0, 26.0, 33.0, 34.0, 41.0, 42.0, 50.0, 51.0, 60.0, 61.0, 80.0, 81.0, 115.0};

Double_t Nch[]       = {5.69, 13.94, 21.42, 28.13, 35.47, 44.62, 53.99, 64.25, 76.92, 98.51};

Double_t Nch2016[]       = {5.63, 14.03, 21.43, 28.12, 35.45, 44.70, 54.22, 64.71, 78.16, 100.74};
Double_t StatNch2016[]   = {0.001, 0.002, 0.002, 0.003, 0.004, 0.006, 0.010, 0.016, 0.045, 0.206};
Double_t SysNch2016[]    = {0.089, 0.038, 0.062, 0.113, 0.215, 0.388, 0.615, 0.896, 1.343, 2.211};

Double_t meanNch         = 14.04;
Double_t statMeanNch     = 0.005;
Double_t sysMeanNch      = 0.038;


void NchComparisonWith2016(){

  TGraphErrors *gNch        = new TGraphErrors(nBins);
  TGraphErrors *gSysNch     = new TGraphErrors(nBins);


  for (int i = 1; i<=nBins; i++){

    gNch->SetPoint(i-1, (MultClassesMaxRef[2*i-2] + MultClassesMaxRef[2*i-1])/2, Nch2016[i]/Nch[i]);
    gNch->SetPointError(i-1, 1, Nch2016[i]/Nch[i]*sqrt(TMath::Power(StatNch2016[i]/Nch2016[i],2) + TMath::Power(statMeanNch/meanNch,2)));
    gSysNch->SetPoint(i-1, (MultClassesMaxRef[2*i-2] + MultClassesMaxRef[2*i-1])/2, Nch2016[i]/Nch[i]);
    gSysNch->SetPointError(i-1, 1, Nch2016[i]/Nch[i]*sqrt(TMath::Power(SysNch2016[i]/Nch2016[i],2) + TMath::Power(sysMeanNch/meanNch,2)));
  }

  TCanvas *cComparison = new TCanvas("cComparison","cComparison",1500,1000);

  gNch->GetYaxis()->SetTitle("(#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>})_{J/#psi 2016}/(#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>})_{This Analysis}");
  gNch->GetXaxis()->SetTitle("N_{trk}^{Corr}");
  gSysNch->GetYaxis()->SetTitle("");
  gSysNch->GetXaxis()->SetTitle("");
  gSysNch->SetFillColor(4);
  gSysNch->SetFillStyle(3001);

  TF1 *Mean = new TF1("Mean","1",0,100);

  cComparison->cd();
  gNch->Draw("AP");
  gSysNch->Draw("E2same");
  Mean->Draw("same");

  return;
}
