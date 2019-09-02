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

Double_t NJPsi[]         = {131927.0, 139000.0, 135933.0, 99213.0, 121177.0, 74624.0, 42770.0, 20622.0, 9561.0, 1123.0};
Double_t StatJPsi[]      = {666.0, 939.0, 899.0, 670.0, 937.0, 729.0, 658.0, 548.0, 284.0, 260.0};
Double_t SysJPsi[]       = {3382.0, 2916.0, 3418.0, 3284.0, 3203.0, 2603.0, 1046.0, 447.0, 347.0, 281.0};

Double_t NJPsiInt          = 768433.0;
Double_t statNJPsiInt      = 1753.0;
Double_t sysNJPsiInt       = 20262.0;

Double_t NJPsiShift[]         = {113779.0, 136141.0, 139130.0, 101724.0, 126560.0, 79957.0, 46684.0, 22278.0, 10583.0, 1240.0};
Double_t StatJPsiShift[]      = {554.0, 939.0, 899.0, 670.0, 937.0, 729.0, 658.0, 548.0, 284.0, 260.0};
Double_t SysJPsiShift[]       = {3382.0, 2916.0, 3418.0, 3284.0, 3203.0, 2603.0, 1046.0, 447.0, 347.0, 281.0};

Double_t NJPsiIntShift          = 771338.0;
Double_t statNJPsiIntShift      = 1753.0;
Double_t sysNJPsiIntShift       = 20262.0;

// Double_t NJPsi2016[]     = {109250.0, 134783.0, 138667.0, 102139.0, 125379.0, 78618.0, 46047.0, 21833.0, 10485.0, 953.0};
// Double_t StatJPsi2016[]  = {435.0, 516.0, 536.0, 484.0, 553.0, 479.0, 373.0, 247.0, 189.0, 59.0};
// Double_t SysJPsi2016[]   = {3065.0, 3693.0, 3908.0, 3034.0, 3601.0, 2287.0, 1348.0, 641.0, 317.0, 24.0};
//
// Double_t NJPsiInt2016          = 771136.0;
// Double_t statNJPsiInt2016      = 1536.0;
// Double_t sysNJPsiInt2016       = 21872.0;

// Values sent by Dhananjaya Thakur (10/07)
Double_t NJPsi2016[]     = {136225.0, 142656.0, 142289.0, 103458.0, 123591.0, 76185.0, 43961.0, 20741.0, 9555.0, 889.0}; //Only mean value modified
Double_t StatJPsi2016[]  = {435.0, 516.0, 536.0, 484.0, 553.0, 479.0, 373.0, 247.0, 189.0, 59.0};
Double_t SysJPsi2016[]   = {3065.0, 3693.0, 3908.0, 3034.0, 3601.0, 2287.0, 1348.0, 641.0, 317.0, 24.0};

Double_t NJPsiInt2016          = 802824.0;  //Only mean value modified
Double_t statNJPsiInt2016      = 1536.0;
Double_t sysNJPsiInt2016       = 21872.0;

Double_t Nch[]       = {5.63, 14.03, 21.43, 28.12, 35.45, 44.70, 54.22, 64.71, 78.16, 100.74};
Double_t StatNch[]   = {0.001, 0.002, 0.002, 0.003, 0.004, 0.006, 0.010, 0.016, 0.045, 0.206};
Double_t SysNch[]    = {0.089, 0.038, 0.062, 0.113, 0.215, 0.388, 0.615, 0.896, 1.343, 2.211};

Double_t meanNch         = 14.02;
Double_t statMeanNch     = 0.005;
Double_t sysMeanNch      = 0.038;


void ComparisonJPsiExtraction2016(){

  Double_t SumBins, SumBins2016;

  TGraphErrors *gResults = new TGraphErrors(nBins);
  TGraphErrors *gSys     = new TGraphErrors(nBins);

  TGraphErrors *gResults2016 = new TGraphErrors(nBins);
  TGraphErrors *gSys2016     = new TGraphErrors(nBins);

  TGraphErrors *gRatio    = new TGraphErrors(nBins);
  TGraphErrors *gSysRatio = new TGraphErrors(nBins);

  TGraphErrors *gRatio2016    = new TGraphErrors(nBins);
  TGraphErrors *gSysRatio2016 = new TGraphErrors(nBins);

  TGraphErrors *gRatio2016Shift    = new TGraphErrors(nBins);
  TGraphErrors *gSysRatio2016Shift = new TGraphErrors(nBins);

  for (int i = 0; i<nBins; i++){
    SumBins += NJPsi[i];
    SumBins2016 += NJPsi2016[i];

    gResults->SetPoint(i, Nch[i]/meanNch, NJPsi[i]/NJPsiInt);
    gResults->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2)), NJPsi[i]/NJPsiInt*sqrt(TMath::Power(StatJPsi[i]/NJPsi[i],2) + TMath::Power(statNJPsiInt/NJPsiInt,2)));
    gSys->SetPoint(i, Nch[i]/meanNch, NJPsi[i]/NJPsiInt);
    gSys->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2)), NJPsi[i]/NJPsiInt*sqrt(TMath::Power(SysJPsi[i]/NJPsi[i],2) + TMath::Power(sysNJPsiInt/NJPsiInt,2)));

    gResults2016->SetPoint(i, Nch[i]/meanNch, NJPsi2016[i]/NJPsiInt2016);
    gResults2016->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2)), NJPsi2016[i]/NJPsiInt2016*sqrt(TMath::Power(StatJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(statNJPsiInt2016/NJPsiInt2016,2)));
    gSys2016->SetPoint(i, Nch[i]/meanNch, NJPsi2016[i]/NJPsiInt2016);
    gSys2016->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2)), NJPsi2016[i]/NJPsiInt2016*sqrt(TMath::Power(SysJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(sysNJPsiInt2016/NJPsiInt2016,2)));

    gRatio->SetPoint(i, Nch[i]/meanNch, (NJPsi[i]/NJPsiInt)/(NJPsi2016[i]/NJPsiInt2016));
    gRatio->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2)), (NJPsi[i]/NJPsiInt)/(NJPsi2016[i]/NJPsiInt2016)*sqrt(TMath::Power(StatJPsi[i]/NJPsi[i],2) + TMath::Power(statNJPsiInt/NJPsiInt,2)));
    gSysRatio->SetPoint(i, Nch[i]/meanNch, (NJPsi[i]/NJPsiInt)/(NJPsi2016[i]/NJPsiInt2016));
    gSysRatio->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2)), (NJPsi[i]/NJPsiInt)/(NJPsi2016[i]/NJPsiInt2016)*sqrt(TMath::Power(SysJPsi[i]/NJPsi[i],2) + TMath::Power(sysNJPsiInt/NJPsiInt,2)));

    gRatio2016->SetPoint(i, Nch[i]/meanNch, (NJPsi2016[i]/NJPsiInt2016)/(NJPsi[i]/NJPsiInt));
    gRatio2016->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2)), (NJPsi2016[i]/NJPsiInt2016)/(NJPsi[i]/NJPsiInt)*sqrt(TMath::Power(StatJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(statNJPsiInt2016/NJPsiInt2016,2)));
    gSysRatio2016->SetPoint(i, Nch[i]/meanNch, (NJPsi2016[i]/NJPsiInt2016)/(NJPsi[i]/NJPsiInt));
    gSysRatio2016->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2)), (NJPsi2016[i]/NJPsiInt2016)/(NJPsi[i]/NJPsiInt)*sqrt(TMath::Power(SysJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(sysNJPsiInt2016/NJPsiInt2016,2)));

    gRatio2016Shift->SetPoint(i, Nch[i]/meanNch, (NJPsi2016[i]/NJPsiInt2016)/(NJPsiShift[i]/NJPsiIntShift));
    gRatio2016Shift->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2)), (NJPsi2016[i]/NJPsiInt2016)/(NJPsiShift[i]/NJPsiIntShift)*sqrt(TMath::Power(StatJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(statNJPsiInt2016/NJPsiInt2016,2)));
    gSysRatio2016Shift->SetPoint(i, Nch[i]/meanNch, (NJPsi2016[i]/NJPsiInt2016)/(NJPsiShift[i]/NJPsiIntShift));
    gSysRatio2016Shift->SetPointError(i, Nch[i]/meanNch*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2)), (NJPsi2016[i]/NJPsiInt2016)/(NJPsiShift[i]/NJPsiIntShift)*sqrt(TMath::Power(SysJPsi2016[i]/NJPsi2016[i],2) + TMath::Power(sysNJPsiInt2016/NJPsiInt2016,2)));
  }

  cout << "Sum over Mult bins NJPsi = " << SumBins << endl;
  cout << "Integrated NJPsi = " << NJPsiInt << endl;
  cout << "Sum over Mult bins NJPsi (2016) = " << SumBins2016 << endl;
  cout << "Integrated NJPsi (2016) = " << NJPsiInt2016 << endl;

  TCanvas *cComparison = new TCanvas("cComparison","cComparison",1500,1000);
  TCanvas *cRatio = new TCanvas("cRatio","cRatio",1500,1000);
  TCanvas *cRatio2016 = new TCanvas("cRatio2016","cRatio2016",1500,1000);
  TCanvas *cRatio2016Shift = new TCanvas("cRatio2016Shift","cRatio2016Shift",1500,1000);

  gResults->GetYaxis()->SetTitle("#frac{N_{J/#psi}}{<N_{J/#psi}>}");
  gResults->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  gSys->GetYaxis()->SetTitle("");
  gSys->GetXaxis()->SetTitle("");
  gSys->SetFillColor(4);
  gSys->SetFillStyle(3001);

  gResults2016->GetYaxis()->SetTitle("");
  gResults2016->GetXaxis()->SetTitle("");
  gSys2016->GetYaxis()->SetTitle("");
  gSys2016->GetXaxis()->SetTitle("");
  gSys2016->SetFillColor(2);
  gSys2016->SetFillStyle(3001);

  gRatio->GetYaxis()->SetTitle("#frac{N_{J/#psi}}{<N_{J/#psi}>}/2016");
  gRatio->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  gSysRatio->GetYaxis()->SetTitle("");
  gSysRatio->GetXaxis()->SetTitle("");
  gSysRatio->SetFillColor(4);
  gSysRatio->SetFillStyle(3001);

  gRatio2016->GetYaxis()->SetTitle("(#frac{N_{J/#psi}}{<N_{J/#psi}>})_{J/#psi 2016}/(#frac{N_{J/#psi}}{<N_{J/#psi}>})_{This Analysis}");
  gRatio2016->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  gSysRatio2016->GetYaxis()->SetTitle("");
  gSysRatio2016->GetXaxis()->SetTitle("");
  gSysRatio2016->SetFillColor(4);
  gSysRatio2016->SetFillStyle(3001);

  gRatio2016Shift->GetYaxis()->SetTitle("(#frac{N_{J/#psi}}{<N_{J/#psi}>})_{J/#psi 2016}/(#frac{N_{J/#psi}}{<N_{J/#psi}>})_{This Analysis (shifted)}");
  gRatio2016Shift->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  gSysRatio2016Shift->GetYaxis()->SetTitle("");
  gSysRatio2016Shift->GetXaxis()->SetTitle("");
  gSysRatio2016Shift->SetFillColor(4);
  gSysRatio2016Shift->SetFillStyle(3001);

  TF1 *Mean = new TF1("Mean","1",0,9);

  cComparison->cd();
  gResults->Draw("AP");
  gSys->Draw("E2same");
  gResults2016->Draw("Psame");
  gSys2016->Draw("E2same");

  cRatio->cd();
  gRatio->Draw("AP");
  gSysRatio->Draw("E2same");
  Mean->Draw("same");

  cRatio2016->cd();
  gRatio2016->Draw("AP");
  gSysRatio2016->Draw("E2same");
  Mean->Draw("same");

  cRatio2016Shift->cd();
  gRatio2016Shift->Draw("AP");
  gSysRatio2016Shift->Draw("E2same");
  Mean->Draw("same");

  return;
}
