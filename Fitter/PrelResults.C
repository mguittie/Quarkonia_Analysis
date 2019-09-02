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

const int nBins           = 9;

Double_t NJPsi[]     = {131927.0, 139000.0, 135933.0, 99213.0, 121177.0, 74624.0, 42770.0, 20622.0, 9561.0, 1123.0};
Double_t StatJPsi[]  = {666.0, 939.0, 899.0, 670.0, 937.0, 729.0, 658.0, 548.0, 284.0, 260.0};
Double_t SysJPsi[]   = {3382.0, 2916.0, 3418.0, 3284.0, 3203.0, 2603.0, 1046.0, 447.0, 347.0, 281.0};

Double_t NPsi2S[]    = {3212.0, 3040.0, 2931.0, 1696.0, 2118.0, 1692.0, 711.0, 354.0, 151.0, 18.0};
Double_t StatPsi2S[] = {159.0, 189.0, 190.0, 169.0, 204.0, 165.0, 133.0, 98.0, 71.0, 28.0};
Double_t SysPsi2S[]  = {276.0, 182.0, 335.0, 293.0, 285.0, 270.0, 100.0, 52.0, 37.0, 8.0};

Double_t Nch[]       = {5.63, 14.03, 21.43, 28.12, 35.45, 44.70, 54.22, 64.71, 78.16, 100.74};
Double_t StatNch[]   = {0.001, 0.002, 0.002, 0.003, 0.004, 0.006, 0.010, 0.016, 0.045, 0.206};
Double_t SysNch[]    = {0.089, 0.038, 0.062, 0.113, 0.215, 0.388, 0.615, 0.896, 1.343, 2.211};

Double_t meanNch         = 14.02;
Double_t statMeanNch     = 0.005;
Double_t sysMeanNch      = 0.038;

// Valeurs Bin 1 (9-14)
// Double_t meanIntegratedRatio     = NPsi2S[1]/NJPsi[1];
// Double_t statMeanIntegratedRatio = NPsi2S[1]/NJPsi[1]*sqrt(TMath::Power(StatPsi2S[1]/NPsi2S[1],2) + TMath::Power(StatJPsi[1]/NJPsi[1],2));
// Double_t sysMeanIntegratedRatio  = NPsi2S[1]/NJPsi[1]*sqrt(TMath::Power(SysPsi2S[1]/NPsi2S[1],2) + TMath::Power(SysJPsi[1]/NJPsi[1],2));

Double_t meanIntegratedRatio     = NPsi2S[1]/NJPsi[1];
Double_t statMeanIntegratedRatio = NPsi2S[1]/NJPsi[1]*sqrt(TMath::Power(StatPsi2S[1]/NPsi2S[1],2) + TMath::Power(StatJPsi[1]/NJPsi[1],2));
Double_t ErrorMeanIntegratedRatio  = sqrt(TMath::Power(TMath::Power(StatPsi2S[1]/NPsi2S[1],2) + TMath::Power(StatJPsi[1]/NJPsi[1],2) + SysPsi2S[1]/NPsi2S[1],2) + TMath::Power(SysJPsi[1]/NJPsi[1],2));

// Valeurs fit Integrated Mult (Voir pourquoi ça ne fit pas)
// Double_t meanIntegratedRatio     = 0;
// Double_t statMeanIntegratedRatio = 0;
// Double_t sysMeanIntegratedRatio  = 0;

void PrelResults(){

  Double_t MeanRatio = 0, statMeanRatio = 0, sysMeanRatio = 0;
  Double_t Ratio[nBins], StatRatio[nBins], SysRatio[nBins], NchRatio[nBins], statNchRatio[nBins], sysNchRatio[nBins];

  TGraphErrors *gResults = new TGraphErrors(nBins);
  TGraphErrors *gSys     = new TGraphErrors(nBins);
  TGraphErrors *gGlobalSys     = new TGraphErrors(1);

  for (int i = 0; i<nBins; i++){
    Ratio[i]       = NPsi2S[i]/NJPsi[i];
    StatRatio[i]   = NPsi2S[i]/NJPsi[i]*sqrt(TMath::Power(StatPsi2S[i]/NPsi2S[i],2) + TMath::Power(StatJPsi[i]/NJPsi[i],2));
    SysRatio[i]    = NPsi2S[i]/NJPsi[i]*sqrt(TMath::Power(SysPsi2S[i]/NPsi2S[i],2) + TMath::Power(SysJPsi[i]/NJPsi[i],2));

    NchRatio[i]     = Nch[i]/meanNch;
    statNchRatio[i] = NchRatio[i]*sqrt(TMath::Power(StatNch[i]/Nch[i],2) + TMath::Power(statMeanNch/meanNch,2));
    sysNchRatio[i]  = NchRatio[i]*sqrt(TMath::Power(SysNch[i]/Nch[i],2) + TMath::Power(sysMeanNch/meanNch,2));

    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "Bin " << i << endl;
    cout << "Ratio = " << Ratio[i] << " +/- " << StatRatio[i] << " (stat) +/- " << SysRatio[i] << " (sys)" << endl;
    cout << "Nch = " << Nch[i] << " +/- " << StatNch[i] << " (stat) +/- " << SysNch[i] << " (sys)" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;
  }

  // cout << "<Psi2S/JPsi> = " << MeanRatio << " +/- " << statMeanRatio << " (stat) +/- " << sysMeanRatio << " (sys)" << endl;

  // for (int i = 0; i<nBins; i++){
  //   gResults->SetPoint(i, NchRatio[i], Ratio[i]/meanIntegratedRatio);
  //   gResults->SetPointError(i, statNchRatio[i], Ratio[i]/meanIntegratedRatio*sqrt(TMath::Power(StatRatio[i]/Ratio[i],2) + TMath::Power(statMeanIntegratedRatio/meanIntegratedRatio,2)));
  //   gSys->SetPoint(i, NchRatio[i], Ratio[i]/meanIntegratedRatio);
  //   gSys->SetPointError(i, sysNchRatio[i], Ratio[i]/meanIntegratedRatio*sqrt(TMath::Power(SysRatio[i]/Ratio[i],2) + TMath::Power(sysMeanIntegratedRatio/meanIntegratedRatio,2)));
  // }

  for (int i = 0; i<nBins; i++){
    gResults->SetPoint(i, NchRatio[i], Ratio[i]/meanIntegratedRatio);
    gResults->SetPointError(i, statNchRatio[i], StatRatio[i]/meanIntegratedRatio);
    gSys->SetPoint(i, NchRatio[i], Ratio[i]/meanIntegratedRatio);
    gSys->SetPointError(i, sysNchRatio[i], SysRatio[i]/meanIntegratedRatio);
  }

  TCanvas *cResults = new TCanvas("cResults","cResults",200,10,700,500);

  gResults->GetYaxis()->SetTitle("#frac{#psi(2S)/J/#psi}{<#psi(2S)/J/#psi>}");
  gResults->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  gSys->GetYaxis()->SetTitle("");
  gSys->GetXaxis()->SetTitle("");
  gSys->SetFillColor(2);
  gSys->SetFillStyle(3001);

  TF1 *Mean = new TF1("Mean","1",0,9);

  gResults->Draw("AP");
  gSys->Draw("E2same");
  Mean->Draw("same");


  gGlobalSys->SetPoint(0, 6., 1);
  gGlobalSys->SetPointError(0, 0.3, ErrorMeanIntegratedRatio);
  gGlobalSys->SetFillColor(2);
  gGlobalSys->SetFillStyle(3001);

  gGlobalSys->Draw("E2same");



  // autre méthode de calcul

  // Double_t MeanJPsi = 0, statMeanJPsi = 0, sysMeanJPsi = 0;
  // Double_t MeanPsi2S = 0, statMeanPsi2S = 0, sysMeanPsi2S = 0;
  // Double_t RatioJPsi[nBins], StatRatioJPsi[nBins], SysRatioJPsi[nBins], RatioPsi2S[nBins], StatRatioPsi2S[nBins], SysRatioPsi2S[nBins];
  //
  // TGraphErrors *gResults2 = new TGraphErrors(nBins);
  // TGraphErrors *gSys2     = new TGraphErrors(nBins);
  //
  // for (int i = 0; i<nBins; i++){
  //   MeanJPsi      += NJPsi[i];
  //   MeanPsi2S     += NPsi2S[i];
  //   statMeanJPsi  += TMath::Power(StatJPsi[i],2);
  //   statMeanPsi2S += TMath::Power(StatPsi2S[i],2);
  //   sysMeanJPsi   += TMath::Power(SysJPsi[i],2);
  //   sysMeanPsi2S  += TMath::Power(SysPsi2S[i],2);
  // }
  // MeanJPsi       = MeanJPsi/nBins;
  // MeanPsi2S      = MeanPsi2S/nBins;
  // statMeanJPsi   = sqrt(statMeanJPsi)/nBins;
  // sysMeanJPsi    = sqrt(sysMeanJPsi)/nBins;
  // statMeanPsi2S  = sqrt(statMeanPsi2S)/nBins;
  // sysMeanPsi2S   = sqrt(sysMeanPsi2S)/nBins;
  //
  // for (int i; i<nBins; i++){
  //   RatioJPsi[i]      = NJPsi[i]/MeanJPsi;
  //   RatioPsi2S[i]     = NPsi2S[i]/MeanPsi2S;
  //   StatRatioJPsi[i]  = RatioJPsi[i]*sqrt(TMath::Power(StatJPsi[i]/NJPsi[i],2) + TMath::Power(statMeanJPsi/MeanJPsi,2));
  //   StatRatioPsi2S[i] = RatioPsi2S[i]*sqrt(TMath::Power(StatPsi2S[i]/NPsi2S[i],2) + TMath::Power(statMeanPsi2S/MeanPsi2S,2));
  //   SysRatioJPsi[i]   = RatioJPsi[i]*sqrt(TMath::Power(SysJPsi[i]/NJPsi[i],2) + TMath::Power(sysMeanJPsi/MeanJPsi,2));
  //   SysRatioPsi2S[i]  = RatioPsi2S[i]*sqrt(TMath::Power(SysPsi2S[i]/NPsi2S[i],2) + TMath::Power(sysMeanPsi2S/MeanPsi2S,2));
  // }
  //
  // for (int i = 0; i<nBins; i++){
  //   gResults2->SetPoint(i, NchRatio[i], RatioPsi2S[i]/RatioJPsi[i]);
  //   gResults2->SetPointError(i, statNchRatio[i], RatioPsi2S[i]/RatioJPsi[i]*sqrt(TMath::Power(StatRatioPsi2S[i]/RatioPsi2S[i],2) + TMath::Power(StatRatioJPsi[i]/RatioJPsi[i],2)));
  //   gSys2->SetPoint(i, NchRatio[i], RatioPsi2S[i]/RatioJPsi[i]);
  //   gSys2->SetPointError(i, sysNchRatio[i], RatioPsi2S[i]/RatioJPsi[i]*sqrt(TMath::Power(SysRatioPsi2S[i]/RatioPsi2S[i],2) + TMath::Power(SysRatioJPsi[i]/RatioJPsi[i],2)));
  // }
  //
  // TCanvas *cResults2 = new TCanvas("cResults2","cResults2",200,10,700,500);
  //
  // gResults2->GetYaxis()->SetTitle("#frac{#psi(2S)/J/#psi}{<#psi(2S)/J/#psi>}");
  // gResults2->GetXaxis()->SetTitle("#frac{dN_{ch}/d#eta}{<dN_{ch}/d#eta>}");
  // gSys2->GetYaxis()->SetTitle("");
  // gSys2->GetXaxis()->SetTitle("");
  // gSys2->SetFillColor(2);
  // gSys2->SetFillStyle(3001);
  //
  // TF1 *Mean2 = new TF1("Mean2","1",0,9);
  //
  // gResults2->Draw("AP");
  // gSys2->Draw("E2same");
  // Mean2->Draw("same");
  return;
}
