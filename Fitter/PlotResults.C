#include "SignalExtraction.C"
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

// const int nBins = 10;

Double_t Nch[]       = {5.63, 14.03, 21.43, 28.12, 35.45, 44.70, 54.22, 64.71, 78.16, 100.74};
Double_t StatNch[]   = {0.001, 0.002, 0.002, 0.003, 0.004, 0.006, 0.010, 0.016, 0.045, 0.206};
Double_t SysNch[]    = {0.089, 0.038, 0.062, 0.113, 0.215, 0.388, 0.615, 0.896, 1.343, 2.211};

Double_t meanNch         = 14.04;
Double_t statMeanNch     = 0.005;
Double_t sysMeanNch      = 0.038;

void DrawDownAndUpLines(Double_t x0, Double_t x1,Double_t downYAt, Double_t upYAt, int color, int style, TLegend *leg, TString strForLineInLegend);

void PlotResults(TString ref = "max", const int nBins = 10){

  std::vector<TString> vectorTests;
  TString testName, rangeName;
  Int_t numberOfTests = 4*nbRanges*nbSig*nbBkg*2; // *4 = Bkg combinations, *2 = nbSig*2 = CB2 with 3 tails + NA60 with 1 tail

  Double_t DoubleRatioTest = 0, StatDoubleRatioTest = 0, StatDoubleRatioTestBin = 0, StatDoubleRatioTestInt = 0;
  Double_t DoubleRatio, StatDoubleRatio, SysDoubleRatio;
  Int_t iTest, numberOfWeights;

  Double_t NchRatio, statNchRatio, sysNchRatio;

  TCanvas* canDoubleRatio = new TCanvas("canDoubleRatio","Psi2S/JPsi/<Psi2S/JPsi>",4000,3000);
  SetCanvasStyle(canDoubleRatio);
  canDoubleRatio->Divide(4,3);

  // TCanvas* canRMS = new TCanvas("canRMS","RMS",4000,3000);
  // SetCanvasStyle(canRMS);
  // canRMS->Divide(4,3);

  TGraphErrors *gResults = new TGraphErrors(nBins);
  TGraphErrors *gSys     = new TGraphErrors(nBins);
  // TGraphErrors *gGlobalSys     = new TGraphErrors(1);

  for (int multclass = 1; multclass<=nBins; multclass ++){
    canDoubleRatio->cd(multclass);

    iTest = 0;
    numberOfWeights = 0;
    Double_t testWeight = 1;
    DoubleRatio = StatDoubleRatio = SysDoubleRatio = 0;
    DoubleRatioTest = StatDoubleRatioTest = StatDoubleRatioTestInt = 0;

    TH1F *histoResults = new TH1F(Form("histoResults_%g_%g",MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]),"",numberOfTests,0,numberOfTests);
    // TH1F *histoRMSResults = new TH1F(Form("histoResultsRMS_%g_%g",MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]),"",1000,0,2.0);


    if (ref == "min"){
      delete histoResults;
      TH1F *histoResults = new TH1F(Form("histoResults_%g_%g",MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]),"",numberOfTests,0,numberOfTests);
    }

    for(int sig=0; sig<nbSig; sig++){
    // sig = 0 : CB2  ;   sig = 1 : NA60

      for(int tail=0; tail<nbTails; tail++){
      // tail = 0 : MC  ;   tail = 1 : pp13TeV  ;  tail = 2 : Free Tail
        if ( sig == 1 && tail == 2) break; // Pas de Free tail pour la NA60
        if ( sig == 1 && tail == 1) break; // Pas de tail pp13TeV pour NA60

        for(int range=0; range<nbRanges; range++){
          for(int bkg=0; bkg<nbBkg; bkg++){
            if (ref == "max"){
              rangeName.Form("%g_%g", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);
            }
            else if (ref == "min"){
              rangeName.Form("%g_%g", MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]);
            }
            else{
              cout << "----- Wrong reference value -----" << endl;
            }

            testName.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data());

            std::vector<double> testFitResults, testFitResults_err;
            GetFitResults(rangeName, testName, testFitResults, testFitResults_err);

            // Loop on bkg functions for integrated Mult
            for(int rangeInt=0; rangeInt<nbRanges; rangeInt++){
              for(int bkgInt=0; bkgInt<nbBkg; bkgInt++){
                if (ref == "max"){
                  rangeName.Form("%g_%g", MultClassesMaxRef[0], MultClassesMaxRef[19]);
                }
                else if (ref == "min"){
                  rangeName.Form("%g_%g", MultClassesMinRef[0], MultClassesMinRef[19]);
                }
                else{
                  cout << "----- Wrong reference value -----" << endl;
                }

                testName.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkgInt]].Data(), Ranges[rangeInt][0], Ranges[rangeInt][1], tailNames[tail].Data());

                std::vector<double> testFitResultsInt, testFitResultsInt_err;
                GetFitResults(rangeName, testName, testFitResultsInt, testFitResultsInt_err);

                DoubleRatioTest        = (testFitResults[1]/testFitResults[0])/(testFitResultsInt[1]/testFitResultsInt[0]);
                StatDoubleRatioTestBin = testFitResults[1]/testFitResults[0]*sqrt(TMath::Power(testFitResults_err[1]/testFitResults[1],2) + TMath::Power(testFitResults_err[0]/testFitResults[0],2));
                StatDoubleRatioTestInt = testFitResultsInt[1]/testFitResultsInt[0]*sqrt(TMath::Power(testFitResultsInt_err[1]/testFitResultsInt[1],2) + TMath::Power(testFitResultsInt_err[0]/testFitResultsInt[0],2));
                StatDoubleRatioTest    = DoubleRatioTest*sqrt(TMath::Power(StatDoubleRatioTestBin/(testFitResults[1]/testFitResults[0]),2) + TMath::Power(StatDoubleRatioTestInt/(testFitResultsInt[1]/testFitResultsInt[0]),2));

                TString doubleTestLabel;
                doubleTestLabel.Form("%s_%s_Range_%1.1f_%1.1f_%s_Int_%s_%1.1f_%1.1f", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data(), arrayFunctionNames[fBackground[bkgInt]].Data(), Ranges[rangeInt][0], Ranges[rangeInt][1]);

                histoResults->SetBinContent(iTest+1,DoubleRatioTest);
                histoResults->SetBinError(iTest+1,StatDoubleRatioTest);
                histoResults->GetXaxis()->SetBinLabel(iTest+1,doubleTestLabel);

                DoubleRatio     += DoubleRatioTest;
                StatDoubleRatio += StatDoubleRatioTest;

                iTest++;
                numberOfWeights += testWeight;

                // histoRMSResults->Fill(DoubleRatioTest);
              }
            }
          }
        }
      }
    }
    DoubleRatio     /= numberOfWeights;
    StatDoubleRatio /= numberOfWeights;

    for (int i=1; i<=iTest; i++){
      SysDoubleRatio += testWeight*TMath::Power( (histoResults->GetBinContent(i,1) - DoubleRatio), 2);
    }

    SysDoubleRatio /= numberOfWeights;
    SysDoubleRatio = sqrt(SysDoubleRatio);

    cout << "Mult range " << multclass << " : Double Ratio = " << DoubleRatio << " +/- " << StatDoubleRatio << " (stat) +/- << " << SysDoubleRatio << " (sys)" << endl;

    histoResults->SetMinimum(DoubleRatio-10*SysDoubleRatio);
    histoResults->SetMaximum(DoubleRatio+10*SysDoubleRatio);
    histoResults->SetMarkerStyle(20);
    histoResults->SetMarkerColor(2);
    histoResults->SetLineColor(1);
    histoResults->GetYaxis()->SetTitle("#frac{#psi(2S)/J/#psi}{<#psi(2S)/J/#psi>}");
    histoResults->GetYaxis()->SetTitleSize(0.05);
    histoResults->GetYaxis()->SetTitleOffset(1);
    histoResults->GetYaxis()->SetLabelSize(0.02);
    histoResults->Draw("p");

    TLine *lineAverage = new TLine(0,DoubleRatio,numberOfTests,DoubleRatio);
    lineAverage->SetLineColor(4);
    lineAverage->SetLineWidth(1);
    lineAverage->Draw("");

    TLegend* legRMSLines = new TLegend(0.142, 0.22, 0.6, 0.32);
    legRMSLines->SetFillStyle(0);
    legRMSLines->SetLineColorAlpha(0,0);
    legRMSLines->SetTextColor(kBlack);
    legRMSLines->SetMargin(0.1);
    DrawDownAndUpLines(0, numberOfTests, DoubleRatio-SysDoubleRatio, DoubleRatio+SysDoubleRatio, 4, 4, legRMSLines,"#pm RMS");
    DrawDownAndUpLines(0, numberOfTests, DoubleRatio-2*SysDoubleRatio, DoubleRatio+2*SysDoubleRatio, 3, 5, legRMSLines,"#pm 2*RMS");
    legRMSLines->Draw();

    TLegend* legPartNumber = new TLegend(0.05, 0.65, 0.8, 0.9);
    legPartNumber->SetFillStyle(0);
    legPartNumber->SetLineColorAlpha(0,0);
    legPartNumber->SetTextColor(kBlack);
    legPartNumber->SetMargin(0.1);

    if (ref == "max"){
      rangeName.Form("%g_%g", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);
    }
    else if (ref == "min"){
      rangeName.Form("%g_%g", MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]);
    }
    else{
      cout << "----- Wrong reference value -----" << endl;
    }

    legPartNumber->AddEntry("NULL",Form("Mult_%s,  #frac{#psi(2S)/J/#psi}{<#psi(2S)/J/#psi>} = %f #pm %f (stat) #pm %f (sys) (%2.2f %%)",rangeName.Data(), DoubleRatio, StatDoubleRatio , SysDoubleRatio, 100*SysDoubleRatio/DoubleRatio),"");
    legPartNumber->Draw();

    gPad->SetBottomMargin(0.2);
    gPad->SetRightMargin(0.1);

    // canRMS->cd(multclass);
    // cout << "RMS = " << histoRMSResults->GetRMS() << endl;
    // histoRMSResults->Draw();

    NchRatio     = Nch[multclass-1]/meanNch;
    statNchRatio = NchRatio*sqrt(TMath::Power(StatNch[multclass-1]/Nch[multclass-1],2) + TMath::Power(statMeanNch/meanNch,2));
    sysNchRatio  = NchRatio*sqrt(TMath::Power(SysNch[multclass-1]/Nch[multclass-1],2) + TMath::Power(sysMeanNch/meanNch,2));

    gResults->SetPoint(multclass, NchRatio, DoubleRatio);
    gResults->SetPointError(multclass, statNchRatio, StatDoubleRatio);
    gSys->SetPoint(multclass, NchRatio, DoubleRatio);
    gSys->SetPointError(multclass, sysNchRatio, SysDoubleRatio);
  }
  canDoubleRatio->SaveAs(Form("Plots/DoubleRatioExtraction_%sRef.pdf", ref.Data()));


  // ---- Final plot ----

  TCanvas *cResults = new TCanvas("cResults","cResults",600,450);

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


  // gGlobalSys->SetPoint(0, 6., 1);
  // gGlobalSys->SetPointError(0, 0.3, ErrorMeanIntegratedRatio);
  // gGlobalSys->SetFillColor(2);
  // gGlobalSys->SetFillStyle(3001);
  //
  // gGlobalSys->Draw("E2same");
  // ---- Final Plot ----

  return;
}

void DrawDownAndUpLines(Double_t x0, Double_t x1,Double_t downYAt, Double_t upYAt, int color, int style, TLegend *leg, TString strForLineInLegend)
{
  TLine *lineUp = new TLine(x0,upYAt,x1,upYAt);
  lineUp->SetLineColor(color);
  lineUp->SetLineWidth(1);
  lineUp->SetLineStyle(style);
  lineUp->Draw("same");

  TLine *lineDown = new TLine(x0,downYAt,x1,downYAt);
  lineDown->SetLineColor(color);
  lineDown->SetLineWidth(1);
  lineDown->SetLineStyle(style);
  lineDown->Draw("same");

  leg->AddEntry(lineUp,strForLineInLegend,"l");
}
