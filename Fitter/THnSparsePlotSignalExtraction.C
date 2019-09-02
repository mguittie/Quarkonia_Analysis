#include "FitOneHisto.C"
#include "THnSparse.h"
#include "TMath.h"
#include "DrawAndRMS.C"

#include <iostream>
#include <string>

const int nbRanges = 2, nbTails = 3, nbSig = 2, nbBkg = 2;

Double_t MultClassesMaxRef[20] = {1.0, 8.0, 9.0, 14.0, 15.0, 20.0, 21.0, 25.0, 26.0, 33.0, 34.0, 41.0, 42.0, 50.0, 51.0, 60.0, 61.0, 80.0, 81.0, 115.0};
Double_t MultClassesMinRef[20] = {1.0, 5.0, 6.0, 10.0, 11.0, 14.0, 15.0, 17.0, 18.0, 23.0, 24.0, 28.0, 29.0, 34.0, 35.0, 42.0, 43.0, 54.0, 55.0, 83.0};

// Signal and Background functions
int fSignal[]         = {kCB21S, kNA601S};
int fBackground[]     = {kVWG, kPol1verPol2};
TString tailNames[] = {"MCTail", "DataTail", "FreeTail"};
// Invariant mass ranges
Double_t Ranges[nbRanges][2] = {{2.0,5.0},{2.2,4.5}};

// Rappel THnSparse structure
// 0: MultClass
// 1: Test Number
// 2: NJpsi
// 3: Error Jpsi
// 4: NPsi2S
// 5: Error Psi2S
// 6: MJPsi
// 7: MJPsi Error
// 8: Sigma JPsi
// 9: Sigma JPsi Error
// 10: Chi2PerNdf

void PlotSignalExtraction(TString fileName = "SignalExtractionResults.root", Double_t fitRangeMin = 2.4, Double_t fitRangeMax = 4.7, int fSig = kCB21S, int fBkgd = kVWG, const int NTests = 16)
{

  TFile *file = new TFile(fileName.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileName.Data());
    return;
  }
  THnSparse *fhSigExt = static_cast<THnSparse*>(file->FindObjectAny("fhSigExt"));
  if (!fhSigExt) return;

  TCanvas* canNJPsi = new TCanvas("canNJPsi","",2000,1500);
  SetCanvasStyle(canNJPsi);
  canNJPsi->Divide(4,3);


  for (int multclass = 1; multclass<=9; multclass ++){
    Int_t test = 1;

    canNJPsi->cd(multclass);
    fhSigExt->GetAxis(0)->SetRangeUser(multclass, multclass+1);

    TH2D *hNJPsi = fhSigExt->Projection(2,1,"a");
    // hNJPsi->SetDirectory(0);
    hNJPsi->SetName(Form("hNJPsi_Mult_%d", multclass));
    hNJPsi->SetMarkerStyle(2);
    hNJPsi->SetMarkerSize(2);
    hNJPsi->SetMarkerColor(2);
    hNJPsi->SetTitle(Form("Number of J/#psi in Mult Class %.1f - %.1f ", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]));
    hNJPsi->GetYaxis()->SetTitle("N_{J/#psi}");

    for(int range=0; range<nbRanges; range++){
      for(int sig=0; sig<nbSig; sig++){
        // sig = 0 : CB2  ;   sig = 1 : NA60
        for(int tail=0; tail<nbTails; tail++){
          // tail = 0 : MC  ;   tail = 1 : pp13TeV  ;  tail = 2 : Free Tail
          if ( sig == 1 && tail == 2) break; // Pas de Free tail pour la NA60
          if ( tail == 1 && sig == 1) break; // Car pas de tail pp13TeV pour NA60

          for(int bkg=0; bkg<nbBkg; bkg++){
            TString testLabel;
            testLabel.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data());
            hNJPsi->GetXaxis()->SetBinLabel(test, testLabel);

            // cout << Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data()) << endl;
            test++;
          }
        }
      }
    }
    hNJPsi->Draw("same");
  }
  canNJPsi->SaveAs("Plots/NJPsiVsTest.pdf");


}
