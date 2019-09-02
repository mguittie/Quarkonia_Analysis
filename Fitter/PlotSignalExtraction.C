// #include "FitOneHisto.C"
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
TString tailNames[]   = {"MCTail", "DataTail", "FreeTail"};

// Invariant mass ranges
Double_t Ranges[nbRanges][2] = {{2.0,5.0},{2.2,4.5}};

void PlotSignalExtraction(TString Opt = "multbins", TString ref = "max")
{
  TCanvas* canNJPsi = new TCanvas("canNJPsi","NJPsi",4000,3000);
  SetCanvasStyle(canNJPsi);

  TCanvas* canNPsi2S = new TCanvas("canNPsi2S","NPsi2S",4000,3000);
  SetCanvasStyle(canNPsi2S);

  TCanvas* canMeanJPsi = new TCanvas("canMeanJPsi","Mean JPsi",4000,3000);
  SetCanvasStyle(canMeanJPsi);

  TCanvas* canWidthJPsi = new TCanvas("canWidthJPsi","Width JPsi",4000,3000);
  SetCanvasStyle(canWidthJPsi);

  TCanvas* canChiSquare = new TCanvas("canChiSquare","Chi Square",4000,3000);
  SetCanvasStyle(canChiSquare);

  TCanvas* canPsi2SOverJPsi = new TCanvas("canPsi2SOverJPsi","Psi2S/JPsi",4000,3000);
  SetCanvasStyle(canPsi2SOverJPsi);

  if (Opt == "multbins"){
    canNJPsi->Divide(4,3);
    canNPsi2S->Divide(4,3);
    canMeanJPsi->Divide(4,3);
    canWidthJPsi->Divide(4,3);
    canChiSquare->Divide(4,3);
    canPsi2SOverJPsi->Divide(4,3);
  }


  std::vector<TString> vectorTests;

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
          vectorTests.push_back(testLabel);
        }
      }
    }
  }

  TString rangeName;

  if (Opt == "multbins"){

    for (int multclass = 1; multclass<=10; multclass ++){

      rangeName.Form("%g_%g", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);

      canNJPsi->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 0);

      canNPsi2S->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 1);

      canMeanJPsi->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 2);

      canWidthJPsi->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 3);

      canChiSquare->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 4);

      canPsi2SOverJPsi->cd(multclass);
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 5);
    }
  }

  if (Opt == "integrated"){

      if (ref == "max"){
        rangeName.Form("%g_%g", MultClassesMaxRef[0], MultClassesMaxRef[19]);
      }
      else if (ref == "min"){
        rangeName.Form("%g_%g", MultClassesMinRef[0], MultClassesMinRef[19]);
      }
      else{
        cout << "----- Wrong reference value -----" << endl;
      }

      canNJPsi->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 0);

      canNPsi2S->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 1);

      canMeanJPsi->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 2);

      canWidthJPsi->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 3);

      canChiSquare->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 4);

      canPsi2SOverJPsi->cd();
      gPad->SetBottomMargin(0.2);
      gPad->SetRightMargin(0.1);
      CreateHistosSigExtAndRMS(rangeName, vectorTests, 5);
  }

  canNJPsi->SaveAs(Form("Plots/NJPsiVsTest_%s.pdf", Opt.Data()));
  canMeanJPsi->SaveAs(Form("Plots/MeanJPsiVsTest_%s.pdf", Opt.Data()));
  canWidthJPsi->SaveAs(Form("Plots/WidthJPsiVsTest_%s.pdf", Opt.Data()));
  canNPsi2S->SaveAs(Form("Plots/NPsi2SVsTest_%s.pdf", Opt.Data()));
  canChiSquare->SaveAs(Form("Plots/Chi2PerNdfVsTest_%s.pdf", Opt.Data()));
  canPsi2SOverJPsi->SaveAs(Form("Plots/Psi2SOverJPsi_%s.pdf", Opt.Data()));
}
