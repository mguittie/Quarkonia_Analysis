#include "FitOneHisto.C"
#include "THnSparse.h"
#include "SetAndGetFitResults.C"

#include <iostream>
#include <string>

using namespace std;

const int nbRanges = 2, nbTails = 3, nbSig = 2, nbBkg = 2;

Double_t MultClassesMaxRef[20] = {1.0, 8.0, 9.0, 14.0, 15.0, 20.0, 21.0, 25.0, 26.0, 33.0, 34.0, 41.0, 42.0, 50.0, 51.0, 60.0, 61.0, 80.0, 81.0, 115.0};
Double_t MultClassesMinRef[20] = {1.0, 5.0, 6.0, 10.0, 11.0, 14.0, 15.0, 17.0, 18.0, 23.0, 24.0, 28.0, 29.0, 34.0, 35.0, 42.0, 43.0, 54.0, 55.0, 83.0};

// Tail parameters (pp@13TeV & MC)
Double_t NA60TailMC[8] = {-1.01, 0.002886, 0.4387, 0.2225, 2.164, 0.1924, 1.284, 0.06384};

Double_t CB2Tailpp[8]  = {0.9743, 7.36, 1.84, 16.06, 0, 0, 0, 0};
Double_t CB2TailMC[8]  = {0.9803, 3.123, 2.198, 2.682, 0, 0, 0, 0};

// Signal and Background functions
int fSignal[]         = {kCB21S, kNA601S};
int fBackground[]     = {kVWG, kPol1verPol2};
TString tailNames[] = {"MCTail", "DataTail", "FreeTail"};

// Invariant mass ranges
Double_t Ranges[nbRanges][2] = {{2.0,5.0},{2.2,4.5}};

// parameter for Free Tails;
Bool_t fixTail = kTRUE;

std::vector<double> GetTailsParam(int iTail, int function);

void SignalExtraction(TString fileName = "AnalysisResults.root", TString Opt = "full", TString ref = "max", Double_t fitRangeMin = 2.4, Double_t fitRangeMax = 4.7, int fSig = kCB21S, int fBkgd = kVWG, Int_t tailType = 0, Double_t multRangeMin = 0.0, Double_t multRangeMax = 115.0)
{

  TFile *file = new TFile(fileName.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileName.Data());
    return;
  }
  THnSparse *fhDimuon = static_cast<THnSparse*>(file->FindObjectAny("CMUL_fhDimuon"));
  if (!fhDimuon) return;

  std::vector<double> tailsParams;
  std::vector<double> fitResults;

  // prepare binning for THnSparse for Signal Extraction
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
  const Int_t nDim = 11;
  Int_t nBins[nDim] = {10, 16, 1000000, 10000, 100000, 10000, 10000, 10000, 10000, 10000, 100000};
  Double_t xMin[nDim] = {0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 2.5, -0.5, -0.5, -0.5, -0.5};
  Double_t xMax[nDim] = {10.5, 16.5, 999999.5, 9999.5, 99999.5, 9999.5, 3.5, 0.5, 0.5, 0.5, 9.5};

  //THnSparse for Dimuon Analysis
  THnSparse *fhSigExt = new THnSparseT<TArrayF>("fhSigExt", "SignalExtraction", nDim, nBins, xMin, xMax);
  fhSigExt->Sumw2();

  Double_t TestInfos[nDim]; // Stockage des résultats du fit avant de remplir le THnSparse
  // Int_t testNumber = 1;

  if( Opt == "full" ){

    for(int multclass=1; multclass<=10; multclass++){
      Int_t testNumber = 1;

      TCanvas* canFit = new TCanvas("canFit","",2000,1500);
      SetCanvasStyle(canFit);
      canFit->Divide(5,4);

      TestInfos[0] = multclass;

      // fhDimuon->GetAxis(4)->SetRangeUser(MultClassesMaxRef[2*multclass-2]-1, MultClassesMaxRef[2*multclass-1]-1);
      fhDimuon->GetAxis(4)->SetRangeUser(MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);
      TH1D *histoInvMass = fhDimuon->Projection(3,"e");
      histoInvMass->Rebin(2);

      if ( ref == "min"){
        fhDimuon->GetAxis(4)->SetRangeUser(MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]);
        histoInvMass = fhDimuon->Projection(3,"e");
        histoInvMass->Rebin(2);
      }

      for(int range=0; range<nbRanges; range++){

        for(int sig=0; sig<nbSig; sig++){
          // sig = 0 : CB2  ;   sig = 1 : NA60

          for(int tail=0; tail<nbTails; tail++){
            // tail = 0 : MC  ;   tail = 1 : pp13TeV  ;  tail = 2 : Free Tail
            if ( tail == 2 ) fixTail = kFALSE;
            if ( sig == 1 && tail == 2) break; // Pas de Free tail pour la NA60
            if ( tail == 1 && sig == 1) break; // Car pas de tail pp13TeV pour NA60

            for(int bkg=0; bkg<nbBkg; bkg++){
              // bkg = 0 : VWG  ;   bkg = 1 : Pol1verPol2

              // Tentative de mettre les paramètres de test correspondant au testNumber sur l'axe des Tests
              TString testLabel;
              testLabel.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data());
              fhSigExt->GetAxis(1)->SetBinLabel(testNumber, testLabel);

              TString rangeName;
              rangeName.Form("%g_%g", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);

              canFit->cd(testNumber);

              cout << Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data()) << endl;

                if( ref == "max" ){
                  fitResults = FitOneHisto(histoInvMass, fSignal[sig], fBackground[bkg], "test", "test", Ranges[range][0], Ranges[range][1], GetTailsParam(tail, fSignal[sig]), fixTail, -1, 1, MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]);
                  cout << "----------------- Fit procedure " << testLabel << " done ----------------" << endl;
                }
                else if( ref == "min" ){
                  fitResults = FitOneHisto(histoInvMass, fSignal[sig], fBackground[bkg], "test", "test", Ranges[range][0], Ranges[range][1], GetTailsParam(tail, fSignal[sig]), fixTail, -1, 1, MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]);
                  cout << "----------------- Fit procedure " << testLabel << " done ----------------" << endl;
                }
                else {
                  cout << "Reference value not declared" << endl;
                }

                TestInfos[1] = testNumber;
                for (int i = 0; i < nDim-2; i++){
                  TestInfos[i+2] = fitResults[i];
                }

                fhSigExt->Fill(TestInfos);
                SetFitResults(rangeName, testLabel, fitResults);
                testNumber++;
              }
            }
          }
        }
        if( ref == "max" ){
          canFit->SaveAs(Form("Plots/SignalExtraction_MaxRef_MultRange_%.1f_%.1f.pdf", MultClassesMaxRef[2*multclass-2], MultClassesMaxRef[2*multclass-1]));
        }
        else if( ref == "min" ){
          canFit->SaveAs(Form("Plots/SignalExtraction_MinRef_MultRange_%.1f_%.1f.pdf", MultClassesMinRef[2*multclass-2], MultClassesMinRef[2*multclass-1]));
        }
        else{
          cout << "Reference value not declared, Canvas not saved" << endl;
        }
        delete histoInvMass;
        delete canFit;
      }
      fhSigExt->SaveAs("SignalExtractionResults.root");
  }

  if( Opt == "singlefit" ){

    fhDimuon->GetAxis(4)->SetRangeUser(multRangeMin,multRangeMax);
    TH1D *histoInvMass = fhDimuon->Projection(3,"e");
    histoInvMass->Rebin(2);

    TString testLabel;
    testLabel.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSig].Data(), arrayFunctionNames[fBkgd].Data(), fitRangeMin, fitRangeMax, tailNames[tailType].Data());

    TString rangeName;
    rangeName.Form("%g_%g", multRangeMin, multRangeMax);

    fitResults = FitOneHisto(histoInvMass, fSig, fBkgd, "test", "test", fitRangeMin, fitRangeMax, GetTailsParam(tailType, fSig), kTRUE,-1, 2, multRangeMin, multRangeMax);

    SetFitResults(rangeName, testLabel, fitResults);

    delete histoInvMass;

    // TString resultsName;
    // resultsName.Form("SignalExtractionResults_Mult_%.1f_%.1f_%s_%s_Range_%.1f_%.1f.root", multRangeMin, multRangeMax, arrayFunctionNames[fSig].Data(), arrayFunctionNames[fBkgd].Data(), fitRangeMin, fitRangeMax);
    // fhSigExt->SaveAs(resultsName);
    // TCanvas* canSigExt = new TCanvas("canSigExt","",2000,1500);
    // SetCanvasStyle(canSigExt);
    // fhSigExt->GetAxis(0)->SetRange(1.,1.);
    // TH2D *hSigExtJpsi = fhSigExt->Projection(2, 1,"e");
    // hSigExtJpsi->SetMarkerStyle(20);
    // hSigExtJpsi->Draw();
  }

  if( Opt == "integrated" ){

      Int_t testNumber = 1;

      TCanvas* canFit = new TCanvas("canFit","",2000,1500);
      SetCanvasStyle(canFit);
      canFit->Divide(5,4);

      if (ref == "min"){
        multRangeMin = MultClassesMinRef[0];
        multRangeMax = MultClassesMinRef[19];
      }

      fhDimuon->GetAxis(4)->SetRangeUser(multRangeMin,multRangeMax);
      TH1D *histoInvMass = fhDimuon->Projection(3,"e");
      // histoInvMass->Rebin(2);

      for(int range=0; range<nbRanges; range++){

        for(int sig=0; sig<nbSig; sig++){
          // sig = 0 : CB2  ;   sig = 1 : NA60

          for(int tail=0; tail<nbTails; tail++){
            // tail = 0 : MC  ;   tail = 1 : pp13TeV  ;  tail = 2 : Free Tail
            if ( tail == 2 ) {
              fixTail = kFALSE;
            }else{
              fixTail = kTRUE;
            }

            if ( sig == 1 && tail == 2) break; // Pas de Free tail pour la NA60
            if ( tail == 1 && sig == 1) break; // Car pas de tail pp13TeV pour NA60

            for(int bkg=0; bkg<nbBkg; bkg++){
              // bkg = 0 : VWG  ;   bkg = 1 : Pol1verPol2
              // if ( sig == 0 && bkg == 1 && range == 0 && tail == 1) break; // fit failed for test CB2_POL1POL2_Range_2.0_5.0_DataTail
              // Tentative de mettre les paramètres de test correspondant au testNumber sur l'axe des Tests
              TString testLabel;
              testLabel.Form("%s_%s_Range_%1.1f_%1.1f_%s", arrayFunctionNames[fSignal[sig]].Data(), arrayFunctionNames[fBackground[bkg]].Data(), Ranges[range][0], Ranges[range][1], tailNames[tail].Data());
              fhSigExt->GetAxis(1)->SetBinLabel(testNumber, testLabel);

              TString rangeName;
              rangeName.Form("%g_%g", multRangeMin, multRangeMax);

              canFit->cd(testNumber);

                  fitResults = FitOneHisto(histoInvMass, fSignal[sig], fBackground[bkg], "test", "test", Ranges[range][0], Ranges[range][1], GetTailsParam(tail, fSignal[sig]), fixTail, -1, 1, multRangeMin, multRangeMax);
                  cout << "----------------- Fit procedure " << testLabel << " done ----------------" << endl;

                SetFitResults(rangeName, testLabel, fitResults);
                testNumber++;
              }
            }
          }
        }
        if( ref == "max" ){
          canFit->SaveAs("Plots/SignalExtraction_refMax_Integrated.pdf");
        }
        else if( ref == "min" ){
          canFit->SaveAs("Plots/SignalExtraction_refMin_Integrated.pdf");
        }
        else{
          cout << "Reference value not declared, Canvas not saved" << endl;
        }
        delete histoInvMass;
        delete canFit;
      }
      fhSigExt->SaveAs("SignalExtractionResults_Integrated.root");
// museau, rideau !
}

std::vector<double> GetTailsParam(int iTail, int function){
  std::vector<double> tailsParams;
  for (int i=0; i<8; i++){
    if( iTail == 0) {
        if (function == kCB21S){
          tailsParams.push_back(CB2TailMC[i]);
        }
        else if (function == kNA601S) {
          tailsParams.push_back(NA60TailMC[i]);
        }
    }
    else if ( iTail == 1 && function == kCB21S){
      tailsParams.push_back(CB2Tailpp[i]);
    }
    else if ( iTail == 2 && function == kCB21S){
      tailsParams.push_back(CB2TailMC[i]);
    }
    else {
      cout << "------- Fail in setting fit parameters -------" << endl;
    }
  }
  return tailsParams;
}
