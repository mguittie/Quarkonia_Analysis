#include "FitOneHisto.C"
#include "THnSparse.h"

void SignalExtraction(TString fileName = "AnalysisResults.root", Double_t fitRangeMin = 2.4, Double_t fitRangeMax = 4.7, int fSig = kCB21S, int fBkgd = kVWG)
{
  // TFile *inputFile = new TFile("ResultsDimuon.root");
  // TH2F *histo2DDimuonInvMassVsPt = ((TH2F*) inputFile->Get("histo2DDimuonInvMassVsPt"));
  // TH1D *histoInvMass = histo2DDimuonInvMassVsPt->ProjectionY("InvMassWeighted",0,400);
  //

  // TFile *inputFile = new TFile("HistosfromTree_LHC18qr.root");
  // TH1D *histoInvMass = ((TH1D*) inputFile->Get("hMassOS_CMUL7_2m"));
  // histoInvMass->Rebin(2);

  Double_t MultClasses[20] = {1.0, 8.0, 9.0, 14.0, 15.0, 20.0, 21.0, 25.0, 26.0, 33.0, 34.0, 41.0, 42.0, 50.0, 51.0, 60.0, 61.0, 80.0, 81.0, 115.0};

  // prepare binning for THnSparse for Signal Extraction
  // 0: MultClass
  // 1: Test Number
  // 2: NJpsi
  // 3: Error Jpsi
  // 4: NPsi2S
  // 5: Error Psi2S

  const Int_t nDim = 6;
  Int_t nBins[nDim] = {10, 4, 100000, 10000, 100000, 10000};
  Double_t xMin[nDim] = {0.5, 0.5, -0.5, -0.5, -0.5, -0.5};
  Double_t xMax[nDim] = {10.5, 4.5, 999999.5, 9999.5, 99999.5, 9999.5};

  //THnSparse for Dimuon Analysis
  THnSparse *fhSigExt = new THnSparseT<TArrayF>("fhSigExt", "SignalExtraction", nDim, nBins, xMin, xMax);
  fhSigExt->Sumw2();

  TFile *file = new TFile(fileName.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileName.Data());
    return;
  }
  THnSparse *fhDimuon = static_cast<THnSparse*>(file->FindObjectAny("CMUL_fhDimuon"));
  if (!fhDimuon) return;

  std::vector<double> tailsParams;
  std::vector<double> fitResults;
  //----------------------------------------------------------------------------------------------//
  //Initialise the roofit variables and import the histogram:
  TCanvas* canFit = new TCanvas("canFit","",2000,1500);
  SetCanvasStyle(canFit);
  canFit->Divide(4,3);
  // canFit->SetLogy();
  //----------------------------------------------------------------------------------------------//
  Double_t TestInfos[nDim] = {0};

  for (int multclass = 1; multclass<=10; multclass ++) {
    canFit->cd(multclass);
    fhDimuon->GetAxis(4)->SetRangeUser(MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    TH1D *histoInvMass = fhDimuon->Projection(3,"e");

    fitResults = FitOneHisto(histoInvMass, kCB21S, kVWG, "test","test", 2.0, 5.0, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    TestInfos[0] = multclass;
    TestInfos[1] = 1;
    for (int i = 0; i <= 3; i++) {
      TestInfos[i+2] = fitResults[i];
    }
    fhSigExt->Fill(TestInfos);

    fitResults = FitOneHisto(histoInvMass, kCB21S, kVWG, "test","test", 2.2, 4.5, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    TestInfos[1] = 2;
    for (int i = 0; i <= 3; i++) {
      TestInfos[i+2] = fitResults[i];
    }
    fhSigExt->Fill(TestInfos);

    // fitResults = FitOneHisto(histoInvMass, kNA601S, kVWG, "test","test", 2.3, 5.5, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // TestInfos[1] = 3;
    // for (int i = 0; i <= 3; i++) {
    //   TestInfos[i+2] = fitResults[i];
    // }
    // fhSigExt->Fill(TestInfos);
    //
    // fitResults = FitOneHisto(histoInvMass, kNA601S, kVWG, "test","test", 2.0, 5.2, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // TestInfos[1] = 4;
    // for (int i = 0; i <= 3; i++) {
    //   TestInfos[i+2] = fitResults[i];
    // }
    // fhSigExt->Fill(TestInfos);
  }

  TCanvas* canSigExt = new TCanvas("canSigExt","",2000,1500);
  SetCanvasStyle(canSigExt);
  fhSigExt->GetAxis(0)->SetRange(1.,1.);
  TH2D *hSigExtJpsi = fhSigExt->Projection(2, 1,"e");
  hSigExtJpsi->SetMarkerStyle(20);
  hSigExtJpsi->Draw();


  canFit->SaveAs(Form("Plots/%s_%s_Range_%2.2f_%2.2f.pdf",arrayFunctionNames[fSig].Data(),arrayFunctionNames[fBkgd].Data(),fitRangeMin,fitRangeMax));


  // fhDimuon->GetAxis(4)->SetRangeUser(81.0, 115.0);
  //
  // TH1D *histoInvMass = fhDimuon->Projection(3,"e");
  // // histoInvMass->SetRangeUser(2.0,5.0);
  //
  // std::vector<double> tailsParams;
  // //FitOneHisto(TH1 *histoInvMass,int fSig, int fBkgd, TString strHistoTitle,TString strRange, Double_t fitRangeMin, Double_t fitRangeMax, std::vector<double> tailsParams, Bool_t fixCBTails, Double_t jpsiMCWidth, int rebiningFator)
  // FitOneHisto(histoInvMass,kCB21S, kVWG, "test","test", 2.4,4.7,tailsParams,kTRUE,-1, 2);

}
