#include "FitOneHisto.C"

void FitTest(TString fileName = "AnalysisResults.root", Double_t fitRangeMin = 2.4, Double_t fitRangeMax = 4.7, int fSig = kCB21S, int fBkgd = kVWG)
{
  // TFile *inputFile = new TFile("ResultsDimuon.root");
  // TH2F *histo2DDimuonInvMassVsPt = ((TH2F*) inputFile->Get("histo2DDimuonInvMassVsPt"));
  // TH1D *histoInvMass = histo2DDimuonInvMassVsPt->ProjectionY("InvMassWeighted",0,400);
  //

  // TFile *inputFile = new TFile("HistosfromTree_LHC18qr.root");
  // TH1D *histoInvMass = ((TH1D*) inputFile->Get("hMassOS_CMUL7_2m"));
  // histoInvMass->Rebin(2);

  Double_t MultClasses[20] = {1.0, 8.0, 9.0, 14.0, 15.0, 20.0, 21.0, 25.0, 26.0, 33.0, 34.0, 41.0, 42.0, 50.0, 51.0, 60.0, 61.0, 80.0, 81.0, 115.0};
  Double_t NA60TailMC[8] = {-1.01, 0.002886, 0.4387, 0.2225, 2.164, 0.1924, 1.284, 0.06384};

  std::vector<double> tailsParams;

  for (int i=0; i<8; i++){
    tailsParams.push_back(NA60TailMC[i]);
  }

  TFile *file = new TFile(fileName.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileName.Data());
    return;
  }
  THnSparse *fhDimuon = static_cast<THnSparse*>(file->FindObjectAny("CMUL_fhDimuon"));
  if (!fhDimuon) return;

  std::vector<double> fitResults;

  //----------------------------------------------------------------------------------------------//
  //Initialise the roofit variables and import the histogram:
  TCanvas* canFit = new TCanvas("canFit","",2000,1500) ;
  SetCanvasStyle(canFit);
  canFit->Divide(4,3);
  // canFit->SetLogy();
  //----------------------------------------------------------------------------------------------//
  for (int multclass = 1; multclass<=1; multclass ++) {
    canFit->cd(multclass);
    fhDimuon->GetAxis(4)->SetRangeUser(MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // fhDimuon->GetAxis(4)->SetRangeUser(1.0, 115.0);
    TH1D *histoInvMass = fhDimuon->Projection(3,"e");

    // FitOneHisto(TH1 *histoInvMass,int fSig, int fBkgd, TString strHistoTitle,TString strRange, Double_t fitRangeMin, Double_t fitRangeMax, std::vector<double> tailsParams, Bool_t fixCBTails, Double_t jpsiMCWidth, int rebiningFator)
    // FitOneHisto(histoInvMass, kCB21S, kVWG, "test","test", 2.4, 4.7, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // fitResults = FitOneHisto(histoInvMass, kCB21S, kVWG, "test","test", 2.3, 5.5, tailsParams,kTRUE,-1, 2, 1.0, 115.0);
    fitResults = FitOneHisto(histoInvMass, fSig, fBkgd, "test","test", fitRangeMin, fitRangeMax, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // fitResults = FitOneHisto(histoInvMass, kCB21S, kVWG, "test","test", 2.0, 5.2, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // fitResults = FitOneHisto(histoInvMass, kNA601S, kVWG, "test","test", 2.3, 5.5, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);
    // fitResults = FitOneHisto(histoInvMass, kNA601S, kVWG, "test","test", 2.0, 5.2, tailsParams,kTRUE,-1, 2, MultClasses[2*multclass-2], MultClasses[2*multclass-1]);

    printf("---------------------------------------------------------------------------------------------------------------\n");
    printf("Nb JPsi = %f +/- %f and Nb Psi2S = %f +/- %f\n", fitResults[0], fitResults[1], fitResults[2], fitResults[3]);
  }

  // canFit->SaveAs(Form("Plots/%s_%s_Range_%2.2f_%2.2f.pdf",arrayFunctionNames[fSig].Data(),arrayFunctionNames[fBkgd].Data(),fitRangeMin,fitRangeMax));


  // fhDimuon->GetAxis(4)->SetRangeUser(81.0, 115.0);
  //
  // TH1D *histoInvMass = fhDimuon->Projection(3,"e");
  // // histoInvMass->SetRangeUser(2.0,5.0);
  //
  // std::vector<double> tailsParams;
  // //FitOneHisto(TH1 *histoInvMass,int fSig, int fBkgd, TString strHistoTitle,TString strRange, Double_t fitRangeMin, Double_t fitRangeMax, std::vector<double> tailsParams, Bool_t fixCBTails, Double_t jpsiMCWidth, int rebiningFator)
  // FitOneHisto(histoInvMass,kCB21S, kVWG, "test","test", 2.4,4.7,tailsParams,kTRUE,-1, 2);

}
