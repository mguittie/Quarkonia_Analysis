/*
An interface to ease getting and setting the fit results/range/test
*/
//needed root classes
#include "Riostream.h"
#include <TFile.h>
#include "TROOT.h"
#include "TH3.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"


const Int_t maxNumberOfRanges = 1000;
const Int_t maxNumberOfTests = 100;
const Int_t numberOfFitVariables = 5;
TString arrayFitVariableNames[numberOfFitVariables] = {"NPsi", "NPsi2s", "Mean", "sigma", "chi2"};

//-------------------------------------------------------------------------------------------------------------------------------------//
//Return the histogram which contains the values.
TH3F *GetFitResultsHisto(Bool_t createIfNotFound = kFALSE, TString rangeName = ""){
  TH3F *histoFitResults;
  if (rangeName == "1_115" || rangeName == "1_83"){
    TFile *outputFile = new TFile("FitResults_Integrated.root");
    histoFitResults =  ((TH3F*) outputFile->Get("histoFitResults"));
  }
  else {
    TFile *outputFile = new TFile("FitResults.root");
    histoFitResults =  ((TH3F*) outputFile->Get("histoFitResults"));
  }

  if(!histoFitResults && createIfNotFound){
    cout<<"fit results file not found: creating a new histogram"<<endl;
    histoFitResults =  new TH3F("histoFitResults","",maxNumberOfRanges,0,maxNumberOfRanges,maxNumberOfTests,0,maxNumberOfTests,numberOfFitVariables,0,numberOfFitVariables);
  }
  return histoFitResults;
}
//-------------------------------------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------------//
void SetFitResults(TString rangeName, TString testName, std::vector<double> fitResults){

  TH3F *histoFitResults = GetFitResultsHisto(kTRUE, rangeName);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  int xBinInHistoFitResults = histoFitResults->GetXaxis()->FindFixBin(rangeName.Data());
  if(xBinInHistoFitResults == -1)
  {
    //This used to work with FindLastBinAbove(), no clue on why it is not working! Replace with a loop over the x-axis labels.
    Int_t nonFilledBin;
    for(nonFilledBin=0;nonFilledBin<maxNumberOfRanges;nonFilledBin++){
      TString strTempoLabel = histoFitResults->GetXaxis()->GetBinLabel(nonFilledBin+1);
      if(strTempoLabel == ""){
        break;
      }
    }
    histoFitResults->GetXaxis()->SetBinLabel(nonFilledBin+1,rangeName);
    xBinInHistoFitResults = nonFilledBin+1;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  int yBinInHistoFitResults = histoFitResults->GetYaxis()->FindFixBin(testName.Data());
  if(yBinInHistoFitResults == -1)
  {
    //This used to work with FindLastBinAbove(), no clue on why it is not working! Replace with a loop over the x-axis labels.
    Int_t nonFilledBin;
    for(nonFilledBin=0;nonFilledBin<maxNumberOfTests;nonFilledBin++){
      TString strTempoLabel = histoFitResults->GetYaxis()->GetBinLabel(nonFilledBin+1);
      if(strTempoLabel == ""){
        break;
      }
    }
    histoFitResults->GetYaxis()->SetBinLabel(nonFilledBin+1,testName);
    yBinInHistoFitResults = nonFilledBin+1;
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


  for(Int_t iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    histoFitResults->SetBinContent(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1, fitResults[iVariable]);
    histoFitResults->SetBinError(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1, fitResults[iVariable+numberOfFitVariables]);
  }

  if (rangeName == "1_115" || rangeName == "1_83"){
    TFile *outputFile = new TFile("FitResults_Integrated.root", "RECREATE");
    histoFitResults->Write();
    outputFile->Close();
  }
  else {
    TFile *outputFile = new TFile("FitResults.root", "RECREATE");
    histoFitResults->Write();
    outputFile->Close();
  }


}
//-------------------------------------------------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------------------------------------------------//
void GetFitResults(TString rangeName, TString testName, std::vector<double> &fitResults, std::vector<double> &fitResults_err){

  TH3F *histoFitResults = GetFitResultsHisto(kFALSE, rangeName);
  if(!histoFitResults){
    cout<<"Fit Results file not found"<<endl;
    return;
  }
  int xBinInHistoFitResults = histoFitResults->GetXaxis()->FindFixBin(rangeName.Data());
  if(xBinInHistoFitResults == -1)
  {
    cout<<"Range not found"<<endl;
    return;
  }
  int yBinInHistoFitResults = histoFitResults->GetYaxis()->FindFixBin(testName.Data());
  if(yBinInHistoFitResults == -1)
  {
    cout<<"test not found"<<endl;
    return;
  }

  for(Int_t iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    fitResults.push_back(histoFitResults->GetBinContent(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1));
    fitResults_err.push_back(histoFitResults->GetBinError(xBinInHistoFitResults, yBinInHistoFitResults, iVariable+1));
  }
}
//-------------------------------------------------------------------------------------------------------------------------------------//
