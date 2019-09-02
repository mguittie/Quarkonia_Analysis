#include "SetAndGetFitResults.C"
#include "Common.C"

std::vector<double> GetMeanAndRMS(TH1F *histo){
  Double_t mean=0;
  Double_t RMS=0;
  int numberOfBins = histo->GetNbinsX();
  for(int iBin=1;iBin<=numberOfBins;iBin++){
    mean += histo->GetBinContent(iBin);
  }
  mean /= numberOfBins;
  for(int iBin=1;iBin<=numberOfBins;iBin++){
    RMS += (mean - histo->GetBinContent(iBin))*(mean - histo->GetBinContent(iBin));
  }
  RMS = TMath::Sqrt(RMS/numberOfBins);

  std::vector<double> vector;
  vector.push_back(mean);
  vector.push_back(RMS);
  return vector;
}

Bool_t IsTestFailed(TH1F *histo,int iTest){
  std::vector<double> chi2MeanAndRMS = GetMeanAndRMS(histo);
  Double_t chi2ForThisTest = histo->GetBinContent(iTest);
  if(chi2ForThisTest > (chi2MeanAndRMS[0]+2*chi2MeanAndRMS[1])) return kTRUE;
  return kFALSE;
}

void CreatePadAndDrawHisto(TH1F *histo,Double_t histoMin,Double_t histoMax, Double_t lineDownAt, Double_t lineUpAt,Double_t bottomMargin, Double_t topMargin,Double_t yMin, Double_t yMax, int color, int markerStyle, const char *name, Double_t numberOfTests, TString strYaxisTitle,bool isGridY){
  TPad *pad = new TPad(Form("padControl_%s",name), Form("padControl_%s",name),0,yMin,1,yMax,0);
  pad->SetBottomMargin(bottomMargin);
  pad->SetTopMargin(topMargin);
  pad->Draw();
  if(isGridY) pad->SetGridy();
  pad->cd();

  if( (histoMin != -1) && (histoMax != -1)){
    histo->SetMinimum(histoMin);
    histo->SetMaximum(histoMax);
  }

  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerSize(1);
  histo->GetYaxis()->SetTitle(strYaxisTitle);
  histo->GetYaxis()->CenterTitle();
  histo->GetYaxis()->SetTitleOffset(0.15);
  histo->GetYaxis()->SetTitleSize(0.25);
  histo->GetYaxis()->SetLabelSize(0.1);
  histo->GetYaxis()->SetTickLength(0.01);
  histo->GetYaxis()->SetNdivisions(4);
  histo->SetMarkerStyle(markerStyle);
  histo->Draw("p");


  TLine *lineUp = new TLine(0,lineDownAt,numberOfTests,lineDownAt);
  lineUp->SetLineColor(color);
  lineUp->SetLineWidth(1);
  lineUp->SetLineStyle(2);


  TLine *lineDown = new TLine(0,lineUpAt,numberOfTests,lineUpAt);
  lineDown->SetLineColor(color);
  lineDown->SetLineWidth(1);
  lineDown->SetLineStyle(2);


  if(!( (lineDownAt == -1) && lineUpAt == -1 )){
    lineDown->Draw("same");
    lineUp->Draw("same");
  }


}

void DrawHisto(TH1F *histo,Double_t histoMin,Double_t histoMax, Double_t lineDownAt, Double_t lineUpAt, int color, int markerStyle, const char *name, Double_t numberOfTests, TString strYaxisTitle){

  if( (histoMin != -1) && (histoMax != -1)){
    histo->SetMinimum(histoMin);
    histo->SetMaximum(histoMax);
  }

  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerSize(1);
  histo->GetYaxis()->SetTitle(strYaxisTitle);
  histo->GetYaxis()->CenterTitle();
  histo->GetYaxis()->SetTitleOffset(0.15);
  histo->GetYaxis()->SetTitleSize(0.25);
  histo->GetYaxis()->SetLabelSize(0.1);
  histo->GetYaxis()->SetTickLength(0.01);
  histo->GetYaxis()->SetNdivisions(4);
  histo->SetMarkerStyle(markerStyle);
  histo->Draw("p");


  TLine *lineUp = new TLine(0,lineDownAt,numberOfTests,lineDownAt);
  lineUp->SetLineColor(color);
  lineUp->SetLineWidth(1);
  lineUp->SetLineStyle(2);


  TLine *lineDown = new TLine(0,lineUpAt,numberOfTests,lineUpAt);
  lineDown->SetLineColor(color);
  lineDown->SetLineWidth(1);
  lineDown->SetLineStyle(2);


  if(!( (lineDownAt == -1) && lineUpAt == -1 )){
    lineDown->Draw("same");
    lineUp->Draw("same");
  }
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

void DrawAndRMS(TString rangeName, std::vector<TString> vectorTests)
{
  ///////////////////////////////////////////
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  ///////////////////////////////////////////


  Int_t numberOfTests = (int)vectorTests.size();


  TH1F *histoFitResults[numberOfFitVariables];

  for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    histoFitResults[iVariable] = new TH1F(Form("histo%s%s",arrayFitVariableNames[iVariable].Data(),rangeName.Data()),"",numberOfTests,0,numberOfTests);
  }

  for(Int_t iTest=0;iTest<numberOfTests;iTest++){
    TString testName = vectorTests[iTest];

    std::vector<double> testFitResults, testFitResults_err;
    GetFitResults(rangeName, testName, testFitResults, testFitResults_err);

    for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++){
      histoFitResults[iVariable]->SetBinContent(iTest+1,testFitResults[iVariable]);
      histoFitResults[iVariable]->SetBinError(iTest+1,testFitResults_err[iVariable]);
      histoFitResults[iVariable]->GetXaxis()->SetBinLabel(iTest+1,testName);
    }

  }

  Double_t numberOfJpsi =0;
  Double_t numberOfJpsi_stat_err =0;
  Double_t numberOfJpsi_sys_err =0;
  Double_t numberOfWeights=0;

  for(int iTest =1;iTest<=numberOfTests;iTest++){
    Double_t testWeight = 1;
    TString strTestLabel = histoFitResults[0]->GetXaxis()->GetBinLabel(iTest);
    if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
    numberOfJpsi += testWeight*histoFitResults[0]->GetBinContent(iTest,1);
    numberOfJpsi_stat_err += testWeight*histoFitResults[0]->GetBinError(iTest,1);
    numberOfWeights += testWeight;
  }

  numberOfJpsi /= numberOfWeights;
  numberOfJpsi_stat_err /= numberOfWeights;


  for(int iTest =1;iTest<=numberOfTests;iTest++){
    Double_t testWeight = 1;
    TString strTestLabel = histoFitResults[0]->GetXaxis()->GetBinLabel(iTest);
    if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
    numberOfJpsi_sys_err += testWeight*TMath::Power( (histoFitResults[0]->GetBinContent(iTest,1) - numberOfJpsi),2);
  }
  numberOfJpsi_sys_err /= numberOfWeights;
  numberOfJpsi_sys_err = TMath::Sqrt(numberOfJpsi_sys_err);


  TCanvas* canFitResult = new TCanvas(Form("canFitResult%s",rangeName.Data()),Form("canFitResult%s",rangeName.Data()),1200,1200) ;

  TPad *padMain = new TPad("padMain", "padMain",0,0.6,1,1,0);
  padMain->SetBottomMargin(0.4);
  padMain->Draw();
  padMain->cd();


  histoFitResults[0]->SetMinimum(numberOfJpsi-10*numberOfJpsi_sys_err);
  histoFitResults[0]->SetMaximum(numberOfJpsi+10*numberOfJpsi_sys_err);
  histoFitResults[0]->SetMarkerStyle(20);
  histoFitResults[0]->SetMarkerColor(4);
  histoFitResults[0]->SetLineColor(1);
  histoFitResults[0]->GetYaxis()->SetTitle("N_{J/#psi}");
  histoFitResults[0]->GetYaxis()->SetTitleSize(0.05);
  histoFitResults[0]->GetYaxis()->SetTitleOffset(0.7);
  histoFitResults[0]->GetYaxis()->SetLabelSize(0.035);
  histoFitResults[0]->Draw("p");


  TLine *lineAverage = new TLine(0,numberOfJpsi,numberOfTests,numberOfJpsi);
  lineAverage->SetLineColor(4);
  lineAverage->SetLineWidth(1);
  lineAverage->Draw("");

  TLegend* legRMSLines = new TLegend(0.142, 0.058, 0.85, 0.20);
  legRMSLines->SetFillStyle(0);
  legRMSLines->SetLineColorAlpha(0,0);
  legRMSLines->SetTextColor(kBlack);
  legRMSLines->SetMargin(0.1);
  DrawDownAndUpLines(0, numberOfTests, numberOfJpsi-numberOfJpsi_sys_err, numberOfJpsi+numberOfJpsi_sys_err, 2, 2, legRMSLines,"#pm RMS");
  DrawDownAndUpLines(0, numberOfTests, numberOfJpsi-2*numberOfJpsi_sys_err, numberOfJpsi+2*numberOfJpsi_sys_err, 3, 5, legRMSLines,"#pm 2*RMS");
  legRMSLines->Draw();


  TLegend* legJpsiNumber = new TLegend(0.07, 0.70, 0.59, 0.84);
  legJpsiNumber->SetFillStyle(0);
  legJpsiNumber->SetLineColorAlpha(0,0);
  legJpsiNumber->SetTextColor(kBlack);
  legJpsiNumber->SetMargin(0.1);
  legJpsiNumber->AddEntry("NULL",Form("%s, N_{J/#psi} = %d #pm %d (stat) #pm %d (sys) (%2.2f %%)",rangeName.Data(), (int)numberOfJpsi,(int)numberOfJpsi_stat_err ,(int)numberOfJpsi_sys_err,100*numberOfJpsi_sys_err/numberOfJpsi),"");
  legJpsiNumber->Draw();


  Double_t padPosition = 0.6;
  Double_t padHeight = 0.125;
  int padCounter =0;

  canFitResult->cd();
  std::vector<double> fitWidthsMeanAndRMS = GetMeanAndRMS(histoFitResults[3]);
  CreatePadAndDrawHisto(histoFitResults[3],-1,-1, fitWidthsMeanAndRMS[0] - fitWidthsMeanAndRMS[1],fitWidthsMeanAndRMS[0] + fitWidthsMeanAndRMS[1], 0,0,padPosition-(padCounter+1)*padHeight,padPosition-(padCounter*padHeight),2,20,"jpsiWidth",numberOfTests,"#sigma_{J/#Psi}",false);
  padCounter++;


  canFitResult->cd();
  std::vector<double> fitMeanPosMeanAndRMS = GetMeanAndRMS(histoFitResults[2]);
  CreatePadAndDrawHisto(histoFitResults[2],-1,-1, fitMeanPosMeanAndRMS[0] - fitMeanPosMeanAndRMS[1],fitMeanPosMeanAndRMS[0] + fitMeanPosMeanAndRMS[1], 0,0,padPosition-(padCounter+1)*padHeight,padPosition-(padCounter*padHeight),3,20,"jpsiMeanPos",numberOfTests,"M_{J/#Psi}",false);
  padCounter++;


  canFitResult->cd();

  std::vector<double> fitChisquareMeanAndRMS = GetMeanAndRMS(histoFitResults[4]);
  CreatePadAndDrawHisto(histoFitResults[4],0.02,3.5, fitChisquareMeanAndRMS[0] - fitChisquareMeanAndRMS[1],fitChisquareMeanAndRMS[0] + fitChisquareMeanAndRMS[1], 0.696,0,0,padPosition-(padCounter*padHeight),7,20,"fitChisquare",numberOfTests,"#chi^{2}/ndf",false);
  histoFitResults[kChi2PerTest]->LabelsOption("V");
  histoFitResults[kChi2PerTest]->GetYaxis()->SetLabelSize(0.05);
  histoFitResults[kChi2PerTest]->GetXaxis()->SetLabelSize(0.05);
  histoFitResults[kChi2PerTest]->GetYaxis()->SetTitleSize(0.07);
  histoFitResults[kChi2PerTest]->GetYaxis()->SetTitleOffset(0.5);


}

void CreateHistosSigExtAndRMS(TString rangeName, std::vector<TString> vectorTests, int variableID)
{
  ///////////////////////////////////////////
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  ///////////////////////////////////////////


  Int_t numberOfTests = (int)vectorTests.size();

  if(variableID == 5){
    TH1F *histoFitResults = new TH1F(Form("histoPsi2SOverJPsi_%s",rangeName.Data()),"",numberOfTests,0,numberOfTests);

    for(Int_t iTest=0;iTest<numberOfTests;iTest++){
      TString testName = vectorTests[iTest];

      std::vector<double> testFitResults, testFitResults_err;
      GetFitResults(rangeName, testName, testFitResults, testFitResults_err);

      histoFitResults->SetBinContent(iTest+1,testFitResults[1]/testFitResults[0]);
      histoFitResults->SetBinError(iTest+1,testFitResults[1]/testFitResults[0]*sqrt(TMath::Power(testFitResults_err[1]/testFitResults[1],2) + TMath::Power(testFitResults_err[0]/testFitResults[0],2)));
      histoFitResults->GetXaxis()->SetBinLabel(iTest+1,testName);
    }
    Double_t varValue =0;
    Double_t varValue_stat_err =0;
    Double_t varValue_sys_err =0;
    Double_t numberOfWeights=0;

    for(int iTest =1;iTest<=numberOfTests;iTest++){
      Double_t testWeight = 1;
      TString strTestLabel = histoFitResults->GetXaxis()->GetBinLabel(iTest);
      if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
      varValue += testWeight*histoFitResults->GetBinContent(iTest,1);
      varValue_stat_err += testWeight*histoFitResults->GetBinError(iTest,1);
      numberOfWeights += testWeight;
    }

    varValue /= numberOfWeights;
    varValue_stat_err /= numberOfWeights;


    for(int iTest =1;iTest<=numberOfTests;iTest++){
      Double_t testWeight = 1;
      TString strTestLabel = histoFitResults->GetXaxis()->GetBinLabel(iTest);
      if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
      varValue_sys_err += testWeight*TMath::Power( (histoFitResults->GetBinContent(iTest,1) - varValue),2);
    }
    varValue_sys_err /= numberOfWeights;
    varValue_sys_err = TMath::Sqrt(varValue_sys_err);


    histoFitResults->SetMinimum(varValue-10*varValue_sys_err);
    histoFitResults->SetMaximum(varValue+10*varValue_sys_err);
    histoFitResults->SetMarkerStyle(20);
    histoFitResults->SetMarkerColor(2);
    histoFitResults->SetLineColor(1);
    histoFitResults->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}");
    histoFitResults->GetYaxis()->SetTitleSize(0.05);
    histoFitResults->GetYaxis()->SetTitleOffset(1);
    histoFitResults->GetYaxis()->SetLabelSize(0.02);
    histoFitResults->Draw("p");


    TLine *lineAverage = new TLine(0,varValue,numberOfTests,varValue);
    lineAverage->SetLineColor(4);
    lineAverage->SetLineWidth(1);
    lineAverage->Draw("");


    TLegend* legRMSLines = new TLegend(0.142, 0.22, 0.6, 0.32);
    legRMSLines->SetFillStyle(0);
    legRMSLines->SetLineColorAlpha(0,0);
    legRMSLines->SetTextColor(kBlack);
    legRMSLines->SetMargin(0.1);
    DrawDownAndUpLines(0, numberOfTests, varValue-varValue_sys_err, varValue+varValue_sys_err, 4, 4, legRMSLines,"#pm RMS");
    DrawDownAndUpLines(0, numberOfTests, varValue-2*varValue_sys_err, varValue+2*varValue_sys_err, 3, 5, legRMSLines,"#pm 2*RMS");
    legRMSLines->Draw();

    TLegend* legPartNumber = new TLegend(0.05, 0.65, 0.8, 0.9);
    legPartNumber->SetFillStyle(0);
    legPartNumber->SetLineColorAlpha(0,0);
    legPartNumber->SetTextColor(kBlack);
    legPartNumber->SetMargin(0.1);
    legPartNumber->AddEntry("NULL",Form("Mult_%s, N_{#psi(2S)}/N_{J/#psi} = %f #pm %f (stat) #pm %f (sys) (%2.2f %%)",rangeName.Data(),varValue,varValue_stat_err ,varValue_sys_err,100*varValue_sys_err/varValue),"");
    legPartNumber->Draw();

    return;
  }

  TH1F *histoFitResults = new TH1F(Form("histo%s%s",arrayFitVariableNames[variableID].Data(),rangeName.Data()),"",numberOfTests,0,numberOfTests);

  // for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++){
  //   histoFitResults[iVariable] = new TH1F(Form("histo%s%s",arrayFitVariableNames[iVariable].Data(),rangeName.Data()),"",numberOfTests,0,numberOfTests);
  // }

  for(Int_t iTest=0;iTest<numberOfTests;iTest++){
    TString testName = vectorTests[iTest];

    std::vector<double> testFitResults, testFitResults_err;
    GetFitResults(rangeName, testName, testFitResults, testFitResults_err);

    // for(int iVariable=0;iVariable<numberOfFitVariables;iVariable++){
    //   histoFitResults[iVariable]->SetBinContent(iTest+1,testFitResults[iVariable]);
    //   histoFitResults[iVariable]->SetBinError(iTest+1,testFitResults_err[iVariable]);
    //   histoFitResults[iVariable]->GetXaxis()->SetBinLabel(iTest+1,testName);
    // }

    histoFitResults->SetBinContent(iTest+1,testFitResults[variableID]);
    histoFitResults->SetBinError(iTest+1,testFitResults_err[variableID]);
    histoFitResults->GetXaxis()->SetBinLabel(iTest+1,testName);

  }

  Double_t varValue =0;
  Double_t varValue_stat_err =0;
  Double_t varValue_sys_err =0;
  Double_t numberOfWeights=0;

  for(int iTest =1;iTest<=numberOfTests;iTest++){
    Double_t testWeight = 1;
    TString strTestLabel = histoFitResults->GetXaxis()->GetBinLabel(iTest);
    if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
    varValue += testWeight*histoFitResults->GetBinContent(iTest,1);
    varValue_stat_err += testWeight*histoFitResults->GetBinError(iTest,1);
    numberOfWeights += testWeight;
  }

  varValue /= numberOfWeights;
  varValue_stat_err /= numberOfWeights;


  for(int iTest =1;iTest<=numberOfTests;iTest++){
    Double_t testWeight = 1;
    TString strTestLabel = histoFitResults->GetXaxis()->GetBinLabel(iTest);
    if(strTestLabel.Contains("DataAt13TeV")) testWeight =2;
    varValue_sys_err += testWeight*TMath::Power( (histoFitResults->GetBinContent(iTest,1) - varValue),2);
  }
  varValue_sys_err /= numberOfWeights;
  varValue_sys_err = TMath::Sqrt(varValue_sys_err);


  histoFitResults->SetMinimum(varValue-10*varValue_sys_err);
  histoFitResults->SetMaximum(varValue+10*varValue_sys_err);
  histoFitResults->SetMarkerStyle(20);
  histoFitResults->SetMarkerColor(2);
  histoFitResults->SetLineColor(1);
  histoFitResults->GetYaxis()->SetTitle(Form("%s",variableNames[variableID]));
  histoFitResults->GetYaxis()->SetTitleSize(0.05);
  histoFitResults->GetYaxis()->SetTitleOffset(1);
  histoFitResults->GetYaxis()->SetLabelSize(0.02);
  histoFitResults->Draw("p");


  TLine *lineAverage = new TLine(0,varValue,numberOfTests,varValue);
  lineAverage->SetLineColor(4);
  lineAverage->SetLineWidth(1);
  lineAverage->Draw("");


  if (variableID <= 1) {
    TLegend* legRMSLines = new TLegend(0.142, 0.22, 0.6, 0.32);
    legRMSLines->SetFillStyle(0);
    legRMSLines->SetLineColorAlpha(0,0);
    legRMSLines->SetTextColor(kBlack);
    legRMSLines->SetMargin(0.1);
    DrawDownAndUpLines(0, numberOfTests, varValue-varValue_sys_err, varValue+varValue_sys_err, 4, 4, legRMSLines,"#pm RMS");
    DrawDownAndUpLines(0, numberOfTests, varValue-2*varValue_sys_err, varValue+2*varValue_sys_err, 3, 5, legRMSLines,"#pm 2*RMS");
    legRMSLines->Draw();

    TLegend* legPartNumber = new TLegend(0.05, 0.65, 0.8, 0.9);
    legPartNumber->SetFillStyle(0);
    legPartNumber->SetLineColorAlpha(0,0);
    legPartNumber->SetTextColor(kBlack);
    legPartNumber->SetMargin(0.1);
    legPartNumber->AddEntry("NULL",Form("Mult_%s, %s = %d #pm %d (stat) #pm %d (sys) (%2.2f %%)",rangeName.Data(), variableNames[variableID],(int)varValue,(int)varValue_stat_err ,(int)varValue_sys_err,100*varValue_sys_err/varValue),"");
    legPartNumber->Draw();
  }
  else {
    TLegend* legPartNumber = new TLegend(0.4, 0.70, 0.6, 0.75);
    legPartNumber->SetFillStyle(0);
    legPartNumber->SetLineColorAlpha(0,0);
    legPartNumber->SetTextColor(kBlack);
    legPartNumber->SetMargin(0.1);
    legPartNumber->AddEntry("NULL",Form("Mult_%s",rangeName.Data()),"");
    legPartNumber->Draw();
  }
}
