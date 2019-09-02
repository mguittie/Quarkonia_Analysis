#include "FitOneHisto.C"
#include "DrawAndRMS.C"

//----------------------------------------------------------------------------------------------//
//IO params
TString inputFileName  = "ResultsDimuon.root";
//----------------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------------------//
//Define here the pt ranges in which you want to do the fit
Double_t ptRanges[][2] = { {7,8} };
// Double_t ptRanges[][2] = { {0,15},{0,2},{2,4},{4,6},{6,15}};
int numberOfPtRanges = sizeof( ptRanges ) / sizeof( ptRanges[0] );
//----------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------//
//here the ingredients of the different tests
Double_t fitRanges[][2] = { {1.7,4.8},{2,4.4} };
Int_t numberOfFitRanges = sizeof( fitRanges ) / sizeof( fitRanges[0] );

Int_t sigFunctions[] = { kCB21S, kNA601S};
Int_t numberOfSigFunctions = sizeof( sigFunctions ) / sizeof( sigFunctions[0] );

Int_t tailsSets[] = { kMC, kData13TeV };
Int_t numberOfTailsSets = sizeof( tailsSets ) / sizeof( tailsSets[0] );

Int_t bkgdFunctions[] = { kPol1OverPol2, kVWG };
Int_t numberOfBkgdFunctions = sizeof( bkgdFunctions ) / sizeof( bkgdFunctions[0] );

Int_t rebiningFator = 3;
//----------------------------------------------------------------------------------------------//


void FitTest()
{

  TFile *inputFile = new TFile(inputFileName);

  std::vector<double> tailsParams;//null
  for(Int_t iPtRange=0;iPtRange<numberOfPtRanges;iPtRange++){

    Double_t pTMin = ptRanges[iPtRange][0];
    Double_t pTMax = ptRanges[iPtRange][1];

    TString rangeName;
    rangeName.Form("Pt%gto%g",pTMin,pTMax);//Update it with rapidity

    TH1F *histoInvMass = ((TH1F*) inputFile->Get(Form("histoDimuonInvMass_Pt%gto%g",pTMin,pTMax)));
    histoInvMass->Rebin(rebiningFator);

    TCanvas* canFit = new TCanvas("CanFit","",2000,1000) ;
    canFit->Divide(4,3);
    SetCanvasStyle(canFit);

    std::vector<TString> vectorTestNames;
    Int_t testCounter=1;
    for(Int_t iSigFunction=0;iSigFunction<numberOfSigFunctions;iSigFunction++){
      Int_t fSig = sigFunctions[iSigFunction];

      for(Int_t iBkgdFunction=0;iBkgdFunction<numberOfBkgdFunctions;iBkgdFunction++){
        Int_t fBkgd = bkgdFunctions[iBkgdFunction];

        for(Int_t iFitRanges=0;iFitRanges<numberOfFitRanges;iFitRanges++){

          Double_t fitRangeMin = fitRanges[iFitRanges][0];
          Double_t fitRangeMax = fitRanges[iFitRanges][1];

          for(Int_t iTailsSet=0;iTailsSet<numberOfTailsSets;iTailsSet++){
            Int_t fTailsSet = tailsSets[iTailsSet];

            if(fTailsSet == kData13TeV && fSig==kNA601S){
              continue;
            }
            std::vector<double> vectorTails = GetTails(fSig,fTailsSet, pTMin , pTMax, -4,-2.5);

            //----------------------------------------------------------------------------------------------//
            //Initialise the roofit variables and import the histogram:
            TString testName;
            testName.Form("sig%s_bkgd%s_tails%s_fitRanges%gto%g",arrayFunctionNames[fSig].Data(),arrayFunctionNames[fBkgd].Data(),arrayTailsSetNames[iTailsSet].Data(),fitRangeMin,fitRangeMax);
            vectorTestNames.push_back(testName);

            canFit->cd(testCounter);
            // canFit->SetLogy();
            //----------------------------------------------------------------------------------------------//

            std::vector<double> testFitResults, testFitResults_err;
            FitOneHisto(histoInvMass, testFitResults, testFitResults_err,fSig, fBkgd, "test","test", fitRangeMin,fitRangeMax,vectorTails,kTRUE,-1);
            SetFitResults(rangeName, testName, testFitResults, testFitResults_err );



            //----------------------------------------------------------------------------------------------//
            testCounter++;
          }//Tails sets
        }//FitRanges
      }//Bckg functions
    }//Signal functions
    DrawAndRMS(rangeName,vectorTestNames);
  }//Pt ranges


}
