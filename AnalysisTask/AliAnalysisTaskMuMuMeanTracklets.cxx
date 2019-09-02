/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <Riostream.h>

// ROOT includes
#include "THnSparse.h"
#include "TMath.h"
#include "TRandom3.h"

// STEER includes
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliAODTracklets.h"
#include "AliCounterCollection.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"

// ANALYSIS includes
#include "AliAnalysisTaskMuMuMeanTracklets.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMuonTrackCuts.h"

ClassImp(AliAnalysisTaskMuMuMeanTracklets)

//________________________________________________________________________
AliAnalysisTaskMuMuMeanTracklets::AliAnalysisTaskMuMuMeanTracklets() :
AliAnalysisTaskSE(),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fhNtrkCorrVsCuts(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fhDimuon(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fUseBinomial(kFALSE),
fRandom(new TRandom3(0)),
fTrigger(""),
fPSTriggerMask(0),
fRejectSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE),
fReject0Tracklet(kFALSE),
fName("")
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskMuMuMeanTracklets::AliAnalysisTaskMuMuMeanTracklets(const char *name) :
AliAnalysisTaskSE(name),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fhNtrkCorrVsCuts(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fhDimuon(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fUseBinomial(kFALSE),
fRandom(new TRandom3(0)),
fTrigger(""),
fPSTriggerMask(0),
fRejectSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE),
fReject0Tracklet(kFALSE),
fName(name)
{
  /// Constructor

  // Output slot #1 writes into a AliCounterCollection
  DefineOutput(1,AliCounterCollection::Class());
  // Output slot #2 writes into a THnSparse
  DefineOutput(2,THnSparse::Class());
  // Output slot #3 writes into a THnSparse
  DefineOutput(3,THnSparse::Class());
  // Output slot #4 writes into a THnSparse
  DefineOutput(4,THnSparse::Class());
  // Output slot #5 writes into a TProfile
  DefineOutput(5,TProfile::Class());
  // Output slot #6 writes into a TProfile
  DefineOutput(6,TProfile::Class());
  // Output slot #7 writes into a TProfile
  DefineOutput(7,THnSparse::Class());

}

//________________________________________________________________________
AliAnalysisTaskMuMuMeanTracklets::~AliAnalysisTaskMuMuMeanTracklets()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fEvents;
    delete fhNtrk;
    delete fhNtrkCorr;
    delete fhNtrkCorrVsCuts;
    delete fpMeanNtrkVsZvtx;
    delete fpMeanNtrkVsZvtxCorr;
    delete fhDimuon;
  }
  delete fpMeanNtrkVsZvtxRef;
  delete fRandom;
}

//________________________________________________________________________
void AliAnalysisTaskMuMuMeanTracklets::NotifyRun()
{
  /// Set run number for cuts
  if ( fMuonTrackCuts ) fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskMuMuMeanTracklets::UserCreateOutputObjects()
{
  /// Create output objects

  printf("\nseed = %u\n\n", fRandom->GetSeed());

  // events analyzed
  fEvents = new AliCounterCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fEvents->AddRubric("event", "any/selected");
  fEvents->AddRubric("run", 10000);
  fEvents->Init();

  // prepare binning for THnSparse for N SPD tracklets studies
  // 0: N SPD tracklets
  // 1: zVtx
  // 2: run
  // 3: N MC particles
  const Int_t nDims = 4;
  Int_t nBins[nDims] = {350, 80, 300000, 350};
  Double_t xMin[nDims] = {-0.5, -10., 99999.5, -0.5};
  Double_t xMax[nDims] = {349.5, 10., 399999.5, 349.5};

  // prepare binning for THnSparse for cut efficiency studies
  // 0: Corrected (if correction profile is provided) N SPD tracklets
  // 1: N MC particles
  // 2: run
  // 3: flag: physics selection (with or without pile-up rejection)
  // 4: flag: SPD pile-up
  // 5: flag: Vtx QA
  // 6: flag: reconstructed |zVtx| < 10cm
  // 7: flag: MC |zVtx| < 10cm
  // 8: flag: NSD
  // 9: flag: INEL>0
  const Int_t nDims2 = 10;
  Int_t nBins2[nDims2] = {350, 350, 300000, 2, 2, 2, 2, 2, 2, 2};
  Double_t xMin2[nDims2] = {-0.5, -0.5, 99999.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
  Double_t xMax2[nDims2] = {349.5, 349.5, 399999.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

  // prepare binning for THnSparse for Dimuon analysis
  // 0: Run number
  // 1: Pt
  // 2: Rapidity
  // 3: Mass
  // 4: Ntracklets
  const Int_t nDimMuMu = 5;
  Int_t nBinsMuMu[nDimMuMu] = {300000, 300 , 12, 2000, 350};
  Double_t xMinMuMu[nDimMuMu] = {99999.5, -0.5, -5.5, -0.5, -0.5};
  Double_t xMaxMuMu[nDimMuMu] = {399999.5, 29.5, -1.5, 19.5, 349.5};

  // create histogram
  fhNtrk = new THnSparseT<TArrayF>(Form("%s_hNtrk",fName.Data()), "N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrk->Sumw2();
  fhNtrkCorr = new THnSparseT<TArrayF>(Form("%s_hNtrkCorr",fName.Data()), "Corrected N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrkCorr->Sumw2();
  fhNtrkCorrVsCuts = new THnSparseT<TArrayF>(Form("%s_hNtrkCorrVsCuts",fName.Data()), "Corrected N SPD tracklets versus cuts", nDims2, nBins2, xMin2, xMax2);
  fhNtrkCorrVsCuts->Sumw2();

  fpMeanNtrkVsZvtx = new TProfile(Form("%s_fpMeanNtrkVsZvtx",fName.Data()), "<Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);
  fpMeanNtrkVsZvtxCorr = new TProfile(Form("%s_fpMeanNtrkVsZvtxCorr",fName.Data()), "corrected <Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);

  //------------------------------------------------------------------------------------------------//
  //THnSparse for Dimuon Analysis
  fhDimuon = new THnSparseT<TArrayF>(Form("%s_fhDimuon",fName.Data()), "Dimuons", nDimMuMu, nBinsMuMu, xMinMuMu, xMaxMuMu);
  fhDimuon->Sumw2();
  //------------------------------------------------------------------------------------------------//

  //------------------------------------------------------------------------------------------------//
  //To make the muon selection easier, this class can be used. Make sure to call ->SetRun() function in the notifyRun() of this task.
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts","StandardMuonTrackCuts");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  //Add the cuts to be used; the nominal values are already set in the AliMuonTrackCuts task itself.
  fMuonTrackCuts->SetFilterMask (AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuPdca);
  //------------------------------------------------------------------------------------------------//

  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fhNtrkCorrVsCuts);
  PostData(5, fpMeanNtrkVsZvtx);
  PostData(6, fpMeanNtrkVsZvtxCorr);
  PostData(7, fhDimuon);

}

//________________________________________________________________________
void AliAnalysisTaskMuMuMeanTracklets::UserExec(Option_t *)
{
  /// Called for each event

  // get the input event
  AliVEvent *evt = InputEvent();

  // cast it to AOD and make sure it is actually AOD
  if (!dynamic_cast<AliAODEvent*>(evt)) return;

  // fill event counters for all events
  fEvents->Count(Form("event:any/run:%d", fCurrentRunNumber));

  // select a specific trigger (select all by default)
  Int_t triggerFlag = 1;
  if (!fTrigger.IsNull()) {
    TString trigger = evt->GetFiredTriggerClasses();
    if (fTrigger == "MB" && !(trigger.Contains("V0L") && trigger.Contains("V0R"))) triggerFlag = 0;
    else if (!trigger.Contains(fTrigger.Data())) triggerFlag = 0;
  }

  // apply physics selection
  Int_t psFlag = (fInputHandler->IsEventSelected() & fPSTriggerMask) ? 1 : 0;

  // reject events with pile-up from SPD
  Int_t spdPUFlag = evt->IsPileupFromSPDInMultBins() ? 1 : 0;

  // apply vertex QA selection
  Int_t vtxQAFlag = 0;
  const AliVVertex* vtx = evt->GetPrimaryVertexSPD();
  if (vtx && vtx->GetNContributors() > 0) {
    Double_t cov[6]={0};
    vtx->GetCovarianceMatrix(cov);
    if (TMath::Sqrt(cov[5]) <= 0.25) vtxQAFlag = 1;
  }

  // get and select on vertex position
  Double_t zVtx = (vtxQAFlag == 1) ? vtx->GetZ() : 0.;
  Int_t zVtxRangeFlag = (zVtx > -10. && zVtx < 10.) ? 1 : 0;

  // get the number of SPD tracklets (WARNING: not reliable without a good SPD vertex)
  AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(evt)->GetTracklets();
  Int_t nTrk = tracklets->GetNumberOfTracklets();
  Int_t nTrkInEtaRange(0);
  for (Int_t i = 0; i < nTrk; i++) {
    Double_t eta = -TMath::Log(TMath::Tan(tracklets->GetTheta(i)/2.));
    if ( eta < -1. || eta > 1. ) continue;
    ++nTrkInEtaRange;
  }

  // get the corrected number of SPD tracklets (WARNING: not reliable without a good SPD vertex)
  Int_t nTrkInEtaRangeCorr = fpMeanNtrkVsZvtxRef ? GetCorrectedNtrk(nTrkInEtaRange,zVtx) : GetCorrectedNtrkFromMultSel();
  Int_t nMCpartInEtaRange, zVtxMCRangeFlag, nsdFlag, inelPosFlag = 0;
  // fill histograms for cut efficiency studies
  Double_t EventInfoVsCuts[10];
  EventInfoVsCuts[0] = fpMeanNtrkVsZvtxRef ? nTrkInEtaRangeCorr : nTrkInEtaRange;
  EventInfoVsCuts[1] = nMCpartInEtaRange;
  EventInfoVsCuts[2] = fCurrentRunNumber;
  EventInfoVsCuts[3] = psFlag;
  EventInfoVsCuts[4] = spdPUFlag;
  EventInfoVsCuts[5] = vtxQAFlag;
  EventInfoVsCuts[6] = zVtxRangeFlag;
  EventInfoVsCuts[7] = zVtxMCRangeFlag;
  EventInfoVsCuts[8] = nsdFlag;
  EventInfoVsCuts[9] = inelPosFlag;
  fhNtrkCorrVsCuts->Fill(EventInfoVsCuts);

  // apply all events selections
  if (fRejectSD && nsdFlag == 0) return;
  if (triggerFlag == 0) return;
  if (fPSTriggerMask != 0 && psFlag == 0) return;
  if (fRejectPUFromSPD && spdPUFlag == 1) return;
  if (fSelectSPDVtxQA && vtxQAFlag == 0) return;
  if (zVtxRangeFlag == 0) return;
  if (fReject0Tracklet && nTrkInEtaRange == 0) return;

  // fill event counters for selected events
  fEvents->Count(Form("event:selected/run:%d", fCurrentRunNumber));

  // fill histograms for N SPD tracklets studies
  Double_t EventInfo[4];
  EventInfo[0] = nTrkInEtaRange;
  EventInfo[1] = zVtx;
  EventInfo[2] = fCurrentRunNumber;
  EventInfo[3] = nMCpartInEtaRange;
  fhNtrk->Fill(EventInfo);
  fpMeanNtrkVsZvtx->Fill(zVtx,nTrkInEtaRange);
  if (nTrkInEtaRangeCorr > 0 || (nTrkInEtaRangeCorr == 0 && !fReject0Tracklet)) {
    EventInfo[0] = nTrkInEtaRangeCorr;
    fhNtrkCorr->Fill(EventInfo);
    fpMeanNtrkVsZvtxCorr->Fill(zVtx,nTrkInEtaRangeCorr);
  }

  //------------------------------------------------------------------------------------------------//
  // Dimuon analysis
  // AliVEvent* aodesd = NULL;
  // aodesd = static_cast<AliVEvent *>(InputEvent());
  TLorentzVector lvMuon1,lvMuon2,lvDimuon;
  int nTracks = AliAnalysisMuonUtility::GetNTracks(evt);
  // Histograms to fill for Dimuon analysis
  Double_t EventInfoMuMu[5];
  EventInfoMuMu[0] = fCurrentRunNumber;
  EventInfoMuMu[4] = nTrkInEtaRangeCorr;

  for (Int_t iMuon1 = 0; iMuon1 < nTracks; iMuon1++) {
    AliVParticle *muonTrack1 = 0;
    muonTrack1 = AliAnalysisMuonUtility::GetTrack(iMuon1,evt);
    if ( !muonTrack1 ) {
      AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon1));
      continue;
    }
    if ( ! fMuonTrackCuts->IsSelected(muonTrack1) ) continue;//include cuts on pDCA, Eta, Rabs.
    TrackToLorentzVector(muonTrack1,lvMuon1);

    for (Int_t iMuon2 = iMuon1+1; iMuon2 < nTracks; iMuon2++) {
      AliVParticle *muonTrack2 = 0;
      muonTrack2 = AliAnalysisMuonUtility::GetTrack(iMuon2,evt);
      if ( !muonTrack2 ) {
        AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon2));
        continue;
      }
      if ( ! fMuonTrackCuts->IsSelected(muonTrack2) ) continue;//include cuts on pDCA, Eta, Rabs.
      TrackToLorentzVector(muonTrack2,lvMuon2);

      lvDimuon = lvMuon1+lvMuon2;
      if(muonTrack1->Charge() == muonTrack2->Charge()) continue;
      // Fill THnSparse with corresponding observables
      EventInfoMuMu[1] = lvDimuon.Pt();
      EventInfoMuMu[2] = lvDimuon.Rapidity();
      EventInfoMuMu[3] = lvDimuon.M();
      fhDimuon->Fill(EventInfoMuMu);
    }
  }
  //------------------------------------------------------------------------------------------------//

  // Post data
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fhNtrkCorrVsCuts);
  PostData(5, fpMeanNtrkVsZvtx);
  PostData(6, fpMeanNtrkVsZvtxCorr);
  PostData(7, fhDimuon);

}

//________________________________________________________________________
Int_t AliAnalysisTaskMuMuMeanTracklets::GetCorrectedNtrk(Int_t nTrkInEtaRange, Double_t zVtx)
{
  /// return the number of SPD tracklet corrected for the Zvtx dependence of <Ntrk>.
  /// return -999 in case of invalid correction

  Double_t meanNtrk = fpMeanNtrkVsZvtxRef->GetBinContent(fpMeanNtrkVsZvtxRef->FindBin(zVtx));
  if (meanNtrk < 1.e-6) return -999;

  if (fUseBinomial && fMeanNtrkRef <= meanNtrk) {

    return fRandom->Binomial(nTrkInEtaRange, fMeanNtrkRef/meanNtrk);

  } else {

    Double_t dN = nTrkInEtaRange*fMeanNtrkRef/meanNtrk - nTrkInEtaRange;
    Int_t sign = (dN > 0.) ? 1 : -1;
    /*
     Int_t nTrkInEtaRangeCorr = -1;
     do {
     nTrkInEtaRangeCorr = nTrkInEtaRange + sign*fRandom->Poisson(TMath::Abs(dN));
     } while (nTrkInEtaRangeCorr < 0);

     return nTrkInEtaRangeCorr;
     */
    return TMath::Max(nTrkInEtaRange + sign*fRandom->Poisson(TMath::Abs(dN)), 0);

  }

}

//________________________________________________________________________
Int_t AliAnalysisTaskMuMuMeanTracklets::GetCorrectedNtrkFromMultSel()
{
  /// return the corrected number of SPD tracklet from the multiplicity framework
  /// return -999 if not set

  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if (!multSelection) return -999;

  return static_cast<Int_t>(multSelection->GetEstimator("SPDTracklets")->GetValue());

}

//------------------------------------------------------------------------------------------------//
//Dummy function to get a TLorentzVector from a track
void AliAnalysisTaskMuMuMeanTracklets::TrackToLorentzVector(AliVParticle *muonTrack, TLorentzVector &lvMuon){
  Double_t muonMass2 = AliAnalysisMuonUtility::MuonMass2();
  Double_t energy = TMath::Sqrt(muonTrack->P()*muonTrack->P() + muonMass2);
  lvMuon.SetPxPyPzE(muonTrack->Px(),muonTrack->Py(),muonTrack->Pz(),energy);
}
