#ifndef AliAnalysisTaskMuMuMeanTracklets_H
#define AliAnalysisTaskMuMuMeanTracklets_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muon
/// \class AliAnalysisTaskMuMuMeanTracklets
/// \brief task to study the mean number of SPD tracklets and correct for the z dependence
//Author: Philippe Pillot - SUBATECH Nantes

#include "TProfile.h"
#include "AliAnalysisTaskSE.h"

class THnSparse;
class TRandom3;
class AliCounterCollection;
class AliMuonTrackCuts;


class AliAnalysisTaskMuMuMeanTracklets : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskMuMuMeanTracklets();
  AliAnalysisTaskMuMuMeanTracklets(const char *name);
  virtual ~AliAnalysisTaskMuMuMeanTracklets();

  virtual void   NotifyRun();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *){}
  void           TrackToLorentzVector(AliVParticle *muonTrack, TLorentzVector &lvMuon);


  void SetMeanNtrkVsZvtxRef(TProfile &meanNtrkVsZvtxRef, Double_t meanNtrkRef = -1.);

  /// use a binomial instead of poissonian distribution to correct Ntrk when fMeanNtrkRef < <Ntrk>
  void UseBinomial(Bool_t flag = kTRUE) {fUseBinomial = flag;}

  /// select a specific trigger class
  void SelectTrigger(TString trigger) {fTrigger = trigger;}

  /// apply physics selection with the given trigger mask
  void ApplyPhysicsSelection(UInt_t triggerMask) {fPSTriggerMask = triggerMask;}

  /// reject single-diffractive MC events
  void RejectSD() {fRejectSD = kTRUE;}

  /// reject pile-up from SPD
  void RejectPUFromSPD() {fRejectPUFromSPD = kTRUE;}

  /// disable SPD vertex QA selection
  void DisableSPDVtxQA() {fSelectSPDVtxQA = kFALSE;}

  /// reject events with 0 (corrected) tracklets in -1 < eta < 1
  void Reject0Tracklet() {fReject0Tracklet = kTRUE;}

private:

  /// Not implemented
  AliAnalysisTaskMuMuMeanTracklets(const AliAnalysisTaskMuMuMeanTracklets& rhs);
  /// Not implemented
  AliAnalysisTaskMuMuMeanTracklets& operator = (const AliAnalysisTaskMuMuMeanTracklets& rhs);

  Int_t GetCorrectedNtrk(Int_t NtrkInEtaRange, Double_t zVtx);
  Int_t GetCorrectedNtrkFromMultSel();

private:

  AliMuonTrackCuts *fMuonTrackCuts;//!

  AliCounterCollection* fEvents;  //!< number of analyzed events
  THnSparse *fhNtrk;              //!< output histogram for number of SPD tracklets
  THnSparse *fhNtrkCorr;          //!< output histogram for corrected number of SPD tracklets
  THnSparse *fhNtrkCorrVsCuts;    //!< output histogram for corrected number of SPD tracklets versus cuts
  TProfile *fpMeanNtrkVsZvtx;     //!< <Ntrk> vs Z profile before correction
  TProfile *fpMeanNtrkVsZvtxCorr; //!< <Ntrk> vs Z profile after correction

  THnSparse *fhDimuon;              //!< output histogram for dimuons

  TProfile *fpMeanNtrkVsZvtxRef;  /// <Ntrk> vs Z profile used to correct Ntrk
  Double_t fMeanNtrkRef;          /// <Ntrk> value used as a reference
  Bool_t   fUseBinomial;          /// use a binomial distribution to correct Ntrk when fMeanNtrkRef < <Ntrk>
  TRandom3 *fRandom;              //!< random number generator

  TString fTrigger;               /// select a specific trigger class
  UInt_t fPSTriggerMask;          /// apply physics selection with the given trigger mask
  Bool_t fRejectSD;               /// reject single-diffractive MC events
  Bool_t fRejectPUFromSPD;        /// reject pile-up from SPD
  Bool_t fSelectSPDVtxQA;         /// select events with a good SPD vertex
  Bool_t fReject0Tracklet;        /// reject events with 0 (corrected) tracklets in -1 < eta < 1
  TString fName;                  /// Trigger name


  ClassDef(AliAnalysisTaskMuMuMeanTracklets, 3);
};

//________________________________________________________________________
inline void AliAnalysisTaskMuMuMeanTracklets::SetMeanNtrkVsZvtxRef(TProfile &meanNtrkVsZvtxRef, Double_t meanNtrkRef)
{
  /// set the <Ntrk> vs Z profile used to correct Ntrk
  /// set the reference <Ntrk> to the given value or to the minimum if not provided
  delete fpMeanNtrkVsZvtxRef;
  fpMeanNtrkVsZvtxRef = new TProfile(meanNtrkVsZvtxRef);
  fpMeanNtrkVsZvtxRef->SetDirectory(0);
  fMeanNtrkRef = (meanNtrkRef > 0.) ? meanNtrkRef : fpMeanNtrkVsZvtxRef->GetMinimum();
  if (fMeanNtrkRef < 1.e-6) AliFatal(Form("Invalid <Ntrk> reference value: %f",fMeanNtrkRef));
}

#endif
