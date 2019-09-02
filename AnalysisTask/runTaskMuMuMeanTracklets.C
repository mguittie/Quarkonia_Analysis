/*
 *  runTaskDimuon.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 10/04/17.
 *  Copyright 2017 SUBATECH. All rights reserved.
 *
 */

 #if !defined(__runGenTuner__)
 #define __runGenTuner__

 #include "/home/alidock/Quarkonia_Analysis/analysis/Facilities/runTaskFacilities.C"

//______________________________________________________________________________
void runTaskMuMuMeanTracklets(TString smode = "local", TString inputFileName = "runListLHC16h.txt",
                          Bool_t MuMu = kFALSE, Bool_t applyPS = kTRUE, Bool_t applyPileupCuts = kTRUE,
                          Bool_t isMC = kTRUE, TString refInput = "")
{
  /// Run the baseline task to test the framework

  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20190201_ROOT6-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskMuMuMeanTracklets";
  TString extraPkgs="";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("/home/alidock/Quarkonia_Analysis/AnalysisTask"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runTaskMuMuMeanTracklets.C"));
  fileList.Add(new TObjString("AddTaskMuMuMeanTracklets.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuMuMeanTracklets.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuMuMeanTracklets.h"));

  // --- grid specific setup ---
  TString dataDir = "/alice/data/2016/LHC16k";
  TString dataPattern = "/pass1/AOD/*AliAOD.Muons.root";
  TString runFormat = "%09d";
  TString outDir = "OutputQuarkonia/LHC16kMuMu";
  TString analysisMacroName = "MuMuMeanTracklets";

  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 40;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 1;

  // --- saf3 specific setup ---
  Bool_t splitDataset = kFALSE;

  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  if (!refInput.IsNull()) {
    CopyInputFileLocally(refInput.Data(), "ReferenceResults.root", 'a');
    fileList.Add(new TObjString("ReferenceResults.root"));
  }

  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {

    if (!RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset)) return;

  } else {

    gSystem->Exec(TString::Format("cp %s __runTask__.C", __FILE__));
    gROOT->LoadMacro("__runTask__.C");
    gSystem->Exec("rm __runTask__.C");
    gROOT->ProcessLineSync(TString::Format("CreateAnalysisTrain(%d, %d, %d, \"%s\")",applyPS, applyPileupCuts, isMC, refInput.Data()));

    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);

    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);

  }

}

#else

Bool_t SetMeanNtrkVsZvtxRef(TObject* meanTracklets);

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t applyPS, Bool_t applyPileupCuts, Bool_t isMC, TString refInput)
{
  /// create the analysis train and configure it

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuMuMeanTrackletsAnalysis");
  //mgr->SetDebugLevel(3);

  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);

  if (applyPS) {
    // event selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = reinterpret_cast<AliPhysicsSelectionTask*>(gROOT->ProcessLineSync(TString::Format("AddTaskPhysicsSelection(%d, %d)", isMC,applyPileupCuts)));
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return;
    }
  }
  /*
  // multiplicity/centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *mult = AddTaskMultSelection(kFALSE);
  if (applyPS) mult->SelectCollisionCandidates(AliVEvent::kINT7inMUON | AliVEvent::kINT7);
  //  if (applyPS) mult->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt7);
  */
  // MeanTracklets task
  gROOT->LoadMacro("AddTaskMuMuMeanTracklets.C");
  AliAnalysisTaskMuMuMeanTracklets* MuMumeanTrackletsCINT = reinterpret_cast<AliAnalysisTaskMuMuMeanTracklets*>(gROOT->ProcessLineSync("AddTaskMuMuMeanTracklets(\"CINT\")"));
  if(!MuMumeanTrackletsCINT) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuMuMeanTracklets not created!");
    return;
  }
  // if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kINT7inMUON | AliVEvent::kINT7);
  if (applyPS) MuMumeanTrackletsCINT->ApplyPhysicsSelection(AliVEvent::kINT7inMUON);
  // if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kMuonUnlikeLowPt7);
  //  meanTracklets->SelectTrigger("MB");
  //  meanTracklets->SelectTrigger("CINT7-B-NOPF-");
  //  meanTracklets->SelectTrigger("CMUL7-B-NOPF-MUFAST");
   // if (isMC) meanTracklets->RejectSD();
  //  meanTracklets->RejectPUFromSPD();
  //  meanTracklets->DisableSPDVtxQA();
  //  meanTracklets->Reject0Tracklet();
  if (!refInput.IsNull() && !SetMeanNtrkVsZvtxRef(MuMumeanTrackletsCINT)) return;

  if (!refInput.IsNull()) {
    // gROOT->LoadMacro("AddTaskMuMuMeanTracklets.C");
    AliAnalysisTaskMuMuMeanTracklets* MuMumeanTrackletsCMUL = reinterpret_cast<AliAnalysisTaskMuMuMeanTracklets*>(gROOT->ProcessLineSync("AddTaskMuMuMeanTracklets()"));
    if(!MuMumeanTrackletsCMUL) {
      Error("CreateAnalysisTrain","AliAnalysisTaskMuMuMeanTracklets not created!");
      return;
    }
    // if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kINT7inMUON | AliVEvent::kINT7);
    // if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kINT7inMUON);
    if (applyPS) MuMumeanTrackletsCMUL->ApplyPhysicsSelection(AliVEvent::kMuonUnlikeLowPt7);

    if (!refInput.IsNull() && !SetMeanNtrkVsZvtxRef(MuMumeanTrackletsCMUL)) return;
  }
    //  meanTracklets->UseBinomial();
}

//______________________________________________________________________________
Bool_t SetMeanNtrkVsZvtxRef(TObject* meanTracklets)
{
  /// set the reference <Ntrk> vs Z profile

  TFile *file = new TFile("ReferenceResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ReferenceResults.root\n");
    return kFALSE;
  }

  TProfile *p = static_cast<TProfile*>(file->FindObjectAny("CINT_fpMeanNtrkVsZvtx"));
  if (!p) {
    Error("CreateAnalysisTrain","cannot find the reference <Ntrk> vs Z profile!");
    return kFALSE;
  }

  AliAnalysisTaskMuMuMeanTracklets *mT = static_cast<AliAnalysisTaskMuMuMeanTracklets*>(meanTracklets);
//  mT->SetMeanNtrkVsZvtxRef(*p);
  mT->SetMeanNtrkVsZvtxRef(*p, p->GetMaximum());

  file->Close();
  delete file;

  return kTRUE;

}

#endif
