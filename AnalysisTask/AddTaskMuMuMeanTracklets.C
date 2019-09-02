AliAnalysisTaskMuMuMeanTracklets* AddTaskMuMuMeanTracklets(TString name = "CMUL")
{
  /// Add AliAnalysisTaskMuMuMeanTracklets to the train (Philippe Pillot)

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskMuMuMeanTracklets","AliAnalysisManager not set!");
    return NULL;
  }

  // This task run on AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("AOD")) {
    Error("AddTaskMuMuMeanTracklets", "AOD input handler needed!");
    return NULL;
  }

  // Create and configure task
  AliAnalysisTaskMuMuMeanTracklets *task = new AliAnalysisTaskMuMuMeanTracklets(name.Data());
  if (!task) {
    Error("AddTaskMuMuMeanTracklets", "MuMuMeanTracklets task cannot be created!");
    return NULL;
  }

  // Add task to analysis manager
  mgr->AddTask(task);

  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuMuMeanTracklets", "Common output file is not defined!");
    return NULL;
  }

  // Create and connect output containers
  AliAnalysisDataContainer *events = mgr->CreateContainer(Form("%s_events",name.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hNtrk = mgr->CreateContainer(Form("%s_hNtrk",name.Data()), THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hNtrkCorr = mgr->CreateContainer(Form("%s_hNtrkCorr",name.Data()), THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hNtrkCorrVsCuts = mgr->CreateContainer(Form("%s_hNtrkCorrVsCuts",name.Data()), THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *pMeanNtrkVsZvtx = mgr->CreateContainer(Form("%s_pMeanNtrkVsZvtx",name.Data()), TProfile::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *pMeanNtrkVsZvtxCorr = mgr->CreateContainer(Form("%s_pMeanNtrkVsZvtxCorr",name.Data()), TProfile::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hDimuon = mgr->CreateContainer(Form("%s_hDimuon",name.Data()), THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, events);
  mgr->ConnectOutput(task, 2, hNtrk);
  mgr->ConnectOutput(task, 3, hNtrkCorr);
  mgr->ConnectOutput(task, 4, hNtrkCorrVsCuts);
  mgr->ConnectOutput(task, 5, pMeanNtrkVsZvtx);
  mgr->ConnectOutput(task, 6, pMeanNtrkVsZvtxCorr);
  mgr->ConnectOutput(task, 7, hDimuon);

  return task;
}
