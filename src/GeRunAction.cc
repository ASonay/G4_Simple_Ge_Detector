#include "GeRunAction.hh"

#include "GeRecorder.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"


GeRunAction::GeRunAction(GeRecorder *rec)
  : recorder(rec)
{
}


GeRunAction::~GeRunAction()
{
}


void GeRunAction::BeginOfRunAction(const G4Run *aRun)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  recorder->BeginOfRecordingRun(aRun);
}


void GeRunAction::EndOfRunAction(const G4Run *aRun)
{
  recorder->EndOfRecordingRun(aRun);
}

