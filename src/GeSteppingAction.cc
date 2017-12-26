
#include "GeSteppingAction.hh"

#include "GeRecorder.hh"
#include "G4Step.hh"


GeSteppingAction::GeSteppingAction(GeRecorder *rec)
    : recorder(rec)
{
}


GeSteppingAction::~GeSteppingAction()
{
}


void GeSteppingAction::UserSteppingAction(const G4Step *aStep)
{
    recorder->BeginOfRecordingStep(aStep);
    recorder->EndOfRecordingStep(aStep);
}
