
#include "GeEventAction.hh"

#include "GeRecorder.hh"
#include "G4Event.hh"
#include "Randomize.hh"


GeEventAction::GeEventAction(GeRecorder *rec)
    : recorder(rec)
{
}


GeEventAction::~GeEventAction()
{
}


void GeEventAction::BeginOfEventAction(const G4Event *evt)
{
    recorder->BeginOfRecordingEvent(evt);
}


void GeEventAction::EndOfEventAction(const G4Event *evt)
{
    recorder->EndOfRecordingEvent(evt);
}
