#include "GeRecorder.hh"
#include "GeRecorderMessenger.hh"


//#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

GeRecorderMessenger::GeRecorderMessenger(GeRecorder *rec)
  : Recorder(rec)
{

  RecMsg = new G4UIdirectory("/Ge/RecMsg/");
  RecMsg->SetGuidance("UI commands specific to GeRecorder");

  FileName = new G4UIcmdWithAString("/Ge/RecMsg/FileName", this);
  FileName->SetGuidance("Specify the output ROOT file name.");
  FileName->SetParameterName("fName", false);
  FileName->AvailableForStates(G4State_PreInit, G4State_Idle);
}


GeRecorderMessenger::~GeRecorderMessenger()
{
  delete RecMsg;
  delete FileName;
}

void GeRecorderMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{	
  if(command == FileName)
    Recorder->SetFileName(newValue);
}
