
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "Randomize.hh"


//User Scripts
//------------------------------------------
#include "GeDetectorConstruction.hh"
#include "GePrimaryGeneratorAction.hh"
#include "GeRunAction.hh"
#include "GeEventAction.hh"
#include "GeSteppingAction.hh"
#include "GeRecorder.hh"
//------------------------------------------

//Physic List
#include "PhysicsList.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char **argv)
{
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif 
 
       
  //------------------------
  runManager->SetUserInitialization(new GeDetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);

  G4VUserPrimaryGeneratorAction *gen_action = new GePrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  GeRecorder *recorder = new GeRecorder();
  
  G4UserRunAction *run_action = new GeRunAction(recorder);
  runManager->SetUserAction(run_action);
  G4UserEventAction *event_action = new GeEventAction(recorder);
  runManager->SetUserAction(event_action);
  G4UserSteppingAction *stepping_action = new GeSteppingAction(recorder);
  runManager->SetUserAction(stepping_action);
  //------------------------
  
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( ! ui ) { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}
