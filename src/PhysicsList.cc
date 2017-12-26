#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "GammaPhysics.hh"
#include "NeutronHPphysics.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{
  G4int verb = 1;
  SetVerboseLevel(verb);
  
  //add new units for radioActive decays
  //
  new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  // 
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);
          
  
  // Decay
  RegisterPhysics(new G4DecayPhysics(verb));

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics(verb));
  
  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(verb) );
  
  // EM physics
  RegisterPhysics( new G4EmStandardPhysics(verb) );

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(verb) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(verb) );
  
  // Gamma-Nuclear Physics
  RegisterPhysics( new GammaPhysics("gamma") );

  // Neutron HP Physics
  RegisterPhysics( new NeutronHPphysics("neutronHP"));  
}


PhysicsList::~PhysicsList()
{ }


void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


void PhysicsList::SetCuts()
{


  if (verboseLevel >1)
    G4cout << "DMXPhysicsList::SetCuts:";
  
  if (verboseLevel>0){
    G4cout << "DMXPhysicsList::SetCuts:";
    G4cout << "CutLength : " 
	   << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  SetCutValue(1*nm, "neutron");
  SetCutValue(100*nm, "alpha");
  SetCutValue(100*nm, "proton");
  SetCutValue(100*nm, "e-");
  SetCutValue(100*nm, "e+");
  SetCutValue(100*nm, "gamma");      
}

