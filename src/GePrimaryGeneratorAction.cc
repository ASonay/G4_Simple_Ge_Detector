#include "GePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


GePrimaryGeneratorAction::GePrimaryGeneratorAction()
{ generator = new G4GeneralParticleSource; }

GePrimaryGeneratorAction::~GePrimaryGeneratorAction()
{ delete generator; }


//call GeneratePrimaryVertex from GeneratePrimaries method
void GePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ generator->GeneratePrimaryVertex(anEvent); }
