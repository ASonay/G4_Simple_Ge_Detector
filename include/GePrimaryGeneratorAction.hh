
#ifndef GePrimaryGeneratorAction_h
#define GePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"


class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;

class GePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GePrimaryGeneratorAction();
  virtual ~GePrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event *);
  
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

private:
  G4GeneralParticleSource *generator;
  G4ParticleGun *fParticleGun;
};

#endif
