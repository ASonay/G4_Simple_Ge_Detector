#ifndef NeutronHPMessenger_h
#define NeutronHPMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NeutronHPphysics;
class G4UIdirectory;
class G4UIcmdWithABool;


class NeutronHPMessenger: public G4UImessenger
{
  public:
    NeutronHPMessenger(NeutronHPphysics*);
   ~NeutronHPMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:    
    NeutronHPphysics*  fNeutronPhysics;
    
    G4UIdirectory*     fPhysDir;      
    G4UIcmdWithABool*  fThermalCmd;
};


#endif
