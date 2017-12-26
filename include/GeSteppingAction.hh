
#ifndef GeSteppingAction_h
#define GeSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class GeRecorder;

class GeSteppingAction : public G4UserSteppingAction
{
public:
    GeSteppingAction(GeRecorder *);
    ~GeSteppingAction();

    void UserSteppingAction(const G4Step *);

private:
    GeRecorder *recorder;
};

#endif
