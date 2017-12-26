#ifndef GeRunAction_h
#define GeRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class GeRecorder;


class GeRunAction : public G4UserRunAction
{
public:
    GeRunAction(GeRecorder *);
    ~GeRunAction();

    void BeginOfRunAction(const G4Run *);
    void   EndOfRunAction(const G4Run *);

private:
    GeRecorder *recorder;
};


#endif
