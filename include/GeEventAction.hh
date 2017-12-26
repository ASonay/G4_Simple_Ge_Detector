
#ifndef GeEventAction_h
#define GeEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class GeRecorder;


class GeEventAction : public G4UserEventAction
{
public:
    GeEventAction(GeRecorder *);
    ~GeEventAction();

    void BeginOfEventAction(const G4Event *);
    void EndOfEventAction(const G4Event *);

private:
    GeRecorder *recorder;
};


#endif
