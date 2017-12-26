#ifndef GeRecorderMessenger_hh
#define GeRecorderMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"


class GeRecorder;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class GeRecorderMessenger:public G4UImessenger
{
public:
    GeRecorderMessenger(GeRecorder *);
    ~GeRecorderMessenger();

    void SetNewValue(G4UIcommand *, G4String);

private:
    GeRecorder             *Recorder;

    G4UIdirectory             *RecMsg;
    G4UIcmdWithAString 		  *FileName;
};

#endif
