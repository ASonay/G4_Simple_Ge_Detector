#ifndef GeRecorder_h
#define GeRecorder_h 1

#include "GeRecorderMessenger.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"

#include "Rtypes.h"  // Basic types used by ROOT.
#include "TROOT.h"   // Entry point to the system.
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNtuple.h"


class G4Run;
class G4Event;
class G4Step;

class GeTrajectory;


class GeRecorder
{
public:
  GeRecorder();
  ~GeRecorder();

  void BeginOfRecordingRun(const G4Run *);
  void EndOfRecordingRun(const G4Run *); 
    
  void BeginOfRecordingEvent(const G4Event *);
  void EndOfRecordingEvent(const G4Event *);

  void BeginOfRecordingStep(const G4Step *);
  void EndOfRecordingStep(const G4Step *);

  void ResetAllMembers();
  void SetFileName(const G4String &);

private:
  TFile *rootFile;
  TTree *tr;
  
  //---
  G4double ratio_R;
  G4int evtTot,evtNb;
  G4int mul_evt;
  //---
  
  Char_t parName[100];
  Char_t procName[100];
  Char_t volName[100];
  G4int hit_Ge,hit_Ge_b;
  G4int ge_trig;
  G4int evt_id;
  G4int idx;
  G4int ntc,ntc_p;
  G4int track_par,track_elastic,track_inelastic;
  G4int index_ctrl;
  G4int countPar;
  G4int countParDec[500];
  G4int getTrackDec[500];
  G4int nindex_tc[500];
  G4int trackID,parentID;
  G4int pretrackID0;
  G4int elastic_onoff,inelastic_onoff,rad_onoff;
  G4double total_time[500];
  G4double edep;
  G4double x, y, z;
  G4double x1, y1, z1, r1, h1;
  G4double dr, dh;
  G4double epson,k,g,f;
  G4double hrec_tot,hrec_alpha,hrec_tot_lindh,hrec_inelastic_tot,hrec_elastic_tot;
  G4double edep_tot;
  G4double Etot,Egun,Er,Ea,Eqr,Eqr_lindh,Eqr_elastic,Eqr_inelastic,Eqr_elastic_lindh,Eqr_inelastic_lindh,Er_elastic,Er_inelastic;
  G4double Etot_tc[100];
  G4double time_inf[500];
  G4String rootFileName;
  G4String particleName;
  G4String parNameDec[500];
  G4String volumeName;
  G4String processName;
  G4String nul;
  G4double locT, globT, propT;
  G4double tdc,livT,time_int,time_fin,pre_time;
  GeRecorderMessenger *recMessenger;
};


#endif
