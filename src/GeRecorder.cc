//---------------------------------------------------------------------//
//---------------------------------------------------------------------//

#include "GeRecorder.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"

#include "G4RunManager.hh"
#include "G4VProcess.hh"
#include "G4TrajectoryContainer.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Track.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"



#define SIGNAL_COUT

//---------------------------------------------------------------------//

GeRecorder::GeRecorder()
{
  ge_trig=0;
  rootFileName = "knownFile.root";
  recMessenger = new  GeRecorderMessenger(this);  
}

//---------------------------------------------------------------------//

GeRecorder::~GeRecorder()
{
  delete rootFile;
  delete recMessenger;
}

//---------------------------------------------------------------------//


double qin_trim(double Tm)
{
  const float qin0 = 0.23236; 
  const float qin1 = 0.053429;
  const float qin2 = 0.0091422;
  const float qin3 = 0.0037246; 
  const float qin4 = 0.0007226; 
  double result;
  double logT;

  logT = log10(1000*Tm);

  result = qin0 + qin1*logT + qin2*pow(logT,2.0) + qin3*pow(logT,3.0) + qin4*pow(logT,4.0);

  return result;
}


void GeRecorder::BeginOfRecordingRun(const G4Run *aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << " Last Event: " << aRun->GetNumberOfEventToBeProcessed() <<G4endl;
  
  G4cout << "Create ROOT file to record data .........." << G4endl;

  rootFile = new TFile(rootFileName, "recreate");

  if (rootFile->IsZombie()) {
    G4cout << "Error opening file" << G4endl;
    exit(-1);
  }
  else {
    tr = new TTree("tr", "Ge Detector Simulation");

    // ID
    tr->Branch("evt_id", &evt_id, "evt_id/D");
    
    tr->Branch("locT", &locT, "locT/D");
    tr->Branch("globT", &globT, "globT/D");
    tr->Branch("propT", &propT, "propT/D");

    tr->Branch("parName", &parName, "parName[20]/C");
    tr->Branch("volName", &volName, "volName[20]/C");
    tr->Branch("procName", &procName, "procName[20]/C");
	
    // position
    tr->Branch("x", &x, "x/D");
    tr->Branch("y", &y, "y/D");
    tr->Branch("z", &z, "z/D");
    
    tr->Branch("deltah", &dh, "deltah/D");
    tr->Branch("deltar", &dr, "deltar/D");
	
    // energy deposition
    tr->Branch("Egun", &Egun, "Egun/D");
    tr->Branch("Etot", &Etot, "Etot/D");
    tr->Branch("Ea", &Ea, "Ea/D");
    tr->Branch("Er", &Er, "Er/D");
    tr->Branch("Er_elastic", &Er_elastic, "Er_elastic/D");
    tr->Branch("Er_inelastic", &Er_inelastic, "Er_inelastic/D");
    tr->Branch("Eqr", &Eqr, "Eqr/D");
    tr->Branch("Eqr_lindh", &Eqr_lindh, "Eqr_lindh/D");
    tr->Branch("Eqr_elastic", &Eqr_elastic, "Eqr_elastic/D");
    tr->Branch("Eqr_inelastic", &Eqr_inelastic, "Eqr_inelastic/D");
    tr->Branch("Eqr_elastic_lindh", &Eqr_elastic_lindh, "Eqr_elastic_lindh/D");
    tr->Branch("Eqr_inelastic_lindh", &Eqr_inelastic_lindh, "Eqr_inelastic_lindh/D");
    
    tr->Branch("Etot_tc", Etot_tc, "Etot_tc[100]/D");
  
    tr->Branch("tdc", &tdc, "tdc/D");
    tr->Branch("livT", &livT, "livT/D");
    tr->Branch("bulk", &hit_Ge_b, "bulk/I");
  }

  mul_evt = aRun->GetNumberOfEventToBeProcessed();
}

//---------------------------------------------------------------------//
void GeRecorder::EndOfRecordingRun(const G4Run *aRun)
{
  
  evtNb=aRun->GetNumberOfEvent();
  ratio_R=100;
  
  G4cout <<"%"<< ratio_R <<" has Done." <<" Total Triger: "<<ge_trig<<G4endl;
   
  
  G4cout << "" << G4endl;
  G4cout << "====================================================================" << G4endl;
  
  rootFile->Write();
  rootFile->Close();

  G4cout << ".................... Close ROOT file" << G4endl;

  G4cout << std::setw(17) << "  Name:"
	 << std::setw(14) << "  NoE:"
	 << std::setw(17) << "  Total Time:" <<G4endl;
  for (int i=0;i<countPar;i++)
    {      
      G4cout << "  " << std::setw(17) <<parNameDec[i]
	     << "  " << std::setw(9) <<countParDec[i]
	     << "  " << std::setw(14) <<total_time[i]<<G4endl;  
    }
}

//---------------------------------------------------------------------//

void GeRecorder::BeginOfRecordingEvent(const G4Event *anEvent)
{
  ResetAllMembers();

  evt_id = anEvent->GetEventID();
  G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex();
  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  Egun = primaryParticle->GetKineticEnergy();
}

//---------------------------------------------------------------------//

void GeRecorder::EndOfRecordingEvent(const G4Event *anEvent)
{
  // information==========================================================//
  evtTot=mul_evt;
  evtNb = anEvent->GetEventID();
  ratio_R = 100.*((float)evtNb/(float)evtTot);

  if(evtNb%1000==0)
    G4cout <<"%"<< ratio_R <<" has Done." <<" Total Triger: "<<ge_trig<<G4endl;
  //======================================================================//
  nul = "no hit Ge";
  evt_id = (Int_t)evt_id;

  // Chang units
  x = x / cm;
  y = y / cm;
  z = z / cm;

  edep_tot = edep_tot / MeV;
  hrec_tot = hrec_tot / MeV;
  hrec_tot_lindh = hrec_tot / MeV;
  hrec_elastic_tot = hrec_elastic_tot / MeV;
  hrec_inelastic_tot = hrec_inelastic_tot / MeV;

  Er = hrec_tot / keV;
  Er_elastic = hrec_elastic_tot / keV;
  Er_inelastic = hrec_inelastic_tot / keV;
  Ea = hrec_alpha / keV;

  epson = 11.5*hrec_tot_lindh*pow(32,-7/3);
  k = 0.133*pow(32,2/3)*pow(72.64,1/2);
  g = 3*pow(epson,0.2)+0.7*pow(epson,0.6)+epson;
  f = k*g/(1+k*g);
  Eqr_lindh = hrec_tot_lindh*f;
  
  epson = 11.5*hrec_elastic_tot*pow(32,-7/3);
  k = 0.133*pow(32,2/3)*pow(72.64,1/2);
  g = 3*pow(epson,0.2)+0.7*pow(epson,0.6)+epson;
  f = k*g/(1+k*g);
  Eqr_elastic_lindh = hrec_elastic_tot*f;
  
  epson = 11.5*hrec_inelastic_tot*pow(32,-7/3);
  k = 0.133*pow(32,2/3)*pow(72.64,1/2);
  g = 3*pow(epson,0.2)+0.7*pow(epson,0.6)+epson;
  f = k*g/(1+k*g);
  Eqr_inelastic_lindh = hrec_inelastic_tot*f;

 if (hrec_tot > 0)
    hrec_tot = hrec_tot*qin_trim(hrec_tot);
  if (hrec_elastic_tot > 0)
    hrec_elastic_tot = hrec_elastic_tot*qin_trim(hrec_elastic_tot);
  if (hrec_inelastic_tot > 0)
    hrec_inelastic_tot = hrec_inelastic_tot*qin_trim(hrec_inelastic_tot);

  
  Etot = edep_tot / keV;
  
  Eqr = hrec_tot / keV;
  Eqr_elastic = hrec_elastic_tot / keV;
  Eqr_inelastic = hrec_inelastic_tot / keV;

  Eqr_lindh = Eqr_lindh / keV;  
  Eqr_elastic_lindh = Eqr_elastic_lindh / keV;
  Eqr_inelastic_lindh = Eqr_inelastic_lindh / keV;

  dr = 30-r1 / mm;
  dh = 71-h1 / mm;
  

  if (hit_Ge==0)
    {
      strcpy(parName, nul.c_str());
      strcpy(volName, nul.c_str());
      strcpy(procName, nul.c_str());
      x = 0;
      y = 0;
      z = 0;
    }

  //if(Etot>0||Eqr>0){
    ge_trig++;
    tr->Fill();
    //}
}

//---------------------------------------------------------------------//

void GeRecorder::BeginOfRecordingStep(const G4Step *aStep)
{
  edep = aStep->GetTotalEnergyDeposit();

  G4Track *fTrack = aStep->GetTrack();   
  trackID = fTrack->GetTrackID();
  parentID = fTrack->GetParentID();
  //G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  if (postStepPoint->GetProcessDefinedStep() != 0)
    processName = postStepPoint->GetProcessDefinedStep()->GetProcessName();

  particleName = fTrack->GetDefinition()->GetParticleName();
  volumeName = fTrack->GetVolume()->GetName();


  if (hit_Ge !=1){
    locT = fTrack->GetLocalTime() / ns;
    globT = fTrack->GetGlobalTime() / ns;
    propT = fTrack->GetProperTime() / ns;
    strcpy(parName, particleName.c_str());
    strcpy(volName, volumeName.c_str());
    strcpy(procName, processName.c_str());
    x = fTrack->GetPosition().x();
    y = fTrack->GetPosition().y();
    z = fTrack->GetPosition().z();
    tdc=fTrack->GetGlobalTime() / ns;
  }
  
  if ((volumeName == "GeDet"||volumeName == "GeCover")){
    hit_Ge=1;
    x1 = fTrack->GetPosition().x();
    y1 = fTrack->GetPosition().y();
    z1 = fTrack->GetPosition().z();
    h1 = z1;
    r1 = sqrt(x1*x1+y1*y1);

  }
  
  if (volumeName == "GeDet")
    hit_Ge_b=1;

  if ((volumeName == "GeDet" || volumeName == "GeCover") && processName == "neutronInelastic"){
    track_inelastic  = fTrack->GetTrackID();
    elastic_onoff = 0;
    inelastic_onoff = 1;
  }
  
  if ((volumeName == "GeDet" || volumeName == "GeCover") && processName == "hadElastic"){
    track_elastic  = fTrack->GetTrackID();
    elastic_onoff = 1;
    inelastic_onoff = 0;
  }


  if ((volumeName == "GeDet" || volumeName == "GeCover") && edep > 0  && (processName == "hIoni" || processName == "ionIoni" || processName == "hadElastic" || processName == "neutronInelastic" || processName == "nCapture" || processName == "nFission") && particleName != "alpha" && fTrack->GetGlobalTime()<100000 )
    hrec_tot += edep;
  
  if ((volumeName == "GeDet" || volumeName == "GeCover") && edep > 0 && particleName == "alpha" && fTrack->GetGlobalTime()<1000 )
    hrec_alpha += edep;
  
  if ((volumeName == "GeDet" || volumeName == "GeCover") && edep > 0  && processName == "ionIoni" && particleName != "alpha" && fTrack->GetGlobalTime()<1000 && elastic_onoff == 1)
    hrec_elastic_tot += edep;
  
  if ((volumeName == "GeDet" || volumeName == "GeCover") && edep > 0  && processName == "ionIoni" && particleName != "alpha" && fTrack->GetGlobalTime()<1000 && inelastic_onoff == 1)
    hrec_inelastic_tot += edep;
	
  if ((volumeName == "GeDet" || volumeName == "GeCover") && edep > 0  && !(processName == "hIoni" || processName == "ionIoni" || processName == "hadElastic" || processName == "neutronInelastic" || processName == "nCapture" || processName == "nFission") && fTrack->GetGlobalTime()<100000 )
    edep_tot += edep;

  //GATHER DECAY CHANNELS INTO DIFFERENT PARTS OF ENERGY#####################################
  if ((volumeName == "GeDet" || volumeName == "GeCover") && fTrack->GetGlobalTime()>100000 && (processName == "RadioactiveDecay")){
    time_inf[ntc] = fTrack->GetGlobalTime();
    nindex_tc[ntc] = ntc;
    ntc++;
  }

  if ((volumeName == "GeDet" || volumeName == "GeCover") && fTrack->GetGlobalTime()>100000){

    for (int i=0;i<ntc+1;i++){
      if (abs(fTrack->GetGlobalTime()-time_inf[i])<100000){
	ntc_p = nindex_tc[i];
	break;
      }
    }
    
    Etot_tc[ntc_p] += edep;
  }
  //GATHER DECAY CHANNELS INTO DIFFERENT PARTS OF ENERGY#####################################

  //ALL PROCESSES############################################################################
  index_ctrl=0;
  if ((volumeName == "GeDet"||volumeName == "GeCover")){
    for (int i=0;i<500;i++){
      if (parNameDec[i]==particleName){
	getTrackDec[i]=fTrack->GetTrackID();
	if (processName == "RadioactiveDecay")
	  total_time[i]=fTrack->GetLocalTime()*pow(10,-9)/86400.;
	if (trackID!=pretrackID0&&total_time[i]<1e+32)
	  {countParDec[i]++; pretrackID0=trackID;}
	index_ctrl=1;
      }
    }
    if (index_ctrl==0)
      {parNameDec[countPar]=particleName; countPar++;}
  }

  //ALL PROCESSES############################################################################
  
}

//---------------------------------------------------------------------//

void GeRecorder::EndOfRecordingStep(const G4Step *aStep)
{
  livT+=aStep->GetDeltaTime();
}

//---------------------------------------------------------------------//

void GeRecorder::ResetAllMembers()
{
  // Reset previous events' data
  evt_id = -999;
    
  x = -999.;
  y = -999.;
  z = -999.;
    
  x1 = -999.;
  y1 = -999.;
  z1 = -999.;
  r1 = -999.; h1 = -999.;

  dr = -999.; dh = -999.;
    
  hit_Ge=0;
  hit_Ge_b=0;

  edep_tot=0;
  hrec_tot=0;
  hrec_alpha=0;
  hrec_tot_lindh=0;
  hrec_elastic_tot=0;
  hrec_inelastic_tot=0;
  
  Etot = 0.;
  Er = 0.;
  Ea = 0.;
  Eqr = 0.;
  Eqr_lindh = 0.;
  Eqr_elastic = 0.;
  Eqr_inelastic = 0.;
  Eqr_elastic_lindh = 0.;
  Eqr_inelastic_lindh = 0.;

  for (int i=0;i<100;i++)
    Etot_tc[i] = 0.0;

  ntc = 0;

  edep_tot = 0.;
  hrec_tot = 0.;

  track_par = -1;
  track_inelastic = -1;
  track_elastic = -1;
  
  elastic_onoff = 0;
  inelastic_onoff = 0;

  locT=0;
  globT=0;
  propT=0;

  tdc=0;
  
  
}

//---------------------------------------------------------------------//

void GeRecorder::SetFileName(const G4String &fName)
{
  rootFileName = fName;
}

//---------------------------------------------------------------------//
