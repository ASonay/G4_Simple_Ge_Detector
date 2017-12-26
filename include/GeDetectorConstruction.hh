
#ifndef GeDetectorConstruction_h
#define GeDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;


class GeDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  GeDetectorConstruction();
  ~GeDetectorConstruction();

  G4VPhysicalVolume *Construct();

private:
  G4bool  fCheckOverlaps;
  // Logical volumes & Physical volumes
  G4Box             *solidWorld;
  G4LogicalVolume   *logicWorld;
  G4VPhysicalVolume *physiWorld;

  //------------------------------------------------------------
  //Ge----------------------------------------------------------
  G4Tubs            *solidCuHolder, *solidCuHolderBase;
  G4LogicalVolume   *logicCuHolder, *logicCuHolderBase;
  G4VPhysicalVolume *physiCuHolder, *physiCuHolderBase;
  
  G4Tubs            *solidGeCover, *solidGeDet;
  G4LogicalVolume   *logicGeCover, *logicGeDet;
  G4VPhysicalVolume *physiGeCover, *physiGeDet;

  G4Tubs            *solidCuDisc, *solidCuPillar;
  G4LogicalVolume   *logicCuDisc, *logicCuPillar;
  G4VPhysicalVolume *physiCuDisc, *physiCuPillar;
};

#endif
