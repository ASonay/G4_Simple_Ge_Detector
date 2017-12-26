
#include "GeDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"


GeDetectorConstruction::GeDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
}


GeDetectorConstruction::~GeDetectorConstruction()
{
}


G4VPhysicalVolume* GeDetectorConstruction::Construct()
{
  //=========================================================================
  // Define Materials && Elements
  //=========================================================================
	
  G4String symbol;		// a = mass of a mole;
  G4double a, z, density;	// z = mean number of protons;
  G4int ncomponents;	// n = number of nucleons in an isotope;
  G4double fractionmass;

  G4Element *N  = new G4Element("Nitrogen", symbol="N",  z=7.,  a=14.01*g/mole);
  G4Element *O  = new G4Element("Oxigen",   symbol="O",  z=8.,  a=16.00*g/mole);
  G4Element *Pb = new G4Element("Lead",     symbol="Pb", z=82., a=207.2*g/mole);

  G4Material *Ge = new G4Material("Germanium", z=32., a=72.64*g/mole,  density=5.323*g/cm3);
  G4Material *Cu = new G4Material("Copper",    z=29., a=63.55*g/mole,  density=8.96*g/cm3);
  
  G4Material *Air = new G4Material("Air", density=1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material *Lead = new G4Material("Lead", density= 11.34*g/cm3, ncomponents=1);
  Lead->AddElement(Pb, fractionmass=1.0);
  
  /* THIS COULD BE USE FOR ENRICHED DETECTORS
  //GERMANIUM################################################################
  //Define isotopes
  G4int natoms, nIsotopes;	// n = number of nucleons in an isotope;
  G4double abundance;
  G4Isotope* Ge67 = new G4Isotope("Ge67",  32, 67, 66.93*g/mole);
  G4Isotope* Ge70 = new G4Isotope("Ge70",  32, 70, 69.92*g/mole);
  G4Isotope* Ge71 = new G4Isotope("Ge71",  32, 71, 70.92*g/mole);
  G4Isotope* Ge72 = new G4Isotope("Ge72",  32, 72, 71.92*g/mole);
  G4Isotope* Ge73 = new G4Isotope("Ge73",  32, 73, 73.0*g/mole);
  G4Isotope* Ge74 = new G4Isotope("Ge74",  32, 74, 74.0*g/mole);
  G4Isotope* Ge76 = new G4Isotope("Ge76",  32, 76, 76.0*g/mole);
   
  //Define element from isotopes
  G4Element* elGe = new G4Element("elGermanium",symbol="elGe",nIsotopes=7);
  elGe->AddIsotope(Ge67,abundance= 30.02*perCent);
  elGe->AddIsotope(Ge70,abundance= 20.83*perCent);
  elGe->AddIsotope(Ge71,abundance= 0.01*perCent);
  elGe->AddIsotope(Ge72,abundance= 27.54*perCent);
  elGe->AddIsotope(Ge73,abundance= 7.72*perCent);
  elGe->AddIsotope(Ge74,abundance= 6.28*perCent);
  elGe->AddIsotope(Ge76,abundance= 7.6*perCent);
  
  //Define material
  G4Material* Ge = new G4Material("Ge", density=5.32*g/cm3,1);
  Ge->AddElement(elGe,natoms=1);
  //GERMANIUM################################################################
  */

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //=========================================================================
  // World
  //=========================================================================
  
  solidWorld = new G4Box("World", 150.*cm, 150.*cm, 150.*cm);
  logicWorld = new G4LogicalVolume(solidWorld, Air, "World");
  physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, fCheckOverlaps);
  
  //=========================================================================
  // Ge Detector
  //=========================================================================

  float r1 = 3.15*cm;
  float r2 = 3.4*cm;
  float Ge_r1 = 3.0*cm;
  float Ge_r2 = 2.93*cm;
  float Ge_height = 3.0*cm;
  float Ge_diff   = 0.035*cm;


  solidCuHolder = new G4Tubs("CuHolder", 0.*cm, r2, 3.6*cm, 0.*deg, 360.*deg);//3.15 3.4
  logicCuHolder = new G4LogicalVolume(solidCuHolder, Cu, "CuHolder");
  physiCuHolder = new G4PVPlacement(0, G4ThreeVector(0., 0., 4.1*cm),
				    logicCuHolder, "CuHolder", logicWorld,
				    false, 0, fCheckOverlaps);


  solidGeCover = new G4Tubs("GeCover", 0.*cm, Ge_r1, Ge_height, 0.*deg, 360.*deg);
  logicGeCover = new G4LogicalVolume(solidGeCover, Ge, "GeCover");
  physiGeCover = new G4PVPlacement(0, G4ThreeVector(0., 0., 0), 
				   logicGeCover, "GeCover", logicCuHolder, 
				   false, 0, fCheckOverlaps);

  solidGeDet = new G4Tubs("GeDet", 0.*cm, Ge_r2, Ge_height - Ge_diff, 0.*deg, 360.*deg);
  logicGeDet = new G4LogicalVolume(solidGeDet, Ge, "GeDet");
  physiGeDet = new G4PVPlacement(0, G4ThreeVector(0., 0., - Ge_diff/2.), 
			      logicGeDet, "GeDet", logicGeCover, 
			      false, 0, fCheckOverlaps);
  
  solidCuHolderBase = new G4Tubs("CuHolderBase", 0.*cm, r1, 0.5*cm, 0.*deg, 360.*deg);
  logicCuHolderBase = new G4LogicalVolume(solidCuHolderBase, Cu, "CuHolderBase");
  physiCuHolderBase = new G4PVPlacement(0, G4ThreeVector(0., 0., 0),
					logicCuHolderBase, "CuHolderBase", logicWorld,
					false, 0, fCheckOverlaps);

  solidCuDisc = new G4Tubs("CuDisc", 0.*cm, 3.75*cm, 1.6*cm, 0.*deg, 360.*deg);
  logicCuDisc = new G4LogicalVolume(solidCuDisc, Cu, "CuDisc");
  physiCuDisc	= new G4PVPlacement(0, G4ThreeVector(0., 0., -2.1*cm),
				    logicCuDisc, "CuDisc", logicWorld,
				    false, 0, fCheckOverlaps);

  solidCuPillar = new G4Tubs("CuPillar", 0.*cm, 0.6*cm, 4.1*cm, 0.*deg, 360.*deg);
  logicCuPillar = new G4LogicalVolume(solidCuPillar, Cu, "CuPillar");
  physiCuPillar = new G4PVPlacement(0, G4ThreeVector(0., 0., -7.8*cm),
				    logicCuPillar, "CuPillar", logicWorld,
				    false, 0, fCheckOverlaps);

  
  //=========================================================================
  // Visualization attributes
  //=========================================================================
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *CuHolderVisAtt = new G4VisAttributes(G4Colour(1., 1., 0., 0.1));
  CuHolderVisAtt->SetForceSolid(true);
  logicCuHolder->SetVisAttributes(CuHolderVisAtt);
  logicCuHolderBase->SetVisAttributes(CuHolderVisAtt);
  logicCuDisc->SetVisAttributes(CuHolderVisAtt);
  logicCuPillar->SetVisAttributes(CuHolderVisAtt);

  G4VisAttributes *GeCoverVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8));
  GeCoverVisAtt->SetForceSolid(true);
  logicGeCover->SetVisAttributes(GeCoverVisAtt);

  G4VisAttributes *GeVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  GeVisAtt->SetForceSolid(true);
  logicGeDet->SetVisAttributes(GeVisAtt);

  return physiWorld;
}

