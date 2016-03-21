//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, 2012-08-01 $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// Developer: Chao Peng
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "CalorimeterParameterisation.hh"
#include "LeadGlassPartParameterisation.hh"
#include "CalorimeterSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:AbsorberMaterial(0),Absorber2Material(0),TargetMaterial(0),CellMaterial(0),
 NeckMaterial(0),ColMaterial(0),VacBoxMaterial(0),defaultMaterial(0),
 solidWorld(0),logicWorld(0),physiWorld(0),
 solidNeck(0),logicNeck(0),physiNeck(0),
 solidCell(0),logicCell(0),physiCell(0),
 solidTarget(0),logicTarget(0),physiTarget(0),
 solidCellWin(0),logicCellWin(0),physiWinIn(0),physiWinOut(0),
 solidVacBoxEnd(0),logicVacBoxEnd(0),physiVacBoxEnd(0),
 solidFlange(0),logicFlange(0),physiFlange(0),
 logicGEMFrame(0),physiGEMFrame1(0),physiGEMFrame2(0),
 logicGEM(0),physiGEM1(0),physiGEM2(0),
 solidCalor(0),logicCalor(0),physiCalor(0),
 solidOuterCalor(0),logicOuterCalor(0),physiOuterCalor(0),
 solidAbsorber(0),logicAbsorber(0),physiAbsorber(0),
 solidAbsorber2(0),logicAbsorber2(0),physiAbsorber2(0),
 magField(0)
{
  // default parameter values of the calorimeter
  AbsorberThickness = 180.*mm;
  SurfaceDiff       = 6.5*cm;
  CalorSizeXY       = 34.85*cm;
  WorldSizeXY       = 150.*cm;
  WorldSizeZ        = 400.*cm;

  //parameter for target cell
  CellR = 1.2*cm;
  CellThickness = 0.003*cm;
  CellHalfL = 2.0*cm;
  NeckHalfL = 2.0*cm;
  ApertureR = 0.2*cm;


  //parameters for the vacuum box window
  fCylinderD = 163.5*cm;
  ArcEndR = 164.9*cm;
  ArcDistance = 24.91*cm;
  WinThickness = 2*mm;


  //parameter for the distance between Vacuum Box and central surface of HyCal
  //it's for the end of the whole vacuum box
  VacBoxtoHyCal = 55.*cm; //should be at least 35
  HyCalCenter = 265.*cm;
  TargetCenter = -300.*cm;


  // materials
  DefineMaterials();
  SetAbsorberMaterial("PbWO4");
  SetAbsorber2Material("Lead Glass");
  SetTargetMaterial("H2gas");
  SetCellMaterial("Kapton");
  SetNeckMaterial("IronMetal");
  SetColMaterial("Ma_Tungsten");
  SetVacBoxMaterial("Ma_Aluminum");
  SetGEMMaterial("NemaG10");

  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials

  G4String symbol;               //a=mass of a mole;
  G4double a, z, density;        //z=mean number of protons;  
  //G4int iz, n;                 //iz=number of protons  in an isotope; 
                                 // n=number of nucleons in an isotope;
  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Fe = new G4Element("Iron"    ,symbol="Fe", z= 26., a= 55.845*g/mole);
  G4Element* Pb = new G4Element("elLead"  ,symbol="Pb", z= 82., a= 207.20*g/mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Element* W  = new G4Element("Tungsten",symbol="W" , z= 74., a= 183.84*g/mole);
  G4Element* Ar = new G4Element("Argon", symbol="Ar", z= 18., a= 39.95*g/mole);
  G4Element* Al = new G4Element("Aluminum", symbol="Al", z= 13., a= 26.98*g/mole);


  //
  // define a material from elements.
  //
  G4Material* PbWO4 = new G4Material("PbWO4", density= 8.300*g/cm3, ncomponents=3);
  PbWO4->AddElement(Pb, natoms=1);
  PbWO4->AddElement(W, natoms=1);
  PbWO4->AddElement(O, natoms=4);

  G4Material* Alma = new G4Material("Ma_Aluminum", density = 2.700*g/cm3, ncomponents=1);
  Alma->AddElement(Al, natoms=1);

  G4Material* Tung = new G4Material("Ma_Tungsten", density = 19.25*g/cm3, ncomponents = 1);
  Tung->AddElement(W, natoms=1);

  G4Material* IronMetal = new G4Material("IronMetal", density = 7.874*g/cm3, ncomponents = 1);
  IronMetal->AddElement(Fe, natoms=1);


  G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  //Torlon4203L
  G4Material* Torlon = new G4Material("Torlon", density= 1.412*g/cm3, ncomponents=5);
  Torlon->AddElement(C, natoms=9);
  Torlon->AddElement(H, natoms=4);
  Torlon->AddElement(N, natoms=2);
  Torlon->AddElement(O, natoms=3);
  Torlon->AddElement(Ar, natoms=1);

  //scintillator
  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material("Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);

  //C2H6 gas
  G4Material* C2H6 = new G4Material("C2H6", density = 1.26e-3*g/cm3, ncomponents=2, kStateGas, 298*kelvin, 1.*atmosphere);
  C2H6->AddElement(C, natoms = 2);
  C2H6->AddElement(H, natoms = 6);

  //Hydrogen gas
  density = 10000* 4.192263e-7*g/cm3;
  G4Material* H2 =  new G4Material("H2gas", density, ncomponents=1, kStateGas, 25*kelvin, 83.02*pascal);
  H2->AddElement(H, natoms=2);

  G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  // overwrite computed meanExcitationEnergy with ICRU recommended value
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  G4Material* G10 = new G4Material("NemaG10", density= 1.700*g/cm3, ncomponents=4);
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //
  G4Material* Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass= 0.7);
  Air->AddElement(O, fractionmass= 0.3);

  G4Material* PbGl = new G4Material("Lead Glass", density= 3.85*g/cm3, ncomponents=2);
  PbGl->AddElement(Pb, fractionmass= 0.5316);
  PbGl->AddMaterial(SiO2, fractionmass= 0.4684);

  //Kapton
  G4Material* Kapton = new G4Material("Kapton", density= 1.42*g/cm3, ncomponents=4);
  Kapton->AddElement(H, fractionmass= 0.0273);
  Kapton->AddElement(C, fractionmass= 0.7213);
  Kapton->AddElement(N, fractionmass= 0.0765);
  Kapton->AddElement(O, fractionmass= 0.1749);

  //
  // examples of vacuum
  //

  G4Material* Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                                      kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Material* beam = new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                                    kStateGas, STP_Temperature, 2.e-2*bar);
  beam->AddMaterial(Air, fractionmass=1.);

  //
  // or use G4-NIST materials data base
  //
  G4NistManager* man = G4NistManager::Instance();
  man->FindOrBuildMaterial("G4_NE_102");//("G4_SODIUM_IODIDE");

  //
  // print table
  //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //default materials of the World
  defaultMaterial  = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // World
  //
  solidWorld = new G4Box("World", WorldSizeXY, WorldSizeXY, WorldSizeZ);
  logicWorld = new G4LogicalVolume(solidWorld,
                                   defaultMaterial,
                                   "World");
  physiWorld = new G4PVPlacement(0,
                                 G4ThreeVector(0., 0. ,0.),
                                 logicWorld,
                                 "World",
                                 0,
                                 false,
                                 0);
  //
  // Target
  //
/*
  solidTarget = new G4Tubs("Cell", 0.*cm, 0.4*cm, 2.*cm, 0, twopi);
  logicTarget = new G4LogicalVolume(solidTarget, defaultMaterial, "H2gas");
  physiTarget = new G4PVPlacement(0, G4ThreeVector(0.,0.,-250.*cm), logicTarget, "Cell", logicWorld, false, 0);
*/


  //
  // Cell
  //

  G4RotationMatrix rm;
  rm.rotateX(-90.*deg);
  G4double CellOR = CellR + CellThickness;

  solidCelltube = new G4Tubs("Horizontal Tube", 0.*cm, CellOR, CellHalfL, 0, twopi);
  solidNecktube = new G4Tubs("Inlet Tube", 0.3*cm, 0.3075*cm, NeckHalfL, 0, twopi);
  solidNeck = new G4SubtractionSolid("Gas inlet", solidNecktube, solidCelltube, G4Transform3D(rm, G4ThreeVector(0., 0., -(NeckHalfL-0.2*cm) - CellOR)));  
  logicNeck = new G4LogicalVolume(solidNeck, CellMaterial, CellMaterial->GetName());
  physiNeck = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(0., (NeckHalfL-0.2*cm) + CellOR, TargetCenter)), logicNeck, "Neck", logicWorld, false, 0);

  solidCell = new G4Tubs("Tube", CellR, CellOR, CellHalfL, 0, twopi);
  logicCell = new G4LogicalVolume(solidCell, CellMaterial, CellMaterial->GetName());
  physiCell = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter), logicCell, "Tube", logicWorld, true, 0);

  //
  // Windows
  //

  solidCellWin  = new G4Tubs("Windows", ApertureR, CellR, 5*um, 0, twopi);
  logicCellWin  = new G4LogicalVolume(solidCellWin, CellMaterial, CellMaterial->GetName());
  physiWinIn  = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter - CellHalfL + 5*um), logicCellWin, "upstream Window", logicWorld, true, 0);
  physiWinOut  = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter + CellHalfL - 5*um), logicCellWin, "downstream Window", logicWorld, true, 0);
/*
  // Chamber Window
  G4Tubs* solidChamberWin = new G4Tubs("chamber windows", 0.3*cm, 17.5*cm, 3.75*um, 0, twopi);
  G4LogicalVolume* logicChamberWin = new G4LogicalVolume(solidChamberWin, CellMaterial, CellMaterial->GetName());
  G4VPhysicalVolume* physiChamberWin = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - VacBoxtoHyCal - 458.*cm), logicChamberWin, "Chamber Windows", logicWorld, false, 0);
*/

/*
  //position test
  G4Box* solidPos = new G4Box("position", 35.*cm, 35.*cm, 0.01*mm);
  G4LogicalVolume* logicPos  = new G4LogicalVolume(solidPos, defaultMaterial, "position");
  G4VPhysicalVolume* physiPos = new G4PVPlacement(0, G4ThreeVector(0., 0., 240.89*cm), logicPos, "position", logicWorld, false, 0);
*/



  //
  //Vacuum Box End, Flange, and Tube
  //
  G4RotationMatrix rmVB;
  rmVB.rotateX(180.*deg);
//  G4double EndDistance = sqrt(ArcEndR*ArcEndR - fCylinderD*fCylinderD/4.);
  G4double EndDistance = ArcEndR - ArcDistance;
  solidVacBoxSph = new G4Sphere("Subtraction Sphere", ArcEndR - WinThickness, ArcEndR, 0.*deg, 360.*deg, 0.*deg, 90.*deg);
  solidVacBoxBox = new G4Box("Subtraction Box", ArcEndR + 1.*mm, ArcEndR + 1.*mm, EndDistance);
  solidVacBoxTub = new G4Tubs("Subtraction Tube for hole", 0.*cm, 3.*cm, ArcEndR + 1.*mm, 0, twopi);
  solidVacBoxShl = new G4SubtractionSolid("Shell shape", solidVacBoxSph, solidVacBoxBox);
  solidVacBoxEnd = new G4SubtractionSolid("Aluminum Cover", solidVacBoxShl, solidVacBoxTub);
  logicVacBoxEnd = new G4LogicalVolume(solidVacBoxEnd, defaultMaterial, VacBoxMaterial->GetName());
  physiVacBoxEnd = new G4PVPlacement(G4Transform3D(rmVB, G4ThreeVector(0., 0., HyCalCenter - 9.*cm - VacBoxtoHyCal + EndDistance)), logicVacBoxEnd, "Vacuum Box End", logicWorld, false, 0);

  solidFlange = new G4Tubs("Flange", 1.9*cm, 3.*cm, 0.5*cm, 0, twopi);
  logicFlange = new G4LogicalVolume(solidFlange, VacBoxMaterial, VacBoxMaterial->GetName());
  physiFlange = new G4PVPlacement(0, G4ThreeVector(0., 0., (HyCalCenter - 9.*cm - VacBoxtoHyCal) - (ArcEndR - EndDistance)), logicFlange, "Vacuum Box Flange", logicWorld, false, 0);


  solidVacTube = new G4Tubs("solidVacTube", 1.8*cm, 1.9*cm, VacBoxtoHyCal + ArcDistance + 9.*cm, 0, twopi);
  logicVacTube = new G4LogicalVolume(solidVacTube, VacBoxMaterial, VacBoxMaterial->GetName());
  physiVacTube = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter), logicVacTube, "Vacuum Tube", logicWorld, false, 0);


  G4Tubs* solidVacBox1 = new G4Tubs("solid_VacBox1", 17.5*cm, 50.*cm, 1.*cm, 0, twopi);
  G4LogicalVolume* logicVacBox1 = new G4LogicalVolume(solidVacBox1, VacBoxMaterial, VacBoxMaterial->GetName());
  G4VPhysicalVolume* physiVacBox1 = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - VacBoxtoHyCal - 458.*cm), logicVacBox1, "Vacuum Box 1", logicWorld, false, 0);

  G4Tubs* solidVacBox2 = new G4Tubs("solid_VacBox2", 49.*cm, 50.*cm, 102.*cm, 0, twopi);
  G4LogicalVolume* logicVacBox2 = new G4LogicalVolume(solidVacBox2, VacBoxMaterial, VacBoxMaterial->GetName());
  G4VPhysicalVolume* physiVacBox2 = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - VacBoxtoHyCal - 355.*cm), logicVacBox2, "Vacuum Box 2", logicWorld, false, 0);

  G4Tubs* solidVacBox3 = new G4Tubs("solid_VacBox3", 50.*cm, 82.*cm, 1.*cm, 0, twopi);
  G4LogicalVolume* logicVacBox3 = new G4LogicalVolume(solidVacBox3, VacBoxMaterial, VacBoxMaterial->GetName());
  G4VPhysicalVolume* physiVacBox3 = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - VacBoxtoHyCal - 252.*cm), logicVacBox3, "Vacuum Box 3", logicWorld, false, 0);

  G4Tubs* solidVacBox4 = new G4Tubs("solid_VacBox4", 81.*cm, 82.*cm, 120.*cm, 0, twopi);
  G4LogicalVolume* logicVacBox4 = new G4LogicalVolume(solidVacBox4, VacBoxMaterial, VacBoxMaterial->GetName());
  G4VPhysicalVolume* physiVacBox4 = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - 9.*cm - VacBoxtoHyCal - 2.*cm - 120.*cm), logicVacBox4, "Vacuum Box 4", logicWorld, false, 0);


  G4Box *solidGEMFrame1 = new G4Box("GEMFrame1", 332.5*mm, 699.9*mm, 12.*mm);
  G4Box *solidGEMPiece1 = new G4Box("GEMPiece1", 275.*mm, 674.4*mm, 12.1*mm);
  G4Box *solidGEMPiece2 = new G4Box("GEMPiece2", 37.*mm, 29.5*mm, 12.2*mm);
  G4Tubs *solidGEMPipeHole = new G4Tubs("GEMPipeHole", 0., 22.*mm, 12.3*mm, 0, twopi);
  G4SubtractionSolid *solidGEMPiece = new G4SubtractionSolid("GEMPiece", solidGEMPiece1, solidGEMPiece2, 0, G4ThreeVector(-245.5*mm,0.,0.));
  G4SubtractionSolid *solidGEMFrame2 = new G4SubtractionSolid("GEMFrame2", solidGEMFrame1, solidGEMPiece);
  G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid("GEM_Frame", solidGEMFrame2, solidGEMPipeHole, 0, G4ThreeVector(-253.*mm,0.,0.));
  logicGEMFrame = new G4LogicalVolume(solidGEMFrame, GEMMaterial, GEMMaterial->GetName());
  physiGEMFrame1 = new G4PVPlacement(0, G4ThreeVector(25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 3.*cm), logicGEMFrame, "GEM_Frame", logicWorld, false, 0);
  G4RotationMatrix rm2;
  rm2.rotateZ(180.*deg);
  physiGEMFrame2 = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector(-25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 7.*cm)), logicGEMFrame, "GEM_Frame", logicWorld, false,0);
  //Position Detectors

  G4Box *solidGEMPiece3 = new G4Box("GEMPiece3", 275.*mm, 674.4*mm, 3.*mm);
  G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid("GEM_Foil", solidGEMPiece3, solidGEMPiece2, 0, G4ThreeVector(-245.5*mm,0.,0.));
  logicGEM = new G4LogicalVolume(solidGEMFoil, GEMMaterial, GEMMaterial->GetName());
  physiGEM1 = new G4PVPlacement(0, G4ThreeVector(25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 3.*cm), logicGEM, "GEM_Foil", logicWorld, false, 0);
  physiGEM2 = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector(-25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 7.*cm)), logicGEM, "GEM_Foil", logicWorld, false,0);


  //Collimator around the central hole
  solidColl = 0;
  logicColl = 0;
  for(int i = 0; i < 12; ++i) {
    physiColl[i] = 0;
  }


  solidColl = new G4Box("solidColBox", 1.025*cm, 1.025*cm, 5.*cm);
  logicColl = new G4LogicalVolume(solidColl, ColMaterial, ColMaterial->GetName());
  for(int i = 0; i < 4; ++i) {
    physiColl[i] = new G4PVPlacement(0, G4ThreeVector((-2.05 - 1.025 + 2.05*(double)i)*cm, (2.05 + 1.025)*cm, HyCalCenter - 9.*cm - 5.1*cm), logicColl, "Collimators", logicWorld, false, 0);
  }

  for(int i = 4; i < 6; ++i) {
    physiColl[i] = new G4PVPlacement(0, G4ThreeVector((-2.05 - 1.025 + 3*2.05*(double)(i-4))*cm, 1.025*cm, HyCalCenter - 9.*cm - 5.1*cm), logicColl, "Collimators", logicWorld, false, 0);
  }

  for(int i = 6; i < 8; ++i) {
    physiColl[i] = new G4PVPlacement(0, G4ThreeVector((-2.05 - 1.025 + 3*2.05*(double)(i-6))*cm, -1.025*cm, HyCalCenter - 9.*cm - 5.1*cm), logicColl, "Collimators", logicWorld, false, 0);
  }

  for(int i = 8; i < 12; ++i) {
    physiColl[i] = new G4PVPlacement(0, G4ThreeVector((-2.05 - 1.025 + 2.05*(double)(i-8))*cm, -(2.05 + 1.025)*cm, HyCalCenter - 9.*cm - 5.1*cm), logicColl, "Collimators", logicWorld, false, 0);
  }


  //HyCal box

  G4Box *solidBox1 = new G4Box("Box1", 70.*cm, 70.*cm, 60*cm);
  G4Box *solidBox2 = new G4Box("Box2", 66*cm, 66*cm, 59.6*cm);
  G4Tubs *solidBoxHole = new G4Tubs("BoxHole", 0., 25.*mm, 60.5*cm, 0, twopi);
  G4SubtractionSolid *solidBox3 = new G4SubtractionSolid("Box3", solidBox1, solidBox2);
  solidHyCalBox = new G4SubtractionSolid("HyCal_Box", solidBox3, solidBoxHole);
  logicHyCalBox = new G4LogicalVolume(solidHyCalBox, NeckMaterial, NeckMaterial->GetName());

  physiHyCalBox = new G4PVPlacement(0,                 //no rotation
                                   G4ThreeVector(0.*cm, 0.*cm, HyCalCenter -9.*cm + 30.*cm),      //at (0,0,0)
                                   logicHyCalBox,        //its logical volume
                                   "HyCal_Box",       //its name
                                   logicWorld,        //its mother  volume
                                   false,             //no boolean operation
                                   0);                //copy number


  //
  //HyCal_central part
  //
  G4Box *solidCentral = new G4Box("Central",             //its name
                         CalorSizeXY, CalorSizeXY, 90.*mm);//size

  G4Box *solidHole = new G4Box("Hole", 2.05*cm, 2.05*cm, 91.*mm);

  solidCalor = new G4SubtractionSolid("Central part with hole", solidCentral, solidHole);
  logicCalor = new G4LogicalVolume(solidCalor,      //its solid
                                   defaultMaterial, //its material
                                   "Central part");  //its name
  physiCalor = new G4PVPlacement(0,                 //no rotation
                                 G4ThreeVector(0.,0., HyCalCenter),      //at (0,0,0)
                                 logicCalor,        //its logical volume
                                 "Central Frame",     //its name
                                 logicWorld,        //its mother  volume
                                 false,             //no boolean operation
                                 0);                //copy number



  solidAbsorber = new G4Box ("Crystal Block", 1.025*cm, 1.025*cm, 90.*mm);
  logicAbsorber = new G4LogicalVolume(solidAbsorber, AbsorberMaterial, AbsorberMaterial->GetName());

  CalorimeterParameter = new CalorimeterParameterisation(34,        //NbofBlocks in a line
                                                         34,        //NbofBlocks in a column
                                                         CalorSizeXY, //Calorimeter Size
                             G4ThreeVector(0.*cm, 0.*cm, 0.*cm), //center position
                             1.025*cm,  //dimension x
                             1.025*cm); //dimension y
  physiAbsorber = new G4PVParameterised("HyCal_Crystal" ,
                                        logicAbsorber,
                                        logicCalor,
                                        kXAxis,
                                        1152,                      //NbofBlocks = 34*34-4
                                        CalorimeterParameter,
                                        false);

 //
  // Lead Glass Part
  //


  G4Box *solidOuterBox = new G4Box("OuterBox", 60.*cm, 60.*cm, 225.*mm);
  G4Box *solidCentralBox = new G4Box("CentralBox", CalorSizeXY, CalorSizeXY, 226.*mm);

  solidOuterCalor = new G4SubtractionSolid("Outer part with hole", solidOuterBox, solidCentralBox);

  logicOuterCalor = new G4LogicalVolume(solidOuterCalor, defaultMaterial, "Outer part");

  physiOuterCalor = new G4PVPlacement(0,                 //no rotation
                                 G4ThreeVector(0.*cm, 0.*cm, HyCalCenter + 13.5*cm - SurfaceDiff),      //at (0,0,0)
                                 logicOuterCalor,        //its logical volume
                                 "Upper_Outerpart",     //its name
                                 logicWorld,        //its mother  volume
                                 false,             //no boolean operation
                                 0);                //copy number

  solidAbsorber2 = new G4Box ("Crystal Block", 1.91*cm, 1.91*cm, 225.*mm);
  logicAbsorber2 = new G4LogicalVolume(solidAbsorber2, Absorber2Material, Absorber2Material->GetName());

  LeadGlassPartParameter = new LeadGlassPartParameterisation(24,        //NbofBlocks in a line
                                                             6,        //NbofBlocks in a column
                                                             G4ThreeVector(0.*cm, 0.*cm, 0.*cm), //center position
                                                             1.91*cm,  //dimension x
                                                             1.91*cm); //dimension y

  physiAbsorber2 = new G4PVParameterised("HyCal_Leadglass" ,
                                         logicAbsorber2,
                                         logicOuterCalor,
                                         kXAxis,
                                         576,                      //NbofBlocks = 34*34-4
                                         LeadGlassPartParameter,
                                         false);




  //
  //Sensitive detectors
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String trackerChamberSDname = "eps/CalorimeterSD";
  CalorimeterSD* aTrackerSD = new CalorimeterSD(trackerChamberSDname);
  SDman->AddNewDetector( aTrackerSD );
  logicAbsorber->SetSensitiveDetector(aTrackerSD);
  logicAbsorber2->SetSensitiveDetector(aTrackerSD);
  logicGEM->SetSensitiveDetector(aTrackerSD);



  //
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  logicCalor->SetVisAttributes (G4VisAttributes::Invisible);
  logicOuterCalor->SetVisAttributes (G4VisAttributes::Invisible);
/*
  logicFlange->SetVisAttributes (G4VisAttributes::Invisible);
  logicGEM->SetVisAttributes (G4VisAttributes::Invisible);
  logicColl->SetVisAttributes (G4VisAttributes::Invisible);
  logicAbsorber->SetVisAttributes (G4VisAttributes::Invisible);
  logicAbsorber2->SetVisAttributes (G4VisAttributes::Invisible);
  logicVacTube->SetVisAttributes (G4VisAttributes::Invisible);
  logicVacBoxEnd->SetVisAttributes (G4VisAttributes::Invisible);
*/
//  logicHyCalBox->SetVisAttributes (G4VisAttributes::Invisible);


  G4VisAttributes* KaptonVisAtt = new G4VisAttributes(G4Colour(0.79, 0.53, 0.));
  KaptonVisAtt->SetVisibility(true);
  logicNeck->SetVisAttributes(KaptonVisAtt);
  logicCell->SetVisAttributes(KaptonVisAtt);
  logicCellWin->SetVisAttributes(KaptonVisAtt);
  //logicChamberWin->SetVisAttributes(KaptonVisAtt);

  G4VisAttributes* CrystalVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.35));
  CrystalVisAtt->SetVisibility(true);
  logicAbsorber->SetVisAttributes(CrystalVisAtt);

  G4VisAttributes* LeadGlassVisAtt= new G4VisAttributes(G4Colour(0.2,0.0,1.0,0.3));
  LeadGlassVisAtt->SetVisibility(true);
  logicAbsorber2->SetVisAttributes(LeadGlassVisAtt);

  G4VisAttributes* HyCalBoxVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0,0.3));
  HyCalBoxVisAtt->SetVisibility(true);
  logicHyCalBox->SetVisAttributes(HyCalBoxVisAtt);

  G4VisAttributes* GEMFrameVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.63));
  GEMFrameVisAtt->SetVisibility(true);
  logicGEMFrame->SetVisAttributes(GEMFrameVisAtt);

  G4VisAttributes* GEMFoilVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.1,0.3));
  logicGEM->SetVisAttributes(GEMFoilVisAtt);

  G4VisAttributes* FlangeVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  FlangeVisAtt->SetVisibility(true);
  logicFlange->SetVisAttributes(FlangeVisAtt);

  G4VisAttributes* TungstenVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  TungstenVisAtt->SetVisibility(true);
  logicColl->SetVisAttributes(TungstenVisAtt);

  G4VisAttributes* AluminumVisAtt= new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  AluminumVisAtt->SetVisibility(true);
  logicVacBoxEnd->SetVisAttributes(AluminumVisAtt);
  logicVacTube->SetVisAttributes(AluminumVisAtt);

  logicVacBox1->SetVisAttributes(AluminumVisAtt);
  logicVacBox2->SetVisAttributes(AluminumVisAtt);
  logicVacBox3->SetVisAttributes(AluminumVisAtt);
  logicVacBox4->SetVisAttributes(AluminumVisAtt);


/*
  // Below are vis attributes that permits someone to test / play 
  // with the interactive expansion / contraction geometry system of the
  // vis/OpenInventor driver :
 {G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  simpleBoxVisAtt->SetVisibility(true);
  delete logicCalor->GetVisAttributes();
  logicCalor->SetVisAttributes(simpleBoxVisAtt);}

 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  logicLayer->SetVisAttributes(atb);}
  
 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  atb->SetForceSolid(true);
  logicAbsorber->SetVisAttributes(atb);}
  
 {//Set opacity = 0.2 then transparency = 1 - 0.2 = 0.8
  G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.2));
  atb->SetForceSolid(true);
  logicGap->SetVisAttributes(atb);}
*/

  //
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) AbsorberMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorber2Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) Absorber2Material = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) TargetMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCellMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) CellMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNeckMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) NeckMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetColMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) ColMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetVacBoxMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) VacBoxMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGEMMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) GEMMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeXY(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  CalorSizeXY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField; //delete the existing magn field

  if(fieldValue!=0.)    // create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
