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
// $Id: DetectorConstruction.cc, 2016-03-29 $
// GEANT4 tag $Name: geant4.10.02.p01 $
// Developer: Chao Peng
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Digitization.hh"
#include "CalorimeterSD.hh"
#include "GasElectronMultiplierSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
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
: solidWorld(0),logicWorld(0),physiWorld(0), magField(0)
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
    SetDefaultMaterials();

    detectorMessenger = new DetectorMessenger(this);
    daq_system = new Digitization();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete daq_system;
    delete detectorMessenger;
}

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

    // define Elements
    G4Element* H  = new G4Element("Ele_Hydrogen", symbol= "H" , z= 1. , a= 1.01*g/mole);
    G4Element* C  = new G4Element("Ele_Carbon"  , symbol= "C" , z= 6. , a= 12.01*g/mole);
    G4Element* N  = new G4Element("Ele_Nitrogen", symbol= "N" , z= 7. , a= 14.01*g/mole);
    G4Element* O  = new G4Element("Ele_Oxygen"  , symbol= "O" , z= 8. , a= 16.00*g/mole);
    G4Element* Fe = new G4Element("Ele_Iron"    , symbol= "Fe", z= 26., a= 55.845*g/mole);
    G4Element* Pb = new G4Element("Ele_Lead"    , symbol= "Pb", z= 82., a= 207.20*g/mole);
    G4Element* Si = new G4Element("Ele_Silicon" , symbol= "Si", z= 14., a= 28.09*g/mole);
    G4Element* W  = new G4Element("Ele_Tungsten", symbol= "W" , z= 74., a= 183.84*g/mole);
    G4Element* Ar = new G4Element("Ele_Argon"   , symbol= "Ar", z= 18., a= 39.95*g/mole);
    G4Element* Al = new G4Element("Ele_Aluminum", symbol= "Al", z= 13., a= 26.98*g/mole);

    // define a material from elements.
    G4Material* PbWO4 = new G4Material("PbWO4", density= 8.300*g/cm3, ncomponents= 3);
    PbWO4->AddElement(Pb, natoms= 1);
    PbWO4->AddElement(W , natoms= 1);
    PbWO4->AddElement(O , natoms= 4);

    G4Material* Alum = new G4Material("Aluminum", density= 2.700*g/cm3, ncomponents= 1);
    Alum->AddElement(Al, natoms= 1);

    G4Material* Tung = new G4Material("Tungsten", density= 19.25*g/cm3, ncomponents= 1);
    Tung->AddElement(W, natoms= 1);

    G4Material* Iron = new G4Material("Iron", density= 7.874*g/cm3, ncomponents= 1);
    Iron->AddElement(Fe, natoms= 1);


    G4Material* SiO2 = new G4Material("quartz",density= 2.200*g/cm3, ncomponents= 2);
    SiO2->AddElement(Si, natoms= 1);
    SiO2->AddElement(O , natoms= 2);

    //Torlon4203L
    G4Material* Torlon = new G4Material("Torlon", density= 1.412*g/cm3, ncomponents= 5);
    Torlon->AddElement(C , natoms= 9);
    Torlon->AddElement(H , natoms= 4);
    Torlon->AddElement(N , natoms= 2);
    Torlon->AddElement(O , natoms= 3);
    Torlon->AddElement(Ar, natoms= 1);

    //scintillator
    G4Material* C9H10 = new G4Material("C9H10", density= 1.032*g/cm3, ncomponents= 2);
    C9H10->AddElement(C, natoms= 9);
    C9H10->AddElement(H, natoms= 10);

    //C2H6 gas
    G4Material* C2H6 = new G4Material("C2H6", density= 1.26e-3*g/cm3, ncomponents= 2, kStateGas, 298*kelvin, 1.*atmosphere);
    C2H6->AddElement(C, natoms= 2);
    C2H6->AddElement(H, natoms= 6);

    //Hydrogen gas
    G4Material* H2 =  new G4Material("H2 Gas", density= 8.988e-5*g/cm3, ncomponents=1, kStateGas, 25*kelvin, 83.02*pascal);
    H2->AddElement(H, natoms= 2);

    G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents= 2);
    H2O->AddElement(H, natoms= 2);
    H2O->AddElement(O, natoms= 1);
    //overwrite computed meanExcitationEnergy with ICRU recommended value
    H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

    //GEM Frame G10
    G4Material* G10 = new G4Material("NemaG10", density= 1.700*g/cm3, ncomponents= 4);
    G10->AddElement(Si, natoms= 1);
    G10->AddElement(O , natoms= 2);
    G10->AddElement(C , natoms= 3);
    G10->AddElement(H , natoms= 3);

    //Air
    G4Material* Air = new G4Material("Air"  , density= 1.225*mg/cm3, ncomponents= 2);
    Air->AddElement(N, fractionmass= 0.7);
    Air->AddElement(O, fractionmass= 0.3);

    //Lead Glass
    G4Material* PbGl = new G4Material("Lead Glass", density= 3.85*g/cm3, ncomponents= 2);
    PbGl->AddElement(Pb, fractionmass= 0.5316);
    PbGl->AddMaterial(SiO2, fractionmass= 0.4684);

    //Kapton
    G4Material* Kapton = new G4Material("Kapton", density= 1.42*g/cm3, ncomponents= 4);
    Kapton->AddElement(H, fractionmass= 0.0273);
    Kapton->AddElement(C, fractionmass= 0.7213);
    Kapton->AddElement(N, fractionmass= 0.0765);
    Kapton->AddElement(O, fractionmass= 0.1749);

    //Vacuum
    G4Material* Vacuum = new G4Material("Galactic", z= 1., a= 1.01*g/mole, density= universe_mean_density,
                                        kStateGas, 2.73*kelvin, 3.e-18*pascal);

    //Beam line vacuum
    G4Material* beam = new G4Material("Beam Line", density= 1.6e-11*g/cm3, ncomponents= 1,
                                      kStateGas, STP_Temperature, 1.33e-3*pascal);
    beam->AddMaterial(Air, fractionmass=1.);

    G4NistManager* man = G4NistManager::Instance();
    man->FindOrBuildMaterial("G4_NE_102");//("G4_SODIUM_IODIDE");

    // Print out material table
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    defaultMaterial = Vacuum;
}

void DetectorConstruction::SetDefaultMaterials()
{
    //default materials of the World
    CenterHyCalMaterial = G4Material::GetMaterial("PbWO4");
    OuterHyCalMaterial = G4Material::GetMaterial("Lead Glass");
    TargetMaterial = G4Material::GetMaterial("H2 Gas");
    CellMaterial = G4Material::GetMaterial("Kapton");
    HyCalBoxMaterial = G4Material::GetMaterial("Torlon");
    CollimatorMaterial = G4Material::GetMaterial("Tungsten");
    VacuumBoxMaterial = G4Material::GetMaterial("Aluminum");
    GEMFrameMaterial = G4Material::GetMaterial("NemaG10");
    GEMFoilMaterial = G4Material::GetMaterial("Galactic");
    GEMGasMaterial = G4Material::GetMaterial("C2H6");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Clean old geometry, if any
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // World
    solidWorld = new G4Box("World", WorldSizeXY, WorldSizeXY, WorldSizeZ);
    logicWorld = new G4LogicalVolume(solidWorld, defaultMaterial, "World");
    physiWorld = new G4PVPlacement(0, G4ThreeVector(0., 0. ,0.), logicWorld,
                                   "World", 0, false, 0);

    // Target
    /* disabled, events directly come from event generator
    solidTarget = new G4Tubs("Cell", 0.*cm, 0.4*cm, 2.*cm, 0, twopi);
    logicTarget = new G4LogicalVolume(solidTarget, TargetMaterial, TargetMaterial->GetName());
    physiTarget = new G4PVPlacement(0, G4ThreeVector(0.,0.,-250.*cm), logicTarget, "Cell", logicWorld, false, 0);
    */

    // Target Cell
    G4RotationMatrix rm;
    rm.rotateX(-90.*deg);
    G4double CellOR = CellR + CellThickness;

    G4Tubs *cellTube = new G4Tubs("Horizontal Tube", 0.*cm, CellOR, CellHalfL, 0, twopi);
    G4Tubs *neckTube = new G4Tubs("Inlet Tube", 0.3*cm, 0.3075*cm, NeckHalfL, 0, twopi);
    solidCellNeck = new G4SubtractionSolid("Gas inlet", neckTube, cellTube,
                                           G4Transform3D(rm, G4ThreeVector(0., 0., -(NeckHalfL-0.2*cm) - CellOR)));
    logicCellNeck = new G4LogicalVolume(solidCellNeck, CellMaterial, CellMaterial->GetName());
    physiCellNeck = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(0., (NeckHalfL-0.2*cm) + CellOR, TargetCenter)),
                                      logicCellNeck, "Neck", logicWorld, false, 0);

    solidCell = new G4Tubs("Tube", CellR, CellOR, CellHalfL, 0, twopi);
    logicCell = new G4LogicalVolume(solidCell, CellMaterial, CellMaterial->GetName());
    physiCell = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter), logicCell, "Tube", logicWorld, true, 0);

    // Windows
    solidCellWin = new G4Tubs("Windows", ApertureR, CellR, 5*um, 0, twopi);
    logicCellWin = new G4LogicalVolume(solidCellWin, CellMaterial, CellMaterial->GetName());
    physiWinIn = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter - CellHalfL + 5*um),
                                   logicCellWin, "Upstream Target Window", logicWorld, true, 0);
    physiWinOut = new G4PVPlacement(0, G4ThreeVector(0., 0., TargetCenter + CellHalfL - 5*um),
                                    logicCellWin, "Downstream Target Window", logicWorld, true, 0);

/*  TODO a detailed target chamber
    // Chamber Window
    G4Tubs* solidChamberWin = new G4Tubs("chamber windows", 0.3*cm, 17.5*cm, 3.75*um, 0, twopi);
    7G4LogicalVolume* logicChamberWin = new G4LogicalVolume(solidChamberWin, CellMaterial, CellMaterial->GetName());
    G4VPhysicalVolume* physiChamberWin = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - VacBoxtoHyCal - 458.*cm),
                                                           logicChamberWin, "Chamber Windows", logicWorld, false, 0);
*/

    //Vacuum Box Window, Flange, and Tube
    G4RotationMatrix rmVB;
    rmVB.rotateX(180.*deg);
    G4double EndDistance = ArcEndR - ArcDistance;
    G4Sphere *vacSphere = new G4Sphere("Sphere", ArcEndR - WinThickness, ArcEndR, 0.*deg, 360.*deg, 0.*deg, 90.*deg);
    G4Box *vacCut = new G4Box("Box", ArcEndR + 1.*mm, ArcEndR + 1.*mm, EndDistance);
    G4Tubs *vacHole = new G4Tubs("Hole", 0.*cm, 3.*cm, ArcEndR + 1.*mm, 0, twopi);
    G4SubtractionSolid *vacShell = new G4SubtractionSolid("Shell", vacSphere, vacCut);

    solidVacBoxWin = new G4SubtractionSolid("Aluminum Cover", vacShell, vacHole);
    logicVacBoxWin = new G4LogicalVolume(solidVacBoxWin, VacuumBoxMaterial, VacuumBoxMaterial->GetName());
    physiVacBoxWin = new G4PVPlacement(G4Transform3D(rmVB, G4ThreeVector(0., 0., HyCalCenter - 9.*cm - VacBoxtoHyCal + EndDistance)),
                                       logicVacBoxWin, "Vacuum Box End", logicWorld, false, 0);

    solidFlange = new G4Tubs("Flange", 1.9*cm, 3.*cm, 0.5*cm, 0, twopi);
    logicFlange = new G4LogicalVolume(solidFlange, VacuumBoxMaterial, VacuumBoxMaterial->GetName());
    physiFlange = new G4PVPlacement(0, G4ThreeVector(0., 0., (HyCalCenter - 9.*cm - VacBoxtoHyCal) - (ArcEndR - EndDistance)),
                                    logicFlange, "Vacuum Box Flange", logicWorld, false, 0);

    solidVacTube = new G4Tubs("solidVacTube", 1.8*cm, 1.9*cm, VacBoxtoHyCal + ArcDistance + 9.*cm, 0, twopi);
    logicVacTube = new G4LogicalVolume(solidVacTube, VacuumBoxMaterial, VacuumBoxMaterial->GetName());
    physiVacTube = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter), logicVacTube, "Vacuum Tube", logicWorld, false, 0);

    G4double rInner[] = {16.5*cm, 49.*cm, 49.*cm, 81.*cm, 81.*cm};
    G4double rOuter[] = {17.5*cm, 50.*cm, 50.*cm, 82.*cm, 82.*cm};
    G4double zPlane[] = {0.*cm, 5.*cm, 206.*cm, 211.*cm, 448.*cm};
    solidVacBox = new G4Polycone("Vacuum Box", 0, twopi, 5, zPlane, rInner, rOuter);
    logicVacBox = new G4LogicalVolume(solidVacBox, VacuumBoxMaterial, VacuumBoxMaterial->GetName());
    physiVacBox = new G4PVPlacement(0, G4ThreeVector(0., 0., HyCalCenter - 9.*cm - VacBoxtoHyCal - 2.*cm - 450.*cm),
                                                        logicVacBox, "Vacuum Box", logicWorld, false, 0);

    //GEM Frame
    G4Box *solidGEMFrame1 = new G4Box("GEMFrame1", 332.5*mm, 699.9*mm, 12.*mm);
    G4Box *solidGEMPiece1 = new G4Box("GEMPiece1", 275.*mm, 674.4*mm, 12.1*mm);
    G4Box *solidGEMPiece2 = new G4Box("GEMPiece2", 37.*mm, 29.5*mm, 12.2*mm);
    G4Tubs *solidGEMPipeHole = new G4Tubs("GEMPipeHole", 0., 22.*mm, 12.3*mm, 0, twopi);
    G4SubtractionSolid *solidGEMPiece = new G4SubtractionSolid("GEMPiece", solidGEMPiece1, solidGEMPiece2,
                                                               0, G4ThreeVector(-245.5*mm,0.,0.));
    G4SubtractionSolid *solidGEMFrame2 = new G4SubtractionSolid("GEMFrame2", solidGEMFrame1, solidGEMPiece);

    G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid("GEM_Frame", solidGEMFrame2, solidGEMPipeHole,
                                                               0, G4ThreeVector(-253.*mm,0.,0.));
    logicGEMFrame = new G4LogicalVolume(solidGEMFrame, defaultMaterial, GEMFrameMaterial->GetName());
    physiGEMFrame1 = new G4PVPlacement(0, G4ThreeVector(25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 3.*cm),
                                       logicGEMFrame, "GEM_Frame", logicWorld, false, 0);
    G4RotationMatrix rm2;
    rm2.rotateZ(180.*deg);
    physiGEMFrame2 = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector(-25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 7.*cm)),
                                       logicGEMFrame, "GEM_Frame", logicWorld, false,0);

    //GEM Foil
    G4Box *solidGEMPiece3 = new G4Box("GEMPiece3", 275.0*mm, 674.4*mm, 3.*mm);
    G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid("GEM_Foil", solidGEMPiece3, solidGEMPiece2,
                                                              0, G4ThreeVector(-245.5*mm,0.,0.));
    logicGEM = new G4LogicalVolume(solidGEMFoil, GEMFoilMaterial, GEMFoilMaterial->GetName());
    physiGEM1 = new G4PVPlacement(0, G4ThreeVector(25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 3.*cm),
                                                   logicGEM, "GEM_Foil", logicWorld, false, 0);
    physiGEM2 = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector(-25.3*cm, 0., HyCalCenter - VacBoxtoHyCal + 7.*cm)),
                                  logicGEM, "GEM_Foil", logicWorld, false,0);



    //HyCal box
    G4Box *solidBox1 = new G4Box("Box1", 70.*cm, 70.*cm, 60*cm);
    G4Box *solidBox2 = new G4Box("Box2", 66*cm, 66*cm, 59.6*cm);
    G4Tubs *solidBoxHole = new G4Tubs("BoxHole", 0., 25.*mm, 60.5*cm, 0, twopi);
    G4SubtractionSolid *solidBox3 = new G4SubtractionSolid("Box3", solidBox1, solidBox2);
    solidHyCalBox = new G4SubtractionSolid("HyCal_Box", solidBox3, solidBoxHole);
    logicHyCalBox = new G4LogicalVolume(solidHyCalBox, HyCalBoxMaterial, HyCalBoxMaterial->GetName());

    physiHyCalBox = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, HyCalCenter -9.*cm + 30.*cm),
                                      logicHyCalBox, "HyCal_Box", logicWorld, false, 0);


    //Collimators around the central hole
    solidColl = new G4Box("solidColBox", 1.025*cm, 1.025*cm, 5.*cm);
    logicColl = new G4LogicalVolume(solidColl, CollimatorMaterial, CollimatorMaterial->GetName());
    double pos_x[12] = {-3.075,-1.025, 1.025, 3.075,-3.075, 3.075,-3.075,3.075,-3.075,-1.025,1.025,3.075};
    double pos_y[12] = {-3.075,-3.075,-3.075,-3.075,-1.025,-1.025, 1.025,1.025, 3.075, 3.075,3.075,3.075};
    for(int i = 0; i < 12; ++i)
    {
        physiColl[i] = new G4PVPlacement(0, G4ThreeVector(pos_x[i]*cm, pos_y[i]*cm, HyCalCenter - 9.*cm - 5.1*cm),
                                         logicColl, "HyCal_Collimator", logicWorld, false, 0);
    }

    //HyCal
    G4Box *solidCentral = new G4Box("Container", 60.*cm, 60.*cm, 23.*cm);
    G4Box *solidHole = new G4Box("Hole", 2.*cm, 2.*cm, 23.1*cm);

    solidCalor = new G4SubtractionSolid("HyCal Container", solidCentral, solidHole);
    logicCalor = new G4LogicalVolume(solidCalor, defaultMaterial, "Virtual Container");
    physiCalor = new G4PVPlacement(0, G4ThreeVector(0.,0., HyCalCenter),
                                   logicCalor, "HyCal Container", logicWorld, false, 0);

    solidAbsorber = new G4Box ("Crystal Block", 1.025*cm, 1.025*cm, 90.*mm);
    logicAbsorber = new G4LogicalVolume(solidAbsorber, CenterHyCalMaterial, CenterHyCalMaterial->GetName());

    HyCalParameterisation *param = new HyCalParameterisation("config/module_list.txt");
    physiAbsorber = new G4PVParameterised("HyCal_Crystal", logicAbsorber, logicCalor,
                                          kUndefined, param->GetNumber(), param, false);

    //Sensitive detectors
    daq_system->RegisterModules(param);
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    CalorimeterSD* HyCalSD = new CalorimeterSD("eps/CalorimeterSD", daq_system);
    GasElectronMultiplierSD *GEMSD = new GasElectronMultiplierSD("eps/GasElectronMultiplierSD", daq_system);
    SDman->AddNewDetector(HyCalSD);
    SDman->AddNewDetector(GEMSD);
    logicAbsorber->SetSensitiveDetector(HyCalSD);
    logicGEM->SetSensitiveDetector(GEMSD);

    // Visualization attributes
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
    logicCalor->SetVisAttributes(G4VisAttributes::Invisible);

    G4VisAttributes* KaptonVisAtt = new G4VisAttributes(G4Colour(0.79, 0.53, 0.));
    KaptonVisAtt->SetVisibility(true);
    logicCell->SetVisAttributes(KaptonVisAtt);
    logicCellNeck->SetVisAttributes(KaptonVisAtt);
    logicCellWin->SetVisAttributes(KaptonVisAtt);
    //logicChamberWin->SetVisAttributes(KaptonVisAtt);

    G4VisAttributes* CrystalVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,1.0,0.1));
    CrystalVisAtt->SetVisibility(true);
    logicAbsorber->SetVisAttributes(CrystalVisAtt);

    G4VisAttributes* HyCalBoxVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0,0.3));
    HyCalBoxVisAtt->SetVisibility(false);
    logicHyCalBox->SetVisAttributes(HyCalBoxVisAtt);

    G4VisAttributes* GEMFrameVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.63));
    GEMFrameVisAtt->SetVisibility(false);
    logicGEMFrame->SetVisAttributes(GEMFrameVisAtt);

    G4VisAttributes* GEMFoilVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.1,0.3));
    GEMFoilVisAtt->SetVisibility(false);
    logicGEM->SetVisAttributes(GEMFoilVisAtt);

    G4VisAttributes* FlangeVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0));
    FlangeVisAtt->SetVisibility(false);
    logicFlange->SetVisAttributes(FlangeVisAtt);

    G4VisAttributes* TungstenVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    TungstenVisAtt->SetVisibility(false);
    logicColl->SetVisAttributes(TungstenVisAtt);

    G4VisAttributes* AluminumVisAtt= new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    AluminumVisAtt->SetVisibility(true);
    logicVacBoxWin->SetVisAttributes(AluminumVisAtt);
    logicVacTube->SetVisAttributes(AluminumVisAtt);
    logicVacBox->SetVisAttributes(AluminumVisAtt);

    //always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) TargetMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
    //apply a global uniform magnetic field along Z axis
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    if(magField) delete magField; //delete the existing magn field

    if(fieldValue != 0.) {
        magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
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
