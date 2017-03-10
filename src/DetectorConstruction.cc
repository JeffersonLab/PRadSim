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
// DetectorConstruction.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "HyCalParameterisation.hh"
#include "StandardDetectorSD.hh"

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"

#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4VisAttributes.hh"

#include <cmath>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4String conf) : G4VUserDetectorConstruction(), fConfig(conf)
{
    if (fConfig != "prad" && fConfig != "drad")
        fConfig = "prad";

    fTargetCenter = -300.0 * cm;
    fTargetR = 15.0 * cm;
    fTargetHalfL = 3.0 * cm;

    fRecoilDetNSeg = 36;
    fRecoilDetHalfL = 2.0 * cm;
    fRecoilDetThickness = 2.0 * mm;

    fGEM1Center = 200.0 * cm;
    fGEM2Center = 250.0 * cm;

    fSciPlaneCenter = 265.0 * cm;

    fCrystalSurf = 300.0 * cm;

    detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    // Material
    std::map<G4String, G4VisAttributes *> mVisAtt;

    G4int z, n;
    G4double a;
    G4double density;
    G4int ncomponents, natoms;
    G4double fractionmass;

    G4NistManager *pNM = G4NistManager::Instance();

    // Define elements
    G4Element *H  = pNM->FindOrBuildElement(z = 1);
    G4Element *C  = pNM->FindOrBuildElement(z = 6);
    G4Element *N  = pNM->FindOrBuildElement(z = 7);
    G4Element *O  = pNM->FindOrBuildElement(z = 8);
    G4Element *Al = pNM->FindOrBuildElement(z = 13);
    G4Element *Si = pNM->FindOrBuildElement(z = 14);
    G4Element *Ar = pNM->FindOrBuildElement(z = 18);
    G4Element *Cu = pNM->FindOrBuildElement(z = 29);
    G4Element *W  = pNM->FindOrBuildElement(z = 74);
    G4Element *Pb = pNM->FindOrBuildElement(z = 86);
    G4Isotope *H2 = new G4Isotope("H2", z = 1, n = 2, a = 2.0141 * g / mole);
    G4Element *D = new G4Element("Deuterium", "D", ncomponents = 1);
    D->AddIsotope(H2, 1.0);

    // Space Vacuum
    G4Material *Galaxy = new G4Material("Galaxy", density = universe_mean_density, ncomponents = 1, kStateGas, 0.1 * kelvin, 1.0e-19 * pascal);
    Galaxy->AddElement(H, fractionmass = 1.0);

    // Air
    G4Material *Air = new G4Material("Air", density = 1.292 * mg / cm3, ncomponents = 2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);

    // Air vacuum of 1.e-6 torr at room temperature, 1 atmosphere = 760 torr
    G4Material *Vacuum = new G4Material("Vacuum", density = 1.0e-6 / 760.0 * 1.292 * mg / cm3, ncomponents = 1, kStateGas, STP_Temperature, 1.0e-6 / 760.0 * atmosphere);
    Vacuum->AddMaterial(Air, fractionmass = 1.0);

    // Hydrogen Gas (T=19.5K, P=470mTorr)
    G4Material *H2Gas =  new G4Material("H2 Gas", density = 0.47 / 760.0 * 273.15 / 19.5 * 0.08988 * mg / cm3, ncomponents = 1, kStateGas, 25.0 * kelvin, 0.6 / 760.0 * atmosphere);
    H2Gas->AddElement(H, natoms = 2);

    // Deuteron Gas 
    G4Material *D2Gas =  new G4Material("D2 Gas", density = 0.47 / 760.0 * 273.15 / 19.5 * 0.1796 * mg / cm3, ncomponents = 1, kStateGas, 25.0 * kelvin, 0.6 / 760.0 * atmosphere);
    D2Gas->AddElement(D, natoms = 2);

    // Copper C101
    G4Material *Copper = new G4Material("Copper", density = 8.92 * g / cm3, ncomponents = 1);
    Copper->AddElement(Cu, natoms = 1);

    // Kapton
    G4Material *Kapton = new G4Material("Kapton", density = 1.42 * g / cm3, ncomponents = 4);
    Kapton->AddElement(H, fractionmass = 0.0273);
    Kapton->AddElement(C, fractionmass = 0.7213);
    Kapton->AddElement(N, fractionmass = 0.0765);
    Kapton->AddElement(O, fractionmass = 0.1749);

    // Si
    G4Material *Silicon = new G4Material("Silicon", density = 2.329 * g / cm3, ncomponents = 1);
    Silicon->AddElement(Si, natoms = 1);

    // Aluminum
    G4Material *Aluminum = new G4Material("Aluminum", density = 2.700 * g / cm3, ncomponents = 1);
    Aluminum->AddElement(Al, natoms = 1);

    // GEM Frame G10
    G4Material *NemaG10 = new G4Material("NemaG10", density = 1.700 * g / cm3, ncomponents = 4);
    NemaG10->AddElement(Si, natoms = 1);
    NemaG10->AddElement(O, natoms = 2);
    NemaG10->AddElement(C, natoms = 3);
    NemaG10->AddElement(H, natoms = 3);

    // CO2 Gas
    G4Material *CO2 = new G4Material("CO2", density = 1.842e-3 * g / cm3, ncomponents = 2);
    CO2->AddElement(C, natoms = 1);
    CO2->AddElement(O, natoms = 2);

    // Ar/CO2 Gas
    G4Material *ArCO2 = new G4Material("ArCO2", density = 1.715e-3 * g / cm3, ncomponents = 2);
    ArCO2->AddElement(Ar, fractionmass = 0.7);
    ArCO2->AddMaterial(CO2, fractionmass = 0.3);

    // Scintillator EJ204
    G4Material *EJ204 = new G4Material("EJ204", density = 1.032 * g / cm3, ncomponents = 2);
    EJ204->AddElement(H, natoms = 521);
    EJ204->AddElement(C, natoms = 474);

    // Torlon4203L
    G4Material *Torlon = new G4Material("Torlon", density = 1.412 * g / cm3, ncomponents = 5);
    Torlon->AddElement(C, natoms = 9);
    Torlon->AddElement(H, natoms = 4);
    Torlon->AddElement(N, natoms = 2);
    Torlon->AddElement(O, natoms = 3);
    Torlon->AddElement(Ar, natoms = 1);

    // Tungsten
    G4Material *Tungsten = new G4Material("Tungsten", density = 19.25 * g / cm3, ncomponents = 1);
    Tungsten->AddElement(W, natoms = 1);

    // PbWO4 Crystal
    G4Material *PbWO4 = new G4Material("PbWO4", density = 8.300 * g / cm3, ncomponents = 3);
    PbWO4->AddElement(Pb, natoms = 1);
    PbWO4->AddElement(W, natoms = 1);
    PbWO4->AddElement(O, natoms = 4);

    // Lead Glass
    G4Material *SiO2 = new G4Material("Quartz", density = 2.200 * g / cm3, ncomponents = 2);
    SiO2->AddElement(Si, natoms = 1);
    SiO2->AddElement(O, natoms = 2);
    G4Material *PbGlass = new G4Material("Lead Glass", density = 3.85 * g / cm3, ncomponents = 2);
    PbGlass->AddElement(Pb, fractionmass = 0.5316);
    PbGlass->AddMaterial(SiO2, fractionmass = 0.4684);

    // Print out material table
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    mVisAtt[Vacuum->GetName()] = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5)); // Vacuum;
    mVisAtt[H2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());
    mVisAtt[D2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());
    mVisAtt[Copper->GetName()] = new G4VisAttributes(G4Colour::Brown());
    mVisAtt[Kapton->GetName()] = new G4VisAttributes(G4Colour::Brown());
    mVisAtt[Silicon->GetName()] = new G4VisAttributes(G4Colour::Green());
    mVisAtt[Aluminum->GetName()] = new G4VisAttributes(G4Colour::Grey());
    mVisAtt[NemaG10->GetName()] = new G4VisAttributes(G4Colour::Magenta());
    mVisAtt[ArCO2->GetName()] = new G4VisAttributes(G4Colour::Yellow());
    mVisAtt[EJ204->GetName()] = new G4VisAttributes(G4Colour::Green());
    mVisAtt[Torlon->GetName()] = new G4VisAttributes(G4Colour::Grey());
    mVisAtt[Tungsten->GetName()] = new G4VisAttributes(G4Colour::Black());
    mVisAtt[PbWO4->GetName()] = new G4VisAttributes(G4Colour::Blue());
    mVisAtt[PbGlass->GetName()] = new G4VisAttributes(G4Colour::Blue());
    mVisAtt[Galaxy->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    G4Material *VacuumMaterial = Vacuum;
    G4Material *TargetMaterial = H2Gas;
    G4Material *TarCelMaterial = Copper;
    G4Material *WindowsMaterial = Kapton;
    G4Material *RecoilDetMaterial = Silicon;
    G4Material *ChamberMaterial = Aluminum;
    G4Material *GEMFrameMaterial = NemaG10;
    G4Material *GEMFoilMaterial = Kapton;
    G4Material *GEMGasMaterial = ArCO2;
    G4Material *SciPlaneMaterial = EJ204;
    G4Material *HyCalBoxMaterial = Torlon;
    G4Material *CollimatorMaterial = Tungsten;
    G4Material *CenterHyCalMaterial = PbWO4;
    //G4Material *OuterHyCalMaterial = PbGlass;
    G4Material *defaultMaterial = Galaxy;

    if (fConfig == "drad")
        TargetMaterial = D2Gas;

    // Sensitive detector manager
    G4SDManager *SDman = G4SDManager::GetSDMpointer();

    // World
    G4double WorldSizeXY = 150.0 * cm;
    G4double WorldSizeZ = 400.0 * cm;
    G4VSolid *solidWorld = new G4Box("Solid World", WorldSizeXY, WorldSizeXY, WorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, defaultMaterial, "Logical World");
    physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    if (fConfig == "prad") {
        // Target
        G4double TargetCenter = -300.0 * cm + 88.9 * mm;
        G4VSolid *solidTargetCon = new G4Box("Solid Target Container", 3.5 * cm, 3.5 * cm, 2.1 * cm);
        G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, defaultMaterial, "Logical Target Container");
        new G4PVPlacement(0, G4ThreeVector(0, 0, TargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

        // Target material
        G4double TargetR = 25.0 * mm;
        G4double TargetHalfL = 20.0 * mm;
        G4VSolid *solidTarget = new G4Tubs("Solid Target Material", 0, TargetR, TargetHalfL, 0, twopi);
        G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetMaterial, "Logical Target Material");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

        // Target cell
        G4double CellXY = 3.5 * cm;
        G4Box *CellBox = new G4Box("Cell Box", CellXY, CellXY, TargetHalfL);
        G4Tubs *CellTube = new G4Tubs("Cell Tube", 0, TargetR, TargetHalfL + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidCell = new G4SubtractionSolid("Solid Target Cell", CellBox, CellTube);
        G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TarCelMaterial, "Logical Target Cell");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

        // Windows
        G4double CellApertureR = 2.0 * mm;
        G4double CellWinThickness = 7.5 * um;
        G4Box *CellWinBox = new G4Box("Cell Window Box", CellXY, CellXY, CellWinThickness / 2.0);
        G4Tubs *CellWinTube = new G4Tubs("Cell Window Tube", 0, CellApertureR, CellWinThickness + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidCellWin = new G4SubtractionSolid("Solid Target Cell", CellWinBox, CellWinTube);
        G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, WindowsMaterial, "Logical Target Window");
        new G4PVPlacement(0, G4ThreeVector(0, 0, -TargetHalfL - CellWinThickness / 2.0), logicCellWin, "Upstream Target Window", logicTargetCon, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0, 0, +TargetHalfL + CellWinThickness / 2.0), logicCellWin, "Downstream Target Window", logicTargetCon, false, 0);
    } else if (fConfig == "drad") {
        // Target
        G4VSolid *solidTargetCon = new G4Tubs("Solid Target Container", 0, fTargetR + 0.1 * cm, fTargetHalfL + 0.1 * cm, 0, twopi);
        G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, defaultMaterial, "Logical Target Container");
        new G4PVPlacement(0, G4ThreeVector(0, 0, fTargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

        // Target material
        G4VSolid *solidTarget = new G4Tubs("Solid Target Material", 0, fTargetR, fTargetHalfL, 0, twopi);
        G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetMaterial, "Logical Target Material");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

        // Target cell
        G4double CellThickness = 0.5 * mm;
        G4VSolid *solidCell = new G4Tubs("Solid Target Cell", fTargetR, fTargetR + CellThickness, fTargetHalfL, 0, twopi);
        G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TarCelMaterial, "Logical Target Cell");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

        // Windows
        G4double CellApertureR = 2.0 * mm;
        G4double CellWinThickness = 7.5 * um;
        G4VSolid *solidCellWin = new G4Tubs("Solid Cell Window", CellApertureR, fTargetR + CellThickness, CellWinThickness / 2.0, 0, twopi);
        G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, WindowsMaterial, "Logical Target Window");
        new G4PVPlacement(0, G4ThreeVector(0, 0, -fTargetHalfL - CellWinThickness / 2.0), logicCellWin, "Upstream Target Window", logicTargetCon, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0, 0, +fTargetHalfL + CellWinThickness / 2.0), logicCellWin, "Downstream Target Window", logicTargetCon, false, 0);

        // Recoil detector
        G4double RecoilDetAng = twopi / fRecoilDetNSeg;
        G4double RecoilDetOR = fTargetR * cos(RecoilDetAng / 2.0) - 0.5 * mm;
        G4double RecoilDetIR = RecoilDetOR - fRecoilDetThickness;
        G4double rInnerRD[] = {RecoilDetIR, RecoilDetIR};
        G4double rOuterRD[] = {RecoilDetOR, RecoilDetOR};
        G4double zPlaneRD[] = {-fRecoilDetHalfL, fRecoilDetHalfL};
        G4VSolid *solidRecoilDet = new G4Polyhedra("Solid Recoil Detector", 0, twopi, fRecoilDetNSeg, 2, zPlaneRD, rInnerRD, rOuterRD);
        G4LogicalVolume *logicRecoilDet = new G4LogicalVolume(solidRecoilDet, RecoilDetMaterial, "Logical Recoil Detector");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicRecoilDet, "Recoil Detector", logicTarget, false, 0);
        StandardDetectorSD *RecoilDetSD = new StandardDetectorSD("pradsim/RecoilDetector", "RD");
        SDman->AddNewDetector(RecoilDetSD);
        logicRecoilDet->SetSensitiveDetector(RecoilDetSD);
    }

    if (fConfig == "prad") {
        // Target chamber
        // For now, only built the downstream chamber with window
        // The downstream chamber window should locate at -3000.0 + 88.9 + 74.0  = -2837.1 mm
        // The length of the downstream chamber is 381.7 mm
        // The total length of the downstream chamber and the tube in total is 710.0 mm
        // Here the downstream chamber and the tube are built together to be the now down stream chamber.
        // So the center of this geometry should be at -2837.1 + 710.0 / 2 = -2482.1 mm
        G4double DownChamberCenter = -248.21 * cm;
        G4double DownChamberHalfL = 71.0 / 2.0 * cm;
        G4double DownChamberUR = 8.00 * cm;

        // Downstream chamber
        G4double rInnerDC[] = {7.56 * cm, 7.56 * cm, 7.56 * cm, 7.56 * cm, 17.30 * cm, 17.30 * cm};
        G4double rOuterDC[] = {8.00 * cm, 8.00 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm};
        G4double zPlaneDC[] = {0, 32.83 * cm, 32.83 * cm, 35.37 * cm, 35.37 * cm, 71.00 * cm};
        G4VSolid *solidDownChamber = new G4Polycone("Solid Downstream Chamber", 0, twopi, 6, zPlaneDC, rInnerDC, rOuterDC);
        G4LogicalVolume *logicDownChamber = new G4LogicalVolume(solidDownChamber, ChamberMaterial, "Logical Downstream Chamber");
        new G4PVPlacement(0, G4ThreeVector(0, 0, DownChamberCenter - DownChamberHalfL), logicDownChamber, "Downstream Chamber", logicWorld, false, 0);

        // Windows
        G4double DownChamberApertureR = 22.8 * mm;
        G4double DownChamberWinThickness = 7.5 * um;
        G4Tubs *solidDownChamberWin = new G4Tubs("Solid Downstream Chamber Window", DownChamberApertureR, DownChamberUR, DownChamberWinThickness / 2.0, 0, twopi);
        G4LogicalVolume *logicDownChamberWin = new G4LogicalVolume(solidDownChamberWin, WindowsMaterial, "Logical Downstream Chamber Window");
        new G4PVPlacement(0, G4ThreeVector(0, 0, DownChamberCenter - DownChamberHalfL - DownChamberWinThickness / 2.0), logicDownChamberWin, "Downstream Chamber Window", logicWorld, false, 0);

        // Vacuum box
        // The length of the vacuum box is 4250.0 mm
        // So the center of this geometry should be at -3000.0 + 88.9 + 74.0 + 710.0 + 2125.0 = -2.1 mm
        G4double VacBoxCenter = -0.21 * cm; // !!! defined twice below for vacuum tube
        G4double VacBoxHalfL = 212.5 * cm; // !!! defined twice below for vacuum tube
        G4double VacBoxMaxR = 78.11 * cm;
        G4double rInner2[] = {17.30 * cm, 17.30 * cm, 50.17 * cm, 50.17 * cm, 78.11 * cm, 78.11 * cm};
        G4double rOuter2[] = {17.78 * cm, 17.78 * cm, 50.80 * cm, 50.80 * cm, 78.74 * cm, 78.74 * cm};
        G4double zPlane2[] = {0, 6.8 * cm, 17.6 * cm, 215.3 * cm, 229.5 * cm, 425.00 * cm};
        G4VSolid *solidVacBox = new G4Polycone("Solid Vacuum Box", 0, twopi, 6, zPlane2, rInner2, rOuter2);
        G4LogicalVolume *logicVacBox = new G4LogicalVolume(solidVacBox, ChamberMaterial, "Logical Vacuum Box");
        new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter - VacBoxHalfL), logicVacBox, "Vacuum Box", logicWorld, false, 0);

        // Window and flange
        G4double ArcDistance = 5.59 * cm; // !!! defined twice below for vacuum tube
        G4double ArcEndR = (ArcDistance * ArcDistance + VacBoxMaxR * VacBoxMaxR) / (2 * ArcDistance);
        G4double ArcEndThickness = 1.6 * mm;
        G4double FlangeIR = 1.9 * cm;
        G4double FlangeOR = 3.0 * cm;
        G4double FlangeHalfL = 0.5 * cm;
        G4Sphere *VacSphere = new G4Sphere("Vac Sphere", ArcEndR - ArcEndThickness, ArcEndR, 0, twopi, pi - asin(VacBoxMaxR / ArcEndR), pi);
        G4Tubs *VacHole = new G4Tubs("Vac Hole", 0, FlangeOR, ArcEndR + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidVacBoxWin = new G4SubtractionSolid("Solid Vacuum Box Window", VacSphere, VacHole);
        G4LogicalVolume *logicVacBoxWin = new G4LogicalVolume(solidVacBoxWin, ChamberMaterial, "Logical Vacuum Box Window");
        new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter + VacBoxHalfL + ArcEndR - ArcDistance), logicVacBoxWin, "Vacuum Box Window", logicWorld, false, 0);
        G4VSolid *solidFlange = new G4Tubs("Solid Vacuum Flange", FlangeIR, FlangeOR, FlangeHalfL, 0, twopi);
        G4LogicalVolume *logicFlange = new G4LogicalVolume(solidFlange, ChamberMaterial, "Logical Vacuum Flange");
        new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter + VacBoxHalfL - ArcDistance + FlangeHalfL), logicFlange, "Vacuum Flange", logicWorld, false, 0);

        // Vacuum Tube
        G4double VacTubeIR = 1.8 * cm;
        G4double VacTubeOR = 1.9 * cm;
        G4double VacTubeL = WorldSizeZ - 10.0 * cm - VacBoxCenter - VacBoxHalfL + ArcDistance;
        G4VSolid *solidVacTube = new G4Tubs("Solid Vacuum Tube", VacTubeIR, VacTubeOR, VacTubeL / 2.0, 0, twopi);
        G4LogicalVolume *logicVacTube = new G4LogicalVolume(solidVacTube, ChamberMaterial, "Logical Vacuum Tube");
        new G4PVPlacement(0, G4ThreeVector(0, 0, WorldSizeZ - 10.0 * cm - VacTubeL / 2.0), logicVacTube, "Vacuum Tube", logicWorld, false, 0);
    } else if (fConfig == "drad") {
        // no chamber for drad configuration now
    }

    // GEM
    // Center of two GEM should be at -3000.0 + 88.9 + (5222.0 + 5183.0) / 2 = 2291.4 mm // (5222.0 + 5183.0) / 2 from Weizhi
    G4double GEMCenter = 229.14 * cm;
    G4double GEMGap = 4.0 * cm; // Gap between two GEM
    G4RotationMatrix rmGEM2;
    rmGEM2.rotateZ(180.0 * deg);
    G4Box *GEMOuterFrame = new G4Box("GEM OuterFrame", 332.5 * mm, 699.9 * mm, 6.0 * mm);
    G4Tubs *GEMPipeHole = new G4Tubs("GEM Pipe Hole", 0, 22.0 * mm, 6.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMCon = new G4SubtractionSolid("Solid GEM Container", GEMOuterFrame, GEMPipeHole, 0, G4ThreeVector(-253.0 * mm, 0, 0));
    G4LogicalVolume *logicGEMCon = new G4LogicalVolume(solidGEMCon, defaultMaterial, "Logical GEM Container");

    if (fConfig == "prad") {
        new G4PVPlacement(0, G4ThreeVector(25.3 * cm, 0, GEMCenter - GEMGap / 2.0), logicGEMCon, "GEM Container", logicWorld, false, 0);
        new G4PVPlacement(G4Transform3D(rmGEM2, G4ThreeVector(-25.3 * cm, 0, GEMCenter + GEMGap / 2.0)), logicGEMCon, "GEM Container", logicWorld, false, 0);
    } else if (fConfig == "drad") {
        new G4PVPlacement(0, G4ThreeVector(25.3 * cm, 0, fGEM1Center - GEMGap / 2.0), logicGEMCon, "GEM Container", logicWorld, false, 0);
        new G4PVPlacement(G4Transform3D(rmGEM2, G4ThreeVector(-25.3 * cm, 0, fGEM1Center + GEMGap / 2.0)), logicGEMCon, "GEM Container", logicWorld, false, 0);
        new G4PVPlacement(0, G4ThreeVector(25.3 * cm, 0, fGEM2Center - GEMGap / 2.0), logicGEMCon, "GEM Container", logicWorld, false, 1);
        new G4PVPlacement(G4Transform3D(rmGEM2, G4ThreeVector(-25.3 * cm, 0, fGEM2Center + GEMGap / 2.0)), logicGEMCon, "GEM Container", logicWorld, false, 1);
    }

    // GEM Gas
    G4Box *GEMGasPiece1 = new G4Box("GEM Gas Piece 1", 275.0 * mm, 674.4 * mm, 6.0 * mm);
    G4Box *GEMGasPiece2 = new G4Box("GEM Gas Piece 2", 37.0 * mm, 29.5 * mm, 6.2 * mm);
    G4SubtractionSolid *solidGEMGas = new G4SubtractionSolid("Solid GEM Gas", GEMGasPiece1, GEMGasPiece2, 0, G4ThreeVector(-245.5 * mm, 0, 0));
    G4LogicalVolume *logicGEMGas = new G4LogicalVolume(solidGEMGas, GEMGasMaterial, "Logical GEM Gas");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMGas, "GEM Gas", logicGEMCon, false, 0);

    // GEM Frame
    G4Box *GEMFramePiece1 = new G4Box("GEM Frame Piece 1", 275.0 * mm, 674.4 * mm, 6.1 * mm);
    G4SubtractionSolid *GEMFramePiece2 = new G4SubtractionSolid("GEM Frame Piece 2", GEMFramePiece1, GEMGasPiece2, 0, G4ThreeVector(-245.5 * mm, 0, 0));
    G4SubtractionSolid *GEMFrameNoHole = new G4SubtractionSolid("GEM Frame No Hole", GEMOuterFrame, GEMFramePiece2);
    G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid("Solid GEM Frame", GEMFrameNoHole, GEMPipeHole, 0, G4ThreeVector(-253.0 * mm, 0, 0));
    G4LogicalVolume *logicGEMFrame = new G4LogicalVolume(solidGEMFrame, GEMFrameMaterial, "Logical GEM Frame");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMFrame, "GEM Frame", logicGEMCon, false, 0);

    // GEM Foil
    G4double GEMFoilThickness = 55.0 * um;
    G4double GEMCoverThickness = 25.0 * um;
    G4Box *GEMCoverBox = new G4Box("GEM Cover Box", 275.0 * mm, 674.4 * mm, GEMCoverThickness / 2.0);
    G4SubtractionSolid *solidGEMCover = new G4SubtractionSolid("Solid GEM Cover", GEMCoverBox, GEMGasPiece2, 0, G4ThreeVector(-245.5 * mm, 0, 0));
    G4LogicalVolume *logicGEMCover = new G4LogicalVolume(solidGEMCover, GEMFoilMaterial, "Logical GEM Cover");
    G4Box *GEMFoilBox = new G4Box("GEM Foil Box", 275.0 * mm, 674.4 * mm, GEMFoilThickness / 2.0);
    G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid("Solid GEM Foil", GEMFoilBox, GEMGasPiece2, 0, G4ThreeVector(-245.5 * mm, 0, 0));
    G4LogicalVolume *logicGEMFoil = new G4LogicalVolume(solidGEMFoil, GEMFoilMaterial, "Logical GEM Foil");
    G4LogicalVolume *logicGEMReadout = new G4LogicalVolume(solidGEMFoil, GEMFoilMaterial, "Logical GEM Readout");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -6.0 * mm + GEMCoverThickness / 2.0), logicGEMCover, "GEM Cover Foil", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -3.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMFoil, "GEM Foil", logicGEMGas, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 2.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 2);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 4.0 * mm), logicGEMReadout, "GEM Readout Foil", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 6.0 * mm - GEMCoverThickness / 2.0), logicGEMCover, "GEM CoverFoil", logicGEMGas, false, 1);
    StandardDetectorSD *GEMSD = new StandardDetectorSD("pradsim/GasElectronMultiplier", "GEM");
    SDman->AddNewDetector(GEMSD);
    logicGEMReadout->SetSensitiveDetector(GEMSD);

    if (fConfig == "drad") {
        // Scintillator Plane
        G4double SciPlaneThickness = 2.0 * mm;
        G4double SciPlaneHalfX = 75.0 * cm;
        G4double SciPlaneHalfY = 75.0 * cm;

        G4VSolid *solidSciPlane = new G4Box("Solid Scintillator Plane", SciPlaneHalfX, SciPlaneHalfY, SciPlaneThickness / 2.0);
        G4LogicalVolume *logicSciPlane = new G4LogicalVolume(solidSciPlane, SciPlaneMaterial, "Logical Scintillator Plane");
        new G4PVPlacement(0, G4ThreeVector(0, 0, fSciPlaneCenter), logicSciPlane, "Scintillator Plane", logicWorld, false, 0);
        StandardDetectorSD *SciPlaneSD = new StandardDetectorSD("pradsim/ScintillatorPlane", "SP");
        SDman->AddNewDetector(SciPlaneSD);
        logicSciPlane->SetSensitiveDetector(SciPlaneSD);
    }

    // HyCal
    // The crystal surface should be at -3000.0 + 88.9 + 5640.0 = 2728.9 mm // 5640.0 from Weizhi
    G4double PbGlassL = 45.0 * cm;
    //G4double CrystalL = 18.0 * cm;
    G4double CrystalDiffL = 10.12 * cm;
    G4double CrystalSurf = 272.89 * cm; // Surface of the PWO

    if (fConfig == "prad")
        fCrystalSurf = CrystalSurf;

    G4double HyCalCenter = fCrystalSurf - CrystalDiffL + PbGlassL / 2.0;

    // HyCal box
    G4double HyCalBoxCenter = HyCalCenter - 9.0 * cm + 30.0 * cm; // Check
    G4Box *HyCalBoxOuter = new G4Box("HyCal Box Outer", 70.0 * cm, 70.0 * cm, 60.0 * cm);
    G4Box *HyCalBoxInner = new G4Box("HyCal Box Inner", 66.0 * cm, 66.0 * cm, 59.6 * cm);
    G4SubtractionSolid *HyCalBoxNoHole = new G4SubtractionSolid("HyCal Box No Hole", HyCalBoxOuter, HyCalBoxInner);
    G4Tubs *HyCalBoxHole = new G4Tubs("HyCal Box Hole", 0, 25.0 * mm, 60.5 * cm, 0, twopi);
    G4SubtractionSolid *solidHyCalBox = new G4SubtractionSolid("Solid HyCal Box", HyCalBoxNoHole, HyCalBoxHole);
    G4LogicalVolume *logicHyCalBox = new G4LogicalVolume(solidHyCalBox, HyCalBoxMaterial, "Logical HyCal Box");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalBoxCenter), logicHyCalBox, "HyCal Box", logicWorld, false, 0);

    // HyCal container
    G4Box *HyCalConPiece1 = new G4Box("HyCal Container Piece 1", 58.21 * cm, 58.17 * cm, PbGlassL / 2.0);
    G4Box *HyCalConPiece2 = new G4Box("HyCal Container Piece 2", 35.30 * cm, 35.27 * cm, CrystalDiffL / 2.0 + 0.5 * mm);
    G4SubtractionSolid *HyCalConBox = new G4SubtractionSolid("HyCal Container Box", HyCalConPiece1, HyCalConPiece2, 0, G4ThreeVector(0, 0, (CrystalDiffL - PbGlassL) / 2.0 - 0.5 * mm));
    G4Box *HyCalConHole = new G4Box("HyCal Container Hole", 2.0 * cm, 2.0 * cm, PbGlassL / 2.0 + 1.0 * mm);
    G4SubtractionSolid *solidHyCalCon = new G4SubtractionSolid("Solid HyCal Container", HyCalConBox, HyCalConHole);
    G4LogicalVolume *logicHyCalCon = new G4LogicalVolume(solidHyCalCon, defaultMaterial, "Logical HyCal Container");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter), logicHyCalCon, "HyCal Container", logicWorld, false, 0);

    // HyCal modules
    G4VSolid *solidAbsorber = new G4Box("Solid Crystal Block", 1.025 * cm, 1.025 * cm, 90.0 * mm);
    G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber, CenterHyCalMaterial, "Logical Crystal Block");
    HyCalParameterisation *param = new HyCalParameterisation("config/hycal.conf");
    new G4PVParameterised("HyCal Crystal", logicAbsorber, logicHyCalCon, kUndefined, param->GetNumber(), param, false);
    StandardDetectorSD *HyCalSD = new StandardDetectorSD("pradsim/HybridCalorimeter", "HC");
    SDman->AddNewDetector(HyCalSD);
    logicAbsorber->SetSensitiveDetector(HyCalSD);

    // Collimators around the central hole
    G4VSolid *CollConBox = new G4Box("HyCal Collimator Container Box", 4.1 * cm, 4.1 * cm, 5.0 * cm);
    G4SubtractionSolid *solidCollCon = new G4SubtractionSolid("Solid HyCal Collimator Container", CollConBox, HyCalConHole);
    G4LogicalVolume *logicCollCon = new G4LogicalVolume(solidCollCon, defaultMaterial, "Logical HyCal Collimator Container");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter - PbGlassL / 2.0 + CrystalDiffL - 5.1 * cm), logicCollCon, "HyCal Collimator Container", logicWorld, false, 0);
    G4VSolid *solidColl = new G4Box("Solid HyCal Collimator", 1.025 * cm, 1.025 * cm, 5.0 * cm);
    G4LogicalVolume *logicColl = new G4LogicalVolume(solidColl, CollimatorMaterial, "Logical HyCal Collimator");
    double pos_x[12] = { -3.075, -1.025, 1.025, 3.075, -3.075, 3.075, -3.075, 3.075, -3.075, -1.025, 1.025, 3.075};
    double pos_y[12] = { -3.075, -3.075, -3.075, -3.075, -1.025, -1.025, 1.025, 1.025, 3.075, 3.075, 3.075, 3.075};

    for (int i = 0; i < 12; ++i)
        new G4PVPlacement(0, G4ThreeVector(pos_x[i] * cm, pos_y[i] * cm, 0), logicColl, "HyCal Collimator", logicCollCon, false, i);

    // Virtual Detector
    G4RotationMatrix rmVD;
    rmVD.rotateZ(90.0 * deg);
    G4Box *VDPiece1 = new G4Box("Virtual Detector Piece 1", 58.21 * cm, 58.17 * cm, 0.5 * mm);
    G4Box *VDPiece2 = new G4Box("Virtual Detector Piece 2", 35.30 * cm, 35.27 * cm, 0.6 * mm);
    G4SubtractionSolid *VDPiece3 = new G4SubtractionSolid("Virtual Detector Piece 3", VDPiece1, VDPiece2);
    G4Box *VDPiece4 = new G4Box("Virtual Detector Piece 4", 1.97 * cm / 2.0, 22.86 * cm / 2.0, 0.6 * mm);
    G4SubtractionSolid *VDPiece5 = new G4SubtractionSolid("Virtual Detector Piece 5", VDPiece3, VDPiece4, 0, G4ThreeVector((58.21 - 1.97 / 2.0 + 0.001) * cm, (58.17 - 22.86 / 2.0 + 0.001) * cm, 0));
    G4SubtractionSolid *VDPiece6 = new G4SubtractionSolid("Virtual Detector Piece 6", VDPiece5, VDPiece4, G4Transform3D(rmVD, G4ThreeVector(-(58.21 - 22.86 / 2.0 + 0.001) * cm, (58.17 - 1.97 / 2.0 + 0.001) * cm, 0)));
    G4SubtractionSolid *VDPiece7 = new G4SubtractionSolid("Virtual Detector Piece 7", VDPiece6, VDPiece4, 0, G4ThreeVector(-(58.21 - 1.97 / 2.0 + 0.001) * cm, -(58.17 - 22.86 / 2.0 + 0.001) * cm, 0));
    G4SubtractionSolid *solidVD1 = new G4SubtractionSolid("Solid Virtual Detector 1", VDPiece7, VDPiece4, G4Transform3D(rmVD, G4ThreeVector((58.21 - 22.86 / 2.0 + 0.001) * cm, -(58.17 - 1.97 / 2.0 + 0.001) * cm, 0)));
    G4Box *VDHole = new G4Box("Virtual Detector Hole", 2.075 * cm, 2.075 * cm, 0.6 * mm);
    G4Box *VDPiece10 = new G4Box("Virtual Detector Piece 10", 35.30 * cm, 35.27 * cm, 0.5 * mm);
    G4SubtractionSolid *solidVD2 = new G4SubtractionSolid("Solid Virtual Detector 2", VDPiece10, VDHole);
    G4LogicalVolume *logicVD1 = new G4LogicalVolume(solidVD1, VacuumMaterial, "Logical Virtual Detector 1");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter - PbGlassL / 2.0 - 0.5 * mm), logicVD1, "Virtual Detector", logicWorld, false, 0);
    G4LogicalVolume *logicVD2 = new G4LogicalVolume(solidVD2, VacuumMaterial, "Logical Virtual Detector 2");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter - PbGlassL / 2.0 + CrystalDiffL - 0.5 * mm), logicVD2, "Virtual Detector", logicWorld, false, 0);
    StandardDetectorSD *VirtualSD = new StandardDetectorSD("pradsim/VirtualDetector", "VD");
    SDman->AddNewDetector(VirtualSD);
    logicVD1->SetSensitiveDetector(VirtualSD);
    logicVD2->SetSensitiveDetector(VirtualSD);

    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(mVisAtt[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    // Always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
