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

#include "CalorimeterSD.hh"
#include "DetectorMessenger.hh"
#include "HyCalParameterisation.hh"
#include "StandardDetectorSD.hh"
#include "TrackingDetectorSD.hh"

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
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VUserDetectorConstruction.hh"

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

    fVisAtts.clear();

    fTargetCenter = -300.0 * cm;
    fTargetR = 10.0 * cm;
    fTargetHalfL = 3.0 * cm;
    fTargetMat = "D2Gas";

    fRecoilDetNSeg = 72;
    fRecoilDetCenter = -300.0 * cm;
    fRecoilDetHalfL = 2.0 * cm;
    fRecoilDetThickness = 2.0 * mm;

    fGEM1Center = 160.0 * cm;
    fGEM2Center = 200.0 * cm;

    fSciPlaneCenter = 240.0 * cm;

    fCrystalSurf = 270.0 * cm;

    if (fConfig == "drad") {
        fRecoilDetSDOn = true;
        fGEMSDOn = true;
        fSciPlaneSDOn = true;
        fHyCalSDOn = 1;
    } else {
        fRecoilDetSDOn = false;
        fGEMSDOn = true;
        fSciPlaneSDOn = false;
        fHyCalSDOn = 1;
    }

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
    // Define materials
    DefineMaterials();

    // Define volumes
    if (fConfig == "drad")
        return DefineDRadVolumes();
    else
        return DefinePRadVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if (fConfig == "drad")
        DefineDRadSDs();
    else
        DefinePRadSDs();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4String symbol;
    G4int z, n;
    G4double a;
    G4double density;
    G4int ncomponents, natoms;
    G4double fractionmass;

    G4NistManager *pNM = G4NistManager::Instance();

    // Define elements from NIST material table
    G4Element *H  = pNM->FindOrBuildElement(z = 1);
    G4Element *He = pNM->FindOrBuildElement(z = 2);
    G4Element *C  = pNM->FindOrBuildElement(z = 6);
    G4Element *N  = pNM->FindOrBuildElement(z = 7);
    G4Element *O  = pNM->FindOrBuildElement(z = 8);
    G4Element *F  = pNM->FindOrBuildElement(z = 9);
    G4Element *Al = pNM->FindOrBuildElement(z = 13);
    G4Element *Si = pNM->FindOrBuildElement(z = 14);
    G4Element *P  = pNM->FindOrBuildElement(z = 15);
    G4Element *S  = pNM->FindOrBuildElement(z = 16);
    G4Element *Ar = pNM->FindOrBuildElement(z = 18);
    G4Element *Cr = pNM->FindOrBuildElement(z = 24);
    G4Element *Mn = pNM->FindOrBuildElement(z = 25);
    G4Element *Fe = pNM->FindOrBuildElement(z = 26);
    G4Element *Ni = pNM->FindOrBuildElement(z = 28);
    G4Element *Cu = pNM->FindOrBuildElement(z = 29);
    G4Element *W  = pNM->FindOrBuildElement(z = 74);
    G4Element *Pb = pNM->FindOrBuildElement(z = 86);

    // Define isotopes
    G4Isotope *H2 = new G4Isotope("H2", z = 1, n = 2, a = 2.0141 * g / mole);
    G4Element *D = new G4Element("Deuterium", symbol = "D", ncomponents = 1);
    D->AddIsotope(H2, 1.0);

    // Define materials
    // Space Vacuum
    G4Material *Galaxy = new G4Material("Galaxy", density = universe_mean_density, ncomponents = 1, kStateGas, 0.1 * kelvin, 1.0e-19 * pascal);
    Galaxy->AddElement(H, fractionmass = 1.0);
    fVisAtts[Galaxy->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Air
    G4Material *Air = new G4Material("Air", density = 1.292 * mg / cm3, ncomponents = 2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);
    fVisAtts[Air->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Air vacuum of 1.e-6 torr at room temperature, 1 atmosphere = 760 torr
    G4Material *Vacuum = new G4Material("Vacuum", density = 1.0e-6 / 760.0 * 1.292 * mg / cm3, ncomponents = 1, kStateGas, STP_Temperature, 1.0e-6 / 760.0 * atmosphere);
    Vacuum->AddMaterial(Air, fractionmass = 1.0);
    fVisAtts[Vacuum->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Hydrogen Gas (T = 19.5 K, P = 470 mTorr)
    G4Material *H2Gas = new G4Material("H2Gas", density = 0.47 / 760.0 * 273.15 / 19.5 * 0.08988 * mg / cm3, ncomponents = 1, kStateGas, 19.5 * kelvin, 0.47 / 760.0 * atmosphere);
    H2Gas->AddElement(H, natoms = 2);
    fVisAtts[H2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());

    // Deuteron Gas
    G4Material *D2Gas = new G4Material("D2Gas", density = 0.47 / 760.0 * 273.15 / 19.5 * 0.1796 * mg / cm3, ncomponents = 1, kStateGas, 19.5 * kelvin, 0.47 / 760.0 * atmosphere);
    D2Gas->AddElement(D, natoms = 2);
    fVisAtts[D2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());

    // Copper C101
    G4Material *Copper = new G4Material("Copper", density = 8.92 * g / cm3, ncomponents = 1);
    Copper->AddElement(Cu, natoms = 1);
    fVisAtts[Copper->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Kapton
    G4Material *Kapton = new G4Material("Kapton", density = 1.42 * g / cm3, ncomponents = 4);
    Kapton->AddElement(H, fractionmass = 0.0273);
    Kapton->AddElement(C, fractionmass = 0.7213);
    Kapton->AddElement(N, fractionmass = 0.0765);
    Kapton->AddElement(O, fractionmass = 0.1749);
    fVisAtts[Kapton->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Si
    G4Material *Silicon = new G4Material("Silicon", density = 2.329 * g / cm3, ncomponents = 1);
    Silicon->AddElement(Si, natoms = 1);
    fVisAtts[Silicon->GetName()] = new G4VisAttributes(G4Colour::Green());

    // Aluminum
    G4Material *Aluminum = new G4Material("Aluminum", density = 2.700 * g / cm3, ncomponents = 1);
    Aluminum->AddElement(Al, natoms = 1);
    fVisAtts[Aluminum->GetName()] = new G4VisAttributes(G4Colour::Grey());

    // Tedlar
    G4Material *Tedlar = new G4Material("Tedlar", density = 1.545 * g / cm3, ncomponents = 3);
    Tedlar->AddElement(H, natoms = 3);
    Tedlar->AddElement(C, natoms = 2);
    Tedlar->AddElement(F, natoms = 1);

    // Stainless Steel
    G4Material *SSteel = new G4Material("SSteel", density = 7.9 * g / cm3, ncomponents = 9);
    SSteel->AddElement(C, fractionmass = 0.0007);
    SSteel->AddElement(Si, fractionmass = 0.01);
    SSteel->AddElement(Mn, fractionmass = 0.02);
    SSteel->AddElement(Ni, fractionmass = 0.09);
    SSteel->AddElement(P, fractionmass = 0.00045);
    SSteel->AddElement(S, fractionmass = 0.00015);
    SSteel->AddElement(Cr, fractionmass = 0.18);
    SSteel->AddElement(N, fractionmass = 0.0011);
    SSteel->AddElement(Fe, fractionmass = 0.6976);
    fVisAtts[SSteel->GetName()] = new G4VisAttributes(G4Colour::Grey());

    // GEM Frame G10
    G4Material *NemaG10 = new G4Material("NemaG10", density = 1.700 * g / cm3, ncomponents = 4);
    NemaG10->AddElement(Si, natoms = 1);
    NemaG10->AddElement(O, natoms = 2);
    NemaG10->AddElement(C, natoms = 3);
    NemaG10->AddElement(H, natoms = 3);
    fVisAtts[NemaG10->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Ar/CO2 Gas
    G4Material *CO2 = new G4Material("CO2", density = 1.842e-3 * g / cm3, ncomponents = 2);
    CO2->AddElement(C, natoms = 1);
    CO2->AddElement(O, natoms = 2);
    G4Material *ArCO2 = new G4Material("ArCO2", density = 1.715e-3 * g / cm3, ncomponents = 2);
    ArCO2->AddElement(Ar, fractionmass = 0.7);
    ArCO2->AddMaterial(CO2, fractionmass = 0.3);
    fVisAtts[ArCO2->GetName()] = new G4VisAttributes(G4Colour::Yellow());

    // He Gas
    G4Material *HeGas = new G4Material("HeGas", density = 0.1786e-3 * g / cm3, ncomponents = 1);
    HeGas->AddElement(He, natoms = 1);
    fVisAtts[HeGas->GetName()] = new G4VisAttributes(G4Colour::Cyan());

    // Scintillator EJ204
    G4Material *EJ204 = new G4Material("EJ204", density = 1.032 * g / cm3, ncomponents = 2);
    EJ204->AddElement(H, natoms = 521);
    EJ204->AddElement(C, natoms = 474);
    fVisAtts[EJ204->GetName()] = new G4VisAttributes(G4Colour::Green());

    // Torlon4203L
    G4Material *Torlon = new G4Material("Torlon", density = 1.412 * g / cm3, ncomponents = 5);
    Torlon->AddElement(C, natoms = 9);
    Torlon->AddElement(H, natoms = 4);
    Torlon->AddElement(N, natoms = 2);
    Torlon->AddElement(O, natoms = 3);
    Torlon->AddElement(Ar, natoms = 1);
    fVisAtts[Torlon->GetName()] = new G4VisAttributes(G4Colour::Grey());

    // Tungsten
    G4Material *Tungsten = new G4Material("Tungsten", density = 19.25 * g / cm3, ncomponents = 1);
    Tungsten->AddElement(W, natoms = 1);
    fVisAtts[Tungsten->GetName()] = new G4VisAttributes(G4Colour::Black());

    // PbWO4 Crystal
    G4Material *PbWO4 = new G4Material("PbWO4", density = 8.300 * g / cm3, ncomponents = 3);
    PbWO4->AddElement(Pb, natoms = 1);
    PbWO4->AddElement(W, natoms = 1);
    PbWO4->AddElement(O, natoms = 4);
    fVisAtts[PbWO4->GetName()] = new G4VisAttributes(G4Colour::Blue());

    // Lead Glass
    G4Material *SiO2 = new G4Material("SiO2", density = 2.200 * g / cm3, ncomponents = 2);
    SiO2->AddElement(Si, natoms = 1);
    SiO2->AddElement(O, natoms = 2);
    G4Material *PbGlass = new G4Material("PbGlass", density = 3.85 * g / cm3, ncomponents = 2);
    PbGlass->AddElement(Pb, fractionmass = 0.5316);
    PbGlass->AddMaterial(SiO2, fractionmass = 0.4684);
    fVisAtts[PbGlass->GetName()] = new G4VisAttributes(G4Colour::Blue());

    // Print out material table
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefinePRadVolumes()
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *TargetM = G4Material::GetMaterial("H2Gas");
    G4Material *TargetCellM = G4Material::GetMaterial("Copper");
    G4Material *TargetWindowM = G4Material::GetMaterial("Kapton");
    G4Material *ChamberM = G4Material::GetMaterial("Aluminum");
    G4Material *ChamberWindowM = G4Material::GetMaterial("Kapton");
    G4Material *VacuumBoxM = G4Material::GetMaterial("Aluminum");
    G4Material *VacuumTubeM = G4Material::GetMaterial("SSteel");
    G4Material *GEMFrameM = G4Material::GetMaterial("NemaG10");
    G4Material *GEMGasM = G4Material::GetMaterial("ArCO2");
    G4Material *GEMFoilM = G4Material::GetMaterial("Kapton");
    G4Material *GEMCuFoilM = G4Material::GetMaterial("Copper");
    G4Material *HyCalBoxM = G4Material::GetMaterial("Torlon");
    G4Material *CollimatorM = G4Material::GetMaterial("Tungsten");
    G4Material *HyCalModuleM = G4Material::GetMaterial("PbWO4");

    // World
    G4double WorldSizeXY = 150.0 * cm;
    G4double WorldSizeZ = 400.0 * cm;
    G4VSolid *solidWorld = new G4Box("WorldS", WorldSizeXY, WorldSizeXY, WorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, DefaultM, "WorldLV");
    G4VPhysicalVolume *physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    // Target
    G4double TargetCenter = -300.0 * cm + 88.9 * mm;

    // Target Container
    G4VSolid *solidTargetCon = new G4Box("TargetContainerS", 3.5 * cm, 3.5 * cm, 2.1 * cm);
    G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, DefaultM, "TargetContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, TargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

    // Target material
    G4double TargetR = 25.0 * mm;
    G4double TargetHalfL = 20.0 * mm;
    G4VSolid *solidTarget = new G4Tubs("TargetS", 0, TargetR, TargetHalfL, 0, twopi);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetM, "TargetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

    // Target cell
    G4double CellXY = 3.5 * cm;
    G4Box *CellBox = new G4Box("CellBox", CellXY, CellXY, TargetHalfL);
    G4Tubs *CellTube = new G4Tubs("CellTube", 0, TargetR, TargetHalfL + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCell = new G4SubtractionSolid("TargetCellS", CellBox, CellTube);
    G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

    // Target cell windows
    G4double CellApertureR = 2.0 * mm;
    G4double CellWinThickness = 7.5 * um;
    G4Box *CellWinBox = new G4Box("CellWinBox", CellXY, CellXY, CellWinThickness / 2.0);
    G4Tubs *CellWinTube = new G4Tubs("CellWinTube", 0, CellApertureR, CellWinThickness + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCellWin = new G4SubtractionSolid("TargetWindowS", CellWinBox, CellWinTube);
    G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, TargetWindowM, "TargetWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -TargetHalfL - CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +TargetHalfL + CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 1);

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
    G4VSolid *solidDownChamber = new G4Polycone("DownstreamChamberS", 0, twopi, 6, zPlaneDC, rInnerDC, rOuterDC);
    G4LogicalVolume *logicDownChamber = new G4LogicalVolume(solidDownChamber, ChamberM, "DownstreamChamberLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, DownChamberCenter - DownChamberHalfL), logicDownChamber, "Downstream Chamber", logicWorld, false, 0);

    // Downstream chamber window
    G4double DownChamberApertureR = 22.8 * mm;
    G4double DownChamberWinThickness = 7.5 * um;
    G4Tubs *solidDownChamberWin = new G4Tubs("DownstreamChamberWindowS", DownChamberApertureR, DownChamberUR, DownChamberWinThickness / 2.0, 0, twopi);
    G4LogicalVolume *logicDownChamberWin = new G4LogicalVolume(solidDownChamberWin, ChamberWindowM, "DownstreamChamberWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, DownChamberCenter - DownChamberHalfL - DownChamberWinThickness / 2.0), logicDownChamberWin, "Downstream Chamber Window", logicWorld, false, 0);

    // Vacuum box
    // The length of the vacuum box is 4250.0 mm
    // So the center of this geometry should be at -3000.0 + 88.9 + 74.0 + 710.0 + 2125.0 = -2.1 mm
    G4double VacBoxCenter = -0.21 * cm;
    G4double VacBoxHalfL = 212.5 * cm;
    G4double VacBoxMaxR = 78.11 * cm;
    G4double rInner2[] = {17.30 * cm, 17.30 * cm, 50.17 * cm, 50.17 * cm, 78.11 * cm, 78.11 * cm};
    G4double rOuter2[] = {17.78 * cm, 17.78 * cm, 50.80 * cm, 50.80 * cm, 78.74 * cm, 78.74 * cm};
    G4double zPlane2[] = {0, 6.8 * cm, 17.6 * cm, 215.3 * cm, 229.5 * cm, 425.00 * cm};
    G4VSolid *solidVacBox = new G4Polycone("VacuumBoxS", 0, twopi, 6, zPlane2, rInner2, rOuter2);
    G4LogicalVolume *logicVacBox = new G4LogicalVolume(solidVacBox, VacuumBoxM, "VacuumBoxLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter - VacBoxHalfL), logicVacBox, "Vacuum Box", logicWorld, false, 0);

    // Vacuum box window
    G4double ArcDistance = 5.59 * cm;
    G4double ArcEndR = (ArcDistance * ArcDistance + VacBoxMaxR * VacBoxMaxR) / (2 * ArcDistance);
    G4double ArcEndThickness = 1.6 * mm;
    G4double VacBoxWinApertureR = 3.0 * cm;
    G4VSolid *solidVacBoxWin = new G4Sphere("VacuumBoxWindowS", ArcEndR - ArcEndThickness, ArcEndR, 0, twopi, pi - asin(VacBoxMaxR / ArcEndR), asin(VacBoxMaxR / ArcEndR) - asin((VacBoxWinApertureR + 0.1 * mm) / ArcEndR));
    G4LogicalVolume *logicVacBoxWin = new G4LogicalVolume(solidVacBoxWin, VacuumBoxM, "VacuumBoxWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter + VacBoxHalfL + ArcEndR - ArcDistance), logicVacBoxWin, "Vacuum Box Window", logicWorld, false, 0);

    // Vacuum Tube
    G4double VacTubeOR = 1.9 * cm;
    G4double VacTubeIR = VacTubeOR - 0.12446 * cm; // 0.049 in = 0.12446 cm from Eugene
    G4double VacTubeL = WorldSizeZ - 10.0 * cm - VacBoxCenter - VacBoxHalfL + ArcDistance;
    G4VSolid *solidVacTube = new G4Tubs("VacuumTubeS", VacTubeIR, VacTubeOR, VacTubeL / 2.0, 0, twopi);
    G4LogicalVolume *logicVacTube = new G4LogicalVolume(solidVacTube, VacuumTubeM, "VacuumTubeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, WorldSizeZ - 10.0 * cm - VacTubeL / 2.0), logicVacTube, "Vacuum Tube", logicWorld, false, 0);

    // Flange on vacuum tube
    G4double FlangeOR = VacBoxWinApertureR;
    G4double FlangeIR = VacTubeOR;
    G4double FlangeHalfL = 0.5 * cm;
    G4VSolid *solidFlange = new G4Tubs("FlangeS", FlangeIR, FlangeOR, FlangeHalfL, 0, twopi);
    G4LogicalVolume *logicFlange = new G4LogicalVolume(solidFlange, VacuumTubeM, "FlangeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, VacBoxCenter + VacBoxHalfL - ArcDistance + FlangeHalfL), logicFlange, "Flange", logicWorld, false, 0);

    // GEM
    // Center of two GEM should be at -3000.0 + 88.9 + (5222.0 + 5183.0) / 2 = 2291.4 mm // (5222.0 + 5183.0) / 2 from Weizhi
    G4double GEMCenter = 229.14 * cm;
    G4double GEMGap = 4.0 * cm; // Gap between two GEM
    G4double GEMHalfX = 55.04 * cm / 2.0;
    G4double GEMHalfY = 122.88 * cm / 2.0;
    G4double GEMHalfThickness = 7.0 * mm;
    G4double GEMHoleR = 2.2 * cm;
    G4double GEMCenterHalfXY = 7.4 * cm / 2.0;
    G4double GEMFrameWidth = 1.5 * cm;
    G4double GEMCenterOffset = GEMHalfX + GEMFrameWidth - GEMCenterHalfXY;

    // GEM Container
    G4Box *GEMConBox = new G4Box("GEMConBox", 1.0 * m, 1.0 * m, (GEMGap + 2.0 * GEMHalfThickness + 1.0 * mm) / 2.0);
    G4Tubs *GEMConTube = new G4Tubs("GEMConTube", 0, GEMHoleR, (GEMGap + 2.0 * GEMHalfThickness + 1.0 * mm) / 2.0 + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMCon = new G4SubtractionSolid("GEMContainerS", GEMConBox, GEMConTube);
    G4LogicalVolume *logicGEMCon = new G4LogicalVolume(solidGEMCon, DefaultM, "GEMContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMCenter), logicGEMCon, "GEM Container", logicWorld, false, 0);

    // GEM
    G4Box *GEMBox = new G4Box("GEMBox", GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfThickness);
    G4Tubs *GEMTube = new G4Tubs("GEMTube", 0, GEMHoleR, GEMHalfThickness + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEM = new G4SubtractionSolid("GEMS", GEMBox, GEMTube, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEM = new G4LogicalVolume(solidGEM, DefaultM, "GEMLV");
    new G4PVPlacement(0, G4ThreeVector(GEMCenterOffset, 0, GEMGap / 2.0), logicGEM, "GEM L", logicGEMCon, false, 0);
    G4RotationMatrix rmGEM;
    rmGEM.rotateZ(180.0 * deg);
    new G4PVPlacement(G4Transform3D(rmGEM, G4ThreeVector(-GEMCenterOffset, 0, -GEMGap / 2.0)), logicGEM, "GEM R", logicGEMCon, false, 1);

    // GEM Gas
    G4Box *GEMGasBox1 = new G4Box("GEMGasBox1", GEMHalfX, GEMHalfY, GEMHalfThickness);
    G4Box *GEMGasBox2 = new G4Box("GEMGasBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMGas = new G4SubtractionSolid("GEMGasS", GEMGasBox1, GEMGasBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMGas = new G4LogicalVolume(solidGEMGas, GEMGasM, "GEMGasLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMGas, "GEM Gas", logicGEM, false, 0);

    // GEM Frame
    G4Box *GEMFrameBox1 = new G4Box("GEMFrameBox1", GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfThickness);
    G4Box *GEMFrameBox2 = new G4Box("GEMFrameBox2", GEMHalfX, GEMHalfY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid("GEMFrameS", GEMFrameBox1, GEMFrameBox2);
    G4LogicalVolume *logicGEMFrame = new G4LogicalVolume(solidGEMFrame, GEMFrameM, "GEMFrameLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMFrame, "GEM Frame", logicGEM, false, 0);
    G4Box *GEMPipeBox = new G4Box("GEMPipeBox", GEMCenterHalfXY - GEMFrameWidth / 2.0, GEMCenterHalfXY, GEMHalfThickness);
    G4Tubs *GEMPipeTube = new G4Tubs("GEMPipeTube", 0, GEMHoleR, GEMHalfThickness + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMPipe = new G4SubtractionSolid("GEMPipeS", GEMPipeBox, GEMPipeTube, 0, G4ThreeVector(-GEMFrameWidth / 2.0, 0, 0));
    G4LogicalVolume *logicGEMPipe = new G4LogicalVolume(solidGEMPipe, GEMFrameM, "GEMPipeLV");
    new G4PVPlacement(0, G4ThreeVector(-GEMCenterOffset + GEMFrameWidth / 2.0, 0, 0), logicGEMPipe, "GEM Pipe", logicGEM, false, 0);

    // GEM Foil
    G4double GEMFoilThickness = 50.0 * um;
    G4double GEMCuFoilThickness = 5.0 * um;
    G4double GEMWinThickness = 25.0 * um;
    G4Box *GEMWinBox1 = new G4Box("GEMWinBox1", GEMHalfX, GEMHalfY, GEMWinThickness / 2.0);
    G4Box *GEMWinBox2 = new G4Box("GEMWinBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMWin = new G4SubtractionSolid("GEMWinS", GEMWinBox1, GEMWinBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMWin = new G4LogicalVolume(solidGEMWin, GEMFoilM, "GEMWinLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -GEMHalfThickness + GEMWinThickness / 2.0), logicGEMWin, "GEM Window", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - GEMWinThickness / 2.0), logicGEMWin, "GEM Window", logicGEMGas, false, 1);

    G4Box *GEMFoilBox1 = new G4Box("GEMFoilBox1", GEMHalfX, GEMHalfY, GEMFoilThickness / 2.0);
    G4Box *GEMFoilBox2 = new G4Box("GEMFoilBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid("GEMFoilS", GEMFoilBox1, GEMFoilBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMFoil = new G4LogicalVolume(solidGEMFoil, GEMFoilM, "GEMFoilLV");
    G4LogicalVolume *logicGEMCathode = new G4LogicalVolume(solidGEMFoil, GEMFoilM, "GEMCathodeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 11.0 * mm), logicGEMCathode, "GEM Cathode", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 2);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 2.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 3);

    G4Box *GEMCuFoilBox1 = new G4Box("GEMCuFoilBox1", GEMHalfX, GEMHalfY, GEMCuFoilThickness / 2.0);
    G4Box *GEMCuFoilBox2 = new G4Box("GEMCuFoilBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMCuFoil = new G4SubtractionSolid("GEMCuFoilS", GEMCuFoilBox1, GEMCuFoilBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMCuFoil = new G4LogicalVolume(solidGEMCuFoil, GEMCuFoilM, "GEMCuFoilLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 11.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 2);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 3);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 4);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 5);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 6);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 2.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 7);

    // HyCal
    // The crystal surface should be at -3000.0 + 88.9 + 5640.0 = 2728.9 mm // 5640.0 from Weizhi
    G4double PbGlassL = 45.0 * cm;
    //G4double CrystalL = 18.0 * cm;
    G4double CrystalDiffL = 10.12 * cm;
    G4double CrystalSurf = 272.89 * cm; // Surface of the PWO
    G4double HyCalCenter = CrystalSurf - CrystalDiffL + PbGlassL / 2.0;

    // HyCal box
    G4double HyCalBoxCenter = HyCalCenter - 9.0 * cm + 30.0 * cm; // Check
    G4Box *HyCalBoxOuter = new G4Box("HyCalBoxOuter", 70.0 * cm, 70.0 * cm, 60.0 * cm);
    G4Box *HyCalBoxInner = new G4Box("HyCalBoxInner", 66.0 * cm, 66.0 * cm, 59.6 * cm);
    G4SubtractionSolid *HyCalBoxNoHole = new G4SubtractionSolid("HyCalBoxNoHole", HyCalBoxOuter, HyCalBoxInner);
    G4Tubs *HyCalBoxHole = new G4Tubs("HyCalBoxHole", 0, 25.0 * mm, 60.5 * cm, 0, twopi);
    G4SubtractionSolid *solidHyCalBox = new G4SubtractionSolid("HyCalBoxS", HyCalBoxNoHole, HyCalBoxHole);
    G4LogicalVolume *logicHyCalBox = new G4LogicalVolume(solidHyCalBox, HyCalBoxM, "HyCalBoxLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalBoxCenter), logicHyCalBox, "HyCal Box", logicWorld, false, 0);

    // HyCal container
    G4Box *HyCalConPiece1 = new G4Box("HyCalConPiece1", 58.21 * cm, 58.17 * cm, PbGlassL / 2.0);
    G4Box *HyCalConPiece2 = new G4Box("HyCalConPiece2", 35.30 * cm, 35.27 * cm, CrystalDiffL / 2.0 + 0.5 * mm);
    G4SubtractionSolid *HyCalConBox = new G4SubtractionSolid("HyCalConBox", HyCalConPiece1, HyCalConPiece2, 0, G4ThreeVector(0, 0, (CrystalDiffL - PbGlassL) / 2.0 - 0.5 * mm));
    G4Box *HyCalConHole = new G4Box("HyCalConHole", 2.0 * cm, 2.0 * cm, PbGlassL / 2.0 + 1.0 * mm);
    G4SubtractionSolid *solidHyCalCon = new G4SubtractionSolid("HyCalContainerS", HyCalConBox, HyCalConHole);
    G4LogicalVolume *logicHyCalCon = new G4LogicalVolume(solidHyCalCon, DefaultM, "HyCalContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter), logicHyCalCon, "HyCal Container", logicWorld, false, 0);

    // HyCal modules
    G4VSolid *solidAbsorber = new G4Box("HyCalModuleS", 1.025 * cm, 1.025 * cm, 90.0 * mm);
    G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber, HyCalModuleM, "HyCalModuleLV");
    HyCalParameterisation *param = new HyCalParameterisation("config/hycal.conf");
    new G4PVParameterised("HyCal Module", logicAbsorber, logicHyCalCon, kUndefined, param->GetNumber(), param, false);

    // Collimator container
    G4VSolid *CollConBox = new G4Box("CollConBox", 4.1 * cm, 4.1 * cm, 5.0 * cm);
    G4SubtractionSolid *solidCollCon = new G4SubtractionSolid("CollimatorContainerS", CollConBox, HyCalConHole);
    G4LogicalVolume *logicCollCon = new G4LogicalVolume(solidCollCon, DefaultM, "CollimatorContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter - PbGlassL / 2.0 + CrystalDiffL - 5.1 * cm), logicCollCon, "Collimator Container", logicWorld, false, 0);

    // Collimators
    G4VSolid *solidColl = new G4Box("CollimatorS", 1.025 * cm, 1.025 * cm, 5.0 * cm);
    G4LogicalVolume *logicColl = new G4LogicalVolume(solidColl, CollimatorM, "CollimatorLV");
    double pos_x[12] = { -3.075, -1.025, 1.025, 3.075, -3.075, 3.075, -3.075, 3.075, -3.075, -1.025, 1.025, 3.075};
    double pos_y[12] = { -3.075, -3.075, -3.075, -3.075, -1.025, -1.025, 1.025, 1.025, 3.075, 3.075, 3.075, 3.075};

    for (int i = 0; i < 12; ++i)
        new G4PVPlacement(0, G4ThreeVector(pos_x[i] * cm, pos_y[i] * cm, 0), logicColl, "Collimator", logicCollCon, false, i);

    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(fVisAtts[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    // Always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefinePRadSDs()
{
    if (fGEMSDOn) {
        TrackingDetectorSD *GEMSD = new TrackingDetectorSD("GEMSD", "GEM");
        G4SDManager::GetSDMpointer()->AddNewDetector(GEMSD);
        SetSensitiveDetector("GEMCathodeLV", GEMSD);
    }

    if (fHyCalSDOn != 0) {
        CalorimeterSD *HyCalSD = new CalorimeterSD("HyCalSD", "HC");
        G4SDManager::GetSDMpointer()->AddNewDetector(HyCalSD);
        SetSensitiveDetector("HyCalModuleLV", HyCalSD);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineDRadVolumes()
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *TargetM = G4Material::GetMaterial("D2Gas");

    if (fTargetMat == "H2Gas" || fTargetMat == "hydrogen")
        TargetM = G4Material::GetMaterial("H2Gas");

    G4Material *TargetCellM = G4Material::GetMaterial("Kapton");
    G4Material *TargetWindowM = G4Material::GetMaterial("Kapton");
    G4Material *RecoilDetectorM = G4Material::GetMaterial("Silicon");
    G4Material *KaptonWindowM = G4Material::GetMaterial("Kapton");
    G4Material *AluminumWindowM = G4Material::GetMaterial("Aluminum");
    G4Material *GEMFrameM = G4Material::GetMaterial("NemaG10");
    G4Material *GEMGasM = G4Material::GetMaterial("ArCO2");
    G4Material *GEMFoilM = G4Material::GetMaterial("Kapton");
    G4Material *GEMCuFoilM = G4Material::GetMaterial("Copper");
    G4Material *HeBagM = G4Material::GetMaterial("HeGas");
    G4Material *ScintillatorPlaneM = G4Material::GetMaterial("EJ204");
    G4Material *HyCalBoxM = G4Material::GetMaterial("Torlon");
    G4Material *CollimatorM = G4Material::GetMaterial("Tungsten");
    G4Material *HyCalModuleM = G4Material::GetMaterial("PbWO4");

    // World
    G4double WorldSizeXY = 150.0 * cm;
    G4double WorldSizeZ = 400.0 * cm;
    G4VSolid *solidWorld = new G4Box("WorldS", WorldSizeXY, WorldSizeXY, WorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, DefaultM, "WorldLV");
    G4VPhysicalVolume *physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    // Target
    // Target container
    G4VSolid *solidTargetCon = new G4Tubs("TargetContainerS", 0, fTargetR + 0.1 * cm, fTargetHalfL + 0.1 * cm, 0, twopi);
    G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, DefaultM, "TargetContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fTargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

    // Target material
    G4VSolid *solidTarget = new G4Tubs("TargetS", 0, fTargetR, fTargetHalfL, 0, twopi);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetM, "TargetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

    // Target cell
    G4double CellThickness = 0.5 * mm;
    G4VSolid *solidCell = new G4Tubs("TargetCellS", fTargetR, fTargetR + CellThickness, fTargetHalfL, 0, twopi);
    G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

    // Windows
    G4double CellApertureR = 2.0 * mm;
    G4double CellWinThickness = 7.5 * um;
    G4VSolid *solidCellWin = new G4Tubs("TargetWindowS", CellApertureR, fTargetR + CellThickness, CellWinThickness / 2.0, 0, twopi);
    G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, TargetWindowM, "TargetWindowS");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -fTargetHalfL - CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +fTargetHalfL + CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 1);

    // Recoil detector
    G4double RecoilDetCenter = fRecoilDetCenter - fTargetCenter;
    G4double RecoilDetAng = twopi / fRecoilDetNSeg;
    G4double RecoilDetOR = fTargetR * cos(RecoilDetAng / 2.0) - 0.5 * mm;
    G4double RecoilDetIR = RecoilDetOR - fRecoilDetThickness;
    G4double rInnerRD[] = {RecoilDetIR, RecoilDetIR};
    G4double rOuterRD[] = {RecoilDetOR, RecoilDetOR};
    G4double zPlaneRD[] = { -fRecoilDetHalfL, fRecoilDetHalfL};
    G4VSolid *solidRecoilDet = new G4Polyhedra("RecoilDetectorS", 0, twopi, fRecoilDetNSeg, 2, zPlaneRD, rInnerRD, rOuterRD);
    G4LogicalVolume *logicRecoilDet = new G4LogicalVolume(solidRecoilDet, RecoilDetectorM, "RecoilDetectorLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, RecoilDetCenter), logicRecoilDet, "Recoil Detector", logicTarget, false, 0);

    // Additional window (Downstream chamber window)
    G4double KaptonWinCenter = fTargetCenter + 74.0 * mm;
    G4double KaptonWinApertureR = 22.8 * mm;
    G4double KaptonWinOR = 8.00 * cm;
    G4double KaptonWinThickness = 7.5 * um;
    G4Tubs *solidKaptonWin = new G4Tubs("KaptonWindowS", KaptonWinApertureR, KaptonWinOR, KaptonWinThickness / 2.0, 0, twopi);
    G4LogicalVolume *logicKaptonWin = new G4LogicalVolume(solidKaptonWin, KaptonWindowM, "KaptonWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, KaptonWinCenter), logicKaptonWin, "Kapton Window", logicWorld, false, 0);

    // Additional window (Vacuum box window)
    G4double AlWinCenter = fGEM1Center - 10.0 * cm;
    G4double AlWinMaxR = 78.11 * cm;
    G4double ArcDistance = 5.59 * cm;
    G4double ArcEndR = (ArcDistance * ArcDistance + AlWinMaxR * AlWinMaxR) / (2 * ArcDistance);
    G4double ArcEndThickness = 1.6 * mm;
    G4double AlWinApertureR = 3.0 * cm;
    G4VSolid *solidAlWin = new G4Sphere("AluminumWindowS", ArcEndR - ArcEndThickness, ArcEndR, 0, twopi, pi - asin(AlWinMaxR / ArcEndR), asin(AlWinMaxR / ArcEndR) - asin((AlWinApertureR + 0.1 * mm) / ArcEndR));
    G4LogicalVolume *logicAlWin = new G4LogicalVolume(solidAlWin, AluminumWindowM, "AluminumWindowS");
    new G4PVPlacement(0, G4ThreeVector(0, 0,  AlWinCenter + ArcEndR - ArcDistance), logicAlWin, "Aluminum Window", logicWorld, false, 0);

    // GEM
    G4double GEMGap = 4.0 * cm; // Gap between two GEM
    G4double GEMHalfX = 55.04 * cm / 2.0;
    G4double GEMHalfY = 122.88 * cm / 2.0;
    G4double GEMHalfThickness = 7.0 * mm;
    G4double GEMHoleR = 2.2 * cm;
    G4double GEMCenterHalfXY = 7.4 * cm / 2.0;
    G4double GEMFrameWidth = 1.5 * cm;
    G4double GEMCenterOffset = GEMHalfX + GEMFrameWidth - GEMCenterHalfXY;

    // GEM Container
    G4Box *GEMConBox = new G4Box("GEMConBox", 1.0 * m, 1.0 * m, (GEMGap + 2.0 * GEMHalfThickness + 1.0 * mm) / 2.0);
    G4Tubs *GEMConTube = new G4Tubs("GEMConTube", 0, GEMHoleR, (GEMGap + 2.0 * GEMHalfThickness + 1.0 * mm) / 2.0 + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMCon = new G4SubtractionSolid("GEMContainerS", GEMConBox, GEMConTube);
    G4LogicalVolume *logicGEMCon = new G4LogicalVolume(solidGEMCon, DefaultM, "GEMContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fGEM1Center), logicGEMCon, "GEM Container", logicWorld, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, fGEM2Center), logicGEMCon, "GEM Container", logicWorld, false, 2);

    // GEM
    G4Box *GEMBox = new G4Box("GEMBox", GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfThickness);
    G4Tubs *GEMTube = new G4Tubs("GEMTube", 0, GEMHoleR, GEMHalfThickness + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEM = new G4SubtractionSolid("GEMS", GEMBox, GEMTube, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEM = new G4LogicalVolume(solidGEM, DefaultM, "GEMLV");
    new G4PVPlacement(0, G4ThreeVector(GEMCenterOffset, 0, GEMGap / 2.0), logicGEM, "GEM L", logicGEMCon, false, 0);
    G4RotationMatrix rmGEM;
    rmGEM.rotateZ(180.0 * deg);
    new G4PVPlacement(G4Transform3D(rmGEM, G4ThreeVector(-GEMCenterOffset, 0, -GEMGap / 2.0)), logicGEM, "GEM R", logicGEMCon, false, 1);

    // GEM Gas
    G4Box *GEMGasBox1 = new G4Box("GEMGasBox1", GEMHalfX, GEMHalfY, GEMHalfThickness);
    G4Box *GEMGasBox2 = new G4Box("GEMGasBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMGas = new G4SubtractionSolid("GEMGasS", GEMGasBox1, GEMGasBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMGas = new G4LogicalVolume(solidGEMGas, GEMGasM, "GEMGasLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMGas, "GEM Gas", logicGEM, false, 0);

    // GEM Frame
    G4Box *GEMFrameBox1 = new G4Box("GEMFrameBox1", GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfThickness);
    G4Box *GEMFrameBox2 = new G4Box("GEMFrameBox2", GEMHalfX, GEMHalfY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid("GEMFrameS", GEMFrameBox1, GEMFrameBox2);
    G4LogicalVolume *logicGEMFrame = new G4LogicalVolume(solidGEMFrame, GEMFrameM, "GEMFrameLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMFrame, "GEM Frame", logicGEM, false, 0);
    G4Box *GEMPipeBox = new G4Box("GEMPipeBox", GEMCenterHalfXY - GEMFrameWidth / 2.0, GEMCenterHalfXY, GEMHalfThickness);
    G4Tubs *GEMPipeTube = new G4Tubs("GEMPipeTube", 0, GEMHoleR, GEMHalfThickness + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMPipe = new G4SubtractionSolid("GEMPipeS", GEMPipeBox, GEMPipeTube, 0, G4ThreeVector(-GEMFrameWidth / 2.0, 0, 0));
    G4LogicalVolume *logicGEMPipe = new G4LogicalVolume(solidGEMPipe, GEMFrameM, "GEMPipeLV");
    new G4PVPlacement(0, G4ThreeVector(-GEMCenterOffset + GEMFrameWidth / 2.0, 0, 0), logicGEMPipe, "GEM Pipe", logicGEM, false, 0);

    // GEM Foil
    G4double GEMFoilThickness = 50.0 * um;
    G4double GEMCuFoilThickness = 5.0 * um;
    G4double GEMWinThickness = 25.0 * um;
    G4Box *GEMWinBox1 = new G4Box("GEMWinBox1", GEMHalfX, GEMHalfY, GEMWinThickness / 2.0);
    G4Box *GEMWinBox2 = new G4Box("GEMWinBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMWin = new G4SubtractionSolid("GEMWinS", GEMWinBox1, GEMWinBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMWin = new G4LogicalVolume(solidGEMWin, GEMFoilM, "GEMWinLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -GEMHalfThickness + GEMWinThickness / 2.0), logicGEMWin, "GEM Window", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - GEMWinThickness / 2.0), logicGEMWin, "GEM Window", logicGEMGas, false, 1);

    G4Box *GEMFoilBox1 = new G4Box("GEMFoilBox1", GEMHalfX, GEMHalfY, GEMFoilThickness / 2.0);
    G4Box *GEMFoilBox2 = new G4Box("GEMFoilBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid("GEMFoilS", GEMFoilBox1, GEMFoilBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMFoil = new G4LogicalVolume(solidGEMFoil, GEMFoilM, "GEMFoilLV");
    G4LogicalVolume *logicGEMCathode = new G4LogicalVolume(solidGEMFoil, GEMFoilM, "GEMCathodeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 11.0 * mm), logicGEMCathode, "GEM Cathode", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 2);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 2.0 * mm), logicGEMFoil, "GEM Foil", logicGEMGas, false, 3);

    G4Box *GEMCuFoilBox1 = new G4Box("GEMCuFoilBox1", GEMHalfX, GEMHalfY, GEMCuFoilThickness / 2.0);
    G4Box *GEMCuFoilBox2 = new G4Box("GEMCuFoilBox2", GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfThickness + 0.1 * mm);
    G4SubtractionSolid *solidGEMCuFoil = new G4SubtractionSolid("GEMCuFoilS", GEMCuFoilBox1, GEMCuFoilBox2, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMCuFoil = new G4LogicalVolume(solidGEMCuFoil, GEMCuFoilM, "GEMCuFoilLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 11.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 8.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 2);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 3);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 6.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 4);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 5);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 4.0 * mm + GEMFoilThickness / 2.0 + GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 6);
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfThickness - 2.0 * mm - GEMFoilThickness / 2.0 - GEMCuFoilThickness / 2.0), logicGEMCuFoil, "GEM Copper", logicGEMGas, false, 7);

    // He bag (Only He gas for now)
    G4Box *HeBagBox = new G4Box("HeBagBox", 1.0 * m, 1.0 * m, (fGEM2Center - fGEM1Center - 5.6 * cm) / 2.0);
    G4Tubs *HeBagTube = new G4Tubs("HeBagTube", 0, 22.0 * mm, (fGEM2Center - fGEM1Center - 5.6 * cm + 1.0 * mm) / 2.0, 0, twopi);
    G4SubtractionSolid *solidHeBag = new G4SubtractionSolid("HeBagS", HeBagBox, HeBagTube);
    G4LogicalVolume *logicHeBag = new G4LogicalVolume(solidHeBag, HeBagM, "HeBagLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, (fGEM1Center + fGEM2Center) / 2.0), logicHeBag, "He Bag", logicWorld, false, 0);

    // Scintillator plane
    G4double SciPlaneThickness = 5.0 * mm;
    G4double SciPlaneHalfX = 75.0 * cm;
    G4double SciPlaneHalfY = 75.0 * cm;
    G4VSolid *solidSciPlane = new G4Box("ScintillatorPlaneS", SciPlaneHalfX, SciPlaneHalfY, SciPlaneThickness / 2.0);
    G4LogicalVolume *logicSciPlane = new G4LogicalVolume(solidSciPlane, ScintillatorPlaneM, "ScintillatorPlaneLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fSciPlaneCenter), logicSciPlane, "Scintillator Plane", logicWorld, false, 0);

    // HyCal
    G4double PbGlassL = 45.0 * cm;
    //G4double CrystalL = 18.0 * cm;
    G4double CrystalDiffL = 10.12 * cm;
    G4double HyCalCenter = fCrystalSurf - CrystalDiffL + PbGlassL / 2.0;

    // HyCal box
    G4double HyCalBoxCenter = HyCalCenter - 9.0 * cm + 30.0 * cm; // Check
    G4Box *HyCalBoxOuter = new G4Box("HyCalBoxOuter", 70.0 * cm, 70.0 * cm, 60.0 * cm);
    G4Box *HyCalBoxInner = new G4Box("HyCalBoxInner", 66.0 * cm, 66.0 * cm, 59.6 * cm);
    G4SubtractionSolid *HyCalBoxNoHole = new G4SubtractionSolid("HyCalBoxNoHole", HyCalBoxOuter, HyCalBoxInner);
    G4Tubs *HyCalBoxHole = new G4Tubs("HyCalBoxHole", 0, 25.0 * mm, 60.5 * cm, 0, twopi);
    G4SubtractionSolid *solidHyCalBox = new G4SubtractionSolid("HyCalBoxS", HyCalBoxNoHole, HyCalBoxHole);
    G4LogicalVolume *logicHyCalBox = new G4LogicalVolume(solidHyCalBox, HyCalBoxM, "HyCalBoxLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalBoxCenter), logicHyCalBox, "HyCal Box", logicWorld, false, 0);

    // HyCal container
    G4Box *HyCalConPiece1 = new G4Box("HyCalConPiece1", 58.21 * cm, 58.17 * cm, PbGlassL / 2.0);
    G4Box *HyCalConPiece2 = new G4Box("HyCalConPiece2", 35.30 * cm, 35.27 * cm, CrystalDiffL / 2.0 + 0.5 * mm);
    G4SubtractionSolid *HyCalConBox = new G4SubtractionSolid("HyCalConBox", HyCalConPiece1, HyCalConPiece2, 0, G4ThreeVector(0, 0, (CrystalDiffL - PbGlassL) / 2.0 - 0.5 * mm));
    G4Box *HyCalConHole = new G4Box("HyCalConHole", 2.0 * cm, 2.0 * cm, PbGlassL / 2.0 + 1.0 * mm);
    G4SubtractionSolid *solidHyCalCon = new G4SubtractionSolid("HyCalContainerS", HyCalConBox, HyCalConHole);
    G4LogicalVolume *logicHyCalCon = new G4LogicalVolume(solidHyCalCon, DefaultM, "HyCalContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter), logicHyCalCon, "HyCal Container", logicWorld, false, 0);

    // HyCal modules
    G4VSolid *solidAbsorber = new G4Box("HyCalModuleS", 1.025 * cm, 1.025 * cm, 90.0 * mm);
    G4LogicalVolume *logicAbsorber = new G4LogicalVolume(solidAbsorber, HyCalModuleM, "HyCalModuleLV");
    HyCalParameterisation *param = new HyCalParameterisation("config/hycal.conf");
    new G4PVParameterised("HyCal Module", logicAbsorber, logicHyCalCon, kUndefined, param->GetNumber(), param, false);

    // Collimator container
    G4VSolid *CollConBox = new G4Box("CollConBox", 4.1 * cm, 4.1 * cm, 5.0 * cm);
    G4SubtractionSolid *solidCollCon = new G4SubtractionSolid("CollimatorContainerS", CollConBox, HyCalConHole);
    G4LogicalVolume *logicCollCon = new G4LogicalVolume(solidCollCon, DefaultM, "CollimatorContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalCenter - PbGlassL / 2.0 + CrystalDiffL - 5.1 * cm), logicCollCon, "Collimator Container", logicWorld, false, 0);

    // Collimators
    G4VSolid *solidColl = new G4Box("CollimatorS", 1.025 * cm, 1.025 * cm, 5.0 * cm);
    G4LogicalVolume *logicColl = new G4LogicalVolume(solidColl, CollimatorM, "CollimatorLV");
    double pos_x[12] = { -3.075, -1.025, 1.025, 3.075, -3.075, 3.075, -3.075, 3.075, -3.075, -1.025, 1.025, 3.075};
    double pos_y[12] = { -3.075, -3.075, -3.075, -3.075, -1.025, -1.025, 1.025, 1.025, 3.075, 3.075, 3.075, 3.075};

    for (int i = 0; i < 12; ++i)
        new G4PVPlacement(0, G4ThreeVector(pos_x[i] * cm, pos_y[i] * cm, 0), logicColl, "Collimator", logicCollCon, false, i);

    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(fVisAtts[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    // Always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineDRadSDs()
{
    if (fRecoilDetSDOn) {
        StandardDetectorSD *RecoilDetSD = new StandardDetectorSD("RecoilDetectorSD", "RD");
        G4SDManager::GetSDMpointer()->AddNewDetector(RecoilDetSD);
        SetSensitiveDetector("RecoilDetectorLV", RecoilDetSD);
    }

    if (fGEMSDOn) {
        TrackingDetectorSD *GEMSD = new TrackingDetectorSD("GEMSD", "GEM");
        G4SDManager::GetSDMpointer()->AddNewDetector(GEMSD);
        SetSensitiveDetector("GEMCathodeLV", GEMSD);
    }

    if (fSciPlaneSDOn) {
        StandardDetectorSD *SciPlaneSD = new StandardDetectorSD("ScintillatorPlaneSD", "SP");
        G4SDManager::GetSDMpointer()->AddNewDetector(SciPlaneSD);
        SetSensitiveDetector("ScintillatorPlaneLV", SciPlaneSD);
    }

    if (fHyCalSDOn == 1) {
        CalorimeterSD *HyCalSD = new CalorimeterSD("HyCalSD", "HC");
        G4SDManager::GetSDMpointer()->AddNewDetector(HyCalSD);
        SetSensitiveDetector("HyCalModuleLV", HyCalSD);
    } else if (fHyCalSDOn == 2) {
        StandardDetectorSD *HyCalSD = new StandardDetectorSD("HyCalSD", "HC");
        G4SDManager::GetSDMpointer()->AddNewDetector(HyCalSD);
        SetSensitiveDetector("HyCalModuleLV", HyCalSD);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
