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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Tubs;
class G4Box;
class G4Sphere;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4VPVParameterisation;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

    DetectorConstruction();
   ~DetectorConstruction();

public:
     G4VPhysicalVolume* Construct();
     void SetAbsorberMaterial (G4String);
     void SetAbsorber2Material (G4String);
     void SetAbsorberThickness(G4double);
     void SetTargetMaterial (G4String);
     void SetCellMaterial (G4String);
     void SetNeckMaterial (G4String);
     void SetColMaterial (G4String);
     void SetVacBoxMaterial (G4String);
     void SetGEMMaterial (G4String);
     void SetCalorSizeXY(G4double);
     void SetMagField(G4double);
     void UpdateGeometry();

public:
     void PrintCalorParameters();
     G4double GetWorldSizeZ() {return WorldSizeZ;};
     G4double GetWorldSizeXY() {return WorldSizeXY;};
     G4double GetCalorThickness() {return CalorThickness;};
     G4double GetCalorSizeXY() {return CalorSizeXY;};

     G4Material* GetAbsorberMaterial() {return AbsorberMaterial;};
     G4Material* GetAbsorber2Material() {return Absorber2Material;};
     G4double GetAbsorberThickness() {return AbsorberThickness;};
     G4Material* GetTargetMaterial() {return TargetMaterial;};
     G4Material* GetCellMaterial() {return CellMaterial;};
     G4Material* GetNeckMaterial() {return NeckMaterial;};
     G4Material* GetColMaterial() {return ColMaterial;};
     G4Material* GetVacBoxMaterial() {return VacBoxMaterial;};
     G4Material* GetGEMMaterial() {return GEMMaterial;};

     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};
     const G4VPhysicalVolume* GetAbsorber() {return physiAbsorber;};
     const G4VPhysicalVolume* GetCell() {return physiCell;};

private:
     G4Material* AbsorberMaterial;
     G4Material* Absorber2Material;
     G4Material* TargetMaterial;
     G4Material* CellMaterial;
     G4Material* NeckMaterial;
     G4Material* ColMaterial;
     G4Material* VacBoxMaterial;
     G4Material* GEMMaterial;
     G4Material* defaultMaterial;

     G4double CalorSizeXY;
     G4double CalorThickness;
     G4double AbsorberThickness;
     G4double SurfaceDiff;

     G4double CellR;
     G4double CellThickness;
     G4double CellHalfL;
     G4double NeckHalfL;
     G4double ApertureR;

     G4double fCylinderD;
     G4double ArcEndR;
     G4double ArcDistance;
     G4double WinThickness;

     G4double VacBoxtoHyCal;
     G4double HyCalCenter;
     G4double TargetCenter;

     G4double WorldSizeXY;
     G4double WorldSizeZ;

     G4Box*                 solidWorld;
     G4LogicalVolume*       logicWorld;
     G4VPhysicalVolume*     physiWorld;

     G4Tubs*                solidCelltube;
     G4Tubs*                solidNecktube;
     G4SubtractionSolid*    solidNeck;
     G4LogicalVolume*       logicNeck;
     G4VPhysicalVolume*     physiNeck;

     G4Tubs*                solidCell;
     G4LogicalVolume*       logicCell;
     G4VPhysicalVolume*     physiCell;

     G4Tubs*                solidTarget;
     G4LogicalVolume*       logicTarget;
     G4VPhysicalVolume*     physiTarget;

     G4Tubs*                solidCellWin;
     G4LogicalVolume*       logicCellWin;
     G4VPhysicalVolume*     physiWinIn;
     G4VPhysicalVolume*     physiWinOut;

     G4Sphere*              solidVacBoxSph;
     G4Box*                 solidVacBoxBox;
     G4SubtractionSolid*    solidVacBoxShl;
     G4Tubs*                solidVacBoxTub;
     G4SubtractionSolid*    solidVacBoxEnd;
     G4LogicalVolume*       logicVacBoxEnd;
     G4VPhysicalVolume*     physiVacBoxEnd;

     G4Tubs*                solidFlange;
     G4LogicalVolume*       logicFlange;
     G4VPhysicalVolume*     physiFlange;

     G4Tubs*                solidVacTube;
     G4LogicalVolume*       logicVacTube;
     G4VPhysicalVolume*     physiVacTube;

     G4Box*                 solidColl;
     G4LogicalVolume*       logicColl;
     G4VPhysicalVolume*     physiColl[12];

     G4LogicalVolume*       logicGEMFrame;
     G4VPhysicalVolume*     physiGEMFrame1;
     G4VPhysicalVolume*     physiGEMFrame2;

     G4LogicalVolume*       logicGEM;
     G4VPhysicalVolume*     physiGEM1;
     G4VPhysicalVolume*     physiGEM2;

     G4SubtractionSolid*    solidHyCalBox;
     G4LogicalVolume*       logicHyCalBox;
     G4VPhysicalVolume*     physiHyCalBox;

     G4SubtractionSolid*    solidCalor;
     G4LogicalVolume*       logicCalor;    //pointer to the logical Calor
     G4VPhysicalVolume*     physiCalor;    //pointer to the physical Calor

     G4SubtractionSolid*    solidOuterCalor;
     G4LogicalVolume*       logicOuterCalor;    //pointer to the logical Calor
     G4VPhysicalVolume*     physiOuterCalor;    //pointer to the physical Calor

     G4Box*                 solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*       logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume*     physiAbsorber; //pointer to the physical Absorber

     G4Box*                 solidAbsorber2; //pointer to the solid Absorber
     G4LogicalVolume*       logicAbsorber2; //pointer to the logical Absorber
     G4VPhysicalVolume*     physiAbsorber2; //pointer to the physical Absorber

     G4VPVParameterisation* CalorimeterParameter;
     G4VPVParameterisation* LeadGlassPartParameter;

     G4UniformMagField*     magField;      //pointer to the magnetic field
     DetectorMessenger*     detectorMessenger;  //pointer to the Messenge

private:
     void DefineMaterials();
};


#endif

