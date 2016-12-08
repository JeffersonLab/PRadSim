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

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <map>

class Digitization;
class RootTree;
class DetectorMessenger;
class G4String;
class G4Material;
class G4VPhysicalVolume;
class G4VisAttributes;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

public:
    G4VPhysicalVolume *Construct();

    void SetTargetMaterial(G4String);
    void UpdateGeometry();

    const G4VPhysicalVolume *GetPhysiWorld() {return physiWorld;}

private:
    G4VPhysicalVolume *physiWorld;

    G4Material *VacuumMaterial;
    G4Material *TargetMaterial;
    G4Material *CellMaterial;
    G4Material *VacuumBoxMaterial;
    G4Material *GEMFrameMaterial;
    G4Material *GEMFoilMaterial;
    G4Material *GEMGasMaterial;
    G4Material *HyCalBoxMaterial;
    G4Material *CollimatorMaterial;
    G4Material *CenterHyCalMaterial;
    G4Material *OuterHyCalMaterial;
    G4Material *defaultMaterial;

    std::map<G4String, G4VisAttributes *> mVisAtt;

    Digitization *daq_system;
    RootTree *otree;

private:
    void DefineMaterials();

private:
    DetectorMessenger *detectorMessenger; //pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
