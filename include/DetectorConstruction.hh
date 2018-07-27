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
// DetectorConstruction.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4String.hh"

#include <map>

class DetectorMessenger;

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(G4String conf);
    virtual ~DetectorConstruction();

public:
    G4VPhysicalVolume *Construct();
    void ConstructSDandField();

    inline void SetTargetPos(G4double z);
    inline void SetRecoilDetectorPos(G4double z);
    inline void SetGEMPos(G4double z1, G4double z2);
    inline void SetScitillatorPlanePos(G4double z);
    inline void SetHyCalPos(G4double z);

    inline void SetTarget(G4double ir, G4double l);
    inline void SetRecoilDetector(G4int n, G4double ir, G4double l, G4double t1, G4double t2);

    inline void SetTargetMaterial(G4String val);
    inline void SetTargetDensityRatio(G4double r);

    inline void SetExtDensityRatio(G4double r);

    inline void EnableSD(G4String detname);
    inline void DisableSD(G4String detname);

    inline void SetAttenuationLG(G4double val);
    inline void SetReflectanceLG(G4double val);

private:
    void DefineMaterials();

    G4VPhysicalVolume *DefinePRadVolumes();
    void DefinePRadSDs();

    G4VPhysicalVolume *DefineDRadVolumes();
    void DefineDRadSDs();

    G4VPhysicalVolume *DefineTestVolumes();
    void DefineTestSDs();

    void AddVaccumBox(G4LogicalVolume *mother);
    void AddGEM(G4LogicalVolume *mother, int layerid, bool culess);
    void AddHyCal(G4LogicalVolume *mother);

    G4String fConfig;

    std::map<G4String, G4VisAttributes *> fVisAtts;

    G4double fWorldSizeXY;
    G4double fWorldSizeZ;

    G4double fTargetCenter;
    G4double fTargetR;
    G4double fTargetHalfL;
    G4String fTargetMat;
    G4double fTargetDensityRatio;

    G4int fRecoilDetNSeg;
    G4double fRecoilDetCenter;
    G4double fRecoilDetR;
    G4double fRecoilDetHalfL;
    G4double fRecoilDetL1Thickness;
    G4double fRecoilDetL2Thickness;

    G4double fDownChamberCenter;
    G4double fVacBoxCenter;

    G4double fGEMCenter[10];

    G4double fSciPlaneCenter;

    G4double fCrystalSurf;

    G4double fExtDensityRatio;

    G4bool fTargetSDOn;
    G4bool fRecoilDetSDOn;
    G4bool fGEMSDOn;
    G4bool fSciPlaneSDOn;
    G4bool fHyCalSDOn;
    G4bool fVirtualSDOn;

    G4double fAttenuationLG;
    G4double fReflectanceLG;

    DetectorMessenger *detectorMessenger; // pointer to the messenger
};

inline void DetectorConstruction::SetTargetPos(G4double z)
{
    fTargetCenter = z;
}

inline void DetectorConstruction::SetRecoilDetectorPos(G4double z)
{
    fRecoilDetCenter = z;
}

inline void DetectorConstruction::SetGEMPos(G4double z1, G4double z2)
{
    if (z1 > -9999) fGEMCenter[0] = z1;

    if (z2 > -9999) fGEMCenter[1] = z2;
}

inline void DetectorConstruction::SetScitillatorPlanePos(G4double z)
{
    fSciPlaneCenter = z;
}

inline void DetectorConstruction::SetHyCalPos(G4double z)
{
    fCrystalSurf = z;
}

inline void DetectorConstruction::SetTarget(G4double ir, G4double l)
{
    if (ir > -9999) fTargetR = ir;

    if (l > -9999) fTargetHalfL = l;
}

inline void DetectorConstruction::SetRecoilDetector(G4int n, G4double ir, G4double l, G4double t1, G4double t2)
{
    if (n > -9999) fRecoilDetNSeg = n;

    if (ir > -9999) fRecoilDetR = ir;

    if (l > -9999) fRecoilDetHalfL = l;

    if (t1 > -9999) fRecoilDetL1Thickness = t1;

    if (t2 > -9999) fRecoilDetL2Thickness = t2;
}

inline void DetectorConstruction::SetTargetMaterial(G4String val)
{
    fTargetMat = val;
}

inline void DetectorConstruction::SetTargetDensityRatio(G4double r)
{
    fTargetDensityRatio = r;
}

inline void DetectorConstruction::SetExtDensityRatio(G4double r)
{
    fExtDensityRatio = r;
}

inline void DetectorConstruction::EnableSD(G4String detname)
{
    if (detname == "Target") fTargetSDOn = true;

    if (detname == "Recoil Detector") fRecoilDetSDOn = true;

    if (detname == "GEM") fGEMSDOn = true;

    if (detname == "Scintillator Plane") fSciPlaneSDOn = true;

    if (detname == "HyCal") fHyCalSDOn = true;

    if (detname == "Virtual Detector") fVirtualSDOn = true;
}

inline void DetectorConstruction::DisableSD(G4String detname)
{
    if (detname == "Target") fTargetSDOn = false;

    if (detname == "Recoil Detector") fRecoilDetSDOn = false;

    if (detname == "GEM") fGEMSDOn = false;

    if (detname == "Scintillator Plane") fSciPlaneSDOn = false;

    if (detname == "HyCal") fHyCalSDOn = false;

    if (detname == "Virtual Detector") fVirtualSDOn = false;
}

inline void DetectorConstruction::SetAttenuationLG(G4double val)
{
    fAttenuationLG = val;
}

inline void DetectorConstruction::SetReflectanceLG(G4double val)
{
    fReflectanceLG = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
