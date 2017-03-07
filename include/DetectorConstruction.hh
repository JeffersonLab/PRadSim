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

class DetectorMessenger;
class G4String;
class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(G4String conf);
    virtual ~DetectorConstruction();

public:
    G4VPhysicalVolume *Construct();

    inline void SetTargetPos(G4double z);
    inline void SetRecoilDetectorPos(G4double z);
    inline void SetGEMPos(G4double z1, G4double z2);
    inline void SetScitillatorPlanePos(G4double z);
    inline void SetHyCalPos(G4double z);

    inline void SetRecoilDetector(G4int n, G4double ir, G4double l, G4double t);

    void UpdateGeometry();

    inline const G4VPhysicalVolume *GetPhysiWorld();

private:
    G4String fConfig;

    G4double fTargetCenter;

    G4double fRecoilDetCenter;
    G4int fRecoilDetNSeg;
    G4double fRecoilDetIR;
    G4double fRecoilDetHalfL;
    G4double fRecoilDetThickness;

    G4double fGEM1Center;
    G4double fGEM2Center;

    G4double fSciPlaneCenter;

    G4double fCrystalSurf;

    G4VPhysicalVolume *physiWorld;

private:
    DetectorMessenger *detectorMessenger; //pointer to the Messenger
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
    if (z1 > -9999) fGEM1Center = z1;

    if (z2 > -9999) fGEM2Center = z2;
}

inline void DetectorConstruction::SetScitillatorPlanePos(G4double z)
{
    fSciPlaneCenter = z;
}

inline void DetectorConstruction::SetHyCalPos(G4double z)
{
    fCrystalSurf = z;
}

inline void DetectorConstruction::SetRecoilDetector(G4int n, G4double ir, G4double l, G4double t)
{
    if (n > -9999) fRecoilDetNSeg = n;

    if (ir > -9999) fRecoilDetIR = ir;

    if (l > -9999) fRecoilDetHalfL = l;

    if (t > -9999) fRecoilDetThickness = t;
}

inline const G4VPhysicalVolume *DetectorConstruction::GetPhysiWorld()
{
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
