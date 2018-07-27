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
// CalorimeterSD.hh
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterSD_h
#define CalorimeterSD_h 1

#include "StandardDetectorSD.hh"

#include "CalorimeterHit.hh"

#include "Math/Interpolator.h"
#include "G4String.hh"

#define InterpolPoints 181
#define InterpolType ROOT::Math::Interpolation::kCSPLINE

#define NModules 1728

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterSD: public StandardDetectorSD
{
public:
    CalorimeterSD(G4String name, G4String abbrev, G4String pwo_filename);
    virtual ~CalorimeterSD();

    virtual void Initialize(G4HCofThisEvent *);
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
    virtual void EndOfEvent(G4HCofThisEvent *);

    inline void SetAttenuationLG(G4double val);
    inline void SetReflectanceLG(G4double val);

protected:
    virtual void Register(TTree *);
    virtual void Clear();

    CalorimeterHitsCollection *fCalorHitsCollection;

    double fAttenuationLG;
    double fReflectanceLG;

    double fTotalEdep;
    double fTotalTrackL;
    double fModuleEdep[NModules];
    double fModuleTrackL[NModules];

    ROOT::Math::Interpolator *fInterpolator;
};

inline void CalorimeterSD::SetAttenuationLG(G4double val)
{
    fAttenuationLG = val;
}

inline void CalorimeterSD::SetReflectanceLG(G4double val)
{
    fReflectanceLG = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
