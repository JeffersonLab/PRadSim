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
// StandardDetectorSD.hh
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StandardDetectorSD_h
#define StandardDetectorSD_h 1

#include "G4VSensitiveDetector.hh"

#include "StandardHit.hh"

#include "G4String.hh"

#define MaxNHits 300

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StandardDetectorSD : public G4VSensitiveDetector
{
public:
    StandardDetectorSD(G4String name, G4String abbrev);
    virtual ~StandardDetectorSD();

    virtual void Initialize(G4HCofThisEvent *);
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
    virtual void EndOfEvent(G4HCofThisEvent *);

protected:
    virtual void Register(TTree *);
    virtual void Clear();

    G4int fID;
    G4String fAbbrev;

    StandardHitsCollection *fHitsCollection;

    bool fRegistered;

    int fN;
    int fPID[MaxNHits]; // Particle ID
    int fTID[MaxNHits]; // Track ID
    int fPTID[MaxNHits]; // Parent Track ID
    int fDID[MaxNHits];
    double fX[MaxNHits];
    double fY[MaxNHits];
    double fZ[MaxNHits];
    double fMomentum[MaxNHits];
    double fTheta[MaxNHits];
    double fPhi[MaxNHits];
    double fTime[MaxNHits];
    double fEdep[MaxNHits];
    double fTrackL[MaxNHits];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
