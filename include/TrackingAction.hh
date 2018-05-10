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
// TrackingAction.hh
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

#include "globals.hh"

#define MaxTracks 3000

class TrackingMessenger;

class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction
{
public:
    TrackingAction();
    virtual ~TrackingAction();

    void PreUserTrackingAction(const G4Track *);
    void PostUserTrackingAction(const G4Track *);

    inline void SetNoSecondary(G4bool val);
    void SetSaveTrackInfo(G4bool val);

    void Clear();

private:
    void Register(TTree *);

    G4bool fNoSecondary;
    G4bool fSaveTrackInfo;
    G4bool fRegistered;

    int fN;
    int fPID[MaxTracks]; // Particle ID
    int fTID[MaxTracks]; // Track ID
    int fPTID[MaxTracks]; // Parent Track ID
    double fX[MaxTracks];
    double fY[MaxTracks];
    double fZ[MaxTracks];
    int fProcessID[MaxTracks];

    TrackingMessenger *trackingMessenger;
};

inline void TrackingAction::SetNoSecondary(G4bool val)
{
    fNoSecondary = val;
}

inline void TrackingAction::SetSaveTrackInfo(G4bool val)
{
    fSaveTrackInfo = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
