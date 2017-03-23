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
// TrackInformation.cc
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackInformation.hh"

#include "G4Allocator.hh"
#include "G4Track.hh"
#include "G4VUserTrackInformation.hh"

#include <unordered_map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Allocator<TrackInformation> TrackInformationAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackInformation::TrackInformation() : G4VUserTrackInformation()
{
    fOriTrackID = 0;
    fOriParticle = NULL;
    fOriPosition = G4ThreeVector(0, 0, 0);
    fOriMomentum = G4ThreeVector(0, 0, 0);
    fOriEnergy = 0;
    fOriTime = 0;
    fAncestorMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackInformation::~TrackInformation()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackInformation::TrackInformation(const G4Track *aTrack)
{
    fOriTrackID = aTrack->GetTrackID();
    fOriParticle = (G4ParticleDefinition *)(aTrack->GetParticleDefinition());
    fOriPosition = aTrack->GetPosition();
    fOriMomentum = aTrack->GetMomentum();
    fOriEnergy = aTrack->GetTotalEnergy();
    fOriTime = aTrack->GetGlobalTime();
    fAncestorMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackInformation::TrackInformation(const TrackInformation *aTrackInfo)
{
    fOriTrackID = aTrackInfo->fOriTrackID;
    fOriParticle = aTrackInfo->fOriParticle;
    fOriPosition = aTrackInfo->fOriPosition;
    fOriMomentum = aTrackInfo->fOriMomentum;
    fOriEnergy = aTrackInfo->fOriEnergy;
    fOriTime = aTrackInfo->fOriTime;
    fAncestorMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TrackInformation::GetAncestor(G4int aDetectorID) const
{
    if (fAncestorMap.count(aDetectorID) > 0)
        return fAncestorMap.at(aDetectorID);

    return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackInformation::SetAncestor(G4int aDetectorID, G4int aTrackID)
{
    fAncestorMap.insert(std::pair<G4int, G4int>(aDetectorID, aTrackID));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
