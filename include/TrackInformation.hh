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
// TrackInformation.hh
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackInformation_h
#define TrackInformation_h 1

#include "G4VUserTrackInformation.hh"

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"

#include "G4ThreeVector.hh"

#include <unordered_map>

class G4Track;

typedef std::unordered_map<G4int, G4int> AncestorMap;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackInformation : public G4VUserTrackInformation
{
public:
    TrackInformation();
    TrackInformation(const G4Track *);
    TrackInformation(const TrackInformation *);
    virtual ~TrackInformation();

    inline void *operator new (size_t);
    inline void operator delete (void *);

    inline bool operator ==(const TrackInformation &) const;

    inline G4int GetOriginalTrackID() const;
    inline G4ParticleDefinition *GetOriginalParticle() const;
    inline G4ThreeVector GetOriginalPosition() const;
    inline G4ThreeVector GetOriginalMomentum() const;
    inline G4double GetOriginalEnergy() const;
    inline G4double GetOriginalTime() const;

    G4int GetAncestor(G4int aDetectorID) const;
    void SetAncestor(G4int aDetectorID, G4int aTrackID);

private:
    G4int                 fOriTrackID;
    G4ParticleDefinition *fOriParticle;
    G4ThreeVector         fOriPosition;
    G4ThreeVector         fOriMomentum;
    G4double              fOriEnergy;
    G4double              fOriTime;

    AncestorMap           fAncestorMap;
};

extern G4Allocator<TrackInformation> TrackInformationAllocator;

inline void *TrackInformation::operator new (size_t)
{
    void *aTrackInfo;
    aTrackInfo = (void *)TrackInformationAllocator.MallocSingle();
    return aTrackInfo;
}

inline void TrackInformation::operator delete (void *aTrackInfo)
{
    TrackInformationAllocator.FreeSingle((TrackInformation *)aTrackInfo);
}

inline G4bool TrackInformation::operator ==(const TrackInformation &right) const
{
    return (this == &right);
}

inline G4int TrackInformation::GetOriginalTrackID() const
{
    return fOriTrackID;
}

inline G4ParticleDefinition *TrackInformation::GetOriginalParticle() const
{
    return fOriParticle;
}

inline G4ThreeVector TrackInformation::GetOriginalPosition() const
{
    return fOriPosition;
}

inline G4ThreeVector TrackInformation::GetOriginalMomentum() const
{
    return fOriMomentum;
}

inline G4double TrackInformation::GetOriginalEnergy() const
{
    return fOriEnergy;
}

inline G4double TrackInformation::GetOriginalTime() const
{
    return fOriTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
