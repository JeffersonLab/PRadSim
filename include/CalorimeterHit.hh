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
// CalorimeterHit.hh
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterHit_h
#define CalorimeterHit_h 1

#include "G4VHit.hh"

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterHit : public G4VHit
{
public:
    CalorimeterHit();
    virtual ~CalorimeterHit();

    inline void *operator new (size_t);
    inline void operator delete (void *);

    bool operator ==(const CalorimeterHit &) const;

    void Print();

    G4double GetEdep() const;
    G4double GetTrackLength() const;

    inline void Add(G4double &de, G4double &dl);

private:
    G4double fEdep;
    G4double fTrackLen;
};

typedef G4THitsCollection<CalorimeterHit> CalorimeterHitsCollection;

extern G4Allocator<CalorimeterHit> CalorimeterHitAllocator;

inline void *CalorimeterHit::operator new (size_t)
{
    void *aHit;
    aHit = (void *)CalorimeterHitAllocator.MallocSingle();
    return aHit;
}

inline void CalorimeterHit::operator delete (void *aHit)
{
    CalorimeterHitAllocator.FreeSingle((CalorimeterHit *)aHit);
}

inline G4double CalorimeterHit::GetEdep() const
{
    return fEdep;
}

inline G4double CalorimeterHit::GetTrackLength() const
{
    return fTrackLen;
}

inline void CalorimeterHit::Add(G4double &de, G4double &dl)
{
    fEdep += de;
    fTrackLen += dl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
