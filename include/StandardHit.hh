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
// StandardHit.hh
// Developer : Jixie Zhang(original), Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StandardHit_h
#define StandardHit_h 1

#include "G4VHit.hh"

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StandardHit : public G4VHit
{
public:
    StandardHit();
    virtual ~StandardHit();

    StandardHit(const StandardHit &right);
    const StandardHit &operator=(const StandardHit &right);

    G4bool operator==(const StandardHit &right) const;

    inline void *operator new (size_t);
    inline void operator delete (void *aHit);

    void Print();
    void Clear();

public:
    inline G4int GetPID() const;
    inline G4int GetTrackID() const;
    inline G4int GetParentTrackID() const;
    inline G4ThreeVector GetInPos() const;
    inline G4ThreeVector GetOutPos() const;
    inline G4ThreeVector GetInMom() const;
    inline G4ThreeVector GetOutMom() const;
    inline G4double GetTime() const;
    inline G4double GetEdep() const;
    inline const G4VPhysicalVolume *GetPhysV() const;
    inline G4int GetCopyNo() const;

    inline void SetPID(G4int &val);
    inline void SetTrackID(G4int &val);
    inline void SetParentTrackID(G4int &val);
    inline void SetInPos(G4ThreeVector &xyz);
    inline void SetOutPos(G4ThreeVector &xyz);
    inline void SetInMom(G4ThreeVector &pxpypz);
    inline void SetOutMom(G4ThreeVector &pxpypz);
    inline void SetTime(G4double &val);
    inline void SetEdep(G4double &val);
    inline void SetPhysV(G4VPhysicalVolume *val);
    inline void SetCopyNo(G4int &val);

    inline void AddEdep(G4double &val);

private:
    G4int         fPID;
    G4int         fTrackID;
    G4int         fPTrackID;
    G4ThreeVector fInPos;
    G4ThreeVector fOutPos;
    G4ThreeVector fInMom;
    G4ThreeVector fOutMom;
    G4double      fTime;
    G4double      fEdep;
    const G4VPhysicalVolume *fPhysV;
    G4int         fCopyNo;
};

typedef G4THitsCollection<StandardHit> StandardHitsCollection;

extern G4Allocator<StandardHit> StandardHitAllocator;

inline void *StandardHit::operator new (size_t)
{
    void *aHit;
    aHit = (void *)StandardHitAllocator.MallocSingle();
    return aHit;
}

inline void StandardHit::operator delete (void *aHit)
{
    StandardHitAllocator.FreeSingle((StandardHit *)aHit);
}

inline G4int StandardHit::GetPID() const
{
    return fPID;
}

inline G4int StandardHit::GetTrackID() const
{
    return fTrackID;
}

inline G4int StandardHit::GetParentTrackID() const
{
    return fPTrackID;
}

inline G4ThreeVector StandardHit::GetInPos() const
{
    return fInPos;
}

inline G4ThreeVector StandardHit::GetOutPos() const
{
    return fOutPos;
}

inline G4ThreeVector StandardHit::GetInMom() const
{
    return fInMom;
}

inline G4ThreeVector StandardHit::GetOutMom() const
{
    return fOutMom;
}

inline G4double StandardHit::GetTime() const
{
    return fTime;
}

inline G4double StandardHit::GetEdep() const
{
    return fEdep;
}

inline const G4VPhysicalVolume *StandardHit::GetPhysV() const
{
    return fPhysV;
}

inline G4int StandardHit::GetCopyNo() const
{
    return fCopyNo;
}

inline void StandardHit::SetPID(G4int &val)
{
    fPID = val;
}

inline void StandardHit::SetTrackID(G4int &val)
{
    fTrackID = val;
}

inline void StandardHit::SetParentTrackID(G4int &val)
{
    fPTrackID = val;
}

inline void StandardHit::SetInPos(G4ThreeVector &xyz)
{
    fInPos = xyz;
}

inline void StandardHit::SetOutPos(G4ThreeVector &xyz)
{
    fOutPos = xyz;
}

inline void StandardHit::SetInMom(G4ThreeVector &pxpypz)
{
    fInMom = pxpypz;
}

inline void StandardHit::SetOutMom(G4ThreeVector &pxpypz)
{
    fOutMom = pxpypz;
}

inline void StandardHit::SetTime(G4double &val)
{
    fTime = val;
}

inline void StandardHit::SetEdep(G4double &val)
{
    fEdep = val;
}

inline void StandardHit::SetPhysV(G4VPhysicalVolume *val)
{
    fPhysV = val;
}

inline void StandardHit::SetCopyNo(G4int &val)
{
    fCopyNo = val;
}

inline void StandardHit::AddEdep(G4double &val)
{
    fEdep += val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
