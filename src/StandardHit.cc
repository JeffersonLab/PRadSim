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
// StandardHit.cc
// Developer : Jixie Zhang(original), Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StandardHit.hh"

#include "G4Allocator.hh"
#include "G4VHit.hh"

#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Allocator<StandardHit> StandardHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardHit::StandardHit() : G4VHit()
{
    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardHit::~StandardHit()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool StandardHit::operator ==(const StandardHit &right) const
{
    return ((fPID == right.fPID) && (fTrackID == right.fTrackID) && (fDetID == right.fDetID) && (fPhysV == right.fPhysV) && (fCopyNo == right.fCopyNo));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardHit::Print()
{
    G4int prec = G4cout.precision(3);
    G4cout << std::setw(10) << fPID << " " << std::setw(5) << fTrackID << " " << std::setw(5) << fPTrackID << " ";
    G4cout << std::setw(5) << G4BestUnit(fInPos.x(), "Length") << " " << std::setw(5) << G4BestUnit(fInPos.y(), "Length") << " " << std::setw(5) << G4BestUnit(fInPos.z(), "Length") << " ";
    G4cout << std::setw(5) << G4BestUnit(fOutPos.x(), "Length") << " " << std::setw(5) << G4BestUnit(fOutPos.y(), "Length") << " " << std::setw(5) << G4BestUnit(fOutPos.z(), "Length") << " ";
    G4cout << std::setw(5) << G4BestUnit(fInMom.mag(), "Energy") << " " << std::setw(8) << G4BestUnit(fOutMom.mag(), "Energy") << " " << std::setw(8) << G4BestUnit(fEdep, "Energy") << " ";
    G4cout << std::setw(5) << G4BestUnit(fTime, "Time") << " ";

    if (fPhysV != 0) G4cout << std::setw(18) << fPhysV->GetName() << " " << fCopyNo;

    G4cout << G4endl;
    G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardHit::Clear()
{
    fPID = 0;
    fTrackID = 0;
    fPTrackID = 0;
    fDetID = 0;
    fInPos.set(0, 0, 0);
    fOutPos.set(0, 0, 0);
    fInMom.set(0, 0, 0);
    fOutMom.set(0, 0, 0);
    fTime = 0;
    fEdep = 0;
    fTrackLen = 0;
    fPhysV = NULL;
    fCopyNo = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
