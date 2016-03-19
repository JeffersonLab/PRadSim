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
//
// $Id: CalorimeterHit.hh, 2012-08-03 $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterHit_h
#define CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterHit : public G4VHit
{
  public:

      CalorimeterHit();
     ~CalorimeterHit();
      CalorimeterHit(const CalorimeterHit&);
      const CalorimeterHit& operator=(const CalorimeterHit&);
      G4int operator==(const CalorimeterHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Printf(std::ostream &os);

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetStepE    (G4double se)      { stepe = se; };  
      void SetEdep     (G4double de)      { edep = de; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      void SetDirection(G4ThreeVector vec){ dir = vec; };
      void SetParName  (G4String name)    { parnam = name; };
      void SetVertex   (G4ThreeVector abc){ ver = abc; };
      void SetVolNumber(G4int number)     { voln = number; };

      G4int GetTrackID()    { return trackID; };
      G4double GetStepE()   { return stepe; };
      G4double GetEdep()    { return edep; };      
      G4ThreeVector GetPos(){ return pos; };
      G4ThreeVector GetDirection() {return dir; };
      G4String GetParName() {return parnam; };
      G4ThreeVector GetVertex() {return ver; }; 
      G4int GetVolNumber()  {return voln; };      
         
  private:
  
      G4int         trackID;
      G4double      stepe;
      G4double      edep;
      G4ThreeVector pos;
      G4ThreeVector dir;
      G4String      parnam;
      G4ThreeVector ver;
      G4int         voln;
      G4int flag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream &operator<<(std::ostream &os,CalorimeterHit &var);

typedef G4THitsCollection<CalorimeterHit> CalorimeterHitsCollection;

extern G4Allocator<CalorimeterHit> CalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CalorimeterHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) CalorimeterHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CalorimeterHit::operator delete(void *aHit)
{
  CalorimeterHitAllocator.FreeSingle((CalorimeterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
