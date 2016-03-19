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
// $Id: CalorimeterHit.cc, 2012-08-03 $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"


G4Allocator<CalorimeterHit> CalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::CalorimeterHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::~CalorimeterHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::CalorimeterHit(const CalorimeterHit& right)
  : G4VHit()
{
  parnam    = right.parnam;
  trackID   = right.trackID;
  stepe     = right.stepe;
  edep      = right.edep;
  pos       = right.pos;
  dir       = right.dir;
  ver       = right.ver;
  voln      = right.voln;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const CalorimeterHit& CalorimeterHit::operator=(const CalorimeterHit& right)
{
  parnam    = right.parnam;
  trackID   = right.trackID;
  stepe     = right.stepe;
  edep      = right.edep;
  pos       = right.pos;
  dir       = right.dir;
  ver       = right.ver;
  voln      = right.voln;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int CalorimeterHit::operator==(const CalorimeterHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterHit::Printf(std::ostream & os)
{
/*
 if( (ver.x()*ver.x()+ver.y()*ver.y())/0.4*0.4*mm*mm >1)
 {flag = 1;}
 else
 {flag = 0;}
*/
/*
 BolckNo = voln;
 if (voln >= 560) BlockNo = voln+2;
 if (voln >= 592) BlockNo = voln+4;
*/
/*
os << pos.x() <<" "<< pos.y()<<" " << pos.z()<<" ";
os << stepe <<"     ";
os << dir.x() <<" "<< dir.y()<<" " << dir.z()<<" ";
os << voln;
os << G4endl;
*/
}

std::ostream &operator<<(std::ostream &os, CalorimeterHit &hit)
{
  hit.Printf(os);
  return os;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

