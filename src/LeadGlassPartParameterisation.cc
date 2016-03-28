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
// $Id: LeadGlassPartParameterisation.cc,2014-04-22 $
// GEANT4 tag $Name: geant-4-9-4-patch-02 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LeadGlassPartParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LeadGlassPartParameterisation::LeadGlassPartParameterisation(
        G4int    NoBlocksX,
        G4int    NoBlocksY,
        G4ThreeVector CenterPosition,          //  center of the first 
        G4double HalfsizeX,
        G4double HalfsizeY)
{
   fNoBlocksX =  NoBlocksX;
   fNoBlocksY =  NoBlocksY;
   fCenterPosition = CenterPosition;
   fHalfsizeX  =  HalfsizeX;
   fHalfsizeY  =  HalfsizeY;
/*
   if( NoChambersX > 0 && NoChambersY > 0 ){
     G4double widthCalorX = fNoChambersX * fsizeX;
     G4double widthCalorY = fNoChambersY * fsizeY;
     if ( (fspacingXY < widthCalorX) || (fspacingXY < widthCalorY) ) {
       G4Exception("CalorimeterParameterisation::CalorimeterParameterisation()",
                   "InvalidSetup", FatalException,
                   "Width>Spacing");
     }
   }
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LeadGlassPartParameterisation::~LeadGlassPartParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LeadGlassPartParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

  G4int BlockIndex, GroupIndex;
  G4int XIndex, YIndex;
  G4double Xposition = 0., Yposition = 0.;
  G4double XOffset, YOffset;

  GroupIndex = copyNo/144;
  BlockIndex = copyNo%144;

  XOffset = fNoBlocksX*fHalfsizeX - 34.85*cm;
  YOffset = fNoBlocksY*fHalfsizeY + 34.85*cm;
 
  YIndex = BlockIndex/24;
  XIndex = BlockIndex%24;

  if(GroupIndex == 0) {
    Xposition = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX + XOffset;
    Yposition = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY + YOffset;
  }

  if(GroupIndex == 2) {
    Xposition = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX - XOffset;
    Yposition = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY - YOffset;
  }

  if(GroupIndex == 1) {
    Yposition = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX + XOffset;
    Xposition = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY - YOffset;
  }

  if(GroupIndex == 3) {
    Yposition = -(fNoBlocksX - 1)*fHalfsizeX + XIndex*2.*fHalfsizeX - XOffset;
    Xposition = (fNoBlocksY - 1)*fHalfsizeY - YIndex*2.*fHalfsizeY + YOffset;
  }

  G4ThreeVector origin(Xposition, Yposition, 0.*cm);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void LeadGlassPartParameterisation::ComputeDimensions(G4Box& /*CalBlock*/,
                                                      const G4int /*copyNo*/,
                                                      const G4VPhysicalVolume* /*volume*/)
const
{
//Nothing
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
