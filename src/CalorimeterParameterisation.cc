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
// $Id: CalorimeterParameterisation.cc,2012-08-01 $
// GEANT4 tag $Name: geant-4-9-4-patch-02 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CalorimeterParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterParameterisation::CalorimeterParameterisation(  
        G4int    NoChambersX,
        G4int    NoChambersY, 
        G4double widthCalor,
        G4ThreeVector CenterPos,          //  center of the first 
        G4double sizeX,
        G4double sizeY) 
{
   fNoChambersX =  NoChambersX; 
   fNoChambersY =  NoChambersY; 
   fspacingXY = widthCalor;
   fCenterPos = CenterPos;
   fsizeX  =  sizeX;
   fsizeY  =  sizeY;
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

CalorimeterParameterisation::~CalorimeterParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

  G4int BlockIndex;
  G4int XIndex, YIndex;
  G4double Xposition, Yposition;

  //
  //Rule out the 4 blocks at the center
  //
    BlockIndex = copyNo;
    if (copyNo >= 560) BlockIndex = copyNo+2;
    if (copyNo >= 592) BlockIndex = copyNo+4;

  //
  //Calculate X position
  //
    XIndex = BlockIndex % fNoChambersX;
    Xposition = fCenterPos.x() + (XIndex - 16)*2.*fsizeX - fsizeX;

  //
  //Calculate Y position
  //
    YIndex = (BlockIndex - XIndex)/fNoChambersX;
    Yposition = fCenterPos.y() + (16 - YIndex)*2.*fsizeY + fsizeY;


//  G4cout << copyNo << "  " << Xposition << "  " << Yposition << G4endl;
  G4ThreeVector origin(Xposition, Yposition, 0.*cm);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CalorimeterParameterisation::ComputeDimensions
(G4Box& CalBlock, const G4int copyNo, const G4VPhysicalVolume*) const
{
/*
  if(copyNo >= 1152 ) {
    CalBlock.SetXHalfLength(1.9*cm);
    CalBlock.SetYHalfLength(1.9*cm);
  }
*/
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
