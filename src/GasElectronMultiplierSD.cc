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
// $Id: GasElectronMultiplierSD.cc,2016-03-29 $
// GEANT4 tag $Name: geant4.10.02.p01 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "GasElectronMultiplierSD.hh"
#include "Digitization.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GasElectronMultiplierSD::GasElectronMultiplierSD(G4String name, Digitization *pdaq)
:G4VSensitiveDetector(name), daq_system(pdaq)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GasElectronMultiplierSD::~GasElectronMultiplierSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GasElectronMultiplierSD::Initialize(G4HCofThisEvent* /*HCE*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool GasElectronMultiplierSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.)
        return false;

    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* PhysVol = theTouchable->GetVolume();

    G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
    double hitx = G4RandGauss::shoot(position.x()/mm, 0.1);
    double hity = G4RandGauss::shoot(position.y()/mm, 0.1);
    double hitz = PhysVol->GetTranslation().z()/mm;

    daq_system->GEMHits(hitx, hity, hitz);
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GasElectronMultiplierSD::EndOfEvent(G4HCofThisEvent*)
{
    daq_system->Event(Digitization::GEM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

