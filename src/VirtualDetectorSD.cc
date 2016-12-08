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
// $Id: VirtualDetectorSD.cc,2016-11-20 $
// GEANT4 tag $Name: geant4.10.02.p01 $
// Developer: Chao Peng
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "VirtualDetectorSD.hh"

#include "RootTree.hh"

#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VirtualDetectorSD::VirtualDetectorSD(G4String name, RootTree *ptree) : G4VSensitiveDetector(name), otree(ptree)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VirtualDetectorSD::~VirtualDetectorSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VirtualDetectorSD::Initialize(G4HCofThisEvent * /*HCE*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool VirtualDetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    
    if (edep == 0.) return false;

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHistory *theTouchable = (G4TouchableHistory *)(preStepPoint->GetTouchable());
    G4VPhysicalVolume *PhysVol = theTouchable->GetVolume();

    G4Track *theTrack = aStep->GetTrack();  
    int pid = theTrack->GetParticleDefinition()->GetPDGEncoding();
    int tid = theTrack->GetTrackID();
    int parenttid = theTrack->GetParentID();
    
    G4ThreeVector position = preStepPoint->GetPosition();
    double hitx = position.x() / mm;
    double hity = position.y() / mm;
    //double hitz = PhysVol->GetTranslation().z() / mm;
    double hitz = position.z() / mm;
    
    G4ThreeVector momentum = preStepPoint->GetMomentum();
    double p = momentum.mag() / MeV;
    double theta = momentum.theta();
    double phi = momentum.phi();

    //G4cout << " " << pid << " " << tid << " " << parenttid << " " << hitx << " " << hity << " " << hitz << " " << p << " " << theta << " " << phi << G4endl;
    otree->UpdateValue(pid, tid, parenttid, hitx, hity, hitz, p, theta, phi);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VirtualDetectorSD::EndOfEvent(G4HCofThisEvent * /*HCE*/)
{
    otree->FillTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
