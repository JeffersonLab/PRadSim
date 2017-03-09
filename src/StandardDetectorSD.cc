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
// StandardDetectorSD.cc
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StandardDetectorSD.hh"

#include "Globals.hh"
#include "RootTree.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSensitiveDetector.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDetectorSD::StandardDetectorSD(G4String name, G4String abbrev) : G4VSensitiveDetector(name), fAbbrev(abbrev)
{
    G4String cname = "Coll";
    cname = fAbbrev + cname;
    collectionName.insert(cname);

    fHCID = -1;
    fHitsCollection = 0;

    gRootTree->RegisterSD(fAbbrev.data());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDetectorSD::~StandardDetectorSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDetectorSD::Initialize(G4HCofThisEvent *HCE)
{
    fHitsCollection = new StandardHitsCollection(SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

    HCE->AddHitsCollection(fHCID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StandardDetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fHitsCollection) return false;

    G4double edep = aStep->GetTotalEnergyDeposit();

    if (edep <= 0.) return true;

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHistory *theTouchable = (G4TouchableHistory *)(preStepPoint->GetTouchable());
    G4VPhysicalVolume *thePhysVol = theTouchable->GetVolume();
    G4Track *theTrack = aStep->GetTrack();

    G4int PID = theTrack->GetParticleDefinition()->GetPDGEncoding();
    G4int TrackID = theTrack->GetTrackID();
    G4int PTrackID = theTrack->GetParentID();

    G4ThreeVector InPos = preStepPoint->GetPosition();
    G4ThreeVector InMom = preStepPoint->GetMomentum();

    G4double Time = preStepPoint->GetGlobalTime();

    G4String PhysName = thePhysVol->GetName();
    G4int CopyNo = theTouchable->GetVolume()->GetCopyNo();

    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    G4ThreeVector OutPos = postStepPoint->GetPosition();
    G4ThreeVector OutMom = postStepPoint->GetMomentum();

    StandardHit *aHit = 0;

    for (G4int i = fHitsCollection->entries() - 1; i >= 0; i--) {
        if ((*fHitsCollection)[i]->GetPhysV()->GetName() == PhysName && (*fHitsCollection)[i]->GetCopyNo() == CopyNo && (*fHitsCollection)[i]->GetTrackID() == TrackID) {
            // found an exist hit
            aHit = (*fHitsCollection)[i];
            break;
        }
    }

    // if found an exits hit, refresh it and accumulate the deposited energy
    // if not, create a new hit and push it into the collection
    if (aHit) {
        aHit->AddEdep(edep);

        if (aHit->GetTime() > Time) aHit->SetTime(Time);

        aHit->SetOutPos(OutPos);
        aHit->SetOutMom(OutMom);
    } else {
        // create a new hit
        aHit = new StandardHit();

        aHit->SetPID(PID);
        aHit->SetTrackID(TrackID);
        aHit->SetParentTrackID(PTrackID);
        aHit->SetInPos(InPos);
        aHit->SetInMom(InMom);
        aHit->SetOutPos(OutPos);
        aHit->SetOutMom(OutMom);
        aHit->SetTime(Time);
        aHit->SetEdep(edep);
        aHit->SetPhysV(thePhysVol);
        aHit->SetCopyNo(CopyNo);

        fHitsCollection->insert(aHit);
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDetectorSD::EndOfEvent(G4HCofThisEvent *HCE)
{
    if (!HCE) return; //no hits collection found

    // the above line just to avoid warning of not use HCE

    int NHitC = fHitsCollection->GetSize();

    if (NHitC <= 0) return;

    for (int i = 0; i < NHitC; i++) {
        StandardHit *aHit = (*fHitsCollection)[i];

        gRootTree->UpdateValue(fAbbrev, aHit->GetPID(), aHit->GetTrackID(), aHit->GetParentTrackID(), aHit->GetInPos().x(), aHit->GetInPos().y(), aHit->GetInPos().z(), aHit->GetInMom().mag(), aHit->GetOutPos().x(), aHit->GetOutPos().y(), aHit->GetOutPos().z(), aHit->GetOutMom().mag(), aHit->GetEdep(), aHit->GetTime(), aHit->GetCopyNo());
    }

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
