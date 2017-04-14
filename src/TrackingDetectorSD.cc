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
// TrackingDetectorSD.cc
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingDetectorSD.hh"

#include "Globals.hh"
#include "RootTree.hh"
#include "StandardDetectorSD.hh"
#include "StandardHit.hh"
#include "TrackInformation.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingDetectorSD::TrackingDetectorSD(G4String name, G4String abbrev) : StandardDetectorSD(name, abbrev)
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingDetectorSD::~TrackingDetectorSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackingDetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fHitsCollection) return false;

    G4Track *theTrack = aStep->GetTrack();
    TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());

    G4double Edep = aStep->GetTotalEnergyDeposit();

    G4int AncestorID = theTrackInfo->GetAncestor(fID);

    if (Edep > 0) {
        G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
        G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
        G4VPhysicalVolume *thePhysVol = theTouchable->GetVolume();

        G4int PID = theTrack->GetParticleDefinition()->GetPDGEncoding();
        G4int TrackID = theTrack->GetTrackID();
        G4int ParentTrackID = theTrack->GetParentID();

        G4int DetectorID = 0;

        for (G4int i = 0; i < theTouchable->GetHistoryDepth(); i++)
            DetectorID += theTouchable->GetCopyNumber(i);

        G4ThreeVector InPos = preStepPoint->GetPosition();
        G4ThreeVector InMom = preStepPoint->GetMomentum();

        G4ThreeVector OutPos = postStepPoint->GetPosition();
        G4ThreeVector OutMom = postStepPoint->GetMomentum();

        G4double Time = preStepPoint->GetGlobalTime();

        G4int CopyNo = theTouchable->GetCopyNumber();

        if (AncestorID < 0) AncestorID = TrackID;

        StandardHit *aHit = NULL;

        for (G4int i = fHitsCollection->entries() - 1; i >= 0; i--) {
            if ((*fHitsCollection)[i]->GetDetectorID() == DetectorID && (*fHitsCollection)[i]->GetTrackID() == AncestorID) {
                aHit = (*fHitsCollection)[i];
                break;
            }
        }

        // if found an exits hit, refresh it and accumulate the deposited energy
        // if not, create a new hit and push it into the collection
        if (aHit) {
            aHit->AddEdep(Edep);

            if (aHit->GetTrackID() == TrackID) {
                if (aHit->GetTime() > Time) aHit->SetTime(Time);

                aHit->SetOutPos(OutPos);
                aHit->SetOutMom(OutMom);
            }
        } else {
            // create a new hit
            aHit = new StandardHit();

            aHit->SetPID(PID);
            aHit->SetTrackID(TrackID);
            aHit->SetParentTrackID(ParentTrackID);
            aHit->SetDetectorID(DetectorID);
            aHit->SetInPos(InPos);
            aHit->SetInMom(InMom);
            aHit->SetOutPos(OutPos);
            aHit->SetOutMom(OutMom);
            aHit->SetTime(Time);
            aHit->SetEdep(Edep);
            aHit->SetPhysV(thePhysVol);
            aHit->SetCopyNo(CopyNo);

            fHitsCollection->insert(aHit);
        }
    }

    G4int nSecondaries = aStep->GetNumberOfSecondariesInCurrentStep();

    if (nSecondaries > 0 && AncestorID >= 0) {
        for (auto &aSecondary : * (aStep->GetSecondaryInCurrentStep())) {
            if (aSecondary->GetUserInformation() == 0) {
                TrackInformation *newTrackInfo = new TrackInformation(theTrackInfo);
                newTrackInfo->SetAncestor(fID, AncestorID);
                G4Track *theSecondary = (G4Track *)aSecondary;
                theSecondary->SetUserInformation(newTrackInfo);
            }
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingDetectorSD::Register(TTree *tree)
{
    StandardDetectorSD::Register(tree);

    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.DID", abbr), fDID, Form("%s.DID[%s.N]/I", abbr, abbr));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
