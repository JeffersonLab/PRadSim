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
// CalorimeterSD.cc
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CalorimeterSD.hh"

#include "CalorimeterHit.hh"
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
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(G4String name, G4String abbrev) : StandardDetectorSD(name, abbrev)
{
    G4String cname = "ModuleColl";
    cname = fAbbrev + cname;
    collectionName.insert(cname);

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 1e+38;
        fModuleTrackL[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent *HCE)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    fHitsCollection = new StandardHitsCollection(SensitiveDetectorName, collectionName[0]);
    fCalorHitsCollection = new CalorimeterHitsCollection(SensitiveDetectorName, collectionName[1]);

    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    HCE->AddHitsCollection(HCID, fHitsCollection);
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(fCalorHitsCollection);
    HCE->AddHitsCollection(HCID, fCalorHitsCollection);

    for (G4int i = 0; i < NModules; i++)
        fCalorHitsCollection->insert(new CalorimeterHit());

    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fHitsCollection || !fCalorHitsCollection) return false;

    G4Track *theTrack = aStep->GetTrack();
    TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());

    G4double Edep = aStep->GetTotalEnergyDeposit();

    G4int AncestorID = theTrackInfo->GetAncestor(fID);

    if (Edep > 0) {
        G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
        G4TouchableHistory *theTouchable = (G4TouchableHistory *)(preStepPoint->GetTouchable());
        G4VPhysicalVolume *thePhysVol = theTouchable->GetVolume();

        G4int PID = theTrack->GetParticleDefinition()->GetPDGEncoding();
        G4int TrackID = theTrack->GetTrackID();
        G4int ParentTrackID = theTrack->GetParentID();

        G4ThreeVector InPos = preStepPoint->GetPosition();
        G4ThreeVector InMom = preStepPoint->GetMomentum();

        G4ThreeVector OutPos = postStepPoint->GetPosition();
        G4ThreeVector OutMom = postStepPoint->GetMomentum();

        G4double Time = preStepPoint->GetGlobalTime();

        G4double StepLength = 0;

        if (theTrack->GetParticleDefinition()->GetPDGCharge() != 0.)
            StepLength = aStep->GetStepLength();

        G4int CopyNo = theTouchable->GetVolume()->GetCopyNo();

        if (AncestorID < 0) AncestorID = TrackID;

        StandardHit *aHit = NULL;

        for (G4int i = fHitsCollection->entries() - 1; i >= 0; i--) {
            if ((*fHitsCollection)[i]->GetTrackID() == AncestorID) {
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

        CalorimeterHit *aCalorHit = (*fCalorHitsCollection)[CopyNo];

        aCalorHit->Add(Edep, StepLength);
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

void CalorimeterSD::EndOfEvent(G4HCofThisEvent *HCE)
{
    if (!HCE) return; //no hits collection found

    int NHits = fHitsCollection->GetSize();

    if (NHits <= 0) return;

    fN = NHits;

    if (fN > MaxNHits) {
        G4cout << "WARNING: " << fN << " hits in " << fHitsCollection->GetName() << " exceed " << MaxNHits << G4endl;
        fN = MaxNHits;
    }

    for (int i = 0; i < fN; i++) {
        StandardHit *aHit = (*fHitsCollection)[i];

        fPID[i] = aHit->GetPID();
        fTID[i] = aHit->GetTrackID();
        fPTID[i] = aHit->GetParentTrackID();
        fDID[i] = aHit->GetDetectorID();
        fX[i] = aHit->GetInPos().x();
        fY[i] = aHit->GetInPos().y();
        fZ[i] = aHit->GetInPos().z();
        fMomentum[i] = aHit->GetInMom().mag();
        fTheta[i] = aHit->GetInMom().theta();
        fPhi[i] = aHit->GetInMom().phi();
        fTime[i] = aHit->GetTime();
    }

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        CalorimeterHit *aCalorHit = (*fCalorHitsCollection)[i];

        fModuleEdep[i] = aCalorHit->GetEdep();
        fModuleTrackL[i] = aCalorHit->GetTrackLength();
        fTotalEdep += fModuleEdep[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Register(TTree *tree)
{
    StandardDetectorSD::Register(tree);

    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.TotalEdep", abbr), &fTotalEdep, Form("%s.TotalEdep/D", abbr));
    tree->Branch(Form("%s.ModuleEdep", abbr), fModuleEdep, Form("%s.ModuleEdep[%d]/D", abbr, NModules));
    tree->Branch(Form("%s.ModuleTrackL", abbr), fModuleTrackL, Form("%s.ModuleTrackL[%d]/D", abbr, NModules));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Clear()
{
    StandardDetectorSD::Clear();

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 1e+38;
        fModuleTrackL[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
