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
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StandardDetectorSD.hh"

#include "GlobalVars.hh"
#include "RootTree.hh"
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
#include "G4VSensitiveDetector.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDetectorSD::StandardDetectorSD(G4String name, G4String abbrev) : G4VSensitiveDetector(name), fAbbrev(abbrev), fHitsCollection(NULL), fRegistered(false)
{
    fID = name.hash() % 100000;
    //G4cout << name << "\t" << fAbbrev << "\t" << fID << G4endl;

    G4String cname = "Coll";
    cname = fAbbrev + cname;
    collectionName.insert(cname);

    fN = 0;

    for (int i = 0; i < MaxNHits; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fDID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
        fTime[i] = 1e+38;
        fEdep[i] = 1e+38;
        fTrackL[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDetectorSD::~StandardDetectorSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDetectorSD::Initialize(G4HCofThisEvent *HCE)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    fHitsCollection = new StandardHitsCollection(SensitiveDetectorName, collectionName[0]);

    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    HCE->AddHitsCollection(HCID, fHitsCollection);

    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StandardDetectorSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fHitsCollection) return false;

    G4Track *theTrack = aStep->GetTrack();
    TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());

    G4double Edep = aStep->GetTotalEnergyDeposit();

    G4int AncestorID = theTrackInfo->GetAncestor(fID);

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
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

    G4int CopyNo = theTouchable->GetCopyNumber();

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
        aHit->AddTrackLength(StepLength);

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
        aHit->SetTrackLength(StepLength);
        aHit->SetPhysV(thePhysVol);
        aHit->SetCopyNo(CopyNo);

        fHitsCollection->insert(aHit);
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

void StandardDetectorSD::EndOfEvent(G4HCofThisEvent *HCE)
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
        fEdep[i] = aHit->GetEdep();
        fTrackL[i] = aHit->GetTrackLength();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDetectorSD::Register(TTree *tree)
{
    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.N", abbr), &fN, Form("%s.N/I", abbr));
    tree->Branch(Form("%s.PID", abbr), fPID, Form("%s.PID[%s.N]/I", abbr, abbr));
    tree->Branch(Form("%s.TID", abbr), fTID, Form("%s.TID[%s.N]/I", abbr, abbr));
    tree->Branch(Form("%s.PTID", abbr), fPTID, Form("%s.PTID[%s.N]/I", abbr, abbr));
    tree->Branch(Form("%s.X", abbr), fX, Form("%s.X[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Y", abbr), fY, Form("%s.Y[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Z", abbr), fZ, Form("%s.Z[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.P", abbr), fMomentum, Form("%s.P[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Theta", abbr), fTheta, Form("%s.Theta[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Phi", abbr), fPhi, Form("%s.Phi[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Time", abbr), fTime, Form("%s.Time[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.Edep", abbr), fEdep, Form("%s.Edep[%s.N]/D", abbr, abbr));
    tree->Branch(Form("%s.TrackL", abbr), fTrackL, Form("%s.TrackL[%s.N]/D", abbr, abbr));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDetectorSD::Clear()
{
    for (int i = 0; i < fN; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fDID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
        fTime[i] = 1e+38;
        fEdep[i] = 1e+38;
        fTrackL[i] = 1e+38;
    }

    fN = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
