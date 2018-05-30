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
// StepRecordSD.cc
// Developer : Chao Gu
// History:
//   May 2018, C. Gu, Add for beam energy loss study.
//
// WARNING: do not use it for large density materials!
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepRecordSD.hh"

#include "GlobalVars.hh"
#include "RootTree.hh"
#include "StandardHit.hh"

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
#include "G4VProcess.hh"
#include "G4VSensitiveDetector.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepRecordSD::StepRecordSD(G4String name, G4String abbrev) : G4VSensitiveDetector(name), fAbbrev(abbrev), fRegistered(false)
{
    fID = name.hash() % 100000;
    //G4cout << name << "\t" << fAbbrev << "\t" << fID << G4endl;

    fN = 0;

    fPID.clear();
    fTID.clear();
    fPTID.clear();
    fX.clear();
    fY.clear();
    fZ.clear();
    fMomentum.clear();
    fTheta.clear();
    fPhi.clear();
    fProcessID.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepRecordSD::~StepRecordSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepRecordSD::Initialize(G4HCofThisEvent *)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StepRecordSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4Track *theTrack = aStep->GetTrack();
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    const G4VProcess *theProcess = preStepPoint->GetProcessDefinedStep();

    G4int PID = theTrack->GetParticleDefinition()->GetPDGEncoding();
    G4int TrackID = theTrack->GetTrackID();
    G4int ParentTrackID = theTrack->GetParentID();

    G4ThreeVector InPos = preStepPoint->GetPosition();
    G4ThreeVector InMom = preStepPoint->GetMomentum();

    G4int ProcessID = -1;

    if (theProcess)
        ProcessID = int(theProcess->GetProcessType()) * 1000 + theProcess->GetProcessSubType();

    fPID.push_back(PID);
    fTID.push_back(TrackID);
    fPTID.push_back(ParentTrackID);
    fX.push_back(InPos.x());
    fY.push_back(InPos.y());
    fZ.push_back(InPos.z());
    fMomentum.push_back(InMom.mag());
    fTheta.push_back(InMom.theta());
    fPhi.push_back(InMom.phi());
    fProcessID.push_back(ProcessID);

    fN++;

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepRecordSD::Register(TTree *tree)
{
    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.N", abbr), &fN, Form("%s.N/I", abbr));
    tree->Branch(Form("%s.PID", abbr), &fPID);
    tree->Branch(Form("%s.TID", abbr), &fTID);
    tree->Branch(Form("%s.PTID", abbr), &fPTID);
    tree->Branch(Form("%s.X", abbr), &fX);
    tree->Branch(Form("%s.Y", abbr), &fY);
    tree->Branch(Form("%s.Z", abbr), &fZ);
    tree->Branch(Form("%s.P", abbr), &fMomentum);
    tree->Branch(Form("%s.Theta", abbr), &fTheta);
    tree->Branch(Form("%s.Phi", abbr), &fPhi);
    tree->Branch(Form("%s.ProcessID", abbr), &fProcessID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepRecordSD::Clear()
{
    fPID.clear();
    fTID.clear();
    fPTID.clear();
    fX.clear();
    fY.clear();
    fZ.clear();
    fMomentum.clear();
    fTheta.clear();
    fPhi.clear();
    fProcessID.clear();

    fN = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
