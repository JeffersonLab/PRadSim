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
// TrackingAction.cc
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Rewrite sensitive detectors.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "GlobalVars.hh"
#include "RootTree.hh"
#include "TrackInformation.hh"
#include "TrackingMessenger.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"
#include "G4UserTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() : G4UserTrackingAction(), fNoSecondary(false), fSaveTrackInfo(false), fRegistered(false)
{
    trackingMessenger = new TrackingMessenger(this);

    fN = 0;

    for (int i = 0; i < MaxTracks; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fProcessID[i] = -9999;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track *aTrack)
{
    if (aTrack->GetParentID() != 0 && fNoSecondary) {
        fpTrackingManager->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
        return;
    }

    if (aTrack->GetParentID() == 0 && aTrack->GetUserInformation() == 0) {
        TrackInformation *aTrackInfo = new TrackInformation(aTrack);
        G4Track *theTrack = (G4Track *)aTrack;
        theTrack->SetUserInformation(aTrackInfo);
    }

    if (fSaveTrackInfo) {
        if (!fRegistered) {
            Register(gRootTree->GetTree());
            fRegistered = true;
        }

        const G4LogicalVolume *theLV = aTrack->GetLogicalVolumeAtVertex();
        const G4VProcess *theProcess = aTrack->GetCreatorProcess();

        if (theLV && theProcess) {
            G4String theLVName = theLV->GetName();
            G4ProcessType theProcessType = theProcess->GetProcessType();
            int theProcessSubType = theProcess->GetProcessSubType();
            const G4ThreeVector thePosition = aTrack->GetVertexPosition();

            if (fN < MaxTracks && theLVName != "PbWO4AbsorberLV" && theLVName != "PbGlassAbsorberLV") {
                fPID[fN] = aTrack->GetParticleDefinition()->GetPDGEncoding();
                fTID[fN] = aTrack->GetTrackID();
                fPTID[fN] = aTrack->GetParentID();
                fX[fN] = thePosition.x();
                fY[fN] = thePosition.y();
                fZ[fN] = thePosition.z();
                fProcessID[fN] = (int)theProcessType * 1000 + theProcessSubType;
                fN++;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track *aTrack)
{
    G4TrackVector *theSecondaries = fpTrackingManager->GimmeSecondaries();

    if (theSecondaries) {
        TrackInformation *theTrackInfo = (TrackInformation *)(aTrack->GetUserInformation());
        size_t nSecondaries = theSecondaries->size();

        if (nSecondaries > 0) {
            for (size_t i = 0; i < nSecondaries; i++) {
                if ((*theSecondaries)[i]->GetUserInformation() == 0) {
                    TrackInformation *newTrackInfo = new TrackInformation(theTrackInfo);
                    (*theSecondaries)[i]->SetUserInformation(newTrackInfo);
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::Register(TTree *tree)
{
    tree->Branch("TR.N", &fN, "TR.N/I");
    tree->Branch("TR.PID", fPID, "TR.PID[TR.N]/I");
    tree->Branch("TR.TID", fTID, "TR.TID[TR.N]/I");
    tree->Branch("TR.PTID", fPTID, "TR.PTID[TR.N]/I");
    tree->Branch("TR.X", fX, "TR.X[TR.N]/D");
    tree->Branch("TR.Y", fY, "TR.Y[TR.N]/D");
    tree->Branch("TR.Z", fZ, "TR.Z[TR.N]/D");
    tree->Branch("TR.ProcessID", fProcessID, "TR.ProcessID[TR.N]/I");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::Clear()
{
    for (int i = 0; i < fN; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fProcessID[i] = -9999;
    }

    fN = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
