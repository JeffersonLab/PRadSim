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
#include "GlobalVars.hh"
#include "RootTree.hh"
#include "StandardDetectorSD.hh"
#include "StandardHit.hh"
#include "TrackInformation.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TMath.h"
#include "TTree.h"

#include "Math/Interpolator.h"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static const G4double ZCRFrontSurf = 273.1 * cm;
//static const G4double ZCRBackSurf = ZCRFrontSurf + 18.0 * cm;
static const G4double ZLGFrontSurf = 273.1 * cm - 9.73 * cm;
static const G4double ZLGBackSurf = ZLGFrontSurf + 45.0 * cm;
static const G4double sin_phi_c[100] = {0.01570732, 0.04710645, 0.07845910, 0.10973431, 0.14090123, 0.17192910, 0.20278730, 0.23344536, 0.26387305, 0.29404033, 0.32391742, 0.35347484, 0.38268343, 0.41151436, 0.43993917, 0.46792981, 0.49545867, 0.52249856, 0.54902282, 0.57500525, 0.60042023, 0.62524266, 0.64944805, 0.67301251, 0.69591280, 0.71812630, 0.73963109, 0.76040597, 0.78043041, 0.79968466, 0.81814972, 0.83580736, 0.85264016, 0.86863151, 0.88376563, 0.89802758, 0.91140328, 0.92387953, 0.93544403, 0.94608536, 0.95579301, 0.96455742, 0.97236992, 0.97922281, 0.98510933, 0.99002366, 0.99396096, 0.99691733, 0.99888987, 0.99987663, 0.99987663, 0.99888987, 0.99691733, 0.99396096, 0.99002366, 0.98510933, 0.97922281, 0.97236992, 0.96455742, 0.95579301, 0.94608536, 0.93544403, 0.92387953, 0.91140328, 0.89802758, 0.88376563, 0.86863151, 0.85264016, 0.83580736, 0.81814972, 0.79968466, 0.78043041, 0.76040597, 0.73963109, 0.71812630, 0.69591280, 0.67301251, 0.64944805, 0.62524266, 0.60042023, 0.57500525, 0.54902282, 0.52249856, 0.49545867, 0.46792981, 0.43993917, 0.41151436, 0.38268343, 0.35347484, 0.32391742, 0.29404033, 0.26387305, 0.23344536, 0.20278730, 0.17192910, 0.14090123, 0.10973431, 0.07845910, 0.04710645, 0.01570732};
static const G4double cos_phi_c[100] = {0.99987663, 0.99888987, 0.99691733, 0.99396096, 0.99002366, 0.98510933, 0.97922281, 0.97236992, 0.96455742, 0.95579301, 0.94608536, 0.93544403, 0.92387953, 0.91140328, 0.89802758, 0.88376563, 0.86863151, 0.85264016, 0.83580736, 0.81814972, 0.79968466, 0.78043041, 0.76040597, 0.73963109, 0.71812630, 0.69591280, 0.67301251, 0.64944805, 0.62524266, 0.60042023, 0.57500525, 0.54902282, 0.52249856, 0.49545867, 0.46792981, 0.43993917, 0.41151436, 0.38268343, 0.35347484, 0.32391742, 0.29404033, 0.26387305, 0.23344536, 0.20278730, 0.17192910, 0.14090123, 0.10973431, 0.07845910, 0.04710645, 0.01570732, -0.01570732, -0.04710645, -0.07845910, -0.10973431, -0.14090123, -0.17192910, -0.20278730, -0.23344536, -0.26387305, -0.29404033, -0.32391742, -0.35347484, -0.38268343, -0.41151436, -0.43993917, -0.46792981, -0.49545867, -0.52249856, -0.54902282, -0.57500525, -0.60042023, -0.62524266, -0.64944805, -0.67301251, -0.69591280, -0.71812630, -0.73963109, -0.76040597, -0.78043041, -0.79968466, -0.81814972, -0.83580736, -0.85264016, -0.86863151, -0.88376563, -0.89802758, -0.91140328, -0.92387953, -0.93544403, -0.94608536, -0.95579301, -0.96455742, -0.97236992, -0.97922281, -0.98510933, -0.99002366, -0.99396096, -0.99691733, -0.99888987, -0.99987663};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(G4String name, G4String abbrev, G4String pwo_filename) : StandardDetectorSD(name, abbrev)
{
    G4String cname = "ModuleColl";
    cname = fAbbrev + cname;
    collectionName.insert(cname);

    fAttenuationLG = 0.0;
    fReflectanceLG = 1.0;

    fTotalEdep = 0;
    fTotalTrackL = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 1e+38;
        fModuleTrackL[i] = 1e+38;
    }

    fInterpolator = new ROOT::Math::Interpolator(InterpolPoints, InterpolType);

    std::ifstream attenuation_file;
    attenuation_file.open(pwo_filename.c_str());
    G4double d[InterpolPoints], a[InterpolPoints];

    for (int i = 0; i < InterpolPoints - 1; i++)
        attenuation_file >> d[i] >> a[i];

    d[InterpolPoints - 1] = d[InterpolPoints - 2] + 1;
    a[InterpolPoints - 1] = a[InterpolPoints - 2];

    fInterpolator->SetData(InterpolPoints, d, a);
    attenuation_file.close();
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

    G4double InZ = (InPos.z() + OutPos.z()) / 2.0;
    G4double InBeta = preStepPoint->GetBeta();

    G4double Time = preStepPoint->GetGlobalTime();

    G4String theMat = thePhysVol->GetLogicalVolume()->GetMaterial()->GetName();

    G4double StepLength = 0;

    if (theTrack->GetParticleDefinition()->GetPDGCharge() != 0.) {
        if (theMat == "PbGlass") {
            // Edep = aStep->GetTotalEnergyDeposit();
            if (InBeta > 1.0 / 1.65) { // 1.65 is the index of Pb-glass
                //G4double theta_0 = InMom.theta() / 2.0;
                //G4double theta_c = acos(1.0 / (InBeta * 1.65));
                G4double cos_theta_c = 1.0 / (InBeta * 1.65);
                G4double sin_theta_c = sqrt(1.0 - cos_theta_c * cos_theta_c);

                G4ThreeVector direction_0 = InMom.unit();

                G4double distance = 0;
                G4int nreflect = 0;

                G4double factor = 0;

                for (G4int i = 0; i < 100; i++) {
                    G4ThreeVector direction_c(sin_theta_c * cos_phi_c[i], sin_theta_c * sin_phi_c[i], cos_theta_c);
                    direction_c.rotateUz(direction_0);

                    if (direction_c.z() > 0) {
                        distance = (ZLGBackSurf - InZ) / direction_c.z();
                        nreflect = int((fabs(distance * direction_c.x()) + 19.0) / 38.0) + int((fabs(distance * direction_c.y()) + 19.0) / 38.0);
                    } else if (direction_c.z() < 0) {
                        distance = -(InZ - ZLGFrontSurf + 450.0) / direction_c.z();
                        nreflect = int((fabs(distance * direction_c.x()) + 19.0) / 38.0) + int((fabs(distance * direction_c.y()) + 19.0) / 38.0) + 1;
                    } else
                        continue;

                    factor = factor + exp(-distance * fAttenuationLG) * pow(fReflectanceLG, nreflect);
                }

                factor = factor / 100.0;

                StepLength = aStep->GetStepLength() * factor * sin_theta_c * sin_theta_c;
            } // else beta < 1 / n, no Cherenkov light
        } else { // PbWO4
            G4double depth = InZ - ZCRFrontSurf;
            G4double factor = fInterpolator->Eval(depth);

            if (std::isnan(factor))
                G4cout << factor << " " << depth << G4endl;

            Edep = Edep * factor;
            StepLength = aStep->GetStepLength();
        }
    }

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
        aHit->SetDetectorID(DetectorID);
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

    CalorimeterHit *aCalorHit = (*fCalorHitsCollection)[DetectorID];

    aCalorHit->Add(Edep, StepLength);

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
    fTotalTrackL = 0;

    for (int i = 0; i < NModules; i++) {
        CalorimeterHit *aCalorHit = (*fCalorHitsCollection)[i];

        fModuleEdep[i] = aCalorHit->GetEdep();
        fModuleTrackL[i] = aCalorHit->GetTrackLength();
        fTotalEdep += fModuleEdep[i];
        fTotalTrackL += fModuleTrackL[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Register(TTree *tree)
{
    StandardDetectorSD::Register(tree);

    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.TotalEdep", abbr), &fTotalEdep, Form("%s.TotalEdep/D", abbr));
    tree->Branch(Form("%s.TotalTrackL", abbr), &fTotalTrackL, Form("%s.TotalTrackL/D", abbr));
    tree->Branch(Form("%s.ModuleEdep", abbr), fModuleEdep, Form("%s.ModuleEdep[%d]/D", abbr, NModules));
    tree->Branch(Form("%s.ModuleTrackL", abbr), fModuleTrackL, Form("%s.ModuleTrackL[%d]/D", abbr, NModules));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Clear()
{
    StandardDetectorSD::Clear();

    fTotalEdep = 0;
    fTotalTrackL = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 1e+38;
        fModuleTrackL[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
