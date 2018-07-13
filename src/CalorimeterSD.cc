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
#include "TTree.h"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//static const G4double ZCRFrontSurf = 272.5 * cm;
//static const G4double ZCRBackSurf = ZCRFrontSurf + 18.0 * cm;
static const G4double ZLGFrontSurf = 272.5 * cm - 9.73 * cm;
static const G4double ZLGBackSurf = ZLGFrontSurf + 45.0 * cm;
static const G4double sin_phi_c[100] = {
    0.015707317311820675, 0.04710645070964266, 0.07845909572784494, 0.10973431109104527, 0.14090123193758267,
    0.17192910027940952, 0.2027872953565125, 0.23344536385590536, 0.2638730499653729, 0.29404032523230395,
    0.3239174181981494, 0.3534748437792571, 0.3826834323650898, 0.41151435860510877, 0.43993916985591514,
    0.46792981426057334, 0.4954586684324076, 0.5224985647159488, 0.5490228179981318, 0.5750052520432786,
    0.6004202253258839, 0.6252426563357051, 0.6494480483301837, 0.6730125135097733, 0.6959127965923143,
    0.7181262977631888, 0.7396310949786097, 0.7604059656000309, 0.7804304073383297, 0.7996846584870905,
    0.8181497174250234, 0.8358073613682702, 0.8526401643540922, 0.8686315144381912, 0.8837656300886935,
    0.8980275757606156, 0.9114032766354452, 0.9238795325112867, 0.9354440308298673, 0.9460853588275453,
    0.9557930147983301, 0.9645574184577981, 0.9723699203976766, 0.9792228106217657, 0.985109326154774,
    0.9900236577165575, 0.9939609554551797, 0.996917333733128, 0.99888987496197, 0.9998766324816606,
    0.9998766324816606, 0.99888987496197, 0.996917333733128, 0.9939609554551797, 0.9900236577165576,
    0.9851093261547739, 0.9792228106217657, 0.9723699203976766, 0.9645574184577981, 0.9557930147983302,
    0.9460853588275453, 0.9354440308298674, 0.9238795325112868, 0.9114032766354453, 0.8980275757606156,
    0.8837656300886935, 0.8686315144381912, 0.8526401643540923, 0.8358073613682702, 0.8181497174250234,
    0.7996846584870907, 0.7804304073383297, 0.760405965600031, 0.7396310949786096, 0.7181262977631888,
    0.6959127965923144, 0.6730125135097733, 0.6494480483301838, 0.6252426563357054, 0.600420225325884,
    0.5750052520432787, 0.549022817998132, 0.5224985647159494, 0.4954586684324074, 0.4679298142605734,
    0.4399391698559153, 0.4115143586051091, 0.3826834323650903, 0.353474843779257, 0.3239174181981494,
    0.2940403252323041, 0.26387304996537325, 0.2334453638559055, 0.20278729535651233, 0.17192910027940955,
    0.14090123193758286, 0.10973431109104566, 0.07845909572784507, 0.047106450709642964, 0.015707317311820707
};
static const G4double cos_phi_c[100] = {
    0.9998766324816606, 0.99888987496197, 0.996917333733128, 0.9939609554551797, 0.9900236577165575,
    0.9851093261547739, 0.9792228106217657, 0.9723699203976766, 0.9645574184577981, 0.9557930147983301,
    0.9460853588275453, 0.9354440308298674, 0.9238795325112867, 0.9114032766354453, 0.8980275757606156,
    0.8837656300886935, 0.8686315144381912, 0.8526401643540923, 0.8358073613682702, 0.8181497174250234,
    0.7996846584870906, 0.7804304073383298, 0.760405965600031, 0.7396310949786097, 0.7181262977631888,
    0.6959127965923144, 0.6730125135097733, 0.6494480483301838, 0.6252426563357053, 0.6004202253258841,
    0.5750052520432786, 0.5490228179981318, 0.5224985647159489, 0.4954586684324076, 0.46792981426057334,
    0.4399391698559151, 0.4115143586051089, 0.38268343236508984, 0.35347484377925714, 0.32391741819814934,
    0.29404032523230406, 0.26387304996537275, 0.23344536385590545, 0.20278729535651271, 0.1719291002794095,
    0.1409012319375828, 0.10973431109104537, 0.078459095727845, 0.04710645070964268, 0.015707317311820648,
    -0.015707317311820523, -0.047106450709642554, -0.07845909572784489, -0.10973431109104526, -0.14090123193758247,
    -0.1719291002794096, -0.20278729535651238, -0.23344536385590556, -0.26387304996537286, -0.2940403252323038,
    -0.32391741819814945, -0.35347484377925703, -0.3826834323650895, -0.41151435860510877, -0.43993916985591514,
    -0.46792981426057323, -0.4954586684324076, -0.5224985647159488, -0.5490228179981318, -0.5750052520432786,
    -0.6004202253258839, -0.6252426563357052, -0.6494480483301835, -0.6730125135097734, -0.6959127965923143,
    -0.7181262977631887, -0.7396310949786098, -0.7604059656000309, -0.7804304073383296, -0.7996846584870906,
    -0.8181497174250233, -0.83580736136827, -0.852640164354092, -0.8686315144381913, -0.8837656300886934,
    -0.8980275757606155, -0.9114032766354451, -0.9238795325112865, -0.9354440308298674, -0.9460853588275453,
    -0.9557930147983301, -0.964557418457798, -0.9723699203976766, -0.9792228106217659, -0.9851093261547739,
    -0.9900236577165575, -0.9939609554551796, -0.996917333733128, -0.9988898749619699, -0.9998766324816606
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(G4String name, G4String abbrev) : StandardDetectorSD(name, abbrev)
{
    G4String cname = "ModuleColl";
    cname = fAbbrev + cname;
    collectionName.insert(cname);

    fAttenuationCR = 1.0e-10;
    fAttenuationLG = 1.0e-10;

    fReflectance = 1.0;

    fTotalEdep = 0;
    fTotalTrackL = 0;

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
                    G4double cos_theta = direction_c.cosTheta();
                    G4double sin_theta = sqrt(1 - cos_theta * cos_theta);

                    if (cos_theta > 0) {
                        distance = (ZLGBackSurf - InZ) / cos_theta;
                        nreflect = int(distance * sin_theta / 38.0);
                    } else if (cos_theta < 0) {
                        distance = -(InZ - ZLGFrontSurf + 450.0) / cos_theta;
                        nreflect = int(distance * sin_theta / 38.0) + 1;
                    } else
                        continue;

                    factor = factor + exp(-distance * fAttenuationLG) * pow(fReflectance, nreflect);
                }

                factor = factor / 100.0;

                StepLength = aStep->GetStepLength() * factor * sin_theta_c * sin_theta_c;
            } // else beta < 1 / n, no Cherenkov light
        } else { // PbWO4
            // Edep = aStep->GetTotalEnergyDeposit();

            // G4double distance = 0;
            // G4int nreflect = 0;

            // G4double factor = 0;

            // for (G4int i = 0; i < 101; i++) {
            //     G4double cos_theta = 1.0 - 2.0 * i / 100.0;
            //     G4double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            //     if (cos_theta > 0) {
            //         distance = (ZCRBackSurf - InZ) / cos_theta;
            //         nreflect = int(distance * sin_theta / 20.5);
            //     } else if (cos_theta < 0) {
            //         distance = -(InZ - ZCRFrontSurf + 180.0) / cos_theta;
            //         nreflect = int(distance * sin_theta / 20.5) + 1;
            //     } else
            //         continue;

            //     factor = factor + exp(-distance * fAttenuationCR) * pow(fReflectance, nreflect);
            // }

            // factor = factor / 101.0;

            // Edep = Edep * factor;

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
