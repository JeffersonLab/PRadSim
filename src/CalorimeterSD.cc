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
// $Id: CalorimeterSD.cc,2012-08-01 $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "CalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(G4String name)
:G4VSensitiveDetector(name)
{
  HyCaldata.open("EnergyDeposit.dat");
  GEMdata.open("GEMPosition.dat");
  G4String HCname;
  collectionName.insert(HCname="HitsCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD(){
  HyCaldata.close();
  GEMdata.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
  HitsCollection = new CalorimeterHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID < 0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, HitsCollection ); 

  for(G4int i = 0; i < 1728; ++i) {EnergyDeposit[i] = 0.;}
  for(G4int i = 0; i < 100; ++i) {GEM_x[i] = 0.; GEM_y[i] = 0.; GEM_n = 0;}// TID[i] = 0; GEM_z[i] = 0.;}
  TrackID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool CalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* PhysVol = theTouchable->GetVolume();

  G4String Name = PhysVol->GetName();
//  G4cout << Name << G4endl;
  if(Name == "GEM_Foil"){
    position = aStep->GetPreStepPoint()->GetPosition();
//    position = aStep->GetTrack()->GetPosition();
    G4double KinEne = aStep->GetTrack()->GetKineticEnergy();
//    G4int TrackNow = aStep->GetTrack()->GetTrackID();
//    if(edep/eV > 11.5 && TrackNow != TrackID && GEM_n < 100 && (position.z() == 2220. || position.z() == 2240.)) {
   if(KinEne > 0.) {
      GEM_x[GEM_n] = G4RandGauss::shoot(position.x()/mm, 0.1)*mm;
      GEM_y[GEM_n] = G4RandGauss::shoot(position.y()/mm, 0.1)*mm;
      GEM_z[GEM_n] = position.z();
      GEM_E[GEM_n] = KinEne;
      GEM_n += 1;
//      GEMdata << position.x()/mm << "  " << position.y()/mm << "  " << position.z()/mm << "  " << KinEne << G4endl;
//      TID[GEM_n] = aStep->GetTrack()->GetParentID();
//      TrackID = TrackNow;
    }
  }
  if(Name == "HyCal_Leadglass"){
    G4int ModuleID = PhysVol->GetCopyNo() + 1152;
    EnergyDeposit[ModuleID] = EnergyDeposit[ModuleID] + edep;
  }
  if(Name == "HyCal_Crystal"){
    G4int ModuleID = PhysVol->GetCopyNo();
    EnergyDeposit[ModuleID] = EnergyDeposit[ModuleID] + edep;
  }
//     newHit->Draw();
  return true;
} 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
//  datafile.open("EnergyDeposit.dat");

  G4double TotEne = 0.;

  for(G4int i = 0; i < 1728; ++i) {
//    TotEne += EnergyDeposit[i];
    if(EnergyDeposit[i]/MeV > 1.) {
//      datafile << i+1 << "  " << EnergyDeposit[i]/MeV << G4endl;
      TotEne += EnergyDeposit[i];
    }
  }

  if(TotEne/MeV > 500.)
  {
    for(G4int i = 0; i < 1728; ++i) {
     if(EnergyDeposit[i]/MeV > 1.) HyCaldata << i+1 << "  " << EnergyDeposit[i]/MeV << G4endl;
    }
    HyCaldata << 0 << "  " << TotEne << G4endl;
//    G4cout << GEM_n << G4endl;
    for(G4int i = 0; i < GEM_n; ++i) {
      GEMdata << GEM_x[i] << "  " << GEM_y[i] << "  " << GEM_z[i] << G4endl; //"  " << GEM_z[i] << "  "<< TID[i] << G4endl;
    }
    GEMdata << 0 << "  " << GEM_n << "  " << 0 << G4endl;
  }
//  datafile.close();

//   if(GEM_n > 0) GEMdata << 0 << "  " << 0 << "  " << 0 << "  " << GEM_n << G4endl;

  if (verboseLevel>0) {
     G4int NbHits = HitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
            << " hits in the Calorimeter: " << G4endl;
     for (G4int i=0;i<NbHits;i++) (*HitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

