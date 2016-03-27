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
#include "Digitization.hh"
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
    G4String HCname;
    collectionName.insert(HCname="HitsCollection");
    gem_hits.reserve(100);
    daq_system = new Digitization();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
    HitsCollection = new CalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
    static G4int HCID = -1;
    if(HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    HCE->AddHitsCollection( HCID, HitsCollection );

    for(int i = 0; i < MAX_MODULE; ++i) ModuleEnergy[i] = 0;
    gem_hits.clear();
    TotalEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool CalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.)
        return false;

    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* PhysVol = theTouchable->GetVolume();

    G4String Name = PhysVol->GetName();
    if(Name == "GEM_Foil") {
        G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
//        if(edep/eV > 26) {
        // currently there is no gem hits reconstruction
        // thus the gem foil material is set to be vacuum
        // to avoid multiple hits count here
            double hitx = G4RandGauss::shoot(position.x()/mm, 0.1);
            double hity = G4RandGauss::shoot(position.y()/mm, 0.1);
            double hitz = PhysVol->GetTranslation().z()/mm;
            gem_hits.push_back(Digitization::GEM_Hit(hitx, hity, hitz));
//        }
    } else if(Name == "HyCal_Leadglass") {
        int ModuleID = PhysVol->GetCopyNo() + 1152;
        ModuleEnergy[ModuleID] += edep/MeV;
        TotalEnergy += edep/MeV;
    } else if(Name == "HyCal_Crystal") {
        int ModuleID = PhysVol->GetCopyNo();
        ModuleEnergy[ModuleID] += edep/MeV;
        TotalEnergy += edep/MeV;
    }

    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
    if(TotalEnergy < 500)
        return;

    daq_system->Event(ModuleEnergy, gem_hits);

    if(verboseLevel > 0) {
        G4int NbHits = HitsCollection->entries();
        G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
               << " hits in the Calorimeter: " << G4endl;
        for(int i = 0; i < NbHits; ++i)
            (*HitsCollection)[i]->Print();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

