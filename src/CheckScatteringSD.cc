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
// CheckScatteringSD.cc
// Developer : Chao Gu
// History:
//   May 2018, C. Gu, Add for beam energy loss study.
//
// WARNING: do not use it for large density materials!
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CheckScatteringSD.hh"

#include "GlobalVars.hh"
#include "RootTree.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ProcessType.hh"
#include "G4VProcess.hh"
#include "G4VSensitiveDetector.hh"

#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CheckScatteringSD::CheckScatteringSD(G4String name, G4String abbrev) : G4VSensitiveDetector(name), fAbbrev(abbrev), fRegistered(false)
{
    fID = name.hash() % 100000;
    //G4cout << name << "\t" << fAbbrev << "\t" << fID << G4endl;

    fCoulomb = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CheckScatteringSD::~CheckScatteringSD()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CheckScatteringSD::Initialize(G4HCofThisEvent *)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CheckScatteringSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    const G4VProcess *Process = postStepPoint->GetProcessDefinedStep();

    if (Process) {
        if (Process->GetProcessType() == G4ProcessType::fElectromagnetic && Process->GetProcessSubType() == 1)
            fCoulomb = true;
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CheckScatteringSD::Register(TTree *tree)
{
    const char *abbr = fAbbrev.data();

    tree->Branch(Form("%s.Coulomb", abbr), &fCoulomb, Form("%s.Coulomb/O", abbr));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CheckScatteringSD::Clear()
{
    fCoulomb = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
