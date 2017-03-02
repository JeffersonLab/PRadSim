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
// PrimaryGeneratorMessenger.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *Gun) : Action(Gun)
{
    gunDir = new G4UIdirectory("/pradsim/gun/");
    gunDir->SetGuidance("PrimaryGenerator control");

    RandCmd = new G4UIcmdWithAString("/pradsim/gun/random", this);
    RandCmd->SetGuidance("Shoot randomly the incident particle.");
    RandCmd->SetGuidance("  Choice : on(default), off");
    RandCmd->SetParameterName("random", true);
    RandCmd->SetDefaultValue("on");
    RandCmd->SetCandidates("on off");
    RandCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    GunTypeCmd = new G4UIcmdWithAString("/pradsim/gun/type", this);
    GunTypeCmd->SetGuidance("Choose a type of event generator.");
    GunTypeCmd->SetGuidance("  Choice : ring (default), elastic, moller");
    GunTypeCmd->SetParameterName("generator", true);
    GunTypeCmd->SetDefaultValue("ring");
    GunTypeCmd->SetCandidates("ring elastic moller");
    GunTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    StartEventCmd = new G4UIcmdWithAnInteger("/pradsim/gun/startEvent", this);
    StartEventCmd->SetGuidance("Set start point for event gun file.");
    StartEventCmd->SetParameterName("nstart", true);
    StartEventCmd->SetDefaultValue(0);
    StartEventCmd->SetRange("nstart>=0");
    StartEventCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete RandCmd;
    delete GunTypeCmd;
    delete StartEventCmd;
    delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
    G4UIcommand *command, G4String newValue)
{
    if (command == RandCmd)
        Action->SetRandFlag(newValue);

    if (command == GunTypeCmd)
        Action->SetGunType(newValue);

    if (command == StartEventCmd)
        Action->SetStartEvent(StartEventCmd->GetNewIntValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
