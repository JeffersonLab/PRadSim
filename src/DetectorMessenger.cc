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
//
// $Id: DetectorMessenger.cc, 2012-08-01 $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// Developer: Chao Peng
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction *Det) : Detector(Det)
{
    pradsimDir = new G4UIdirectory("/pradsim/");
    pradsimDir->SetGuidance("UI commands of this example");

    detDir = new G4UIdirectory("/pradsim/det/");
    detDir->SetGuidance("detector control");

    TargetMaterCmd = new G4UIcmdWithAString("/pradsim/det/setTargetMat", this);
    TargetMaterCmd->SetGuidance("Select Material of the Target.");
    TargetMaterCmd->SetParameterName("choice", false);
    TargetMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    UpdateCmd = new G4UIcmdWithoutParameter("/pradsim/det/update", this);
    UpdateCmd->SetGuidance("Update calorimeter geometry.");
    UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    UpdateCmd->SetGuidance("if you changed geometrical value(s).");
    UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
    delete TargetMaterCmd;
    delete UpdateCmd;
    delete detDir;
    delete pradsimDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == TargetMaterCmd)
        Detector->SetTargetMaterial(newValue);

    if (command == UpdateCmd)
        Detector->UpdateGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
