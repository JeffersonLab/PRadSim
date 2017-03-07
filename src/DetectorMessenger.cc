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
// DetectorMessenger.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction *Det) : Detector(Det)
{
    PRadSimDir = new G4UIdirectory("/pradsim/");
    PRadSimDir->SetGuidance("UI commands of this example");

    DetDir = new G4UIdirectory("/pradsim/det/");
    DetDir->SetGuidance("Detector control");

    ZDir = new G4UIdirectory("/pradsim/det/z");
    ZDir->SetGuidance("Detector z coords control");

    TargetZCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/z/target", this);
    TargetZCmd->SetGuidance("Set fTargetCenter");
    TargetZCmd->SetParameterName("targetz", false);
    TargetZCmd->SetDefaultUnit("cm");

    GEM1ZCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/z/gem1", this);
    GEM1ZCmd->SetGuidance("Set fGEM1Center");
    GEM1ZCmd->SetParameterName("gem1z", false);
    GEM1ZCmd->SetDefaultUnit("cm");

    GEM2ZCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/z/gem2", this);
    GEM2ZCmd->SetGuidance("Set fGEM2Center");
    GEM2ZCmd->SetParameterName("gem2z", false);
    GEM2ZCmd->SetDefaultUnit("cm");

    SciPlaneZCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/z/plane", this);
    SciPlaneZCmd->SetGuidance("Set fSciPlaneCenter");
    SciPlaneZCmd->SetParameterName("planez", false);
    SciPlaneZCmd->SetDefaultUnit("cm");

    HyCalZCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/z/hycal", this);
    HyCalZCmd->SetGuidance("Set fCrystalSurf");
    HyCalZCmd->SetParameterName("hycalz", false);
    HyCalZCmd->SetDefaultUnit("cm");

    TargetDir = new G4UIdirectory("/pradsim/det/target");
    TargetDir->SetGuidance("Target control");

    TargetRCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/target/radius", this);
    TargetRCmd->SetGuidance("Set fTargetR");
    TargetRCmd->SetParameterName("targetr", false);
    TargetRCmd->SetDefaultUnit("mm");

    TargetHalfLCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/target/halfl", this);
    TargetHalfLCmd->SetGuidance("Set fTargetHalfL");
    TargetHalfLCmd->SetParameterName("targetl", false);
    TargetHalfLCmd->SetDefaultUnit("mm");

    RecoilDetDir = new G4UIdirectory("/pradsim/det/recoil");
    RecoilDetDir->SetGuidance("Recoil detector control");

    RecoilDetNSegCmd = new G4UIcmdWithAnInteger("/pradsim/det/recoil/nseg", this);
    RecoilDetNSegCmd->SetGuidance("Set fRecoilDetNSeg");
    RecoilDetNSegCmd->SetParameterName("recoiln", true);
    RecoilDetNSegCmd->SetDefaultValue(8);

    RecoilDetHalfLCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/recoil/halfl", this);
    RecoilDetHalfLCmd->SetGuidance("Set fRecoilDetHalfL");
    RecoilDetHalfLCmd->SetParameterName("recoill", false);
    RecoilDetHalfLCmd->SetDefaultUnit("mm");

    RecoilDetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/det/recoil/thick", this);
    RecoilDetThicknessCmd->SetGuidance("Set fRecoilDetThickness");
    RecoilDetThicknessCmd->SetParameterName("recoilt", false);
    RecoilDetThicknessCmd->SetDefaultUnit("mm");

    UpdateCmd = new G4UIcmdWithoutParameter("/pradsim/det/update", this);
    UpdateCmd->SetGuidance("Update calorimeter geometry.");
    UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    UpdateCmd->SetGuidance("if you changed geometrical value(s).");
    UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
    delete UpdateCmd;
    delete RecoilDetNSegCmd;
    delete RecoilDetHalfLCmd;
    delete RecoilDetThicknessCmd;
    delete RecoilDetDir;
    delete TargetDir;
    delete TargetRCmd;
    delete TargetHalfLCmd;
    delete TargetZCmd;
    delete GEM1ZCmd;
    delete GEM2ZCmd;
    delete SciPlaneZCmd;
    delete HyCalZCmd;
    delete DetDir;
    delete PRadSimDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == TargetZCmd)
        Detector->SetTargetPos(TargetZCmd->GetNewDoubleValue(newValue));

    if (command == GEM1ZCmd)
        Detector->SetGEMPos(GEM1ZCmd->GetNewDoubleValue(newValue), -10000);

    if (command == GEM2ZCmd)
        Detector->SetGEMPos(-10000, GEM2ZCmd->GetNewDoubleValue(newValue));

    if (command == SciPlaneZCmd)
        Detector->SetScitillatorPlanePos(SciPlaneZCmd->GetNewDoubleValue(newValue));

    if (command == HyCalZCmd)
        Detector->SetHyCalPos(HyCalZCmd->GetNewDoubleValue(newValue));

    if (command == TargetRCmd)
        Detector->SetTarget(TargetRCmd->GetNewDoubleValue(newValue), -10000);

    if (command == TargetHalfLCmd)
        Detector->SetTarget(-10000, TargetHalfLCmd->GetNewDoubleValue(newValue));

    if (command == RecoilDetNSegCmd)
        Detector->SetRecoilDetector(RecoilDetNSegCmd->GetNewIntValue(newValue), -10000, -10000);

    if (command == RecoilDetHalfLCmd)
        Detector->SetRecoilDetector(-10000, RecoilDetHalfLCmd->GetNewDoubleValue(newValue), -10000);

    if (command == RecoilDetThicknessCmd)
        Detector->SetRecoilDetector(-10000, -10000, RecoilDetThicknessCmd->GetNewDoubleValue(newValue));

    if (command == UpdateCmd)
        Detector->UpdateGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
