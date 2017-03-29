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

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *act) : G4UImessenger(), Action(act)
{
    GunDir = new G4UIdirectory("/pradsim/gun/");
    GunDir->SetGuidance("Primary generator control");

    GunTypeCmd = new G4UIcmdWithAString("/pradsim/gun/type", this);
    GunTypeCmd->SetGuidance("Choose a type of event generator.");
    GunTypeCmd->SetGuidance("  Choice : ring, elastic, moller");
    GunTypeCmd->SetParameterName("generator", false);
    GunTypeCmd->SetCandidates("ring elastic moller");

    RecoilCmd = new G4UIcmdWithAString("/pradsim/gun/recoil", this);
    RecoilCmd->SetGuidance("Choose a type of recoil particle.");
    RecoilCmd->SetGuidance("  Choice : none, proton, deuteron");
    RecoilCmd->SetParameterName("recoil", false);
    RecoilCmd->SetCandidates("none proton deuteron");

    EBeamCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/ebeam", this);
    EBeamCmd->SetGuidance("Set fE");
    EBeamCmd->SetParameterName("ebeam", false);
    EBeamCmd->SetDefaultUnit("MeV");
    
    ThetaDir = new G4UIdirectory("/pradsim/gun/theta/");
    ThetaDir->SetGuidance("Scattering angle control");
    
    ThetaLowCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/theta/low", this);
    ThetaLowCmd->SetGuidance("Set fThetaLo");
    ThetaLowCmd->SetParameterName("thetalo", false);
    ThetaLowCmd->SetDefaultUnit("deg");
    
    ThetaHighCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/theta/high", this);
    ThetaHighCmd->SetGuidance("Set fThetaHi");
    ThetaHighCmd->SetParameterName("thetahi", false);
    ThetaHighCmd->SetDefaultUnit("deg");

    EvTypeCmd = new G4UIcmdWithAString("/pradsim/gun/evtype", this);
    EvTypeCmd->SetGuidance("Choose a type of event.");
    EvTypeCmd->SetGuidance("  Choice : ep, moller");
    EvTypeCmd->SetParameterName("evtype", false);
    EvTypeCmd->SetCandidates("ep moller");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete ThetaLowCmd;
    delete ThetaHighCmd;
    delete ThetaDir;
    delete GunTypeCmd;
    delete RecoilCmd;
    delete EBeamCmd;
    delete GunDir;
    delete EvTypeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == GunTypeCmd)
        Action->SetGunType(newValue);

    if (command == RecoilCmd)
        Action->SetRecoilParticle(newValue);
    
    if (command == EBeamCmd)
        Action->SetBeamEnergy(EBeamCmd->GetNewDoubleValue(newValue));
    
    if (command == ThetaLowCmd)
        Action->SetThetaRange(ThetaLowCmd->GetNewDoubleValue(newValue), -10000);
    
    if (command == ThetaHighCmd)
        Action->SetThetaRange(-10000, ThetaHighCmd->GetNewDoubleValue(newValue));

    if (command == EvTypeCmd)
        Action->SetEventType(newValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
