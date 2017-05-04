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
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *act) : G4UImessenger(), Action(act)
{
    GunDir = new G4UIdirectory("/pradsim/gun/");
    GunDir->SetGuidance("Primary generator control");

    GunTypeCmd = new G4UIcmdWithAString("/pradsim/gun/type", this);
    GunTypeCmd->SetGuidance("Choose a type of event generator.");
    GunTypeCmd->SetGuidance("  Choice : point, ring, file");
    GunTypeCmd->SetParameterName("guntype", false);
    GunTypeCmd->SetCandidates("point ring file");

    EventTypeCmd = new G4UIcmdWithAString("/pradsim/gun/evtype", this);
    EventTypeCmd->SetGuidance("Choose a type of model.");
    EventTypeCmd->SetGuidance("  Choice : elastic, moller");
    EventTypeCmd->SetParameterName("evtype", false);
    EventTypeCmd->SetCandidates("elastic moller");

    RecoilCmd = new G4UIcmdWithAString("/pradsim/gun/recoil", this);
    RecoilCmd->SetGuidance("Choose a type of recoil particle.");
    RecoilCmd->SetGuidance("  Choice : none, proton, deuteron");
    RecoilCmd->SetParameterName("recoil", false);
    RecoilCmd->SetCandidates("none proton deuteron");

    EBeamCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/ebeam", this);
    EBeamCmd->SetGuidance("Set fE");
    EBeamCmd->SetParameterName("ebeam", false);
    EBeamCmd->SetDefaultUnit("MeV");

    PosCmd = new G4UIcmdWith3VectorAndUnit("/pradsim/gun/pos", this);
    PosCmd->SetGuidance("Set fX,fY,fZ");
    PosCmd->SetParameterName("x", "y", "z", false);
    PosCmd->SetDefaultUnit("mm");

    ThetaCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/theta", this);
    ThetaCmd->SetGuidance("Set fTheta");
    ThetaCmd->SetParameterName("theta", false);
    ThetaCmd->SetDefaultUnit("deg");

    PhiCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/phi", this);
    PhiCmd->SetGuidance("Set fPhi");
    PhiCmd->SetParameterName("phi", false);
    PhiCmd->SetDefaultUnit("deg");

    ThetaLowCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/thetalow", this);
    ThetaLowCmd->SetGuidance("Set fThetaLo");
    ThetaLowCmd->SetParameterName("thetalo", false);
    ThetaLowCmd->SetDefaultUnit("deg");

    ThetaHighCmd = new G4UIcmdWithADoubleAndUnit("/pradsim/gun/thetahigh", this);
    ThetaHighCmd->SetGuidance("Set fThetaHi");
    ThetaHighCmd->SetParameterName("thetahi", false);
    ThetaHighCmd->SetDefaultUnit("deg");

    EventFileCmd = new G4UIcmdWithAString("/pradsim/gun/path", this);
    EventFileCmd->SetGuidance("Choose path of event file");
    EventFileCmd->SetParameterName("path", false);

    TargetProfileCmd = new G4UIcmdWithAString("/pradsim/gun/target", this);
    TargetProfileCmd->SetGuidance("Choose path of target profile");
    TargetProfileCmd->SetParameterName("target", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete ThetaLowCmd;
    delete ThetaHighCmd;
    delete GunTypeCmd;
    delete EventTypeCmd;
    delete RecoilCmd;
    delete EBeamCmd;
    delete PosCmd;
    delete ThetaCmd;
    delete PhiCmd;
    delete EventFileCmd;
    delete TargetProfileCmd;
    delete GunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == GunTypeCmd)
        Action->SetGunType(newValue);

    if (command == EventTypeCmd)
        Action->SetEventType(newValue);

    if (command == RecoilCmd)
        Action->SetRecoilParticle(newValue);

    if (command == EBeamCmd)
        Action->SetBeamEnergy(EBeamCmd->GetNewDoubleValue(newValue));

    if (command == PosCmd)
        Action->SetPosition(PosCmd->GetNew3VectorValue(newValue));

    if (command == ThetaCmd)
        Action->SetTheta(ThetaCmd->GetNewDoubleValue(newValue));

    if (command == PhiCmd)
        Action->SetPhi(PhiCmd->GetNewDoubleValue(newValue));

    if (command == ThetaLowCmd)
        Action->SetThetaRange(ThetaLowCmd->GetNewDoubleValue(newValue), -10000);

    if (command == ThetaHighCmd)
        Action->SetThetaRange(-10000, ThetaHighCmd->GetNewDoubleValue(newValue));

    if (command == EventFileCmd)
        Action->SetEventFile(newValue);

    if (command == TargetProfileCmd)
        Action->SetTargetProfile(newValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
