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
// $Id: PrimaryGeneratorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                          PrimaryGeneratorAction* pga)
:pgaction(pga)
{
  genDir = new G4UIdirectory("/genevent/");
  genDir->SetGuidance("PrimaryGenerator control");
   
  eseppCmd = new G4UIcmdWithAString("/genevent/esepp",this);
  eseppCmd->SetGuidance("Generate event with radiative correction");
  eseppCmd->SetGuidance("  Choice : on, off(default)");
  eseppCmd->SetParameterName("choice",true);
  eseppCmd->SetDefaultValue("off");
  eseppCmd->SetCandidates("on off");

  mollerCmd = new G4UIcmdWithAString("/genevent/moller",this);
  mollerCmd->SetGuidance("Add moller events to generator");
  mollerCmd->SetGuidance("  Choice : on, off(default), only");
  mollerCmd->SetParameterName("choice",true);
  mollerCmd->SetDefaultValue("off");
  mollerCmd->SetCandidates("on off only");

  rosenCmd = new G4UIcmdWithAString("/genevent/rosen",this);
  rosenCmd->SetGuidance("Add rosen events to generator");
  rosenCmd->SetGuidance("  Choice : on, off(default)");
  rosenCmd->SetParameterName("choice",true);
  rosenCmd->SetDefaultValue("off");
  rosenCmd->SetCandidates("on off");

  quickCmd = new G4UIcmdWithAString("/genevent/quick",this);
  quickCmd->SetGuidance("Quick initialization");
  quickCmd->SetGuidance("  Choice : on, off(default)");
  quickCmd->SetParameterName("choice",true);
  quickCmd->SetDefaultValue("off");
  quickCmd->SetCandidates("on off");

  targetCmd = new G4UIcmdWithAString("/genevent/target",this);
  targetCmd->SetGuidance("Target");
  targetCmd->SetGuidance("  Choice : on->full target, off->punctual(default)");
  targetCmd->SetParameterName("choice",true);
  targetCmd->SetDefaultValue("off");
  targetCmd->SetCandidates("on off");

  leptonCmd = new G4UIcmdWithAnInteger("/genevent/lepton",this);
  leptonCmd->SetGuidance("Choose lepton");
  leptonCmd->SetGuidance("  Choice : e- (1)(default), e+ (2), e- and e+ (3), mu- (4), mu+ (5), mu- and mu+ (6)");
  leptonCmd->SetParameterName("choice",true);
  leptonCmd->SetDefaultValue(1);

  modeCmd = new G4UIcmdWithAnInteger("/genevent/mode",this);
  modeCmd->SetGuidance("Choose mode");
  modeCmd->SetGuidance("  Choice : cf esepp code");
  modeCmd->SetParameterName("choice",true);
  modeCmd->SetDefaultValue(12);

  structCmd = new G4UIcmdWithAnInteger("/genevent/struct",this);
  structCmd->SetGuidance("Choose struct");
  structCmd->SetGuidance("  Choice : cf esepp code");
  structCmd->SetParameterName("choice",true);
  structCmd->SetDefaultValue(5);

  tpeCmd = new G4UIcmdWithAnInteger("/genevent/tpe",this);
  tpeCmd->SetGuidance("Choose tpe");
  tpeCmd->SetGuidance("  Choice : cf esepp code");
  tpeCmd->SetParameterName("choice",true);
  tpeCmd->SetDefaultValue(3);

  vpolCmd = new G4UIcmdWithAnInteger("/genevent/vpol",this);
  vpolCmd->SetGuidance("Choose vpol");
  vpolCmd->SetGuidance("  Choice : cf esepp code");
  vpolCmd->SetParameterName("choice",true);
  vpolCmd->SetDefaultValue(3);

  scangleCmd = new G4UIcmdWithADouble("/genevent/scangle",this);
  scangleCmd->SetGuidance("Choose scattered angle (no RC) in degree");
  scangleCmd->SetParameterName("sc_angle",true);
  scangleCmd->SetDefaultValue(0.8);

  eliCmd = new G4UIcmdWithADouble("/genevent/eli",this);
  eliCmd->SetGuidance("Choose incident lepton energy in GeV");
  eliCmd->SetParameterName("E_li",true);
  eliCmd->SetDefaultValue(1.1);

  egcutCmd = new G4UIcmdWithADouble("/genevent/egcut",this);
  egcutCmd->SetGuidance("Choose minimal photon energy in GeV");
  egcutCmd->SetParameterName("E_g_cut",true);
  egcutCmd->SetDefaultValue(0.0001);

  egmaxCmd = new G4UIcmdWithADouble("/genevent/egmax",this);
  egmaxCmd->SetGuidance("Choose maximal photon energy in GeV");
  egmaxCmd->SetParameterName("E_g_max",true);
  egmaxCmd->SetDefaultValue(0.7);

  thetaminCmd = new G4UIcmdWithADouble("/genevent/thetamin",this);
  thetaminCmd->SetGuidance("Choose minimal angle in degree");
  thetaminCmd->SetParameterName("theta_min",true);
  thetaminCmd->SetDefaultValue(0.1);

  thetamaxCmd = new G4UIcmdWithADouble("/genevent/thetamax",this);
  thetamaxCmd->SetGuidance("Choose maximal angle in degree");
  thetamaxCmd->SetParameterName("theta_max",true);
  thetamaxCmd->SetDefaultValue(10.);

  phiminCmd = new G4UIcmdWithADouble("/genevent/phimin",this);
  phiminCmd->SetGuidance("Choose minimal azimuthal angle in degree");
  phiminCmd->SetParameterName("phi_min",true);
  phiminCmd->SetDefaultValue(-180.);

  phimaxCmd = new G4UIcmdWithADouble("/genevent/phimax",this);
  phimaxCmd->SetGuidance("Choose maximal azimuthal angle in degree");
  phimaxCmd->SetParameterName("phi_max",true);
  phimaxCmd->SetDefaultValue(180.);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete eseppCmd;
  delete mollerCmd;
  delete rosenCmd;
  delete quickCmd;
  delete targetCmd;
  delete leptonCmd;
  delete modeCmd;
  delete structCmd;
  delete tpeCmd;
  delete vpolCmd;
  delete scangleCmd;
  delete eliCmd;
  delete egcutCmd;
  delete egmaxCmd;
  delete thetaminCmd;
  delete thetamaxCmd;
  delete phiminCmd;
  delete phimaxCmd;
  delete genDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == eseppCmd) pgaction->SetEseppFlag(newValue);
  if( command == mollerCmd) pgaction->SetMollerFlag(newValue);
  if( command == rosenCmd) pgaction->SetRosenFlag(newValue);
  if( command == quickCmd) pgaction->SetQuickFlag(newValue);
  if( command == targetCmd) pgaction->SetTargetFlag(newValue);
  if( command == leptonCmd) pgaction->SetLeptonFlag(leptonCmd->GetNewIntValue(newValue));
  if( command == modeCmd) pgaction->SetModeFlag(modeCmd->GetNewIntValue(newValue));
  if( command == structCmd) pgaction->SetStructFlag(structCmd->GetNewIntValue(newValue));
  if( command == tpeCmd) pgaction->SetTpeFlag(tpeCmd->GetNewIntValue(newValue));
  if( command == vpolCmd) pgaction->SetVpolFlag(vpolCmd->GetNewIntValue(newValue));
  if( command == scangleCmd) pgaction->SetScangle(scangleCmd->GetNewDoubleValue(newValue));
  if( command == egcutCmd) pgaction->SetEgcut(egcutCmd->GetNewDoubleValue(newValue));
  if( command == egmaxCmd) pgaction->SetEgmax(egmaxCmd->GetNewDoubleValue(newValue));
  if( command == eliCmd) pgaction->SetEli(eliCmd->GetNewDoubleValue(newValue));
  if( command == thetaminCmd) pgaction->SetThetamin(thetaminCmd->GetNewDoubleValue(newValue));
  if( command == thetamaxCmd) pgaction->SetThetamax(thetamaxCmd->GetNewDoubleValue(newValue));
  if( command == phiminCmd) pgaction->SetPhimin(phiminCmd->GetNewDoubleValue(newValue));
  if( command == phimaxCmd) pgaction->SetPhimax(phimaxCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

