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
// $Id: RunAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "EventAction.hh"
#include "g4root.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* evA)
  : evAction(evA)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("test_prad");

  //event info
  analysisManager->CreateNtuple("events", "events");
  analysisManager->CreateNtupleIColumn("eventID"); //0
  analysisManager->CreateNtupleSColumn("eventType"); //1
  //incoming particle info
  analysisManager->CreateNtupleDColumn("mass_i"); //2
  analysisManager->CreateNtupleDColumn("E_i"); //3
  analysisManager->CreateNtupleDColumn("px_i"); //4
  analysisManager->CreateNtupleDColumn("py_i"); //5
  analysisManager->CreateNtupleDColumn("pz_i"); //6
  analysisManager->CreateNtupleDColumn("theta_i"); //7
  analysisManager->CreateNtupleDColumn("phi_i"); //8
  //vertex info
  analysisManager->CreateNtupleDColumn("t_vertex"); //9
  analysisManager->CreateNtupleDColumn("x_vertex"); //10 
  analysisManager->CreateNtupleDColumn("y_vertex"); //11
  analysisManager->CreateNtupleDColumn("z_vertex"); //12
  analysisManager->CreateNtupleIColumn("nOut"); //13
  analysisManager->CreateNtupleDColumn("Q2"); //14
  analysisManager->CreateNtupleDColumn("nu"); //15
  analysisManager->CreateNtupleDColumn("y"); //16
  //outgoing particles info
  analysisManager->CreateNtupleIColumn("charge_o", evAction->GetCharge_o()); //17
  analysisManager->CreateNtupleDColumn("mass_o", evAction->GetMass_o()); //18
  analysisManager->CreateNtupleDColumn("E_o", evAction->GetE_o()); //19
  analysisManager->CreateNtupleDColumn("px_o", evAction->GetPx_o()); //20
  analysisManager->CreateNtupleDColumn("py_o", evAction->GetPy_o()); //21
  analysisManager->CreateNtupleDColumn("pz_o", evAction->GetPz_o()); //22
  analysisManager->CreateNtupleDColumn("theta_o", evAction->GetTheta_o()); //23
  analysisManager->CreateNtupleDColumn("phi_o", evAction->GetPhi_o()); //24
   
  analysisManager->FinishNtuple();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
