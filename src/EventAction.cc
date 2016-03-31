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
// $Id: EventAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "EventActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <cmath>
#include <vector>
#include "g4root.hh"
#include "TLorentzVector.h"

#include "Randomize.hh"
#include <iomanip>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PrimaryGeneratorAction* pgAction)
  : pga(pgAction),eventMessenger(0),printModulo(100000)
{
  eventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    CLHEP::HepRandom::showEngineStatus();
  }
  
  charge_o.clear();
  mass_o.clear();
  px_o.clear();
  py_o.clear();
  pz_o.clear();
  E_o.clear();
  theta_o.clear();
  phi_o.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Get vertex info
  G4int eventID = evt->GetEventID();
  G4int nOut = evt->GetNumberOfPrimaryVertex();

  G4String event_type = pga->GetEventType();
  G4double m_l = pga->GetMli();
  TLorentzVector coord = pga->GetCoord();
  TLorentzVector v_l1i = pga->GetV_l1i();
  TLorentzVector v_l2i = pga->GetV_l2i();
  TLorentzVector v_l1f = pga->GetV_l1f();
  TLorentzVector v_l2f = pga->GetV_l2f();
  TLorentzVector v_kf = pga->GetV_kf();
  vector<TLorentzVector> v_out = {v_l1f,v_l2f,v_kf};

  analysisManager->FillNtupleIColumn(0,eventID);
  analysisManager->FillNtupleSColumn(1,event_type);
  analysisManager->FillNtupleDColumn(2,m_l);
  analysisManager->FillNtupleDColumn(3,v_l1i.E());
  analysisManager->FillNtupleDColumn(4,v_l1i.Px());
  analysisManager->FillNtupleDColumn(5,v_l1i.Py());
  analysisManager->FillNtupleDColumn(6,v_l1i.Pz());
  analysisManager->FillNtupleDColumn(7,v_l1i.Theta());
  analysisManager->FillNtupleDColumn(8,v_l1i.Phi());
  analysisManager->FillNtupleDColumn(9,coord.T());
  analysisManager->FillNtupleDColumn(10,coord.X());
  analysisManager->FillNtupleDColumn(11,coord.Y());
  analysisManager->FillNtupleDColumn(12,coord.Y());
  analysisManager->FillNtupleIColumn(13,nOut);
  analysisManager->FillNtupleDColumn(14,-(v_l1i-v_l1f).Mag2());
  analysisManager->FillNtupleDColumn(15,(v_l1i-v_l1f).E());
  analysisManager->FillNtupleDColumn(16,(v_l1i-v_l1f).E()/v_l1i.E());
 
  for (G4int i=0; i<nOut; i++) {
    G4PrimaryVertex* pv = evt->GetPrimaryVertex(i);
    G4PrimaryParticle* par_o = pv->GetPrimary();
    charge_o.push_back(par_o->GetCharge());
    mass_o.push_back(par_o->GetMass());
    E_o.push_back(v_out[i].E());
    px_o.push_back(v_out[i].Px());
    py_o.push_back(v_out[i].Py());
    pz_o.push_back(v_out[i].Pz());
    theta_o.push_back(v_out[i].Theta());
    phi_o.push_back(v_out[i].Phi());
  }
  
  analysisManager->AddNtupleRow(0);
  
}
