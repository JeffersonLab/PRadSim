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
// $Id: EventAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

class EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(PrimaryGeneratorAction*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);

  std::vector<G4String>& GetType_o(){return type_o;}
  std::vector<G4int>& GetCharge_o(){return charge_o;}
  std::vector<G4double>& GetMass_o(){return mass_o;}
  std::vector<G4double>& GetE_o(){return E_o;}
  std::vector<G4double>& GetPx_o(){return px_o;}
  std::vector<G4double>& GetPy_o(){return py_o;}
  std::vector<G4double>& GetPz_o(){return pz_o;}
  std::vector<G4double>& GetTheta_o(){return theta_o;}
  std::vector<G4double>& GetPhi_o(){return phi_o;}

  void SetPrintModulo(G4int val) {printModulo = val;};
  
private:
  
  EventActionMessenger*  eventMessenger;
  PrimaryGeneratorAction* pga;
  G4int printModulo;
  std::vector<G4String> type_o;
  std::vector<G4int> charge_o;
  std::vector<G4double> mass_o;
  std::vector<G4double> E_o;
  std::vector<G4double> px_o;
  std::vector<G4double> py_o;
  std::vector<G4double> pz_o;
  std::vector<G4double> theta_o;
  std::vector<G4double> phi_o;
  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
