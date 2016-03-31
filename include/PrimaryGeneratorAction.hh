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
// $Id: PrimaryGeneratorAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "TMath.h"
#include "TLorentzVector.h"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*);    
  virtual ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  
  void Initialize_flags();
  
  void SetEseppFlag(G4String st) {fl_esepp =  (st=="on") ? true : false;}
  void SetMollerFlag(G4String st) {fl_moller = (st=="on") ? true : false;
    fl_onlymoller = (st=="only") ? true : false;}
  void SetRosenFlag(G4String st) {fl_rosen =  (st=="on") ? true : false;}
  void SetQuickFlag(G4String st) {fl_quick =  (st=="on") ? true : false;}
  void SetTargetFlag(G4String st) {fl_target =  (st=="on") ? true : false;}
  void SetLeptonFlag(G4int fl) {fl_lepton = fl;}
  void SetModeFlag(G4int fl) {fl_mode = fl;}
  void SetStructFlag(G4int fl) {fl_struct = fl;}
  void SetTpeFlag(G4int fl) {fl_tpe = fl;}
  void SetVpolFlag(G4int fl) {fl_vpol = fl;}
  void SetScangle(G4double fl) {sc_angle = fl;}
  void SetEli(G4double fl) {E_li_pga = fl;}
  void SetEgcut(G4double fl) {E_g_cut_pga = fl;}
  void SetEgmax(G4double fl) {E_g_max_pga = fl;}
  void SetThetamin(G4double fl) {theta_min_pga = fl*TMath::Pi()/180.;}
  void SetThetamax(G4double fl) {theta_max_pga = fl*TMath::Pi()/180.;}
  void SetPhimin(G4double fl) {phi_min_pga = fl*TMath::Pi()/180.;}
  void SetPhimax(G4double fl) {phi_max_pga = fl*TMath::Pi()/180.;}

  TLorentzVector GetV_l1i() {if (event_type_pga->contains("m")) return *l1i_pga;
    else return *v_li_pga;}
  TLorentzVector GetV_l2i() {if (event_type_pga->contains("m")) return *l2i_pga;
    else return *v_pi_pga;}
  TLorentzVector GetV_l1f() {if (event_type_pga->contains("m")) return *l1f_pga;
    else return *v_lf_pga;}
  TLorentzVector GetV_l2f() {if (event_type_pga->contains("m")) return *l2f_pga;
    else return *v_pf_pga;}
  TLorentzVector GetV_kf() {if (event_type_pga->contains("m")) return *vkf_pga;
    else return *v_kf_pga;}
  TLorentzVector GetCoord() {return *coord_pga;}
  G4double GetMli() {return *m_l_pga;}
  G4String GetEventType() {return *event_type_pga;}
   
private:
  G4ParticleGun* particleGun;	 //pointer a to G4  class
  DetectorConstruction* detector;     //pointer to the geometry
    
  PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class

  G4bool fl_esepp, fl_quick, fl_target, fl_rosen, fl_moller, fl_onlymoller;
  G4int fl_lepton, fl_mode, fl_struct, fl_tpe, fl_vpol;
  
  G4double sc_angle;
  G4double E_li_pga;
  G4double E_g_cut_pga;
  G4double E_g_max_pga;
  G4double theta_min_pga;
  G4double theta_max_pga;
  G4double phi_min_pga;
  G4double phi_max_pga;

  G4double* m_l_pga;
  G4String* event_type_pga;
  TLorentzVector *coord_pga;
  TLorentzVector *v_li_pga, *v_pi_pga, *v_lf_pga, *v_pf_pga, *v_kf_pga;
  TLorentzVector *l1i_pga, *l2i_pga, *l1f_pga, *l2f_pga, *vkf_pga;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


