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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

#include <cmath>
#include "GenEvent.hh"
#include "TLorentzVector.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
  :detector(DC),
   fl_lepton(1), fl_mode(12), fl_struct(2), fl_tpe(1), fl_vpol(3), fl_quick(false), fl_target(false), fl_rosen(false), fl_esepp(false), fl_moller(false), fl_onlymoller(true),
   sc_angle(0.8*degrad), E_li_pga(1.), E_g_cut_pga(0.0001), E_g_max_pga(0.7), theta_min_pga(0.2*degrad), theta_max_pga(10.*degrad), phi_min_pga(-180.*degrad), phi_max_pga(180.*degrad) 
   
{
  particleGun  = new G4ParticleGun(1);
  gunMessenger = new PrimaryGeneratorMessenger(this);
  
  m_l_pga = &m_l;
  event_type_pga = &event_type;
  coord_pga = &coord;
  v_li_pga = &v_li;
  v_pi_pga = &v_pi;
  l1i_pga = &l1i;
  l2i_pga = &l2i;
  v_lf_pga = &v_lf;
  v_pf_pga = &v_pf;
  v_kf_pga = &v_kf;
  l1f_pga = &l1f;
  l2f_pga = &l2f;
  vkf_pga = &vkf;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::Initialize_flags()
{
  flag_esepp = fl_esepp;
  flag_lepton = fl_lepton;
  flag_mode = fl_mode;
  flag_struct = fl_struct;
  flag_tpe = fl_tpe;
  flag_vpol = fl_vpol;
  flag_quick = fl_quick;
  flag_target = fl_target;
  flag_rosen = fl_rosen;
  flag_moller = fl_moller;
  flag_onlymoller = fl_onlymoller;
  
  E_li = E_li_pga;
  E_g_cut = E_g_cut_pga;
  E_g_max = E_g_max_pga;
  theta_min = theta_min_pga;
  theta_max = theta_max_pga;
  phi_min = phi_min_pga;
  phi_max = phi_max_pga;

  coord = TLorentzVector(0.,0.,-300.*cm,0.);
  cell_xmin = coord.X()-2.*mm;
  cell_xmax = coord.X()+2.*mm;
  cell_ymin = coord.Y()-2.*mm;
  cell_ymax = coord.Y()+2.*mm;
  cell_zmin = coord.Z()-2.*cm;
  cell_zmax = coord.Z()+2.*cm;
  
  if (flag_lepton<4) m_l = m_e;
  else m_l = m_mu;
  
  m_l2 = Pow2(m_l); m_l4 = Pow4(m_l);
  
  if (flag_esepp) Initialize_genevent();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // initialize flags and kinematics
  if (anEvent->GetEventID()==0) {
    nevents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    Initialize_flags();
  }

  // generation of vertex position
  if (flag_target) coord = scell();
  particleGun->SetParticlePosition(G4ThreeVector(coord.X(), coord.Y(), coord.Z()));

  // generation of e- without radiative correction at a given angle
  if (!flag_esepp) {
    particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
    G4double phi = 2.*Pi*G4UniformRand();
    G4double kx_lf = sin(sc_angle)*cos(phi);
    G4double ky_lf = sin(sc_angle)*sin(phi);
    G4double kz_lf = cos(sc_angle);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(kx_lf, ky_lf, kz_lf));
    particleGun->SetParticleEnergy(E_li*GeV);
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  // generation of events by the RC generator
  else {
    
    Loop_genevent();
    // moller events
    if (event_type=="elm" || event_type=="brm") {
      
      particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(l1f.Px(),l1f.Py(),l1f.Pz()));
      particleGun->SetParticleEnergy(l1f.E()*GeV);
      particleGun->GeneratePrimaryVertex(anEvent);
      
      particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(l2f.Px(),l2f.Py(),l2f.Pz()));
      particleGun->SetParticleEnergy(l2f.E()*GeV);
      particleGun->GeneratePrimaryVertex(anEvent);

      if (event_type=="brm") {
	particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(vkf.Px(),vkf.Py(),vkf.Pz()));
	particleGun->SetParticleEnergy(vkf.E()*GeV);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
      
    }
    // proton scattering events
    else {
      particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(v_lf.Px(),v_lf.Py(),v_lf.Pz()));
      particleGun->SetParticleEnergy(v_lf.E()*GeV);
      particleGun->GeneratePrimaryVertex(anEvent);
      
      particleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(v_pf.Px(),v_pf.Py(),v_pf.Pz()));
      particleGun->SetParticleEnergy(v_pf.E()*GeV);
      particleGun->GeneratePrimaryVertex(anEvent);

      if (event_type=="bre1" || event_type=="bre2") {
	particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(v_kf.Px(),v_kf.Py(),v_kf.Pz()));
	particleGun->SetParticleEnergy(v_kf.E()*GeV);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

