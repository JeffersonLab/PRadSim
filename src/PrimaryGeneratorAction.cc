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
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh" 
#include "G4SystemOfUnits.hh"

#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(
                                             DetectorConstruction* DC)
:Detector(DC),rndmFlag("on")
{
  n_particle = 1;
  particles.open("ep_out_e-.dat");
  particleGun  = new G4ParticleGun(n_particle);
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(1100.*MeV);
//  G4double position = -0.5*(Detector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm,-280.*cm));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  particles.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 

  if (rndmFlag = "on")
  {
    G4double x0 = 0., y0 = 0., z0 = -330.;
    G4double Ene = 0.;
    G4double kx = 0., ky = 0., kz = 1.;
    G4double theta, phi;
//    G4double tmp1[3], tmp2[3], tmp3[3];
//    G4ParticleDefinition* particle;
//    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;

/*
    G4double halosig = 0.12;
    G4double r0 = G4RandGauss::shoot(0, halosig);
    theta = 3.141592653589793*2.*G4UniformRand();
    x0 = sin(theta)*r0;
    y0 = cos(theta)*r0;

    particleGun->SetParticlePosition(G4ThreeVector(x0*cm,y0*cm,z0*cm));
    particleGun->SetParticleEnergy(1100.*MeV);
    particleGun->GeneratePrimaryVertex(anEvent);
*/

//    G4double tmp1, tmp2, tmp3, PID;
/*
    G4ParticleDefinition* particle
                = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    particleGun->SetParticleDefinition(particle);
*/

    // generate a theta ring
    theta = 0.8;
    Ene = 1097*MeV;
    x0 = 0.*cm;
    y0 = 0.*cm;
    z0 = -300.*cm + 4.*(0.5 - G4UniformRand())*cm;
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    theta = theta/180.*3.14159265358979;
    phi = 2.*3.14159265358979*G4UniformRand();
    kx = sin(theta)*cos(phi);
    ky = sin(theta)*sin(phi);
    kz = cos(theta);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
    particleGun->SetParticleEnergy(Ene);
    particleGun->GeneratePrimaryVertex(anEvent);



  //RCEP PART
/*
    particles >> tmp1[0] >> tmp2[0] >> tmp3[0] >> tmp1[1] >> tmp2[1] >> tmp3[1] >> tmp1[2] >> tmp2[2] >> tmp3[2];

    particle = particleTable->FindParticle(particleName="e-");
    particleGun->SetParticleDefinition(particle);
    x0 = G4RandGauss::shoot(0., 0.008)*cm;
    y0 = G4RandGauss::shoot(0., 0.008)*cm;
    z0 = -300.*cm + 4.*(0.5 - G4UniformRand())*cm;
    Ene = tmp1[0];
    theta = tmp2[0];
    phi = tmp3[0];
//    G4cout << Ene << "  " << theta/3.14159265358979*180. << G4endl;
    kx = sin(theta)*cos(phi);
    ky = sin(theta)*sin(phi);
    kz = cos(theta);
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
    particleGun->SetParticleEnergy(Ene*MeV);
    particleGun->GeneratePrimaryVertex(anEvent);

    if(tmp1[2] > 0.00) {
      particle = particleTable->FindParticle(particleName="gamma");
      particleGun->SetParticleDefinition(particle);
      Ene = tmp1[2];
      theta = tmp2[2];
      phi = tmp3[2];
      kx = sin(theta)*cos(phi);
      ky = sin(theta)*sin(phi);
      kz = cos(theta);
      particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
      particleGun->SetParticleEnergy(Ene*MeV);
      particleGun->GeneratePrimaryVertex(anEvent);
    }
*/
//    particleGun->GeneratePrimaryVertex(anEvent);


/*
  //RCEE PART
    particles >> tmp1[0] >> tmp2[0] >> tmp3[0] >> tmp1[1] >> tmp2[1] >> tmp3[1] >> tmp1[2] >> tmp2[2] >> tmp3[2];
  //1st e
    particle = particleTable->FindParticle(particleName="e-");
    particleGun->SetParticleDefinition(particle);
    x0 = 0.*cm;
    y0 = 0.*cm;
    z0 = -300.*cm + 4.*(0.5 - G4UniformRand())*cm;
    Ene = tmp1[0];
    theta = tmp2[0]; 
    phi = tmp3[0];
//    G4cout << Ene << "  " << theta/3.14159265358979*180. << G4endl;
    kx = sin(theta)*cos(phi);
    ky = sin(theta)*sin(phi);
    kz = cos(theta);
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
    particleGun->SetParticleEnergy(Ene*MeV);
    particleGun->GeneratePrimaryVertex(anEvent);

  //2nd e
    Ene = tmp1[1];
    theta = tmp2[1];
    phi = tmp3[1];
    kx = sin(theta)*cos(phi);
    ky = sin(theta)*sin(phi);
    kz = cos(theta);
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
    particleGun->SetParticleEnergy(Ene*MeV);
    particleGun->GeneratePrimaryVertex(anEvent);

    if(tmp1[2] > 0.00) {
      particle = particleTable->FindParticle(particleName="gamma");
      particleGun->SetParticleDefinition(particle);
      Ene = tmp1[2];
      theta = tmp2[2];
      phi = tmp3[2];
      kx = sin(theta)*cos(phi);
      ky = sin(theta)*sin(phi);
      kz = cos(theta);
      particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(kx, ky, kz));
      particleGun->SetParticleEnergy(Ene*MeV);
      particleGun->GeneratePrimaryVertex(anEvent);
    }
*/


/*
    G4cout << x0 <<"  " << y0 << "  " << z0 << G4endl;
    G4cout << kx << "  " << ky << "  " << kz << G4endl;
    G4cout << Ene << G4endl;
    G4cout << particleGun->GetParticleDefinition()->GetParticleName() << G4endl;  
*/
//    n_particle += 1;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

