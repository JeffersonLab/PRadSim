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
// PrimaryGenerator.cc
// Developer : Chao Gu
// History:
//   Mar 2017, C. Gu, Add for DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "ConfigParser.h"
#include "Globals.hh"
#include "RootTree.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPrimaryGenerator.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <cmath>
#include <fstream>

static double me = 0.510998928 * MeV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator(G4String type, G4double e, G4double thlo, G4double thhi, G4bool rec, G4String par) : G4VPrimaryGenerator(), fRegistered(false), fEventType(type), fRecoilOn(rec), fRecoilParticle(par), fTargetInfo(false), fTargetCenter(-300 * cm), fTargetHalfL(0), fEBeam(e), fThetaLo(thlo), fThetaHi(thhi), fTargetMass(0)
{
    if (fRecoilParticle != "proton" && fRecoilParticle != "deuteron")
        fRecoilParticle = "proton";

    if (fEventType != "elastic" && fEventType != "moller")
        fEventType = "elastic";

    fN = 0;

    for (int i = 0; i < MaxN; i++) {
        fPID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fE[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event *anEvent)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    Clear();

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    if (!fTargetInfo) {
        G4VPhysicalVolume *physiTargetCon = G4PhysicalVolumeStore::GetInstance()->GetVolume("Target Container");
        G4LogicalVolume *logicTarget = G4LogicalVolumeStore::GetInstance()->GetVolume("TargetLV");

        if (physiTargetCon) fTargetCenter = physiTargetCon->GetObjectTranslation().z();

        G4Tubs *solidTarget = NULL;

        if (logicTarget) solidTarget = dynamic_cast<G4Tubs *>(logicTarget->GetSolid());

        if (solidTarget)
            fTargetHalfL = solidTarget->GetZHalfLength();
        else
            G4cout << "WARNING: target volume not found" << G4endl;

        fTargetMass = particleTable->FindParticle(fRecoilParticle)->GetPDGMass();
    }

    double x = G4RandGauss::shoot(0, 0.08) * mm;
    double y = G4RandGauss::shoot(0, 0.08) * mm;
    double z = fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());
    double theta_l = fThetaLo + (fThetaHi - fThetaLo) * G4UniformRand();
    double phi_l = twopi * G4UniformRand();

    double M = fTargetMass;
    double E = fEBeam;
    double P = sqrt(E * E - me * me);
    double cosang = cos(theta_l);
    double p_l = 0, e_l = 0, a = 0;

    if (fEventType == "elastic") { // me is not ignored
        p_l = (P * M / (E + M - P * cosang)) * (((E + M) * sqrt(1 - (me / M) * (me / M) * (1 - cosang * cosang)) + (E + (me / M) * me) * cosang) / (E + M + P * cosang));
        e_l = sqrt(p_l * p_l + me * me);
    } else if (fEventType == "moller") {
        a = (E - me) / (E + me);
        e_l = me * (1 + a * cosang * cosang) / (1 - a * cosang * cosang);
        p_l = sqrt(e_l * e_l - me * me);
    }

    if (e_l > E) G4cout << "WARNING: super-elastic event found" << G4endl;

    G4PrimaryVertex *vertexL = new G4PrimaryVertex(x, y, z, 0);
    G4PrimaryParticle *particleL = new G4PrimaryParticle(particleTable->FindParticle("e-"));
    double kx_l = sin(theta_l) * cos(phi_l);
    double ky_l = sin(theta_l) * sin(phi_l);
    double kz_l = cos(theta_l);
    particleL->SetMomentumDirection(G4ThreeVector(kx_l, ky_l, kz_l));
    particleL->SetTotalEnergy(e_l);
    vertexL->SetPrimary(particleL);

    anEvent->AddPrimaryVertex(vertexL);

    fPID[fN] = particleL->GetPDGcode();
    fX[fN] = x;
    fY[fN] = y;
    fZ[fN] = z;
    fE[fN] = e_l;
    fMomentum[fN] = p_l;
    fTheta[fN] = theta_l;
    fPhi[fN] = phi_l;
    fN++;

    if (fEventType == "moller") {
        double e_l2 = E - e_l;
        double p_l2 = sqrt(e_l2 * e_l2  - me * me);
        double theta_l2 = asin(p_l * sin(theta_l) / p_l2);
        double phi_l2 = phi_l + pi;

        G4PrimaryVertex *vertexL2 = new G4PrimaryVertex(x, y, z, 0);
        G4PrimaryParticle *particleL2 = new G4PrimaryParticle(particleTable->FindParticle("e-"));
        double kx_l2 = sin(theta_l2) * cos(phi_l2);
        double ky_l2 = sin(theta_l2) * sin(phi_l2);
        double kz_l2 = cos(theta_l2);
        particleL2->SetMomentumDirection(G4ThreeVector(kx_l2, ky_l2, kz_l2));
        particleL2->SetTotalEnergy(e_l2);
        vertexL2->SetPrimary(particleL2);

        anEvent->AddPrimaryVertex(vertexL2);

        fPID[fN] = particleL2->GetPDGcode();
        fX[fN] = x;
        fY[fN] = y;
        fZ[fN] = z;
        fE[fN] = e_l2;
        fMomentum[fN] = p_l2;
        fTheta[fN] = theta_l2;
        fPhi[fN] = phi_l2;
        fN++;
    } else if (fRecoilOn) {
        double theta_h = atan(p_l * sin(theta_l) / (P - p_l * cosang));

        if (theta_h < 0) theta_h = theta_h + twopi;

        if (theta_h > twopi) theta_h = theta_h - twopi;

        double phi_h = phi_l + pi;

        if (phi_h > twopi) phi_h = phi_h - twopi;

        double p_h = sqrt(P * P - 2 * P * p_l * cosang + p_l * p_l);
        double e_h = sqrt(p_h * p_h + M * M);

        G4PrimaryVertex *vertexH = new G4PrimaryVertex(x, y, z, 0);
        G4PrimaryParticle *particleH = new G4PrimaryParticle(particleTable->FindParticle(fRecoilParticle));
        double kx_h = sin(theta_h) * cos(phi_h);
        double ky_h = sin(theta_h) * sin(phi_h);
        double kz_h = cos(theta_h);
        particleH->SetMomentumDirection(G4ThreeVector(kx_h, ky_h, kz_h));
        particleH->SetTotalEnergy(e_h);
        vertexH->SetPrimary(particleH);

        anEvent->AddPrimaryVertex(vertexH);

        fPID[fN] = particleH->GetPDGcode();
        fX[fN] = x;
        fY[fN] = y;
        fZ[fN] = z;
        fE[fN] = e_h;
        fMomentum[fN] = p_h;
        fTheta[fN] = theta_h;
        fPhi[fN] = phi_h;
        fN++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::Register(TTree *tree)
{
    tree->Branch("GUN.N", &fN, "GUN.N/I");
    tree->Branch("GUN.PID", fPID, "GUN.PID[GUN.N]/I");
    tree->Branch("GUN.X", fX, "GUN.X[GUN.N]/D");
    tree->Branch("GUN.Y", fY, "GUN.Y[GUN.N]/D");
    tree->Branch("GUN.Z", fZ, "GUN.Z[GUN.N]/D");
    tree->Branch("GUN.E", fE, "GUN.E[GUN.N]/D");
    tree->Branch("GUN.P", fMomentum, "GUN.P[GUN.N]/D");
    tree->Branch("GUN.Theta", fTheta, "GUN.Theta[GUN.N]/D");
    tree->Branch("GUN.Phi", fPhi, "GUN.Phi[GUN.N]/D");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::Print() const
{
    G4int prec = G4cout.precision(3);

    for (int i = 0; i < fN; i++) {
        G4cout << std::setw(10) << fPID[i] << " ";
        G4cout << std::setw(5) << G4BestUnit(fX[i], "Length") << " " << std::setw(5) << G4BestUnit(fY[i], "Length") << " " << std::setw(5) << G4BestUnit(fZ[i], "Length") << " ";
        G4cout << std::setw(5) << G4BestUnit(fE[i], "Energy") << " " << std::setw(8) << G4BestUnit(fMomentum[i], "Energy") << " ";
        G4cout << std::setw(5) << fTheta[i] / pi * 180 << " deg ";
        G4cout << G4endl;
    }

    G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::Clear()
{
    for (int i = 0; i < fN; i++) {
        fPID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fE[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
    }

    fN = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::PRadPrimaryGenerator(G4String type, G4bool rec, G4String par) : PrimaryGenerator(type, 0, 0, 0, rec, par)
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::PRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path) : PrimaryGenerator(type, 0, 0, 0, rec, par)
{
    if (path.empty()) {
        if (fEventType == "elastic")
            path = "epelastic.dat";
        else if (fEventType == "moller")
            path = "moller.dat";
    }

    if (!fParser.ReadFile(path)) {
        G4cout << "ERROR: failed to read event file " << "\"" << path << "\"" << G4endl;
        exit(1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::~PRadPrimaryGenerator()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadPrimaryGenerator::GeneratePrimaryVertex(G4Event *anEvent)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    Clear();

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    if (!fTargetInfo) {
        G4VPhysicalVolume *physiTargetCon = G4PhysicalVolumeStore::GetInstance()->GetVolume("Target Container");
        G4LogicalVolume *logicTarget = G4LogicalVolumeStore::GetInstance()->GetVolume("TargetLV");

        if (physiTargetCon) fTargetCenter = physiTargetCon->GetObjectTranslation().z();

        G4Tubs *solidTarget = NULL;

        if (logicTarget)
            solidTarget = dynamic_cast<G4Tubs *>(logicTarget->GetSolid());

        if (solidTarget)
            fTargetHalfL = solidTarget->GetZHalfLength();
        else
            G4cout << "WARNING: target volume not found" << G4endl;
    }

    double e_l = 0, theta_l = 0, phi_l = 0;
    double e_h = 0, theta_h = 0, phi_h = 0;
    double e_p = 0, theta_p = 0, phi_p = 0;

    while (fParser.ParseLine()) {
        if (!fParser.CheckElements(9))
            continue;
        else {
            fParser >> e_l >> theta_l >> phi_l >> e_h >> theta_h >> phi_h >> e_p >> theta_p >> phi_p;
            break;
        }
    }

    double x = G4RandGauss::shoot(0, 0.08) * mm;
    double y = G4RandGauss::shoot(0, 0.08) * mm;
    double z = fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());

    G4PrimaryVertex *vertexL = new G4PrimaryVertex(x, y, z, 0);
    G4PrimaryParticle *particleL = new G4PrimaryParticle(particleTable->FindParticle("e-"));
    double m_l = particleL->GetParticleDefinition()->GetPDGMass();
    double p_l = sqrt(e_l * e_l - m_l * m_l);
    double kx_l = sin(theta_l) * cos(phi_l);
    double ky_l = sin(theta_l) * sin(phi_l);
    double kz_l = cos(theta_l);
    particleL->SetMomentumDirection(G4ThreeVector(kx_l, ky_l, kz_l));
    particleL->SetTotalEnergy(e_l);
    vertexL->SetPrimary(particleL);

    anEvent->AddPrimaryVertex(vertexL);

    fPID[fN] = particleL->GetPDGcode();
    fX[fN] = x;
    fY[fN] = y;
    fZ[fN] = z;
    fE[fN] = e_l;
    fMomentum[fN] = p_l;
    fTheta[fN] = theta_l;
    fPhi[fN] = phi_l;
    fN++;

    if (fRecoilOn || fEventType == "moller") {
        G4PrimaryVertex *vertexH = new G4PrimaryVertex(x, y, z, 0);
        G4PrimaryParticle *particleH = NULL;

        if (fEventType == "moller")
            particleH = new G4PrimaryParticle(particleTable->FindParticle("e-"));
        else
            particleH = new G4PrimaryParticle(particleTable->FindParticle(fRecoilParticle));

        double m_h = particleH->GetParticleDefinition()->GetPDGMass();
        double p_h = sqrt(e_h * e_h - m_h * m_h);
        double kx_h = sin(theta_h) * cos(phi_h);
        double ky_h = sin(theta_h) * sin(phi_h);
        double kz_h = cos(theta_h);
        particleH->SetMomentumDirection(G4ThreeVector(kx_h, ky_h, kz_h));
        particleH->SetTotalEnergy(e_h);
        vertexH->SetPrimary(particleH);

        anEvent->AddPrimaryVertex(vertexH);

        fPID[fN] = particleH->GetPDGcode();
        fX[fN] = x;
        fY[fN] = y;
        fZ[fN] = z;
        fE[fN] = e_h;
        fMomentum[fN] = p_h;
        fTheta[fN] = theta_h;
        fPhi[fN] = phi_h;
        fN++;
    }

    if (e_p > 0) {
        G4PrimaryVertex *vertexP = new G4PrimaryVertex(x, y, z, 0);
        G4PrimaryParticle *particleP = new G4PrimaryParticle(particleTable->FindParticle("gamma"));
        double kx_p = sin(theta_p) * cos(phi_p);
        double ky_p = sin(theta_p) * sin(phi_p);
        double kz_p = cos(theta_p);
        particleP->SetMomentumDirection(G4ThreeVector(kx_p, ky_p, kz_p));
        particleP->SetTotalEnergy(e_p);
        vertexP->SetPrimary(particleP);

        anEvent->AddPrimaryVertex(vertexP);

        fPID[fN] = particleP->GetPDGcode();
        fX[fN] = x;
        fY[fN] = y;
        fZ[fN] = z;
        fE[fN] = e_p;
        fMomentum[fN] = e_p;
        fTheta[fN] = theta_p;
        fPhi[fN] = phi_p;
        fN++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRadPrimaryGenerator::DRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path) : PRadPrimaryGenerator(type, rec, par)
{
    if (path.empty()) {
        if (fEventType == "elastic")
            path = "edelastic.dat";

        else if (fEventType == "moller")
            path = "moller.dat";
    }

    if (!fParser.ReadFile(path)) {
        G4cout << "ERROR: failed to read event file " << "\"" << path << "\"" << G4endl;
        exit(1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRadPrimaryGenerator::~DRadPrimaryGenerator()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
