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
// Developer : Chao Gu, Weizhi Xiong
// History:
//   Mar 2017, C. Gu, Add for DRad configuration.
//   Apr 2017, W. Xiong, Add target thickness profile.
//   May 2017, C. Gu, Add Deuteron disintegration.
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
#include "TMath.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom2.h"
#include "TTree.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"

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
#include <vector>

static double me = 0.510998928 * MeV;
static double mmu = 105.6583745 * MeV;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator() : G4VPrimaryGenerator(), fRegistered(false), fTargetInfo(false), fTargetCenter(0), fTargetHalfL(0), fEventType(""), fRecoilOn(false), fRecoilParticle(""), fEBeam(0), fReactX(0), fReactY(0), fReactZ(0), fReactTheta(0), fReactPhi(0), fReactThetaLo(0), fReactThetaHi(0), fTargetMass(0)
{
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

PrimaryGenerator::PrimaryGenerator(G4String type, G4double e, G4double x, G4double y, G4double z, G4double theta, G4double phi, G4bool rec, G4String par) : G4VPrimaryGenerator(), fRegistered(false), fTargetInfo(false), fTargetCenter(0), fTargetHalfL(0), fEventType(type), fRecoilOn(rec), fRecoilParticle(par), fEBeam(e), fReactX(x), fReactY(y), fReactZ(z), fReactTheta(theta), fReactPhi(phi), fReactThetaLo(-1e5), fReactThetaHi(-1e5), fTargetMass(0)
{
    if (fRecoilParticle != "proton" && fRecoilParticle != "deuteron")
        fRecoilParticle = "proton";

    if (fEventType != "point" && fEventType != "elastic" && fEventType != "moller" && fEventType != "inelastic")
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

PrimaryGenerator::PrimaryGenerator(G4String type, G4double e, G4double thlo, G4double thhi, G4bool rec, G4String par) : G4VPrimaryGenerator(), fRegistered(false), fTargetInfo(false), fTargetCenter(0), fTargetHalfL(0), fEventType(type), fRecoilOn(rec), fRecoilParticle(par), fEBeam(e), fReactX(-1e5), fReactY(-1e5), fReactZ(-1e5), fReactTheta(-1e5), fReactPhi(-1e5), fReactThetaLo(thlo), fReactThetaHi(thhi), fTargetMass(0)
{
    if (fRecoilParticle != "proton" && fRecoilParticle != "deuteron")
        fRecoilParticle = "proton";

    if (fEventType != "point" && fEventType != "elastic" && fEventType != "moller" && fEventType != "inelastic")
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

        fTargetInfo = true;
    }

    double x, y, z, theta_l, phi_l;

    if (fEventType == "point") {
        x = fReactX;
        y = fReactY;
        z = fReactZ;
        theta_l = fReactTheta;
        phi_l = fReactPhi;
    } else {
        x = G4RandGauss::shoot(0, 0.08) * mm;
        y = G4RandGauss::shoot(0, 0.08) * mm;
        z = fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());
        theta_l = fReactThetaLo + (fReactThetaHi - fReactThetaLo) * G4UniformRand();
        phi_l = twopi * G4UniformRand();
    }

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
    } else if (fEventType == "point") {
        e_l = E;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::PRadPrimaryGenerator(G4String type, G4bool rec, G4String par) : PrimaryGenerator(), fEventType(type), fRecoilOn(rec), fRecoilParticle(par), fTargetProfile(NULL), fZGenerator(NULL), fPseRan(NULL), fFoamI(NULL)
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::PRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path, G4String profile): PrimaryGenerator(), fEventType(type), fRecoilOn(rec), fRecoilParticle(par), fTargetProfile(NULL), fZGenerator(NULL), fPseRan(NULL), fFoamI(NULL)
{
    if (!profile.empty())
        LoadTargetProfile(profile);

    if (path.empty()) {
        if (fEventType == "elastic")
            path = "epelastic.dat";
        else if (fEventType == "moller")
            path = "moller.dat";
    }

    // only open file, do not read the whole file into memory
    if (!fParser.OpenFile(path)) {
        G4cout << "ERROR: failed to read event file " << "\"" << path << "\"" << G4endl;
        exit(1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadPrimaryGenerator::~PRadPrimaryGenerator()
{
    if (fZGenerator)
        delete fZGenerator;

    if (fPseRan)
        delete fPseRan;

    if (fFoamI)
        delete fFoamI;
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

        fTargetInfo = true;
    }

    if (fEventType == "inelastic") {
        int pid[4];
        double p[4][3];

        while (fParser.ParseLine()) {
            if (!fParser.CheckElements(16))  continue;
            else {
                fParser >> pid[0] >> p[0][0] >> p[0][1] >> p[0][2] >> pid[1] >> p[1][0] >> p[1][1] >> p[1][2] >> pid[2] >> p[2][0] >> p[2][1] >> p[2][2] >> pid[3] >> p[3][0] >> p[3][1] >> p[3][2];
                break;
            }
        }

        double x = G4RandGauss::shoot(0, 0.08) * mm;
        double y = G4RandGauss::shoot(0, 0.08) * mm;
        double z = GenerateZ();
        G4PrimaryVertex *vertexL = new G4PrimaryVertex(x, y, z, 0);

        for (int i = 0; i < 4; i++) {
            if (p[i][2] <= 0.) continue;

            G4PrimaryParticle *particleL = new G4PrimaryParticle(pid[i], p[i][0], p[i][1], p[i][2]);
            G4PrimaryVertex *vertexL = new G4PrimaryVertex(x, y, z, 0);
            vertexL->SetPrimary(particleL);
            anEvent->AddPrimaryVertex(vertexL);

            fPID[fN] = pid[i];
            fX[fN] = x;
            fY[fN] = y;
            fZ[fN] = z;
            fE[fN] = particleL->GetTotalEnergy();
            fMomentum[fN] = particleL->GetTotalMomentum();
            fTheta[fN] = particleL->GetMomentum().theta();
            fPhi[fN] = particleL->GetMomentum().phi();
            fN++;
        }

        return;
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
    double z = GenerateZ();

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

void PRadPrimaryGenerator::LoadTargetProfile(const std::string &path)
{
    if (path.empty()) return;

    ConfigParser c_parser;
    c_parser.ReadFile(path);

    std::vector<double> z, density;
    z.clear();
    density.clear();
    double tempz, tempd;

    while (c_parser.ParseLine()) {
        if (!c_parser.CheckElements(2))
            continue;

        c_parser >> tempz >> tempd;

        if (!z.empty() && tempz * 10.0 <= z.back()) continue;

        // convert z to mm, and convert density to H_atom/cm^3
        z.push_back(tempz * 10.0);
        density.push_back(tempd * 2.0);
    }

    fTargetProfile = new ROOT::Math::Interpolator(z, density, ROOT::Math::Interpolation::kLINEAR);
    fZMin = z[0];
    fZMax = z.back();

    fPseRan = new TRandom2(0);
    fFoamI = new TargetProfileIntegrand(this);

    fZGenerator = new TFoam("Z Generator");
    fZGenerator->SetkDim(1);
    fZGenerator->SetnCells(10000); // Set number of cells
    fZGenerator->SetnSampl(500); // Set number of samples
    fZGenerator->SetOptRej(1); // Unweighted events in MC generation
    fZGenerator->SetRho(fFoamI); // Set distribution function
    fZGenerator->SetPseRan(fPseRan); // Set random number generator
    fZGenerator->SetChat(0); // Set "chat level" in the standard output
    fZGenerator->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PRadPrimaryGenerator::GenerateZ()
{
    if (fZGenerator) {
        double rvect[1];
        fZGenerator->MakeEvent();
        fZGenerator->GetMCvect(rvect);
        return fTargetCenter + fZMin + (fZMax - fZMin) * rvect[0];
    } else
        return fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TargetProfileIntegrand::TargetProfileIntegrand(PRadPrimaryGenerator *gen)
{
    fTargetProfile = gen->fTargetProfile;
    fZMin = gen->fZMin;
    fZMax = gen->fZMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double TargetProfileIntegrand::Density(int, double *arg)
{
    double z = fZMin + (fZMax - fZMin) * arg[0];

    return fTargetProfile->Eval(z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRadPrimaryGenerator::DRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path) : PRadPrimaryGenerator(type, rec, par)
{
    if (path.empty()) {
        if (fEventType == "elastic")
            path = "edelastic.dat";
        else if (fEventType == "moller")
            path = "moller.dat";
        else if (fEventType == "disintegration")
            path = "edepn.dat";
    }

    // OpenFile doesn't read the whole file into memory
    if (!fParser.OpenFile(path)) {
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

void DRadPrimaryGenerator::GeneratePrimaryVertex(G4Event *anEvent)
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

        fTargetInfo = true;
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
    double z = GenerateZ();

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

    if (fRecoilOn || fEventType == "moller" || fEventType == "disintegration") {
        G4PrimaryVertex *vertexH = new G4PrimaryVertex(x, y, z, 0);
        G4PrimaryParticle *particleH = NULL;

        if (fEventType == "moller")
            particleH = new G4PrimaryParticle(particleTable->FindParticle("e-"));
        else if (fEventType == "disintegration")
            particleH = new G4PrimaryParticle(particleTable->FindParticle("proton"));
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
        G4PrimaryParticle *particleP = NULL;

        if (fEventType == "disintegration")
            //particleP = new G4PrimaryParticle(particleTable->FindParticle("neutron"));
            return;
        else
            particleP = new G4PrimaryParticle(particleTable->FindParticle("gamma"));

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

double DRadPrimaryGenerator::GenerateZ()
{
    return fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DeuteronDisintegration::DeuteronDisintegration(G4double e, G4double eflo, G4double efhi, G4double thlo, G4double thhi) : PrimaryGenerator(), fEBeam(e), fEnpLo(eflo), fEnpHi(efhi), fReactThetaLo(thlo), fReactThetaHi(thhi)
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DeuteronDisintegration::~DeuteronDisintegration()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DeuteronDisintegration::GeneratePrimaryVertex(G4Event *anEvent)
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

        fTargetInfo = true;
    }

    double Md = particleTable->FindParticle("deuteron")->GetPDGMass();
    double Mp = particleTable->FindParticle("proton")->GetPDGMass();
    double Mn = particleTable->FindParticle("neutron")->GetPDGMass();

    double x, y, z, theta_l, phi_l;

    x = G4RandGauss::shoot(0, 0.08) * mm;
    y = G4RandGauss::shoot(0, 0.08) * mm;
    z = fTargetCenter + fTargetHalfL * 2 * (0.5 - G4UniformRand());
    theta_l = fReactThetaLo + (fReactThetaHi - fReactThetaLo) * G4UniformRand();
    phi_l = twopi * G4UniformRand();

    double E = fEBeam;
    double P = sqrt(E * E - me * me);

    double enp = fEnpLo + (fEnpHi - fEnpLo) * G4UniformRand();

    double w = enp + Mn + Mp;
    double w2 = w * w;
    double sinth22 = sin(theta_l / 2.0) * sin(theta_l / 2.0);
    double e_l = (2 * Md * E + Md * Md - w2) / (2 * Md + 4 * E * sinth22);
    double p_l = sqrt(e_l * e_l - me * me);

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

    G4ThreeVector vi_l(0, 0, P);
    G4ThreeVector vf_l(p_l * sin(theta_l), 0, p_l * cos(theta_l));
    G4ThreeVector vq = vi_l - vf_l;
    double vq2 = vq.mag2();

    double Mp2 = Mp * Mp;
    double Mn2 = Mn * Mn;
    double ecm_p = (w2 + Mp2 - Mn2) / (2 * w);
    double pcm_p = sqrt(ecm_p * ecm_p - Mp2);

    // Boost variables
    double ee = sqrt(w2 + vq2);
    double pp = sqrt(vq2);
    double gamma = ee / w;
    double v = pp / ee;

    double theta_cm = acos(-1.0 + 2.0 * G4UniformRand());
    //double theta_cm = pi * G4UniformRand();

    G4ThreeVector vf_p(pcm_p * sin(theta_cm), 0, gamma * (pcm_p * cos(theta_cm) + v * ecm_p)); // NOTE: not lab system
    double dtheta = vf_p.theta(); // angle between proton momentum and virtual photon momentum

    if (vq.theta() - dtheta > 0) {
        vf_p.setTheta(vq.theta() - dtheta);
        vf_p.setPhi(vq.phi());
    } else {
        vf_p.setTheta(dtheta - vq.theta());
        vf_p.setPhi(vq.phi() + pi);
    }

    double phi = twopi * G4UniformRand(); // angle between scattering plane and interaction plane

    vf_p.rotate(vq, phi);
    vf_p.rotateZ(phi_l);

    double p_p = vf_p.mag();
    double theta_p = vf_p.theta();
    double phi_p = vf_p.phi();
    double e_p = sqrt(Mp * Mp + vf_p.mag2());

    G4PrimaryVertex *vertexP = new G4PrimaryVertex(x, y, z, 0);
    G4PrimaryParticle *particleP = new G4PrimaryParticle(particleTable->FindParticle("proton"));
    particleP->SetMomentumDirection(vf_p.unit());
    particleP->SetTotalEnergy(e_p);
    vertexP->SetPrimary(particleP);

    anEvent->AddPrimaryVertex(vertexP);

    fPID[fN] = particleP->GetPDGcode();
    fX[fN] = x;
    fY[fN] = y;
    fZ[fN] = z;
    fE[fN] = e_p;
    fMomentum[fN] = p_p;
    fTheta[fN] = theta_p;
    fPhi[fN] = phi_p;
    fN++;

    /* Neutron
    G4ThreeVector vf_n = vq - vf_p;

    double p_n = vf_n.mag();
    double theta_n = vf_n.theta();
    double phi_n = vf_n.phi();
    double e_n = sqrt(Mn * Mn + vf_n.mag2());

    G4PrimaryVertex *vertexN = new G4PrimaryVertex(x, y, z, 0);
    G4PrimaryParticle *particleN = new G4PrimaryParticle(particleTable->FindParticle("neutron"));
    particleN->SetMomentumDirection(vf_n.unit());
    particleN->SetTotalEnergy(e_n);
    vertexN->SetPrimary(particleN);

    anEvent->AddPrimaryVertex(vertexN);

    fPID[fN] = particleN->GetPDGcode();
    fX[fN] = x;
    fY[fN] = y;
    fZ[fN] = z;
    fE[fN] = e_n;
    fMomentum[fN] = p_n;
    fTheta[fN] = theta_n;
    fPhi[fN] = phi_n;
    fN++;
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CosmicsGenerator::CosmicsGenerator() : G4VPrimaryGenerator(), fRegistered(false)
{
    fEMin = 0.5 * GeV;
    fEMax = 100.0 * GeV;
    fZenithMin = 0.0 * deg;
    fZenithMax = 120.0 * deg;

    fPseRan = new TRandom2();
    fFoamI = new CosmicsIntegrand(this, 4.28, 854.0, 174.0, 3.09);

    fETGenerator = new TFoam("ET Generator");
    fETGenerator->SetkDim(2);
    fETGenerator->SetnCells(20000); // Set number of cells
    fETGenerator->SetnSampl(1000); // Set number of samples
    fETGenerator->SetOptRej(1); // Unweighted events in MC generation
    fETGenerator->SetRho(fFoamI); // Set distribution function
    fETGenerator->SetPseRan(fPseRan); // Set random number generator
    fETGenerator->SetChat(0); // Set "chat level" in the standard output
    fETGenerator->Initialize();

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

CosmicsGenerator::~CosmicsGenerator()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CosmicsGenerator::GeneratePrimaryVertex(G4Event *anEvent)
{
    if (!fRegistered) {
        Register(gRootTree->GetTree());
        fRegistered = true;
    }

    Clear();

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    double x, y, z, theta_l, phi_l;
    double p_l, e_l;
    double rvect[2];

    double tempr, tempt, tempp;

    tempr = 2.8 * m;
    tempt = acos(G4UniformRand());
    tempp = twopi * G4UniformRand();

    x = tempr * sin(tempt) * sin(tempp);
    y = -1.4 * m + tempr * cos(tempt);
    z = 2.85 * m + tempr * sin(tempt) * cos(tempp);

    fETGenerator->MakeEvent();
    fETGenerator->GetMCvect(rvect);
    double E = fEMin + (fEMax - fEMin) * rvect[0];
    double zenith = fZenithMin + (fZenithMax - fZenithMin) * rvect[1];

    theta_l = 180.0 * deg - zenith;
    phi_l = twopi * G4UniformRand();

    e_l = E;
    p_l = sqrt(E * E - mmu * mmu);

    G4PrimaryVertex *vertexL = new G4PrimaryVertex(x, y, z, 0);
    G4PrimaryParticle *particleL = new G4PrimaryParticle(particleTable->FindParticle("mu-"));
    double kx_l = sin(theta_l) * sin(phi_l);
    double ky_l = cos(theta_l);
    double kz_l = sin(theta_l) * cos(phi_l);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CosmicsGenerator::Register(TTree *tree)
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

void CosmicsGenerator::Print() const
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

void CosmicsGenerator::Clear()
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

CosmicsIntegrand::CosmicsIntegrand(CosmicsGenerator *gen, double e0, double eps, double rd, double nn) :  E0(e0), epsilon(eps), Rd(rd), n(nn)
{
    fEMin = gen->fEMin;
    fEMax = gen->fEMax;
    fZenithMin = gen->fZenithMin;
    fZenithMax = gen->fZenithMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double CosmicsIntegrand::Density(int, double *arg)
{
    double E = (fEMin + (fEMax - fEMin) * arg[0]) / GeV;
    double Theta = fZenithMin + (fZenithMax - fZenithMin) * arg[1];

    double cosTheta = cos(Theta);

    return TMath::Power(E0 + E, -n) / (1 + E / epsilon) * TMath::Power(sqrt(Rd * Rd * cosTheta * cosTheta + 2 * Rd + 1) - Rd * cosTheta, -(n - 1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
