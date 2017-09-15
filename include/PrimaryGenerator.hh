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
// PrimaryGenerator.hh
// Developer : Chao Gu, Weizhi Xiong
// History:
//   Mar 2017, C. Gu, Add for DRad configuration.
//   Apr 2017, W. Xiong, Add target thickness profile.
//   May 2017, C. Gu, Add Deuteron disintegration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "ConfigParser.h"

#include "TFoamIntegrand.h"
#include "Math/Interpolator.h"

#include "G4VPrimaryGenerator.hh"

#include "G4String.hh"

#define MaxN 10

class G4Event;

class TTree;
class TFoam;
class TRandom2;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGenerator : public G4VPrimaryGenerator
{
public:
    PrimaryGenerator();
    PrimaryGenerator(G4String type, G4double e, G4double x, G4double y, G4double z, G4double theta, G4double phi, G4bool rec, G4String par);
    PrimaryGenerator(G4String type, G4double e, G4double thlo, G4double thhi, G4bool rec, G4String par);
    virtual ~PrimaryGenerator();

    virtual void GeneratePrimaryVertex(G4Event *);

protected:
    void Register(TTree *);

    void Print() const;
    void Clear();

    bool fRegistered;

    int fN;
    int fPID[MaxN];
    double fX[MaxN], fY[MaxN], fZ[MaxN];
    double fE[MaxN], fMomentum[MaxN];
    double fTheta[MaxN], fPhi[MaxN];

    G4bool fTargetInfo;
    G4double fTargetCenter, fTargetHalfL;

private:
    G4String fEventType;

    G4bool fRecoilOn;
    G4String fRecoilParticle;

    G4double fEBeam;

    G4double fReactX, fReactY, fReactZ;
    G4double fReactTheta, fReactPhi;
    G4double fReactThetaLo, fReactThetaHi;

    G4double fTargetMass;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PRadPrimaryGenerator;

class TargetProfileIntegrand : public TFoamIntegrand
{
public:
    TargetProfileIntegrand(PRadPrimaryGenerator *gen);

    double Density(int nDim, double *arg);

    ROOT::Math::Interpolator *fTargetProfile;
    double fZMin, fZMax;
};

class PRadPrimaryGenerator : public PrimaryGenerator
{
    friend class TargetProfileIntegrand;

public:
    PRadPrimaryGenerator(G4String type, G4bool rec, G4String par); // DRadPrimaryGenerator uses this
    PRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path, G4String profile);
    virtual ~PRadPrimaryGenerator();

    virtual void GeneratePrimaryVertex(G4Event *);

protected:
    void LoadTargetProfile(const std::string &path);

    virtual double GenerateZ();

    G4String fEventType;

    G4bool fRecoilOn;
    G4String fRecoilParticle;

    ROOT::Math::Interpolator *fTargetProfile;
    double fZMin, fZMax;

    TFoam *fZGenerator;
    TRandom2 *fPseRan;
    TargetProfileIntegrand *fFoamI;

    ConfigParser fParser;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRadPrimaryGenerator : public PRadPrimaryGenerator
{
public:
    DRadPrimaryGenerator(G4String type, G4bool rec, G4String par, G4String path);
    virtual ~DRadPrimaryGenerator();
    
    virtual void GeneratePrimaryVertex(G4Event *);

protected:
    virtual double GenerateZ();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DeuteronDisintegration : public PrimaryGenerator
{
public:
    DeuteronDisintegration(G4double e, G4double enplo, G4double enphi, G4double thlo, G4double thhi);
    virtual ~DeuteronDisintegration();

    virtual void GeneratePrimaryVertex(G4Event *);

private:
    G4double fEBeam;

    G4double fEnpLo, fEnpHi;
    G4double fReactThetaLo, fReactThetaHi;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CosmicsGenerator;

class CosmicsIntegrand : public TFoamIntegrand
{
public:
    CosmicsIntegrand(CosmicsGenerator *gen, double e0, double eps, double rd, double nn);

    double Density(int nDim, double *arg);

    double E0, epsilon, Rd, n;
    
    double fEMin, fEMax;
    double fZenithMin, fZenithMax;
};

class CosmicsGenerator : public G4VPrimaryGenerator
{
    friend class CosmicsIntegrand;
    
public:
    CosmicsGenerator();
    virtual ~CosmicsGenerator();

    virtual void GeneratePrimaryVertex(G4Event *);

protected:
    void Register(TTree *);

    void Print() const;
    void Clear();

    bool fRegistered;

    int fN;
    int fPID[MaxN];
    double fX[MaxN], fY[MaxN], fZ[MaxN];
    double fE[MaxN], fMomentum[MaxN];
    double fTheta[MaxN], fPhi[MaxN];
         
    double fEMin, fEMax;
    double fZenithMin, fZenithMax;
    
    TFoam *fETGenerator;
    TRandom2 *fPseRan;
    CosmicsIntegrand *fFoamI;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
