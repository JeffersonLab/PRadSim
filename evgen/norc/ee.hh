//
// ee.h
// Developer : Chao Gu
// History:
//   Apr 2017, C. Gu, moller event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EE_h
#define EE_h 1

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <iostream>

#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.00001

#define InterpolPoints 50000
#define InterpolType ROOT::Math::Interpolation::kCSPLINE

#define Abs   TMath::Abs
#define Exp   TMath::Exp
#define Log   TMath::Log
#define DiLog TMath::DiLog
#define Sqrt  TMath::Sqrt
#define Sin   TMath::Sin
#define Cos   TMath::Cos
#define Tan   TMath::Tan
#define ASin  TMath::ASin
#define ACos  TMath::ACos
#define ATan  TMath::ATan
#define ATan2 TMath::ATan2

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const double Pi = TMath::Pi();
const double deg = Pi / 180.; // Degree to radian conversion
const double m = 0.51099893e-3; // Mass of the electron/positron (in GeV)
const double m2 = TMath::Power(m, 2);
const double M = m;
const double alpha = 1. / 137.036; // Fine-structure constant
const double e = Sqrt(4. * Pi *alpha);   // Electron charge magnitude
const double fm = 0.197327; // GeV^{-1} to fm conversion
const double mkb = 389.379404; // GeV^{-2} to mkbarn conversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ee.dat";
const char *ifilename = "elastic_ee.info";

double Ei_1;

double theta_min, theta_max;
double phi_min, phi_max;

int fselect;

double Ef_1, Ef_2;
double theta_1, theta_2;
double phi_1, phi_2;

double omega;

double theta[InterpolPoints];
double xs_sin[InterpolPoints];

TLorentzVector vi_1, vi_2;
TLorentzVector vf_1, vf_2;

double ElasticXS_Sin(double theta);

ROOT::Math::Functor1D Func_XS_Sin(&ElasticXS_Sin);
ROOT::Math::GSLIntegrator Integrator_XS_Sin(IntOpt);
ROOT::Math::Interpolator Interpolator_XS_Sin(InterpolPoints, InterpolType);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline double Pow2(double arg) // arg^2
{
    return TMath::Power(arg, 2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return m * (Ei_1 + m + (Ei_1 - m) * Pow2(Cos(theta))) / (Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS()
{
    double s = (vi_1 + vi_2) * (vi_1 + vi_2);
    double t = (vf_1 - vi_1) * (vf_1 - vi_1);
    double u = (vf_1 - vi_2) * (vf_1 - vi_2);

    double Adir = Pow2(s - 2. * m2) + Pow2(u - 2. * m2) + 4. * m2 * t;
    double Aex = Pow2(s - 2. * m2) + Pow2(t - 2. * m2) + 4. * m2 * u;
    double Aint = -(s - 2. * m2) * (s - 6. * m2);

    return 2. * Pow2(alpha) * (Cos(theta_1) / Pow2(Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta_1)))) * (Adir / Pow2(t) + Aex / Pow2(u) - 2. * Aint / (t * u));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetFinalFourMomenta()
{
    vf_1.SetPxPyPzE(Sqrt(Pow2(Ef_1) - m2) * Sin(theta_1), 0., Sqrt(Pow2(Ef_1) - m2) * Cos(theta_1), Ef_1);
    vf_2 = vi_1 + vi_2 - vf_1;

    // Checking the kinematics:
    //if (Abs(M2 - vf_2 * vf_2) > 1.e-8)
    //    std::cout << "Warning: bad kinematics! M^2 - vf^2 = " << M2 - vf_2 * vf_2 << " GeV^2" << std::endl;

    // Kinematic parameters of the final proton:
    Ef_2 = vf_2.E();
    theta_2 = vf_2.Theta();
    phi_2 = vf_2.Phi();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    theta_1 = theta; // Theta angle for the lepton

    Ef_1 = ElasticEnergy(theta_1); // Full energy of the scattered lepton
    SetFinalFourMomenta();

    return ElasticXS() * Sin(theta_1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1 = theta_min + arg[0] * (theta_max - theta_min);
        Ef_1 = ElasticEnergy(theta_1);
        SetFinalFourMomenta();
        return (phi_max - phi_min) * (theta_max - theta_min) * Interpolator_XS_Sin.Eval(theta_1) * mkb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
