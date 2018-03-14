//
// ed.h
// Developer : Chao Gu
// History:
//   Apr 2017, C. Gu, d(e,e')d elastic event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ED_h
#define ED_h 1

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
const double M = 1875.612928e-3; // Mass of the deuteron (in GeV)
const double M2 = TMath::Power(M, 2);
const double alpha = 1. / 137.036; // Fine-structure constant
const double e = Sqrt(4. * Pi *alpha);   // Electron charge magnitude
const double fm = 0.197327; // GeV^{-1} to fm conversion
const double mkb = 389.379404; // GeV^{-2} to mkbarn conversion

const double rd = 1.0 * 2.130 / fm; // GeV^{-1}
const double mepsilon = 0.936 * 0.0022246 / 0.197 / 0.197;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ed.dat";
const char *ifilename = "elastic_ed.info";

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

inline double Pow4(double arg) // arg^4
{
    return TMath::Power(arg, 4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return ((Ei_1 + M) * (M * Ei_1 + m2) + Sqrt(M2 - Pow2(m * Sin(theta))) * (Pow2(Ei_1) - m2) * Cos(theta)) / (Pow2(Ei_1 + M) - (Pow2(Ei_1) - m2) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GCd(double q2)
{
    if (fselect == 1)
        return 1. / Pow2(1. + q2 / 12. * Pow2(rd));
    else
        return Exp(-q2 / 3.5) / (1. + q2 / mepsilon) / Pow2(1. + q2 / 0.71);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GQd(double q2)
{
    if (fselect == 1)
        return 0;
    else
        return 25.8298 / 1.01 * (Exp(-q2) + 0.01 * Exp(-q2 / 100.)) / (1. + q2 / mepsilon) / Pow2(1. + q2 / 0.71);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GMd(double q2)
{
    if (fselect == 1)
        return 0;
    else
        return 1.7487 * Exp(-q2 / 2.5) / (1. + q2 / mepsilon) / Pow2(1. + q2 / 0.71);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS()
{
    double q2 = -2. * M * (Ef_1 - Ei_1);
    double eta = q2 / (4. * M2);
    double cost2 = Cos(theta_1 / 2.0);
    double mott = Pow2(alpha * cost2 / (2. * Ei_1 * (1. - cost2 * cost2)));

    double GCd2 = Pow2(GCd(q2));
    double GMd2 = Pow2(GMd(q2));
    double GQd2 = Pow2(GQd(q2));

    double A = GCd2 + 8. / 9. * Pow2(eta) * GQd2 + 2. / 3. * eta * GMd2;
    double B = 4. / 3. * eta * (1 + eta) * GMd2;

    return mott * Ef_1 / Ei_1 * (A + B * Pow2(Tan(theta_1 / 2.)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetFinalFourMomenta()
{
    vf_1.SetPxPyPzE(Sqrt(Pow2(Ef_1) - m2) * Sin(theta_1), 0., Sqrt(Pow2(Ef_1) - m2) * Cos(theta_1), Ef_1);
    vf_2 = vi_1 + vi_2 - vf_1;

    // Checking the kinematics:
    if (Abs(M2 - vf_2 * vf_2) > 1.e-8)
        std::cout << "Warning: bad kinematics! M^2 - vf^2 = " << M2 - vf_2 *vf_2 << " GeV^2" << std::endl;

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
