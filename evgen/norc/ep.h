//
// ep.h
// Developer : Chao Gu
// History:
//   Apr 2017, C. Gu, p(e,e')p elastic event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EP_h
#define EP_h 1

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <iostream>

#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.00001

#define InterpolPoints 10000
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
const double deg = TMath::Pi()/180.; // Degree to radian conversion
const double m = 0.51099893e-3; // Mass of the electron/positron (in GeV)
const double m2 = TMath::Power(m, 2);
const double M = 938.272046e-3; // Mass of the proton (in GeV)
const double M2 = TMath::Power(M, 2);
const double alpha = 1. / 137.036; // Fine-structure constant
const double e = Sqrt(4. * Pi * alpha); // Electron charge magnitude
const double mu = 2.79284736; // Magnetic moment of the proton
const double fm = 0.197327; // GeV^{-1} to fm conversion
const double mkb = 389.379404; // GeV^{-2} to mkbarn conversion

const double rp = 0.876 / fm; // GeV^{-1}
// Proton form factor
// S. Venkat, J. Arrington, G. A. Miller and X. Zhan, Phys. Rev. C, 83(2011)015203
// Electric FF
const double a11 = 2.90966;
const double a12 = -1.11542229;
const double a13 = 3.866171e-2;
const double b11 = 14.5187212;
const double b12 = 40.88333;
const double b13 = 99.999998;
const double b14 = 4.579e-5;
const double b15 = 10.3580447;
// Magnetic FF
const double a21 = -1.43573;
const double a22 = 1.19052066;
const double a23 = 2.5455841e-1;
const double b21 = 9.70703681;
const double b22 = 3.7357e-4;
const double b23 = 6.0e-8;
const double b24 = 9.9527277;
const double b25 = 12.7977739;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char* filename = "elastic_ep.dat";

double E_li; // Full energy of the initial lepton

double theta_min; // Minimum theta angle for the lepton
double theta_max; // Maximum theta angle for the lepton

double phi_min; // Minimum phi angle for the lepton
double phi_max; // Maximum phi angle for the lepton

int select;

double E_lf; // Full energy of the final lepton
double E_p; // Full energy of the final proton

double theta_l; // Theta angle of the final lepton
double theta_p; // Theta angle of the final proton

double phi_l; // Phi angle of the final lepton
double phi_p; // Phi angle of the final proton

double omega; // Solid angle (steradian)

double theta[InterpolPoints];
double xs_sin[InterpolPoints];

TLorentzVector v_li; // Initial lepton four-momentum
TLorentzVector v_pi; // Initial proton four-momentum
TLorentzVector v_lf; // Final lepton four-momentum
TLorentzVector v_pf; // Final proton four-momentum

TRandom3 *PseRan = new TRandom3();

double ElasticXS_Sin(double theta);

ROOT::Math::Functor1D Func_XS_Sin(&ElasticXS_Sin);
ROOT::Math::GSLIntegrator Integrator_XS_Sin(IntOpt);
ROOT::Math::Interpolator Interpolator_XS_Sin(InterpolPoints, InterpolType);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline double Pow2(double arg) // arg^2
{
    return TMath::Power(arg, 2);
}

inline double Pow3(double arg) // arg^3
{
    return TMath::Power(arg, 3);
}

inline double Pow4(double arg) // arg^4
{
    return TMath::Power(arg, 4);
}

inline double Pow5(double arg) // arg^5
{
    return TMath::Power(arg, 5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return ((E_li + M) * (M * E_li + m2) + Sqrt(M2 - Pow2(m * Sin(theta))) * (Pow2(E_li) - m2) * Cos(theta)) / (Pow2(E_li + M) - (Pow2(E_li) - m2) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GEp(double qq)
{
    double t = Abs(qq) / (4. * M * M); // Tau

    if (select == 1)
        return 1. - (-qq) / 6. * Pow2(rp);
    else
        return (1. + a11 * t + a12 * Pow2(t) + a13 * Pow3(t)) / (1. + b11 * t + b12 * Pow2(t) + b13 * Pow3(t) + b14 * Pow4(t) + b15 * Pow5(t));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GMp(double qq)
{
    double t = Abs(qq) / (4. * M * M); // Tau

    if (select == 1)
        return 0.;
    else
        return mu * (1. + a21 * t + a22 * Pow2(t) + a23 * Pow3(t)) / (1. + b21 * t + b22 * Pow2(t) + b23 * Pow3(t)) + b24 * Pow4(t) + b25 * Pow5(t);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS()
{
    double qq = 2. * M * (E_lf - E_li); // Four-momentum transfer squared
    double tau = -qq / (4. * M * M); // Tau
    //double eps = 1. / (1. + 2. * (1. + tau) * Pow2(Tan(theta_l / 2.))); // Epsilon
    // To calculate the Rosenbluth cross section without neglecting the lepton mass:
    double myeps = 1. / (1. - 2. * (1. + tau) * (qq + 2. * m2) / (4. * E_li * E_lf + qq)); // Modified epsilon
    double d = (E_lf / E_li) * Sqrt((Pow2(E_li) - m2) / (Pow2(E_lf) - m2));

    // NB: the lepton mass isn't neglected here, see arXiv:1401.2959
    return Pow2(alpha / (2. * E_li)) * ((1. + qq / (4. * E_li * E_lf)) / Pow2(qq / (4. * E_li * E_lf))) * (1. / d) * (M * (Pow2(E_lf) - m2) / (M * E_li * E_lf + m2 * (E_lf - E_li - M))) * (1. / (myeps * (1. + tau))) * (myeps * Pow2(GEp(qq)) + tau * Pow2(GMp(qq)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetFinalFourMomenta()
{
    v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - m2) * Sin(theta_l), 0., Sqrt(Pow2(E_lf) - m2) * Cos(theta_l), E_lf);
    v_pf = v_li + v_pi - v_lf;

    // Checking the kinematics:
    if (Abs(M2 - v_pf * v_pf) > 1.e-8) {
        std::cout << std::endl << "Warning: bad kinematics! M^2 - v_pf^2 = " << M2 - v_pf * v_pf << " GeV^2" << std::endl;
    }

    // Kinematic parameters of the final proton:
    E_p = v_pf.E();
    theta_p = v_pf.Theta();
    phi_p = v_pf.Phi();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    theta_l = theta; // Theta angle for the lepton

    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
    SetFinalFourMomenta();

    return ElasticXS() * Sin(theta_l);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg) {
        theta_l = theta_min + arg[0] * (theta_max - theta_min);
        E_lf = ElasticEnergy(theta_l);
        SetFinalFourMomenta();
        return (phi_max - phi_min) * (theta_max - theta_min) * Interpolator_XS_Sin.Eval(theta_l) * mkb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
