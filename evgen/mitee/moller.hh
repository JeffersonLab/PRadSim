//
// moller.hh
// Developer : Chao Gu
// Beased on Chao Peng's RC moller cross section code in PRadAnalyzer
// History:
//   Apr 2017, C. Gu, moller event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MOLLER_h
#define MOLLER_h 1

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TF1.h"
#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

#include <iostream>
#include <cmath>

#define InterpolPoints 1001
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

const double pi = TMath::Pi();
const double pi2 = pi * pi;
const double deg = pi / 180.0;
const double m = 0.51099893; // MeV
const double m2 = m * m;
const double m4 = m2 * m2;
const double mmu = 105.6583745; // MeV
const double alpha = 1.0 / 137.035999139;
const double e = Sqrt(4.0 * pi *alpha);
const double alp_pi = alpha / pi;
const double alp2 = alpha * alpha;
const double alp3 = alp2 * alpha;
const double mkb = 38937.9323 * 1e4; // MeV^{-2} to mkbarn conversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ee.dat";
const char *ifilename = "elastic_ee.info";

double Ei_1;

double theta_min, theta_max;

double Ef_1, theta_1, phi_1;
double Ef_2, theta_2, phi_2;

TLorentzVector vi_1, vi_2;
TLorentzVector vf_1, vf_2;

double E_cm;

TLorentzVector cm;

double E_cmp, P_cmp;

double theta_min_cm, theta_max_cm;

double Ef_1_cm, theta_1_cm, phi_1_cm;
double Ef_2_cm, theta_2_cm, phi_2_cm;

TLorentzVector vi_1_cm, vi_2_cm;
TLorentzVector vf_1_cm, vf_2_cm;
TVector3 v3f_1_cm;

double omega;

double E_g_min, E_g_max;
double E_g_frac_min, E_g_frac_max;

double E_g, theta_g, phi_g;
TLorentzVector v_g;

double E_g_cm, theta_g_cm, phi_g_cm;
TLorentzVector v_g_cm;

double theta_cm[InterpolPoints];
double xs_elastic_sin[InterpolPoints];

double SSD, TTD, UUD;

TRandom *PseRan;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS(double s, double t, double u);
double NonRadXS(double s, double t, double u);

double BornXS_Sin(double theta);
double ElasticXS_Sin(double theta);

ROOT::Math::Interpolator Interpolator_ElasticXS_Sin(InterpolPoints, InterpolType);

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

inline double Pow6(double arg) // arg^6
{
    return TMath::Power(arg, 6);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double symWeight(double th)
{
    double symFactor;

    if ((th < theta_max_cm) && (th > theta_min_cm)) symFactor = 2.0;
    else symFactor = 2.0;

    return symFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return m * (Ei_1 + m + (Ei_1 - m) * Pow2(Cos(theta))) / (Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Mandelstam(double theta, double &s, double &t, double &u)
{
    vf_1_cm.SetPxPyPzE(P_cmp * Sin(theta), 0.0, P_cmp * Cos(theta), E_cmp);
    vf_2_cm.SetPxPyPzE(-P_cmp * Sin(theta), 0.0, -P_cmp * Cos(theta), E_cmp);
    vf_1 = vf_1_cm;
    vf_2 = vf_2_cm;
    vf_1.Boost(cm.BoostVector());
    vf_2.Boost(cm.BoostVector());

    s = (vi_1_cm + vi_2_cm) * (vi_1_cm + vi_2_cm);
    t = (vf_1_cm - vi_1_cm) * (vf_1_cm - vi_1_cm);
    u = (vf_1_cm - vi_2_cm) * (vf_1_cm - vi_2_cm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double M2(double s, double t, double u)
{
    return 64.0 * pi2 * alp2 * (m4 / (t * t) * ((s * s + u * u) / (2.0 * m4) + 4.0 * t / m2 - 4.0) + m4 / (u * u) * ((s * s + t * t) / (2.0 * m4) + 4.0 * u / m2 - 4.0) + m4 / (u * t) * (s / m2 - 2.0) * (s / m2 - 6.0));
}

double BornXS(double s, double t, double u)
{
    // CM system
    return 0.5 * M2(s, t, u) / (64.0 * pi2 * s);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LoopCorrections.hh"

Double_t SoftPhoton_Moller_Integrand(Double_t *x, Double_t *par)
{
    double var = x[0];

    return ((Sqrt(SSD) * (-2.0 + TTD) * Log((Sqrt(SSD) + Sqrt(-4.0 + SSD + 4.0 * TTD * (1.0 - var) * var)) / (Sqrt(SSD) - Sqrt(-4.0 + SSD + 4.0 * TTD * (1.0 - var) * var)))) / ((1.0 - TTD * (1.0 - var) * var) * Sqrt(-4.0 + SSD + 4.0 * TTD * (1.0 - var) * var)) + (Sqrt(SSD) * (-2.0 + UUD) * Log((Sqrt(SSD) + Sqrt(-4.0 + SSD + 4.0 * UUD * (1.0 - var) * var)) / (Sqrt(SSD) - Sqrt(-4.0 + SSD + 4.0 * UUD * (1.0 - var) * var)))) / ((1.0 - UUD * (1.0 - var) * var) * Sqrt(-4.0 + SSD + 4.0 * UUD * (1.0 - var) * var)));
}

double SoftPhoton_Moller_Integral()
{
    TF1 f("Integrand", &SoftPhoton_Moller_Integrand, 0, 0.5, 0);
    ROOT::Math::WrappedTF1 wf1(f);
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
    // Set parameters of the integration
    ig.SetFunction(wf1);
    ig.SetRelTolerance(0.0001);
    return ig.Integral(0, 0.5);
}

double NonRadXS(double s, double t, double u)
{
    SSD = s / m2;
    UUD = u / m2;
    TTD = t / m2;

    double III = SoftPhoton_Moller_Integral();

    double dE = E_g_min;

    double delta = MollerLoop(SSD, TTD, UUD) + ((2.0 * alpha * (III + (2.0 * Sqrt(SSD) * Log((Sqrt(-4.0 + SSD) + Sqrt(SSD)) / 2.0)) / Sqrt(-4.0 + SSD) + 4.0 * Log(m / (2.0 * dE)) * (0.5 + ((-2.0 + SSD) * Log((Sqrt(-4.0 + SSD) + Sqrt(SSD)) / 2.0)) / (Sqrt(-4.0 + SSD) * Sqrt(SSD)) + ((-2 + TTD) * Log((Sqrt(4.0 - TTD) + Sqrt(-TTD)) / 2.0)) / (Sqrt(4.0 - TTD) * Sqrt(-TTD)) + ((-2.0 + UUD) * Log((Sqrt(4.0 - UUD) + Sqrt(-UUD)) / 2.0)) / (Sqrt(4.0 - UUD) * Sqrt(-UUD))) + ((-4.0 + 2.0 * SSD) * (pi2 / 6.0 + ((4.0 - SSD) * Pow2(Log((Sqrt(-4.0 + SSD) + Sqrt(SSD)) / 2.0))) / (-4.0 + SSD) + Log((Sqrt(-4.0 + SSD) + Sqrt(SSD)) / 2.0) * Log(-4.0 + SSD) - DiLog((-2.0 + SSD - Sqrt(-4.0 * SSD + Pow2(SSD))) / 2.0))) / Sqrt(-4 * SSD + Pow2(SSD)))) / pi);

    return BornXS(s, t, u) * (1 + delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS_Sin(double theta)
{
    // CM system
    double s, t, u;
    Mandelstam(theta, s, t, u);

    return BornXS(s, t, u) * symWeight(vf_2_cm.Theta()) * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    // CM system
    double s, t, u;
    Mandelstam(theta, s, t, u);

    return NonRadXS(s, t, u) * symWeight(vf_2_cm.Theta()) * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1_cm = theta_min_cm + arg[0] * (theta_max_cm - theta_min_cm);

        return 2.0 * pi * (theta_max_cm - theta_min_cm) * Interpolator_ElasticXS_Sin.Eval(theta_1_cm);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MSqBrem.hh"

double bremXS(double xi)
{
    double bCS = E_g_cm / (2.0 * E_cmp * P_cmp * xi) * Pow2(vf_1_cm.P());

    return 0.5 * bCS * Mh2(vi_1_cm, vi_2_cm, vf_1_cm, v_g_cm) / (32.0 * Pow5(2 * pi) * m2);
}

class BremsIntegrandM: public TFoamIntegrand
{
public:
    BremsIntegrandM() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1_cm = theta_min_cm + arg[0] * (theta_max_cm - theta_min_cm);
        E_g_cm = E_g_min + arg[1] * (E_g_max - E_g_min);
        theta_g_cm = arg[2] * pi;
        phi_g_cm = arg[3] * 2.0 * pi;

        double E0_g_cm = 2.0 * E_cmp * (E_cmp - m) / (2.0 * E_cmp - m);

        v3f_1_cm.SetXYZ(Sin(theta_1_cm), 0.0, Cos(theta_1_cm));
        v_g_cm.SetPxPyPzE(E_g_cm * Sin(theta_g_cm) * Cos(phi_g_cm), E_g_cm * Sin(theta_g_cm) * Sin(phi_g_cm), E_g_cm * Cos(theta_g_cm), E_g_cm);

        double cosa = Cos(v3f_1_cm.Angle(v_g_cm.Vect()));
        double xi = Sqrt(4.0 * Pow2(E_cmp * (E_cmp - E_g_cm)) / m2 - Pow2(2.0 * E_cmp - E_g_cm) + Pow2(E_g_cm * cosa)) / m;

        if ((E_g_cm > E0_g_cm) && (cosa > -Sqrt(Pow2(2.0 * E_cmp - E_g_cm) - 4.0 * Pow2(E_cmp * (E_cmp - E_g_cm)) / m2) / E_g_cm)) return 0.0;

        Ef_1_cm = (2.0 * E_cmp * (E_cmp - E_g_cm) * (2.0 * E_cmp - E_g_cm) - m2 * E_g_cm * xi * cosa) / (Pow2(2.0 * E_cmp - E_g_cm) - Pow2(E_g_cm * cosa));
        vf_1_cm.SetPxPyPzE(Sqrt(Pow2(Ef_1_cm) - m2) * Sin(theta_1_cm), 0.0, Sqrt(Pow2(Ef_1_cm) - m2) * Cos(theta_1_cm), Ef_1_cm);

        vf_2_cm = vi_1_cm + vi_2_cm - vf_1_cm - v_g_cm;

        double aFlag = 1.0;

        if (E_g_cm > E0_g_cm) aFlag = 2.0;

        double xs = bremXS(xi) * symWeight(vf_2_cm.Theta()) * aFlag;

        return 2.0 * pi * (theta_max_cm - theta_min_cm) * 2.0 * pi * pi * (E_g_max - E_g_min) * Sin(theta_1_cm) * Sin(theta_g_cm) * xs * mkb;
    }
};

class BremsIntegrandP: public TFoamIntegrand
{
public:
    BremsIntegrandP() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1_cm = theta_min_cm + arg[0] * (theta_max_cm - theta_min_cm);
        E_g_cm = E_g_min + arg[1] * (E_g_max - E_g_min);
        theta_g_cm = arg[2] * pi;
        phi_g_cm = arg[3] * 2.0 * pi;

        double E0_g_cm = 2.0 * E_cmp * (E_cmp - m) / (2.0 * E_cmp - m);

        if (E_g_cm < E0_g_cm) return 0.0;

        v3f_1_cm.SetXYZ(Sin(theta_1_cm), 0.0, Cos(theta_1_cm));
        v_g_cm.SetPxPyPzE(E_g_cm * Sin(theta_g_cm) * Cos(phi_g_cm), E_g_cm * Sin(theta_g_cm) * Sin(phi_g_cm), E_g_cm * Cos(theta_g_cm), E_g_cm);

        double cosa = Cos(v3f_1_cm.Angle(v_g_cm.Vect()));
        double xi = Sqrt(4.0 * Pow2(E_cmp * (E_cmp - E_g_cm)) / m2 - Pow2(2.0 * E_cmp - E_g_cm) + Pow2(E_g_cm * cosa)) / m;

        if (cosa > -Sqrt(Pow2(2.0 * E_cmp - E_g_cm) - 4.0 * Pow2(E_cmp * (E_cmp - E_g_cm)) / m2) / E_g_cm) return 0.0;

        Ef_1_cm = (2.0 * E_cmp * (E_cmp - E_g_cm) * (2.0 * E_cmp - E_g_cm) + m2 * E_g_cm * xi * cosa) / (Pow2(2.0 * E_cmp - E_g_cm) - Pow2(E_g_cm * cosa));
        vf_1_cm.SetPxPyPzE(Sqrt(Pow2(Ef_1_cm) - m2) * Sin(theta_1_cm), 0.0, Sqrt(Pow2(Ef_1_cm) - m2) * Cos(theta_1_cm), Ef_1_cm);

        vf_2_cm = vi_1_cm + vi_2_cm - vf_1_cm - v_g_cm;

        double aFlag = 1.0;

        if (E_g_cm > E0_g_cm) aFlag = 2.0;

        double xs = bremXS(xi) * symWeight(vf_2_cm.Theta()) * aFlag;

        return 2.0 * pi * (theta_max_cm - theta_min_cm) * 2.0 * pi * pi * (E_g_max - E_g_min) * Sin(theta_1_cm) * Sin(theta_g_cm) * xs * mkb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
