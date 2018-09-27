//
// newep.hh
// Developer : Chao Gu
// Based on Eur. Phys. J. A 51(2015)1
// History:
//   Apr 2017, C. Gu, ep event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NEWEP_h
#define NEWEP_h 1

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

#include "Math/Interpolator.h"

#include <gsl/gsl_integration.h>

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

const double pi = 3.14159265358979323846;
const double pi2 = pi * pi;
const double deg = pi / 180.0;
const double m = 0.5109989461e-3; // GeV
const double m2 = m * m;
const double m4 = m2 * m2;
const double M = 938.272046e-3; // GeV
const double M2 = M * M;
const double M4 = M2 * M2;
const double mmu = 105.6583745e-3; // GeV
const double mtau = 1776.82e-3; // GeV
const double alpha = 0.72973525664e-2;
const double alpha_pi = alpha / pi;
const double alpha2 = alpha * alpha;
const double alpha3 = alpha2 * alpha;
const double mu = 2.792782;
const double e = Sqrt(4.0 * pi *alpha);
const double mkb = 389.379; // MeV^{-2} to mkbarn conversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ep.dat";
const char *ifilename = "elastic_ep.info";

int Use_TPE;

double Ei_e;

double theta_min, theta_max;

double Ef_e, theta_e, phi_e;
double Ef_p, theta_p, phi_p;

TLorentzVector vi_e, vi_p;
TLorentzVector vf_e, vf_p;

double omega;

double E_g_min, E_g_cut;
double v_min, v_cut;

double E_g, theta_g, phi_g;
TLorentzVector v_g;

double theta[InterpolPoints];
double xs_elastic_sin[InterpolPoints];
double xs_born_sin[InterpolPoints];

TRandom *PseRan;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS(double s, double q2);
double NonRadXS(double s, double q2);

double BornXS_Sin(double theta);
double ElasticXS_Sin(double theta);

ROOT::Math::Interpolator Interpolator_ElasticXS_Sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator Interpolator_BornXS_Sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator TPE_Feshbach(32, InterpolType);
ROOT::Math::Interpolator TPE_Oleksandr(32, InterpolType);

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

double GEp(double q2)
{
    const double a11 = 2.90966;
    const double a12 = -1.11542229;
    const double a13 = 3.866171e-2;
    const double b11 = 14.5187212;
    const double b12 = 40.88333;
    const double b13 = 99.999998;
    const double b14 = 4.579e-5;
    const double b15 = 10.3580447;

    double t = q2 / (4.0 * M2); // Tau

    return (1.0 + a11 * t + a12 * Pow2(t) + a13 * Pow3(t)) / (1.0 + b11 * t + b12 * Pow2(t) + b13 * Pow3(t) + b14 * Pow4(t) + b15 * Pow5(t));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GMp(double q2)
{
    const double a21 = -1.43573;
    const double a22 = 1.19052066;
    const double a23 = 2.5455841e-1;
    const double b21 = 9.70703681;
    const double b22 = 3.7357e-4;
    const double b23 = 6.0e-8;
    const double b24 = 9.9527277;
    const double b25 = 12.7977739;

    double t = q2 / (4.0 * M2); // Tau

    return mu * (1.0 + a21 * t + a22 * Pow2(t) + a23 * Pow3(t)) / (1.0 + b21 * t + b22 * Pow2(t) + b23 * Pow3(t) + b24 * Pow4(t) + b25 * Pow5(t));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return ((Ei_e + M) * (M * Ei_e + m2) + Sqrt(M2 - Pow2(m * Sin(theta))) * (Pow2(Ei_e) - m2) * Cos(theta)) / (Pow2(Ei_e + M) - (Pow2(Ei_e) - m2) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS(double s, double q2)
{
    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (4)

    double theta_1B = q2 - 2.0 * m2; // eq. (7)
    double theta_2B = (s * (s - q2) - M2 * q2) / (2.0 * M2); // eq. (8)

    double t = q2 / 4.0 / M2;
    double f1 = 4.0 * t * M2 * Pow2(GMp(q2));
    double f2 = 4.0 * M2 * (Pow2(GEp(q2)) + t * Pow2(GMp(q2))) / (1.0 + t);

    double sig_born = alpha2 / lambda_s / Pow2(q2) * (theta_1B * f1 + theta_2B * f2); // eq. (3)

    return sig_born;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Finite part of the bremsstrahlung cross-section integrated over phik
double sig_rad_1(double s, double q2, double v, double tau)
{
    double x = s - q2 - v;

    double lambda_s = s * s - 4.0 * m2 * M2;
    double lambda_y = Pow2(s - x) + 4.0 * M2 * q2;

    double b1 = (-lambda_y * tau - (Pow2(s) - Pow2(x)) * tau - 2.0 * (s + x) * q2) / 2.0;
    double b2 = (-lambda_y * tau + (Pow2(s) - Pow2(x)) * tau + 2.0 * (s + x) * q2) / 2.0;
    double c1 = -(4.0 * (M2 * Pow2(tau) - (s - x) * tau - q2) * m2 - Pow2(s * tau + q2));
    double c2 = -(4.0 * (M2 * Pow2(tau) - (s - x) * tau - q2) * m2 - Pow2(x * tau - q2));
    double sc1 = Sqrt(c1);
    double sc2 = Sqrt(c2);

    double F = 1.0 / Sqrt(lambda_y);
    double F_d = ((s + x) * ((s - x) * tau + 2.0 * q2)) / (sc1 * sc2 * (sc1 + sc2));
    double F_1p = 1.0 / sc1 + 1.0 / sc2;
    double F_2p = m2 * (b2 / sc2 / c2 - b1 / sc1 / c1);
    double F_2m = m2 * (b2 / sc2 / c2 + b1 / sc1 / c1);
    double F_IR = F_2p - (q2 + 2.0 * m2) * F_d;

    double theta_11 = 4.0 * (q2 - 2.0 * m2) * F_IR;
    double theta_12 = 4.0 * F_IR * tau;
    double theta_13 = -4.0 * F - 2.0 * Pow2(tau) * F_d;
    double theta_21 = 2.0 * (s * x - M2 * q2) * F_IR / M2;
    double theta_22 = (2.0 * (s + x) * F_2m + (Pow2(s) - Pow2(x)) * F_1p + 2.0 * (s - x - 2.0 * M2 * tau) * F_IR - tau * Pow2(s + x) * F_d) / 2.0 / M2;
    double theta_23 = (4.0 * M2 * F + (4.0 * m2 + 2.0 * M2 * Pow2(tau) - (s - x) * tau) * F_d - (s + x) * F_1p) / 2.0 / M2;

    double factor = 1.0 + tau;
    double r = v / factor;

    double t = q2 + r * tau;
    double tau_t = t / 4.0 / M2;
    double F1 = 4.0 * tau_t * M2 * Pow2(GMp(t));
    double F2 = 4.0 * M2 * (Pow2(GEp(t)) + tau_t * Pow2(GMp(t))) / (1.0 + tau_t);

    double theta_1 = theta_11 / r + theta_12 + theta_13 * r;
    double theta_2 = theta_21 / r + theta_22 + theta_23 * r;

    double result = -alpha2 * alpha_pi / 4.0 / lambda_s * (theta_1 * F1 + theta_2 * F2) / t / t / factor + alpha_pi * F_IR / r * BornXS(s, q2) / factor;

    return result;
}

double func1(double x, void *params)
{
    double *p = (double *)params;

    double F1 = sig_rad_1(p[0], p[1], p[2], x);
    return F1;
}

// Finite part of the bremsstrahlung cross-section integrated over tau and phik
double sig_rad_2(double s, double q2, double v)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    double result, error;
    double par[3] = {s, q2, v};

    gsl_function F;
    F.function = &func1;
    F.params = &par;

    double x = s - q2 - v;
    double lambda_y = Pow2(s - x) + 4.0 * M2 * q2;

    double tau_max = (s - x + Sqrt(lambda_y)) / 2.0 / M2;
    double tau_min = (s - x - Sqrt(lambda_y)) / 2.0 / M2;

    gsl_integration_qags(&F, tau_min, tau_max, 0, 1e-6, 10000, w, &result, &error);

    gsl_integration_workspace_free(w);

    return result;
}

double func2(double x, void *params)
{
    double *p = (double *)params;

    double F2 = sig_rad_2(p[0], p[1], x);
    return F2;
}

// Finite part of the bremsstrahlung cross-section integrated over v, tau and phik
double sig_rad_3(double s, double q2)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    double result, error;
    double par[2] = {s, q2};

    gsl_function F;
    F.function = &func2;
    F.params = &par;

    gsl_integration_qags(&F, 0, v_min, 0, 1e-6, 10000, w, &result, &error);

    gsl_integration_workspace_free(w);

    return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double NonRadXS(double s, double q2)
{
    double sig_0 = BornXS(s, q2);

    double x = s - q2;
    double q2_m = q2 + 2.0 * m2; // eq. (27)

    double lambda_m = q2 * q2 + 4.0 * m2 * q2; // eq. (29)
    double slambda_m = Sqrt(lambda_m);
    double l_m = 1.0 / slambda_m * Log((slambda_m + q2) / (slambda_m - q2)); // eq. (28)

    double lambda_s = s * s - 4.0 * m2 * M2; // eq. (4)
    double slambda_s = Sqrt(lambda_s);
    double l_s = 1.0 / slambda_s * Log((s + slambda_s) / (s - slambda_s)); // eq. (30)

    double lambda_0x = x * x - 4.0 * m2 * M2; // eq. (32)
    double slambda_0x = Sqrt(lambda_0x);
    double l_0x = 1.0 / slambda_0x * Log((x + slambda_0x) / (x - slambda_0x)); // eq. (31)

    double a = (s * x - 2.0 * M2 * (q2 - 2.0 * m2)) / 2.0 / M2; // eq. (33)
    double b = (q2 * (s * x - M2 * q2) - m2 * q2 * (q2 + 4.0 * M2)) / M2; // eq. (34)

    // eq. (35)
    // NOTICE: typo in the expression of gamma1,2
    auto S_phi = [](const double & s, const double & lambda, const double & a, const double & b) {
        double sb = Sqrt(b);
        double slambda = Sqrt(lambda);

        double delta[5] = {0.0, 1.0, 1.0, -1.0, -1.0};
        double theta[5] = {1.0, -1.0, 1.0, -1.0, 1.0};

        double D = (s + a) * (lambda * a - s * b) + 0.25 * (lambda + b) * (lambda + b);
        double sD = Sqrt(D);
        double gamma_u = (Sqrt(b + lambda) - sb) / slambda;

        double result = 0.0;

        for (int j = 1; j <= 4; j++) {
            double a_j = s - delta[j] * slambda;
            double tau_j = -a * slambda + 0.5 * delta[j] * (b - lambda) + theta[j] * sD;

            for (int i = 1; i <= 2; i++) {
                double gamma_i = theta[i] * Pow2(sb + theta[i] * slambda) / (b - lambda);
                double gamma_1 = - (sb - slambda) * (sb - slambda) / (b - lambda);

                for (int k = 1; k <= 2; k++) {
                    double gamma_jk = -(a_j * sb - theta[k] * Sqrt(b * a_j * a_j + tau_j * tau_j)) / tau_j;
                    result += theta[i] * delta[j] * (DiLog((gamma_i - gamma_u) / (gamma_i - gamma_jk)) + DiLog((gamma_u + theta[i]) / (gamma_jk + theta[i])) - DiLog((gamma_i - gamma_1) / (gamma_i - gamma_jk)) - DiLog((gamma_1 + theta[i]) / (gamma_jk + theta[i])));
                }
            }
        }

        return result * s / 2.0 / slambda;
    };

    double v_limit = 0.99 * 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    double q2mlm = q2_m * l_m - 1.0;

    double delta_VR = 2.0 * q2mlm * Log(v_max / m / M) + 0.5 * (s * l_s + x * l_0x) + S_phi(q2_m, lambda_m, a, b) + (1.5 * q2 + 4.0 * m2) * l_m - 2.0 - q2_m / slambda_m * (0.5 * lambda_m * l_m * l_m + 2.0 * DiLog(2.0 * slambda_m / (q2 + slambda_m)) - pi2 / 2.0); // eq. (40)

    // eq. (41)
    double delta_vac = 0.0;
    double lepton_mass[3] = {m, mmu, mtau};

    for (auto &vac_m : lepton_mass) {
        double vac_m2 = vac_m * vac_m;
        double vac_slambda_m = Sqrt(q2 * q2 + 4.0 * vac_m2 * q2);
        double vac_l_m = 1.0 / vac_slambda_m * Log((vac_slambda_m + q2) / (vac_slambda_m - q2));
        delta_vac += 2.0 / 3.0 * (q2 + 2.0 * vac_m2) * vac_l_m - 10.0 / 9.0 + 8.0 / 3.0 * vac_m2 / q2 * (1.0 - 2.0 * vac_m2 * vac_l_m);
    }

    // Note: in ELRADGEN, the delta_inf (delta_infx in the code) is written
    // as (Log(q2 / m2) - 1.0) * Log(v_max * v_max / s / x), which is actually
    // the form with in URA
    double delta_inf = q2mlm * Log(v_max * v_max / s / x); // eq. (42)

    double delta_Fs = -2.0 * alpha_pi * q2mlm * Log(v_max / v_min);

    double sig_Fs = sig_rad_3(s, q2);

    double t = q2 / 4.0 / M2;
    double f1 = 4.0 * t * M2 * Pow2(GMp(q2));
    double f2 = 4.0 * M2 * (Pow2(GEp(q2)) + t * Pow2(GMp(q2))) / (1.0 + t);
    double sigma_AMM = alpha3 * m2 * l_m * (12.0 * M2 * f1 - (q2 + 4.0 * M2) * f2) / (4.0 * pi * M2 * q2 * lambda_s);

    // soft part of the TPE contribution
    double delta_2g = 0.0;

    if (Use_TPE != 0)
        delta_2g = -(Log(Ei_e / Ef_e) * Log(q2 * q2 / (4.0 * M2 * Ei_e * Ef_e)) + 2.0 * DiLog(1.0 - 0.5 * M / Ei_e) - 2.0 * DiLog(1.0 - 0.5 * M / Ef_e));

    double result = (1.0 + alpha_pi * (delta_VR + delta_vac + delta_2g - delta_inf)) * sig_0 * Exp(alpha_pi * delta_inf) + sigma_AMM + sig_0 * delta_Fs + sig_Fs;

    if (Use_TPE == 1)
        result = result * (1 + TPE_Feshbach.Eval(q2));
    else if (Use_TPE == 2)
        result = result * (1 + TPE_Oleksandr.Eval(q2));

    return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalSQ2(double theta, double &s, double &q2)
{
    Ef_e = ElasticEnergy(theta);

    vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta), Ef_e);
    vf_p = vi_e + vi_p - vf_e;

    s = 2.0 * vi_e * vi_p;
    q2 = -(vi_e - vf_e) * (vi_e - vf_e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS_Sin(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double jacob = 2.0 * Pow2(Ei_e) / Pow2(1.0 + Ei_e / M * (1.0 - Cos(theta)));

    return BornXS(s, q2) * Sin(theta) * jacob * mkb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double jacob = 2.0 * Pow2(Ei_e) / Pow2(1.0 + Ei_e / M * (1.0 - Cos(theta)));

    return NonRadXS(s, q2) * Sin(theta) * jacob * mkb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);

        return (theta_max - theta_min) * Interpolator_ElasticXS_Sin.Eval(theta_e);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double matrix_r(double s, double q2, double v, double tau, double phi)
{
    double x = s - q2 - v;
    double lambda_y = Pow2(s - x) + 4.0 * M2 * q2;

    double tau_max = (s - x + Sqrt(lambda_y)) / 2.0 / M2;
    double tau_min = (s - x - Sqrt(lambda_y)) / 2.0 / M2;

    if (tau <= tau_min || tau >= tau_max) return 0.0;

    double lambda_z = (tau - tau_min) * (tau_max - tau) * (s * x * q2 - M2 * Pow2(q2) - m2 * lambda_y);
    double slambda_z = (lambda_z > 0.0) ? Sqrt(lambda_z) : 0.0;

    double z1 = (q2 * (s + x) + tau * (s * (s - x) + 2.0 * M2 * q2) - 2.0 * M * slambda_z * Cos(phi)) / lambda_y;
    double z2 = (q2 * (s + x) + tau * (x * (s - x) - 2.0 * M2 * q2) - 2.0 * M * slambda_z * Cos(phi)) / lambda_y;

    double F = 0.5 / pi / Sqrt(lambda_y);
    double F_d = F / z1 / z2;
    double F_1p = F / z1 + F / z2;
    double F_2p = F * m2 * (1.0 / Pow2(z2) + 1.0 / Pow2(z1));
    double F_2m = F * m2 * (1.0 / Pow2(z2) - 1.0 / Pow2(z1));
    double F_IR = F_2p - (q2 + 2.0 * m2) * F_d;

    double theta_11 = 4.0 * (q2 - 2.0 * m2) * F_IR;
    double theta_12 = 4.0 * F_IR * tau;
    double theta_13 = -4.0 * F - 2.0 * Pow2(tau) * F_d;
    double theta_21 = 2.0 * (s * x - M2 * q2) * F_IR / M2;
    double theta_22 = (2.0 * (s + x) * F_2m + (Pow2(s) - Pow2(x)) * F_1p + 2.0 * (s - x - 2.0 * M2 * tau) * F_IR - tau * Pow2(s + x) * F_d) / 2.0 / M2;
    double theta_23 = (4.0 * M2 * F + (4.0 * m2 + 2.0 * M2 * Pow2(tau) - (s - x) * tau) * F_d - (s + x) * F_1p) / 2.0 / M2;

    double factor = 1.0 + tau;
    double r = v / factor;

    double t = q2 + r * tau;
    double tau_t = t / 4.0 / M2;
    double F1 = 4.0 * tau_t * M2 * Pow2(GMp(t));
    double F2 = 4.0 * M2 * (Pow2(GEp(t)) + tau_t * Pow2(GMp(t))) / (1.0 + tau_t);

    double theta_1 = theta_11 / r + theta_12 + theta_13 * r;
    double theta_2 = theta_21 / r + theta_22 + theta_23 * r;

    double result = Pow6(e) / Pow2(t) * (-4.0 * pi) * Sqrt(lambda_y) * (theta_1 * F1 + theta_2 * F2) / r; // matrix

    return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BremsIntegrand: public TFoamIntegrand
{
public:
    BremsIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);
        E_g = E_g_min + arg[1] * (E_g_cut - E_g_min);
        theta_g = arg[2] * pi;
        phi_g = arg[3] * 2.0 * pi;

        if (E_g > M * (Ei_e - m) / (M + Ei_e - Sqrt(Pow2(Ei_e) - m2) * Cos(theta_g))) return 0.0;

        double A = Sqrt(Pow2(Ei_e) - m2) * Cos(theta_e) - E_g * (Cos(theta_e) * Cos(theta_g) + Sin(theta_e) * Sin(theta_g) * Cos(phi_g));
        double B = Ei_e + M - E_g;
        double C = E_g * (Ei_e + M - Sqrt(Pow2(Ei_e) - m2) * Cos(theta_g)) - M * Ei_e - m2;

        if (m2 * (Pow2(A) - Pow2(B)) + Pow2(C) < 0.0) return 0.0;

        if (Abs(Pow2(A) - Pow2(B)) < 1.0e-12) return 0.0;

        Ef_e = (B * C - A * Sqrt(m2 * (Pow2(A) - Pow2(B)) + Pow2(C))) / (Pow2(A) - Pow2(B));

        if (Abs(A * Sqrt(Pow2(Ef_e) - m2) - B * Ef_e - C) > 1.0e-9) return 0.0;

        if (Ef_e < m || Ef_e > Ei_e - E_g) return 0.0;

        v_g.SetPxPyPzE(E_g * Sin(theta_g) * Cos(phi_g), E_g * Sin(theta_g) * Sin(phi_g), E_g * Cos(theta_g), E_g);
        vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta_e), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta_e), Ef_e);
        vf_p = vi_e + vi_p - vf_e - v_g;

        TVector3 k1 = vi_e.Vect(), k2 = vf_e.Vect();
        TVector3 q = (vi_e - vf_e).Vect(), k = v_g.Vect();
        TVector3 k1k2 = k1.Cross(k2);
        TVector3 qk = q.Cross(k);

        double s = 2.0 * Ei_e * M;
        double q2 = -(vi_e - vf_e) * (vi_e - vf_e);
        double v = (vi_e + vi_p - vf_e) * (vi_e + vi_p - vf_e) - M2;
        double tau = v_g * (vi_e - vf_e) / (v_g * vi_p);
        double phik = qk.Angle(k1k2);

        double matrix = matrix_r(s, q2, v, tau, phik);

        double xs = matrix * (1.0 / (1024.0 * Pow5(pi) * Sqrt(Pow2(Ei_e) - m2) * M)) * E_g * ((Pow2(Ef_e) - m2) / Abs(A * Ef_e - B * Sqrt(Pow2(Ef_e) - m2)));

        return (theta_max - theta_min) * (E_g_cut - E_g_min) * pi * (2.0 * pi) * (2.0 * pi) * xs * Sin(theta_e) * Sin(theta_g) * mkb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
