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

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <iostream>
#include <cmath>

#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.001

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
const double M = 938.272046e-3; // GeV
const double M2 = M * M;
const double mmu = 105.6583745e-3; // GeV
const double mtau = 1776.82e-3; // GeV
const double alpha = 0.72973525664e-2;
const double alpha_pi = alpha / pi;
const double alpha2 = alpha * alpha;
const double alpha3 = alpha2 * alpha;
const double mu = 2.792782;
const double mkb = 389.379; // MeV^{-2} to mkbarn conversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ep.dat";
const char *ifilename = "elastic_ep.info";

double Ei_e;

double theta_min, theta_max;

double v_min, v_cut;

double Ef_e, theta_e, phi_e;
double Ef_p, theta_p, phi_p;

TLorentzVector vi_e, vi_p;
TLorentzVector vf_e, vf_p;

double omega;

double E_g, theta_g, phi_g;
TLorentzVector v_g;

double theta[InterpolPoints];
double xs_sin[InterpolPoints];
double xs_elastic_sin[InterpolPoints];
double xs_brems_sin[InterpolPoints];

TRandom *PseRan;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern "C"
{
    void elrad_init(double Elab, double vmin);
    double elrad_sigfs(double tmin, double tmax, double q2);
    double elrad_sigfh(double tmin, double tmax, double q2, double vcut);
    void elrad_gentvp(double s, double q2, double vmin, double vmax, double *tgen, double *vgen, double *pgen);
};

double BornXS(double s, double q2);
double NonRadXS(double s, double q2);
double RadXS(double s, double q2);

double BornXS_Sin(double theta);
double ElasticXS_Sin(double theta);
double BremsXS_Sin(double theta);

double EPXS_Sin(double theta);

ROOT::Math::Functor1D Func_EPXS_Sin(&EPXS_Sin);
ROOT::Math::GSLIntegrator Integrator_EPXS_Sin(IntOpt);

ROOT::Math::Functor1D Func_ElasticXS_Sin(&ElasticXS_Sin);
ROOT::Math::GSLIntegrator Integrator_ElasticXS_Sin(IntOpt);

ROOT::Math::Interpolator Interpolator_EPXS_Sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator Interpolator_ElasticXS_Sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator Interpolator_BremsXS_Sin(InterpolPoints, InterpolType);

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

    double tau = q2 / 4.0 / M2;
    double f1 = 4.0 * tau * M2 * Pow2(GMp(q2));
    double f2 = 4.0 * M2 * (Pow2(GEp(q2)) + tau * Pow2(GMp(q2))) / (1.0 + tau);

    double sig_born = alpha2 / lambda_s / Pow2(q2) * (theta_1B * f1 + theta_2B * f2); // eq. (3)
 
    return sig_born;
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
    auto S_phi = [](const double &s, const double &lambda, const double &a, const double &b) {
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

    double v_limit = 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
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
    // as Log(q2 / m2) - 1.0) * Log(v_max * v_max / s / x), which is actually
    // the form with in URA
    double delta_inf = q2mlm * Log(v_max * v_max / s / x); // eq. (42)

    double tmin = (2.0 * M2 * q2 + v_min * (q2 + v_min - Sqrt(Pow2(q2 + v_min) + 4.0 * M2 * q2))) / 2.0 / (M2 + v_min);
    double tmax = (2.0 * M2 * q2 + v_min * (q2 + v_min + Sqrt(Pow2(q2 + v_min) + 4.0 * M2 * q2))) / 2.0 / (M2 + v_min);
    double sig_Fs = elrad_sigfs(tmin, tmax, q2);

    double delta_Fs = -2.0 * alpha_pi * q2mlm * Log(v_max / v_min);

    // Ignore sigma_AMM (eq. (38)) since it is always very small
    double result = (1.0 + alpha_pi * (delta_VR + delta_vac - delta_inf)) * sig_0 * Exp(alpha_pi * delta_inf) + sig_0 * delta_Fs + sig_Fs;

    return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double RadXS(double s, double q2)
{
    double lambda_s = s * s - 4.0 * m2 * M2;
    
    double v_limit = 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    double tmin = (2.0 * M2 * q2 + v_max * (q2 + v_max - Sqrt(Pow2(q2 + v_max) + 4.0 * M2 * q2))) / 2.0 / (M2 + v_max);
    double tmax = (2.0 * M2 * q2 + v_max * (q2 + v_max + Sqrt(Pow2(q2 + v_max) + 4.0 * M2 * q2))) / 2.0 / (M2 + v_max);
    double sig_Fh = elrad_sigfh(tmin, tmax, q2, v_max);

    return sig_Fh;
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

    return BornXS(s, q2) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double jacob = 2.0 * Pow2(Ei_e) / Pow2(1.0 + Ei_e / M * (1.0 - Cos(theta)));

    return NonRadXS(s, q2) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BremsXS_Sin(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double jacob = 2.0 * Pow2(Ei_e) / Pow2(1.0 + Ei_e / M * (1.0 - Cos(theta)));

    return RadXS(s, q2) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double EPXS_Sin(double theta)
{
    double sigma_elastic_sin = ElasticXS_Sin(theta);
    double sigma_brems_sin = BremsXS_Sin(theta);

    return sigma_elastic_sin + sigma_brems_sin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RecBremsKins(double theta)
{
    double s, q2;
    CalSQ2(theta, s, q2);

    double lambda_s = s * s - 4.0 * m2 * M2;

    double v_limit = 2.0 * q2 * (lambda_s - q2 * (s + m2 + M2)) / (q2 * (s + 2.0 * m2) + Sqrt(q2 * lambda_s * (q2 + 4.0 * m2)));
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    double t = 0.0, v = 0.0, phi = 0.0;
    elrad_gentvp(s, q2, v_min, v_max, &t, &v, &phi);

    double lambda_1 = s * (q2 + v) + 2.0 * M2 * q2;
    double lambda_2 = t * (q2 + v) + 2.0 * M2 * (q2 + t);
    double lambda_3 = t * v * (q2 - t + v) - M2 * Pow2(q2 - t);
    double lambda_4 = q2 * s * (v_limit - v);

    double lambda_q = Pow2(q2 + v) + 4.0 * M2 * q2;

    while (std::isnan(t) || std::isnan(v) || std::isnan(phi) || lambda_3 < 0 || lambda_4 < 0 || lambda_q < 0) {
        elrad_gentvp(s, q2, v_min, v_max, &t, &v, &phi);

        lambda_1 = s * (q2 + v) + 2.0 * M2 * q2;
        lambda_2 = t * (q2 + v) + 2.0 * M2 * (q2 + t);
        lambda_3 = t * v * (q2 - t + v) - M2 * Pow2(q2 - t);
        lambda_4 = q2 * s * (v_limit - v);

        lambda_q = Pow2(q2 + v) + 4.0 * M2 * q2;
    }
    
    double slambda_3 = Sqrt(lambda_3);
    double slambda_4 = Sqrt(lambda_4);
    double slambda_q = Sqrt(lambda_q);

    double cosphi = Cos(phi);
    double sinphi = Sin(phi);

    if (PseRan->Rndm() > 0.5) slambda_3 *= -1.0;

    vf_e.SetPxPyPzE(slambda_4 / s, 0.0, (s * s - lambda_1) / (2.0 * M * s), (s - q2 - v) / 2.0 / M);
    vf_p.SetPxPyPzE(
        (slambda_3 * lambda_1 * cosphi + lambda_2 * slambda_4) / (lambda_q * s),
        (slambda_3 * slambda_q * s * sinphi) / (lambda_q * s),
        (lambda_1 * lambda_2 - 4 * M2 * slambda_3 * slambda_4 * cosphi) / (2.0 * lambda_q * M * s),
        (t + 2.0 * M2) / 2.0 / M);
    v_g.SetPxPyPzE(
        (-slambda_3 * lambda_1 * cosphi + (lambda_q - lambda_2) * slambda_4) / (lambda_q * s),
        (-slambda_3 * slambda_q * s * sinphi) / (lambda_q * s),
        (lambda_1 * (lambda_q - lambda_2) + 4.0 * M2 * slambda_3 * slambda_4 * cosphi) / (2.0 * lambda_q * M * s),
        (q2 + v - t) / 2.0 / M);

    if (std::isnan(vf_p[0]) || std::isnan(vf_p[1]) || std::isnan(vf_p[2]) || std::isnan(vf_p[3])) {
        std::cout << vf_e[0] << std::endl;
        std::cout << vf_p[0] << " " << vf_p[1] << " " << vf_p[2] << " " << vf_p[3] << std::endl;
        std::cout << lambda_1 << " " << lambda_2 << " " << lambda_3 << " " << lambda_4 << " " << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);

        return 2.0 * pi * (theta_max - theta_min) * Interpolator_ElasticXS_Sin.Eval(theta_e);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BremsIntegrand: public TFoamIntegrand
{
public:
    BremsIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_e = theta_min + arg[0] * (theta_max - theta_min);

        return 2.0 * pi * (theta_max - theta_min) * Interpolator_BremsXS_Sin.Eval(theta_e);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
