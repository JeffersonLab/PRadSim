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

#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <iostream>

#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.0001

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
const double mmu = 105.6583745; // MeV
const double mtau = 1776.82; // MeV
const double alpha = 1.0 / 137.035999139;
const double alp_pi = alpha / pi;
const double alp2 = alpha * alpha;
const double alp3 = alp2 * alpha;
const double mkb = 38937.9323 * 1e4; // MeV^{-2} to mkbarn conversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char *filename = "elastic_ee.dat";
const char *ifilename = "elastic_ee.info";

double Ei_1;

double theta_min, theta_max;

double v_min, v_cut;

double Ef_1, Ef_2;
double theta_1, theta_2;
double phi_1, phi_2;

TLorentzVector vi_1, vi_2;
TLorentzVector vf_1, vf_2;

double omega;

double E_g, theta_g, phi_g;
TLorentzVector v_g;

double temps, tempt, tempv, tempt1, tempz;

double theta[InterpolPoints];
double xs_elastic_sin[InterpolPoints];
double xs_brems_sin[InterpolPoints];

TRandom *PseRan;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern "C"
{
    void merad_init(double Elab);
    double merad_sigfs(double vmin, double t, double pl);
    double merad_sigfh(double vmin, double vmax, double t, double pl);
    double merad_sigr(double t, double t1, double v, double z, double pl, int ikey);
    void merad_genvt1z(double t, double vmin, double vmax, double *vgen, double *t1gen, double *zgen);
};

double BornXS(double s, double t, double u0);
double NonRadXS(double s, double t, double u0);
double RadXS(double s, double t, double u0);

double SigR(double t, double t1, double v, double z, int ikey);

void BornLevel(double s, double t, double u0, double &sig_0);
void VirtualPhoton(double s, double t, double u0, double &sig_0, double &sig_S, double &sig_vert, double &sig_B);
void InfraredDivergent(double s, double t, double u0, double &delta_1H, double &delta_1S, double &delta_1inf);

double BornXS_Sin(double theta);
double ElasticXS_Sin(double theta);
double BremsXS_Sin(double theta);

double MollerXS_Sin(double theta);

void RecBremsKins(double theta);

ROOT::Math::Functor1D Func_BornXS_Sin(&BornXS_Sin);
ROOT::Math::GSLIntegrator Integrator_BornXS_Sin(IntOpt);

ROOT::Math::Functor1D Func_MollerXS_Sin(&MollerXS_Sin);
ROOT::Math::GSLIntegrator Integrator_MollerXS_Sin(IntOpt);

ROOT::Math::Functor1D Func_ElasticXS_Sin(&ElasticXS_Sin);
ROOT::Math::GSLIntegrator Integrator_ElasticXS_Sin(IntOpt);

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

double ElasticEnergy(double theta)
{
    return m * (Ei_1 + m + (Ei_1 - m) * Pow2(Cos(theta))) / (Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS(double s, double t, double u0)
{
    double sig_0t, sig_0u;
    // t channel
    BornLevel(s, t, u0, sig_0t);
    // u0 channel
    BornLevel(s, u0, t, sig_0u);

    return (sig_0t + sig_0u) * (s - 2.0 * m2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double NonRadXS(double s, double t, double u0)
{
    // virtual photon part of Moller cross section
    // the t and u channels should be calculated together
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    VirtualPhoton(s, t, u0, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    VirtualPhoton(s, u0, t, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // infrared divergent part of real photon emission
    // the t and u channels are not separated
    double delta_1H, delta_1S, delta_1inf;
    InfraredDivergent(s, t, u0, delta_1H, delta_1S, delta_1inf);

    // the "soft" Bremsstrahlung part of the radiative cross section
    // below v_min, photon emission is not detectable
    double sig_Fs = merad_sigfs(v_min, t, 0.0);

    double result = (1.0 + alp_pi * (delta_1H + delta_1S + delta_1inf)) * (sig_0t + sig_0u) + sig_St + sig_Su + sig_vertt + sig_vertu + sig_Bt + sig_Bu;

    return result * (s - 2.0 * m2) + sig_Fs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double RadXS(double s, double t, double u0)
{
    double v_limit = 0.99 * ((s * t + Sqrt(s * (s - 4.0 * m2) * t * (t - 4.0 * m2))) / 2.0 / m2);
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;
    // the "hard" Bremsstrahlung part of the radiative cross section
    double sig_Fh = merad_sigfh(v_min, v_max, t, 0.0);

    return sig_Fh;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double SigR(double t, double t1, double v, double z, int ikey)
{
    double sigr = merad_sigr(t, t1, v, z, 0.0, ikey);

    return sigr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BornLevel(double s, double t, double u0, double &sig_0)
{
    double s2t = s * s * t, st2 = s * t * t, s3 = Pow3(s);
    double u02t = u0 * u0 * t, u0t2 = u0 * t * t, u03 = Pow3(u0), t3 = Pow3(t);

    // frequently used variables
    double xi_s = Sqrt(1.0 - 4.0 * m2 / s);
    double xi_t = Sqrt(1.0 - 4.0 * m2 / t);
    double xi_u0 = Sqrt(1.0 - 4.0 * m2 / u0);
    double xi_s2 = xi_s * xi_s, xi_s4 = xi_s2 * xi_s2;
    double xi_t2 = xi_t * xi_t, xi_t4 = xi_t2 * xi_t2;
    double xi_u02 = xi_u0 * xi_u0, xi_u04 = xi_u02 * xi_u02;

    // equation (49), Born Level
    sig_0 = (u0 * u0 / xi_s2 / 4.0 / s * (4.0 * xi_u04 - Pow2(1.0 - xi_u02) * (2.0 + t / u0)) - s * s * xi_s4 / u0) * 2.0 * pi * alp2 / t / t / s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VirtualPhoton(double s, double t, double u0, double &sig_0, double &sig_S, double &sig_vert, double &sig_B)
{
    BornLevel(s, t, u0, sig_0);

    double s2t = s * s * t, st2 = s * t * t, s3 = Pow3(s);
    double u02t = u0 * u0 * t, u0t2 = u0 * t * t, u03 = Pow3(u0), t3 = Pow3(t);

    // frequently used variables
    double xi_s = Sqrt(1.0 - 4.0 * m2 / s);
    double xi_t = Sqrt(1.0 - 4.0 * m2 / t);
    double xi_u0 = Sqrt(1.0 - 4.0 * m2 / u0);
    double xi_s2 = xi_s * xi_s, xi_s4 = xi_s2 * xi_s2;
    double xi_t2 = xi_t * xi_t, xi_t4 = xi_t2 * xi_t2;
    double xi_u02 = xi_u0 * xi_u0, xi_u04 = xi_u02 * xi_u02;

    // singularity term, appears in delta_ver and delta_box
    // we only need the divergence free part of the sigma_ver and simga_box,
    // which are obtained by substituting lambda = m so log(lambda/m) = 0
    double log_m = 0.0; // log(lambda/m), where lambda is the infinitesimal photon mass

    // other frequently used variables
    // Q^2 (-t) related, equation (27) - (29)
    double Q2_m = -t + 2.*m2;
    double lambda_m = t * t - 4.*m2 * t;
    double slambda_m = Sqrt(lambda_m);
    double L_m = 1. / slambda_m * Log((slambda_m - t) / (slambda_m + t));

    // s related
    double log_s = Log((1.0 + xi_s) / (1.0 - xi_s));
    double log_2s = Log((1.0 + xi_s) / 2.0 / xi_s);
    double Li2_s = DiLog((xi_s - 1.0) / 2.0 / xi_s);
    double Li2_sp = DiLog(2.0 * xi_s / (xi_s + 1.0));

    // t related
    double log_t = Log((1.0 + xi_t) / (xi_t - 1.0));
    double log_2t = Log((1.0 + xi_t) / 2.0);
    double log_4t = Log((xi_t2 - 1.0) / 4.0);
    double Li2_t = DiLog((1.0 - xi_t) / 2.0);

    // u0 related
    double log_u0 = Log((1.0 + xi_u0) / (xi_u0 - 1.0));
    double log_2u0 = Log((xi_u0 - 1.0) / 2.0 / xi_u0);
    double log_2u0p = Log((xi_u0 + 1.0) / 2.0 / xi_u0);
    double Li2_u0 = DiLog((1.0 + xi_u0) / 2.0 / xi_u0);

    // vacuum polarization for all leptons, factorized part
    // equation (41) with Q^2 -> -t
    double delta_vac = 0.0;
    double lepton_mass[3] = {m, mmu, mtau};

    for (auto &vac_m : lepton_mass) {
        double vac_m2 = vac_m * vac_m;
        double vac_slambda_m = Sqrt(t * t - 4.0 * vac_m2 * t);
        double vac_L_m = 1.0 / vac_slambda_m * Log((vac_slambda_m - t) / (vac_slambda_m + t));
        delta_vac += 2.0 / 3.0 * (-t + 2 * vac_m2) * vac_L_m - 10.0 / 9.0 - 8.0 / 3.0 * vac_m2 / t * (1.0 - 2.0 * vac_m2 * vac_L_m);
    }

    // equation (50)
    sig_S = alp_pi * delta_vac * sig_0;

    // vertex correction, factorized part
    // equation (36) with Q^2 -> -t
    double delta_vert = 2.0 * (Q2_m * L_m - 1.0) * log_m + (4.0 * m2 - 3.0 / 2.0 * t) * L_m - 2.0 - Q2_m / slambda_m * (lambda_m * L_m * L_m / 2.0 + 2.0 * DiLog(2.0 * slambda_m / (slambda_m - t)) - pi2 / 2.0);

    // vertex correction, non-factorized part, anomalous magnetic moment
    // euqation (52)
    double sig_AMM = - 4.0 * alp3 / st2 / xi_t * m2 * log_t * (3.0 * (s - 2.0 * m2) / u0 + (10.0 * m2 - 3.0 * u0) / (s - 4.0 * m2));

    // equation (51)
    sig_vert = 2.0 * alp_pi * delta_vert * sig_0 + sig_AMM;

    // box diagram, factorized part
    // equation (54)
    double delta_box = (1.0 + xi_s2) / xi_s * (-4.0 * log_s * log_m + log_s * log_s - 2 * pi2 + 4.0 * Li2_sp) + (1 + xi_u02) / xi_u0 * (4.0 * log_u0 * log_m - log_u0 * log_u0 + 2.0 * log_2u0p * log_2u0p - pi2 / 3.0 + 4 * Li2_u0);

    // box diagram, non-factorized part
    // equation (A.1) and (A.2)
    double sig_B1_t1 = 1.0 / 12.0 / xi_s / t * ((xi_s2 + 1.0) * (xi_s4 - 6.0 * xi_s2 - 3.0) * s2t - 2.0 * xi_s2 * Pow3(xi_s2 + 1.0) * s3 - 12.0 * xi_s2 * st2 - 4.0 * t3) * (4.0 * pi2 + 3.0 * log_s * log_s - 6.0 * log_2s * log_2s - 12.0 * Li2_s - 6.0 * log_s * Log(-xi_s2 * s / t));
    double sig_B1_t2 = 1.0 / 12.0 / Pow3(xi_t) / t * (xi_t2 * (xi_t2 - 3.0) * (3.0 * xi_t2 + 1.0) * t3 - 2.0 * u0t2 * (3.0 * Pow6(xi_t) + 2.0 * xi_t4 + 10.0 * xi_t2 - 1.0) - 4.0 * u02t * (5.0 * xi_t4 + 4.0 * xi_t2 - 1.0) - 16.0 * xi_t2 * u03) * (4.0 * pi2 - 6.0 * log_2t * log_2t + 3.0 * log_t * log_t - 12.0 * Li2_t);
    double sig_B1_t3 = 1.0 / xi_s * log_s * (xi_s2 * s + t) * ((xi_s4 + xi_s2) * s - 2.0 * (xi_s2 - 2.0) * t);
    double sig_B1_t4 = log_4t * (2.0 * t * t - (xi_s4 + xi_s2) * s * s + (3.0 * xi_s2 - 1.0) * s * t - 2.0 * s * (t + 2.0 * u0) / xi_t2);
    double sig_B1 = sig_B1_t1 + sig_B1_t2 + sig_B1_t3 + sig_B1_t4;

    double sig_B2_t1 = 1.0 / 12.0 / xi_u0 / t * (4.0 * t3 - 2.0 * u0t2 * (xi_u04 - 6.0 * xi_u02 - 1.0) + u02t * (-Pow6(xi_u0) + xi_u04 + 9.0 * xi_u02 + 7.0) + 2.0 * Pow3(xi_u02 + 1.0) * u03) * (3.0 * log_u0 * log_u0 - 6.0 * log_2u0 * log_2u0 - 12.0 * Li2_u0 + 6.0 * log_u0 * Log(xi_u02 * u0 / t) + pi2);
    double sig_B2_t2 = 1.0 / 12.0 / xi_t2 / xi_t / t * (xi_t2 * (-xi_t4 + 2.0 * xi_t2 + 3.0) * t3 + 2.0 * (Pow6(xi_t2) - 4.0 * xi_t4 + 8.0 * xi_t2 + 1.0) * u0t2 + 4.0 * (3.0 * xi_t4 + 1.0) * u02t + 16.0 * xi_t2 * u03) * (-6.0 * log_2t * log_2t + 3.0 * log_t * log_t - 12.0 * Li2_t + 4.0 * pi2);
    double sig_B2_t3 = log_4t * (2.0 * u0 / xi_t2 * (xi_t2 * t + t + 2.0 * u0) + (t - u0) * (2.0 * t + xi_u02 * u0 + u0));
    double sig_B2_t4 = - 1.0 / xi_u0 * log_u0 * (xi_u02 * (t - u0) - 2.0 * t) * (2.0 * t + xi_u02 * u0 + u0);
    double sig_B2 = sig_B2_t1 + sig_B2_t2 + sig_B2_t3 + sig_B2_t4;

    // equation (53)
    sig_B = alp_pi / 2.0 * delta_box * sig_0 + alp3 / xi_s2 / s2t / u0 * (sig_B1 + sig_B2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void InfraredDivergent(double s, double t, double u0, double &delta_1H, double &delta_1S, double &delta_1inf)
{
    // frequently used variables
    double xi_s = Sqrt(1.0 - 4.0 * m2 / s);
    double xi_t = Sqrt(1.0 - 4.0 * m2 / t);
    double xi_u0 = Sqrt(1.0 - 4.0 * m2 / u0);
    double xi_s2 = xi_s * xi_s;
    double xi_t2 = xi_t * xi_t;
    double xi_u02 = xi_u0 * xi_u0;
    double log_s = Log((1.0 + xi_s) / (1.0 - xi_s));
    double log_t = Log((1.0 + xi_t) / (xi_t - 1.0));
    double log_u0 = Log((1.0 + xi_u0) / (xi_u0 - 1.0));

    // equation (A.5) - (A.13)
    double v_limit = 0.99 * ((s * t + Sqrt(s * (s - 4.0 * m2) * t * (t - 4.0 * m2))) / 2.0 / m2); // equation (A.6)
    double v_max = (v_limit > v_min) ? v_min : v_limit;

    double z_u1 = Sqrt((xi_u02 * (v_max + u0) - v_max) / u0) / xi_u0;
    double z_u2 = Sqrt((v_max + xi_u02 * u0) / (v_max + u0)) / xi_u0;

    auto H = [v_max](const double & ch) {
        double xi_ch = Sqrt(1.0 - 4.0 * m2 / ch), xi_ch2 = xi_ch * xi_ch;
        double z_ch = xi_ch / v_max * (Sqrt(xi_ch2 * ch * ch - 2.0 * ch * v_max + v_max * v_max) - xi_ch * ch) + 1.0; // equation (A.9)
        double z_1 = 1.0 + xi_ch;
        double z_2 = Pow2(1.0 + xi_ch) / (1.0 - xi_ch);
        double z_3 = 1.0 - xi_ch;
        double z_4 = Pow2(1.0 - xi_ch) / (1.0 + xi_ch);
        double Li2_z1 = DiLog(z_ch / z_1);
        double Li2_z2 = DiLog(z_ch / z_2);
        double Li2_z3 = DiLog(z_ch / z_3);
        double Li2_z4 = DiLog(z_ch / z_4);
        // NOTICE: fabs is not in the original function
        // however, the term in log can be negative and thus result in
        // undefined behavior for real numbers
        double log_term = Log(Pow2((xi_ch + 1.0) / (xi_ch - 1.0))) * Log(fabs((Pow2(z_ch - 1.0) - xi_ch2) / (1.0 - xi_ch2)));

        return (xi_ch2 + 1.0) / 2.0 / xi_ch * (Li2_z1 + Li2_z2 - Li2_z3 - Li2_z4 - log_term);
    };

    // equation (A.3)
    double delta_1H_t1 = DiLog(4.0 * xi_u0 / Pow2(xi_u0 + 1.0)) - DiLog(-4.0 * xi_u0 / Pow2(xi_u0 - 1.0)) - 2.0 * DiLog(2.0 * xi_u0 / (xi_u0 - 1.0)) + 2.0 * DiLog(2.0 * xi_u0 / (xi_u0 + 1.0));
    double delta_1H_t2 = DiLog(2.0 * (z_u1 - 1.0) * xi_u0 / Pow2(xi_u0 - 1.0)) + DiLog(-2.0 * (z_u1 + 1.0) * xi_u0 / Pow2(xi_u0 - 1.0)) - DiLog(-2.0 * (z_u1 - 1.0) * xi_u0 / Pow2(xi_u0 + 1.0)) - DiLog(2.0 * (z_u1 + 1.0) * xi_u0 / Pow2(xi_u0 + 1.0));
    double delta_1H_t3 = 2.0 * DiLog(-(z_u2 - 1.0) * xi_u0 / (xi_u0 - 1.0)) + 2.0 * DiLog((z_u2 + 1.0) * xi_u0 / (xi_u0 - 1.0)) - 2.0 * DiLog((z_u2 + 1.0) * xi_u0 / (xi_u0 + 1.0)) - 2.0 * DiLog((1.0 - z_u2) * xi_u0 / (xi_u0 + 1.0));
    double delta_1H_t4 = 2.0 * log_u0 * Log((xi_u02 * z_u2 * z_u2 - 1.0) / (xi_u02 - 1.0));
    delta_1H = Log(1.0 + v_max / m2) + H(s) - H(t) + (xi_u02 + 1.0) / 2.0 / xi_u0 * (delta_1H_t1 + delta_1H_t2 + delta_1H_t3 + delta_1H_t4);

    // equation (A.14) - (A.15)
    auto S_phi = [](const double & s_1, const double & s_2, const double & s_3) {
        double lambda_1 = s_1 * s_1 - 16.0 * m2 * m2, slambda_1 = Sqrt(lambda_1);
        double lambda_2 = s_2 * s_2 - 16.0 * m2 * m2, slambda_2 = Sqrt(lambda_2);
        double lambda_3 = s_3 * s_3 - 16.0 * m2 * m2, slambda_3 = Sqrt(lambda_3);
        // z_u and z_d
        double z_ud[2] = {slambda_1 / slambda_2 - 1.0, (s_1 * s_2 - 4.0 * m2 * s_3) / lambda_2 - 1.0};
        // z_1, z_2, z_3, z_4
        double z[4] = {1.0 / slambda_2 *(4.0 * m2 * (s_3 - slambda_3) / (s_2 - slambda_2) - s_1 - slambda_2), 1.0 / slambda_2 *(4.0 * m2 * (s_3 + slambda_3) / (s_2 - slambda_2) - s_1 - slambda_2), 1.0 / slambda_2 *(s_1 - slambda_2 - 4.0 * m2 * (s_3 + slambda_3) / (s_2 + slambda_2)), 1.0 / slambda_2 *(s_1 - slambda_2 - 4.0 * m2 * (s_3 - slambda_3) / (s_2 + slambda_2))};
        // Sj
        double Sj[4] = {1, 1, -1, -1};
        // (-1)^(i + 1), i from 1 to 4 but index is from 0 to 3
        double Si[4] = {1, -1, 1, -1};
        // z_u term - z_d term
        double Sk[2] = {1, -1};

        double result = 0.0;

        for (int k = 0; k < 2; ++k) {
            // TODO, check with authors
            // it is noted in the reference that
            // S_phi(s1, s2, s3) == S_phi(s2, s1, s3)
            // but this part could not satisfy the relation
            double term = Log((s_2 - slambda_2) / (s_2 + slambda_2)) * Log((z_ud[k] - z[0]) * (z_ud[k] - z[2]) / (z_ud[k] - z[1]) / (z_ud[k] - z[3]));

            double sum_term = 0.0;

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    double sign = Sj[j] * Si[i];

                    if (i == j)
                        sum_term += sign * Pow2(Log(Abs(z_ud[k] - z[i]))) / 2.0;

                    else {
                        // the input to log may be negative and thus
                        // the result may be a complex number
                        // TODO check with authors for this part
                        double spence_term = DiLog((z_ud[k] - z[i]) / (z[j] - z[i]));
                        sum_term += sign * (Log(fabs(z_ud[k] - z[i])) * Log(fabs(z[i] - z[j])) - spence_term);
                    }
                }
            }

            result += s_3 / 2.0 / slambda_3 * (term + sum_term) * Sk[k];
        }

        return result;
    };

    // equation (A.4)
    double delta_1S_t1 = (xi_s2 + 1.0) / 2.0 / xi_s * (log_s * log_s + log_s + DiLog(4.0 * xi_s / Pow2(xi_s + 1.0)));
    double delta_1S_t2 = -(xi_t2 + 1.0) / 2.0 / xi_t * (log_t * log_t - log_t + DiLog(4.0 * xi_t / Pow2(xi_t + 1.0)));
    double delta_1S_t3 = -(xi_u02 + 1.0) / 2.0 / xi_u0 * (log_u0 * log_u0 - log_u0 + DiLog(4.0 * xi_u0 / Pow2(xi_u0 + 1.0)));
    double delta_1S_t4 = -S_phi(-(xi_u02 + 1.0) * u0, (xi_s2 + 1.0) * s, -(xi_t2 + 1.0) * t) + S_phi(-(xi_u02 + 1.0) * u0, -(xi_t2 + 1.0) * t, (xi_s2 + 1.0) * s) - S_phi(-(xi_t2 + 1.0) * t, (xi_s2 + 1.0) * s, -(xi_u02 + 1.0) * u0) + 1.0;
    delta_1S = delta_1S_t1 + delta_1S_t2 + delta_1S_t3 + delta_1S_t4;

    // equation (60)
    double J_0 = -2.0 * ((xi_s2 + 1.0) / xi_s * log_s - (xi_t2 + 1.0) / xi_t * log_t - (xi_u02 + 1.0) / xi_u0 * log_u0 + 2.0);

    // equation (66)
    delta_1inf = J_0 * Log(v_max / m2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Mandelstam(double theta, double &s, double &t, double &u0)
{
    Ef_1 = ElasticEnergy(theta);

    vf_1.SetPxPyPzE(Sqrt(Pow2(Ef_1) - m2) * Sin(theta), 0.0, Sqrt(Pow2(Ef_1) - m2) * Cos(theta), Ef_1);
    vf_2 = vi_1 + vi_2 - vf_1;

    s = (vi_1 + vi_2) * (vi_1 + vi_2);
    t = (vf_1 - vi_1) * (vf_1 - vi_1);
    u0 = (vf_1 - vi_2) * (vf_1 - vi_2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BornXS_Sin(double theta)
{
    double s, t, u0;
    Mandelstam(theta, s, t, u0);

    double jacob = 4.0 * m * (Ei_1 - m) / Pow2(Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta))) / 2.0 / pi;

    return BornXS(s, t, u0) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticXS_Sin(double theta)
{
    double s, t, u0;
    Mandelstam(theta, s, t, u0);

    double jacob = 4.0 * m * (Ei_1 - m) / Pow2(Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta))) / 2.0 / pi;

    return NonRadXS(s, t, u0) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double BremsXS_Sin(double theta)
{
    double s, t, u0;
    Mandelstam(theta, s, t, u0);

    double jacob = 4.0 * m * (Ei_1 - m) / Pow2(Ei_1 + m - (Ei_1 - m) * Pow2(Cos(theta))) / 2.0 / pi;

    return RadXS(s, t, u0) * jacob * mkb * Sin(theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MollerXS_Sin(double theta)
{
    double sigma_elastic_sin = ElasticXS_Sin(theta);
    double sigma_brems_sin = BremsXS_Sin(theta);

    return sigma_elastic_sin + sigma_brems_sin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RecBremsKins(double theta)
{
    double s, t, u0;
    Mandelstam(theta, s, t, u0);

    double v_limit = 0.99 * ((s * t + Sqrt(s * (s - 4.0 * m2) * t * (t - 4.0 * m2))) / 2.0 / m2);
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;

    double v = 0.0, t1 = 0.0, z = 0.0;
    merad_genvt1z(t, v_min, v_max, &v, &t1, &z);

    double sqrts = Sqrt(s);
    double lambdas = s * (s - 4.0 * m2);

    double lambda1 = Pow2(s - v) - 4 * s * m2;
    double lambda2 = 2.0 * t + s - v - 4.0 * m2;
    double lambda3 = -s * t * (s + t - v - 4.0 * m2) - m2 * Pow2(v);
    double lambda4 = s * (s - v - 4.0 * m2) - (s + v) * z;
    double lambda5 = v * z * (s - v - z) - m2 * Pow2(v + z);
    //double lambda6 = s * (v - z) - v * (v + z);
    double lambda7 = (s + 2.0 * t1 - z - 4.0 * m2) * lambda1 - lambda2 * lambda4;
    double lambda8 = 16 * lambda3 * lambda5 - lambda7 * lambda7;

    while (std::isnan(v) || std::isnan(t1) || std::isnan(z) || lambda3 < 0.0 || lambdas * lambda1 * lambda8 < 0.0) {
        merad_genvt1z(t, v_min, v_max, &v, &t1, &z);

        lambda1 = Pow2(s - v) - 4 * s * m2;
        lambda2 = 2.0 * t + s - v - 4.0 * m2;
        lambda3 = -s * t * (s + t - v - 4.0 * m2) - m2 * Pow2(v);
        lambda4 = s * (s - v - 4.0 * m2) - (s + v) * z;
        lambda5 = v * z * (s - v - z) - m2 * Pow2(v + z);
        lambda7 = (s + 2.0 * t1 - z - 4.0 * m2) * lambda1 - lambda2 * lambda4;
        lambda8 = 16 * lambda3 * lambda5 - lambda7 * lambda7;
    }

    double lambda34 = lambda3 * lambda4;
    double lambda27 = lambda2 * lambda7;
    double lambda24 = lambda2 * lambda4;

    double slambda3 = Sqrt(lambda3);
    double slambdas = Sqrt(lambdas);
    double slambdas18 = Sqrt(lambdas * lambda1 * lambda8);

    if (PseRan->Rndm() > 0.5) slambdas18 *= -1.0;

    double lambda1_slambdas3 = lambda1 * slambdas * slambda3;

    vf_1.SetPxPyPzE(slambda3 / slambdas, 0.0, sqrts * lambda2 / 2.0 / slambdas, (s - v) / 2.0 / sqrts);
    vf_2.SetPxPyPzE(-(4.0 * lambda34 - s * lambda27) / (4.0 * lambda1_slambdas3), slambdas18 / (4.0 * lambda1_slambdas3), -sqrts * (lambda7 + lambda24) / (2.0 * lambda1 * slambdas), (s - z) / (2.0 * sqrts));

    TVector3 b = (vi_1 + vi_2).BoostVector();
    vf_1.Boost(b);
    vf_2.Boost(b);

    v_g = vi_1 + vi_2 - vf_1 - vf_2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElasticIntegrand: public TFoamIntegrand
{
public:
    ElasticIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1 = theta_min + arg[0] * (theta_max - theta_min);

        return 2.0 * pi * (theta_max - theta_min) * Interpolator_ElasticXS_Sin.Eval(theta_1);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BremsIntegrand: public TFoamIntegrand
{
public:
    BremsIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_1 = theta_min + arg[0] * (theta_max - theta_min);

        return 2.0 * pi * (theta_max - theta_min) * Interpolator_BremsXS_Sin.Eval(theta_1);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

/*
        double s = (vi_1 + vi_2) * (vi_1 + vi_2);

        double t_min = -4.0 * Ei_1 * Ei_1 * Pow2(Cos(theta_max / 2.0));
        double t_max = -m2;
        double t = t_min + arg[0] * (t_max - t_min);

        double v_limit = 0.99 * ((s * t + Sqrt(s * (s - 4.0 * m2) * t * (t - 4.0 * m2))) / 2.0 / m2);
        double v_max = (v_limit > v_cut) ? v_cut : v_limit;
        double v = v_min + arg[1] * (v_max - v_min);

        double tau = v + m2;
        double t1_min = (v * (t - v) + 2.0 * m2 * t - v * Sqrt((t - v) * (t - v) - 4.0 * m2 * t)) / tau / 2.0;
        //double t1_max = (v * (t - v) + 2.0 * m2 * t + v * Sqrt((t - v) * (t - v) - 4.0 * m2 * t)) / tau / 2.0;
        double t1_max = m2 * t * t / tau / t1_min;
        double t1 = t1_min + arg[2] * (t1_max - t1_min);

        double Az = (v - t) * (v - t) - 4.0 * m2 * t;
        double Bz = - (v * (2.0 * m2 * (t + t1) + t1 * (v - t))) + s * (-t * t + t1 * v + t * (t1 + v));
        double Cz = Pow2(s * (t - t1) + t1 * v) - 4.0 * m2 * (s * (t - t1) * (t - t1) + t1 * v * v);
        double z_min = (-Bz - Sqrt(Bz * Bz - Az * Cz)) / Az;
        //double z_max = (-Bz + Sqrt(Bz * Bz - Az * Cz)) / Az;
        double z_max = Cz / Az / z_min;
        double z = z_min + arg[0] * (z_max - z_min);

        double sigr = SigR(t, t1, v, z, 0);


        double lambda1 = (s - v) * (s - v) - 4 * s * m2;
        double lambda2 = 2.0 * t + s - v - 4.0 * m2;
        double lambda3 = -s * t * (s + t - v - 4.0 * m2) - m2 * v * v;
        double lambda4 = s * (s - v - 4.0 * m2) - (s + v) * z;
        double lambda5 = v * z * (s - v - z) - m2 * (v + z) * (v + z);
        double lambda6 = s * (v - z) - v * (v + z);
        double lambda7 = (s + 2.0 * t1 - z - 4.0 * m2) * lambda1 - lambda2 * lambda4;
        double lambda8 = 16 * lambda3 * lambda5 - lambda7 * lambda7;
        double lambdas = s * (s - 4.0 * m2);

        double slambda1 = Sqrt(lambda1);
        double slambda3 = Sqrt(lambda3);
        double slambda8 = Sqrt(lambda8);
        double slambdas = Sqrt(lambdas);
        double sqrts = Sqrt(s);

        vf_1.SetPxPyPzE(slambda3 / slambdas, 0.0, sqrts * lambda2 / 2.0 / slambdas, (s - v) / (2.0 * sqrts));
        double lambda1sss3 = lambda1 * slambdas * slambda3;
        vf_2.SetPxPyPzE(-(4.0 * lambda3 * lambda4 - s * lambda2 * lambda7) / (4.0 * lambda1sss3), slambdas * slambda1 * slambda8 / (4.0 * lambda1sss3), -sqrts * (lambda7 + lambda2 * lambda4) / (2.0 * lambda1 * slambdas), (s - z) / (2.0 * sqrts));

        //if (Abs(vf_1.M() - m) > 1.0e-5 || Abs(vf_2.M() - m) > 1.0e-5) return 0.0;

        TVector3 b = (vi_1 + vi_2).BoostVector();
        vf_1.Boost(b);
        vf_2.Boost(b);
        v_g = vi_1 + vi_2 - vf_1 - vf_2;

        //if (vf_1.E() + vf_2.E() > Ei_1) return 0.0;

        theta_1 = vf_1.Theta();
        if (theta_1 > theta_max || theta_1 < theta_min) return 0.0;

        return sigr;
*/
