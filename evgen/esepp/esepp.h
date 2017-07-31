//----------------------------------------------------------------------------------------------
// This file is part of ESEPP (version 1.4).
//
// ESEPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ESEPP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ESEPP.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (c) Alexander Gramolin, 2014. E-mail: gramolin (at) inp.nsk.su 
// http://gramolin.com/esepp/
//==============================================================================================


#include <Riostream.h>
#include <sstream>
#include <time.h>

#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h" // http://root.cern.ch/root/html/TRandom3.html

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TLorentzVector.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

using namespace std;

// Parameters for numerical integration:
#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.00001

// Parameters for interpolation:
#define InterpolPoints 1000
#define InterpolType ROOT::Math::Interpolation::kCSPLINE

// Short notations for mathematical functions:
#define Abs   TMath::Abs
#define Exp   TMath::Exp
#define Log   TMath::Log
#define DiLog TMath::DiLog
#define Sqrt  TMath::Sqrt
#define Sin   TMath::Sin
#define Cos   TMath::Cos
#define Tan   TMath::Tan
#define ASin  TMath::ASin
#define ATan  TMath::ATan


// The random number generator:
TRandom3 *PseRan = new TRandom3();

time_t starttime, stoptime; // For the timer

TString filename;

FILE *fvpol;    // File "vpol.dat" (for calculation of the vacuum polarization)

FILE *fe;       // Text file for e-/mu- events
FILE *fp;       // Text file for e+/mu+ events
FILE *f0;       // Text file for Rosenbluth events

TNtuple *ntp_e; // Ntuple for e-/mu- events
TNtuple *ntp_p; // Ntuple for e+/mu+ events
TNtuple *ntp_0; // Ntuple for Rosenbluth events

TFile *froot_e; // ROOT file for e-/mu- events
TFile *froot_p; // ROOT file for e+/mu+ events
TFile *froot_0; // ROOT file for Rosenbluth events

int nCells_1D; // Number of cells for 1D case
int nSampl_1D; // Number of samples for 1Dcase
int nCells_4D; // Number of cells for 4D case
int nSampl_4D; // Number of samples for 4D case

// Flags:
int flag_lepton;
int flag_mode;
int flag_struct;
int flag_phi;
int flag_cell;
int flag_tpe;
int flag_vpol;
bool flag_quick;
bool flag_init;
bool flag_info;
bool flag_dat;
bool flag_root;
bool flag_vepp;
bool flag_target;
bool flag_warn;
bool flag_rosen;
bool flag_bint;

double m; // The lepton mass (m = m_e OR m = m_mu)
double m2, m4; // Powers of m (for convenience)

double E_li; // Full energy of the initial lepton

// Energies of the final particles:
double E_lf; // Full energy of the final lepton
double E_p;  // Full energy of the final proton
double E_g;  // Energy of the photon

// Theta angles ("polar") of the final particles:
double theta_l; // Theta angle of the final lepton
double theta_p; // Theta angle of the final proton
double theta_g; // Theta angle of the photon

// Phi angles ("azimuthal") of the final particles:
double phi_l; // Phi angle of the final lepton
double phi_p; // Phi angle of the final proton
double phi_g; // Phi angle of the photon

double theta_min; // Minimum theta angle for the lepton
double theta_max; // Maximum theta angle for the lepton

double phi_min; // Minimum phi angle for the lepton
double phi_max; // Maximum phi angle for the lepton

double omega; // Soid angle (steradian)

// Cut and maximum energies for the photon (E_g_cut < E_g < E_g_max):
double E_g_cut, E_g_max;

// For the evaluation of the final lepton energy:
double A, B, C; // Coefficients of the quadratic equation
double en_sign; // Sign to be changed to get two different roots

// Some kinematic variables:
double qq;     // Four-momentum transfer
double tau, t; // tau = -q^2 / (4*M^2)
double eps;    // Virtual photon polarization
double q_12;   // Four-momentum transfer squared (when the lepton emits the photon)
double q_22;   // Four-momentum transfer squared (when the proton emits the photon)

// Some variables to calculate the cross section of elastic scattering:
double myeps;
double d;

// Form factors of the proton:
double F10, F20; // F1(0) and F2(0)
double F11, F12; // F1(q_12) and F1(q_22)
double F21, F22; // F2(q_12) and F2(q_22)

// Storage cell parameters:
double cell_zmin; // Minimum z-coordinate of the storage cell
double cell_zmax; // Maximum z-coordinate of the storage cell
double cell_zinj; // Z-coordinate of the injection point
double zcoord;    // Z-coordinate of the event vertex

// Some auxiliary variables:
double delta_sum;
double M_sum;

double xi; // A variable to calculate the vacuum polarization

double ret;

double res, res1, res2, res3, res4;

// Some four-momenta:
TLorentzVector v_li; // Initial lepton four-momentum
TLorentzVector v_pi; // Initial proton four-momentum
TLorentzVector v_lf; // Final lepton four-momentum
TLorentzVector v_pf; // Final proton four-momentum
TLorentzVector v_kf; // Photon four-momentum

TLorentzVector p_x; // LorentzVector for numerical integration

// Scalar products and their squares:
double kfli; // Photon, lepton initial
double kflf; // Photon, lepton final
double kfpi; // Photon, proton initial
double kfpf; // Photon, proton final
double lilf; // Lepton initial, lepton final
double lipi; // Lepton initial, proton initial
double lipf; // Lepton initial, proton final
double pipf; // Proton initial, proton final
double lfpi; // Lepton final, proton initial
double lfpf; // Lepton final, proton final

// Variables for integration of the bremsstrahlung cross section:
double bre1, bre2, bre_error, brp1, brp2, brp_error;
double sum_e, sum_p;

// Beam current integral required for VEPP-3 (kilocoulomb);
double kiloc;

long nevents; // Total number of events
long ninit;   // Number of events for initialization
long loop;    // Counter variable
int i;

// Event counters:
long count_ElE = 0; // Elastic scattering, e-/mu-
long count_ElP = 0; // Elastic scattering, e+/mu+
long count_BrE1 = 0; // Bremsstrahlung, e-/mu-, root "-"
long count_BrE2 = 0; // Bremsstrahlung, e-/mu-, root "+"
long count_BrP1 = 0; // Bremsstrahlung, e+/mu+, root "-"
long count_BrP2 = 0; // Bremsstrahlung, e+/mu+, root "+"

// Some functions: -----------------------------------------------------------------------------
double f_li_lf(double x);
double f_li_pi(double x);
double f_li_pf(double x);
double f_lf_pi(double x);
double f_lf_pf(double x);
double f_pi_pf(double x);

double f_el_e(double theta_l);
double f_el_p(double theta_l);

double f_ros(double theta_l);

ROOT::Math::Functor1D func_li_lf(&f_li_lf);
ROOT::Math::Functor1D func_li_pi(&f_li_pi);
ROOT::Math::Functor1D func_li_pf(&f_li_pf);
ROOT::Math::Functor1D func_lf_pi(&f_lf_pi);
ROOT::Math::Functor1D func_lf_pf(&f_lf_pf);
ROOT::Math::Functor1D func_pi_pf(&f_pi_pf);

ROOT::Math::Functor1D func_el_e(&f_el_e);
ROOT::Math::Functor1D func_el_p(&f_el_p);

ROOT::Math::Functor1D func_ros(&f_ros);

// Integrators: --------------------------------------------------------------------------------
ROOT::Math::GSLIntegrator i_li_lf(IntOpt);
ROOT::Math::GSLIntegrator i_li_pi(IntOpt);
ROOT::Math::GSLIntegrator i_li_pf(IntOpt);
ROOT::Math::GSLIntegrator i_lf_pi(IntOpt);
ROOT::Math::GSLIntegrator i_lf_pf(IntOpt);
ROOT::Math::GSLIntegrator i_pi_pf(IntOpt);

ROOT::Math::GSLIntegrator i_el_e(IntOpt);
ROOT::Math::GSLIntegrator i_el_p(IntOpt);

ROOT::Math::GSLIntegrator i_ros(IntOpt);

// Interpolators: ------------------------------------------------------------------------------
ROOT::Math::Interpolator inter_vpol(10000, InterpolType);
ROOT::Math::Interpolator inter_ros_sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_ee(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_ep(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_pp(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_virt(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_prime(InterpolPoints, InterpolType);

// Arrays for interpolation:
double s[10000];
double rep[10000];
double xx[InterpolPoints];
double y_ros_sin[InterpolPoints];
double y_brem_ee[InterpolPoints];
double y_brem_ep[InterpolPoints];
double y_brem_pp[InterpolPoints];
double y_virt[InterpolPoints];
double y_prime[InterpolPoints];


//----------------------------------------------------------------------------------------------
// Raising to different powers:
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
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Form factors of the proton:

// Electric form factor: -----------------------------------------------------------------------
inline double G_E(double qq)
  {
  t = Abs(qq)/(4.*M*M); // Tau

  if (flag_struct == 1) return 1.; // Point-like proton
  
  if (flag_struct == 3) // Kelly parametrization (J.J. Kelly, PRC 70 (2004) 068202)
    {
    return (1. + a11_K*t)/(1. + b11_K*t + b12_K*Pow2(t) + b13_K*Pow3(t));
    }
    
  if (flag_struct == 4) // Puckett parametrization (A.J.R. Puckett, arXiv:1008.0855)
    {
    return (1. + a11_P*t)/(1. + b11_P*t + b12_P*Pow2(t) + b13_P*Pow3(t));
    }
    
  if (flag_struct == 5) // Arbitrary parametrization from the file "const.h"
    {
      return (1. + a11*t + a12*Pow2(t) + a13*Pow3(t))/(1. + b11*t + b12*Pow2(t) + b13*Pow3(t) + b14*Pow4(t) + b15*Pow5(t));
    }

  return 1./Pow2(1. + Abs(qq)/0.71); // Dipole formula
  }

// Magnetic form factor: -----------------------------------------------------------------------
inline double G_M(double qq)
  {
  t = Abs(qq)/(4.*M*M); // Tau

  if (flag_struct == 1) return mu; // Point-like proton

  if (flag_struct == 3) // Kelly parametrization (J.J. Kelly, PRC 70 (2004) 068202)
    {
    return mu*(1. + a21_K*t)/(1. + b21_K*t + b22_K*Pow2(t) + b23_K*Pow3(t));
    }
    
  if (flag_struct == 4) // Puckett parametrization (A.J.R. Puckett, arXiv:1008.0855)
    {
    return mu*(1. + a21_P*t)/(1. + b21_P*t + b22_P*Pow2(t) + b23_P*Pow3(t));
    }
    
  if (flag_struct == 5) // Arbitrary parametrization from the file "const.h"
    {
        return mu*(1. + a21*t + a22*Pow2(t) + a23*Pow3(t))/(1. + b21*t + b22*Pow2(t) + b23*Pow3(t) + b24*Pow4(t) + b25*Pow5(t));
    }
    
  return mu*G_E(qq); // Dipole formula
  }

// Dirac form factor F1: -----------------------------------------------------------------------
inline double F1(double qq)
  {
  if (flag_struct == 1) return 1.; // Point-like proton
  
  t = Abs(qq)/(4.*M*M); // Tau
  return (G_E(qq) + t*G_M(qq))/(1. + t);
  }

// Pauli form factor F2: -----------------------------------------------------------------------
inline double F2(double qq)
  {
  if (flag_struct == 1) return mu - 1.; // Point-like proton
  
  t = Abs(qq)/(4.*M*M); // Tau
  return (G_M(qq) - G_E(qq))/(1. + t);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function EvalKinematicParams (evaluation of some kinematic parameters):
void EvalKinematicParams(double E_gamma)
  {
  qq = 2.*M*(E_lf - E_li + E_gamma); // Four-momentum transfer squared
  tau = -qq/(4.*M*M); // Tau
  eps = 1./(1. + 2.*(1. + tau)*Pow2(Tan(theta_l/2.))); // Epsilon
  // To calculate the Rosenbluth cross section without neglecting the lepton mass:
  myeps = 1./(1. - 2.*(1. + tau)*(qq + 2.*m*m)/(4.*E_li*E_lf + qq)); // Modified epsilon
  d = (E_lf/E_li)*Sqrt((Pow2(E_li) - Pow2(m))/(Pow2(E_lf) - Pow2(m)));
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Soft-photon bremsstrahlung cross section:
double BremSoftLeptonOnly() // Lepton term
  {
  return -(alpha*v_kf.E()/(4.*Pi*Pi))*(Pow2(m)/Pow2(v_kf*v_li) - 2.*(v_li*v_lf)/((v_kf*v_li)*(v_kf*v_lf)) + Pow2(m)/Pow2(v_kf*v_lf));
  }

double BremSoftProtonOnly() // Proton term
  {
  return -(alpha*v_kf.E()/(4.*Pi*Pi))*(Pow2(M)/Pow2(v_kf*v_pi) - 2.*(v_pi*v_pf)/((v_kf*v_pi)*(v_kf*v_pf)) + Pow2(M)/Pow2(v_kf*v_pf));
  }
  
double BremSoftInterference() // Interference term
  {
  return -(alpha*v_kf.E()/(2.*Pi*Pi))*(-(v_li*v_pi)/((v_kf*v_li)*(v_kf*v_pi)) + (v_li*v_pf)/((v_kf*v_li)*(v_kf*v_pf)) + (v_lf*v_pi)/((v_kf*v_lf)*(v_kf*v_pi)) - (v_lf*v_pf)/((v_kf*v_lf)*(v_kf*v_pf)));
  }
//==============================================================================================
 

//----------------------------------------------------------------------------------------------
// Functions for numerical integration (soft-photon bremsstrahlung):
double f_li_lf(double x) // For the calculation of B(v_li, v_lf, E_g_cut)
  {
  p_x = x*v_li + (1. - x)*v_lf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_li_pi(double x) // For the calculation of B(v_li, v_pi, E_g_cut)
  {
  p_x = x*v_li + (1. - x)*v_pi;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_li_pf(double x) // For the calculation of B(v_li, v_pf, E_g_cut)
  {
  p_x = x*v_li + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_lf_pi(double x) // For the calculation of B(v_lf, v_pi, E_g_cut)
  {
  p_x = x*v_lf + (1. - x)*v_pi;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_lf_pf(double x) // For the calculation of B(v_lf, v_pf, E_g_cut)
  {
  p_x = x*v_lf + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_pi_pf(double x) // For the calculation of B(v_pi, v_pf, E_g_cut)
  {
  p_x = x*v_pi + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Integration of the soft-photon bremsstrahlung cross section:
inline double d_li_li() // B(v_li, v_li, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/m) + E_li*Log(m/(E_li + Sqrt(Pow2(E_li) - Pow2(m))))/Sqrt(Pow2(E_li) - Pow2(m)))/Pi;
  }

inline double d_li_lf() // B(v_li, v_lf, E_g_cut)
  {
  return (v_li*v_lf)*i_li_lf.Integral(0., 1.)/(4.*Pi);
  }

inline double d_lf_lf() // B(v_lf, v_lf, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/m) + E_lf*Log(m/(E_lf + Sqrt(Pow2(E_lf) - Pow2(m))))/Sqrt(Pow2(E_lf) - Pow2(m)))/Pi;
  }

inline double d_li_pi() // B(v_li, v_pi, E_g_cut)
  {
  return (v_li*v_pi)*i_li_pi.Integral(0., 1.)/(4.*Pi);
  }

inline double d_li_pf() // B(v_li, v_pf, E_g_cut)
  {
  return (v_li*v_pf)*i_li_pf.Integral(0., 1.)/(4.*Pi);
  }

inline double d_lf_pi() // B(v_lf, v_pi, E_g_cut)
  {
  return (v_lf*v_pi)*i_lf_pi.Integral(0., 1.)/(4.*Pi);
  }

inline double d_lf_pf() // B(v_lf, v_pf, E_g_cut)
  {
  return (v_lf*v_pf)*i_lf_pf.Integral(0., 1.)/(4.*Pi);
  }

inline double d_pi_pi() // B(v_pi, v_pi, E_g_cut)
  {
  return (Log(2.*E_g_cut/M) - 1.)/(2.*Pi);
  }

inline double d_pi_pf() // B(v_pi, v_pf, E_g_cut)
  {
  return (v_pi*v_pf)*i_pi_pf.Integral(0., 1.)/(4.*Pi);
  }

inline double d_pf_pf() // B(v_pf, v_pf, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/M) + E_p*Log(M/(E_p + Sqrt(Pow2(E_p) - Pow2(M))))/Sqrt(Pow2(E_p) - Pow2(M)))/Pi;
  }

// Lepton term:
inline double d_brem_ee()
  {
  return -2.*alpha*(d_li_li() - 2.*d_li_lf() + d_lf_lf());
  }

// Interference term:
inline double d_brem_ep()
  {
  return 4.*alpha*(d_li_pi() - d_li_pf() - d_lf_pi() + d_lf_pf());
  }

// Proton term:
inline double d_brem_pp()
  {
  return -2.*alpha*(d_pi_pi() - 2.*d_pi_pf() + d_pf_pf());
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Virtual-photon corrections:

// Vacuum polarization contribution from electron-positron loops:
double d_vac_e()
  {
  return (2.*alpha/(3.*Pi))*(-5./3. + Log(-qq/Pow2(m_e)));
  }

// Vacuum polarization contribution from muon loops:
double d_vac_mu()
  {
  xi = -(4.*Pow2(m_mu)/qq)/Pow2(1. + Sqrt(1. - 4.*Pow2(m_mu)/qq));
  return (2*alpha/(3.*Pi))*(-5./3. + 4.*xi/Pow2(1. - xi) - ((1. - 4.*xi + Pow2(xi))/Pow2(1. - xi))*((1. + xi)/(1. - xi))*Log(xi));
  }
  
// Vacuum polarization contribution from tau loops:
double d_vac_tau()
  {
  xi = -(4.*Pow2(m_tau)/qq)/Pow2(1. + Sqrt(1. - 4.*Pow2(m_tau)/qq));
  return (2*alpha/(3.*Pi))*(-5./3. + 4.*xi/Pow2(1. - xi) - ((1. - 4.*xi + Pow2(xi))/Pow2(1. - xi))*((1. + xi)/(1. - xi))*Log(xi));
  }
  
// Electron vertex correction:
double d_vertex()
  {
  return (alpha/Pi)*(3.*Log(-qq/Pow2(m))/2. - 2.);
  }

// Two-photon exchange contribution:
double d_prime()
  {
  if (flag_tpe == 1) return 0.; // Approach of Mo & Tsai to the TPE diagrams
    
  // Approach of Maximon & Tjon to the TPE diagrams (see arXiv:1401.2959):    
  return -(alpha/Pi)*(Log(E_li/E_lf)*Log(qq*qq/(4.*M*M*E_li*E_lf)) + 2.*DiLog(1. - 0.5*M/E_li) - 2.*DiLog(1. - 0.5*M/E_lf));
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function ElasticEnergy (final lepton energy in the case of purely elastic scattering):
double ElasticEnergy(double theta)
  {
  return ((E_li + M)*(M*E_li + Pow2(m)) + Sqrt(Pow2(M) - Pow2(m*Sin(theta)))*(Pow2(E_li) - Pow2(m))*Cos(theta))/(Pow2(E_li + M) - (Pow2(E_li) - Pow2(m))*Pow2(Cos(theta)));
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function EvalEnergy (exact evaluation of the final lepton energy):
int EvalEnergy()
  {
  // Calculation of the coefficients A, B and C:
  A = Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g));
  B = E_li + M - E_g;
  C = E_g*(E_li + M - Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_g)) - M*E_li - Pow2(m);
  
  // Checking the radicand:
  if (Pow2(m)*(Pow2(A) - Pow2(B)) + Pow2(C) < 0.) return 0;
  
  // Checking the denominator:
  if (Abs(Pow2(A) - Pow2(B)) < 1.e-12)
    {
    cout << endl << "Warning: too small denominator!" << endl;
    flag_warn = true;
    return 0;
    }
    
  // The final lepton energy (remember that "en_sign" can be "-1" or "+1"):
  E_lf = (B*C + en_sign*A*Sqrt(Pow2(m)*(Pow2(A) - Pow2(B)) + Pow2(C)))/(Pow2(A) - Pow2(B));
  
  // Checking whether E_lf is a physical root:
  if (Abs(A*Sqrt(Pow2(E_lf) - Pow2(m)) - B*E_lf - C) > 1.e-9) return 0;
  
  // Checking the E_lf value:
  if (E_lf < m || E_lf > E_li - E_g) return 0;
  
  return 1; // If everything is okay
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function SetFinalFourMomenta:
void SetFinalFourMomenta()
  {
  v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
  v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m))*Cos(theta_l), E_lf);
  v_pf = v_li + v_pi - v_lf - v_kf;
  
  // Checking the kinematics:
  if (Abs(Pow2(M) - v_pf*v_pf) > 1.e-8)
    {
    cout << endl << "Warning: bad kinematics! M^2 - v_pf^2 = " << Pow2(M) - v_pf*v_pf << " GeV^2" << endl;
    flag_warn = true;
    }

  // Kinematic parameters of the final proton:
  E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function RosenbluthCS (elastic scattering cross section in the leading order):
double RosenbluthCS()
  {
  // NB: the lepton mass isn't neglected here, see arXiv:1401.2959
  return Pow2(alpha/(2.*E_li))*((1. + qq/(4.*E_li*E_lf))/Pow2(qq/(4.*E_li*E_lf)))*(1./d)*(M*(Pow2(E_lf) - Pow2(m))/(M*E_li*E_lf + m*m*(E_lf - E_li - M)))*(1./(myeps*(1. + tau)))*(myeps*Pow2(G_E(qq)) + tau*Pow2(G_M(qq)));
  }

double f_ros(double theta)
  {
  theta_l = theta; // Theta angle for the lepton

  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  
  return RosenbluthCS()*Sin(theta_l);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Functions f_el_e and f_el_p (for integration of the elastic parts of cross sections):
double f_el_e(double theta) // e-/mu-
  {
  theta_l = theta; // Theta angle for the lepton
  
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  SetFinalFourMomenta(); // Set four-momenta for the final particles
  
  switch (flag_mode)
    {
    case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
      {
      return (1. + d_brem_ee() + inter_virt.Eval(theta_l) + d_prime())*RosenbluthCS()*Sin(theta_l);
      }
      
    case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
      {
      return (1. + d_brem_pp() + inter_virt.Eval(theta_l) + d_prime())*RosenbluthCS()*Sin(theta_l);
      }
      
    case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
      {
      return (1. + d_brem_ee() + d_brem_pp() + d_brem_ep() + inter_virt.Eval(theta_l) + d_prime())*RosenbluthCS()*Sin(theta_l);
      }
    }

  return 0.;
  }

double f_el_p(double theta) // e+/mu+
  {
  theta_l = theta; // Theta angle for the lepton
  
  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  SetFinalFourMomenta(); // Set four-momenta for the final particles
  
  switch (flag_mode)
    {
    case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
      {
      return (1. + d_brem_ee() + inter_virt.Eval(theta_l) - d_prime())*RosenbluthCS()*Sin(theta_l);
      }
      
    case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
      {
      return (1. + d_brem_pp() + inter_virt.Eval(theta_l) - d_prime())*RosenbluthCS()*Sin(theta_l);
      }
      
    case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
      {
      return (1. + d_brem_ee() + d_brem_pp() - d_brem_ep() + inter_virt.Eval(theta_l) - d_prime())*RosenbluthCS()*Sin(theta_l);
      }
    }

  return 0.;
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function scell (evaluation of the z-coordinate of the event):
double scell()
  {
  double r = PseRan->Rndm(); // Random number

  if (flag_cell == 1) return cell_zmin + r*(cell_zmax - cell_zmin); // Uniform distribution
  else // Triangle distribution
    {
    if (r < (cell_zinj - cell_zmin)/(cell_zmax - cell_zmin)) return cell_zmin + Sqrt(r*(cell_zmax - cell_zmin)*(cell_zinj - cell_zmin));
    else return cell_zmax - Sqrt((1. - r)*(cell_zmax - cell_zmin)*(cell_zmax - cell_zinj));
    }
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function EvalAllProducts (evaluation of the four-momentum products):
void EvalAllProducts()
  {
  // Scalar products of the four-momenta:
  kfli = v_kf*v_li; // Photon (v_kf) and initial lepton (v_li)
  kflf = v_kf*v_lf; // Photon (v_kf) and final lepton (v_lf)
  kfpi = v_kf*v_pi; // Photon (v_kf) and initial proton (v_pi)
  kfpf = v_kf*v_pf; // Photon (v_kf) and final proton (v_pf)
  lilf = v_li*v_lf; // Initial lepton (v_li) and final lepton (v_lf)
  lipi = v_li*v_pi; // Initial lepton (v_li) and initial proton (v_pi)
  lipf = v_li*v_pf; // Initial lepton (v_li) and final proton (v_pf)
  pipf = v_pi*v_pf; // Initial proton (v_pi) and final proton (v_pf)
  lfpi = v_lf*v_pi; // Final lepton (v_lf) and initial proton (v_pi)
  lfpf = v_lf*v_pf; // Final lepton (v_lf) and final proton (v_pf)

  // Four-momentum transfers squared:
  q_12 = (v_pf - v_pi)*(v_pf - v_pi);
  q_22 = (v_li - v_lf)*(v_li - v_lf);
  
  // Values for the proton form factors:
  F11 = F1(q_12);
  F21 = F2(q_12);

  F12 = F1(q_22);
  F22 = F2(q_22);

  // In the case of zero four-momentum transfer:
  F10 = F1(0.);
  F20 = F2(0.);
  }
//==============================================================================================

