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

#include "TFile.h"
#include "TRandom3.h" // http://root.cern.ch/root/html/TRandom3.html
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

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

using namespace std;


//----------------------------------------------------------------------------------------------
// Some mathematical and physical constants:
const G4double Pi = TMath::Pi();
const G4double m_e = 0.51099893e-3;          // Mass of the electron/positron (in GeV)
const G4double m_mu = 105.658372e-3;         // Mass of the muon/antimuon (in GeV)
const G4double m_tau = 1.77682;              // Mass of the tau/antitau (in GeV)
const G4double M = 938.272046e-3;            // Mass of the proton (in GeV)
const G4double M2 = TMath::Power(M,2);
const G4double M4 = TMath::Power(M,4);
const G4double mu = 2.79284736;              // Magnetic moment of the proton
const G4double alpha = 1./137.036;           // Fine-structure constant
const G4double e = Sqrt(4.*Pi*alpha); // Electron charge magnitude

const G4double degrad = TMath::Pi()/180.; // Degree to radian conversion

const G4double mb = 0.389379338; // GeV^{-2} to mbarn conversion (GeV^{-2} = 0.389379304 millibarn)
const G4double mkb = 389.379338; // GeV^{-2} to mkbarn conversion


//----------------------------------------------------------------------------------------------
// Kelly parametrization of the proton form factors (J.J. Kelly, PRC 70 (2004) 068202):
const G4double a11_K = -0.24; // Electric form factor
const G4double b11_K = 10.98;
const G4double b12_K = 12.82;
const G4double b13_K = 21.97;

const G4double a21_K = 0.12; // Magnetic form factor
const G4double b21_K = 10.97;
const G4double b22_K = 18.86;
const G4double b23_K = 6.55;


//----------------------------------------------------------------------------------------------
// Puckett parametrization of the proton form factors (A.J.R. Puckett, arXiv:1008.0855):
const G4double a11_P = -0.299; // Electric form factor
const G4double b11_P = 11.11;
const G4double b12_P = 14.11;
const G4double b13_P = 15.7;

const G4double a21_P = 0.081; // Magnetic form factor
const G4double b21_P = 11.15;
const G4double b22_P = 18.45;
const G4double b23_P = 5.31;


//----------------------------------------------------------------------------------------------
// You can use your own parametrization of the proton form factors:
const G4double a11 = 0.; // Electric form factor
const G4double b11 = 9.92;
const G4double b12 = 24.6;
const G4double b13 = 0.;

const G4double a21 = 0.; // Magnetic form factor
const G4double b21 = 9.92;
const G4double b22 = 24.6;
const G4double b23 = 0.;

// The random number generator:
TRandom3 *PseRan = new TRandom3();

G4String event_type;

FILE *fvpol;    // File "vpol.dat" (for calculation of the vacuum polarization)

G4int nCells_1D; // Number of cells for 1D case
G4int nSampl_1D; // Number of samples for 1Dcase
G4int nCells_4D; // Number of cells for 4D case
G4int nSampl_4D; // Number of samples for 4D case

// Flags:
G4int flag_lepton;
G4int flag_mode;
G4int flag_struct;
G4int flag_tpe;
G4int flag_vpol;
G4bool flag_quick;
G4bool flag_target;
G4bool flag_rosen;
G4bool flag_esepp;
G4bool flag_moller;
G4bool flag_onlymoller;

G4double m_l; // The lepton mass (m = m_e OR m = m_mu)
G4double m_l2, m_l4; // Powers of m (for convenience)

G4double E_li; // Full energy of the initial lepton

// Energies of the final particles:
G4double E_lf; // Full energy of the final lepton
G4double E_p;  // Full energy of the final proton
G4double E_g;  // Energy of the photon

// Theta angles ("polar") of the final particles:
G4double theta_l; // Theta angle of the final lepton
G4double theta_p; // Theta angle of the final proton
G4double theta_g; // Theta angle of the photon

// Phi angles ("azimuthal") of the final particles:
G4double phi_l; // Phi angle of the final lepton
G4double phi_p; // Phi angle of the final proton
G4double phi_g; // Phi angle of the photon

G4double theta_min; // Minimum theta angle for the lepton
G4double theta_max; // Maximum theta angle for the lepton

G4double phi_min; // Minimum phi angle for the lepton
G4double phi_max; // Maximum phi angle for the lepton

G4double omega; // Soid angle (steradian)

// Cut and maximum energies for the photon (E_g_cut < E_g < E_g_max):
G4double E_g_cut, E_g_max;

// For the evaluation of the final lepton energy:
G4double A, B, C; // Coefficients of the quadratic equation
G4double en_sign; // Sign to be changed to get two different roots

// Some kinematic variables:
G4double qq;     // Four-momentum transfer
G4double qq_M;     // Four-momentum transfer (Moller)
G4double tau, tt; // tau = -q^2 / (4*M^2)
G4double eps;    // Virtual photon polarization
G4double q_12;   // Four-momentum transfer squared (when the lepton emits the photon)
G4double q_22;   // Four-momentum transfer squared (when the proton emits the photon)
G4double E_1, theta_1, phi_1; // For the lepton #1 (Moller)
G4double E_2, theta_2, phi_2; // For the lepton #2 (Moller)
//G4double E_g, theta_g, phi_g; // For the bremsstrahlung photon (Moller)

G4double Adir, Aex, Aint; // (Moller)
G4double cospsi; // (Moller)

// Some variables to calculate the cross section of elastic scattering:
G4double myeps;
G4double d;

G4double sss, ttt, uuu; // Mandelstam variables (Moller)

// Form factors of the proton:
G4double F10, F20; // F1(0) and F2(0)
G4double F11, F12; // F1(q_12) and F1(q_22)
G4double F21, F22; // F2(q_12) and F2(q_22)

// Storage cell parameters:
G4double cell_xmin, cell_xmax; // x-coordinate of the storage cell
G4double cell_ymin, cell_ymax; // x-coordinate of the storage celll
G4double cell_zmin, cell_zmax; // z-coordinate of the storage cell
G4double cell_zinj; // Z-coordinate of the injection point
TLorentzVector coord;    // Z-coordinate of the event vertex

// Some auxiliary variables:
G4double delta_sum;
G4double M_sum;

G4double xi; // A variable to calculate the vacuum polarization

G4double ret;

G4double res, res1, res2, res3, res4;

// Some four-momenta:
TLorentzVector v_li; // Initial lepton four-momentum
TLorentzVector v_pi; // Initial proton four-momentum
TLorentzVector v_lf; // Final lepton four-momentum
TLorentzVector v_pf; // Final proton four-momentum
TLorentzVector v_kf; // Photon four-momentum
TLorentzVector l1i; // Incident lepton four-momentum (Moller)
TLorentzVector l2i; // Target lepton four-momentum (Moller)
TLorentzVector l1f; // Scattered lepton #1 (Moller)
TLorentzVector l2f; // Scattered lepton #2 (Moller)
TLorentzVector vkf; // Bremsstrahlung photon (Moller)

TLorentzVector p_x; // LorentzVector for numerical integration

// Scalar products and their squares:
G4double kfli; // Photon, lepton initial
G4double kflf; // Photon, lepton final
G4double kfpi; // Photon, proton initial
G4double kfpf; // Photon, proton final
G4double lilf; // Lepton initial, lepton final
G4double lipi; // Lepton initial, proton initial
G4double lipf; // Lepton initial, proton final
G4double pipf; // Proton initial, proton final
G4double lfpi; // Lepton final, proton initial
G4double lfpf; // Lepton final, proton final

// Variables for integration of the bremsstrahlung cross section:
G4double bre1, bre2, bre_error, brp1, brp2, brp_error;
G4double sum_all;

G4double mel, mbr, mbr_error; //(moller)

G4double brems_factor; // To get the soft-photon bremsstrahlung cross section (Moller)

G4long nevents; // Total number of events
G4long ninit;   // Number of events for initialization
G4long loop;    // Counter variable
G4long nev_ElE; // Elastic scattering, e-/mu-
G4long nev_ElP; // Elastic scattering, e+/mu+
  
G4long nev_BrE1; // Bremsstrahlung, e-/mu-, root "-"
G4long nev_BrE2; // Bremsstrahlung, e-/mu-, root "+"

G4long nev_BrP1; // Bremsstrahlung, e+/mu+, root "-"
G4long nev_BrP2; // Bremsstrahlung, e+/mu+, root "+"

G4long nev_ElM; // (Moller)
G4long nev_BrM; // (Moller) Bremsstrahlung
  
G4long nev_0; // Number of Rosenbluth events
G4long nev_e; // Total number of events for e-/mu-
G4long nev_p; // Total number of events for e+/mu+
G4long nev_M; // Total number of Moller events

// Event counters:
G4long count_ElE = 0; // Elastic scattering, e-/mu-
G4long count_ElP = 0; // Elastic scattering, e+/mu+
G4long count_BrE1 = 0; // Bremsstrahlung, e-/mu-, root "-"
G4long count_BrE2 = 0; // Bremsstrahlung, e-/mu-, root "+"
G4long count_BrP1 = 0; // Bremsstrahlung, e+/mu+, root "-"
G4long count_BrP2 = 0; // Bremsstrahlung, e+/mu+, root "+"
G4long count_ElM = 0;    // Moller counter variable
G4long count_BrM = 0;    // Moller bremsstrahlung counter variable

G4double elast_e;
G4double elast_p;
G4double rosen;

// Foam simulators
TFoam *Foam_Ros = new TFoam("Foam_Ros"); // Rosenbluth events
TFoam *Foam_ElE = new TFoam("Foam_ElE"); // Elastic scattering, e-/mu-
TFoam *Foam_BrE1 = new TFoam("Foam_BrE1"); // First-order bremsstrahlung, e-/mu-, en_sign = -1
TFoam *Foam_BrE2 = new TFoam("Foam_BrE2"); // First-order bremsstrahlung, e-/mu-, en_sign = +1
TFoam *Foam_ElP = new TFoam("Foam_ElP"); // Elastic scattering, e+/mu+
TFoam *Foam_BrP1 = new TFoam("Foam_BrP1"); // First-order bremsstrahlung, e+/mu+, en_sign = -1
TFoam *Foam_BrP2 = new TFoam("Foam_BrP2"); // First-order bremsstrahlung, e+/mu+, en_sign = +1

TFoam *Foam1 = new TFoam("Foam1"); // Moller elastic scattering
TFoam *Foam2 = new TFoam("Foam2"); // Moller bremsstrahlung


// Some functions: -----------------------------------------------------------------------------
G4double f_li_lf(G4double x);
G4double f_li_pi(G4double x);
G4double f_li_pf(G4double x);
G4double f_lf_pi(G4double x);
G4double f_lf_pf(G4double x);
G4double f_pi_pf(G4double x);

G4double f_el_e(G4double theta_l);
G4double f_el_p(G4double theta_l);

G4double f_ros(G4double theta_l);

G4double f_l1i_l1f(G4double x); // (Moller)
G4double f_l1i_l2i(G4double x); // (Moller)
G4double f_l1i_l2f(G4double x); // (Moller)
G4double f_l1f_l2i(G4double x); // (Moller)
G4double f_l1f_l2f(G4double x); // (Moller) 
G4double f_l2i_l2f(G4double x); // (Moller)

G4double f_mol(G4double theta_l);

ROOT::Math::Functor1D func_li_lf(&f_li_lf);
ROOT::Math::Functor1D func_li_pi(&f_li_pi);
ROOT::Math::Functor1D func_li_pf(&f_li_pf);
ROOT::Math::Functor1D func_lf_pi(&f_lf_pi);
ROOT::Math::Functor1D func_lf_pf(&f_lf_pf);
ROOT::Math::Functor1D func_pi_pf(&f_pi_pf);

ROOT::Math::Functor1D func_el_e(&f_el_e);
ROOT::Math::Functor1D func_el_p(&f_el_p);

ROOT::Math::Functor1D func_ros(&f_ros);

ROOT::Math::Functor1D func_l1i_l1f(&f_l1i_l1f); // (Moller)
ROOT::Math::Functor1D func_l1i_l2i(&f_l1i_l2i); // (Moller)
ROOT::Math::Functor1D func_l1i_l2f(&f_l1i_l2f); // (Moller)
ROOT::Math::Functor1D func_l1f_l2i(&f_l1f_l2i); // (Moller)
ROOT::Math::Functor1D func_l1f_l2f(&f_l1f_l2f); // (Moller)
ROOT::Math::Functor1D func_l2i_l2f(&f_l2i_l2f); // (Moller)

ROOT::Math::Functor1D func_mol(&f_mol); // (Moller)

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

ROOT::Math::GSLIntegrator i_l1i_l1f(IntOpt); // (Moller)
ROOT::Math::GSLIntegrator i_l1i_l2i(IntOpt); // (Moller)
ROOT::Math::GSLIntegrator i_l1i_l2f(IntOpt); // (Moller)
ROOT::Math::GSLIntegrator i_l1f_l2i(IntOpt); // (Moller)
ROOT::Math::GSLIntegrator i_l1f_l2f(IntOpt); // (Moller)
ROOT::Math::GSLIntegrator i_l2i_l2f(IntOpt); // (Moller)

ROOT::Math::GSLIntegrator i_mol(IntOpt); // (Moller)

// Interpolators: ------------------------------------------------------------------------------
ROOT::Math::Interpolator inter_vpol(10000, InterpolType);
ROOT::Math::Interpolator inter_ros_sin(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_ee(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_ep(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_brem_pp(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_virt(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_prime(InterpolPoints, InterpolType);

ROOT::Math::Interpolator inter_vpol_M(10000, InterpolType); // (Moller)
ROOT::Math::Interpolator inter_brems_M(InterpolPoints, InterpolType); // (Moller)
ROOT::Math::Interpolator inter_virt_M(InterpolPoints, InterpolType);  // (Moller)
ROOT::Math::Interpolator inter_Mol(InterpolPoints, InterpolType); // (Moller)

// Arrays for interpolation:
G4double ss[10000];
G4double rep[10000];
G4double xx[InterpolPoints];
G4double xx_M[InterpolPoints];
G4double y_ros_sin[InterpolPoints];
G4double y_brem_ee[InterpolPoints];
G4double y_brem_ep[InterpolPoints];
G4double y_brem_pp[InterpolPoints];
G4double y_virt[InterpolPoints];
G4double y_prime[InterpolPoints];

G4double y_brems_M[InterpolPoints]; // (Moller)
G4double y_virt_M[InterpolPoints]; // (Moller)
G4double y_Mol[InterpolPoints]; // (Moller)


//----------------------------------------------------------------------------------------------
// Raising to different powers:
inline G4double Pow2(G4double arg) // arg^2
{
  return TMath::Power(arg, 2);
}

inline G4double Pow3(G4double arg) // arg^3
{
  return TMath::Power(arg, 3);
}

inline G4double Pow4(G4double arg) // arg^4
{
  return TMath::Power(arg, 4);
}

inline G4double Pow5(G4double arg) // arg^5
{
  return TMath::Power(arg, 5);
}

inline G4double Pow6(G4double arg) // arg^6
{
  return TMath::Power(arg, 6);
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Form factors of the proton:

// Electric form factor: -----------------------------------------------------------------------
inline G4double G_E(G4double qqq)
{
  tt = Abs(qqq)/(4.*M*M); // Tau

  if (flag_struct == 1) return 1.; // Point-like proton
  
  if (flag_struct == 3) // Kelly parametrization (J.J. Kelly, PRC 70 (2004) 068202)
    {
      return (1. + a11_K*tt)/(1. + b11_K*tt + b12_K*Pow2(tt) + b13_K*Pow3(tt));
    }
    
  if (flag_struct == 4) // Puckett parametrization (A.J.R. Puckett, arXiv:1008.0855)
    {
      return (1. + a11_P*tt)/(1. + b11_P*tt + b12_P*Pow2(tt) + b13_P*Pow3(tt));
    }
    
  if (flag_struct == 5) // Arbitrary parametrization from the file "const.h"
    {
      return (1. + a11*tt)/(1. + b11*tt + b12*Pow2(tt) + b13*Pow3(tt));
    }

  return 1./Pow2(1. + Abs(qqq)/0.71); // Dipole formula
}

// Magnetic form factor: -----------------------------------------------------------------------
inline G4double G_M(G4double qqq)
{
  tt = Abs(qqq)/(4.*M*M); // Tau

  if (flag_struct == 1) return mu; // Point-like proton

  if (flag_struct == 3) // Kelly parametrization (J.J. Kelly, PRC 70 (2004) 068202)
    {
      return mu*(1. + a21_K*tt)/(1. + b21_K*tt + b22_K*Pow2(tt) + b23_K*Pow3(tt));
    }
    
  if (flag_struct == 4) // Puckett parametrization (A.J.R. Puckett, arXiv:1008.0855)
    {
      return mu*(1. + a21_P*tt)/(1. + b21_P*tt + b22_P*Pow2(tt) + b23_P*Pow3(tt));
    }
    
  if (flag_struct == 5) // Arbitrary parametrization from the file "const.h"
    {
      return mu*(1. + a21*tt)/(1. + b21*tt + b22*Pow2(tt) + b23*Pow3(tt));
    }
    
  return mu*G_E(qqq); // Dipole formula
}

// Dirac form factor F1: -----------------------------------------------------------------------
inline G4double F1(G4double qqq)
{
  if (flag_struct == 1) return 1.; // Point-like proton
  
  tt = Abs(qqq)/(4.*M*M); // Tau
  return (G_E(qqq) + tt*G_M(qqq))/(1. + tt);
}

// Pauli form factor F2: -----------------------------------------------------------------------
inline G4double F2(G4double qqq)
{
  if (flag_struct == 1) return mu - 1.; // Point-like proton
  
  tt = Abs(qqq)/(4.*M*M); // Tau
  return (G_M(qqq) - G_E(qqq))/(1. + tt);
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function EvalKinematicParams (evaluation of some kinematic parameters):
void EvalKinematicParams(G4double E_gamma)
{
  qq = 2.*M*(E_lf - E_li + E_gamma); // Four-momentum transfer squared
  tau = -qq/(4.*M*M); // Tau
  eps = 1./(1. + 2.*(1. + tau)*Pow2(Tan(theta_l/2.))); // Epsilon
  // To calculate the Rosenbluth cross section without neglecting the lepton mass:
  myeps = 1./(1. - 2.*(1. + tau)*(qq + 2.*m_l*m_l)/(4.*E_li*E_lf + qq)); // Modified epsilon
  d = (E_lf/E_li)*Sqrt((Pow2(E_li) - Pow2(m_l))/(Pow2(E_lf) - Pow2(m_l)));
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Soft-photon bremsstrahlung cross section:
G4double BremSoftLeptonOnly() // Lepton term
{
  return -(alpha*v_kf.E()/(4.*Pi*Pi))*(Pow2(m_l)/Pow2(v_kf*v_li) - 2.*(v_li*v_lf)/((v_kf*v_li)*(v_kf*v_lf)) + Pow2(m_l)/Pow2(v_kf*v_lf));
}

G4double BremSoftProtonOnly() // Proton term
{
  return -(alpha*v_kf.E()/(4.*Pi*Pi))*(Pow2(M)/Pow2(v_kf*v_pi) - 2.*(v_pi*v_pf)/((v_kf*v_pi)*(v_kf*v_pf)) + Pow2(M)/Pow2(v_kf*v_pf));
}
  
G4double BremSoftInterference() // Interference term
{
  return -(alpha*v_kf.E()/(2.*Pi*Pi))*(-(v_li*v_pi)/((v_kf*v_li)*(v_kf*v_pi)) + (v_li*v_pf)/((v_kf*v_li)*(v_kf*v_pf)) + (v_lf*v_pi)/((v_kf*v_lf)*(v_kf*v_pi)) - (v_lf*v_pf)/((v_kf*v_lf)*(v_kf*v_pf)));
}
//==============================================================================================
 

//----------------------------------------------------------------------------------------------
// Functions for numerical integration (soft-photon bremsstrahlung):
G4double f_li_lf(G4double x) // For the calculation of B(v_li, v_lf, E_g_cut)
{
  p_x = x*v_li + (1. - x)*v_lf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_li_pi(G4double x) // For the calculation of B(v_li, v_pi, E_g_cut)
{
  p_x = x*v_li + (1. - x)*v_pi;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_li_pf(G4double x) // For the calculation of B(v_li, v_pf, E_g_cut)
{
  p_x = x*v_li + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_lf_pi(G4double x) // For the calculation of B(v_lf, v_pi, E_g_cut)
{
  p_x = x*v_lf + (1. - x)*v_pi;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_lf_pf(G4double x) // For the calculation of B(v_lf, v_pf, E_g_cut)
{
  p_x = x*v_lf + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_pi_pf(G4double x) // For the calculation of B(v_pi, v_pf, E_g_cut)
{
  p_x = x*v_pi + (1. - x)*v_pf;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Functions for numerical integration (soft-photon bremsstrahlung):
G4double f_l1i_l1f(G4double x) // For the calculation of B(l1i, l1f, E_g_cut)
{
  p_x = x*l1i + (1. - x)*l1f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_l1i_l2i(G4double x) // For the calculation of B(l1i, l2i, E_g_cut)
{
  p_x = x*l1i + (1. - x)*l2i;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_l1i_l2f(G4double x) // For the calculation of B(l1i, l2f, E_g_cut)
{
  p_x = x*l1i + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_l1f_l2i(G4double x) // For the calculation of B(l1f, l2i, E_g_cut)
{
  p_x = x*l1f + (1. - x)*l2i;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_l1f_l2f(G4double x) // For the calculation of B(l1f, l2f, E_g_cut)
{
  p_x = x*l1f + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}

G4double f_l2i_l2f(G4double x) // For the calculation of B(l2i, l2f, E_g_cut)
{
  p_x = x*l2i + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Integration of the soft-photon bremsstrahlung cross section: (Moller)
inline G4double d_li_li() // B(v_li, v_li, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/m_l) + E_li*Log(m_l/(E_li + Sqrt(Pow2(E_li) - Pow2(m_l))))/Sqrt(Pow2(E_li) - Pow2(m_l)))/Pi;
}

inline G4double d_li_lf() // B(v_li, v_lf, E_g_cut)
{
  return (v_li*v_lf)*i_li_lf.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_lf_lf() // B(v_lf, v_lf, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/m_l) + E_lf*Log(m_l/(E_lf + Sqrt(Pow2(E_lf) - Pow2(m_l))))/Sqrt(Pow2(E_lf) - Pow2(m_l)))/Pi;
}

inline G4double d_li_pi() // B(v_li, v_pi, E_g_cut)
{
  return (v_li*v_pi)*i_li_pi.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_li_pf() // B(v_li, v_pf, E_g_cut)
{
  return (v_li*v_pf)*i_li_pf.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_lf_pi() // B(v_lf, v_pi, E_g_cut)
{
  return (v_lf*v_pi)*i_lf_pi.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_lf_pf() // B(v_lf, v_pf, E_g_cut)
{
  return (v_lf*v_pf)*i_lf_pf.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_pi_pi() // B(v_pi, v_pi, E_g_cut)
{
  return (Log(2.*E_g_cut/M) - 1.)/(2.*Pi);
}

inline G4double d_pi_pf() // B(v_pi, v_pf, E_g_cut)
{
  return (v_pi*v_pf)*i_pi_pf.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_pf_pf() // B(v_pf, v_pf, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/M) + E_p*Log(M/(E_p + Sqrt(Pow2(E_p) - Pow2(M))))/Sqrt(Pow2(E_p) - Pow2(M)))/Pi;
}

// Lepton term:
inline G4double d_brem_ee()
{
  return -2.*alpha*(d_li_li() - 2.*d_li_lf() + d_lf_lf());
}

// Interference term:
inline G4double d_brem_ep()
{
  return 4.*alpha*(d_li_pi() - d_li_pf() - d_lf_pi() + d_lf_pf());
}

// Proton term:
inline G4double d_brem_pp()
{
  return -2.*alpha*(d_pi_pi() - 2.*d_pi_pf() + d_pf_pf());
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Integration of the soft-photon bremsstrahlung cross section: (Moller)
inline G4double d_l1i_l1i() // B(l1i, l1i, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/m_l) + E_li*Log(m_l/(E_li + Sqrt(Pow2(E_li) - Pow2(m_l))))/Sqrt(Pow2(E_li) - Pow2(m_l)))/Pi;
}

inline G4double d_l1i_l1f() // B(l1i, l1f, E_g_cut)
{
  return (l1i*l1f)*i_l1i_l1f.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l1f_l1f() // B(l1f, l1f, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/m_l) + E_1*Log(m_l/(E_1 + Sqrt(Pow2(E_1) - Pow2(m_l))))/Sqrt(Pow2(E_1) - Pow2(m_l)))/Pi;
}

inline G4double d_l1i_l2i() // B(l1i, l2i, E_g_cut)
{
  return (l1i*l2i)*i_l1i_l2i.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l1i_l2f() // B(l1i, l2f, E_g_cut)
{
  return (l1i*l2f)*i_l1i_l2f.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l1f_l2i() // B(l1f, l2i, E_g_cut)
{
  return (l1f*l2i)*i_l1f_l2i.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l1f_l2f() // B(l1f, l2f, E_g_cut)
{
  return (l1f*l2f)*i_l1f_l2f.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l2i_l2i() // B(l2i, l2i, E_g_cut)
{
  return (Log(2.*E_g_cut/m_l) - 1.)/(2.*Pi);
}

inline G4double d_l2i_l2f() // B(l2i, l2f, E_g_cut)
{
  return (l2i*l2f)*i_l2i_l2f.Integral(0., 1.)/(4.*Pi);
}

inline G4double d_l2f_l2f() // B(l2f, l2f, E_g_cut)
{
  return 0.5*(Log(2.*E_g_cut/m_l) + E_2*Log(m_l/(E_2 + Sqrt(Pow2(E_2) - Pow2(m_l))))/Sqrt(Pow2(E_2) - Pow2(m_l)))/Pi;
}

// Total bremsstrahlung correction, Moller scattering:
inline G4double d_brems_M()
{
  return -2.*alpha*(d_l1i_l1i() - 2.*d_l1i_l1f() + d_l1f_l1f() + d_l2i_l2i() - 2.*d_l2i_l2f() + d_l2f_l2f() + 2.*d_l1i_l2i() - 2.*d_l1i_l2f() - 2.*d_l1f_l2i() + 2.*d_l1f_l2f());
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Virtual-photon corrections:

// Vacuum polarization contribution from electron-positron loops:
G4double d_vac_e()
{
  return (2.*alpha/(3.*Pi))*(-5./3. + Log(-qq/Pow2(m_e)));
}

// Vacuum polarization contribution from muon loops:
G4double d_vac_mu()
{
  xi = -(4.*Pow2(m_mu)/qq)/Pow2(1. + Sqrt(1. - 4.*Pow2(m_mu)/qq));
  return (2*alpha/(3.*Pi))*(-5./3. + 4.*xi/Pow2(1. - xi) - ((1. - 4.*xi + Pow2(xi))/Pow2(1. - xi))*((1. + xi)/(1. - xi))*Log(xi));
}
  
// Vacuum polarization contribution from tau loops:
G4double d_vac_tau()
{
  xi = -(4.*Pow2(m_tau)/qq)/Pow2(1. + Sqrt(1. - 4.*Pow2(m_tau)/qq));
  return (2*alpha/(3.*Pi))*(-5./3. + 4.*xi/Pow2(1. - xi) - ((1. - 4.*xi + Pow2(xi))/Pow2(1. - xi))*((1. + xi)/(1. - xi))*Log(xi));
}
  
// Electron vertex correction:
G4double d_vertex()
{
  return (alpha/Pi)*(3.*Log(-qq/Pow2(m_l))/2. - 2.);
}

// Vertex correction:
inline double d_vertex_M()
{
  return (2.*alpha/Pi)*(3.*Log(-qq_M/Pow2(m_l))/2. - 2.); //NB: twice larger than for ep scattering
}

// Two-photon exchange contribution:
G4double d_prime()
{
  if (flag_tpe == 1) return 0.; // Approach of Mo & Tsai to the TPE diagrams
    
  // Approach of Maximon & Tjon to the TPE diagrams (see arXiv:1401.2959):    
  return -(alpha/Pi)*(Log(E_li/E_lf)*Log(qq*qq/(4.*M*M*E_li*E_lf)) + 2.*DiLog(1. - 0.5*M/E_li) - 2.*DiLog(1. - 0.5*M/E_lf));
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Function ElasticEnergy (final lepton energy in the case of purely elastic scattering):
G4double ElasticEnergy(G4double theta)
{
  return ((E_li + M)*(M*E_li + Pow2(m_l)) + Sqrt(Pow2(M) - Pow2(m_l*Sin(theta)))*(Pow2(E_li) - Pow2(m_l))*Cos(theta))/(Pow2(E_li + M) - (Pow2(E_li) - Pow2(m_l))*Pow2(Cos(theta)));
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function EvalEnergy (exact evaluation of the final lepton energy):
G4int EvalEnergy()
{
  // Calculation of the coefficients A, B and C:
  A = Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g));
  B = E_li + M - E_g;
  C = E_g*(E_li + M - Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_g)) - M*E_li - Pow2(m_l);
  
  // Checking the radicand:
  if (Pow2(m_l)*(Pow2(A) - Pow2(B)) + Pow2(C) < 0.) return 0;
  
  // Checking the denominator:
  if (Abs(Pow2(A) - Pow2(B)) < 1.e-12)
    {
      cout << endl << "Warning: too small denominator!" << endl;
      return 0;
    }
    
  // The final lepton energy (remember that "en_sign" can be "-1" or "+1"):
  E_lf = (B*C + en_sign*A*Sqrt(Pow2(m_l)*(Pow2(A) - Pow2(B)) + Pow2(C)))/(Pow2(A) - Pow2(B));
  
  // Checking whether E_lf is a physical root:
  if (Abs(A*Sqrt(Pow2(E_lf) - Pow2(m_l)) - B*E_lf - C) > 1.e-9) return 0;
  
  // Checking the E_lf value:
  if (E_lf < m_l || E_lf > E_li - E_g) return 0;
  
  return 1; // If everything is okay
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function SetFinalFourMomenta:
void SetFinalFourMomenta()
{
  v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
  v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m_l))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m_l))*Cos(theta_l), E_lf);
  v_pf = v_li + v_pi - v_lf - v_kf;
  
  // Checking the kinematics:
  if (Abs(Pow2(M) - v_pf*v_pf) > 1.e-8)
    {
      cout << endl << "Warning: bad kinematics! M^2 - v_pf^2 = " << Pow2(M) - v_pf*v_pf << " GeV^2" << endl;
    }

  // Kinematic parameters of the final proton:
  E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
}
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function RosenbluthCS (elastic scattering cross section in the leading order):
G4double RosenbluthCS()
{
  // NB: the lepton mass isn't neglected here, see arXiv:1401.2959
  return Pow2(alpha/(2.*E_li))*((1. + qq/(4.*E_li*E_lf))/Pow2(qq/(4.*E_li*E_lf)))*(1./d)*(M*(Pow2(E_lf) - Pow2(m_l))/(M*E_li*E_lf + m_l*m_l*(E_lf - E_li - M)))*(1./(myeps*(1. + tau)))*(myeps*Pow2(G_E(qq)) + tau*Pow2(G_M(qq)));
}

G4double f_ros(G4double theta)
{
  theta_l = theta; // Theta angle for the lepton

  E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

  EvalKinematicParams(0.); // Evaluation of some kinematic parameters
  
  return RosenbluthCS()*Sin(theta_l);
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Function MCS (Moller cross section in the lowest order): (Moller)
G4double MCS()
{
  // Energy of the scattered lepton l1'
  E_1 = m_l*(E_li + m_l + (E_li - m_l)*Pow2(Cos(theta_1)))/(E_li + m_l - (E_li - m_l)*Pow2(Cos(theta_1)));

  // Four-momenta l1' and l2':
  l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m_l))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m_l))*Cos(theta_1), E_1);
  l2f = l1i + l2i - l1f;

  // Full energy of the lepton l2':
  E_2 = l2f.E();

  // Mandelstam variables t and u:
  ttt = (l1f - l1i)*(l1f - l1i);
  uuu = (l1f - l2i)*(l1f - l2i);

  // Adir, Aex and Aint for Moller scattering:
  Adir = Pow2(sss - 2.*m_l*m_l) + Pow2(uuu - 2.*m_l*m_l) + 4.*m_l*m_l*ttt;
  Aex = Pow2(sss - 2.*m_l*m_l) + Pow2(ttt - 2.*m_l*m_l) + 4.*m_l*m_l*uuu;
  Aint = -(sss - 2.*m_l*m_l)*(sss - 6.*m_l*m_l);

  return 2.*Pow2(alpha)*(Cos(theta_1)/Pow2(E_li + m_l - (E_li - m_l)*Pow2(Cos(theta_1))))*(Adir/Pow2(ttt) + Aex/Pow2(uuu) - 2.*Aint/(ttt*uuu));
}

G4double f_mol (G4double theta)
{
  theta_1 = theta;
  return MCS()*(1. + inter_virt_M.Eval(theta_1) + inter_brems_M.Eval(theta_1))*Sin(theta_1);
}
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Functions f_el_e and f_el_p (for integration of the elastic parts of cross sections):
G4double f_el_e(G4double theta) // e-/mu-
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

G4double f_el_p(G4double theta) // e+/mu+
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
TLorentzVector scell()
{
  G4double rx = PseRan->Rndm(); // Random number
  G4double ry = PseRan->Rndm(); // Random number
  G4double rz = PseRan->Rndm(); // Random number
  return TLorentzVector(cell_xmin + rx*(cell_xmax - cell_xmin),cell_ymin + ry*(cell_ymax - cell_ymin),cell_zmin + rz*(cell_zmax - cell_zmin),0.);
  
/*    if (flag_cell == 1) return cell_zmin + r*(cell_zmax - cell_zmin); // Uniform distribution
  else // Triangle distribution
    {
      if (r < (cell_zinj - cell_zmin)/(cell_zmax - cell_zmin)) return cell_zmin + Sqrt(r*(cell_zmax - cell_zmin)*(cell_zinj - cell_zmin));
      else return cell_zmax - Sqrt((1. - r)*(cell_zmax - cell_zmin)*(cell_zmax - cell_zinj));
      }*/
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



//----------------------------------------------------------------------------------------------
// First-order bremsstrahlung, proton term:
G4double pterm ()
{
  res = -((pow(F22,2)*pow(kfpi,4)*(-lilf + m_l2)*M2*
	   ((-2*pow(F10,2) + pow(F20,2))*kfpf + 2*pow(F10,2)*M2) + 2*pow(F10,2)*pow(kfpf,2)*M4*
	   (pow(F22,2)*pow(kfpf,2)*(-lilf + m_l2) - 
	    4*pow(F12,2)*kflf*lipf*M2 - 
	    4*F12*F22*kflf*lipf*M2 - 
	    pow(F22,2)*kflf*lipf*M2 - 
	    4*F12*F22*lfpf*lipf*M2 - 
	    3*pow(F22,2)*lfpf*lipf*M2 + 
	    4*pow(F12,2)*lfpi*lipf*M2 + 
	    4*F12*F22*lfpi*lipf*M2 + 
	    pow(F22,2)*lfpi*lipf*M2 + 
	    4*F12*F22*kflf*lipi*M2 + 
	    3*pow(F22,2)*kflf*lipi*M2 + 
	    4*pow(F12,2)*lfpf*lipi*M2 + 
	    4*F12*F22*lfpf*lipi*M2 + 
	    pow(F22,2)*lfpf*lipi*M2 - 
	    4*F12*F22*lfpi*lipi*M2 - 
	    3*pow(F22,2)*lfpi*lipi*M2 - 
	    4*pow(F12,2)*lilf*M4 - 4*F12*F22*lilf*M4 - 
	    pow(F22,2)*lilf*M4 + 
	    8*pow(F12,2)*m_l2*M4 + 
	    12*F12*F22*m_l2*M4 + 
	    5*pow(F22,2)*m_l2*M4 - 
	    pow(F22,2)*kflf*lipf*pipf + pow(F22,2)*lfpf*lipf*pipf + 
	    pow(F22,2)*lfpi*lipf*pipf - pow(F22,2)*kflf*lipi*pipf + 
	    pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf + 
	    4*F12*F22*lilf*M2*pipf + 
	    2*pow(F22,2)*lilf*M2*pipf - 
	    4*pow(F12,2)*m_l2*M2*pipf - 
	    12*F12*F22*m_l2*M2*pipf - 
	    6*pow(F22,2)*m_l2*M2*pipf - 
	    pow(F22,2)*lilf*pow(pipf,2) + 
	    pow(F22,2)*m_l2*pow(pipf,2) + 
	    kfpf*(-(pow(F22,2)*lfpf*lipf) - pow(F22,2)*lfpi*lipf - 
		  pow(F22,2)*lfpf*lipi - pow(F22,2)*lfpi*lipi + 
		  pow(F22,2)*kflf*(lipf + lipi) + 
		  4*pow(F12,2)*m_l2*M2 + 
		  12*F12*F22*m_l2*M2 + 
		  6*pow(F22,2)*m_l2*M2 - 
		  2*pow(F22,2)*m_l2*pipf - 
		  2*F22*lilf*((2*F12 + F22)*M2 - F22*pipf)) - 
	    kfli*(F22*kflf*(F22*kfpf + (4*F12 + 3*F22)*M2 - 
			    F22*pipf) - F22*lfpi*(F22*kfpf + (4*F12 + 3*F22)*M2 - F22*pipf) + 
		  lfpf*(-(pow(F22,2)*kfpf) + pow(2*F12 + F22,2)*M2 + 
			pow(F22,2)*pipf))) - 
	   pow(kfpi,3)*(-4*pow(F20,2)*pow(F22,2)*pow(kfpf,3)*
			m_l2 + F22*kfli*(-(F22*kflf*M2*((-2*pow(F10,2) + pow(F20,2))*kfpf + 
							2*pow(F10,2)*M2)) - 
					 F22*lfpi*M2*((-2*pow(F10,2) + pow(F20,2))*kfpf + 
						      2*pow(F10,2)*M2) + 
					 lfpf*(2*pow(F20,2)*F22*pow(kfpf,2) + 
					       (4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
						pow(F20,2)*(2*F12 + F22))*kfpf*M2 - 
					       2*pow(F10,2)*F22*M4)) + 
			F22*pow(kfpf,2)*(2*pow(F20,2)*F22*kflf*lipf + 
					 4*pow(F20,2)*F22*lfpf*lipf + 2*pow(F20,2)*F22*lfpi*lipf + 
					 2*pow(F20,2)*F22*lfpf*lipi - 
					 8*F10*F12*F20*m_l2*M2 - 
					 8*F12*pow(F20,2)*m_l2*M2 + 
					 2*pow(F10,2)*F22*m_l2*M2 - 
					 8*F10*F20*F22*m_l2*M2 - 
					 7*pow(F20,2)*F22*m_l2*M2 + 
					 2*pow(F20,2)*F22*m_l2*pipf - 
					 F22*lilf*((2*pow(F10,2) + pow(F20,2))*M2 + 
						   2*pow(F20,2)*pipf)) + 
			2*pow(F10,2)*M4*(-(pow(F22,2)*lfpf*lipf) - pow(F22,2)*lfpi*lipf - 
					 pow(F22,2)*lfpf*lipi - pow(F22,2)*lfpi*lipi - 
					 pow(F22,2)*kflf*(lipf + lipi) + 
					 4*pow(F12,2)*m_l2*M2 + 
					 12*F12*F22*m_l2*M2 + 
					 6*pow(F22,2)*m_l2*M2 - 
					 2*pow(F22,2)*m_l2*pipf - 
					 2*F22*lilf*((2*F12 + F22)*M2 - F22*pipf)) + 
			kfpf*M2*(-8*pow(F12,2)*pow(F20,2)*lfpf*lipf - 
				 12*F10*F12*F20*F22*lfpf*lipf - 
				 12*F12*pow(F20,2)*F22*lfpf*lipf + 
				 6*pow(F10,2)*pow(F22,2)*lfpf*lipf - 
				 12*F10*F20*pow(F22,2)*lfpf*lipf - 
				 5*pow(F20,2)*pow(F22,2)*lfpf*lipf + 
				 2*F10*F12*F20*F22*lfpi*lipf + 
				 2*F12*pow(F20,2)*F22*lfpi*lipf + 
				 4*pow(F10,2)*pow(F22,2)*lfpi*lipf + 
				 2*F10*F20*pow(F22,2)*lfpi*lipf + 
				 pow(F20,2)*pow(F22,2)*lfpi*lipf + 
				 2*F10*F12*F20*F22*lfpf*lipi + 
				 2*F12*pow(F20,2)*F22*lfpf*lipi + 
				 4*pow(F10,2)*pow(F22,2)*lfpf*lipi + 
				 2*F10*F20*pow(F22,2)*lfpf*lipi + 
				 pow(F20,2)*pow(F22,2)*lfpf*lipi + 
				 2*pow(F10,2)*pow(F22,2)*lfpi*lipi - 
				 pow(F20,2)*pow(F22,2)*lfpi*lipi + 
				 F22*kflf*((4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
					    pow(F20,2)*(2*F12 + F22))*lipf + 
					   (2*pow(F10,2) - pow(F20,2))*F22*lipi) - 
				 8*pow(F10,2)*pow(F12,2)*m_l2*M2 - 
				 8*F10*pow(F12,2)*F20*m_l2*M2 - 
				 8*pow(F10,2)*F12*F22*m_l2*M2 - 
				 4*F10*F12*F20*F22*m_l2*M2 + 
				 4*F12*pow(F20,2)*F22*m_l2*M2 + 
				 2*pow(F10,2)*pow(F22,2)*m_l2*M2 + 
				 4*F10*F20*pow(F22,2)*m_l2*M2 + 
				 3*pow(F20,2)*pow(F22,2)*m_l2*M2 + 
				 6*pow(F10,2)*pow(F22,2)*m_l2*pipf - 
				 pow(F20,2)*pow(F22,2)*m_l2*pipf + 
				 lilf*((4*pow(F12,2)*pow(F20,2) + 
					4*F12*F20*(F10 + F20)*F22 + 
					(-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*pow(F22,2))*M2 + 
				       (-6*pow(F10,2) + pow(F20,2))*pow(F22,2)*pipf))) + 
	   kfpf*kfpi*M2*((2*pow(F10,2) - pow(F20,2))*pow(F22,2)*
			 pow(kfpf,3)*(-lilf + m_l2) - 
			 4*pow(F10,2)*pipf*(-4*F12*F22*lfpf*lipf*M2 - 
					    3*pow(F22,2)*lfpf*lipf*M2 + 
					    4*pow(F12,2)*lfpi*lipf*M2 + 
					    4*F12*F22*lfpi*lipf*M2 + 
					    pow(F22,2)*lfpi*lipf*M2 + 
					    4*pow(F12,2)*lfpf*lipi*M2 + 
					    4*F12*F22*lfpf*lipi*M2 + 
					    pow(F22,2)*lfpf*lipi*M2 - 
					    4*F12*F22*lfpi*lipi*M2 - 
					    3*pow(F22,2)*lfpi*lipi*M2 + 
					    2*pow(F12 + F22,2)*kflf*(-lipf + lipi)*M2 + 
					    8*pow(F12,2)*m_l2*M4 + 
					    12*F12*F22*m_l2*M4 + 
					    5*pow(F22,2)*m_l2*M4 + 
					    pow(F22,2)*lfpf*lipf*pipf + pow(F22,2)*lfpi*lipf*pipf + 
					    pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf - 
					    4*pow(F12,2)*m_l2*M2*pipf - 
					    12*F12*F22*m_l2*M2*pipf - 
					    6*pow(F22,2)*m_l2*M2*pipf + 
					    pow(F22,2)*m_l2*pow(pipf,2) - 
					    lilf*pow((2*F12 + F22)*M2 - F22*pipf,2)) - 
			 pow(kfpf,2)*(2*pow(F10,2)*pow(F22,2)*lfpf*lipf - 
				      pow(F20,2)*pow(F22,2)*lfpf*lipf + 
				      2*F10*F12*F20*F22*lfpi*lipf + 
				      2*F12*pow(F20,2)*F22*lfpi*lipf + 
				      4*pow(F10,2)*pow(F22,2)*lfpi*lipf + 
				      2*F10*F20*pow(F22,2)*lfpi*lipf + 
				      pow(F20,2)*pow(F22,2)*lfpi*lipf + 
				      2*F10*F12*F20*F22*lfpf*lipi + 
				      2*F12*pow(F20,2)*F22*lfpf*lipi + 
				      4*pow(F10,2)*pow(F22,2)*lfpf*lipi + 
				      2*F10*F20*pow(F22,2)*lfpf*lipi + 
				      pow(F20,2)*pow(F22,2)*lfpf*lipi - 
				      8*pow(F12,2)*pow(F20,2)*lfpi*lipi - 
				      12*F10*F12*F20*F22*lfpi*lipi - 
				      12*F12*pow(F20,2)*F22*lfpi*lipi + 
				      6*pow(F10,2)*pow(F22,2)*lfpi*lipi - 
				      12*F10*F20*pow(F22,2)*lfpi*lipi - 
				      5*pow(F20,2)*pow(F22,2)*lfpi*lipi - 
				      F22*kflf*((2*pow(F10,2) - pow(F20,2))*F22*lipf + 
						(4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
						 pow(F20,2)*(2*F12 + F22))*lipi) - 
				      8*pow(F10,2)*pow(F12,2)*m_l2*M2 - 
				      8*F10*pow(F12,2)*F20*m_l2*M2 - 
				      8*pow(F10,2)*F12*F22*m_l2*M2 - 
				      4*F10*F12*F20*F22*m_l2*M2 + 
				      4*F12*pow(F20,2)*F22*m_l2*M2 + 
				      2*pow(F10,2)*pow(F22,2)*m_l2*M2 + 
				      4*F10*F20*pow(F22,2)*m_l2*M2 + 
				      3*pow(F20,2)*pow(F22,2)*m_l2*M2 + 
				      6*pow(F10,2)*pow(F22,2)*m_l2*pipf - 
				      pow(F20,2)*pow(F22,2)*m_l2*pipf + 
				      lilf*((4*pow(F12,2)*pow(F20,2) + 
					     4*F12*F20*(F10 + F20)*F22 + 
					     (-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*
					     pow(F22,2))*M2 + 
					    (-6*pow(F10,2) + pow(F20,2))*pow(F22,2)*pipf)) + 
			 kfli*(F22*pow(kfpf,2)*
			       ((2*pow(F10,2) - pow(F20,2))*F22*lfpf + 
				(4*pow(F10,2)*F22 + 2*F10*F20*(F12 + F22) + 
				 pow(F20,2)*(2*F12 + F22))*lfpi) - 
			       8*pow(F10,2)*pow(F12 + F22,2)*(-lfpf + lfpi)*M2*
			       pipf + 2*kfpf*(lfpi*((2*pow(F10,2)*F22*(2*F12 + F22) + 
						     pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
						     F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
						    (4*pow(F12,2)*pow(F20,2) + 
						     7*F12*F20*(F10 + F20)*F22 + 
						     (3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*
						     pow(F22,2))*pipf) + 
					      lfpf*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
								   pow(F20,2)*(2*F12 + F22))*M2) + 
						    F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
							 pow(F20,2)*(F12 + F22))*pipf)) + 
			       kflf*((-2*pow(F10,2) + pow(F20,2))*pow(F22,2)*pow(kfpf,2) - 
				     2*F22*kfpf*((pow(F10,2)*F22 - 2*F10*F20*(F12 + F22) - 
						  2*pow(F20,2)*(F12 + F22))*M2 + 
						 2*F20*(F10 + F20)*(F12 + F22)*pipf) + 
				     4*(4*F10*F12*(F10 + F20)*(F12 + F22)*M4 + 
					(pow(F10,2)*pow(F22,2) - 
					 4*F10*F20*pow(F12 + F22,2) - 
					 2*pow(F20,2)*pow(F12 + F22,2))*M2*pipf + 
					(2*pow(F12,2)*pow(F20,2) + 
					 4*F12*F20*(F10 + F20)*F22 + 
					 (pow(F10,2) + 4*F10*F20 + 2*pow(F20,2))*
					 pow(F22,2))*pow(pipf,2)))) - 
			 2*kfpf*(kflf*(-(lipi*((2*pow(F10,2)*F22*(2*F12 + F22) + 
						pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
						F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
					       (4*pow(F12,2)*pow(F20,2) + 
						7*F12*F20*(F10 + F20)*F22 + 
						(3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*
						pow(F22,2))*pipf)) + 
				       lipf*((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
							  pow(F20,2)*(2*F12 + F22))*M2 + 
					     F22*(pow(F10,2)*F22 - F10*F20*(F12 + F22) - 
						  pow(F20,2)*(F12 + F22))*pipf)) + 
				 pow(F10,2)*(-(pow(F22,2)*lfpf*lipf*M2) - 
					     4*pow(F12,2)*lfpi*lipf*M2 + 
					     pow(F22,2)*lfpi*lipf*M2 + 
					     4*pow(F12,2)*m_l2*M4 + 
					     12*F12*F22*m_l2*M4 + 
					     6*pow(F22,2)*m_l2*M4 - 
					     2*pow(F22,2)*lfpf*lipf*pipf - 
					     4*pow(F22,2)*lfpi*lipf*pipf + 
					     8*pow(F12,2)*m_l2*M2*pipf + 
					     24*F12*F22*m_l2*M2*pipf + 
					     10*pow(F22,2)*m_l2*M2*pipf - 
					     4*pow(F22,2)*m_l2*pow(pipf,2) - 
					     2*F22*lilf*(M2 + 2*pipf)*
					     ((2*F12 + F22)*M2 - F22*pipf) + 
					     lipi*(lfpi*((-8*pow(F12,2) + 3*pow(F22,2))*
							 M2 - 6*pow(F22,2)*pipf) + 
						   lfpf*((-4*pow(F12,2) + pow(F22,2))*M2 - 
							 4*pow(F22,2)*pipf))))) + 
	   pow(kfpi,2)*(F22*pow(kfpf,3)*
			(2*pow(F20,2)*F22*lfpi*lipf - 2*pow(F20,2)*F22*kflf*lipi + 
			 2*pow(F20,2)*F22*lfpf*lipi + 4*pow(F20,2)*F22*lfpi*lipi - 
			 8*F10*F12*F20*m_l2*M2 - 
			 8*F12*pow(F20,2)*m_l2*M2 + 
			 2*pow(F10,2)*F22*m_l2*M2 - 
			 8*F10*F20*F22*m_l2*M2 - 
			 7*pow(F20,2)*F22*m_l2*M2 + 
			 2*pow(F20,2)*F22*m_l2*pipf - 
			 F22*lilf*((2*pow(F10,2) + pow(F20,2))*M2 + 
				   2*pow(F20,2)*pipf)) - 
			2*pow(kfpf,2)*(-2*F10*F12*F20*F22*lfpf*lipf*M2 - 
				       2*F12*pow(F20,2)*F22*lfpf*lipf*M2 - 
				       4*pow(F10,2)*pow(F22,2)*lfpf*lipf*M2 - 
				       2*F10*F20*pow(F22,2)*lfpf*lipf*M2 - 
				       2*pow(F20,2)*pow(F22,2)*lfpf*lipf*M2 + 
				       4*pow(F12,2)*pow(F20,2)*lfpi*lipf*M2 + 
				       6*F10*F12*F20*F22*lfpi*lipf*M2 + 
				       6*F12*pow(F20,2)*F22*lfpi*lipf*M2 - 
				       4*pow(F10,2)*pow(F22,2)*lfpi*lipf*M2 + 
				       6*F10*F20*pow(F22,2)*lfpi*lipf*M2 + 
				       2*pow(F20,2)*pow(F22,2)*lfpi*lipf*M2 + 
				       4*pow(F12,2)*pow(F20,2)*lfpf*lipi*M2 + 
				       6*F10*F12*F20*F22*lfpf*lipi*M2 + 
				       6*F12*pow(F20,2)*F22*lfpf*lipi*M2 - 
				       4*pow(F10,2)*pow(F22,2)*lfpf*lipi*M2 + 
				       6*F10*F20*pow(F22,2)*lfpf*lipi*M2 + 
				       2*pow(F20,2)*pow(F22,2)*lfpf*lipi*M2 - 
				       2*F10*F12*F20*F22*lfpi*lipi*M2 - 
				       2*F12*pow(F20,2)*F22*lfpi*lipi*M2 - 
				       4*pow(F10,2)*pow(F22,2)*lfpi*lipi*M2 - 
				       2*F10*F20*pow(F22,2)*lfpi*lipi*M2 - 
				       2*pow(F20,2)*pow(F22,2)*lfpi*lipi*M2 - 
				       F22*(pow(F10,2)*F22 - 3*F10*F20*(F12 + F22) - 
					    3*pow(F20,2)*(F12 + F22))*kflf*(-lipf + lipi)*M2 +
				       8*F10*pow(F12,2)*F20*m_l2*M4 + 
				       4*pow(F12,2)*pow(F20,2)*m_l2*M4 + 
				       8*pow(F10,2)*F12*F22*m_l2*M4 + 
				       16*F10*F12*F20*F22*m_l2*M4 + 
				       8*F12*pow(F20,2)*F22*m_l2*M4 + 
				       8*F10*F20*pow(F22,2)*m_l2*M4 + 
				       4*pow(F20,2)*pow(F22,2)*m_l2*M4 + 
				       pow(F20,2)*pow(F22,2)*lfpf*lipf*pipf + 
				       pow(F20,2)*pow(F22,2)*lfpi*lipf*pipf + 
				       pow(F20,2)*pow(F22,2)*lfpf*lipi*pipf + 
				       pow(F20,2)*pow(F22,2)*lfpi*lipi*pipf - 
				       4*pow(F12,2)*pow(F20,2)*m_l2*M2*pipf - 
				       12*F10*F12*F20*F22*m_l2*M2*pipf - 
				       12*F12*pow(F20,2)*F22*m_l2*M2*pipf - 
				       6*pow(F10,2)*pow(F22,2)*m_l2*M2*pipf - 
				       12*F10*F20*pow(F22,2)*m_l2*M2*pipf - 
				       7*pow(F20,2)*pow(F22,2)*m_l2*M2*pipf + 
				       pow(F20,2)*pow(F22,2)*m_l2*pow(pipf,2) + 
				       lilf*(4*pow(F10,2)*pow(F22,2)*M4 - 
					     (4*pow(F12,2)*pow(F20,2) + 
					      4*F12*F20*(F10 + F20)*F22 + 
					      (-2*pow(F10,2) + 4*F10*F20 + pow(F20,2))*
					      pow(F22,2))*M2*pipf - 
					     pow(F20,2)*pow(F22,2)*pow(pipf,2))) + 
			2*pow(F10,2)*M4*
			(-4*F12*F22*lfpf*lipf*M2 - 
			 3*pow(F22,2)*lfpf*lipf*M2 + 
			 4*pow(F12,2)*lfpi*lipf*M2 + 
			 4*F12*F22*lfpi*lipf*M2 + 
			 pow(F22,2)*lfpi*lipf*M2 + 
			 4*pow(F12,2)*lfpf*lipi*M2 + 
			 4*F12*F22*lfpf*lipi*M2 + 
			 pow(F22,2)*lfpf*lipi*M2 - 
			 4*F12*F22*lfpi*lipi*M2 - 
			 3*pow(F22,2)*lfpi*lipi*M2 + 
			 8*pow(F12,2)*m_l2*M4 + 
			 12*F12*F22*m_l2*M4 + 
			 5*pow(F22,2)*m_l2*M4 + 
			 pow(F22,2)*lfpf*lipf*pipf + pow(F22,2)*lfpi*lipf*pipf + 
			 pow(F22,2)*lfpf*lipi*pipf + pow(F22,2)*lfpi*lipi*pipf - 
			 4*pow(F12,2)*m_l2*M2*pipf - 
			 12*F12*F22*m_l2*M2*pipf - 
			 6*pow(F22,2)*m_l2*M2*pipf + 
			 pow(F22,2)*m_l2*pow(pipf,2) - 
			 lilf*pow((2*F12 + F22)*M2 - F22*pipf,2) + 
			 kflf*(F22*lipf*(-((4*F12 + 3*F22)*M2) + F22*pipf) + 
			       lipi*(pow(2*F12 + F22,2)*M2 + pow(F22,2)*pipf))) +
			2*kfli*(F22*kflf*(pow(F10,2)*M4*
					  (-((4*F12 + 3*F22)*M2) + F22*pipf) + 
					  pow(kfpf,2)*(-2*pow(F10,2)*F22*M2 + 
						       pow(F20,2)*F22*pipf) + 
					  kfpf*M2*((pow(F10,2)*F22 - 2*F10*F20*(F12 + F22) - 
						    2*pow(F20,2)*(F12 + F22))*M2 + 
						   2*F20*(F10 + F20)*(F12 + F22)*pipf)) + 
				lfpf*M2*(F22*(-(pow(F10,2)*F22) + 3*F10*F20*(F12 + F22) + 
					      3*pow(F20,2)*(F12 + F22))*pow(kfpf,2) + 
					 pow(F10,2)*F22*M2*
					 (-((4*F12 + 3*F22)*M2) + F22*pipf) + 
					 kfpf*((2*pow(F10,2)*F22*(2*F12 + F22) + 
						pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
						F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*M2 - 
					       (4*pow(F12,2)*pow(F20,2) + 
						7*F12*F20*(F10 + F20)*F22 + 
						(3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*pow(F22,2))*pipf)) + 
				lfpi*(-(pow(F20,2)*pow(F22,2)*pow(kfpf,3)) + 
				      F22*(pow(F10,2)*F22 - 3*F10*F20*(F12 + F22) - 
					   3*pow(F20,2)*(F12 + F22))*pow(kfpf,2)*M2 + 
				      pow(F10,2)*M4*
				      (pow(2*F12 + F22,2)*M2 + pow(F22,2)*pipf) + 
				      kfpf*M2*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
							      pow(F20,2)*(2*F12 + F22))*M2) + 
					       F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
						    pow(F20,2)*(F12 + F22))*pipf))) + 
			2*kfpf*M2*(kflf*
				   (lipf*((2*pow(F10,2)*F22*(2*F12 + F22) + 
					   pow(F20,2)*(2*pow(F12,2) + 3*F12*F22 + pow(F22,2)) + 
					   F10*F20*(4*pow(F12,2) + 7*F12*F22 + 3*pow(F22,2)))*
					  M2 - (4*pow(F12,2)*pow(F20,2) + 
						7*F12*F20*(F10 + F20)*F22 + 
						(3*pow(F10,2) + 7*F10*F20 + 3*pow(F20,2))*pow(F22,2))*pipf) + 
				    lipi*(-((F12 + F22)*(4*pow(F10,2)*F12 + F10*F20*(4*F12 - 3*F22) - 
							 pow(F20,2)*(2*F12 + F22))*M2) + 
					  F22*(-(pow(F10,2)*F22) + F10*F20*(F12 + F22) + 
					       pow(F20,2)*(F12 + F22))*pipf)) + 
				   pow(F10,2)*(-8*pow(F12,2)*lfpf*lipf*M2 + 
					       3*pow(F22,2)*lfpf*lipf*M2 - 
					       4*pow(F12,2)*lfpi*lipf*M2 + 
					       pow(F22,2)*lfpi*lipf*M2 + 
					       4*pow(F12,2)*m_l2*M4 + 
					       12*F12*F22*m_l2*M4 + 
					       6*pow(F22,2)*m_l2*M4 - 
					       6*pow(F22,2)*lfpf*lipf*pipf - 
					       4*pow(F22,2)*lfpi*lipf*pipf + 
					       8*pow(F12,2)*m_l2*M2*pipf + 
					       24*F12*F22*m_l2*M2*pipf + 
					       10*pow(F22,2)*m_l2*M2*pipf - 
					       4*pow(F22,2)*m_l2*pow(pipf,2) - 
					       2*F22*lilf*(M2 + 2*pipf)*
					       ((2*F12 + F22)*M2 - F22*pipf) - 
					       lipi*(pow(F22,2)*lfpi*(M2 + 2*pipf) + 
						     lfpf*((4*pow(F12,2) - pow(F22,2))*M2 + 
							   4*pow(F22,2)*pipf))))))/(pow(kfpf,2)*pow(kfpi,2)*M4));

  return (Pow6(e)/Pow2(q_22))*res;
}

//----------------------------------------------------------------------------------------------
// First-order bremsstrahlung, lepton term:
G4double lterm ()
{
  res = (2*(pow(kfli,3)*(-kflf + m_l2)*pow((2*F11 + F21)*M2 - F21*pipf,2) - 
	    pow(kflf,2)*m_l2*(-4*pow(F11,2)*kfpf*lfpi*M2 - 
			      4*F11*F21*kfpf*lfpi*M2 - 
			      pow(F21,2)*kfpf*lfpi*M2 + 
			      4*F11*F21*kfpi*lfpi*M2 + 
			      3*pow(F21,2)*kfpi*lfpi*M2 + 
			      4*pow(F11,2)*lfpi*lipf*M2 + 
			      4*F11*F21*lfpi*lipf*M2 + 
			      pow(F21,2)*lfpi*lipf*M2 - 
			      4*F11*F21*lfpi*lipi*M2 - 
			      3*pow(F21,2)*lfpi*lipi*M2 + 
			      (8*pow(F11,2) + 12*F11*F21 + 5*pow(F21,2))*m_l2*
			      M4 - pow(F21,2)*kfpf*lfpi*pipf - 
			      pow(F21,2)*kfpi*lfpi*pipf + pow(F21,2)*lfpi*lipf*pipf + 
			      pow(F21,2)*lfpi*lipi*pipf - 
			      4*pow(F11,2)*m_l2*M2*pipf - 
			      12*F11*F21*m_l2*M2*pipf - 
			      6*pow(F21,2)*m_l2*M2*pipf + 
			      pow(F21,2)*m_l2*pow(pipf,2) + 
			      kflf*pow((2*F11 + F21)*M2 - F21*pipf,2) - 
			      lilf*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
			      lfpf*((F21*(4*F11 + 3*F21)*kfpf - pow(2*F11 + F21,2)*kfpi - 
				     F21*(4*F11 + 3*F21)*lipf + pow(2*F11 + F21,2)*lipi)*
				    M2 + pow(F21,2)*(-kfpf - kfpi + lipf + lipi)*pipf))
	    + pow(kfli,2)*(-(m_l2*
			     (4*pow(F11,2)*lfpi*lipf*M2 + 
			      4*F11*F21*lfpi*lipf*M2 + 
			      pow(F21,2)*lfpi*lipf*M2 - 
			      4*F11*F21*lfpi*lipi*M2 - 
			      3*pow(F21,2)*lfpi*lipi*M2 - 
			      4*pow(F11,2)*lilf*M4 - 
			      4*F11*F21*lilf*M4 - pow(F21,2)*lilf*M4 + 
			      (8*pow(F11,2) + 12*F11*F21 + 5*pow(F21,2))*m_l2*
			      M4 + pow(F21,2)*lfpi*lipf*pipf + 
			      pow(F21,2)*lfpi*lipi*pipf + 
			      4*F11*F21*lilf*M2*pipf + 
			      2*pow(F21,2)*lilf*M2*pipf - 
			      4*pow(F11,2)*m_l2*M2*pipf - 
			      12*F11*F21*m_l2*M2*pipf - 
			      6*pow(F21,2)*m_l2*M2*pipf - 
			      pow(F21,2)*lilf*pow(pipf,2) + 
			      pow(F21,2)*m_l2*pow(pipf,2) + 
			      kfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
				     pow(2*F11 + F21,2)*lipi)*M2 + 
				    pow(F21,2)*(lipf + lipi)*pipf) + 
			      lfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
				     pow(2*F11 + F21,2)*lipi)*M2 + 
				    pow(F21,2)*(lipf + lipi)*pipf) + 
			      kfpi*((pow(2*F11 + F21,2)*lipf - 
				     F21*(4*F11 + 3*F21)*lipi)*M2 + 
				    pow(F21,2)*(lipf + lipi)*pipf))) + 
			   kflf*(-4*F11*F21*pow(lfpi,2)*M2 - 
				 3*pow(F21,2)*pow(lfpi,2)*M2 + 
				 4*pow(F11,2)*lfpi*lipf*M2 + 
				 4*F11*F21*lfpi*lipf*M2 + 
				 pow(F21,2)*lfpi*lipf*M2 - 
				 4*F11*F21*lfpi*lipi*M2 - 
				 3*pow(F21,2)*lfpi*lipi*M2 - 
				 8*pow(F11,2)*lilf*M4 - 8*F11*F21*lilf*M4 - 
				 2*pow(F21,2)*lilf*M4 - 
				 pow(2*F11 + F21,2)*m_l2*M4 + 
				 pow(F21,2)*pow(lfpi,2)*pipf + 
				 pow(F21,2)*lfpi*lipf*pipf + pow(F21,2)*lfpi*lipi*pipf + 
				 8*F11*F21*lilf*M2*pipf + 
				 4*pow(F21,2)*lilf*M2*pipf + 
				 4*F11*F21*m_l2*M2*pipf + 
				 2*pow(F21,2)*m_l2*M2*pipf - 
				 2*pow(F21,2)*lilf*pow(pipf,2) - 
				 pow(F21,2)*m_l2*pow(pipf,2) + 
				 F21*pow(lfpf,2)*(-((4*F11 + 3*F21)*M2) + F21*pipf) + 
				 kfpf*((-(F21*(4*F11 + 3*F21)*lipf) + 
					pow(2*F11 + F21,2)*lipi)*M2 + 
				       pow(F21,2)*(lipf + lipi)*pipf) + 
				 kfpi*((pow(2*F11 + F21,2)*lipf - F21*(4*F11 + 3*F21)*lipi)*
				       M2 + pow(F21,2)*(lipf + lipi)*pipf) + 
				 lfpf*((2*pow(2*F11 + F21,2)*lfpi - 
					F21*(4*F11 + 3*F21)*lipf + pow(2*F11 + F21,2)*lipi)*
				       M2 + pow(F21,2)*(2*lfpi + lipf + lipi)*pipf))) + 
	    kflf*kfli*(-4*pow(F11,2)*kfpf*lfpi*lilf*M2 - 
		       4*F11*F21*kfpf*lfpi*lilf*M2 - 
		       pow(F21,2)*kfpf*lfpi*lilf*M2 - 
		       4*F11*F21*kfpf*lilf*lipf*M2 - 
		       3*pow(F21,2)*kfpf*lilf*lipf*M2 + 
		       8*pow(F11,2)*lfpi*lilf*lipf*M2 + 
		       8*F11*F21*lfpi*lilf*lipf*M2 + 
		       2*pow(F21,2)*lfpi*lilf*lipf*M2 + 
		       4*pow(F11,2)*kfpf*lilf*lipi*M2 + 
		       4*F11*F21*kfpf*lilf*lipi*M2 + 
		       pow(F21,2)*kfpf*lilf*lipi*M2 - 
		       8*F11*F21*lfpi*lilf*lipi*M2 - 
		       6*pow(F21,2)*lfpi*lilf*lipi*M2 + 
		       4*F11*F21*pow(kfpf,2)*m_l2*M2 + 
		       3*pow(F21,2)*pow(kfpf,2)*m_l2*M2 - 
		       8*pow(F11,2)*pow(lilf,2)*M4 - 
		       8*F11*F21*pow(lilf,2)*M4 - 
		       2*pow(F21,2)*pow(lilf,2)*M4 + 
		       16*pow(F11,2)*lilf*m_l2*M4 + 
		       24*F11*F21*lilf*m_l2*M4 + 
		       10*pow(F21,2)*lilf*m_l2*M4 - 
		       pow(F21,2)*kfpf*lfpi*lilf*pipf + 
		       pow(F21,2)*kfpf*lilf*lipf*pipf + 
		       2*pow(F21,2)*lfpi*lilf*lipf*pipf + 
		       pow(F21,2)*kfpf*lilf*lipi*pipf + 
		       2*pow(F21,2)*lfpi*lilf*lipi*pipf - 
		       pow(F21,2)*pow(kfpf,2)*m_l2*pipf + 
		       8*F11*F21*pow(lilf,2)*M2*pipf + 
		       4*pow(F21,2)*pow(lilf,2)*M2*pipf - 
		       8*pow(F11,2)*lilf*m_l2*M2*pipf - 
		       24*F11*F21*lilf*m_l2*M2*pipf - 
		       12*pow(F21,2)*lilf*m_l2*M2*pipf - 
		       2*pow(F21,2)*pow(lilf,2)*pow(pipf,2) + 
		       2*pow(F21,2)*lilf*m_l2*pow(pipf,2) - 
		       pow(kflf,2)*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
		       F21*pow(kfpi,2)*m_l2*((4*F11 + 3*F21)*M2 - F21*pipf) + 
		       lfpf*lilf*((F21*(4*F11 + 3*F21)*kfpf - 
				   2*F21*(4*F11 + 3*F21)*lipf + 2*pow(2*F11 + F21,2)*lipi)*
				  M2 + pow(F21,2)*(-kfpf + 2*(lipf + lipi))*pipf) + 
		       kfpi*(-4*F11*F21*kflf*lfpi*M2 - 
			     3*pow(F21,2)*kflf*lfpi*M2 + 
			     4*F11*F21*lfpi*lilf*M2 + 
			     3*pow(F21,2)*lfpi*lilf*M2 + 
			     4*pow(F11,2)*lilf*lipf*M2 + 
			     4*F11*F21*lilf*lipf*M2 + 
			     pow(F21,2)*lilf*lipf*M2 + 
			     pow(F21,2)*kflf*lfpi*pipf - pow(F21,2)*lfpi*lilf*pipf + 
			     pow(F21,2)*lilf*lipf*pipf + 
			     F21*lilf*lipi*(-((4*F11 + 3*F21)*M2) + F21*pipf) + 
			     lfpf*(kflf - lilf)*(pow(2*F11 + F21,2)*M2 + pow(F21,2)*pipf) - 
			     2*kfpf*m_l2*(pow(2*F11 + F21,2)*M2 + 
					  pow(F21,2)*pipf)) + 
		       kflf*(4*pow(F11,2)*kfpf*lfpi*M2 + 
			     4*F11*F21*kfpf*lfpi*M2 + 
			     pow(F21,2)*kfpf*lfpi*M2 - 
			     4*pow(F11,2)*lfpi*lipf*M2 - 
			     4*F11*F21*lfpi*lipf*M2 - 
			     pow(F21,2)*lfpi*lipf*M2 + 
			     4*F11*F21*pow(lipf,2)*M2 + 
			     3*pow(F21,2)*pow(lipf,2)*M2 + 
			     pow(2*F11 + F21,2)*m_l2*M4 + 
			     pow(F21,2)*kfpf*lfpi*pipf - pow(F21,2)*lfpi*lipf*pipf - 
			     pow(F21,2)*pow(lipf,2)*pipf - 
			     4*F11*F21*m_l2*M2*pipf - 
			     2*pow(F21,2)*m_l2*M2*pipf + 
			     pow(F21,2)*m_l2*pow(pipf,2) + 
			     2*lilf*pow((2*F11 + F21)*M2 - F21*pipf,2) + 
			     F21*lfpf*(-kfpf + lipf)*((4*F11 + 3*F21)*M2 - F21*pipf) + 
			     F21*pow(lipi,2)*((4*F11 + 3*F21)*M2 - F21*pipf) + 
			     lipi*(F21*lfpi*((4*F11 + 3*F21)*M2 - F21*pipf) + 
				   2*lipf*(-(pow(2*F11 + F21,2)*M2) - 
					   pow(F21,2)*pipf) - 
				   lfpf*(pow(2*F11 + F21,2)*M2 + pow(F21,2)*pipf))))))/(pow(kflf,2)*pow(kfli,2)*M2);

  return (Pow6(e)/Pow2(q_12))*res;
}

//----------------------------------------------------------------------------------------------
// First-order bremsstrahlung, interference term:
G4double lpterm ()
{
  res1 = (F20*F21*F22*kflf*pow(kfpf,2)*lilf + 
	  2*F20*F21*F22*pow(kfpf,2)*lfpf*lilf - 
	  4*F11*F20*F22*pow(kfpf,2)*lfpi*lilf + 
	  2*F10*F21*F22*pow(kfpf,2)*lfpi*lilf - 
	  2*F20*F21*F22*pow(kfpf,2)*lfpi*lilf - 
	  F20*F21*F22*pow(kflf,2)*kfpf*lipf - 
	  4*F20*F21*F22*kflf*pow(kfpf,2)*lipf - 
	  4*F20*F21*F22*kflf*kfpf*lfpf*lipf + 2*F12*F20*F21*kflf*kfpf*lfpi*lipf - 
	  2*F11*F20*F22*kflf*kfpf*lfpi*lipf - 2*F10*F21*F22*kflf*kfpf*lfpi*lipf - 
	  4*F20*F21*F22*kflf*kfpf*lfpi*lipf + 
	  4*F12*F20*F21*pow(kfpf,2)*lfpi*lipf - 
	  4*F11*F20*F22*pow(kfpf,2)*lfpi*lipf - 
	  2*F10*F21*F22*pow(kfpf,2)*lfpi*lipf + 
	  8*F12*F20*F21*kfpf*lfpf*lfpi*lipf - 8*F11*F20*F22*kfpf*lfpf*lfpi*lipf + 
	  4*F12*F20*F21*kfpf*pow(lfpi,2)*lipf - 
	  4*F11*F20*F22*kfpf*pow(lfpi,2)*lipf + 
	  2*F11*F20*F22*pow(kflf,2)*kfpf*lipi + 
	  F20*F21*F22*pow(kflf,2)*kfpf*lipi + 
	  12*F11*F20*F22*kflf*pow(kfpf,2)*lipi - 
	  2*F10*F21*F22*kflf*pow(kfpf,2)*lipi + 
	  8*F20*F21*F22*kflf*pow(kfpf,2)*lipi + 
	  2*F12*F20*F21*kflf*kfpf*lfpf*lipi + 
	  10*F11*F20*F22*kflf*kfpf*lfpf*lipi + 
	  8*F20*F21*F22*kflf*kfpf*lfpf*lipi - 
	  4*F10*F21*F22*pow(kfpf,2)*lfpf*lipi - 
	  4*F10*F21*F22*kfpf*pow(lfpf,2)*lipi + 
	  4*F12*F20*F21*kflf*kfpf*lfpi*lipi + 4*F11*F20*F22*kflf*kfpf*lfpi*lipi - 
	  2*F10*F21*F22*kflf*kfpf*lfpi*lipi + 4*F20*F21*F22*kflf*kfpf*lfpi*lipi - 
	  4*F11*F20*F22*pow(kfpf,2)*lfpi*lipi - 
	  6*F10*F21*F22*pow(kfpf,2)*lfpi*lipi - 
	  4*F20*F21*F22*pow(kfpf,2)*lfpi*lipi + 
	  4*F12*F20*F21*kfpf*lfpf*lfpi*lipi - 4*F11*F20*F22*kfpf*lfpf*lfpi*lipi - 
	  8*F10*F21*F22*kfpf*lfpf*lfpi*lipi - 
	  4*F10*F21*F22*kfpf*pow(lfpi,2)*lipi - 
	  F20*F21*F22*kflf*pow(kfpf,2)*m_l2 + 
	  2*F20*F21*F22*pow(kfpf,3)*m_l2 + 
	  4*F11*F20*F22*pow(kfpf,2)*lfpi*m_l2 - 
	  2*F10*F21*F22*pow(kfpf,2)*lfpi*m_l2 + 
	  4*F20*F21*F22*pow(kfpf,2)*lfpi*m_l2 + 
	  2*F20*F21*F22*pow(kfpf,2)*lipf*m_l2 - 
	  8*F11*F20*F22*pow(kfpf,2)*lipi*m_l2 - 
	  6*F20*F21*F22*pow(kfpf,2)*lipi*m_l2 + 
	  8*F11*F12*F20*kflf*kfpf*lilf*M2 + 
	  4*F12*F20*F21*kflf*kfpf*lilf*M2 + 
	  4*F10*F11*F22*kflf*kfpf*lilf*M2 + 
	  4*F11*F20*F22*kflf*kfpf*lilf*M2 + 
	  2*F10*F21*F22*kflf*kfpf*lilf*M2 + 
	  2*F20*F21*F22*kflf*kfpf*lilf*M2 - 
	  8*F11*F12*F20*pow(kfpf,2)*lilf*M2 - 
	  6*F12*F20*F21*pow(kfpf,2)*lilf*M2 - 
	  12*F10*F11*F22*pow(kfpf,2)*lilf*M2 - 
	  6*F11*F20*F22*pow(kfpf,2)*lilf*M2 - 
	  8*F10*F21*F22*pow(kfpf,2)*lilf*M2 - 
	  4*F20*F21*F22*pow(kfpf,2)*lilf*M2 + 
	  4*F10*F12*F21*kfpf*lfpf*lilf*M2 - 
	  4*F10*F11*F22*kfpf*lfpf*lilf*M2 - 
	  16*F10*F11*F12*kfpf*lfpi*lilf*M2 - 
	  12*F10*F12*F21*kfpf*lfpi*lilf*M2 - 
	  4*F10*F11*F22*kfpf*lfpi*lilf*M2 + 
	  2*F12*F20*F21*pow(kflf,2)*lipf*M2 - 
	  4*F10*F11*F22*pow(kflf,2)*lipf*M2 + 
	  2*F11*F20*F22*pow(kflf,2)*lipf*M2 - 
	  2*F10*F21*F22*pow(kflf,2)*lipf*M2 + 
	  4*F20*F21*F22*pow(kflf,2)*lipf*M2 - 
	  8*F11*F12*F20*kflf*kfpf*lipf*M2 - 
	  8*F10*F12*F21*kflf*kfpf*lipf*M2 - 
	  10*F12*F20*F21*kflf*kfpf*lipf*M2 - 
	  4*F10*F11*F22*kflf*kfpf*lipf*M2 - 
	  6*F11*F20*F22*kflf*kfpf*lipf*M2 - 
	  6*F10*F21*F22*kflf*kfpf*lipf*M2 - 
	  8*F20*F21*F22*kflf*kfpf*lipf*M2 - 
	  12*F10*F12*F21*kflf*lfpf*lipf*M2 - 
	  20*F10*F11*F22*kflf*lfpf*lipf*M2 - 
	  24*F10*F21*F22*kflf*lfpf*lipf*M2 - 
	  16*F10*F12*F21*kfpf*lfpf*lipf*M2 - 
	  16*F10*F11*F22*kfpf*lfpf*lipf*M2 - 
	  24*F10*F21*F22*kfpf*lfpf*lipf*M2 - 
	  16*F10*F12*F21*pow(lfpf,2)*lipf*M2 - 
	  16*F10*F11*F22*pow(lfpf,2)*lipf*M2 - 
	  24*F10*F21*F22*pow(lfpf,2)*lipf*M2 + 
	  16*F10*F11*F12*kflf*lfpi*lipf*M2 + 
	  4*F10*F12*F21*kflf*lfpi*lipf*M2 - 
	  4*F12*F20*F21*kflf*lfpi*lipf*M2 + 
	  12*F10*F11*F22*kflf*lfpi*lipf*M2 + 
	  4*F11*F20*F22*kflf*lfpi*lipf*M2 + 
	  4*F10*F21*F22*kflf*lfpi*lipf*M2 + 
	  16*F10*F11*F12*kfpf*lfpi*lipf*M2 + 
	  8*F10*F12*F21*kfpf*lfpi*lipf*M2 + 
	  8*F10*F11*F22*kfpf*lfpi*lipf*M2 + 
	  4*F10*F21*F22*kfpf*lfpi*lipf*M2 + 
	  32*F10*F11*F12*lfpf*lfpi*lipf*M2 + 
	  16*F10*F12*F21*lfpf*lfpi*lipf*M2 + 
	  16*F10*F11*F22*lfpf*lfpi*lipf*M2 + 
	  8*F10*F21*F22*lfpf*lfpi*lipf*M2 - 
	  8*F11*F12*F20*pow(kflf,2)*lipi*M2 - 
	  6*F12*F20*F21*pow(kflf,2)*lipi*M2 - 
	  4*F10*F11*F22*pow(kflf,2)*lipi*M2 - 
	  6*F11*F20*F22*pow(kflf,2)*lipi*M2 - 
	  2*F10*F21*F22*pow(kflf,2)*lipi*M2 - 
	  4*F20*F21*F22*pow(kflf,2)*lipi*M2 + 
	  32*F10*F11*F12*kflf*kfpf*lipi*M2 + 
	  8*F11*F12*F20*kflf*kfpf*lipi*M2 + 
	  24*F10*F12*F21*kflf*kfpf*lipi*M2 + 
	  6*F12*F20*F21*kflf*kfpf*lipi*M2 + 
	  4*F10*F11*F22*kflf*kfpf*lipi*M2 + 
	  6*F11*F20*F22*kflf*kfpf*lipi*M2 + 
	  4*F10*F21*F22*kflf*kfpf*lipi*M2 + 
	  4*F20*F21*F22*kflf*kfpf*lipi*M2 + 
	  32*F10*F11*F12*kflf*lfpf*lipi*M2 + 
	  20*F10*F12*F21*kflf*lfpf*lipi*M2 + 
	  12*F10*F11*F22*kflf*lfpf*lipi*M2 + 
	  8*F10*F21*F22*kflf*lfpf*lipi*M2 + 
	  32*F10*F11*F12*kfpf*lfpf*lipi*M2 + 
	  16*F10*F12*F21*kfpf*lfpf*lipi*M2 + 
	  16*F10*F11*F22*kfpf*lfpf*lipi*M2 + 
	  8*F10*F21*F22*kfpf*lfpf*lipi*M2 + 
	  32*F10*F11*F12*pow(lfpf,2)*lipi*M2 + 
	  16*F10*F12*F21*pow(lfpf,2)*lipi*M2 + 
	  16*F10*F11*F22*pow(lfpf,2)*lipi*M2 + 
	  8*F10*F21*F22*pow(lfpf,2)*lipi*M2 - 
	  12*F10*F12*F21*kflf*lfpi*lipi*M2 - 
	  4*F12*F20*F21*kflf*lfpi*lipi*M2 - 
	  4*F10*F11*F22*kflf*lfpi*lipi*M2 + 
	  4*F11*F20*F22*kflf*lfpi*lipi*M2 - 
	  12*F10*F21*F22*kflf*lfpi*lipi*M2 - 
	  8*F10*F12*F21*kfpf*lfpi*lipi*M2 - 
	  8*F10*F11*F22*kfpf*lfpi*lipi*M2 - 
	  12*F10*F21*F22*kfpf*lfpi*lipi*M2 - 
	  16*F10*F12*F21*lfpf*lfpi*lipi*M2 - 
	  16*F10*F11*F22*lfpf*lfpi*lipi*M2 - 
	  24*F10*F21*F22*lfpf*lfpi*lipi*M2 - 
	  16*F11*F12*F20*kflf*kfpf*m_l2*M2 - 
	  12*F12*F20*F21*kflf*kfpf*m_l2*M2 - 
	  12*F11*F20*F22*kflf*kfpf*m_l2*M2 + 
	  2*F10*F21*F22*kflf*kfpf*m_l2*M2 - 
	  10*F20*F21*F22*kflf*kfpf*m_l2*M2 + 
	  16*F11*F12*F20*pow(kfpf,2)*m_l2*M2 + 
	  4*F10*F12*F21*pow(kfpf,2)*m_l2*M2 + 
	  14*F12*F20*F21*pow(kfpf,2)*m_l2*M2 + 
	  28*F10*F11*F22*pow(kfpf,2)*m_l2*M2 + 
	  14*F11*F20*F22*pow(kfpf,2)*m_l2*M2 + 
	  24*F10*F21*F22*pow(kfpf,2)*m_l2*M2 + 
	  12*F20*F21*F22*pow(kfpf,2)*m_l2*M2 + 
	  24*F10*F11*F22*kfpf*lfpf*m_l2*M2 + 
	  20*F10*F21*F22*kfpf*lfpf*m_l2*M2 + 
	  16*F10*F11*F12*kfpf*lfpi*m_l2*M2 + 
	  16*F10*F12*F21*kfpf*lfpi*m_l2*M2 + 
	  8*F10*F11*F22*kfpf*lfpi*m_l2*M2 + 
	  4*F10*F21*F22*kfpf*lfpi*m_l2*M2 + 
	  4*F10*F12*F21*kfpf*lipf*m_l2*M2 + 
	  4*F10*F11*F22*kfpf*lipf*m_l2*M2 + 
	  4*F10*F21*F22*kfpf*lipf*m_l2*M2 - 
	  16*F10*F11*F12*kfpf*lipi*m_l2*M2 - 
	  12*F10*F12*F21*kfpf*lipi*m_l2*M2 + 
	  4*F10*F11*F22*kfpf*lipi*m_l2*M2 + 
	  4*F10*F21*F22*kfpf*lipi*m_l2*M2 - 
	  16*F10*F11*F12*kflf*lilf*M4 - 
	  8*F10*F12*F21*kflf*lilf*M4 - 
	  8*F10*F11*F22*kflf*lilf*M4 - 
	  4*F10*F21*F22*kflf*lilf*M4 - 
	  16*F10*F11*F12*kfpf*lilf*M4 - 
	  8*F10*F12*F21*kfpf*lilf*M4 - 
	  8*F10*F11*F22*kfpf*lilf*M4 - 
	  4*F10*F21*F22*kfpf*lilf*M4 - 
	  32*F10*F11*F12*lfpf*lilf*M4 - 
	  16*F10*F12*F21*lfpf*lilf*M4 - 
	  16*F10*F11*F22*lfpf*lilf*M4 - 
	  8*F10*F21*F22*lfpf*lilf*M4 + 
	  16*F10*F11*F12*kflf*lipf*M4 + 
	  16*F10*F12*F21*kflf*lipf*M4 + 
	  16*F10*F11*F22*kflf*lipf*M4 + 
	  16*F10*F21*F22*kflf*lipf*M4 - 
	  16*F10*F11*F12*kflf*lipi*M4 - 
	  8*F10*F12*F21*kflf*lipi*M4 - 
	  8*F10*F11*F22*kflf*lipi*M4 - 
	  4*F10*F21*F22*kflf*lipi*M4 + 
	  32*F10*F11*F12*kflf*m_l2*M4 + 
	  24*F10*F12*F21*kflf*m_l2*M4 + 
	  24*F10*F11*F22*kflf*m_l2*M4 + 
	  20*F10*F21*F22*kflf*m_l2*M4 + 
	  32*F10*F11*F12*kfpf*m_l2*M4 + 
	  24*F10*F12*F21*kfpf*m_l2*M4 + 
	  24*F10*F11*F22*kfpf*m_l2*M4 + 
	  20*F10*F21*F22*kfpf*m_l2*M4 + 
	  64*F10*F11*F12*lfpf*m_l2*M4 + 
	  48*F10*F12*F21*lfpf*m_l2*M4 + 
	  48*F10*F11*F22*lfpf*m_l2*M4 + 
	  40*F10*F21*F22*lfpf*m_l2*M4 + 
	  pow(kfpi,2)*(F21*F22*kflf*(F20*lilf + 2*F10*lipf - F20*m_l2) + 
		       2*(-2*F12*F20*F21*lfpf*lipf + F10*F21*F22*lfpf*lipf - 
			  2*F20*F21*F22*lfpf*lipf + F10*F21*F22*lfpf*lipi + 
			  F20*F21*F22*kfpf*m_l2 + F10*F21*F22*lfpf*m_l2 + 
			  2*F10*F12*F21*m_l2*M2 + 
			  F12*F20*F21*m_l2*M2 + 
			  2*F10*F11*F22*m_l2*M2 + 
			  F11*F20*F22*m_l2*M2 + 
			  4*F10*F21*F22*m_l2*M2 + 
			  2*F20*F21*F22*m_l2*M2 + 
			  lilf*(-(F10*F21*F22*lfpf) + 
				(F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2))) -
	  4*F12*F20*F21*kflf*kfpf*lilf*pipf - 2*F10*F21*F22*kflf*kfpf*lilf*pipf - 
	  2*F20*F21*F22*kflf*kfpf*lilf*pipf + 
	  4*F11*F20*F22*pow(kfpf,2)*lilf*pipf + 
	  4*F10*F21*F22*pow(kfpf,2)*lilf*pipf + 
	  4*F20*F21*F22*pow(kfpf,2)*lilf*pipf - 
	  4*F12*F20*F21*kfpf*lfpf*lilf*pipf + 4*F11*F20*F22*kfpf*lfpf*lilf*pipf - 
	  4*F12*F20*F21*kfpf*lfpi*lilf*pipf + 4*F11*F20*F22*kfpf*lfpi*lilf*pipf + 
	  2*F12*F20*F21*pow(kflf,2)*lipf*pipf - 
	  2*F11*F20*F22*pow(kflf,2)*lipf*pipf + 
	  2*F10*F21*F22*pow(kflf,2)*lipf*pipf + 
	  4*F12*F20*F21*kflf*kfpf*lipf*pipf + 6*F10*F21*F22*kflf*kfpf*lipf*pipf + 
	  4*F20*F21*F22*kflf*kfpf*lipf*pipf + 4*F12*F20*F21*kflf*lfpf*lipf*pipf - 
	  4*F11*F20*F22*kflf*lfpf*lipf*pipf + 8*F10*F21*F22*kflf*lfpf*lipf*pipf + 
	  8*F10*F21*F22*kfpf*lfpf*lipf*pipf + 
	  8*F10*F21*F22*pow(lfpf,2)*lipf*pipf + 
	  4*F10*F21*F22*kflf*lfpi*lipf*pipf + 4*F10*F21*F22*kfpf*lfpi*lipf*pipf + 
	  8*F10*F21*F22*lfpf*lfpi*lipf*pipf + 
	  2*F12*F20*F21*pow(kflf,2)*lipi*pipf - 
	  2*F11*F20*F22*pow(kflf,2)*lipi*pipf + 
	  2*F10*F21*F22*pow(kflf,2)*lipi*pipf + 
	  4*F12*F20*F21*kflf*kfpf*lipi*pipf - 4*F11*F20*F22*kflf*kfpf*lipi*pipf + 
	  8*F10*F21*F22*kflf*kfpf*lipi*pipf + 4*F12*F20*F21*kflf*lfpf*lipi*pipf - 
	  4*F11*F20*F22*kflf*lfpf*lipi*pipf + 8*F10*F21*F22*kflf*lfpf*lipi*pipf + 
	  8*F10*F21*F22*kfpf*lfpf*lipi*pipf + 
	  8*F10*F21*F22*pow(lfpf,2)*lipi*pipf + 
	  4*F10*F21*F22*kflf*lfpi*lipi*pipf + 4*F10*F21*F22*kfpf*lfpi*lipi*pipf + 
	  8*F10*F21*F22*lfpf*lfpi*lipi*pipf + 
	  4*F12*F20*F21*kflf*kfpf*m_l2*pipf + 
	  4*F11*F20*F22*kflf*kfpf*m_l2*pipf + 
	  2*F10*F21*F22*kflf*kfpf*m_l2*pipf + 
	  6*F20*F21*F22*kflf*kfpf*m_l2*pipf - 
	  4*F12*F20*F21*pow(kfpf,2)*m_l2*pipf - 
	  4*F11*F20*F22*pow(kfpf,2)*m_l2*pipf - 
	  8*F10*F21*F22*pow(kfpf,2)*m_l2*pipf - 
	  8*F20*F21*F22*pow(kfpf,2)*m_l2*pipf - 
	  4*F10*F21*F22*kfpf*lfpf*m_l2*pipf - 
	  4*F10*F21*F22*kfpf*lfpi*m_l2*pipf - 
	  4*F12*F20*F21*kfpf*lipf*m_l2*pipf + 
	  4*F11*F20*F22*kfpf*lipf*m_l2*pipf - 
	  4*F10*F21*F22*kfpf*lipf*m_l2*pipf - 
	  4*F12*F20*F21*kfpf*lipi*m_l2*pipf + 
	  4*F11*F20*F22*kfpf*lipi*m_l2*pipf - 
	  4*F10*F21*F22*kfpf*lipi*m_l2*pipf + 
	  8*F10*F12*F21*kflf*lilf*M2*pipf + 
	  8*F10*F11*F22*kflf*lilf*M2*pipf + 
	  8*F10*F21*F22*kflf*lilf*M2*pipf + 
	  8*F10*F12*F21*kfpf*lilf*M2*pipf + 
	  8*F10*F11*F22*kfpf*lilf*M2*pipf + 
	  8*F10*F21*F22*kfpf*lilf*M2*pipf + 
	  16*F10*F12*F21*lfpf*lilf*M2*pipf + 
	  16*F10*F11*F22*lfpf*lilf*M2*pipf + 
	  16*F10*F21*F22*lfpf*lilf*M2*pipf - 
	  16*F10*F11*F12*kflf*lipf*M2*pipf - 
	  16*F10*F12*F21*kflf*lipf*M2*pipf - 
	  16*F10*F11*F22*kflf*lipf*M2*pipf - 
	  16*F10*F21*F22*kflf*lipf*M2*pipf + 
	  8*F10*F12*F21*kflf*lipi*M2*pipf + 
	  8*F10*F11*F22*kflf*lipi*M2*pipf + 
	  8*F10*F21*F22*kflf*lipi*M2*pipf - 
	  16*F10*F11*F12*kflf*m_l2*M2*pipf - 
	  24*F10*F12*F21*kflf*m_l2*M2*pipf - 
	  24*F10*F11*F22*kflf*m_l2*M2*pipf - 
	  24*F10*F21*F22*kflf*m_l2*M2*pipf - 
	  16*F10*F11*F12*kfpf*m_l2*M2*pipf - 
	  24*F10*F12*F21*kfpf*m_l2*M2*pipf - 
	  24*F10*F11*F22*kfpf*m_l2*M2*pipf - 
	  24*F10*F21*F22*kfpf*m_l2*M2*pipf - 
	  32*F10*F11*F12*lfpf*m_l2*M2*pipf - 
	  48*F10*F12*F21*lfpf*m_l2*M2*pipf - 
	  48*F10*F11*F22*lfpf*m_l2*M2*pipf - 
	  48*F10*F21*F22*lfpf*m_l2*M2*pipf - 
	  4*F10*F21*F22*kflf*lilf*pow(pipf,2) - 
	  4*F10*F21*F22*kfpf*lilf*pow(pipf,2) - 
	  8*F10*F21*F22*lfpf*lilf*pow(pipf,2) - 
	  4*F10*F21*F22*kflf*lipi*pow(pipf,2) + 
	  4*F10*F21*F22*kflf*m_l2*pow(pipf,2) + 
	  4*F10*F21*F22*kfpf*m_l2*pow(pipf,2) + 
	  8*F10*F21*F22*lfpf*m_l2*pow(pipf,2) + 
	  kfli*(4*F20*(F11 + F21)*F22*pow(kflf,2)*(M2 - pipf) + 
		kfpi*(F22*kflf*(2*F11*F20*lfpf + F20*F21*lfpf - F20*F21*lfpi + 
				8*F10*F11*M2 + 4*F11*F20*M2 + 
				6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) +
		      2*((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lfpf,2) - 
			 2*(F10 + F20)*(F11 + F21)*F22*m_l2*M2 + 
			 lfpi*(F10*F21*F22*lfpf + 
			       (-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*M2) + 
			 lfpf*((F11*(-4*F12*F20 + 2*F10*F22 - 3*F20*F22) + 
				F21*(-3*F12*F20 + F10*F22 - 2*F20*F22))*M2 + 
			       F21*(2*F12*F20 - F10*F22 + 2*F20*F22)*pipf))) - 
		kflf*(F22*kfpf*(F20*F21*lfpf - F20*(2*F11 + F21)*lfpi + 
				8*F10*F11*M2 + 8*F11*F20*M2 + 
				6*F10*F21*M2 + 4*F20*F21*M2 + 
				(4*F11*F20 - 2*F10*F21)*pipf) + 
		      2*(M2*(-((F12*F20*F21 + 
				(-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*lfpf) +
			     2*(2*F11 + F21)*(F20*(F12 + F22) + F10*(2*F12 + F22))*M2) - 
			 ((F12*F20*F21 + (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*
			  lfpf + 4*((F10 + F20)*F21*(F12 + F22) + 
				    F11*(F12*F20 + (F10 + F20)*F22))*M2)*pipf + 
			 2*F21*(F12*F20 + (F10 + F20)*F22)*pow(pipf,2) + 
			 lfpi*((F21*(3*F12*F20 - 3*F10*F22 - 2*F20*F22) + 
				F11*(4*F12*F20 - (2*F10 + F20)*F22))*M2 - 
			       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) - 
		2*(kfpf*(F20*F21*F22*pow(lfpf,2) + 
			 (F10 + F20)*F21*F22*pow(lfpi,2) - 
			 4*F11*F12*F20*lfpf*M2 - 
			 3*F12*F20*F21*lfpf*M2 + 
			 2*F10*F11*F22*lfpf*M2 - 
			 3*F11*F20*F22*lfpf*M2 + 
			 2*F10*F21*F22*lfpf*M2 - 
			 2*F20*F21*F22*lfpf*M2 - 
			 2*F10*F11*F22*m_l2*M2 - 
			 2*F10*F21*F22*m_l2*M2 - 
			 2*F20*(F11 + F21)*F22*(-lfpf + m_l2)*pipf + 
			 lfpi*((-2*F11*F20*F22 + F10*F21*F22)*lfpf + 
			       (F12*F20*F21 + 2*F10*F11*F22 - F11*F20*F22 + 
				F10*F21*F22)*M2 - F10*F21*F22*pipf)) + 
		   2*((F10 + F20)*(F12*F21 - F11*F22)*pow(lfpi,2)*M2 - 
		      lfpf*lfpi*((-(F12*F20*F21) + F11*F20*F22 + 
				  F10*(2*F11 + F21)*(2*F12 + F22))*M2 + 
				 (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
		      lfpf*(F10*M2*((F12*F21 + 3*(F11 + F21)*F22)*lfpf + 
				    (2*F11 + F21)*(2*F12 + F22)*M2) - 
			    ((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf + 
			     2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf + 
			    F10*F21*F22*pow(pipf,2))))) + 
	  kfpi*(F20*F22*pow(kflf,2)*((2*F11 + F21)*lipf - F21*lipi) + 
		2*kflf*(-2*F12*F20*F21*lfpf*lipf + 4*F11*F20*F22*lfpf*lipf + 
			F10*F21*F22*lfpf*lipf + F20*F21*F22*lfpf*lipf - 
			F12*F20*F21*lfpi*lipf + F11*F20*F22*lfpi*lipf - 
			F20*F21*F22*lfpi*lipf + 
			F22*kfpf*(-(F20*(2*F11 + F21)*lilf) + 2*F11*F20*lipf + 
				  F10*F21*lipf - (F10 + 2*F20)*F21*lipi + 
				  2*F11*F20*m_l2 + F20*F21*m_l2) + 
			2*F10*F11*F22*lilf*M2 - 2*F11*F20*F22*lilf*M2 + 
			F10*F21*F22*lilf*M2 - F20*F21*F22*lilf*M2 - 
			4*F10*F12*F21*lipf*M2 - F12*F20*F21*lipf*M2 - 
			2*F10*F11*F22*lipf*M2 + F11*F20*F22*lipf*M2 - 
			4*F10*F21*F22*lipf*M2 + 
			4*F11*F12*F20*m_l2*M2 + 
			4*F12*F20*F21*m_l2*M2 - 
			4*F10*F11*F22*m_l2*M2 + 
			4*F11*F20*F22*m_l2*M2 - 
			3*F10*F21*F22*m_l2*M2 + 
			3*F20*F21*F22*m_l2*M2 - F10*F21*F22*lilf*pipf + 
			F20*F21*F22*lilf*pipf + 2*F12*F20*F21*lipf*pipf + 
			2*F10*F21*F22*lipf*pipf + 2*F20*F21*F22*lipf*pipf + 
			F10*F21*F22*m_l2*pipf - F20*F21*F22*m_l2*pipf - 
			lipi*(F12*F20*F21*lfpf - F11*F20*F22*lfpf - F10*F21*F22*lfpf + 
			      F20*F21*F22*lfpf + F20*F21*F22*lfpi + 
			      4*F10*F12*F21*M2 + 3*F12*F20*F21*M2 - 
			      2*F10*F11*F22*M2 - F11*F20*F22*M2 + 
			      3*F10*F21*F22*M2 + 2*F20*F21*F22*M2 + 
			      F10*F21*F22*pipf)) - 
		2*(2*F20*(2*F11 + F21)*F22*pow(kfpf,2)*m_l2 + 
		   kfpf*(2*F12*F20*F21*lfpf*lipf - 2*F11*F20*F22*lfpf*lipf - 
			 3*F10*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf - 
			 F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf + 
			 2*F11*F20*F22*lfpf*m_l2 - F10*F21*F22*lfpf*m_l2 + 
			 2*F20*F21*F22*lfpf*m_l2 + F10*F21*F22*lfpi*m_l2 - 
			 F20*F21*F22*lipf*m_l2 - 
			 F22*lipi*((2*F11*F20 + F10*F21 + 2*F20*F21)*lfpf - 
				   F10*F21*lfpi + F20*F21*m_l2) + 
			 8*F10*F11*F12*m_l2*M2 + 
			 8*F11*F12*F20*m_l2*M2 + 
			 4*F10*F12*F21*m_l2*M2 + 
			 6*F12*F20*F21*m_l2*M2 + 
			 16*F10*F11*F22*m_l2*M2 + 
			 6*F11*F20*F22*m_l2*M2 + 
			 12*F10*F21*F22*m_l2*M2 + 
			 4*F20*F21*F22*m_l2*M2 - 
			 lilf*(-2*F11*F20*F22*lfpf - F10*F21*F22*lfpf - 
			       F20*F21*F22*lfpf + (F10 + F20)*F21*F22*lfpi + 
			       4*F11*F12*F20*M2 + 4*F12*F20*F21*M2 + 
			       8*F10*F11*F22*M2 + 2*F11*F20*F22*M2 + 
			       6*F10*F21*F22*M2 + 2*F20*F21*F22*M2 - 
			       2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) - 
		   2*(F10*F12*F21*lfpf*lilf*M2 + 
		      F12*F20*F21*lfpf*lilf*M2 + 
		      3*F10*F11*F22*lfpf*lilf*M2 - 
		      F11*F20*F22*lfpf*lilf*M2 + 
		      2*F10*F21*F22*lfpf*lilf*M2 + 
		      F10*F12*F21*lfpi*lilf*M2 + 
		      F12*F20*F21*lfpi*lilf*M2 - 
		      F10*F11*F22*lfpi*lilf*M2 - 
		      F11*F20*F22*lfpi*lilf*M2 - 
		      4*F10*F11*F12*lfpf*m_l2*M2 - 
		      4*F10*F12*F21*lfpf*m_l2*M2 - 
		      8*F10*F11*F22*lfpf*m_l2*M2 - 
		      6*F10*F21*F22*lfpf*m_l2*M2 - 
		      2*F10*F21*F22*lfpf*lilf*pipf + 
		      2*F10*F21*F22*lfpf*m_l2*pipf + 
		      lipi*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*
			    pow(lfpf,2) + (F10 + F20)*(F12*F21 - F11*F22)*m_l2*M2 + 
			    F10*lfpf*(F21*F22*lfpi - 
				      (2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 + 
				      F21*F22*pipf)) + 
		      lipf*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*pow(lfpf,2) + 
			    (F10 + F20)*(F12*F21 - F11*F22)*m_l2*M2 + 
			    lfpf*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpi + 
				  F10*((2*F11 + F21)*(2*F12 + F22)*M2 + 
				       F21*F22*pipf)))))))/(2.*kfpf*M2);

  res2 = (2*F20*F21*F22*pow(kfpi,3)*m_l2 + 
	  pow(kfpi,2)*(F22*kflf*(-(F20*F21*lilf) + 12*F11*F20*lipf - 
				 2*F10*F21*lipf + 8*F20*F21*lipf - 4*F20*F21*lipi + 
				 F20*F21*m_l2) - 
		       2*(2*F11*F20*F22*lfpf*lilf - F10*F21*F22*lfpf*lilf + 
			  F20*F21*F22*lfpf*lilf - F20*F21*F22*lfpi*lilf - 
			  2*F11*F20*F22*lfpf*lipf - 3*F10*F21*F22*lfpf*lipf - 
			  2*F20*F21*F22*lfpf*lipf - 2*F10*F21*F22*lfpi*lipf + 
			  2*F20*(2*F11 + F21)*F22*kfpf*m_l2 - 
			  2*F11*F20*F22*lfpf*m_l2 + F10*F21*F22*lfpf*m_l2 - 
			  2*F20*F21*F22*lfpf*m_l2 + 4*F11*F20*F22*lipf*m_l2 + 
			  3*F20*F21*F22*lipf*m_l2 - 
			  lipi*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*lfpf + 
				F20*F21*F22*m_l2) - 4*F11*F12*F20*lilf*M2 - 
			  3*F12*F20*F21*lilf*M2 - 6*F10*F11*F22*lilf*M2 - 
			  3*F11*F20*F22*lilf*M2 - 4*F10*F21*F22*lilf*M2 - 
			  2*F20*F21*F22*lilf*M2 + 
			  8*F11*F12*F20*m_l2*M2 + 
			  2*F10*F12*F21*m_l2*M2 + 
			  7*F12*F20*F21*m_l2*M2 + 
			  14*F10*F11*F22*m_l2*M2 + 
			  7*F11*F20*F22*m_l2*M2 + 
			  12*F10*F21*F22*m_l2*M2 + 
			  6*F20*F21*F22*m_l2*M2 + 2*F11*F20*F22*lilf*pipf + 
			  2*F10*F21*F22*lilf*pipf + 2*F20*F21*F22*lilf*pipf - 
			  2*F12*F20*F21*m_l2*pipf - 2*F11*F20*F22*m_l2*pipf - 
			  4*F10*F21*F22*m_l2*pipf - 4*F20*F21*F22*m_l2*pipf)) - 
	  pow(kflf,2)*(4*F20*(F11 + F21)*F22*kfli*(M2 - pipf) + 
		       lipf*(-(F20*F21*F22*kfpf) + 
			     2*((4*F11*F12*F20 + 3*F12*F20*F21 + 2*F10*F11*F22 + 
				 3*F11*F20*F22 + F10*F21*F22 + 2*F20*F21*F22)*M2 - 
				(F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)) + 
		       lipi*(F20*(2*F11 + F21)*F22*kfpf - 
			     2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
				(F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
	  kflf*(F21*F22*pow(kfpf,2)*(2*F10*lipi + F20*(-lilf + m_l2)) + 
		2*kfpf*(-(F20*F21*F22*lfpf*lipf) - F12*F20*F21*lfpi*lipf + 
			F11*F20*F22*lfpi*lipf + F10*F21*F22*lfpi*lipf - 
			F20*F21*F22*lfpi*lipf + 4*F10*F12*F21*lipf*M2 + 
			3*F12*F20*F21*lipf*M2 - 2*F10*F11*F22*lipf*M2 - 
			F11*F20*F22*lipf*M2 + 3*F10*F21*F22*lipf*M2 + 
			2*F20*F21*F22*lipf*M2 + 
			4*F11*F12*F20*m_l2*M2 + 
			4*F12*F20*F21*m_l2*M2 - 
			4*F10*F11*F22*m_l2*M2 + 
			4*F11*F20*F22*m_l2*M2 - 
			3*F10*F21*F22*m_l2*M2 + 
			3*F20*F21*F22*m_l2*M2 + F10*F21*F22*lipf*pipf + 
			F10*F21*F22*m_l2*pipf - F20*F21*F22*m_l2*pipf + 
			(F10 - F20)*F22*lilf*((2*F11 + F21)*M2 - F21*pipf) + 
			lipi*(-(F12*F20*F21*lfpf) + F11*F20*F22*lfpf - 
			      F20*F21*F22*lfpf + 
			      (-2*F12*F20*F21 + (4*F11*F20 + (F10 + F20)*F21)*F22)*lfpi + 
			      4*F10*F12*F21*M2 + F12*F20*F21*M2 + 
			      2*F10*F11*F22*M2 - F11*F20*F22*M2 + 
			      4*F10*F21*F22*M2 - 
			      2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) - 
		4*(-3*F10*F12*F21*lfpf*lipf*M2 - 
		   F12*F20*F21*lfpf*lipf*M2 - 
		   F10*F11*F22*lfpf*lipf*M2 + 
		   F11*F20*F22*lfpf*lipf*M2 - 
		   3*F10*F21*F22*lfpf*lipf*M2 + 
		   8*F10*F11*F12*lfpi*lipf*M2 + 
		   5*F10*F12*F21*lfpi*lipf*M2 + 
		   3*F10*F11*F22*lfpi*lipf*M2 + 
		   2*F10*F21*F22*lfpi*lipf*M2 - 
		   4*F10*F11*F12*lilf*M4 - 2*F10*F12*F21*lilf*M4 - 
		   2*F10*F11*F22*lilf*M4 - F10*F21*F22*lilf*M4 + 
		   4*F10*F11*F12*lipf*M4 + 2*F10*F12*F21*lipf*M4 + 
		   2*F10*F11*F22*lipf*M4 + F10*F21*F22*lipf*M4 + 
		   8*F10*F11*F12*m_l2*M4 + 
		   6*F10*F12*F21*m_l2*M4 + 
		   6*F10*F11*F22*m_l2*M4 + 
		   5*F10*F21*F22*m_l2*M4 + 
		   F10*F21*F22*lfpf*lipf*pipf + F12*F20*F21*lfpi*lipf*pipf - 
		   F11*F20*F22*lfpi*lipf*pipf + 2*F10*F21*F22*lfpi*lipf*pipf + 
		   2*F10*F12*F21*lilf*M2*pipf + 
		   2*F10*F11*F22*lilf*M2*pipf + 
		   2*F10*F21*F22*lilf*M2*pipf - 
		   2*F10*F12*F21*lipf*M2*pipf - 
		   2*F10*F11*F22*lipf*M2*pipf - 
		   2*F10*F21*F22*lipf*M2*pipf - 
		   4*F10*F11*F12*m_l2*M2*pipf - 
		   6*F10*F12*F21*m_l2*M2*pipf - 
		   6*F10*F11*F22*m_l2*M2*pipf - 
		   6*F10*F21*F22*m_l2*M2*pipf - 
		   F10*F21*F22*lilf*pow(pipf,2) + 
		   F10*F21*F22*lipf*pow(pipf,2) + 
		   F10*F21*F22*m_l2*pow(pipf,2) + 
		   lipi*(M2*((-(F12*F20*F21) + F11*F20*F22 + 
			      F10*(F21*(F12 + F22) + F11*(4*F12 + 3*F22)))*lfpf - 
			     4*F10*(F11 + F21)*(F12 + F22)*M2) + 
			 F10*(F21*F22*lfpf + 4*(F11 + F21)*(F12 + F22)*M2)*
			 pipf - lfpi*(F10*(3*F12*F21 + 5*F11*F22 + 6*F21*F22)*
				      M2 + (-(F12*F20*F21) + F11*F20*F22 - 2*F10*F21*F22)*pipf))) + 
		kfli*(F22*kfpf*(F20*F21*lfpf - F20*(2*F11 + F21)*lfpi + 
				8*F10*F11*M2 + 4*F11*F20*M2 + 
				6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) +
		      2*(M2*((F21*(-3*F12*F20 + 3*F10*F22 + 2*F20*F22) + 
			      F11*(-4*F12*F20 + (2*F10 + F20)*F22))*lfpf + 
			     2*(2*F11 + F21)*(F20*(F12 + F22) + F10*(2*F12 + F22))*M2) - 
			 (-((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf) + 
			  4*((F10 + F20)*F21*(F12 + F22) + 
			     F11*(F12*F20 + (F10 + F20)*F22))*M2)*pipf + 
			 2*F21*(F12*F20 + (F10 + F20)*F22)*pow(pipf,2) + 
			 lfpi*((F12*F20*F21 + 
				(-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + (F12*F20*F21 + 
											  (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) + 
	  2*(-(pow(kfpf,2)*(F21*lfpi*((-2*F12*F20 + (F10 - 2*F20)*F22)*lipi + 
				      F10*F22*(lipf - m_l2)) + 
			    (2*F10 + F20)*(F12*F21 + (F11 + 2*F21)*F22)*m_l2*
			    M2 + lilf*(F10*F21*F22*lfpi + 
				       (F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2))) + 
	     4*F10*lfpi*(-2*F12*F21*lfpf*lipf*M2 - 
			 2*F11*F22*lfpf*lipf*M2 - 
			 3*F21*F22*lfpf*lipf*M2 + 
			 4*F11*F12*lfpi*lipf*M2 + 
			 2*F12*F21*lfpi*lipf*M2 + 
			 2*F11*F22*lfpi*lipf*M2 + F21*F22*lfpi*lipf*M2 + 
			 8*F11*F12*m_l2*M4 + 
			 6*F12*F21*m_l2*M4 + 
			 6*F11*F22*m_l2*M4 + 
			 5*F21*F22*m_l2*M4 + F21*F22*lfpf*lipf*pipf + 
			 F21*F22*lfpi*lipf*pipf - 4*F11*F12*m_l2*M2*pipf - 
			 6*F12*F21*m_l2*M2*pipf - 
			 6*F11*F22*m_l2*M2*pipf - 
			 6*F21*F22*m_l2*M2*pipf + 
			 F21*F22*m_l2*pow(pipf,2) - 
			 lilf*((2*F11 + F21)*M2 - F21*pipf)*
			 ((2*F12 + F22)*M2 - F22*pipf) + 
			 lipi*(lfpf*((2*F11 + F21)*(2*F12 + F22)*M2 + 
				     F21*F22*pipf) + lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
							   F21*F22*pipf))) - 2*kfpf*(F10*F12*F21*lfpf*lilf*M2 + 
										     F12*F20*F21*lfpf*lilf*M2 - 
										     F10*F11*F22*lfpf*lilf*M2 - 
										     F11*F20*F22*lfpf*lilf*M2 + 
										     F10*F12*F21*lfpi*lilf*M2 + 
										     F12*F20*F21*lfpi*lilf*M2 + 
										     3*F10*F11*F22*lfpi*lilf*M2 - 
										     F11*F20*F22*lfpi*lilf*M2 + 
										     2*F10*F21*F22*lfpi*lilf*M2 - 
										     4*F10*F11*F12*lfpi*m_l2*M2 - 
										     4*F10*F12*F21*lfpi*m_l2*M2 - 
										     8*F10*F11*F22*lfpi*m_l2*M2 - 
										     6*F10*F21*F22*lfpi*m_l2*M2 - 
										     2*F10*F21*F22*lfpi*lilf*pipf + 
										     2*F10*F21*F22*lfpi*m_l2*pipf + 
										     lipf*((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*pow(lfpi,2) + 
											   (F10 + F20)*(F12*F21 - F11*F22)*m_l2*M2 + 
											   F10*lfpi*(F21*F22*lfpf + 
												     (2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - F21*F22*pipf)) + 
										     lipi*((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*pow(lfpi,2) + 
											   (F10 + F20)*(F12*F21 - F11*F22)*m_l2*M2 - 
											   lfpi*((F12*F20*F21 - (F11*F20 + F10*F21)*F22)*lfpf + 
												 F10*(2*F11 + F21)*(2*F12 + F22)*M2 + 
												 F10*F21*F22*pipf))) - 
	     kfli*(kfpf*(-((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lfpi,2)) + 
			 ((-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*
			  lfpf + 2*(F10 + F20)*(F11 + F21)*F22*m_l2)*M2 -
			 lfpi*(F10*F21*F22*lfpf + 
			       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
				3*F11*F20*F22 - F10*F21*F22 + 2*F20*F21*F22)*M2 +
			       F21*(-2*F12*F20 + (F10 - 2*F20)*F22)*pipf)) +
		   2*(-((F10 + F20)*(F12*F21 - F11*F22)*pow(lfpf,2)*M2) - 
		      pow(lfpi,2)*(F10*(F12*F21 + 3*(F11 + F21)*F22)*M2 - 
				   (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
		      lfpi*(M2*((-(F12*F20*F21) + F11*F20*F22 + 
				 F10*(2*F11 + F21)*(2*F12 + F22))*lfpf + 
				F10*(2*F11 + F21)*(2*F12 + F22)*M2) + 
			    ((F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*lfpf - 
			     2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf + 
			    F10*F21*F22*pow(pipf,2))))) + 
	  kfpi*(F20*F22*pow(kflf,2)*(-((2*F11 + F21)*lipf) + F21*lipi) - 
		2*kflf*(-2*F12*F20*F21*lfpf*lipf - 2*F11*F20*F22*lfpf*lipf + 
			F10*F21*F22*lfpf*lipf - 2*F20*F21*F22*lfpf*lipf - 
			F12*F20*F21*lfpi*lipf - 5*F11*F20*F22*lfpi*lipf - 
			4*F20*F21*F22*lfpi*lipf + 
			F22*kfpf*(-(F20*(2*F11 + F21)*lilf) + F10*F21*lipf + 
				  2*F20*F21*lipf - (2*F11*F20 + F10*F21)*lipi + 
				  2*F11*F20*m_l2 + F20*F21*m_l2) - 
			4*F11*F12*F20*lilf*M2 - 2*F12*F20*F21*lilf*M2 - 
			2*F10*F11*F22*lilf*M2 - 2*F11*F20*F22*lilf*M2 - 
			F10*F21*F22*lilf*M2 - F20*F21*F22*lilf*M2 + 
			16*F10*F11*F12*lipf*M2 + 4*F11*F12*F20*lipf*M2 + 
			12*F10*F12*F21*lipf*M2 + 3*F12*F20*F21*lipf*M2 + 
			2*F10*F11*F22*lipf*M2 + 3*F11*F20*F22*lipf*M2 + 
			2*F10*F21*F22*lipf*M2 + 2*F20*F21*F22*lipf*M2 + 
			8*F11*F12*F20*m_l2*M2 + 
			6*F12*F20*F21*m_l2*M2 + 
			6*F11*F20*F22*m_l2*M2 - 
			F10*F21*F22*m_l2*M2 + 
			5*F20*F21*F22*m_l2*M2 + 2*F12*F20*F21*lilf*pipf + 
			F10*F21*F22*lilf*pipf + F20*F21*F22*lilf*pipf + 
			2*F12*F20*F21*lipf*pipf - 2*F11*F20*F22*lipf*pipf + 
			4*F10*F21*F22*lipf*pipf - 2*F12*F20*F21*m_l2*pipf - 
			2*F11*F20*F22*m_l2*pipf - F10*F21*F22*m_l2*pipf - 
			3*F20*F21*F22*m_l2*pipf + 
			lipi*(-(F12*F20*F21*lfpf) + F11*F20*F22*lfpf + 
			      F10*F21*F22*lfpf + 2*F20*F21*F22*lfpf + 2*F20*F21*F22*lfpi - 
			      4*F11*F12*F20*M2 - 4*F10*F12*F21*M2 - 
			      5*F12*F20*F21*M2 - 2*F10*F11*F22*M2 - 
			      3*F11*F20*F22*M2 - 3*F10*F21*F22*M2 - 
			      4*F20*F21*F22*M2 + 
			      F21*(2*F12*F20 + 3*F10*F22 + 2*F20*F22)*pipf)) + 
		kfli*(-(F22*kflf*(2*F11*F20*lfpf + F20*F21*lfpf - F20*F21*lfpi + 
				  8*F10*F11*M2 + 8*F11*F20*M2 + 
				  6*F10*F21*M2 + 4*F20*F21*M2 + 
				  (4*F11*F20 - 2*F10*F21)*pipf)) + 
		      2*(-(F10*F21*F22*pow(lfpf,2)) - F20*F21*F22*pow(lfpf,2) - 
			 F20*F21*F22*pow(lfpi,2) + F12*F20*F21*lfpf*M2 + 
			 2*F10*F11*F22*lfpf*M2 - F11*F20*F22*lfpf*M2 + 
			 F10*F21*F22*lfpf*M2 + 
			 2*F10*F11*F22*m_l2*M2 + 
			 2*F10*F21*F22*m_l2*M2 - F10*F21*F22*lfpf*pipf + 
			 2*F11*F20*F22*m_l2*pipf + 
			 2*F20*F21*F22*m_l2*pipf + 
			 lfpi*((2*F11*F20 - F10*F21)*F22*lfpf - 
			       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
				3*F11*F20*F22 - 2*F10*F21*F22 + 2*F20*F21*F22)*
			       M2 + 2*F20*(F11 + F21)*F22*pipf))) + 
		2*(F20*F21*F22*pow(kfpf,2)*m_l2 + 
		   kfpf*(F10*F21*F22*lfpf*lipf - 2*F11*F20*F22*lfpi*lipf - 
			 F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf - 
			 F10*F21*F22*lfpf*m_l2 - 2*F11*F20*F22*lfpi*m_l2 + 
			 F10*F21*F22*lfpi*m_l2 - 2*F20*F21*F22*lfpi*m_l2 + 
			 F20*F21*F22*lipf*m_l2 + 
			 lipi*((2*F12*F20*F21 - 2*F11*F20*F22 - 3*F10*F21*F22)*lfpi + 
			       F21*(-((2*F12*F20 + F10*F22 + 2*F20*F22)*lfpf) + 
				    F20*F22*m_l2)) + 
			 8*F10*F11*F12*m_l2*M2 + 
			 8*F11*F12*F20*m_l2*M2 + 
			 4*F10*F12*F21*m_l2*M2 + 
			 6*F12*F20*F21*m_l2*M2 + 
			 16*F10*F11*F22*m_l2*M2 + 
			 6*F11*F20*F22*m_l2*M2 + 
			 12*F10*F21*F22*m_l2*M2 + 
			 4*F20*F21*F22*m_l2*M2 - 
			 lilf*(-(F10*F21*F22*lfpf) - F20*F21*F22*lfpf + 
			       (2*F11*F20 + (F10 + F20)*F21)*F22*lfpi + 
			       4*F11*F12*F20*M2 + 4*F12*F20*F21*M2 + 
			       8*F10*F11*F22*M2 + 2*F11*F20*F22*M2 + 
			       6*F10*F21*F22*M2 + 2*F20*F21*F22*M2 - 
			       2*F21*(F12*F20 + (F10 + F20)*F22)*pipf)) + 
		   2*(F10*F21*F22*pow(lfpf,2)*lipf - F12*F20*F21*lfpf*lfpi*lipf + 
		      F11*F20*F22*lfpf*lfpi*lipf + 2*F10*F21*F22*lfpf*lfpi*lipf + 
		      F10*F21*F22*pow(lfpi,2)*lipf - 
		      2*F10*F12*F21*lfpf*lipf*M2 - 
		      2*F10*F11*F22*lfpf*lipf*M2 - 
		      3*F10*F21*F22*lfpf*lipf*M2 + 
		      8*F10*F11*F12*lfpi*lipf*M2 + 
		      4*F10*F12*F21*lfpi*lipf*M2 + 
		      4*F10*F11*F22*lfpi*lipf*M2 + 
		      2*F10*F21*F22*lfpi*lipf*M2 - 
		      4*F10*F11*F12*lfpf*m_l2*M2 - 
		      4*F10*F12*F21*lfpf*m_l2*M2 - 
		      2*F10*F11*F22*lfpf*m_l2*M2 - 
		      F10*F21*F22*lfpf*m_l2*M2 - 
		      6*F10*F11*F22*lfpi*m_l2*M2 - 
		      5*F10*F21*F22*lfpi*m_l2*M2 + 
		      4*F10*F11*F12*lipf*m_l2*M2 + 
		      3*F10*F12*F21*lipf*m_l2*M2 - 
		      F10*F11*F22*lipf*m_l2*M2 - 
		      F10*F21*F22*lipf*m_l2*M2 + 
		      8*F10*F11*F12*m_l2*M4 + 
		      6*F10*F12*F21*m_l2*M4 + 
		      6*F10*F11*F22*m_l2*M4 + 
		      5*F10*F21*F22*m_l2*M4 + 
		      F10*F21*F22*lfpf*lipf*pipf + 2*F10*F21*F22*lfpi*lipf*pipf + 
		      F10*F21*F22*lfpf*m_l2*pipf + 
		      F10*F21*F22*lfpi*m_l2*pipf + 
		      F12*F20*F21*lipf*m_l2*pipf - 
		      F11*F20*F22*lipf*m_l2*pipf + 
		      F10*F21*F22*lipf*m_l2*pipf - 
		      4*F10*F11*F12*m_l2*M2*pipf - 
		      6*F10*F12*F21*m_l2*M2*pipf - 
		      6*F10*F11*F22*m_l2*M2*pipf - 
		      6*F10*F21*F22*m_l2*M2*pipf + 
		      F10*F21*F22*m_l2*pow(pipf,2) + 
		      lilf*(F10*M2*((4*F11*F12 + 3*F12*F21 + F11*F22)*lfpf - 
				    (2*F11 + F21)*(2*F12 + F22)*M2) + 
			    (F20*(F12*F21 - F11*F22)*lfpf + 
			     2*F10*(F12*F21 + (F11 + F21)*F22)*M2)*pipf - 
			    F10*F21*F22*pow(pipf,2) - 
			    (F12*F21 - F11*F22)*lfpi*(F10*M2 - F20*pipf)) + 
		      lipi*(-(F12*F20*F21*pow(lfpf,2)) + 
			    F11*F20*F22*pow(lfpf,2) + 4*F10*F11*F12*lfpf*M2 + 
			    2*F10*F12*F21*lfpf*M2 + 
			    2*F10*F11*F22*lfpf*M2 + 
			    F10*F21*F22*lfpf*M2 - 
			    F10*F12*F21*m_l2*M2 - 
			    F10*F11*F22*m_l2*M2 - 
			    F10*F21*F22*m_l2*M2 + 
			    (F10*F21*F22*lfpf + 
			     (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*m_l2)*
			    pipf - 2*lfpi*(F20*(F12*F21 - F11*F22)*lfpf + 
					   F10*(2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - 
					   F10*F21*F22*pipf))))))/(2.*kfpi*M2);

  res3 = -(pow(kfli,2)*(4*F20*(F11 + F21)*F22*kflf*(M2 - pipf) - 
			lfpi*(-(F20*(2*F11 + F21)*F22*kfpf) + F20*F21*F22*kfpi + 
			      2*((4*F11*F12*F20 + 3*F12*F20*F21 + 2*F10*F11*F22 + 
				  3*F11*F20*F22 + F10*F21*F22 + 2*F20*F21*F22)*M2 - 
				 (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)) + 
			lfpf*(-(F20*F21*F22*kfpf) + F20*(2*F11 + F21)*F22*kfpi + 
			      2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
				 (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
	   kfli*(4*F20*F21*F22*pow(kfpf,2)*lfpf - 
		 12*F11*F20*F22*pow(kfpf,2)*lfpi + 
		 2*F10*F21*F22*pow(kfpf,2)*lfpi - 
		 8*F20*F21*F22*pow(kfpf,2)*lfpi + 
		 F20*F21*F22*pow(kfpf,2)*lilf - 4*F20*F21*F22*kfpf*lfpf*lipf + 
		 2*F12*F20*F21*kfpf*lfpi*lipf + 10*F11*F20*F22*kfpf*lfpi*lipf + 
		 8*F20*F21*F22*kfpf*lfpi*lipf + 2*F12*F20*F21*kfpf*lfpf*lipi - 
		 2*F11*F20*F22*kfpf*lfpf*lipi - 2*F10*F21*F22*kfpf*lfpf*lipi - 
		 4*F20*F21*F22*kfpf*lfpf*lipi + 4*F12*F20*F21*kfpf*lfpi*lipi + 
		 4*F11*F20*F22*kfpf*lfpi*lipi - 2*F10*F21*F22*kfpf*lfpi*lipi + 
		 4*F20*F21*F22*kfpf*lfpi*lipi - 
		 F20*F21*F22*pow(kfpf,2)*m_l2 + 
		 F21*F22*pow(kfpi,2)*(-2*F10*lfpf + F20*lilf - F20*m_l2) + 
		 8*F11*F12*F20*kfpf*lfpf*M2 + 
		 8*F10*F12*F21*kfpf*lfpf*M2 + 
		 10*F12*F20*F21*kfpf*lfpf*M2 + 
		 4*F10*F11*F22*kfpf*lfpf*M2 + 
		 6*F11*F20*F22*kfpf*lfpf*M2 + 
		 6*F10*F21*F22*kfpf*lfpf*M2 + 
		 8*F20*F21*F22*kfpf*lfpf*M2 - 
		 32*F10*F11*F12*kfpf*lfpi*M2 - 
		 8*F11*F12*F20*kfpf*lfpi*M2 - 
		 24*F10*F12*F21*kfpf*lfpi*M2 - 
		 6*F12*F20*F21*kfpf*lfpi*M2 - 
		 4*F10*F11*F22*kfpf*lfpi*M2 - 
		 6*F11*F20*F22*kfpf*lfpi*M2 - 
		 4*F10*F21*F22*kfpf*lfpi*M2 - 
		 4*F20*F21*F22*kfpf*lfpi*M2 + 
		 8*F11*F12*F20*kfpf*lilf*M2 + 
		 4*F12*F20*F21*kfpf*lilf*M2 + 
		 4*F10*F11*F22*kfpf*lilf*M2 + 
		 4*F11*F20*F22*kfpf*lilf*M2 + 
		 2*F10*F21*F22*kfpf*lilf*M2 + 
		 2*F20*F21*F22*kfpf*lilf*M2 - 
		 12*F10*F12*F21*lfpf*lipf*M2 - 
		 20*F10*F11*F22*lfpf*lipf*M2 - 
		 24*F10*F21*F22*lfpf*lipf*M2 + 
		 32*F10*F11*F12*lfpi*lipf*M2 + 
		 20*F10*F12*F21*lfpi*lipf*M2 + 
		 12*F10*F11*F22*lfpi*lipf*M2 + 
		 8*F10*F21*F22*lfpi*lipf*M2 + 
		 16*F10*F11*F12*lfpf*lipi*M2 + 
		 4*F10*F12*F21*lfpf*lipi*M2 - 
		 4*F12*F20*F21*lfpf*lipi*M2 + 
		 12*F10*F11*F22*lfpf*lipi*M2 + 
		 4*F11*F20*F22*lfpf*lipi*M2 + 
		 4*F10*F21*F22*lfpf*lipi*M2 - 
		 12*F10*F12*F21*lfpi*lipi*M2 - 
		 4*F12*F20*F21*lfpi*lipi*M2 - 
		 4*F10*F11*F22*lfpi*lipi*M2 + 
		 4*F11*F20*F22*lfpi*lipi*M2 - 
		 12*F10*F21*F22*lfpi*lipi*M2 - 
		 16*F11*F12*F20*kfpf*m_l2*M2 - 
		 12*F12*F20*F21*kfpf*m_l2*M2 - 
		 12*F11*F20*F22*kfpf*m_l2*M2 + 
		 2*F10*F21*F22*kfpf*m_l2*M2 - 
		 10*F20*F21*F22*kfpf*m_l2*M2 - 
		 16*F10*F11*F12*lfpf*M4 - 16*F10*F12*F21*lfpf*M4 - 
		 16*F10*F11*F22*lfpf*M4 - 16*F10*F21*F22*lfpf*M4 + 
		 16*F10*F11*F12*lfpi*M4 + 8*F10*F12*F21*lfpi*M4 + 
		 8*F10*F11*F22*lfpi*M4 + 4*F10*F21*F22*lfpi*M4 - 
		 16*F10*F11*F12*lilf*M4 - 8*F10*F12*F21*lilf*M4 - 
		 8*F10*F11*F22*lilf*M4 - 4*F10*F21*F22*lilf*M4 + 
		 32*F10*F11*F12*m_l2*M4 + 
		 24*F10*F12*F21*m_l2*M4 + 
		 24*F10*F11*F22*m_l2*M4 + 
		 20*F10*F21*F22*m_l2*M4 - 
		 4*F12*F20*F21*kfpf*lfpf*pipf - 6*F10*F21*F22*kfpf*lfpf*pipf - 
		 4*F20*F21*F22*kfpf*lfpf*pipf - 4*F12*F20*F21*kfpf*lfpi*pipf + 
		 4*F11*F20*F22*kfpf*lfpi*pipf - 8*F10*F21*F22*kfpf*lfpi*pipf - 
		 4*F12*F20*F21*kfpf*lilf*pipf - 2*F10*F21*F22*kfpf*lilf*pipf - 
		 2*F20*F21*F22*kfpf*lilf*pipf + 4*F12*F20*F21*lfpf*lipf*pipf - 
		 4*F11*F20*F22*lfpf*lipf*pipf + 8*F10*F21*F22*lfpf*lipf*pipf + 
		 4*F12*F20*F21*lfpi*lipf*pipf - 4*F11*F20*F22*lfpi*lipf*pipf + 
		 8*F10*F21*F22*lfpi*lipf*pipf + 4*F10*F21*F22*lfpf*lipi*pipf + 
		 4*F10*F21*F22*lfpi*lipi*pipf + 
		 4*F12*F20*F21*kfpf*m_l2*pipf + 
		 4*F11*F20*F22*kfpf*m_l2*pipf + 
		 2*F10*F21*F22*kfpf*m_l2*pipf + 
		 6*F20*F21*F22*kfpf*m_l2*pipf + 
		 16*F10*F11*F12*lfpf*M2*pipf + 
		 16*F10*F12*F21*lfpf*M2*pipf + 
		 16*F10*F11*F22*lfpf*M2*pipf + 
		 16*F10*F21*F22*lfpf*M2*pipf - 
		 8*F10*F12*F21*lfpi*M2*pipf - 
		 8*F10*F11*F22*lfpi*M2*pipf - 
		 8*F10*F21*F22*lfpi*M2*pipf + 
		 8*F10*F12*F21*lilf*M2*pipf + 
		 8*F10*F11*F22*lilf*M2*pipf + 
		 8*F10*F21*F22*lilf*M2*pipf - 
		 16*F10*F11*F12*m_l2*M2*pipf - 
		 24*F10*F12*F21*m_l2*M2*pipf - 
		 24*F10*F11*F22*m_l2*M2*pipf - 
		 24*F10*F21*F22*m_l2*M2*pipf + 
		 4*F10*F21*F22*lfpi*pow(pipf,2) - 
		 4*F10*F21*F22*lilf*pow(pipf,2) + 
		 4*F10*F21*F22*m_l2*pow(pipf,2) + 
		 kfpi*(-(F22*kflf*(-(F20*(2*F11 + F21)*lipf) + F20*F21*lipi + 
				   8*F10*F11*M2 + 4*F11*F20*M2 + 
				   6*F10*F21*M2 + 4*F20*F21*M2 - 
				   2*F10*F21*pipf)) + 
		       2*(-2*F12*F20*F21*lfpf*lipf + 4*F11*F20*F22*lfpf*lipf + 
			  F10*F21*F22*lfpf*lipf + F20*F21*F22*lfpf*lipf - 
			  F12*F20*F21*lfpi*lipf + F11*F20*F22*lfpi*lipf + 
			  F10*F21*F22*lfpi*lipf - F20*F21*F22*lfpi*lipf - 
			  F12*F20*F21*lfpf*lipi + F11*F20*F22*lfpf*lipi - 
			  F20*F21*F22*lfpf*lipi - F20*F21*F22*lfpi*lipi + 
			  F22*kfpf*(-2*F11*F20*lfpf - F10*F21*lfpf + 
				    (F10 + 2*F20)*F21*lfpi - F20*(2*F11 + F21)*lilf + 
				    2*F11*F20*m_l2 + F20*F21*m_l2) + 
			  4*F10*F12*F21*lfpf*M2 + 
			  F12*F20*F21*lfpf*M2 + 
			  2*F10*F11*F22*lfpf*M2 - 
			  F11*F20*F22*lfpf*M2 + 
			  4*F10*F21*F22*lfpf*M2 + 
			  4*F10*F12*F21*lfpi*M2 + 
			  3*F12*F20*F21*lfpi*M2 - 
			  2*F10*F11*F22*lfpi*M2 - 
			  F11*F20*F22*lfpi*M2 + 
			  3*F10*F21*F22*lfpi*M2 + 
			  2*F20*F21*F22*lfpi*M2 + 
			  4*F11*F12*F20*m_l2*M2 + 
			  4*F12*F20*F21*m_l2*M2 - 
			  4*F10*F11*F22*m_l2*M2 + 
			  4*F11*F20*F22*m_l2*M2 - 
			  3*F10*F21*F22*m_l2*M2 + 
			  3*F20*F21*F22*m_l2*M2 - 
			  2*F12*F20*F21*lfpf*pipf - 2*F10*F21*F22*lfpf*pipf - 
			  2*F20*F21*F22*lfpf*pipf + F10*F21*F22*lfpi*pipf + 
			  F10*F21*F22*m_l2*pipf - F20*F21*F22*m_l2*pipf + 
			  (F10 - F20)*F22*lilf*((2*F11 + F21)*M2 - F21*pipf))) + 
		 kflf*(F22*kfpf*(-(F20*F21*lipf) + F20*(2*F11 + F21)*lipi + 
				 8*F10*F11*M2 + 8*F11*F20*M2 + 
				 6*F10*F21*M2 + 4*F20*F21*M2 + 
				 4*F11*F20*pipf - 2*F10*F21*pipf) + 
		       2*(2*((2*F11 + F21)*M2 - F21*pipf)*
			  ((F20*(F12 + F22) + F10*(2*F12 + F22))*M2 - 
			   (F12*F20 + (F10 + F20)*F22)*pipf) + 
			  lipi*((F11*(-4*F12*F20 + 2*F10*F22 + F20*F22) + 
				 F21*(-3*F12*F20 + 3*F10*F22 + 2*F20*F22))*M2 + 
				(F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
			  lipf*((F12*F20*F21 + 
				 (-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + 
				(F12*F20*F21 + (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) - 
	   2*(F20*F21*F22*pow(kfpf,3)*m_l2 + 
	      pow(kfpi,2)*(-2*F12*F20*F21*lfpf*lipf + F10*F21*F22*lfpf*lipf - 
			   2*F20*F21*F22*lfpf*lipf + F10*F21*F22*lfpi*lipf + 
			   F20*F21*F22*kfpf*m_l2 - F10*F21*F22*lipf*m_l2 + 
			   2*F10*F12*F21*m_l2*M2 + 
			   F12*F20*F21*m_l2*M2 + 
			   2*F10*F11*F22*m_l2*M2 + 
			   F11*F20*F22*m_l2*M2 + 
			   4*F10*F21*F22*m_l2*M2 + 
			   2*F20*F21*F22*m_l2*M2 + 
			   lilf*(F10*F21*F22*lipf + 
				 (F12*F20*F21 - (F11*F20 + 2*F10*(F11 + F21))*F22)*M2)) + 
	      pow(kfpf,2)*(-2*F10*F21*F22*lfpi*lipf - 
			   F20*F21*F22*lfpf*m_l2 + 4*F11*F20*F22*lfpi*m_l2 + 
			   3*F20*F21*F22*lfpi*m_l2 - 
			   lipi*(-2*F12*F20*F21*lfpf + 2*F11*F20*F22*lfpf + 
				 F10*F21*F22*lfpf + 
				 (2*F11*F20 + 3*F10*F21 + 2*F20*F21)*F22*lfpi + 
				 (-2*F11*F20 + (F10 - F20)*F21)*F22*lilf + 
				 2*F11*F20*F22*m_l2 - F10*F21*F22*m_l2 + 
				 2*F20*F21*F22*m_l2) + 
			   8*F11*F12*F20*m_l2*M2 + 
			   2*F10*F12*F21*m_l2*M2 + 
			   7*F12*F20*F21*m_l2*M2 + 
			   14*F10*F11*F22*m_l2*M2 + 
			   7*F11*F20*F22*m_l2*M2 + 
			   12*F10*F21*F22*m_l2*M2 + 
			   6*F20*F21*F22*m_l2*M2 - 
			   2*F12*F20*F21*m_l2*pipf - 2*F11*F20*F22*m_l2*pipf - 
			   4*F10*F21*F22*m_l2*pipf - 4*F20*F21*F22*m_l2*pipf - 
			   lilf*(F20*F21*F22*lipf + 
				 (F21*(3*F12*F20 + 4*F10*F22 + 2*F20*F22) + 
				  F11*(4*F12*F20 + 6*F10*F22 + 3*F20*F22))*M2 - 
				 2*(F11*F20 + (F10 + F20)*F21)*F22*pipf)) - 
	      2*(2*F10*lipf*(-2*F12*F21*lfpf*lipf*M2 - 
			     2*F11*F22*lfpf*lipf*M2 - 
			     3*F21*F22*lfpf*lipf*M2 + 
			     4*F11*F12*lfpi*lipf*M2 + 
			     2*F12*F21*lfpi*lipf*M2 + 
			     2*F11*F22*lfpi*lipf*M2 + 
			     F21*F22*lfpi*lipf*M2 + 
			     8*F11*F12*m_l2*M4 + 
			     6*F12*F21*m_l2*M4 + 
			     6*F11*F22*m_l2*M4 + 
			     5*F21*F22*m_l2*M4 + F21*F22*lfpf*lipf*pipf + 
			     F21*F22*lfpi*lipf*pipf - 
			     4*F11*F12*m_l2*M2*pipf - 
			     6*F12*F21*m_l2*M2*pipf - 
			     6*F11*F22*m_l2*M2*pipf - 
			     6*F21*F22*m_l2*M2*pipf + 
			     F21*F22*m_l2*pow(pipf,2) - 
			     lilf*((2*F11 + F21)*M2 - F21*pipf)*
			     ((2*F12 + F22)*M2 - F22*pipf) + 
			     lipi*(lfpf*((2*F11 + F21)*(2*F12 + F22)*M2 + 
					 F21*F22*pipf) + 
				   lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*
					   M2) + F21*F22*pipf))) + 
		 kflf*(-((F10 + F20)*(F12*F21 - F11*F22)*pow(lipi,2)*M2) + 
		       lipf*lipi*((-(F12*F20*F21) + F11*F20*F22 + 
				   F10*(2*F11 + F21)*(2*F12 + F22))*M2 + 
				  (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) + 
		       lipf*(F10*((2*F11 + F21)*M2 - F21*pipf)*
			     ((2*F12 + F22)*M2 - F22*pipf) - 
			     lipf*(F10*(F12*F21 + 3*(F11 + F21)*F22)*M2 - 
				   (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf)))) + 
	      kfpf*(kflf*(F20*F21*F22*pow(lipf,2) + 
			  (F10 + F20)*F21*F22*pow(lipi,2) - 
			  2*(F11 + F21)*F22*m_l2*(F10*M2 + F20*pipf) + 
			  lipi*((-2*F11*F20*F22 + F10*F21*F22)*lipf - 
				(F12*F20*F21 + 2*F10*F11*F22 - F11*F20*F22 + 
				 F10*F21*F22)*M2 + F10*F21*F22*pipf) + 
			  lipf*((F21*(3*F12*F20 - 2*F10*F22 + 2*F20*F22) + 
				 F11*(4*F12*F20 - 2*F10*F22 + 3*F20*F22))*M2 - 
				2*F20*(F11 + F21)*F22*pipf)) + 
		    2*(F10*F21*F22*lfpi*pow(lipf,2) + 
		       (F20*(-(F12*F21) + F11*F22)*lfpf + F10*F21*F22*lfpi)*
		       pow(lipi,2) - 4*F10*F12*F21*lfpf*lipf*M2 - 
		       4*F10*F11*F22*lfpf*lipf*M2 - 
		       6*F10*F21*F22*lfpf*lipf*M2 + 
		       8*F10*F11*F12*lfpi*lipf*M2 + 
		       4*F10*F12*F21*lfpi*lipf*M2 + 
		       4*F10*F11*F22*lfpi*lipf*M2 + 
		       2*F10*F21*F22*lfpi*lipf*M2 - 
		       F10*F12*F21*lfpf*m_l2*M2 - 
		       F10*F11*F22*lfpf*m_l2*M2 - 
		       F10*F21*F22*lfpf*m_l2*M2 + 
		       4*F10*F11*F12*lfpi*m_l2*M2 + 
		       3*F10*F12*F21*lfpi*m_l2*M2 - 
		       F10*F11*F22*lfpi*m_l2*M2 - 
		       F10*F21*F22*lfpi*m_l2*M2 - 
		       6*F10*F11*F22*lipf*m_l2*M2 - 
		       5*F10*F21*F22*lipf*m_l2*M2 + 
		       8*F10*F11*F12*m_l2*M4 + 
		       6*F10*F12*F21*m_l2*M4 + 
		       6*F10*F11*F22*m_l2*M4 + 
		       5*F10*F21*F22*m_l2*M4 + 
		       2*F10*F21*F22*lfpf*lipf*pipf + 
		       2*F10*F21*F22*lfpi*lipf*pipf + 
		       F12*F20*F21*lfpf*m_l2*pipf - 
		       F11*F20*F22*lfpf*m_l2*pipf + 
		       F10*F21*F22*lfpf*m_l2*pipf + 
		       F12*F20*F21*lfpi*m_l2*pipf - 
		       F11*F20*F22*lfpi*m_l2*pipf + 
		       F10*F21*F22*lfpi*m_l2*pipf + 
		       F10*F21*F22*lipf*m_l2*pipf - 
		       4*F10*F11*F12*m_l2*M2*pipf - 
		       6*F10*F12*F21*m_l2*M2*pipf - 
		       6*F10*F11*F22*m_l2*M2*pipf - 
		       6*F10*F21*F22*m_l2*M2*pipf + 
		       F10*F21*F22*m_l2*pow(pipf,2) - 
		       lilf*((F12*F21 - F11*F22)*lipf*(F10*M2 - F20*pipf) + 
			     F10*((2*F11 + F21)*M2 - F21*pipf)*
			     ((2*F12 + F22)*M2 - F22*pipf)) + 
		       lipi*(-2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf + 
			     4*F10*F11*F12*lfpf*M2 + 
			     2*F10*F12*F21*lfpf*M2 + 
			     2*F10*F11*F22*lfpf*M2 + 
			     F10*F21*F22*lfpf*M2 - 
			     4*F10*F11*F12*m_l2*M2 - 
			     4*F10*F12*F21*m_l2*M2 - 
			     2*F10*F11*F22*m_l2*M2 - 
			     F10*F21*F22*m_l2*M2 + 
			     F10*F21*F22*lfpf*pipf + F10*F21*F22*m_l2*pipf + 
			     lilf*(F10*(4*F11*F12 + 3*F12*F21 + F11*F22)*M2 + 
				   F20*(F12*F21 - F11*F22)*pipf) + 
			     lfpi*((-(F12*F20*F21) + F11*F20*F22 + 2*F10*F21*F22)*lipf + 
				   F10*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
					F21*F22*pipf))))) + 
	      kfpi*(-2*F20*(2*F11 + F21)*F22*pow(kfpf,2)*m_l2 + 
		    kflf*(-((2*F11*F20 + (F10 + 2*F20)*F21)*F22*pow(lipf,2)) + 
			  2*(F10 + F20)*(F11 + F21)*F22*m_l2*M2 + 
			  lipi*(-(F10*F21*F22*lipf) + 
				(-(F12*F20*F21) + (F11*F20 + 2*F10*(F11 + F21))*F22)*M2) + 
			  lipf*((F11*(-4*F12*F20 + 2*F10*F22 - 3*F20*F22) + 
				 F21*(-3*F12*F20 + F10*F22 - 2*F20*F22))*M2 + 
				F21*(2*F12*F20 - F10*F22 + 2*F20*F22)*pipf)) + 
		    kfpf*(-2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf + 
			  3*F10*F21*F22*lfpf*lipf + 2*F11*F20*F22*lfpi*lipf + 
			  F10*F21*F22*lfpi*lipf + 2*F20*F21*F22*lfpi*lipf - 
			  F20*F21*F22*lfpf*m_l2 - F20*F21*F22*lfpi*m_l2 + 
			  2*F11*F20*F22*lipf*m_l2 - 
			  F10*F21*F22*lipf*m_l2 + 
			  2*F20*F21*F22*lipf*m_l2 + 
			  F21*lipi*(2*F12*F20*lfpf + F10*F22*lfpf + 2*F20*F22*lfpf - 
				    F10*F22*lfpi - (F10 + F20)*F22*lilf + F10*F22*m_l2) - 
			  8*F10*F11*F12*m_l2*M2 - 
			  8*F11*F12*F20*m_l2*M2 - 
			  4*F10*F12*F21*m_l2*M2 - 
			  6*F12*F20*F21*m_l2*M2 - 
			  16*F10*F11*F22*m_l2*M2 - 
			  6*F11*F20*F22*m_l2*M2 - 
			  12*F10*F21*F22*m_l2*M2 - 
			  4*F20*F21*F22*m_l2*M2 + 
			  lilf*((2*F11*F20 + (F10 + F20)*F21)*F22*lipf + 
				2*((F21*(2*F12*F20 + 3*F10*F22 + F20*F22) + 
				    F11*(2*F12*F20 + 4*F10*F22 + F20*F22))*M2 - 
				   F21*(F12*F20 + (F10 + F20)*F22)*pipf))) - 
		    2*(((-2*F12*F20*F21 + 2*F11*F20*F22 + F10*F21*F22)*lfpf + 
			(-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpi)*pow(lipf,2) + 
		       (F10 + F20)*(F12*F21 - F11*F22)*(lfpf + lfpi)*m_l2*M2 + 
		       lipi*(((-(F12*F20*F21) + F11*F20*F22 + F10*F21*F22)*lfpf + 
			      F10*F21*F22*lfpi)*lipf + 
			     (F10 + F20)*(F12*F21 - F11*F22)*lilf*M2) - 
		       lipf*(-(lilf*((F10*F12*F21 + F12*F20*F21 + 3*F10*F11*F22 - 
				      F11*F20*F22 + 2*F10*F21*F22)*M2 - 
				     2*F10*F21*F22*pipf)) + 
			     F10*(((2*F11 + F21)*(2*F12 + F22)*lfpf + 
				   2*(2*F11*F12 + 2*F12*F21 + 4*F11*F22 + 3*F21*F22)*m_l2)*M2 + 
				  F21*F22*(lfpf - 2*m_l2)*pipf + 
				  lfpi*(-((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2) + 
					F21*F22*pipf)))))))/(2.*kfpf*M2);

  res4 = (pow(kfli,2)*(4*F20*(F11 + F21)*F22*kflf*(M2 - pipf) + 
		       lfpf*(-(F20*F21*F22*kfpf) + F20*(2*F11 + F21)*F22*kfpi + 
			     8*F11*F12*F20*M2 + 6*F12*F20*F21*M2 + 
			     4*F10*F11*F22*M2 + 6*F11*F20*F22*M2 + 
			     2*F10*F21*F22*M2 + 4*F20*F21*F22*M2 - 
			     2*F12*F20*F21*pipf + 2*F11*F20*F22*pipf - 2*F10*F21*F22*pipf) + 
		       lfpi*(F20*(2*F11 + F21)*F22*kfpf - F20*F21*F22*kfpi - 
			     2*((F12*F20*F21 + (-(F10*(2*F11 + F21)) + F20*(F11 + 2*F21))*F22)*M2 + 
				(F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf))) + 
	  kfli*(2*F10*F21*F22*pow(kfpf,2)*lfpi + 
		F20*F21*F22*pow(kfpf,2)*lilf + 2*F20*F21*F22*kfpf*lfpf*lipf + 
		2*F12*F20*F21*kfpf*lfpi*lipf - 2*F11*F20*F22*kfpf*lfpi*lipf + 
		2*F20*F21*F22*kfpf*lfpi*lipf + 2*F12*F20*F21*kfpf*lfpf*lipi - 
		2*F11*F20*F22*kfpf*lfpf*lipi - 2*F10*F21*F22*kfpf*lfpf*lipi + 
		2*F20*F21*F22*kfpf*lfpf*lipi + 4*F12*F20*F21*kfpf*lfpi*lipi - 
		8*F11*F20*F22*kfpf*lfpi*lipi - 2*F10*F21*F22*kfpf*lfpi*lipi - 
		2*F20*F21*F22*kfpf*lfpi*lipi - 
		F20*F21*F22*pow(kfpf,2)*m_l2 + 
		F22*pow(kfpi,2)*(12*F11*F20*lfpf - 2*F10*F21*lfpf + 
				 8*F20*F21*lfpf - 4*F20*F21*lfpi + F20*F21*lilf - 
				 F20*F21*m_l2) + 8*F10*F12*F21*kfpf*lfpf*M2 + 
		6*F12*F20*F21*kfpf*lfpf*M2 - 
		4*F10*F11*F22*kfpf*lfpf*M2 - 
		2*F11*F20*F22*kfpf*lfpf*M2 + 
		6*F10*F21*F22*kfpf*lfpf*M2 + 
		4*F20*F21*F22*kfpf*lfpf*M2 + 
		8*F10*F12*F21*kfpf*lfpi*M2 + 
		2*F12*F20*F21*kfpf*lfpi*M2 + 
		4*F10*F11*F22*kfpf*lfpi*M2 - 
		2*F11*F20*F22*kfpf*lfpi*M2 + 
		8*F10*F21*F22*kfpf*lfpi*M2 - 
		4*F10*F11*F22*kfpf*lilf*M2 + 
		4*F11*F20*F22*kfpf*lilf*M2 - 
		2*F10*F21*F22*kfpf*lilf*M2 + 
		2*F20*F21*F22*kfpf*lilf*M2 - 
		12*F10*F12*F21*lfpf*lipf*M2 - 
		4*F12*F20*F21*lfpf*lipf*M2 - 
		4*F10*F11*F22*lfpf*lipf*M2 + 
		4*F11*F20*F22*lfpf*lipf*M2 - 
		12*F10*F21*F22*lfpf*lipf*M2 + 
		16*F10*F11*F12*lfpi*lipf*M2 + 
		4*F10*F12*F21*lfpi*lipf*M2 - 
		4*F12*F20*F21*lfpi*lipf*M2 + 
		12*F10*F11*F22*lfpi*lipf*M2 + 
		4*F11*F20*F22*lfpi*lipf*M2 + 
		4*F10*F21*F22*lfpi*lipf*M2 + 
		32*F10*F11*F12*lfpf*lipi*M2 + 
		20*F10*F12*F21*lfpf*lipi*M2 + 
		12*F10*F11*F22*lfpf*lipi*M2 + 
		8*F10*F21*F22*lfpf*lipi*M2 - 
		12*F10*F12*F21*lfpi*lipi*M2 - 
		20*F10*F11*F22*lfpi*lipi*M2 - 
		24*F10*F21*F22*lfpi*lipi*M2 - 
		8*F11*F12*F20*kfpf*m_l2*M2 - 
		8*F12*F20*F21*kfpf*m_l2*M2 + 
		8*F10*F11*F22*kfpf*m_l2*M2 - 
		8*F11*F20*F22*kfpf*m_l2*M2 + 
		6*F10*F21*F22*kfpf*m_l2*M2 - 
		6*F20*F21*F22*kfpf*m_l2*M2 - 
		16*F10*F11*F12*lfpf*M4 - 8*F10*F12*F21*lfpf*M4 - 
		8*F10*F11*F22*lfpf*M4 - 4*F10*F21*F22*lfpf*M4 + 
		16*F10*F11*F12*lfpi*M4 + 16*F10*F12*F21*lfpi*M4 + 
		16*F10*F11*F22*lfpi*M4 + 16*F10*F21*F22*lfpi*M4 - 
		16*F10*F11*F12*lilf*M4 - 8*F10*F12*F21*lilf*M4 - 
		8*F10*F11*F22*lilf*M4 - 4*F10*F21*F22*lilf*M4 + 
		32*F10*F11*F12*m_l2*M4 + 
		24*F10*F12*F21*m_l2*M4 + 
		24*F10*F11*F22*m_l2*M4 + 
		20*F10*F21*F22*m_l2*M4 + 
		2*F10*F21*F22*kfpf*lfpf*pipf - 4*F12*F20*F21*kfpf*lfpi*pipf - 
		4*F10*F21*F22*kfpf*lfpi*pipf - 4*F20*F21*F22*kfpf*lfpi*pipf + 
		2*F10*F21*F22*kfpf*lilf*pipf - 2*F20*F21*F22*kfpf*lilf*pipf + 
		4*F10*F21*F22*lfpf*lipf*pipf + 4*F10*F21*F22*lfpi*lipf*pipf + 
		4*F12*F20*F21*lfpf*lipi*pipf - 4*F11*F20*F22*lfpf*lipi*pipf + 
		8*F10*F21*F22*lfpf*lipi*pipf + 4*F12*F20*F21*lfpi*lipi*pipf - 
		4*F11*F20*F22*lfpi*lipi*pipf + 8*F10*F21*F22*lfpi*lipi*pipf - 
		2*F10*F21*F22*kfpf*m_l2*pipf + 
		2*F20*F21*F22*kfpf*m_l2*pipf + 
		8*F10*F12*F21*lfpf*M2*pipf + 
		8*F10*F11*F22*lfpf*M2*pipf + 
		8*F10*F21*F22*lfpf*M2*pipf - 
		16*F10*F11*F12*lfpi*M2*pipf - 
		16*F10*F12*F21*lfpi*M2*pipf - 
		16*F10*F11*F22*lfpi*M2*pipf - 
		16*F10*F21*F22*lfpi*M2*pipf + 
		8*F10*F12*F21*lilf*M2*pipf + 
		8*F10*F11*F22*lilf*M2*pipf + 
		8*F10*F21*F22*lilf*M2*pipf - 
		16*F10*F11*F12*m_l2*M2*pipf - 
		24*F10*F12*F21*m_l2*M2*pipf - 
		24*F10*F11*F22*m_l2*M2*pipf - 
		24*F10*F21*F22*m_l2*M2*pipf - 
		4*F10*F21*F22*lfpf*pow(pipf,2) - 
		4*F10*F21*F22*lilf*pow(pipf,2) + 
		4*F10*F21*F22*m_l2*pow(pipf,2) + 
		kfpi*(-(F22*kflf*(-(F20*(2*F11 + F21)*lipf) + F20*F21*lipi + 
				  8*F10*F11*M2 + 8*F11*F20*M2 + 
				  6*F10*F21*M2 + 4*F20*F21*M2 + 
				  4*F11*F20*pipf - 2*F10*F21*pipf)) - 
		      2*(2*F12*F20*F21*lfpf*lipf + 2*F11*F20*F22*lfpf*lipf - 
			 F10*F21*F22*lfpf*lipf + 2*F20*F21*F22*lfpf*lipf + 
			 F12*F20*F21*lfpi*lipf - F11*F20*F22*lfpi*lipf - 
			 F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf + 
			 F12*F20*F21*lfpf*lipi + 5*F11*F20*F22*lfpf*lipi + 
			 4*F20*F21*F22*lfpf*lipi - 2*F20*F21*F22*lfpi*lipi + 
			 F22*kfpf*(F10*F21*lfpf + 2*F20*F21*lfpf - 
				   (2*F11*F20 + F10*F21)*lfpi + F20*(2*F11 + F21)*lilf - 
				   2*F11*F20*m_l2 - F20*F21*m_l2) + 
			 16*F10*F11*F12*lfpf*M2 + 
			 4*F11*F12*F20*lfpf*M2 + 
			 12*F10*F12*F21*lfpf*M2 + 
			 3*F12*F20*F21*lfpf*M2 + 
			 2*F10*F11*F22*lfpf*M2 + 
			 3*F11*F20*F22*lfpf*M2 + 
			 2*F10*F21*F22*lfpf*M2 + 
			 2*F20*F21*F22*lfpf*M2 - 
			 4*F11*F12*F20*lfpi*M2 - 
			 4*F10*F12*F21*lfpi*M2 - 
			 5*F12*F20*F21*lfpi*M2 - 
			 2*F10*F11*F22*lfpi*M2 - 
			 3*F11*F20*F22*lfpi*M2 - 
			 3*F10*F21*F22*lfpi*M2 - 
			 4*F20*F21*F22*lfpi*M2 - 
			 8*F11*F12*F20*m_l2*M2 - 
			 6*F12*F20*F21*m_l2*M2 - 
			 6*F11*F20*F22*m_l2*M2 + 
			 F10*F21*F22*m_l2*M2 - 
			 5*F20*F21*F22*m_l2*M2 + 
			 2*F12*F20*F21*lfpf*pipf - 2*F11*F20*F22*lfpf*pipf + 
			 4*F10*F21*F22*lfpf*pipf + 2*F12*F20*F21*lfpi*pipf + 
			 3*F10*F21*F22*lfpi*pipf + 2*F20*F21*F22*lfpi*pipf + 
			 2*F12*F20*F21*m_l2*pipf + 
			 2*F11*F20*F22*m_l2*pipf + 
			 F10*F21*F22*m_l2*pipf + 
			 3*F20*F21*F22*m_l2*pipf + 
			 (2*F12*F20 + (F10 + F20)*F22)*lilf*((2*F11 + F21)*M2 - F21*pipf))) + 
		kflf*(F22*kfpf*(-(F20*F21*lipf) + F20*(2*F11 + F21)*lipi + 
				8*F10*F11*M2 + 4*F11*F20*M2 + 
				6*F10*F21*M2 + 4*F20*F21*M2 - 2*F10*F21*pipf) + 
		      2*(2*((2*F11 + F21)*M2 - F21*pipf)*
			 ((F20*(F12 + F22) + F10*(2*F12 + F22))*M2 - 
			  (F12*F20 + (F10 + F20)*F22)*pipf) + 
			 lipf*((F21*(3*F12*F20 - 3*F10*F22 - 2*F20*F22) + 
				F11*(4*F12*F20 - (2*F10 + F20)*F22))*M2 - 
			       (F12*F20*F21 - F11*F20*F22 + F10*F21*F22)*pipf) - lipi*((F12*F20*F21 + 
											(-6*F10*F11 + F11*F20 - 5*F10*F21 + 2*F20*F21)*F22)*M2 + (F12*F20*F21 + 
																		  (-5*F11*F20 + F10*F21 - 4*F20*F21)*F22)*pipf)))) + 
	  2*(-(F10*F21*F22*pow(kfpf,2)*lfpf*lipi) + 
	     2*F12*F20*F21*pow(kfpf,2)*lfpi*lipi - 
	     F10*F21*F22*pow(kfpf,2)*lfpi*lipi + 
	     2*F20*F21*F22*pow(kfpf,2)*lfpi*lipi + 
	     F10*F21*F22*pow(kfpf,2)*lilf*lipi - 
	     F10*F21*F22*kflf*kfpf*lipf*lipi + 
	     2*F10*F21*F22*kfpf*lfpf*lipf*lipi - 
	     2*F12*F20*F21*kfpf*lfpi*lipf*lipi + 
	     2*F11*F20*F22*kfpf*lfpi*lipf*lipi + 
	     2*F10*F21*F22*kfpf*lfpi*lipf*lipi - 
	     2*F11*F20*F22*kflf*kfpf*pow(lipi,2) - 
	     F10*F21*F22*kflf*kfpf*pow(lipi,2) - 
	     2*F20*F21*F22*kflf*kfpf*pow(lipi,2) - 
	     2*F12*F20*F21*kfpf*lfpf*pow(lipi,2) + 
	     2*F11*F20*F22*kfpf*lfpf*pow(lipi,2) + 
	     2*F10*F21*F22*kfpf*lfpf*pow(lipi,2) - 
	     4*F12*F20*F21*kfpf*lfpi*pow(lipi,2) + 
	     4*F11*F20*F22*kfpf*lfpi*pow(lipi,2) + 
	     2*F10*F21*F22*kfpf*lfpi*pow(lipi,2) + 
	     F20*F21*F22*pow(kfpi,3)*m_l2 - 
	     F10*F21*F22*pow(kfpf,2)*lipi*m_l2 - 
	     F12*F20*F21*pow(kfpf,2)*lilf*M2 + 
	     2*F10*F11*F22*pow(kfpf,2)*lilf*M2 + 
	     F11*F20*F22*pow(kfpf,2)*lilf*M2 + 
	     2*F10*F21*F22*pow(kfpf,2)*lilf*M2 + 
	     F12*F20*F21*kflf*kfpf*lipf*M2 - 
	     2*F10*F11*F22*kflf*kfpf*lipf*M2 - 
	     F11*F20*F22*kflf*kfpf*lipf*M2 - 
	     2*F10*F21*F22*kflf*kfpf*lipf*M2 + 
	     2*F10*F12*F21*kfpf*lilf*lipf*M2 + 
	     2*F12*F20*F21*kfpf*lilf*lipf*M2 - 
	     2*F10*F11*F22*kfpf*lilf*lipf*M2 - 
	     2*F11*F20*F22*kfpf*lilf*lipf*M2 - 
	     2*F10*F12*F21*kflf*pow(lipf,2)*M2 - 
	     2*F12*F20*F21*kflf*pow(lipf,2)*M2 + 
	     2*F10*F11*F22*kflf*pow(lipf,2)*M2 + 
	     2*F11*F20*F22*kflf*pow(lipf,2)*M2 + 
	     4*F11*F12*F20*kflf*kfpf*lipi*M2 + 
	     3*F12*F20*F21*kflf*kfpf*lipi*M2 - 
	     2*F10*F11*F22*kflf*kfpf*lipi*M2 + 
	     3*F11*F20*F22*kflf*kfpf*lipi*M2 - 
	     F10*F21*F22*kflf*kfpf*lipi*M2 + 
	     2*F20*F21*F22*kflf*kfpf*lipi*M2 - 
	     4*F10*F12*F21*kfpf*lfpf*lipi*M2 - 
	     4*F10*F11*F22*kfpf*lfpf*lipi*M2 - 
	     6*F10*F21*F22*kfpf*lfpf*lipi*M2 + 
	     8*F10*F11*F12*kfpf*lfpi*lipi*M2 + 
	     4*F10*F12*F21*kfpf*lfpi*lipi*M2 + 
	     4*F10*F11*F22*kfpf*lfpi*lipi*M2 + 
	     2*F10*F21*F22*kfpf*lfpi*lipi*M2 + 
	     2*F10*F12*F21*kfpf*lilf*lipi*M2 + 
	     2*F12*F20*F21*kfpf*lilf*lipi*M2 + 
	     6*F10*F11*F22*kfpf*lilf*lipi*M2 - 
	     2*F11*F20*F22*kfpf*lilf*lipi*M2 + 
	     4*F10*F21*F22*kfpf*lilf*lipi*M2 + 
	     8*F10*F11*F12*kflf*lipf*lipi*M2 + 
	     4*F10*F12*F21*kflf*lipf*lipi*M2 - 
	     2*F12*F20*F21*kflf*lipf*lipi*M2 + 
	     4*F10*F11*F22*kflf*lipf*lipi*M2 + 
	     2*F11*F20*F22*kflf*lipf*lipi*M2 + 
	     2*F10*F21*F22*kflf*lipf*lipi*M2 + 
	     8*F10*F12*F21*lfpf*lipf*lipi*M2 + 
	     8*F10*F11*F22*lfpf*lipf*lipi*M2 + 
	     12*F10*F21*F22*lfpf*lipf*lipi*M2 - 
	     16*F10*F11*F12*lfpi*lipf*lipi*M2 - 
	     8*F10*F12*F21*lfpi*lipf*lipi*M2 - 
	     8*F10*F11*F22*lfpi*lipf*lipi*M2 - 
	     4*F10*F21*F22*lfpi*lipf*lipi*M2 - 
	     2*F10*F12*F21*kflf*pow(lipi,2)*M2 - 
	     6*F10*F11*F22*kflf*pow(lipi,2)*M2 - 
	     6*F10*F21*F22*kflf*pow(lipi,2)*M2 - 
	     16*F10*F11*F12*lfpf*pow(lipi,2)*M2 - 
	     8*F10*F12*F21*lfpf*pow(lipi,2)*M2 - 
	     8*F10*F11*F22*lfpf*pow(lipi,2)*M2 - 
	     4*F10*F21*F22*lfpf*pow(lipi,2)*M2 + 
	     8*F10*F12*F21*lfpi*pow(lipi,2)*M2 + 
	     8*F10*F11*F22*lfpi*pow(lipi,2)*M2 + 
	     12*F10*F21*F22*lfpi*pow(lipi,2)*M2 + 
	     2*F10*F11*F22*kflf*kfpf*m_l2*M2 + 
	     2*F11*F20*F22*kflf*kfpf*m_l2*M2 + 
	     2*F10*F21*F22*kflf*kfpf*m_l2*M2 + 
	     2*F20*F21*F22*kflf*kfpf*m_l2*M2 - 
	     2*F10*F12*F21*pow(kfpf,2)*m_l2*M2 - 
	     F12*F20*F21*pow(kfpf,2)*m_l2*M2 - 
	     2*F10*F11*F22*pow(kfpf,2)*m_l2*M2 - 
	     F11*F20*F22*pow(kfpf,2)*m_l2*M2 - 
	     4*F10*F21*F22*pow(kfpf,2)*m_l2*M2 - 
	     2*F20*F21*F22*pow(kfpf,2)*m_l2*M2 + 
	     2*F10*F12*F21*kfpf*lfpf*m_l2*M2 + 
	     2*F12*F20*F21*kfpf*lfpf*m_l2*M2 - 
	     2*F10*F11*F22*kfpf*lfpf*m_l2*M2 - 
	     2*F11*F20*F22*kfpf*lfpf*m_l2*M2 + 
	     2*F10*F12*F21*kfpf*lfpi*m_l2*M2 + 
	     2*F12*F20*F21*kfpf*lfpi*m_l2*M2 - 
	     2*F10*F11*F22*kfpf*lfpi*m_l2*M2 - 
	     2*F11*F20*F22*kfpf*lfpi*m_l2*M2 - 
	     8*F10*F11*F12*kfpf*lipi*m_l2*M2 - 
	     8*F10*F12*F21*kfpf*lipi*m_l2*M2 - 
	     16*F10*F11*F22*kfpf*lipi*m_l2*M2 - 
	     12*F10*F21*F22*kfpf*lipi*m_l2*M2 - 
	     8*F10*F11*F12*kflf*lipi*M4 - 
	     4*F10*F12*F21*kflf*lipi*M4 - 
	     4*F10*F11*F22*kflf*lipi*M4 - 
	     2*F10*F21*F22*kflf*lipi*M4 + 
	     16*F10*F11*F12*lilf*lipi*M4 + 
	     8*F10*F12*F21*lilf*lipi*M4 + 
	     8*F10*F11*F22*lilf*lipi*M4 + 
	     4*F10*F21*F22*lilf*lipi*M4 - 
	     32*F10*F11*F12*lipi*m_l2*M4 - 
	     24*F10*F12*F21*lipi*m_l2*M4 - 
	     24*F10*F11*F22*lipi*m_l2*M4 - 
	     20*F10*F21*F22*lipi*m_l2*M4 - 
	     2*F12*F20*F21*kflf*kfpf*lipi*pipf + 
	     F10*F21*F22*kflf*kfpf*lipi*pipf - 
	     2*F20*F21*F22*kflf*kfpf*lipi*pipf + 
	     2*F10*F21*F22*kfpf*lfpf*lipi*pipf + 
	     2*F10*F21*F22*kfpf*lfpi*lipi*pipf - 
	     4*F10*F21*F22*kfpf*lilf*lipi*pipf + 
	     2*F12*F20*F21*kflf*lipf*lipi*pipf - 
	     2*F11*F20*F22*kflf*lipf*lipi*pipf + 
	     2*F10*F21*F22*kflf*lipf*lipi*pipf - 
	     4*F10*F21*F22*lfpf*lipf*lipi*pipf - 
	     4*F10*F21*F22*lfpi*lipf*lipi*pipf + 
	     2*F12*F20*F21*kflf*pow(lipi,2)*pipf - 
	     2*F11*F20*F22*kflf*pow(lipi,2)*pipf + 
	     2*F10*F21*F22*kflf*pow(lipi,2)*pipf - 
	     4*F10*F21*F22*lfpf*pow(lipi,2)*pipf - 
	     4*F10*F21*F22*lfpi*pow(lipi,2)*pipf + 
	     4*F10*F21*F22*kfpf*lipi*m_l2*pipf + 
	     4*F10*F12*F21*kflf*lipi*M2*pipf + 
	     4*F10*F11*F22*kflf*lipi*M2*pipf + 
	     4*F10*F21*F22*kflf*lipi*M2*pipf - 
	     8*F10*F12*F21*lilf*lipi*M2*pipf - 
	     8*F10*F11*F22*lilf*lipi*M2*pipf - 
	     8*F10*F21*F22*lilf*lipi*M2*pipf + 
	     16*F10*F11*F12*lipi*m_l2*M2*pipf + 
	     24*F10*F12*F21*lipi*m_l2*M2*pipf + 
	     24*F10*F11*F22*lipi*m_l2*M2*pipf + 
	     24*F10*F21*F22*lipi*m_l2*M2*pipf - 
	     2*F10*F21*F22*kflf*lipi*pow(pipf,2) + 
	     4*F10*F21*F22*lilf*lipi*pow(pipf,2) - 
	     4*F10*F21*F22*lipi*m_l2*pow(pipf,2) + 
	     pow(kfpi,2)*(2*F11*F20*F22*lfpf*lipf + 3*F10*F21*F22*lfpf*lipf + 
			  2*F20*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf + 
			  2*F11*F20*F22*lfpi*lipf + F10*F21*F22*lfpi*lipf + 
			  2*F10*F21*F22*lfpf*lipi - 
			  2*F20*(2*F11 + F21)*F22*kfpf*m_l2 + 
			  4*F11*F20*F22*lfpf*m_l2 + 3*F20*F21*F22*lfpf*m_l2 - 
			  F20*F21*F22*lfpi*m_l2 - 2*F11*F20*F22*lipf*m_l2 + 
			  F10*F21*F22*lipf*m_l2 - 2*F20*F21*F22*lipf*m_l2 - 
			  8*F11*F12*F20*m_l2*M2 - 
			  2*F10*F12*F21*m_l2*M2 - 
			  7*F12*F20*F21*m_l2*M2 - 
			  14*F10*F11*F22*m_l2*M2 - 
			  7*F11*F20*F22*m_l2*M2 - 
			  12*F10*F21*F22*m_l2*M2 - 
			  6*F20*F21*F22*m_l2*M2 + 
			  2*F12*F20*F21*m_l2*pipf + 2*F11*F20*F22*m_l2*pipf + 
			  4*F10*F21*F22*m_l2*pipf + 4*F20*F21*F22*m_l2*pipf + 
			  lilf*((2*F11*F20 - F10*F21 + F20*F21)*F22*lipf - 
				F20*F21*F22*lipi + 4*F11*F12*F20*M2 + 
				3*F12*F20*F21*M2 + 6*F10*F11*F22*M2 + 
				3*F11*F20*F22*M2 + 4*F10*F21*F22*M2 + 
				2*F20*F21*F22*M2 - 2*F11*F20*F22*pipf - 
				2*F10*F21*F22*pipf - 2*F20*F21*F22*pipf)) + 
	     kfpi*(F20*F21*F22*pow(kfpf,2)*m_l2 + 
		   kflf*((F10 + F20)*F21*F22*pow(lipf,2) + 
			 F20*F21*F22*pow(lipi,2) - 
			 2*(F11 + F21)*F22*m_l2*(F10*M2 + F20*pipf) + 
			 lipf*((F12*F20*F21 + (2*F10*F11 - F11*F20 + F10*F21)*F22)*
			       M2 - F10*F21*F22*pipf) + 
			 lipi*((-2*F11*F20*F22 + F10*F21*F22)*lipf - 
			       (4*F11*F12*F20 + 3*F12*F20*F21 - 2*F10*F11*F22 + 
				3*F11*F20*F22 - 2*F10*F21*F22 + 2*F20*F21*F22)*M2 + 
			       2*F20*(F11 + F21)*F22*pipf)) + 
		   kfpf*(F10*F21*F22*lfpf*lipf - 2*F12*F20*F21*lfpi*lipf - 
			 F10*F21*F22*lfpi*lipf - 2*F20*F21*F22*lfpi*lipf - 
			 F20*F21*F22*lfpf*m_l2 - F20*F21*F22*lfpi*m_l2 + 
			 F10*F21*F22*lipf*m_l2 + 
			 lipi*((2*F12*F20*F21 - 2*F11*F20*F22 - 3*F10*F21*F22)*lfpi + 
			       (2*F11*F20 + (F10 + F20)*F21)*F22*lilf + 
			       F22*(-((2*F11*F20 + F10*F21 + 2*F20*F21)*lfpf) + 
				    (2*F11*F20 - F10*F21 + 2*F20*F21)*m_l2)) + 
			 8*F10*F11*F12*m_l2*M2 + 
			 8*F11*F12*F20*m_l2*M2 + 
			 4*F10*F12*F21*m_l2*M2 + 
			 6*F12*F20*F21*m_l2*M2 + 
			 16*F10*F11*F22*m_l2*M2 + 
			 6*F11*F20*F22*m_l2*M2 + 
			 12*F10*F21*F22*m_l2*M2 + 
			 4*F20*F21*F22*m_l2*M2 - 
			 lilf*((F10 + F20)*F21*F22*lipf + 
			       2*((F21*(2*F12*F20 + 3*F10*F22 + F20*F22) + 
				   F11*(2*F12*F20 + 4*F10*F22 + F20*F22))*M2 - 
				  F21*(F12*F20 + (F10 + F20)*F22)*pipf))) + 
		   2*(-(F10*F21*F22*lfpf*pow(lipf,2)) + 
		      F12*F20*F21*lfpi*pow(lipf,2) - 
		      F11*F20*F22*lfpi*pow(lipf,2) - 
		      F10*F21*F22*lfpf*pow(lipi,2) - 
		      2*F10*F12*F21*lfpf*lipf*M2 - 
		      2*F10*F11*F22*lfpf*lipf*M2 - 
		      3*F10*F21*F22*lfpf*lipf*M2 + 
		      4*F10*F11*F12*lfpi*lipf*M2 + 
		      2*F10*F12*F21*lfpi*lipf*M2 + 
		      2*F10*F11*F22*lfpi*lipf*M2 + 
		      F10*F21*F22*lfpi*lipf*M2 - 
		      4*F10*F11*F12*lfpf*m_l2*M2 - 
		      3*F10*F12*F21*lfpf*m_l2*M2 + 
		      F10*F11*F22*lfpf*m_l2*M2 + 
		      F10*F21*F22*lfpf*m_l2*M2 + 
		      F10*F12*F21*lfpi*m_l2*M2 + 
		      F10*F11*F22*lfpi*m_l2*M2 + 
		      F10*F21*F22*lfpi*m_l2*M2 + 
		      4*F10*F11*F12*lipf*m_l2*M2 + 
		      4*F10*F12*F21*lipf*m_l2*M2 + 
		      2*F10*F11*F22*lipf*m_l2*M2 + 
		      F10*F21*F22*lipf*m_l2*M2 + 
		      8*F10*F11*F12*m_l2*M4 + 
		      6*F10*F12*F21*m_l2*M4 + 
		      6*F10*F11*F22*m_l2*M4 + 
		      5*F10*F21*F22*m_l2*M4 + 
		      F10*F21*F22*lfpf*lipf*pipf + F10*F21*F22*lfpi*lipf*pipf - 
		      F12*F20*F21*lfpf*m_l2*pipf + 
		      F11*F20*F22*lfpf*m_l2*pipf - 
		      F10*F21*F22*lfpf*m_l2*pipf - 
		      F12*F20*F21*lfpi*m_l2*pipf + 
		      F11*F20*F22*lfpi*m_l2*pipf - 
		      F10*F21*F22*lfpi*m_l2*pipf - 
		      F10*F21*F22*lipf*m_l2*pipf - 
		      4*F10*F11*F12*m_l2*M2*pipf - 
		      6*F10*F12*F21*m_l2*M2*pipf - 
		      6*F10*F11*F22*m_l2*M2*pipf - 
		      6*F10*F21*F22*m_l2*M2*pipf + 
		      F10*F21*F22*m_l2*pow(pipf,2) - 
		      lilf*(F10*((2*F11 + F21)*M2 - F21*pipf)*
			    ((2*F12 + F22)*M2 - F22*pipf) + 
			    lipf*(F10*(4*F11*F12 + 3*F12*F21 + F11*F22)*M2 + 
				  F20*(F12*F21 - F11*F22)*pipf)) + 
		      lipi*(F12*F20*F21*lfpf*lipf - F11*F20*F22*lfpf*lipf - 
			    2*F10*F21*F22*lfpf*lipf + 8*F10*F11*F12*lfpf*M2 + 
			    4*F10*F12*F21*lfpf*M2 + 
			    4*F10*F11*F22*lfpf*M2 + 
			    2*F10*F21*F22*lfpf*M2 + 
			    6*F10*F11*F22*m_l2*M2 + 
			    5*F10*F21*F22*m_l2*M2 + 
			    2*F10*F21*F22*lfpf*pipf - F10*F21*F22*m_l2*pipf + 
			    (F12*F21 - F11*F22)*lilf*(F10*M2 - F20*pipf) - 
			    2*lfpi*((-(F12*F20*F21) + F11*F20*F22)*lipf + 
				    F10*((2*F12*F21 + 2*F11*F22 + 3*F21*F22)*M2 - 
					 F21*F22*pipf)))))))/(2.*kfpi*M2);

  return (Pow6(e)/(q_12*q_22))*((res1 - res2)/kflf - (res4 - res3)/kfli);
}





//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR0 (1D, Rosenbluth events):
class TFDISTR0: public TFoamIntegrand {
public:
  TFDISTR0(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton

    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

    E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
    SetFinalFourMomenta(); // Set four-momenta for the final particles
  
    return (phi_max - phi_min)*(theta_max - theta_min)*inter_ros_sin.Eval(theta_l)*mkb;
  }
};
// The end of TFDISTR0
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1e (1D, elastic scattering, e-/mu-):
class TFDISTR1e: public TFoamIntegrand {
public:
  TFDISTR1e(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton
  
    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

    E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
    SetFinalFourMomenta(); // Set four-momenta for the final particles
  
    // Mode of calculation for bremsstrahlung:
    switch (flag_mode)
      {
      case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
	{
	  delta_sum = inter_brem_ee.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
	  break;
	}
      
      case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
	{
	  delta_sum = inter_brem_pp.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
	  break;
	}
      
      case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
	{
	  delta_sum = inter_brem_ee.Eval(theta_l) + inter_brem_pp.Eval(theta_l) + inter_brem_ep.Eval(theta_l) + inter_virt.Eval(theta_l) + inter_prime.Eval(theta_l);
	  break;
	}
      }
  
    return (phi_max - phi_min)*(theta_max - theta_min)*(1. + delta_sum)*inter_ros_sin.Eval(theta_l)*mkb;
  }
};
// The end of TFDISTR1e
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1p (1D, elastic scatterring, e+/mu+):
class TFDISTR1p: public TFoamIntegrand {
public:
  TFDISTR1p(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    theta_l = theta_min + arg[0]*(theta_max - theta_min); // Theta angle for the lepton

    E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton

    E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
    SetFinalFourMomenta(); // Set four-momenta for the final particles

    // Mode of calculation for bremsstrahlung:
    switch (flag_mode)
      {
      case 1: case 4: case 7: case 10: // Only lepton bremsstrahlung
	{
	  delta_sum = inter_brem_ee.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
	  break;
	}
      
      case 2: case 5: case 8: case 11: // Only proton bremsstrahlung
	{
	  delta_sum = inter_brem_pp.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
	  break;
	}
      
      case 3: case 6: case 9: case 12: // All the bremsstrahlung terms
	{
	  delta_sum = inter_brem_ee.Eval(theta_l) + inter_brem_pp.Eval(theta_l) - inter_brem_ep.Eval(theta_l) + inter_virt.Eval(theta_l) - inter_prime.Eval(theta_l);
	  break;
	}
      }
  
    return (phi_max - phi_min)*(theta_max - theta_min)*(1. + delta_sum)*inter_ros_sin.Eval(theta_l)*mkb;
  }
};
// The end of TFDISTR1p
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2e (4D, bremsstrahlung, e-/mu-):
class TFDISTR2e: public TFoamIntegrand {
public:
  TFDISTR2e(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    // Four arguments (basic kinematic variables):
    theta_g = arg[0]*Pi; // Theta angle for the photon
    phi_g = arg[1]*2.*Pi; // Phi angle for the photon
    E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
    theta_l = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton
  
    if (flag_mode > 3)
      {
	// Checking the photon energy:
	if (E_g > M*(E_li - m_l)/(M + E_li - Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_g))) return 0.;
      
	// Evaluation of the final lepton energy:
	if (EvalEnergy() == 0) return 0.;

	EvalKinematicParams(E_g); // Evaluation of some kinematic parameters
	SetFinalFourMomenta(); // Set four-momenta for the final particles
      }
    else // Primary soft-photon approximation
      {
	E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
      
	// Four-momenta of the final lepton and proton (elastic kinematics):
	v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m_l))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m_l))*Cos(theta_l), E_lf);
	v_pf = v_li + v_pi - v_lf;
	E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
      
	// Four-momentum of the photon:
	v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
      }

    // Mode of calculation for bremsstrahlung:
    switch (flag_mode)
      {
      case 1: // Primary soft-photon approximation, only lepton bremsstrahlung
	{
	  ret =  2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}

      case 2: // Primary soft-photon approximation, only proton bremsstrahlung
	{
	  ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}
      
      case 3: // Primary soft-photon approximation, all the terms
	{
	  ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}

      case 4: // Modified soft-photon approximation, only lepton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}

      case 5: // Modified soft-photon approximation, only proton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}

      case 6: // Modified soft-photon approximation, all the terms
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}

      case 7: // Improved soft-photon approximation, only lepton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}

      case 8: // Improved soft-photon approximation, only proton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}

      case 9: // Improved soft-photon approximation, all the terms
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() + BremSoftInterference())*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}

      case 10: // Accurate QED calculation, only lepton bremsstrahlung
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = lterm();

	  // Checking the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }

	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}

      case 11: // Accurate QED calculation, only proton bremsstrahlung
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = pterm();

	  // Checking the the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }

	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 12: // Accurate QED calculation, all the terms
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = lterm() + pterm() + lpterm();

	  // Checking the the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }
      
	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      }
  
    return 0.;
  }
};
// The end of TFDISTR2e
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2p (4D, bremsstrahlung, e+/mu+):
class TFDISTR2p: public TFoamIntegrand {
public:
  TFDISTR2p(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    // Four arguments (basic kinematic variables):
    theta_g = arg[0]*Pi; // Theta angle for the photon
    phi_g = arg[1]*2.*Pi; // Phi angle for the photon
    E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
    theta_l = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton
  
    if (flag_mode > 3)
      {
	// Checking the photon energy:
	if (E_g > M*(E_li - m_l)/(M + E_li - Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_g))) return 0.;

	// Evaluation of the final lepton energy:
	if (EvalEnergy() == 0) return 0.;

	EvalKinematicParams(E_g); // Evaluation of some kinematic parameters
	SetFinalFourMomenta(); // Set four-momenta for the final particles
      }
    else // Primary soft-photon approximation
      {
	E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
      
	// Four-momenta of the final lepton and proton (elastic kinematics):
	v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - Pow2(m_l))*Sin(theta_l), 0., Sqrt(Pow2(E_lf) - Pow2(m_l))*Cos(theta_l), E_lf);
	v_pf = v_li + v_pi - v_lf;
	E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();
      
	// Four-momentum of the photon:
	v_kf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);
      }

    // Mode of calculation for bremsstrahlung:
    switch (flag_mode)
      {
      case 1: // Primary soft-photon approximation, only lepton bremsstrahlung
	{
	  ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}

      case 2: // Primary soft-photon approximation, only proton bremsstrahlung
	{
	  ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}
      
      case 3: // Primary soft-photon approximation, all the terms
	{
	  ret = 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
      
	  // Back to the elastic kinematics:
	  E_g = 0.; theta_g = 0.; phi_g = 0.;
	  v_kf.SetPxPyPzE(0., 0., 0., 0.);
            
	  return ret;
	}

      case 4: // Modified soft-photon approximation, only lepton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 5: // Modified soft-photon approximation, only proton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 6: // Modified soft-photon approximation, all the terms
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*inter_ros_sin.Eval(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 7: // Improved soft-photon approximation, only lepton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftLeptonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 8: // Improved soft-photon approximation, only proton bremsstrahlung
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*BremSoftProtonOnly()*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 9: // Improved soft-photon approximation, all the terms
	{
	  return 2.*Pow2(Pi)*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*(BremSoftLeptonOnly() + BremSoftProtonOnly() - BremSoftInterference())*RosenbluthCS()*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 10: // Accurate QED calculation, only lepton bremsstrahlung
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = lterm();

	  // Checking the the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }

	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}

      case 11: // Accurate QED calculation, only proton bremsstrahlung
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = pterm();

	  // Checking the the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }

	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      
      case 12: // Accurate QED calculation, all the terms
	{
	  EvalAllProducts(); // Evaluation of the four-momentum products
	
	  // Square of the total amplitude:
	  M_sum = lterm() + pterm() - lpterm();

	  // Checking the the square of the amplitude:
	  if (M_sum < 0.)
	    {
	      cout << endl << "Warning: negative value! M_sum = " << M_sum << endl;
	      return 0.;
	    }

	  return (1./(512.*Pow3(Pi)*Sqrt(Pow2(E_li) - Pow2(m_l))*M))*E_g*((Pow2(E_lf) - Pow2(m_l))/Abs((Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_l) - E_g*(Cos(theta_l)*Cos(theta_g) + Sin(theta_l)*Sin(theta_g)*Cos(phi_g)))*E_lf - (E_li + M - E_g)*Sqrt(Pow2(E_lf) - Pow2(m_l))))*M_sum*(phi_max - phi_min)*(theta_max - theta_min)*(E_g_max - E_g_cut)*Sin(theta_l)*Sin(theta_g)*mkb;
	}
      }
  
    return 0.;
  }
};
// The end of TFDISTR2p
//==============================================================================================

//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1 (1D, Moller elastic scattering):
class TFDISTR1: public TFoamIntegrand {
public:
  TFDISTR1(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    // Theta angle:
    theta_1 = theta_min + (theta_max - theta_min)*arg[0];
    E_1 = m_l*(E_li + m_l + (E_li - m_l)*Pow2(Cos(theta_1)))/(E_li + m_l - (E_li - m_l)*Pow2(Cos(theta_1)));
    E_g = 0.; theta_g = 0.; phi_g = 0.;

    l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m_l))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m_l))*Cos(theta_1), E_1);
    l2f = l1i + l2i - l1f;
    return 2.*Pi*(theta_max - theta_min)*(1. + inter_virt_M.Eval(theta_1) + inter_brems_M.Eval(theta_1))*inter_Mol.Eval(theta_1)*Sin(theta_1)*mkb;
  }
};
// The end of TFDISTR1
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2 (4D, Moller bremsstrahlung in the soft-photon approximation):
class TFDISTR2: public TFoamIntegrand {
public:
  TFDISTR2(){};
  Double_t Density(G4int nDim, Double_t *arg) {

    // Four arguments (basic kinematic variables):
    theta_g = arg[0]*Pi; // Theta angle for the photon
    phi_g = arg[1]*2.*Pi; // Phi angle for the photon
    E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
    theta_1 = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton

    if (E_g > m_l*(E_li - m_l)/(m_l + E_li - Sqrt(Pow2(E_li) - Pow2(m_l))*Cos(theta_g))) return 0.;
  
    E_1 = m_l*(E_li + m_l + (E_li - m_l)*Pow2(Cos(theta_1)))/(E_li + m_l - (E_li - m_l)*Pow2(Cos(theta_1)));
    if (E_1 < m_l || E_1 > E_li) return 0.;

    l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m_l))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m_l))*Cos(theta_1), E_1);
    vkf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);

    l2f = l1i + l2i - l1f;
    if (abs(l2f*l2f - m_l*m_l) > 1.e-5) cout << endl << "Bad kinematics! #1: " << l2f*l2f - m_l*m_l << endl;

    // Soft-photon bremsstrahlung factor:
    brems_factor = (-alpha*E_g/(4.*Pi*Pi))*((((1./(vkf*l1i))*l1i + (1./(vkf*l2i))*l2i - (1./(vkf*l1f))*l1f - (1./(vkf*l2f))*l2f)*((1./(vkf*l1i))*l1i + (1./(vkf*l2i))*l2i - (1./(vkf*l1f))*l1f - (1./(vkf*l2f))*l2f)));

    cospsi = Cos(theta_1)*Cos(theta_g) + Sin(theta_1)*Sin(theta_g)*Cos(phi_g);

    // More realistic kinematics to output:
    E_1 = (m_l*(E_li - E_g) - E_li*E_g*(1. - Cos(theta_g)))/(m_l + E_li*(1. - Cos(theta_1)) - E_g*(1. - cospsi));
    if (E_1 < m_l || E_1 > E_li - E_g) return 0.;

    l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m_l))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m_l))*Cos(theta_1), E_1);

    l2f = l1i + l2i - l1f - vkf;
    if (abs(l2f*l2f - m_l*m_l) > 1.e-2) cout << endl << "Bad kinematics! #2: " << l2f*l2f - m_l*m_l << endl;

    return (2*Pi)*(2.*Pi*Pi)*(E_g_max - E_g_cut)*(theta_max - theta_min)*brems_factor*inter_Mol.Eval(theta_1)*Sin(theta_1)*Sin(theta_g)*mkb;
  }
};
// The end of TFDISTR2
//==============================================================================================








void Initialize_genevent()
{
  omega = (phi_max - phi_min)*(Cos(theta_min) - Cos(theta_max)); // Solid angle (steradian)
  // Initialization of the random number generator:
  PseRan->SetSeed(0);
  
  if (!flag_onlymoller) {
    // Set the initial lepton four-momentum (in the laboratory frame):
    v_li.SetPxPyPzE(0., 0., Sqrt(Pow2(E_li) - Pow2(m_l)), E_li);
    
    // Set the initial proton four-momentum (in the laboratory frame):
    v_pi.SetPxPyPzE(0., 0., 0., M);
    
    // Setting the functions and tolerance for numerical integration:
    i_li_lf.SetFunction(func_li_lf); // To calculate B(v_li, v_lf, E_g_cut)
    i_li_lf.SetRelTolerance(IntTol);
    
    i_li_pi.SetFunction(func_li_pi); // To calculate B(v_li, v_pi, E_g_cut)
    i_li_pi.SetRelTolerance(IntTol);
    
    i_li_pf.SetFunction(func_li_pf); // To calculate B(v_li, v_pf, E_g_cut)
    i_li_pf.SetRelTolerance(IntTol);
    
    i_lf_pi.SetFunction(func_lf_pi); // To calculate B(v_lf, v_pi, E_g_cut)
    i_lf_pi.SetRelTolerance(IntTol);
    
    i_lf_pf.SetFunction(func_lf_pf); // To calculate B(v_lf, v_pf, E_g_cut)
    i_lf_pf.SetRelTolerance(IntTol);
    
    i_pi_pf.SetFunction(func_pi_pf); // To calculate B(v_pi, v_pf, E_g_cut)
    i_pi_pf.SetRelTolerance(IntTol);
    
    // To integrate the "elastic" part of cross section (with E_g < E_g_cut):
    i_el_e.SetFunction(func_el_e); // e-/mu-
    i_el_e.SetRelTolerance(IntTol);
    
    i_el_p.SetFunction(func_el_p); // e+/mu+
    i_el_p.SetRelTolerance(IntTol);
    
    // To integrate the Rosenbluth differential cross section:
    i_ros.SetFunction(func_ros);
    i_ros.SetRelTolerance(IntTol);
  }

  if (flag_moller) {
    // Four momenta of the initial particles:
    l1i.SetPxPyPzE(0., 0., Sqrt(E_li*E_li - m_l*m_l), E_li); // Incident electron/positron
    l2i.SetPxPyPzE(0., 0. , 0., m_l); // Target electron
    
    // Mandelstam variable s:
    sss = (l1i + l2i)*(l1i + l2i);

    // Setting the functions and tolerance for numerical integration: (Moller)
    i_l1i_l1f.SetFunction(func_l1i_l1f); // To calculate B(l1i, l1f, E_g_cut)
    i_l1i_l1f.SetRelTolerance(IntTol);
    
    i_l1i_l2i.SetFunction(func_l1i_l2i); // To calculate B(l1i, l2i, E_g_cut)
    i_l1i_l2i.SetRelTolerance(IntTol);
    
    i_l1i_l2f.SetFunction(func_l1i_l2f); // To calculate B(l1i, l2f, E_g_cut)
    i_l1i_l2f.SetRelTolerance(IntTol);

    i_l1f_l2i.SetFunction(func_l1f_l2i); // To calculate B(l1f, l2i, E_g_cut)
    i_l1f_l2i.SetRelTolerance(IntTol);

    i_l1f_l2f.SetFunction(func_l1f_l2f); // To calculate B(l1f, l2f, E_g_cut)
    i_l1f_l2f.SetRelTolerance(IntTol);

    i_l2i_l2f.SetFunction(func_l2i_l2f); // To calculate B(l2i, l2f, E_g_cut)
    i_l2i_l2f.SetRelTolerance(IntTol);

    i_mol.SetFunction(func_mol);
    i_mol.SetRelTolerance(IntTol);
    
  }

  // Interpolation of bremsstrahlung and virtual-photon corrections:
  if (flag_vpol == 3) // In the case of accurate calculation of the vacuum polarization
    {
      fvpol = fopen("vpol.dat", "r");
      G4int npoints = 0;
      if (fvpol !=NULL) { 
	char str[128];
	  
	while (!feof(fvpol)) // Reading from the file "vpol.dat"
	  {
	    str[0] = 0;
	    fgets(str, 128, fvpol);
	    if (feof(fvpol) || strlen(str) == 0) break; // The end or empty string
	      
	    if (str[0] != '/')
	      {
		sscanf(str, "%lf %lf", &ss[npoints], &rep[npoints]);
		rep[npoints] = 2.*rep[npoints];
		npoints++;
	      }
	  }
	fclose(fvpol); // Closing the file "vpol.dat"
      }
      else {G4cout << "no polarization file" << G4endl;}
      inter_vpol.SetData(npoints, ss, rep); // Interpolating
    }
      
  E_g = 0.; theta_g = 0.; phi_g = 0.; // Elastic scattering
  
  for (G4int i = 0; i < InterpolPoints; i++)
    {
      if (!flag_onlymoller) {
	theta_l = theta_min + i*(theta_max - theta_min)/(InterpolPoints - 1); // Theta angle for the lepton
	E_lf = ElasticEnergy(theta_l); // Full energy of the scattered lepton
	EvalKinematicParams(0.); // Evaluation of some kinematic parameters
	SetFinalFourMomenta();   // Set four-momenta for the final particles
	xx[i] = theta_l;               // Theta angle for the lepton
	y_ros_sin[i] = f_ros(theta_l); // RosenbluthCS()*Sin(theta_l)
	y_brem_ee[i] = d_brem_ee();    // Bremsstrahlung correction from the lepton term
	y_brem_pp[i] = d_brem_pp();    // Bremsstrahlung correction from the proton term
	y_brem_ep[i] = d_brem_ep();    // Bremsstrahlung correction from the interference term
	y_prime[i] = d_prime();        // TPE contribution by Maximon & Tjon
	if (flag_vpol == 1) y_virt[i] = d_vac_e() + d_vertex();
	if (flag_vpol == 2) y_virt[i] = d_vac_e() + d_vac_mu() + d_vac_tau() + d_vertex();
	if (flag_vpol == 3) y_virt[i] = inter_vpol.Eval(-qq) + d_vertex();
      }
      
      if (flag_moller) {
	// moller kinematics
	theta_1 = theta_min + i*(theta_max - theta_min)/(InterpolPoints - 1);
	E_1 = m_l*(E_li + m_l + (E_li - m_l)*Pow2(Cos(theta_1)))/(E_li + m_l - (E_li - m_l)*Pow2(Cos(theta_1)));
	
	l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m_l))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m_l))*Cos(theta_1), E_1);
	l2f = l1i + l2i - l1f;
	  
	E_2 = l2f.E();
	qq_M = (l1i - l1f)*(l1i - l1f);
	xx[i] = theta_1;
	y_virt_M[i] = inter_vpol.Eval(-qq_M) + d_vertex_M();
	y_brems_M[i] = d_brems_M(); // Bremsstrahlung correction from the lepton term
	y_Mol[i] = MCS();
	  
      }
    }
  
  // Interpolating functions:
  if (!flag_onlymoller) {
    inter_ros_sin.SetData(InterpolPoints, xx, y_ros_sin); // RosenbluthCS()*Sin(theta_l)
    inter_brem_ee.SetData(InterpolPoints, xx, y_brem_ee); // Lepton bremsstrahlung term
    inter_brem_pp.SetData(InterpolPoints, xx, y_brem_pp); // Proton bremsstrahlung term
    inter_brem_ep.SetData(InterpolPoints, xx, y_brem_ep); // Interference bremsstrahlung term
    inter_virt.SetData(InterpolPoints, xx, y_virt);       // Virtual-photon correction
    inter_prime.SetData(InterpolPoints, xx, y_prime);     // TPE contribution by Maximon & Tjon
  }

  if (flag_moller) {
    inter_brems_M.SetData(InterpolPoints, xx, y_brems_M); // Moller bremsstrahlung
    inter_virt_M.SetData(InterpolPoints, xx, y_virt_M);       // Virtual-photon correction
    inter_Mol.SetData(InterpolPoints, xx, y_Mol);         // Moller elastic cross section
  }
  //==============================================================================================

  //----------------------------------------------------------------------------------------------
  // mFOAM simulators (S. Jadach, P. Sawicki, arXiv:physics/0506084):

  // Some options for the simulators:
  if (flag_quick == false) // Default options
    {
      nCells_1D = 1000;  // Number of cells for 1D case
      nSampl_1D = 200;   // Number of samples for 1Dcase
      nCells_4D = 30000; // Number of cells for 4D case
      nSampl_4D = 1500;  // Number of samples for 4D case
    }
  else // Options for the case of quick calculation
    {
      nCells_1D = 500;   // Number of cells for 1D case
      nSampl_1D = 100;   // Number of samples for 1Dcase
      nCells_4D = 10000; // Number of cells for 4D case
      nSampl_4D = 600;   // Number of samples for 4D case
    }
  
  if (flag_rosen == true && !flag_onlymoller) // Rosenbluth events
    {
      // Distribution function:
      TFDISTR0 *Rho_Ros = new TFDISTR0();

      // Initialization of the Foam_Ros simulator (Rosenbluth events):
      cout << endl << "Initialization (Rosenbluth scattering):" << endl;
      Foam_Ros->SetkDim(1);           // Set number of dimensions
      Foam_Ros->SetnCells(nCells_1D); // Set number of cells
      Foam_Ros->SetnSampl(nSampl_1D); // Set number os samples
      Foam_Ros->SetOptRej(1);         // Unweighted events in MC generation
      Foam_Ros->SetRho(Rho_Ros);      // Set distribution function
      Foam_Ros->SetPseRan(PseRan);    // Set random number generator
      Foam_Ros->SetChat(1);           // Set "chat level" in the standard output
      Foam_Ros->Initialize();         // Initialization
    }
  
  if (flag_lepton != 2 && flag_lepton != 5)  // For negatively charged leptons
    {
      if (!flag_onlymoller) {
	// Distribution functions:
	TFoamIntegrand *Rho_ElE = new TFDISTR1e(); // Elastic scattering, e-/mu-
	TFoamIntegrand *Rho_BrE = new TFDISTR2e(); // First-order bremsstrahlung, e-/mu-
      
	cout << endl;
      
	// Initialization of the Foam_ElE1 simulator (elastic scattering, e-/mu-):
	if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (elastic scattering, e-):" << endl;
	if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (elastic scattering, mu-):" << endl;
	Foam_ElE->SetkDim(1);           // Set number of dimensions
	Foam_ElE->SetnCells(nCells_1D); // Set number of cells
	Foam_ElE->SetnSampl(nSampl_1D); // Set number os samples
	Foam_ElE->SetOptRej(1);         // Unweighted events in MC generation
	Foam_ElE->SetRho(Rho_ElE);      // Set distribution function
	Foam_ElE->SetPseRan(PseRan);    // Set random number generator
	Foam_ElE->SetChat(1);           // Set "chat level" in the standard output
	Foam_ElE->Initialize();         // Initialization
      
	cout << endl;
  
	en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
  
	// Initialization of the Foam_BrE1 simulator (bremsstrahlung, e-/mu-, root "-"):
	if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e-, root \"-\"):" << endl;
	if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu-, root \"-\"):" << endl;
	Foam_BrE1->SetkDim(4);           // Set number of dimensions
	Foam_BrE1->SetnCells(nCells_4D); // Set number of cells
	Foam_BrE1->SetnSampl(nSampl_4D); // Set number os samples
	Foam_BrE1->SetOptRej(1);         // Unweighted events in MC generation
	Foam_BrE1->SetRho(Rho_BrE);      // Set distribution function
	Foam_BrE1->SetPseRan(PseRan);    // Set random number generator
	Foam_BrE1->SetChat(1);           // Set "chat level" in the standard output
	Foam_BrE1->Initialize();         // Initialization
	
	if (flag_mode > 3)
	  {
	    cout << endl;
	    
	    en_sign = +1.; // The "+" root of the equation for the scattered lepton energy  

	    // Initialization of the Foam_BrE2 simulator (bremsstrahlung, e-/mu-, root "+"):
	    if (flag_lepton == 1 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e-, root \"+\"):" << endl;
	    if (flag_lepton == 4 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu-, root \"+\"):" << endl;
	    Foam_BrE2->SetkDim(4);           // Set number of dimensions
	    Foam_BrE2->SetnCells(nCells_4D); // Set number of cells
	    Foam_BrE2->SetnSampl(nSampl_4D); // Set number os samples
	    Foam_BrE2->SetOptRej(1);         // Unweighted events in MC generation
	    Foam_BrE2->SetRho(Rho_BrE);      // Set distribution function
	    Foam_BrE2->SetPseRan(PseRan);    // Set random number generator
	    Foam_BrE2->SetChat(1);           // Set "chat level" in the standard output
	    Foam_BrE2->Initialize();         // Initialization
	    
	    bre2 = Foam_BrE2->GetPrimary();
	    if (bre2 == 0) cout << "Foam_BrE2 (bremsstrahlung, e-/mu-, root \"+\"): no need to use." << endl;
	  }
	else bre2 = 0;
      }
      
      if (flag_moller)
	{
	  
	  TFoamIntegrand *Rho1 = new TFDISTR1(); // Moller elastic scattering
	  TFoamIntegrand *Rho2 = new TFDISTR2(); // Moller bremsstrahlung

	  // Initialization of the Foam1 simulator (Moller elastic scattering):
	  cout << endl << "Initialization of Foam1 (Moller elastic scattering):" << endl;
	  Foam1->SetkDim(1);        // Set number of dimensions
	  Foam1->SetnCells(nCells_1D);   // Set number of cells
	  Foam1->SetnSampl(nSampl_1D);    // Set number os samples
	  Foam1->SetOptRej(1);      // Unweighted events in MC generation
	  Foam1->SetRho(Rho1);      // Set distribution function
	  Foam1->SetPseRan(PseRan); // Set random number generator
	  Foam1->SetChat(1);        // Set "chat level" in the standard output
	  Foam1->Initialize();      // Initialization
	    
	  // Initialization of the Foam2 simulator (Moller bremsstrahlung):
	  cout << endl << "Initialization of Foam2 (Moller bremsstrahlung):" << endl;
	  Foam2->SetkDim(4);        // Set number of dimensions
	  Foam2->SetnCells(nCells_4D);   // Set number of cells
	  Foam2->SetnSampl(nSampl_4D);    // Set number os samples
	  Foam2->SetOptRej(1);      // Unweighted events in MC generation
	  Foam2->SetRho(Rho2);      // Set distribution function
	  Foam2->SetPseRan(PseRan); // Set random number generator
	  Foam2->SetChat(1);        // Set "chat level" in the standard output
	  Foam2->Initialize();      // Initialization
	    
	  mbr = Foam2->GetPrimary();
	  if (mel == 0) cout << "Foam1 (Moller bremsstrahlung): no need to use." << endl;
	    
	}
      else
	{mel = 0; mbr = 0; mbr_error = 0;}
	
    }
    
  if (flag_lepton != 1 && flag_lepton != 4 && !flag_onlymoller)  // For positively charged leptons
    {
      // Distribution functions:
      TFoamIntegrand *Rho_ElP = new TFDISTR1p(); // Elastic scattering, e+/mu+
      TFoamIntegrand *Rho_BrP = new TFDISTR2p(); // First-order bremsstrahlung, e+/mu+

      cout << endl;
  
      // Initialization of the Foam_ElP simulator (elastic scattering, e+/mu+):
      if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (elastic scattering, e+):" << endl;
      if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (elastic scattering, mu+):" << endl;
      Foam_ElP->SetkDim(1);           // Set number of dimensions
      Foam_ElP->SetnCells(nCells_1D); // Set number of cells
      Foam_ElP->SetnSampl(nSampl_1D); // Set number os samples
      Foam_ElP->SetOptRej(1);         // Unweighted events in MC generation
      Foam_ElP->SetRho(Rho_ElP);      // Set distribution function
      Foam_ElP->SetPseRan(PseRan);    // Set random number generator
      Foam_ElP->SetChat(1);           // Set "chat level" in the standard output
      Foam_ElP->Initialize();         // Initialization

      cout << endl;
  
      en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
  
      // Initialization of the Foam_BrP simulator (bremsstrahlung, e+/mu+, root "-"):
      if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e+, root \"-\"):" << endl;
      if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu+, root \"-\"):" << endl;
      Foam_BrP1->SetkDim(4);           // Set number of dimensions
      Foam_BrP1->SetnCells(nCells_4D); // Set number of cells
      Foam_BrP1->SetnSampl(nSampl_4D); // Set number os samples
      Foam_BrP1->SetOptRej(1);         // Unweighted events in MC generation
      Foam_BrP1->SetRho(Rho_BrP);      // Set distribution function
      Foam_BrP1->SetPseRan(PseRan);    // Set random number generator
      Foam_BrP1->SetChat(1);           // Set "chat level" in the standard output
      Foam_BrP1->Initialize();         // Initialization
  
      if (flag_mode > 3)
	{
	  cout << endl;

	  en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
	  
	  // Initialization of the Foam_BrP simulator (bremsstrahlung, e+/mu+, root "+"):
	  if (flag_lepton == 2 || flag_lepton == 3) cout << "Initialization (bremsstrahlung, e+, root \"+\"):" << endl;
	  if (flag_lepton == 5 || flag_lepton == 6) cout << "Initialization (bremsstrahlung, mu+, root \"+\"):" << endl;
	  Foam_BrP2->SetkDim(4);           // Set number of dimensions
	  Foam_BrP2->SetnCells(nCells_4D); // Set number of cells
	  Foam_BrP2->SetnSampl(nSampl_4D); // Set number os samples
	  Foam_BrP2->SetOptRej(1);         // Unweighted events in MC generation
	  Foam_BrP2->SetRho(Rho_BrP);      // Set distribution function
	  Foam_BrP2->SetPseRan(PseRan);    // Set random number generator
	  Foam_BrP2->SetChat(1);           // Set "chat level" in the standard output
	  Foam_BrP2->Initialize();         // Initialization
	  
	  brp2 = Foam_BrP2->GetPrimary();
	  if (brp2 == 0) cout << "Foam_BrP2 (bremsstrahlung, e+/mu+, root \"+\"): no need to use." << endl;
	}
      else brp2 = 0;
    }
  //==============================================================================================

    
  //----------------------------------------------------------------------------------------------
  // Initialization of the generator (calculation of the bremsstrahlung cross sections):
  
  // The number of events for initialization (default: 10 million):
  if (!flag_quick) ninit = G4long (1e7);
  else ninit = G4long (1e6);
  
  cout << setprecision(3) << endl;
  en_sign = -1.; // The "-" root of the equation for the scattered lepton energy

  // Negatively charged leptons:------------------------------------------------------------------
  if ((flag_lepton == 1 || flag_lepton == 4) && !flag_onlymoller)
    {
      for (loop = 0; loop < ninit; loop++)
	{
	  Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event for e-/mu-, root "-"
    
	  // The progress bar:
	  if (loop % G4long (0.01*(ninit - 1)))
	    {
	      cout << "Initialization 1 (" << ninit << " events): ";
	      cout <<  G4long (loop/(0.01*(ninit - 1)))<< "%" << '\r' << flush;
	    }
	}
      
      cout << endl;  
      Foam_BrE1->GetIntegMC(bre1, bre_error); // Cross section for e-/mu-, root "-"
    }

  // Positively charged leptons:------------------------------------------------------------------
  if ((flag_lepton == 2 || flag_lepton == 5) && !flag_onlymoller)
    {
      for (loop = 0; loop < ninit; loop++)
	{
	  Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event for e+/mu+, root "-"
	  
	  // The progress bar:
	  if (loop % G4long (0.01*(ninit - 1)))
	    {
	      cout << "Initialization 1 (" <<ninit << " events): ";
	      cout <<  G4long (loop/(0.01*(ninit - 1)))<< "%" << '\r' << flush;
	    }
	}
      
      cout << endl;
      Foam_BrP1->GetIntegMC(brp1, brp_error); // Cross section for e+/mu+, root "-"
    }

  // Both negatively and positively charged leptons:----------------------------------------------
  if ((flag_lepton == 3 || flag_lepton == 6) && !flag_onlymoller)
    {
      for (loop = 0; loop < ninit; loop++)
	{
	  Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event for e-/mu-, root "-"
	  Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event for e+/mu+, root "-"
    
	  // The progress status:
	  if (loop % G4long (0.01*(ninit - 1)))
	    {
	      cout << "Initialization 1 (" << ninit << " events): ";
	      cout << G4long (loop/(0.01*(ninit - 1)))<< "%" << '\r' << flush;
	    }
	}

      cout << endl;
      Foam_BrE1->GetIntegMC(bre1, bre_error); // Cross section for e-/mu-, root "-"
      Foam_BrP1->GetIntegMC(brp1, brp_error); // Cross section for e+/mu+, root "-"
    }

  // Initialization 2 (for "+" roots):
  if (bre2 == 0 && brp2 == 0 || flag_onlymoller) cout << "Initialization 2: no need to perform";
  else
    {
      en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
      
      for (loop = 0; loop < G4long (0.05*ninit); loop++)
	{
	  // Generating a bremsstrahlung event for e-/mu-, root "+":
	  if (bre2 > 0) Foam_BrE2->MakeEvent();
	  
	  // Generating a bremsstrahlung event for e+/mu+, root "+":
	  if (brp2 > 0) Foam_BrP2->MakeEvent();
	  
	  // The progress status:
	  if (loop % G4long (0.01*(G4long (0.05*ninit) - 1)))
	    {
	      cout << "Initialization 2 (" << G4long (0.05*ninit) << " events): ";
	      cout << G4long (loop/(0.01*(G4long (0.05*ninit) - 1))) << "%" << '\r' << flush;
	    }
	}
    }

  cout << endl;
  if (bre2 > 0) Foam_BrE2->GetIntegMC(bre2, bre_error); // Cross section for e-/mu-, root "+"
  if (brp2 > 0) Foam_BrP2->GetIntegMC(brp2, brp_error); // Cross section for e+/mu+, root "+"

  // Initialization 3:
  if (mbr == 0 || !flag_moller) cout << "Initialization 3: no need to perform" <<endl;
  else
    {
      for (loop = 0; loop < ninit; loop++)
	{
	  // Moller bremsstrahlung:
	  Foam2->MakeEvent();
	
	  // Progress bar:
	  if (loop % int (0.01*(ninit - 1)))
	    {
	      cout << "Initialization: " << int (loop/((0.01*(ninit - 1)))) << "%" << '\r' << flush;
	    }
	}
      Foam2->GetIntegMC(mbr, mbr_error);
    }
    
  //==============================================================================================
  // Elastic scattering cross sections:
  E_g = 0.; theta_g = 0.; phi_g = 0.;
  if (flag_onlymoller) {
    elast_e = 0.;
    elast_p = 0.;
    brp1 = 0.;
    brp2 = 0.;
    bre1 = 0.;
    bre2 = 0.;
  }
  else if (flag_lepton==1 || flag_lepton==4) {
    elast_e = (phi_max - phi_min)*mkb*i_el_e.Integral(theta_min, theta_max); // e-/mu-
    elast_p = 0.;
    brp1 = 0.;
    brp2 = 0.;
  }
  else if (flag_lepton==2 || flag_lepton==5) {
    elast_e = 0.;
    elast_p = (phi_max - phi_min)*mkb*i_el_p.Integral(theta_min, theta_max); // e+/mu+
    bre1 = 0.;
    bre2 = 0.;
  }
  else if (flag_lepton==3 || flag_lepton==6) {
    elast_e = (phi_max - phi_min)*mkb*i_el_e.Integral(theta_min, theta_max); // e-/mu-
    elast_p = (phi_max - phi_min)*mkb*i_el_p.Integral(theta_min, theta_max); // e+/mu+
  }

  // Rosenbluth cross section:
  if (flag_rosen) rosen = (phi_max - phi_min)*mkb*i_ros.Integral(theta_min, theta_max);
  else rosen = 0.;
  
  // Cross section for Moller elastic part:
  if (flag_moller) mel = (phi_max - phi_min)*mkb*i_mol.Integral(theta_min, theta_max);
  else mel = 0.;

  sum_all = rosen + elast_e + bre1 + bre2 + mel + mbr + elast_p + brp1 + brp2;
  //----------------------------------------------------------------------------------------------
  // The numbers of events to generate:
  nev_ElE = 0; // Elastic scattering, e-/mu-
  nev_ElP = 0; // Elastic scattering, e+/mu+
  
  nev_BrE1 = 0; // Bremsstrahlung, e-/mu-, root "-"
  nev_BrE2 = 0; // Bremsstrahlung, e-/mu-, root "+"

  nev_BrP1 = 0; // Bremsstrahlung, e+/mu+, root "-"
  nev_BrP2 = 0; // Bremsstrahlung, e+/mu+, root "+"

  nev_ElM = 0; // Moller Elastic scattering
  nev_BrM = 0; // Moller Bremsstrahlung
  
  nev_0 = 0; // Number of Rosenbluth events
  nev_e = 0; // Total number of events for e-/mu-
  nev_p = 0; // Total number of events for e+/mu+
  nev_M = 0; // Total number of Moller events
     
  nev_0 = G4long (nevents*rosen/sum_all); // Rosenbluth scattering
      
  nev_ElE = G4long (nevents*elast_e/sum_all); // Elastic scattering e-/mu-
  nev_ElP = G4long (nevents*elast_p/sum_all); // Elastic scattering e+/mu+
  
  nev_BrE2 = G4long (nevents*bre2/sum_all); // Bremsstrahlung, root "+"
  nev_BrE1 = G4long (nevents*bre1/sum_all); // Bremsstrahlung, root "-"

  nev_BrP1 = G4long (nevents*brp1/sum_all);
  nev_BrP2 = G4long (nevents*brp2/sum_all);

  nev_ElM = G4long (nevents*mel/sum_all);
  nev_BrM = G4long (nevents*mbr/sum_all);

  if (!flag_onlymoller) nev_ElE = nevents - (nev_ElP + nev_BrE1 + nev_BrE2 + nev_BrP1 + nev_BrP2 + nev_ElM + nev_BrM);
  else nev_ElM = nevents - nev_BrM;
 
  nev_e = nev_ElE + nev_BrE1 + nev_BrE2; // Total number of events for e-/mu-
  nev_p = nev_ElP + nev_BrP1 + nev_BrP2; // Total number of events for e+/mu+
  nev_M = nev_ElM + nev_BrM;

  G4cout << G4endl;
  G4cout << "List of events to be generated:" << G4endl;
  G4cout << "-------------------------------" << G4endl;
  G4cout << "total number of events: " << nevents << G4endl;
  G4cout << "number of rosenbluth events: " << nev_0 << G4endl;
  G4cout << "total number of l- events: " << nev_e+nev_M << G4endl;
  G4cout << "  number of elastic l- events: " << nev_ElE << G4endl;
  G4cout << "  number of bremsstrahlung l- events: " << nev_BrE1 + nev_BrE2 << G4endl;
  G4cout << "  number of elastic moller events: " << nev_ElM << G4endl;
  G4cout << "  number of bremsstrahlung moller events: " << nev_BrM << G4endl;
  G4cout << "total number of l+ events: " << nev_p << G4endl;
  G4cout << "  number of elastic l+ events: " << nev_ElP << G4endl;
  G4cout << "  number of bremsstrahlung l+ events: " << nev_BrP1 + nev_BrP2 << G4endl;
  G4cout << G4endl;
  
  loop=0;
  
}

void Loop_genevent()
{

  //----------------------------------------------------------------------------------------------  
  // Generating events for the Rosenbluth scattering:
  if (loop<nev_0) {
    event_type = "ros";
    Foam_Ros->MakeEvent(); // Generating a Rosenbluth event
      
      
    // The progress status for Rosenbluth events:
    if ((loop % G4long (0.01*(nev_0 - 1)))==0) {
      cout << "Generation of " << nev_0 << " Rosenbluth events: " <<  G4long (loop/(0.01*(nev_0 - 1)))<< "%" << '\r' << flush;
    }
  }
  
  if (nev_0!=0 && loop==nev_0) cout << endl;
  
  // The end of generating events for the Rosenbluth scattering
  //==============================================================================================

  //----------------------------------------------------------------------------------------------
  // Generating events for the case of e-/mu-:
  
  if (nev_0<=loop && loop<nev_0+nev_e) {
    // "Elastic" scattering events: --------------------------------------------------------------
    if ((nev_ElE - count_ElE) > 0 && PseRan->Rndm() < 1.*(nev_ElE - count_ElE)/(nev_e - count_ElE - count_BrE1 - count_BrE2)) {
      event_type = "ele";
      Foam_ElE->MakeEvent(); // Generating an elastic event, e-/mu-
      count_ElE++; // Counter for elastic scattering events, e-/mu-
    }
    // Bremsstrahlung events: --------------------------------------------------------------------
    else {
      
      if (nev_BrE2 - count_BrE2 > 0 && count_BrE2/nev_BrE2 < count_BrE1/nev_BrE1) {
	event_type = "bre2";
	en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
	Foam_BrE2->MakeEvent(); // Generating a bremsstrahlung event, e-/mu-, root "+"
	count_BrE2++; // Counter for bremsstrahlung events, e-/mu-, root "+"
	     
      }
      else {
	event_type = "bre1";
	en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
	Foam_BrE1->MakeEvent(); // Generating a bremsstrahlung event, e-/mu-, root "-"
	count_BrE1++; // Counter for bremsstrahlung events, e-/mu-, root "-"
      }
    }

    // The progress status for e-/mu- events:
    if (((loop-nev_0) % G4long (0.01*(nev_e - 1)))==0) {
      if (flag_lepton == 1 || flag_lepton == 3) cout << "Generation of " << nev_e << " events for e-: ";
      if (flag_lepton == 4 || flag_lepton == 6) cout << "Generation of " << nev_e << " events for mu-: ";
      cout <<  G4long ((loop-nev_0)/(0.01*(nev_e - 1)))<< "%" << '\r' << flush;
    }
  }
  
  if (nev_e!=0 && loop==nev_0+nev_e) cout << endl;

  // The end of generating events for the case of e-/mu-
  //==============================================================================================
  
  //----------------------------------------------------------------------------------------------
  // Generating events for the case of e+/mu+:
  if (nev_0+nev_e<=loop && loop<nev_0+nev_e+nev_p) {
    // "Elastic" scattering events: --------------------------------------------------------------
    if ((nev_ElP - count_ElP) > 0 && PseRan->Rndm() < 1.*(nev_ElP - count_ElP)/(nev_p - count_ElP - count_BrP1 - count_BrP2)) {
      event_type = "elp";
      Foam_ElP->MakeEvent(); // Generating an elastic event, e+/mu+
      count_ElP++; // Counter for elastic scattering events, e+/mu+
    }
    // Bremsstrahlung events: --------------------------------------------------------------------
    else {
      if (nev_BrP2 - count_BrP2 > 0 && count_BrP2/nev_BrP2 < count_BrP1/nev_BrP1) {
	event_type = "brp2";
	en_sign = +1.; // The "+" root of the equation for the scattered lepton energy
	Foam_BrP2->MakeEvent(); // Generating a bremsstrahlung event, e+/mu+, root "+"
	count_BrP2++; // Counter for bremsstrahlung events, e+/mu+, root "+"
      }
      else {
	event_type = "brp1";
	en_sign = -1.; // The "-" root of the equation for the scattered lepton energy
	Foam_BrP1->MakeEvent(); // Generating a bremsstrahlung event, e+/mu+, root "-"
	count_BrP1++; // Counter for bremsstrahlung events, e+/mu+, root "-"
      }
    }

    // The progress status for e+/mu+ events:
    if (((loop-nev_0-nev_e) % G4long (0.01*(nev_p - 1)))==0)
      {
	if (flag_lepton == 2 || flag_lepton == 3) cout << "Generation of " << nev_p << " events for e+: ";
	if (flag_lepton == 5 || flag_lepton == 6) cout << "Generation of " << nev_p << " events for mu+: ";
	cout << G4long ((loop-nev_0-nev_e)/(0.01*(nev_p - 1))) << "%" << '\r' << flush;
      }
  }
  
  if (nev_p!=0 && loop==nev_0+nev_e+nev_p) cout << endl;
  
  // The end of generating events for the case of e+/mu+
  //==============================================================================================

  //----------------------------------------------------------------------------------------------
  // Elastic Moller scattering:
  if (nev_0+nev_e+nev_p<=loop && loop<nevents) {
    if ((count_ElM < nev_ElM) && PseRan->Rndm() < 1.*(nev_ElM - count_ElM)/(nev_M - count_ElM - count_BrM)) {
      event_type = "elm";
      Foam1->MakeEvent();
      count_ElM++;
    }
    else { // Moller bremsstrahlung:
      event_type = "brm";
      Foam2->MakeEvent();
      count_BrM++;
    }

    // Progress bar:
    if (((loop-nev_0-nev_e-nev_p) % G4long (0.01*(nev_M - 1)))==0) {
      cout << "Generation of " << nev_M << " moller events: " << G4long ((loop-nev_0-nev_e-nev_p)/(0.01*(nev_M - 1))) << "%" << '\r' << flush;
    }

  }
  
  if (loop==nevents) cout << endl;
  
  if (loop<nev_0+nev_e+nev_p) {
    // Generating azimuthal angle (phi) for the lepton:
    phi_l = phi_min + (phi_max - phi_min)*(PseRan->Rndm());
    
    // Rotation the final particle four-momenta around the z-axis:
    v_lf.RotateZ(phi_l); // Four-momentum of the final lepton
    v_pf.RotateZ(phi_l); // Four-momentum of the final proton
    v_kf.RotateZ(phi_l); // Four-momentum of the photon

    // Kinematic parameters of the final particles after rotation:
    E_lf = v_lf.E(); theta_l = v_lf.Theta(); phi_l = v_lf.Phi(); // Final lepton
    E_p = v_pf.E(); theta_p = v_pf.Theta(); phi_p = v_pf.Phi();  // Final proton
    E_g = v_kf.E(); theta_g = v_kf.Theta(); phi_g = v_kf.Phi();  // Photon

  }
  else {
    // Generating azimuthal angle (phi) for the lepton:
    phi_1 = phi_min + (phi_max - phi_min)*(PseRan->Rndm());
    
    // Rotation the final particle four-momenta around the z-axis:
    l1f.RotateZ(phi_1); // Four-momentum of the final lepton 1
    l2f.RotateZ(phi_1); // Four-momentum of the final lepton 2
    vkf.RotateZ(phi_1); // Four-momentum of the bremsstrahlung photon

    // Kinematic parameters to write to the output file:
    E_1 = l1f.E(); theta_1 = l1f.Theta(); phi_1 = l1f.Phi();
    E_2 = l2f.E(); theta_2 = l2f.Theta(); phi_2 = l2f.Phi();
    E_g = vkf.E(); theta_g = vkf.Theta(); phi_g = vkf.Phi();

  }
  
  loop++;;
    
  //==============================================================================================
}
