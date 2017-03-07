//----------------------------------------------------------------------------------------------
// Program MOLLER (radiative Moller scattering in the soft-photon approximation)
//
// (c) Alexander Gramolin, April 2014
//----------------------------------------------------------------------------------------------

// To support large output files:
#define _FILE_OFFSET_BITS 64

#include <Riostream.h>
#include <sstream>

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#define ninit 10000000 // Number of events for initialization

// Parameters for numerical integration:
#define IntOpt ROOT::Math::IntegrationOneDim::kADAPTIVE
#define IntTol 0.00001

// Parameters for interpolation:
#define InterpolPoints 1000
#define InterpolType ROOT::Math::Interpolation::kCSPLINE

using namespace std;

const double Pi = TMath::Pi();
const double degrad = Pi/180.;
const double m = 0.51099906e-3;      // Electron/positron mass (in GeV)
const double alpha = 1./137.0359895; // Fine-structure constant
const double mkb = 389.379404;       // GeV^{-2} to mkbarn conversion


// Some mathematical functions:

double Log (double arg) // Logarithm
  {
  return TMath::Log(arg);
  }

double Sqrt (double arg) // Square root
  {
  return TMath::Sqrt(arg);
  }

double Pow2 (double arg) // Pow2
  {
  return arg*arg;
  }

double Pow4 (double arg) // Pow4
  {
  return arg*arg*arg*arg;
  }

double Sin (double arg) // Sine
  {
  return TMath::Sin(arg);
  }

double Cos (double arg) // Cosine
  {
  return TMath::Cos(arg);
  }


// Incident lepton four-momentum:
TLorentzVector l1i;

// Target lepton four-momentum:
TLorentzVector l2i;

// Four-momenta of the final particles:
TLorentzVector l1f; // Scattered lepton #1
TLorentzVector l2f; // Scattered lepton #2
TLorentzVector vkf; // Bremsstrahlung photon

TLorentzVector p_x;

double s, t, u; // Mandelstam variables

double Adir, Aex, Aint;

double cospsi, omega, qq;

double A, B, C;

double E_li; // Initial lepton energy (in GeV)

double theta_min, theta_max;

double E_g_cut, E_g_max;

int i;
long n;
long nevents;
long nev_el, count_nev_el = 0;

// Kinematic variables:
double E_1, theta_1, phi_1; // For the lepton #1
double E_2, theta_2, phi_2; // For the lepton #2
double E_g, theta_g, phi_g; // For the bremsstrahlung photon

double brems_factor; // To get the soft-photon bremsstrahlung cross section

// The random number generator:
TRandom3 *PseRan = new TRandom3();

time_t starttime, stoptime; // For the timer

char mychar[64];


// Some functions: -----------------------------------------------------------------------------
double f_l1i_l1f(double x);
double f_l1i_l2i(double x);
double f_l1i_l2f(double x);
double f_l1f_l2i(double x);
double f_l1f_l2f(double x);
double f_l2i_l2f(double x);

ROOT::Math::Functor1D func_l1i_l1f(&f_l1i_l1f);
ROOT::Math::Functor1D func_l1i_l2i(&f_l1i_l2i);
ROOT::Math::Functor1D func_l1i_l2f(&f_l1i_l2f);
ROOT::Math::Functor1D func_l1f_l2i(&f_l1f_l2i);
ROOT::Math::Functor1D func_l1f_l2f(&f_l1f_l2f);
ROOT::Math::Functor1D func_l2i_l2f(&f_l2i_l2f);

// Integrators: --------------------------------------------------------------------------------
ROOT::Math::GSLIntegrator i_l1i_l1f(IntOpt);
ROOT::Math::GSLIntegrator i_l1i_l2i(IntOpt);
ROOT::Math::GSLIntegrator i_l1i_l2f(IntOpt);
ROOT::Math::GSLIntegrator i_l1f_l2i(IntOpt);
ROOT::Math::GSLIntegrator i_l1f_l2f(IntOpt);
ROOT::Math::GSLIntegrator i_l2i_l2f(IntOpt);

// Interpolators: ------------------------------------------------------------------------------
ROOT::Math::Interpolator inter_vpol(10000, InterpolType);
ROOT::Math::Interpolator inter_brems_M(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_virt(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_Mol(InterpolPoints, InterpolType);

// Arrays for interpolation:
double mys[10000];
double rep[10000];
double xx[InterpolPoints];
double y_brems_M[InterpolPoints];
double y_virt[InterpolPoints];
double y_Mol[InterpolPoints];


//----------------------------------------------------------------------------------------------
// Vertex correction:
inline double d_vertex()
  {
  return (2.*alpha/Pi)*(3.*Log(-qq/Pow2(m))/2. - 2.); // NB: twice larger than for ep scattering
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Functions for numerical integration (soft-photon bremsstrahlung):
double f_l1i_l1f(double x) // For the calculation of B(l1i, l1f, E_g_cut)
  {
  p_x = x*l1i + (1. - x)*l1f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_l1i_l2i(double x) // For the calculation of B(l1i, l2i, E_g_cut)
  {
  p_x = x*l1i + (1. - x)*l2i;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_l1i_l2f(double x) // For the calculation of B(l1i, l2f, E_g_cut)
  {
  p_x = x*l1i + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_l1f_l2i(double x) // For the calculation of B(l1f, l2i, E_g_cut)
  {
  p_x = x*l1f + (1. - x)*l2i;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_l1f_l2f(double x) // For the calculation of B(l1f, l2f, E_g_cut)
  {
  p_x = x*l1f + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_l2i_l2f(double x) // For the calculation of B(l2i, l2f, E_g_cut)
  {
  p_x = x*l2i + (1. - x)*l2f;
  return (Log(4.*Pow2(E_g_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Integration of the soft-photon bremsstrahlung cross section:
inline double d_l1i_l1i() // B(l1i, l1i, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/m) + E_li*Log(m/(E_li + Sqrt(Pow2(E_li) - Pow2(m))))/Sqrt(Pow2(E_li) - Pow2(m)))/Pi;
  }

inline double d_l1i_l1f() // B(l1i, l1f, E_g_cut)
  {
  return (l1i*l1f)*i_l1i_l1f.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l1f_l1f() // B(l1f, l1f, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/m) + E_1*Log(m/(E_1 + Sqrt(Pow2(E_1) - Pow2(m))))/Sqrt(Pow2(E_1) - Pow2(m)))/Pi;
  }

inline double d_l1i_l2i() // B(l1i, l2i, E_g_cut)
  {
  return (l1i*l2i)*i_l1i_l2i.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l1i_l2f() // B(l1i, l2f, E_g_cut)
  {
  return (l1i*l2f)*i_l1i_l2f.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l1f_l2i() // B(l1f, l2i, E_g_cut)
  {
  return (l1f*l2i)*i_l1f_l2i.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l1f_l2f() // B(l1f, l2f, E_g_cut)
  {
  return (l1f*l2f)*i_l1f_l2f.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l2i_l2i() // B(l2i, l2i, E_g_cut)
  {
  return (Log(2.*E_g_cut/m) - 1.)/(2.*Pi);
  }

inline double d_l2i_l2f() // B(l2i, l2f, E_g_cut)
  {
  return (l2i*l2f)*i_l2i_l2f.Integral(0., 1.)/(4.*Pi);
  }

inline double d_l2f_l2f() // B(l2f, l2f, E_g_cut)
  {
  return 0.5*(Log(2.*E_g_cut/m) + E_2*Log(m/(E_2 + Sqrt(Pow2(E_2) - Pow2(m))))/Sqrt(Pow2(E_2) - Pow2(m)))/Pi;
  }

// Total bremsstrahlung correction, Moller scattering:
inline double d_brems_M()
  {
  return -2.*alpha*(d_l1i_l1i() - 2.*d_l1i_l1f() + d_l1f_l1f() + d_l2i_l2i() - 2.*d_l2i_l2f() + d_l2f_l2f() + 2.*d_l1i_l2i() - 2.*d_l1i_l2f() - 2.*d_l1f_l2i() + 2.*d_l1f_l2f());
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function MCS (Moller cross section in the lowest order):
double MCS()
  {
  // Energy of the scattered lepton l1'
  E_1 = m*(E_li + m + (E_li - m)*Pow2(Cos(theta_1)))/(E_li + m - (E_li - m)*Pow2(Cos(theta_1)));

  // Four-momenta l1' and l2':
  l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m))*Cos(theta_1), E_1);
  l2f = l1i + l2i - l1f;

  // Full energy of the lepton l2':
  E_2 = l2f.E();

  // Mandelstam variables t and u:
  t = (l1f - l1i)*(l1f - l1i);
  u = (l1f - l2i)*(l1f - l2i);

  // Adir, Aex and Aint for Moller scattering:
  Adir = Pow2(s - 2.*m*m) + Pow2(u - 2.*m*m) + 4.*m*m*t;
  Aex = Pow2(s - 2.*m*m) + Pow2(t - 2.*m*m) + 4.*m*m*u;
  Aint = -(s - 2.*m*m)*(s - 6.*m*m);

  return 2.*Pow2(alpha)*(Cos(theta_1)/Pow2(E_li + m - (E_li - m)*Pow2(Cos(theta_1))))*(Adir/Pow2(t) + Aex/Pow2(u) - 2.*Aint/(t*u));
  }

double f_mol (double theta)
  {
  theta_1 = theta;
  return MCS()*(1. + inter_virt.Eval(theta_1) + inter_brems_M.Eval(theta_1))*Sin(theta_1);
  }
  
ROOT::Math::Functor1D func_mol(&f_mol);
ROOT::Math::GSLIntegrator i_mol(IntOpt);
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1 (1D, Moller elastic scattering):
class TFDISTR1: public TFoamIntegrand {
public:
  TFDISTR1(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Theta angle:
  theta_1 = theta_min + (theta_max - theta_min)*arg[0];

  return 2.*Pi*(theta_max - theta_min)*(1. + inter_virt.Eval(theta_1) + inter_brems_M.Eval(theta_1))*inter_Mol.Eval(theta_1)*Sin(theta_1)*mkb;
}
};
// The end of TFDISTR1
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2 (4D, Moller bremsstrahlung in the soft-photon approximation):
class TFDISTR2: public TFoamIntegrand {
public:
  TFDISTR2(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Four arguments (basic kinematic variables):
  theta_g = arg[0]*Pi; // Theta angle for the photon
  phi_g = arg[1]*2.*Pi; // Phi angle for the photon
  E_g = E_g_cut + arg[2]*(E_g_max - E_g_cut); // Energy for the photon
  theta_1 = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the lepton

  if (E_g > m*(E_li - m)/(m + E_li - Sqrt(Pow2(E_li) - Pow2(m))*Cos(theta_g))) return 0.;
  
  E_1 = m*(E_li + m + (E_li - m)*Pow2(Cos(theta_1)))/(E_li + m - (E_li - m)*Pow2(Cos(theta_1)));
  if (E_1 < m || E_1 > E_li) return 0.;

  l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m))*Cos(theta_1), E_1);
  vkf.SetPxPyPzE(E_g*Sin(theta_g)*Cos(phi_g), E_g*Sin(theta_g)*Sin(phi_g), E_g*Cos(theta_g), E_g);

  l2f = l1i + l2i - l1f;
  if (abs(l2f*l2f - m*m) > 1.e-5) cout << endl << "Bad kinematics! #1: " << l2f*l2f - m*m << endl;

  // Soft-photon bremsstrahlung factor:
  brems_factor = (-alpha*E_g/(4.*Pi*Pi))*((((1./(vkf*l1i))*l1i + (1./(vkf*l2i))*l2i - (1./(vkf*l1f))*l1f - (1./(vkf*l2f))*l2f)*((1./(vkf*l1i))*l1i + (1./(vkf*l2i))*l2i - (1./(vkf*l1f))*l1f - (1./(vkf*l2f))*l2f)));

  cospsi = Cos(theta_1)*Cos(theta_g) + Sin(theta_1)*Sin(theta_g)*Cos(phi_g);

  // More realistic kinematics to output:
  E_1 = (m*(E_li - E_g) - E_li*E_g*(1. - Cos(theta_g)))/(m + E_li*(1. - Cos(theta_1)) - E_g*(1. - cospsi));
  if (E_1 < m || E_1 > E_li - E_g) return 0.;

  l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m))*Cos(theta_1), E_1);

  l2f = l1i + l2i - l1f - vkf;
  if (abs(l2f*l2f - m*m) > 1.e-2) cout << endl << "Bad kinematics! #2: " << l2f*l2f - m*m << endl;

  return (2*Pi)*(2.*Pi*Pi)*(E_g_max - E_g_cut)*(theta_max - theta_min)*brems_factor*inter_Mol.Eval(theta_1)*Sin(theta_1)*Sin(theta_g)*mkb;
}
};
// The end of TFDISTR2
//==============================================================================================


int main(int argc, char **argv)
{

// Initial dialog:
cout << endl << "Full energy of the incident lepton (MeV): " << flush;
cin.getline(mychar, 64);
E_li = 0.001*atof(mychar);

cout << "Minimum polar angle of the lepton #1 (degree): " << flush;
cin.getline(mychar, 64);
theta_min = atof(mychar)*degrad;

cout << "Maximum polar angle of the lepton #1 (degree): " << flush;
cin.getline(mychar, 64);
theta_max = atof(mychar)*degrad;

cout << "Cut-off energy of bremsstrahlung photons (MeV): " << flush;
cin.getline(mychar, 64);
E_g_cut = 0.001*atof(mychar);

cout << "Maximum energy of bremsstrahlung photons (MeV): " << flush;
cin.getline(mychar, 64);
E_g_max = 0.001*atof(mychar);

cout << "Number of events to generate: " << flush;
cin.getline(mychar, 64);
nevents = atol(mychar);

time(&starttime); // Starting the timer

// Set seed to the random generator generator:
PseRan->SetSeed(0);

// Four momenta of the initial particles:
l1i.SetPxPyPzE(0., 0., Sqrt(E_li*E_li - m*m), E_li); // Incident electron/positron
l2i.SetPxPyPzE(0., 0. , 0., m); // Target electron

// Mandelstam variable s:
s = (l1i + l2i)*(l1i + l2i);

// Solid angle (steradian):
omega = 2.*Pi*(Cos(theta_min) - Cos(theta_max));

//----------------------------------------------------------------------------------------------
// Setting the functions and tolerance for numerical integration:
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

//----------------------------------------------------------------------------------------------
// Interpolation of bremsstrahlung and virtual-photon corrections:
int npoints = 0;
char str[128];

FILE *fvpol = fopen("vpol.dat", "r"); // Opening the file "vpol.dat"
if (fvpol == NULL) { cout << "Can't open file \"vpol.dat\"!" << endl; return EXIT_FAILURE; }

while (!feof(fvpol)) // Reading from the file "vpol.dat"
  {
  str[0] = 0;
  fgets(str, 128, fvpol);
  if (feof(fvpol) || strlen(str) == 0) break; // The end or empty string
    
  if (str[0] != '/')
    {
    sscanf(str, "%lf %lf", &mys[npoints], &rep[npoints]);
    rep[npoints] = 2.*rep[npoints];
    npoints++;
    }
  }
    
fclose(fvpol); // Closing the file "vpol.dat"
  
inter_vpol.SetData(npoints, mys, rep); // Interpolating

for (i = 0; i < InterpolPoints; i++)
  {
  // Theta angle for the lepton 1:
  theta_1 = theta_min + i*(theta_max - theta_min)/(InterpolPoints - 1);

  // Energy of the scattered electron l1'
  E_1 = m*(E_li + m + (E_li - m)*Pow2(Cos(theta_1)))/(E_li + m - (E_li - m)*Pow2(Cos(theta_1)));

  // Four-momenta l1' and l2':
  l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m))*Cos(theta_1), E_1);
  l2f = l1i + l2i - l1f;

  E_2 = l2f.E();

  // Four-momentum transfer squared:
  qq = (l1i - l1f)*(l1i - l1f);

  xx[i] = theta_1;            // Theta angle for the lepton
  y_brems_M[i] = d_brems_M(); // Bremsstrahlung correction from the lepton term

  // Moller elastic cross section:
  y_Mol[i] = MCS();
  
  // Virtual-photon correction:
  y_virt[i] = inter_vpol.Eval(-qq) + d_vertex();
  }

// Interpolating functions:
inter_brems_M.SetData(InterpolPoints, xx, y_brems_M); // Moller bremsstrahlung
inter_virt.SetData(InterpolPoints, xx, y_virt);       // Virtual-photon correction
inter_Mol.SetData(InterpolPoints, xx, y_Mol);         // Moller elastic cross section
//==============================================================================================


// Text file for Moller scattering events:
FILE *fdat = fopen ("events.dat", "w");
if (fdat == NULL) { cout << "Error opening file events.dat!" << endl; return EXIT_FAILURE; }

// ROOT file and ntuple for Moller scattering events:
TFile *froot = new TFile("events.root", "RECREATE");
if (froot->IsZombie()) { cout << "Error opening file events.root!" << endl; return EXIT_FAILURE; }
TNtuple *ntp = new TNtuple("ntp", "Moller ntuple", "E_1:theta_1:phi_1:E_2:theta_2:phi_2:E_g:theta_g:phi_g");

// Foam simulators
TFoam *Foam1 = new TFoam("Foam1"); // Moller elastic scattering
TFoam *Foam2 = new TFoam("Foam2"); // Moller bremsstrahlung

// Distribution functions:
TFoamIntegrand *Rho1 = new TFDISTR1(); // Moller elastic scattering
TFoamIntegrand *Rho2 = new TFDISTR2(); // Moller bremsstrahlung

// Initialization of the Foam1 simulator (Moller elastic scattering):
cout << endl << "Initialization of Foam1 (Moller elastic scattering):" << endl;
Foam1->SetkDim(1);        // Set number of dimensions
Foam1->SetnCells(1000);   // Set number of cells
Foam1->SetnSampl(200);    // Set number os samples
Foam1->SetOptRej(1);      // Unweighted events in MC generation
Foam1->SetRho(Rho1);      // Set distribution function
Foam1->SetPseRan(PseRan); // Set random number generator
Foam1->SetChat(1);        // Set "chat level" in the standard output
Foam1->Initialize();      // Initialization

// Initialization of the Foam2 simulator (Moller bremsstrahlung):
cout << endl << "Initialization of Foam2 (Moller bremsstrahlung):" << endl;
Foam2->SetkDim(4);        // Set number of dimensions
Foam2->SetnCells(30000);   // Set number of cells
Foam2->SetnSampl(1500);    // Set number os samples
Foam2->SetOptRej(1);      // Unweighted events in MC generation
Foam2->SetRho(Rho2);      // Set distribution function
Foam2->SetPseRan(PseRan); // Set random number generator
Foam2->SetChat(1);        // Set "chat level" in the standard output
Foam2->Initialize();      // Initialization

// Integrated cross sections:
double mel; // Moller elastic scattering
double mbr, mbr_error; // Moller bremsstrahlung


//----------------------------------------------------------------------------------------------
// Initialization:
i = 0; cout << endl;
for (n = 0; n < ninit; n++)
  {
  // Moller bremsstrahlung:
  Foam2->MakeEvent();

  // Progress bar:
  if (n == int (0.01*i*(ninit - 1)))
    {
    cout << "Initialization: " << i << "%" << '\r' << flush;
    i += 1;
    }
  }
//==============================================================================================


// Cross section for Moller bremsstrahlung:
Foam2->GetIntegMC(mbr, mbr_error);

// Cross section for elastic part:
mel = 2.*Pi*mkb*i_mol.Integral(theta_min, theta_max);

// Numbers of events to generate:
nev_el = nevents*mel/(mel + mbr); // Moller elastic events


//----------------------------------------------------------------------------------------------
// Event generation:
i = 0; cout << endl;
for (n = 0; n < nevents; n++)
  {
  // Elastic Moller scattering:
  if ((nev_el - count_nev_el) > 0 && PseRan->Rndm() < 1.*(nev_el - count_nev_el)/(nevents - n))
    {
    Foam1->MakeEvent();
    count_nev_el++;

    // Energy of the scattered lepton l1'
    E_1 = m*(E_li + m + (E_li - m)*Pow2(Cos(theta_1)))/(E_li + m - (E_li - m)*Pow2(Cos(theta_1)));

    // Four-momenta l1' and l2':
    l1f.SetPxPyPzE(Sqrt(Pow2(E_1) - Pow2(m))*Sin(theta_1), 0., Sqrt(Pow2(E_1) - Pow2(m))*Cos(theta_1), E_1);
    l2f = l1i + l2i - l1f;

    // Generating azimuthal angle (phi) for the lepton 1:
    phi_1 = 2.*Pi*(PseRan->Rndm());

    // Rotation the final particle four-momenta around the z-axis:
    l1f.RotateZ(phi_1); // Four-momentum of the final lepton 1
    l2f.RotateZ(phi_1); // Four-momentum of the final lepton 2

    E_g = 0.; theta_g = 0.; phi_g = 0.; // Purely elastic events
    }
  else // Moller bremsstrahlung:
    {
    Foam2->MakeEvent();

    // Generating azimuthal angle (phi) for the lepton 1:
    phi_1 = 2.*Pi*(PseRan->Rndm());

    // Rotation the final particle four-momenta around the z-axis:
    l1f.RotateZ(phi_1); // Four-momentum of the final lepton 1
    l2f.RotateZ(phi_1); // Four-momentum of the final lepton 2
    vkf.RotateZ(phi_1); // Four-momentum of the bremsstrahlung photon

    E_g = 1000.*vkf.E(); theta_g = vkf.Theta(); phi_g = vkf.Phi();
    }

  // Kinematic parameters to write to the output file:
  E_1 = 1000.*l1f.E(); theta_1 = l1f.Theta(); phi_1 = l1f.Phi();
  E_2 = 1000.*l2f.E(); theta_2 = l2f.Theta(); phi_2 = l2f.Phi();

  ntp->Fill(E_1, theta_1, phi_1, E_2, theta_2, phi_2, E_g, theta_g, phi_g);

  fprintf(fdat, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", E_1, theta_1, phi_1, E_2, theta_2, phi_2, E_g, theta_g, phi_g);

  // Progress bar:
  if (n == int (0.01*i*(nevents - 1)))
    {
    cout << "Generation of " << nevents << " events: " << i << "%" << '\r' << flush;
    i += 1;
    }
  }
cout << endl;
//==============================================================================================


// Finalization of simulators:
double x1, x2;
cout << endl << "Finalization of Foam1 (elastic Moller scattering):" << endl;
Foam1->Finalize(x1, x2);
cout << endl << "Finalization of Foam2 (Moller bremsstrahlung):" << endl;
Foam2->Finalize(x1, x2);

// Output files:
froot->cd();
cout << endl;
ntp->Print();
froot->Write();
froot->Close();
fclose(fdat);

time (&stoptime); // Stopping the timer

// Some information useful to know:

cout << endl << "DIFFERENTIAL CROSS SECTION: " << 0.001*(mel + mbr)/omega << " millibarn / steradian" << endl;

cout << endl << "INTEGRATED LUMINOSITY: " << nevents/(mel + mbr) << " inverse microbarn" << endl;

cout << endl << "BEAMTIME (assuming the target thickness of 10^18 at/cm^2" << endl;
cout << "and the beam current of 10 nA): " << nevents/(6.25e-2*(mel + mbr)) << " seconds" << endl;

cout << endl << "It took " << stoptime - starttime << " seconds to generate " << nevents << " events. Bye!" << endl << endl;

return EXIT_SUCCESS;
}

