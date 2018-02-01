//----------------------------------------------------------------------------------------------
// Program EDS (radiative electron-deuteron scattering in the soft-photon approximation)
// http://www.inp.nsk.su/~gramolin/eds/
//
// (c) Alexander Gramolin, July 2017
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
const double m_e = 0.510998928e-3;    // Electron/positron mass (in GeV)
const double M = 1.875612859;         // Deuteron mass (in GeV)
const double alpha = 1/137.035999074; // Fine-structure constant
const double mkb = 389.379338;        // GeV^{-2} to mkbarn conversion


//----------------------------------------------------------------------------------------------
// Some mathematical functions:
double Log(double arg) // Logarithm
  {
  return TMath::Log(arg);
  }

double Sqrt(double arg) // Square root
  {
  return TMath::Sqrt(arg);
  }

double Pow2(double arg) // Pow2
  {
  return arg*arg;
  }

double Pow4(double arg) // Pow4
  {
  return arg*arg*arg*arg;
  }
  
double Power(double arg1, double arg2) // Any power
  {
  return TMath::Power(arg1, arg2);
  }

double Sin(double arg) // Sine
  {
  return TMath::Sin(arg);
  }

double Cos(double arg) // Cosine
  {
  return TMath::Cos(arg);
  }

double Tan(double arg) // Tangent
  {
  return TMath::Tan(arg);
  }

double Abs(double arg) // Absolute value
  {
  return TMath::Abs(arg);
  }
//==============================================================================================

  
//----------------------------------------------------------------------------------------------
// Deuteron form factors:
// (we use Parameterization I from Eur. Phys. J A 7, 421 (2000),
// see http://irfu.cea.fr/Sphn/T20/Parametrisations/)

const double Gc0 = 1., Gq0 = 25.83, Gm0 = 1.714; // Static moments
const double Qc0 = 4.21, Qq0 = 8.1, Qm0 = 7.37; // Nodes, in fm^{-1}
const double ac1 = 6.740e-1, aq1 = 8.796e-1, am1 = 5.804e-1;
const double ac2 = 2.246e-2, aq2 = -5.656e-2, am2 = 8.701e-2;
const double ac3 = 9.806e-3, aq3 = 1.933e-2, am3 = -3.624e-3;
const double ac4 = -2.709e-4, aq4 = -6.734e-4, am4 = 3.448e-4;
const double ac5 = 3.793e-6, aq5 = 9.438e-6, am5 = -2.818e-6;

double Gc(double QQ) // Charge form Factor
  {
  // Conversion from GeV to fm^{-1}:
  double Qfm = Sqrt(QQ)/0.1973269718;
  
  return Gc0*(1. - Pow2(Qfm/Qc0))/(1. + ac1*Pow2(Qfm) + ac2*Pow4(Qfm) + ac3*Power(Qfm, 6) + ac4*Power(Qfm, 8) + ac5*Power(Qfm, 10));
  }

double Gq(double QQ) // Quadrupole form factor
  {
  // Conversion from GeV to fm^{-1}:
  double Qfm = Sqrt(QQ)/0.1973269718;
  
  return Gq0*(1. - Pow2(Qfm/Qq0))/(1. + aq1*Pow2(Qfm) + aq2*Pow4(Qfm) + aq3*Power(Qfm, 6) + aq4*Power(Qfm, 8) + aq5*Power(Qfm, 10));
  }

double Gm(double QQ) // Magnetic form factor
  {
  // Conversion from GeV to fm^{-1}:
  double Qfm = Sqrt(QQ)/0.1973269718;
  
  return Gm0*(1. - Pow2(Qfm/Qm0))/(1. + am1*Pow2(Qfm) + am2*Pow4(Qfm) + am3*Power(Qfm, 6) + am4*Power(Qfm, 8) + am5*Power(Qfm, 10));
  }

double A(double QQ) // Structure function A
  {
  double tau = QQ/Pow2(2.*M);
  
  return Pow2(Gc(QQ)) + 8.*Pow2(tau*Gq(QQ))/9. + 2.*tau*Pow2(Gm(QQ))/3.;
  }

double B(double QQ) // Structure function B
  {
  double tau = QQ/Pow2(2.*M);
  
  return 4.*tau*(1. + tau)*Pow2(Gm(QQ))/3.;
  }
//==============================================================================================


// Four-momenta of the initial particles:
TLorentzVector p1; // Incident electron
TLorentzVector p2; // Target deuteron

// Four-momenta of the final particles:
TLorentzVector p3; // Scattered electron
TLorentzVector p4; // Recoil deuteron
TLorentzVector p5; // Bremsstrahlung photon

TLorentzVector p_x;

double cospsi, omega, QQ, tau;

double E_1; // Full energy of the incident electron (in GeV)

double theta_min, theta_max;

double E_5_cut, E_5_max;

int i;
long n;
long nevents; // Total number of events to generate
long nev_el, count_nev_el = 0;

// Kinematic variables:
double E_3, theta_3, phi_3; // For the scattered electron
double E_4, theta_4, phi_4; // For the recoil deuteron
double E_5, theta_5, phi_5; // For the bremsstrahlung photon

double brems_factor; // To get the soft-photon bremsstrahlung cross section

// The random number generator:
TRandom3 *PseRan = new TRandom3();

time_t starttime, stoptime; // For the timer

char mychar[64];


// Some functions: -----------------------------------------------------------------------------
double f_p1_p2(double x);
double f_p1_p3(double x);
double f_p1_p4(double x);
double f_p2_p3(double x);
double f_p2_p4(double x);
double f_p3_p4(double x);

ROOT::Math::Functor1D func_p1_p2(&f_p1_p2);
ROOT::Math::Functor1D func_p1_p3(&f_p1_p3);
ROOT::Math::Functor1D func_p1_p4(&f_p1_p4);
ROOT::Math::Functor1D func_p2_p3(&f_p2_p3);
ROOT::Math::Functor1D func_p2_p4(&f_p2_p4);
ROOT::Math::Functor1D func_p3_p4(&f_p3_p4);

// Integrators: --------------------------------------------------------------------------------
ROOT::Math::GSLIntegrator i_p1_p2(IntOpt);
ROOT::Math::GSLIntegrator i_p1_p3(IntOpt);
ROOT::Math::GSLIntegrator i_p1_p4(IntOpt);
ROOT::Math::GSLIntegrator i_p2_p3(IntOpt);
ROOT::Math::GSLIntegrator i_p2_p4(IntOpt);
ROOT::Math::GSLIntegrator i_p3_p4(IntOpt);

// Interpolators: ------------------------------------------------------------------------------
ROOT::Math::Interpolator inter_vpol(10000, InterpolType);
ROOT::Math::Interpolator inter_brems(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_virt(InterpolPoints, InterpolType);
ROOT::Math::Interpolator inter_ecs(InterpolPoints, InterpolType);

// Arrays for interpolation:
double mys[10000];
double rep[10000];
double xx[InterpolPoints];
double y_brems[InterpolPoints];
double y_virt[InterpolPoints];
double y_ecs[InterpolPoints];


//----------------------------------------------------------------------------------------------
// Vertex correction:
inline double d_vertex()
  {
  return (alpha/Pi)*(3.*Log(QQ/Pow2(m_e))/2. - 2.);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Functions for numerical integration (soft-photon bremsstrahlung):
double f_p1_p2(double x) // For the calculation of B(p1, p2, E_5_cut)
  {
  p_x = x*p1 + (1. - x)*p2;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_p1_p3(double x) // For the calculation of B(p1, p3, E_5_cut)
  {
  p_x = x*p1 + (1. - x)*p3;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_p1_p4(double x) // For the calculation of B(p1, p4, E_5_cut)
  {
  p_x = x*p1 + (1. - x)*p4;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_p2_p3(double x) // For the calculation of B(p2, p3, E_5_cut)
  {
  p_x = x*p2 + (1. - x)*p3;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_p2_p4(double x) // For the calculation of B(p2, p4, E_5_cut)
  {
  p_x = x*p2 + (1. - x)*p4;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }

double f_p3_p4(double x) // For the calculation of B(p3, p4, E_5_cut)
  {
  p_x = x*p3 + (1. - x)*p4;
  return (Log(4.*Pow2(E_5_cut)/(p_x*p_x)) + (p_x.E())*Log((p_x.E() - (p_x.Vect()).Mag())/(p_x.E() + (p_x.Vect()).Mag()))/((p_x.Vect()).Mag()))/(p_x*p_x);
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Integration of the soft-photon bremsstrahlung cross section:
inline double d_p1_p1() // B(p1, p1, E_5_cut)
  {
  return 0.5*(Log(2.*E_5_cut/m_e) + E_1*Log(m_e/(E_1 + Sqrt(Pow2(E_1) - Pow2(m_e))))/Sqrt(Pow2(E_1) - Pow2(m_e)))/Pi;
  }

inline double d_p1_p2() // B(p1, p2, E_5_cut)
  {
  return (p1*p2)*i_p1_p2.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p1_p3() // B(p1, p3, E_5_cut)
  {
  return (p1*p3)*i_p1_p3.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p1_p4() // B(p1, p4, E_5_cut)
  {
  return (p1*p4)*i_p1_p4.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p2_p2() // B(p2, p2, E_5_cut)
  {
  return 0.5*(Log(2.*E_5_cut/M) - 1.)/Pi;
  }

inline double d_p2_p3() // B(p2, p3, E_5_cut)
  {
  return (p2*p3)*i_p2_p3.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p2_p4() // B(p2, p4, E_5_cut)
  {
  return (p2*p4)*i_p2_p4.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p3_p3() // B(p3, p3, E_5_cut)
  {
  return 0.5*(Log(2.*E_5_cut/m_e) + E_3*Log(m_e/(E_3 + Sqrt(Pow2(E_3) - Pow2(m_e))))/Sqrt(Pow2(E_3) - Pow2(m_e)))/Pi;
  }

inline double d_p3_p4() // B(p3, p4, E_5_cut)
  {
  return (p3*p4)*i_p3_p4.Integral(0., 1.)/(4.*Pi);
  }

inline double d_p4_p4() // B(p4, p4, E_5_cut)
  {
  return 0.5*(Log(2.*E_5_cut/M) + E_4*Log(M/(E_4 + Sqrt(Pow2(E_4) - Pow2(M))))/Sqrt(Pow2(E_4) - Pow2(M)))/Pi;
  }

// Total bremsstrahlung correction:
inline double d_brems()
  {
  return -2.*alpha*(d_p1_p1() - 2.*d_p1_p2() - 2.*d_p1_p3() + 2.*d_p1_p4() + d_p2_p2() + 2.*d_p2_p3() - 2.*d_p2_p4() + d_p3_p3() - 2.*d_p3_p4() + d_p4_p4());
  }
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Function ECS (elastic cross section in the lowest order):
double ECS()
  {
  tau = QQ/Pow2(2.*M);
  
  return Pow2(0.5*alpha/E_1)*Pow2(Cos(0.5*theta_3)/Pow2(Sin(0.5*theta_3)))*(E_3/E_1)*(A(QQ) + B(QQ)*Pow2(Tan(0.5*theta_3)));
  }

double f_ecs(double arg)
  {
  theta_3 = arg;
  E_3 = M*E_1/(M + E_1*(1. - Cos(theta_3)));
  QQ = 2.*M*Pow2(E_1)*(1. - Cos(theta_3))/(M + E_1*(1. - Cos(theta_3)));
  return ECS()*(1. + inter_virt.Eval(theta_3) + inter_brems.Eval(theta_3))*Sin(theta_3);
  }
  
ROOT::Math::Functor1D func_ecs(&f_ecs);
ROOT::Math::GSLIntegrator i_ecs(IntOpt);
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR1 (1D, elastic scattering):
class TFDISTR1: public TFoamIntegrand {
public:
  TFDISTR1(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Theta angle:
  theta_3 = theta_min + (theta_max - theta_min)*arg[0];

  return 2.*Pi*(theta_max - theta_min)*(1. + inter_virt.Eval(theta_3) + inter_brems.Eval(theta_3))*inter_ecs.Eval(theta_3)*Sin(theta_3)*mkb;
}
};
// The end of TFDISTR1
//==============================================================================================


//----------------------------------------------------------------------------------------------
// Distribution function TFDISTR2 (4D, bremsstrahlung in the soft-photon approximation):
class TFDISTR2: public TFoamIntegrand {
public:
  TFDISTR2(){};
  Double_t Density(int nDim, Double_t *arg) {

  // Four arguments (basic kinematic variables):
  theta_5 = arg[0]*Pi; // Theta angle for the photon
  phi_5 = arg[1]*2.*Pi; // Phi angle for the photon
  E_5 = E_5_cut + arg[2]*(E_5_max - E_5_cut); // Energy for the photon
  theta_3 = theta_min + arg[3]*(theta_max - theta_min); // Theta angle for the electron
  
  if (E_5 > M*(E_1 - m_e)/(M + E_1 - Sqrt(Pow2(E_1) - Pow2(m_e))*Cos(theta_5))) return 0.;
  
  cospsi = Cos(theta_3)*Cos(theta_5) + Sin(theta_3)*Sin(theta_5)*Cos(phi_5);
  
  E_3 = (M*(E_1 - E_5) - E_1*E_5*(1. - Cos(theta_5)))/(M + E_1*(1. - Cos(theta_3)) - E_5*(1. - cospsi));
  if (E_3 < m_e || E_3 > E_1 - E_5) return 0.;
  
  p3.SetPxPyPzE(Sqrt(Pow2(E_3) - Pow2(m_e))*Sin(theta_3), 0., Sqrt(Pow2(E_3) - Pow2(m_e))*Cos(theta_3), E_3);
  p5.SetPxPyPzE(E_5*Sin(theta_5)*Cos(phi_5), E_5*Sin(theta_5)*Sin(phi_5), E_5*Cos(theta_5), E_5);
  p4 = p1 + p2 - p3 - p5;
  
  if (Abs(p4*p4 - M*M) > 0.001) cout << endl << "Bad kinematics: " << p4*p4 - M*M << endl;
  
  // Soft-photon bremsstrahlung factor:
  brems_factor = (-alpha*E_5/Pow2(2.*Pi))*(-(1./(p1*p5))*p1 + (1./(p2*p5))*p2 + (1./(p3*p5))*p3 - (1./(p4*p5))*p4)*(-(1./(p1*p5))*p1 + (1./(p2*p5))*p2 + (1./(p3*p5))*p3 - (1./(p4*p5))*p4);
  
  return (2*Pi)*(2.*Pi*Pi)*(E_5_max - E_5_cut)*(theta_max - theta_min)*brems_factor*inter_ecs.Eval(theta_3)*Sin(theta_3)*Sin(theta_5)*mkb;
}
};
// The end of TFDISTR2
//==============================================================================================


int main(int argc, char **argv)
{
cout << endl << "Radiative electron-deuteron scattering in the soft-photon approximation." << endl;

// Initial dialog:
cout << endl << "Full energy of the incident electron (MeV): " << flush;
cin.getline(mychar, 64);
E_1 = 0.001*atof(mychar);

cout << "Minimum polar angle of the scattered electron (degree): " << flush;
cin.getline(mychar, 64);
theta_min = atof(mychar)*degrad;

cout << "Maximum polar angle of the scattered electron (degree): " << flush;
cin.getline(mychar, 64);
theta_max = atof(mychar)*degrad;

cout << "Cut-off energy of bremsstrahlung photons (MeV): " << flush;
cin.getline(mychar, 64);
E_5_cut = 0.001*atof(mychar);

cout << "Maximum energy of bremsstrahlung photons (MeV): " << flush;
cin.getline(mychar, 64);
E_5_max = 0.001*atof(mychar);

cout << "Number of events to generate: " << flush;
cin.getline(mychar, 64);
nevents = atol(mychar);

time(&starttime); // Starting the timer

// Set seed to the random generator generator:
PseRan->SetSeed(0);

// Four-momenta of the initial particles:
p1.SetPxPyPzE(0., 0., Sqrt(Pow2(E_1) - Pow2(m_e)), E_1); // Incident electron
p2.SetPxPyPzE(0., 0., 0., M); // Target deuteron

// Solid angle (steradian):
omega = 2.*Pi*(Cos(theta_min) - Cos(theta_max));


//----------------------------------------------------------------------------------------------
// Setting the functions and tolerance for numerical integration:
i_p1_p2.SetFunction(func_p1_p2); // To calculate B(p1, p2, E_5_cut)
i_p1_p2.SetRelTolerance(IntTol);

i_p1_p3.SetFunction(func_p1_p3); // To calculate B(p1, p3, E_5_cut)
i_p1_p3.SetRelTolerance(IntTol);

i_p1_p4.SetFunction(func_p1_p4); // To calculate B(p1, p4, E_5_cut)
i_p1_p4.SetRelTolerance(IntTol);

i_p2_p3.SetFunction(func_p2_p3); // To calculate B(p2, p3, E_5_cut)
i_p2_p3.SetRelTolerance(IntTol);

i_p2_p4.SetFunction(func_p2_p4); // To calculate B(p2, p4, E_5_cut)
i_p2_p4.SetRelTolerance(IntTol);

i_p3_p4.SetFunction(func_p3_p4); // To calculate B(p3, p4, E_5_cut)
i_p3_p4.SetRelTolerance(IntTol);

i_ecs.SetFunction(func_ecs);
i_ecs.SetRelTolerance(IntTol);


//----------------------------------------------------------------------------------------------
// Interpolation of the bremsstrahlung and virtual-photon corrections:
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
  // Electron scattering angle:
  theta_3 = theta_min + i*(theta_max - theta_min)/(InterpolPoints - 1);

  // Full energy of the scattered electron:
  E_3 = M*E_1/(M + E_1*(1. - Cos(theta_3)));

  // Four-momenta p3 and p4:
  p3.SetPxPyPzE(Sqrt(Pow2(E_3) - Pow2(m_e))*Sin(theta_3), 0., Sqrt(Pow2(E_3) - Pow2(m_e))*Cos(theta_3), E_3);
  p4 = p1 + p2 - p3;

  E_4 = p4.E();

  // Four-momentum transfer squared:
  QQ = 2.*M*Pow2(E_1)*(1. - Cos(theta_3))/(M + E_1*(1. - Cos(theta_3)));

  xx[i] = theta_3;        // Electron scattering angle
  y_brems[i] = d_brems(); // Bremsstrahlung

  // Elastic scattering cross section:
  y_ecs[i] = ECS();
  
  // Virtual-photon correction:
  y_virt[i] = inter_vpol.Eval(QQ) + d_vertex();
  }

// Interpolating functions:
inter_brems.SetData(InterpolPoints, xx, y_brems); // Bremsstrahlung correction
inter_virt.SetData(InterpolPoints, xx, y_virt);   // Virtual-photon correction
inter_ecs.SetData(InterpolPoints, xx, y_ecs);     // Elastic scattering cross section
//==============================================================================================


// Text file to output events:
FILE *fdat = fopen ("events.dat", "w");
if (fdat == NULL) { cout << "Error opening file events.dat!" << endl; return EXIT_FAILURE; }

// ROOT file to output events:
TFile *froot = new TFile("events.root", "RECREATE");
if (froot->IsZombie()) { cout << "Error opening file events.root!" << endl; return EXIT_FAILURE; }
TNtuple *ntp = new TNtuple("ntp", "ed scattering", "E_e:theta_e:phi_e:E_d:theta_d:phi_d:E_g:theta_g:phi_g");

// Foam simulators
TFoam *Foam1 = new TFoam("Foam1"); // Elastic scattering
TFoam *Foam2 = new TFoam("Foam2"); // Bremsstrahlung

// Distribution functions:
TFoamIntegrand *Rho1 = new TFDISTR1(); // Elastic scattering
TFoamIntegrand *Rho2 = new TFDISTR2(); // Bremsstrahlung

// Initialization of the Foam1 simulator (elastic scattering):
cout << endl << "Initialization of Foam1 (elastic scattering):" << endl;
Foam1->SetkDim(1);        // Set number of dimensions
Foam1->SetnCells(1000);   // Set number of cells
Foam1->SetnSampl(200);    // Set number os samples
Foam1->SetOptRej(1);      // Unweighted events in MC generation
Foam1->SetRho(Rho1);      // Set distribution function
Foam1->SetPseRan(PseRan); // Set random number generator
Foam1->SetChat(1);        // Set "chat level" in the standard output
Foam1->Initialize();      // Initialization

// Initialization of the Foam2 simulator (bremsstrahlung):
cout << endl << "Initialization of Foam2 (bremsstrahlung):" << endl;
Foam2->SetkDim(4);        // Set number of dimensions
Foam2->SetnCells(30000);   // Set number of cells
Foam2->SetnSampl(1500);    // Set number os samples
Foam2->SetOptRej(1);      // Unweighted events in MC generation
Foam2->SetRho(Rho2);      // Set distribution function
Foam2->SetPseRan(PseRan); // Set random number generator
Foam2->SetChat(1);        // Set "chat level" in the standard output
Foam2->Initialize();      // Initialization

// Integrated cross sections:
double el; // For elastic scattering
double br, br_error; // For bremsstrahlung


//----------------------------------------------------------------------------------------------
// Initialization:
i = 0; cout << endl;
for (n = 0; n < ninit; n++)
  {
  // Bremsstrahlung:
  Foam2->MakeEvent();

  // Progress bar:
  if (n == int (0.01*i*(ninit - 1)))
    {
    cout << "Initialization: " << i << "%" << '\r' << flush;
    i += 1;
    }
  }
//==============================================================================================


// Cross section for bremsstrahlung (microbarn):
Foam2->GetIntegMC(br, br_error);

// Cross section for the elastic part (microbarn):
el = 2.*Pi*mkb*i_ecs.Integral(theta_min, theta_max);

// Numbers of events to generate:
nev_el = nevents*el/(el + br); // Elastic events


//----------------------------------------------------------------------------------------------
// Event generation:
i = 0; cout << endl;
for (n = 0; n < nevents; n++)
  {
  // Elastic scattering:
  if ((nev_el - count_nev_el) > 0 && PseRan->Rndm() < 1.*(nev_el - count_nev_el)/(nevents - n))
    {
    Foam1->MakeEvent();
    count_nev_el++;

    // Full energy of the scattered electron:
    E_3 = M*E_1/(M + E_1*(1. - Cos(theta_3)));

    // Four-momenta p3 and p4:
    p3.SetPxPyPzE(Sqrt(Pow2(E_3) - Pow2(m_e))*Sin(theta_3), 0., Sqrt(Pow2(E_3) - Pow2(m_e))*Cos(theta_3), E_3);
    p4 = p1 + p2 - p3;

    // Generate the azimuthal angle (phi) for the scattered electron:
    phi_3 = 2.*Pi*(PseRan->Rndm());

    // Rotate the final particle four-momenta around the z-axis:
    p3.RotateZ(phi_3); // Four-momentum of the scattered electron
    p4.RotateZ(phi_3); // Four-momentum of the recoil deuteron

    E_5 = 0.; theta_5 = 0.; phi_5 = 0.; // Purely elastic events
    }
  else // Bremsstrahlung:
    {
    Foam2->MakeEvent();

    // Generate the azimuthal angle (phi) for the scattered electron:
    phi_3 = 2.*Pi*(PseRan->Rndm());

    // Rotate the final particle four-momenta around the z-axis:
    p3.RotateZ(phi_3); // Four-momentum of the scattered electron
    p4.RotateZ(phi_3); // Four-momentum of the recoil deuteron
    p5.RotateZ(phi_3); // Four-momentum of the bremsstrahlung photon

    E_5 = 1000.*p5.E(); theta_5 = p5.Theta(); phi_5 = p5.Phi();
    }

  // Kinematic parameters to write to the output file:
  E_3 = 1000.*p3.E(); theta_3 = p3.Theta(); phi_3 = p3.Phi();
  E_4 = 1000.*p4.E(); theta_4 = p4.Theta(); phi_4 = p4.Phi();

  ntp->Fill(E_3, theta_3, phi_3, E_4, theta_4, phi_4, E_5, theta_5, phi_5);

  fprintf(fdat, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", E_3, theta_3, phi_3, E_4, theta_4, phi_4, E_5, theta_5, phi_5);

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
cout << endl << "Finalization of Foam1 (elastic scattering):" << endl;
Foam1->Finalize(x1, x2);
cout << endl << "Finalization of Foam2 (bremsstrahlung):" << endl;
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

cout << endl << "DIFFERENTIAL CROSS SECTION: " << (el + br)/omega << " microbarn / steradian" << endl;

cout << endl << "INTEGRATED LUMINOSITY: " << nevents/(el + br) << " inverse microbarn" << endl;

cout << endl << "BEAMTIME (assuming the target thickness of 2*10^18 at/cm^2" << endl;
cout << "and the beam current of 10 nA): " << 0.5*nevents/(6.25e-2*(el + br)) << " seconds" << endl;

cout << endl << "It took " << stoptime - starttime << " seconds to generate " << nevents << " events. Bye!" << endl << endl;

return EXIT_SUCCESS;
}

