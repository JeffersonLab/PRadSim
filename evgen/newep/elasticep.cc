//
// elasticep.cxx
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom2.h"

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>

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

const double Pi = TMath::Pi();
const double deg = Pi / 180.0;
const double m = 0.51099893e-3;
const double m2 = TMath::Power(m, 2);
const double M = 0.938272046;
const double M2 = TMath::Power(M, 2);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern "C" {
    void mainf_(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

inline double Pow2(double arg)
{
    return TMath::Power(arg, 2);
}

double ElasticEnergy(double theta);

void RadCorEPXS(double q2);
double RadCorEPXS_Sin(double theta);

void Clear();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double E_li;

double theta_min, theta_max;
double phi_min, phi_max;

double E_lf, theta_l, phi_l;
double E_p, theta_p, phi_p;
double E_g, theta_g, phi_g;

double omega;

int randin, randout;

double ppx = 0.0, ppy = 0.0, ppz = 0.0, pp0 = 0.0, pp = 0.0;
double q2p = 0.0, phip = 0.0;
double epx = 0.0, epy = 0.0, epz = 0.0, ep0 = 0.0, pe = 0.0;
double phx = 0.0, phy = 0.0, phz = 0.0, ph0 = 0.0, pph = 0.0;
double sigmaborn = 0.0, sigmatot = 0.0, sigmarad = 0.0, sigmabsv = 0.0;
double q2 = 0.0, vmin = 0.0, phi = 0.0, ebeam = 0.0, weightsig = 0.0;
double vcut = 0.0, pvgen = 0.0, ptgen = 0.0;
double pichgen = 0.0;

double thetat[InterpolPoints];
double xs_sin[InterpolPoints];

TLorentzVector v_li, v_pi;
TLorentzVector v_lf, v_pf;

TRandom2 *PseRan = new TRandom2();

ROOT::Math::Functor1D Func_XS_Sin(&RadCorEPXS_Sin);
ROOT::Math::GSLIntegrator Integrator_XS_Sin(IntOpt);
ROOT::Math::Interpolator Interpolator_XS_Sin(InterpolPoints, InterpolType);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RadCorEPIntegrand: public TFoamIntegrand
{
public:
    RadCorEPIntegrand() {};

    Double_t Density(int nDim, Double_t *arg)
    {
        theta_l = theta_min + arg[0] * (theta_max - theta_min);

        return (phi_max - phi_min) * (theta_max - theta_min) * Interpolator_XS_Sin.Eval(theta_l);
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main()
{
    char mychar[64];

    std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    std::cin.getline(mychar, 64);
    E_li = 0.001 * atof(mychar); // MeV

    std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_min = atof(mychar) * deg;

    std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_max = atof(mychar) * deg;

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int nevents = atoi(mychar);

    //E_li = 1.1;
    //theta_min = 0.6 * deg;
    //theta_max = 6.0 * deg;
    //int nevents = 50000;

    PseRan->SetSeed(0);

    phi_min = -180.0 * deg;
    phi_max = 180.0 * deg;

    omega = (phi_max - phi_min) * (Cos(theta_min) - Cos(theta_max));

    ebeam = E_li;

    randin = 0;
    vmin = 0.000025335;
    vcut = 0.05; //0.0001;//0.05;

    v_li.SetPxPyPzE(0., 0., Sqrt(Pow2(E_li) - m2), E_li);
    v_pi.SetPxPyPzE(0., 0., 0., M);

    FILE *fp = fopen("xs.dat", "w");

    /*
    for (int i = 0; i < 101; i++) {
        double q2 = Exp((-6 + i * 0.06) * Log(10.0));
        double tt = 0;

        RadCorEPXS(q2);
        fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / Pi, q2, sigmaborn, sigmatot, sigmarad, sigmabsv, sigmatot / sigmaborn - 1, sigmarad / sigmaborn, sigmabsv / sigmaborn);
        printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / Pi, q2, sigmaborn, sigmatot, sigmarad, sigmabsv, sigmatot / sigmaborn - 1, sigmarad / sigmaborn, sigmabsv / sigmaborn);
    }

    fclose(fp);
    exit(0);
    */

    for (int i = 0; i < InterpolPoints; i++) {
        if (i % 1000 == 0 && i != 0) std::cout << i << std::endl;

        theta_l = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
        E_lf = ElasticEnergy(theta_l);
        double q2 = -2.0 * M * (E_lf - E_li);

        thetat[i] = theta_l;
        xs_sin[i] = RadCorEPXS_Sin(theta_l);
        fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_l * 180.0 / Pi, q2, sigmaborn, sigmatot, sigmarad, sigmabsv, sigmatot / sigmaborn - 1, sigmarad / sigmaborn, sigmabsv / sigmaborn);
        printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_l * 180.0 / Pi, q2, sigmaborn, sigmatot, sigmarad, sigmabsv, sigmatot / sigmaborn - 1, sigmarad / sigmaborn, sigmabsv / sigmaborn);
    }

    fclose(fp);

    Integrator_XS_Sin.SetFunction(Func_XS_Sin);
    Integrator_XS_Sin.SetRelTolerance(IntTol);

    Interpolator_XS_Sin.SetData(InterpolPoints, thetat, xs_sin);

    double xsint = (phi_max - phi_min) * Integrator_XS_Sin.Integral(theta_min, theta_max);

    TFoam *FoamX = new TFoam("FoamX");
    TFoamIntegrand *pIntegrand = new RadCorEPIntegrand();
    FoamX->SetkDim(1);
    FoamX->SetnCells(10000); // Set number of cells
    FoamX->SetnSampl(500); // Set number os samples
    FoamX->SetOptRej(1); // Unweighted events in MC generation
    FoamX->SetRho(pIntegrand); // Set distribution function
    FoamX->SetPseRan(PseRan); // Set random number generator
    //FoamX->SetChat(1); // Set "chat level" in the standard output
    FoamX->Initialize();

    fp = fopen("epelastic.dat", "w");

    for (int i = 0; i < nevents; ++i) {
        if (i % 1000 == 0 && i != 0) std::cout << i << std::endl;

        FoamX->MakeEvent();

        RadCorEPXS_Sin(theta_l);
        phi = phi_min + (phi_max - phi_min) * (PseRan->Rndm());

        if (ph0 > 0) {
            E_lf = ep0;
            theta_l = ACos(epz / ep0);
            phi_l = ATan2(epy, epx) + phi;

            if (phi_l > Pi) phi_l = phi_l - 2 * Pi;

            E_p = pp0;
            theta_p = ACos(ppz / pp0);
            phi_p = ATan2(ppy, ppx) + phi;

            if (phi_p > Pi) phi_p = phi_p - 2 * Pi;

            E_g = ph0;
            theta_g = ACos(phz / ph0);
            phi_g = ATan2(phy, phx) + phi;

            if (phi_g > Pi) phi_g = phi_g - 2 * Pi;
        } else {
            phi_l = phi;

            v_lf.SetPxPyPzE(Sqrt(Pow2(E_lf) - m2) * Sin(theta_l), 0., Sqrt(Pow2(E_lf) - m2) * Cos(theta_l), E_lf);
            v_pf = v_li + v_pi - v_lf;

            if (Abs(M2 - v_pf * v_pf) > 1.e-8)
                std::cout << "Warning: bad kinematics! M^2 - v_pf^2 = " << M2 - v_pf *v_pf << " GeV^2" << std::endl;

            E_p = v_pf.E();
            theta_p = v_pf.Theta();

            if (phi_l < 0.) phi_p = phi_l + Pi;
            else phi_p = phi_p = phi_l - Pi;

            E_g = 0.0;
            theta_g = 0.0;
            phi_g = 0.0;
        }

        fprintf(fp, "%9.3lf %8.6lf %8.5lf %9.3lf %8.6lf %8.5lf %9.3lf %8.6lf %8.5lf\n", 1000. * E_lf, theta_l, phi_l, 1000. * E_p, theta_p, phi_p, 1000. * E_g, theta_g, phi_g);
        //printf("%9.3lf %8.6lf %8.5lf %9.3lf %8.6lf %8.5lf %9.3lf %8.6lf %8.5lf\n", 1000. * E_lf, theta_l, phi_l, 1000. * E_p, theta_p, phi_p, 1000. * E_g, theta_g, phi_g);
    }

    fclose(fp);

    fp = fopen("epelastic.info", "w");

    fprintf(fp, "cross section (averaged over the solid angle):\n");
    fprintf(fp, "%lf microbarn / steradian\n", xsint / omega);
    fprintf(fp, "integrated luminosity:\n");
    fprintf(fp, "%lf inverse microbarn\n", nevents / xsint);

    printf("cross section (averaged over the solid angle):\n");
    printf("%lf microbarn / steradian\n", xsint / omega);
    printf("integrated luminosity:\n");
    printf("%lf inverse microbarn\n", nevents / xsint);

    fclose(fp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ElasticEnergy(double theta)
{
    return ((E_li + M) * (M * E_li + m2) + Sqrt(M2 - Pow2(m * Sin(theta))) * (Pow2(E_li) - m2) * Cos(theta)) / (Pow2(E_li + M) - (Pow2(E_li) - m2) * Pow2(Cos(theta)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double RadCorEPXS_Sin(double theta)
{
    Clear();

    theta_l = theta;

    E_lf = ElasticEnergy(theta_l);
    q2 = -2.0 * M * (E_lf - E_li);

    //std::cout << std::endl;
    //std::cout << theta_l * 180.0 / Pi << " " << q2 << " " << std::endl;
    mainf_(&randin, &randout, &ebeam, &q2, &vmin, &vcut, &phi, &ppx, &ppy, &ppz, &pp0, &epx, &epy, &epz, &ep0, &phx, &phy, &phz, &ph0, &pp, &pe, &pph, &q2p, &phip, &sigmaborn, &sigmatot, &sigmarad, &sigmabsv, &weightsig, &ptgen, &pvgen, &pichgen);
    //std::cout << epx << " " << epy << " " << epz << " " << ep0 << std::endl;
    //std::cout << ppx << " " << ppy << " " << ppz << " " << pp0 << std::endl;
    //std::cout << phx << " " << phy << " " << phz << " " << ph0 << std::endl;

    // convert ds/dQ2/dphi to ds/dOmega
    double jacob = 2 * Pow2(E_li) / Pow2(1 + E_li / M * (1 - Cos(theta_l))) * 1e-3;
    sigmaborn *= jacob;
    sigmatot *= jacob;
    sigmarad *= jacob;
    sigmabsv *= jacob;

    return sigmatot * Sin(theta_l);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RadCorEPXS(double qq)
{
    Clear();

    q2 = qq;

    mainf_(&randin, &randout, &ebeam, &q2, &vmin, &vcut, &phi, &ppx, &ppy, &ppz, &pp0, &epx, &epy, &epz, &ep0, &phx, &phy, &phz, &ph0, &pp, &pe, &pph, &q2p, &phip, &sigmaborn, &sigmatot, &sigmarad, &sigmabsv, &weightsig, &ptgen, &pvgen, &pichgen);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Clear()
{
    q2 = 0.0;
    phi = 0.0;
    ppx = 0.0;
    ppy = 0.0;
    ppz = 0.0;
    pp0 = 0.0;
    epx = 0.0;
    epy = 0.0;
    epz = 0.0;
    ep0 = 0.0;
    phx = 0.0;
    phy = 0.0;
    phz = 0.0;
    ph0 = 0.0;
    pp = 0.0;
    pe = 0.0;
    pph = 0.0;
    q2p = 0.0;
    phip = 0.0;
    sigmaborn = 0.0;
    sigmatot = 0.0;
    sigmarad = 0.0;
    sigmabsv = 0.0;
    weightsig = 0.0;
    ptgen = 0.0;
    pvgen = 0.0;
    pichgen = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
