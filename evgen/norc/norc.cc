//
// norc.cxx
// Developer : Chao Gu
// History:
//   Apr 2017, C. Gu, elastic event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"

#ifdef ELASTIC_EE
    #include "ee.hh"
#elif ELASTIC_ED
    #include "ed.hh"
#else
    #include "ep.hh"
#endif

#include <cstdlib>
#include <cstdio>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main()
{
    char mychar[64];

    std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    std::cin.getline(mychar, 64);
    Ei_1 = 0.001 * atof(mychar); // GeV

    std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_min = atof(mychar) * deg;

    std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_max = atof(mychar) * deg;

    std::cout << "Parameterization to be used: " << std::flush;
    std::cin.getline(mychar, 64);
    fselect = atoi(mychar);

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int N = atoi(mychar);

    TRandom2 *PseRan = new TRandom2();
    PseRan->SetSeed((int)(time(NULL)));

    phi_min = -180.0 * deg;
    phi_max = 180.0 * deg;

    omega = (phi_max - phi_min) * (Cos(theta_min) - Cos(theta_max));

    vi_1.SetPxPyPzE(0., 0., Sqrt(Pow2(Ei_1) - m2), Ei_1);
    vi_2.SetPxPyPzE(0., 0., 0., M);

    for (int i = 0; i < InterpolPoints; i++) {
        theta_1 = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
        theta[i] = theta_1;
        xs_sin[i] = ElasticXS_Sin(theta_1);
    }

    Integrator_XS_Sin.SetFunction(Func_XS_Sin);
    Integrator_XS_Sin.SetRelTolerance(IntTol);

    double xsint = (phi_max - phi_min) * mkb * Integrator_XS_Sin.Integral(theta_min, theta_max);

    Interpolator_XS_Sin.SetData(InterpolPoints, theta, xs_sin);

    TFoam *FoamX = new TFoam("FoamX");
    TFoamIntegrand *pFoamI = new ElasticIntegrand();
    FoamX->SetkDim(1);
    FoamX->SetnCells(10000); // Set number of cells
    FoamX->SetnSampl(500); // Set number of samples
    FoamX->SetOptRej(1); // Unweighted events in MC generation
    FoamX->SetRho(pFoamI); // Set distribution function
    FoamX->SetPseRan(PseRan); // Set random number generator
    //FoamX->SetChat(1); // Set "chat level" in the standard output
    FoamX->Initialize();

    FILE *fp = fopen(filename, "w");

    for (int i = 0; i < N; ++i) {
        if (i % 10000 == 0 && i != 0) std::cout << i << std::endl;

        FoamX->MakeEvent();

        phi_1 = phi_min + (phi_max - phi_min) * (PseRan->Rndm());

        if (phi_1 < 0.) phi_2 = phi_1 + Pi;
        else phi_2 = phi_1 - Pi;

        fprintf(fp, "%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", 1000. * Ef_1, theta_1, phi_1, 1000. * Ef_2, theta_2, phi_2, 0., 0., 0.);
        //printf("%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", 1000. * Ef_1, theta_1, phi_1, 1000. * Ef_2, theta_2, phi_2, 0., 0., 0.);
    }

    fclose(fp);

    fp = fopen(ifilename, "w");

    fprintf(fp, "beam energy:\n");
    fprintf(fp, "%lf MeV\n", Ei_1 * 1000);
    fprintf(fp, "polar angle range:\n");
    fprintf(fp, "%lf ~ %lf deg\n", theta_min / deg, theta_max / deg);
    fprintf(fp, "angle acceptance (the solid angle):\n");
    fprintf(fp, "%lf steradian\n", omega);
    fprintf(fp, "cross section (averaged over the solid angle):\n");
    fprintf(fp, "%lf microbarn / steradian\n", xsint / omega);
    fprintf(fp, "integrated luminosity:\n");
    fprintf(fp, "%lf inverse microbarn\n", N / xsint);

    printf("beam energy:\n");
    printf("%lf MeV\n", Ei_1 * 1000);
    printf("polar angle range:\n");
    printf("%lf ~ %lf deg\n", theta_min / deg, theta_max / deg);
    printf("angle acceptance (the solid angle):\n");
    printf("%lf steradian\n", omega);
    printf("cross section (averaged over the solid angle):\n");
    printf("%lf microbarn / steradian\n", xsint / omega);
    printf("integrated luminosity:\n");
    printf("%lf inverse microbarn\n", N / xsint);

    fclose(fp);
}
