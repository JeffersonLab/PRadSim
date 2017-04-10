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

#ifdef ELASTIC_ED
    #include "ed.h"
#else
    #include "ep.h"
#endif

#include <cstdlib>
#include <cstdio>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main()
{
    char mychar[64];
    /*
    std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    std::cin.getline(mychar, 64);
    Ei = 0.001 * atof(mychar); // MeV

    std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    thetamin = atof(mychar) * deg;

    std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    thetamax = atof(mychar) * deg;

    std::cout << "Parameterization to be used: " << std::flush;
    std::cin.getline(mychar, 64);
    select = atoi(mychar);

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int N = atoi(mychar);
    */
    E_li = 1.1;
    theta_min = 0.6 * deg;
    theta_max = 6.0 * deg;
    select = 1;
    int N = 100000;

    PseRan->SetSeed(0);

    phi_min = -180.0 * deg;
    phi_max = 180.0 * deg;

    omega = (phi_max - phi_min) * (Cos(theta_min) - Cos(theta_max));

    v_li.SetPxPyPzE(0., 0., Sqrt(Pow2(E_li) - m2), E_li);
    v_pi.SetPxPyPzE(0., 0., 0., M);

    FILE *fp = fopen(filename, "w");

    for (int i = 0; i < InterpolPoints; i++) {
        theta_l = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
        E_lf = ElasticEnergy(theta_l);

        theta[i] = theta_l;
        xs_sin[i] = ElasticXS_Sin(theta_l);
    }

    Integrator_XS_Sin.SetFunction(Func_XS_Sin);
    Integrator_XS_Sin.SetRelTolerance(IntTol);

    double xsint = (phi_max - phi_min) * mkb * Integrator_XS_Sin.Integral(theta_min, theta_max);

    Interpolator_XS_Sin.SetData(InterpolPoints, theta, xs_sin);

    TFoam *FoamX = new TFoam("FoamX");
    TFoamIntegrand *pFoamI = new ElasticIntegrand();
    FoamX->SetkDim(1);
    FoamX->SetnCells(1000); // Set number of cells
    FoamX->SetnSampl(200); // Set number os samples
    FoamX->SetOptRej(1); // Unweighted events in MC generation
    FoamX->SetRho(pFoamI); // Set distribution function
    FoamX->SetPseRan(PseRan); // Set random number generator
    //FoamX->SetChat(1); // Set "chat level" in the standard output
    FoamX->Initialize();

    for(int i = 0; i < N; ++i) {
        FoamX->MakeEvent();

        phi_l = phi_min + (phi_max - phi_min) * (PseRan->Rndm());

        if (phi_l < 0.) phi_p = phi_l + Pi;
        else phi_p = phi_l - Pi;

        fprintf(fp, "%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", 1000. * E_lf, theta_l, phi_l, 1000. * E_p, theta_p, phi_p, 0., 0., 0.);
        //printf("%8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf %8.2lf %6.4lf %7.4lf\n", 1000. * E_lf, theta_l, phi_l, 1000. * E_p, theta_p, phi_p, 0., 0., 0.);
    }

    std::cout << std::endl;
    std::cout << "cross section (averaged over the solid angle):" << std::endl;
    std::cout << xsint / omega << " microbarn / steradian" << std::endl;
    std::cout << "integrated luminosity:" << std::endl;
    std::cout << N / xsint << " inverse microbarn" << std::endl;

    fclose(fp);
}
