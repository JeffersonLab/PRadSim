 //
// newep.cc
// Developer : Chao Gu
// Based on Eur. Phys. J. A 51(2015)1
// History:
//   Apr 2017, C. Gu, ep event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"

#include "newep.hh"

#include <cstdlib>
#include <cstdio>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main()
{
    char mychar[64];

    std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    std::cin.getline(mychar, 64);
    Ei_e = atof(mychar) * 1e-3;

    std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_min = atof(mychar) * deg;

    std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_max = atof(mychar) * deg;

    std::cout << "Minimum elasticity cut (MeV^2, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempvmin = atof(mychar);

    if (tempvmin <= 0)
        v_min = 0.000025335;
    else
        v_min = tempvmin * 1e-6;

    std::cout << "Maximum elasticity cut (MeV^2, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempvcut = atof(mychar);

    if (tempvcut <= 0)
        v_cut = 0.05;
    else
        v_cut = tempvcut * 1e-6;

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int N = atoi(mychar);
/*
    Ei_e = 2.142;
    theta_min = 0.5 * deg;
    theta_max = 8.0 * deg;
    int N = 100000;
    v_min = 25.335e-6;
    v_cut = Pow2(2000.0) * 1e-6;
    //v_cut = 10000.0e-6;
*/
    PseRan = new TRandom2();
    PseRan->SetSeed((int)(time(NULL)));

    omega = 2.0 * pi * (Cos(theta_min) - Cos(theta_max));

    vi_e.SetPxPyPzE(0.0, 0.0, Sqrt(Pow2(Ei_e) - m2), Ei_e);
    vi_p.SetPxPyPzE(0.0, 0.0, 0.0, M);

    elrad_init(Ei_e, v_min);

    FILE *fp = fopen("xs.dat", "w");
/*
    for (int i = 0; i < 101; i++) {
        double q2 = Exp((-6 + i * 0.06) * Log(10.0));
        double ef = Ei_e - q2 / 2.0 / M;
        double tt = ACos(1.0 - (Ei_e / ef - 1.0) * M / Ei_e);

        double sinth = Sin(tt);

        double sigma_born = BornXS_Sin(tt) / sinth;

        double sigma_elastic = ElasticXS_Sin(tt) / sinth;
        double sigma_brems = BremsXS_Sin(tt) / sinth;
        double sigma_total = sigma_elastic + sigma_brems;

        fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
        printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
    }

    fclose(fp);
    exit(0);
*/
    for (int i = 0; i < InterpolPoints; i++) {
        theta_e = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
        theta[i] = theta_e;

        double sinth = Sin(theta_e);

        xs_elastic_sin[i] = ElasticXS_Sin(theta_e);
        xs_brems_sin[i] = BremsXS_Sin(theta_e);
        xs_sin[i] = xs_elastic_sin[i] + xs_brems_sin[i];

        double sigma_born = BornXS_Sin(theta_e) / sinth;

        double sigma_total = (xs_elastic_sin[i] + xs_brems_sin[i]) / sinth;
        double sigma_elastic = xs_elastic_sin[i] / sinth;
        double sigma_brems = xs_brems_sin[i] / sinth;

        double q2 = - (vf_e - vi_e) * (vf_e - vi_e);

        fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_e * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
        printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_e * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
    }

    fclose(fp);

    Integrator_EPXS_Sin.SetFunction(Func_EPXS_Sin);
    Integrator_EPXS_Sin.SetRelTolerance(IntTol);

    Integrator_ElasticXS_Sin.SetFunction(Func_ElasticXS_Sin);
    Integrator_ElasticXS_Sin.SetRelTolerance(IntTol);

    //double xsint = 2 * pi * Integrator_EPXS_Sin.Integral(theta_min, theta_max);
    //double elxsint = 2 * pi * Integrator_ElasticXS_Sin.Integral(theta_min, theta_max);

    //int n_elastic = int(N * (elxsint / xsint));

    Interpolator_EPXS_Sin.SetData(InterpolPoints, theta, xs_sin);
    Interpolator_ElasticXS_Sin.SetData(InterpolPoints, theta, xs_elastic_sin);
    Interpolator_BremsXS_Sin.SetData(InterpolPoints, theta, xs_brems_sin);

    double xsint = 2 * pi * Interpolator_EPXS_Sin.Integ(theta_min, theta_max);
    double elxsint = 2 * pi * Interpolator_ElasticXS_Sin.Integ(theta_min, theta_max);

    int n_elastic = int(N * (elxsint / xsint));

    TFoam *FoamElastic = new TFoam("FoamElastic");
    TFoamIntegrand *pFoamElastic = new ElasticIntegrand();
    FoamElastic->SetkDim(1);
    FoamElastic->SetnCells(10000); // Set number of cells
    FoamElastic->SetnSampl(2000); // Set number of samples
    FoamElastic->SetOptRej(1); // Unweighted events in MC generation
    FoamElastic->SetRho(pFoamElastic); // Set distribution function
    FoamElastic->SetPseRan(PseRan); // Set random number generator
    //FoamElastic->SetChat(1); // Set "chat level" in the standard output
    FoamElastic->Initialize();

    TFoam *FoamBrems = new TFoam("FoamBrems");
    TFoamIntegrand *pFoamBrems = new BremsIntegrand();
    FoamBrems->SetkDim(1);
    FoamBrems->SetnCells(10000); // Set number of cells
    FoamBrems->SetnSampl(2000); // Set number of samples
    FoamBrems->SetOptRej(1); // Unweighted events in MC generation
    FoamBrems->SetRho(pFoamBrems); // Set distribution function
    FoamBrems->SetPseRan(PseRan); // Set random number generator
    //FoamBrems->SetChat(1); // Set "chat level" in the standard output
    FoamBrems->Initialize();

    fp = fopen(filename, "w");

    int count_elastic = 0, count_brems = 0;

    for (int i = 0; i < N; ++i) {
        if (i % 10000 == 0 && i != 0) std::cout << i << std::endl;

        if ((n_elastic - count_elastic) > 0 && PseRan->Rndm() < 1.0 * (n_elastic - count_elastic) / (N - i)) {
            FoamElastic->MakeEvent();

            Ef_e = ElasticEnergy(theta_e);

            vf_e.SetPxPyPzE(Sqrt(Pow2(Ef_e) - m2) * Sin(theta_e), 0.0, Sqrt(Pow2(Ef_e) - m2) * Cos(theta_e), Ef_e);
            vf_p = vi_e + vi_p - vf_e;

            double phi = 2.0 * pi * PseRan->Rndm();

            vf_e.RotateZ(phi);
            vf_p.RotateZ(phi);

            Ef_e = vf_e.E();
            theta_e = vf_e.Theta();
            phi_e = vf_e.Phi();

            Ef_p = vf_p.E();
            theta_p = vf_p.Theta();
            phi_p = vf_p.Phi();

            E_g = 0.0;
            theta_g = 0.0;
            phi_g = 0.0;

            count_elastic++;
        } else {
            FoamBrems->MakeEvent();

            RecBremsKins(theta_e);

            double phi = 2.0 * pi * PseRan->Rndm();

            vf_e.RotateZ(phi);
            vf_p.RotateZ(phi);
            v_g.RotateZ(phi);

            Ef_e = vf_e.E();
            theta_e = vf_e.Theta();
            phi_e = vf_e.Phi();

            Ef_p = vf_p.E();
            theta_p = vf_p.Theta();
            phi_p = vf_p.Phi();

            E_g = v_g.E();
            theta_g = v_g.Theta();
            phi_g = v_g.Phi();

            count_brems++;
        }

        fprintf(fp, "%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_e * 1000, theta_e, phi_e, Ef_p * 1000, theta_p, phi_p, E_g * 1000, theta_g, phi_g);
        //printf("%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_e * 1000, theta_e, phi_e, Ef_p * 1000, theta_p, phi_p, E_g * 1000, theta_g, phi_g);
    }

    fclose(fp);

    fp = fopen(ifilename, "w");

    fprintf(fp, "beam energy:\n");
    fprintf(fp, "%lf MeV\n", Ei_e * 1000);
    fprintf(fp, "polar angle range:\n");
    fprintf(fp, "%lf ~ %lf deg\n", theta_min / deg, theta_max / deg);
    fprintf(fp, "angle acceptance (the solid angle):\n");
    fprintf(fp, "%lf steradian\n", omega);
    fprintf(fp, "cross section (averaged over the solid angle):\n");
    fprintf(fp, "%lf microbarn / steradian\n", xsint / omega);
    fprintf(fp, "integrated luminosity:\n");
    fprintf(fp, "%lf inverse microbarn\n", N / xsint);

    printf("beam energy:\n");
    printf("%lf MeV\n", Ei_e * 1000);
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
