//
// moller.cc
// Developer : Chao Gu
// Based on Chao Peng's RC moller cross section code in PRadAnalyzer
// History:
//   Apr 2017, C. Gu, moller event generator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"

#include "moller.hh"

#include <cstdlib>
#include <cstdio>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main()
{
    char mychar[64];

    std::cout << "Full energy of the incident lepton (MeV): " << std::flush;
    std::cin.getline(mychar, 64);
    Ei_1 = atof(mychar);

    std::cout << "Minimum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_min = atof(mychar) * deg;

    std::cout << "Maximum polar angle of the electron (degree): " << std::flush;
    std::cin.getline(mychar, 64);
    theta_max = atof(mychar) * deg;

    std::cout << "Minimum elasticity cut (MeV, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempvmin = atof(mychar);

    if (tempvmin <= 0)
        v_min = 4.69;
    else
        v_min = tempvmin;

    std::cout << "Maximum elasticity cut (MeV, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempvcut = atof(mychar);

    if (tempvcut <= 0)
        v_cut = 100;
    else
        v_cut = tempvcut;

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int N = atoi(mychar);

    //Ei_1 = 1100;
    //theta_min = 0.7 * deg;
    //theta_max = 3.8 * deg;
    //int N = 100;
    //v_min = 4.69;
    //v_cut = 100;

    PseRan = new TRandom2();
    PseRan->SetSeed((int)(time(NULL)));

    omega = 2.0 * pi * (Cos(theta_min) - Cos(theta_max));

    vi_1.SetPxPyPzE(0.0, 0.0, Sqrt(Pow2(Ei_1) - m2), Ei_1);
    vi_2.SetPxPyPzE(0.0, 0.0, 0.0, m);

    merad_init(Ei_1);

    FILE *fp = fopen("xs.dat", "w");
    /*
        for (int i = 0; i < 101; i++) {
            double q2 = Exp((i * 0.03) * Log(10.0));
            double ef = Ei_1 - q2 / (2 * m);
            double tt = ACos(Sqrt((Ei_1 + m) * (ef - m) / (Ei_1 - m) / (ef + m)));

            double sinth = Sin(tt);

            double sigma_born = BornXS_Sin(tt) / sinth;

            double sigma_elastic = ElasticXS_Sin(tt) / sinth;
            double sigma_brems = BremsXS_Sin(tt) / sinth;
            double sigma_total = sigma_elastic + sigma_brems;

            fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / pi, q2 / 1e6, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
            printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", tt * 180.0 / pi, q2 / 1e6, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
        }

        fclose(fp);
        exit(0);
    */

    for (int i = 0; i < InterpolPoints; i++) {
        theta_1 = theta_min + i * (theta_max - theta_min) / (InterpolPoints - 1);
        theta[i] = theta_1;

        double sinth = Sin(theta_1);

        xs_elastic_sin[i] = ElasticXS_Sin(theta_1);
        xs_brems_sin[i] = BremsXS_Sin(theta_1);

        double sigma_born = BornXS_Sin(theta_1) / sinth;

        double sigma_total = (xs_elastic_sin[i] + xs_brems_sin[i]) / sinth;
        double sigma_elastic = xs_elastic_sin[i] / sinth;
        double sigma_brems = xs_brems_sin[i] / sinth;

        double q2 = - (vf_1 - vi_1) * (vf_1 - vi_1) / 1e6;

        //printf("%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_1 * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
        fprintf(fp, "%8.6lf %10.4le %10.4le %10.4le %10.4le %10.4le %8.6lf %8.6lf %8.6lf\n", theta_1 * 180.0 / pi, q2, sigma_born, sigma_total, sigma_elastic, sigma_brems, sigma_total / sigma_born - 1, sigma_elastic / sigma_born, sigma_brems / sigma_born);
    }

    fclose(fp);

    Integrator_BornXS_Sin.SetFunction(Func_BornXS_Sin);
    Integrator_BornXS_Sin.SetRelTolerance(IntTol);

    Integrator_MollerXS_Sin.SetFunction(Func_MollerXS_Sin);
    Integrator_MollerXS_Sin.SetRelTolerance(IntTol);

    Integrator_ElasticXS_Sin.SetFunction(Func_ElasticXS_Sin);
    Integrator_ElasticXS_Sin.SetRelTolerance(IntTol);

    double xsint_born = 2 * pi * Integrator_BornXS_Sin.Integral(theta_min, theta_max);
    double xsint = 2 * pi * Integrator_MollerXS_Sin.Integral(theta_min, theta_max);
    double elxsint = 2 * pi * Integrator_ElasticXS_Sin.Integral(theta_min, theta_max);

    int n_elastic = int(N * (elxsint / xsint));

    std::cerr << xsint_born << " " << xsint << std::endl;

    Interpolator_ElasticXS_Sin.SetData(InterpolPoints, theta, xs_elastic_sin);
    Interpolator_BremsXS_Sin.SetData(InterpolPoints, theta, xs_brems_sin);

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

            Ef_1 = ElasticEnergy(theta_1);

            vf_1.SetPxPyPzE(Sqrt(Pow2(Ef_1) - m2) * Sin(theta_1), 0.0, Sqrt(Pow2(Ef_1) - m2) * Cos(theta_1), Ef_1);
            vf_2 = vi_1 + vi_2 - vf_1;

            double phi = 2.0 * pi * PseRan->Rndm();

            vf_1.RotateZ(phi);
            vf_2.RotateZ(phi);

            Ef_1 = vf_1.E();
            theta_1 = vf_1.Theta();
            phi_1 = vf_1.Phi();

            Ef_2 = vf_2.E();
            theta_2 = vf_2.Theta();
            phi_2 = vf_2.Phi();

            E_g = 0.0;
            theta_g = 0.0;
            phi_g = 0.0;

            count_elastic++;
        } else {
            FoamBrems->MakeEvent();

            RecBremsKins(theta_1);

            double phi = 2.0 * pi * PseRan->Rndm();

            vf_1.RotateZ(phi);
            vf_2.RotateZ(phi);
            v_g.RotateZ(phi);

            Ef_1 = vf_1.E();
            theta_1 = vf_1.Theta();
            phi_1 = vf_1.Phi();

            Ef_2 = vf_2.E();
            theta_2 = vf_2.Theta();
            phi_2 = vf_2.Phi();

            E_g = v_g.E();
            theta_g = v_g.Theta();
            phi_g = v_g.Phi();

            count_brems++;
        }

        //printf("%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_1, theta_1, phi_1, Ef_2, theta_2, phi_2, E_g, theta_g, phi_g);
        fprintf(fp, "%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_1, theta_1, phi_1, Ef_2, theta_2, phi_2, E_g, theta_g, phi_g);
    }

    fclose(fp);

    fp = fopen(ifilename, "w");

    fprintf(fp, "beam energy:\n");
    fprintf(fp, "%lf MeV\n", Ei_1);
    fprintf(fp, "polar angle range:\n");
    fprintf(fp, "%lf ~ %lf deg\n", theta_min / deg, theta_max / deg);
    fprintf(fp, "angle acceptance (the solid angle):\n");
    fprintf(fp, "%lf steradian\n", omega);
    fprintf(fp, "cross section (averaged over the solid angle):\n");
    fprintf(fp, "%lf microbarn / steradian\n", xsint / omega);
    fprintf(fp, "integrated luminosity:\n");
    fprintf(fp, "%lf inverse microbarn\n", N / xsint);

    printf("beam energy:\n");
    printf("%lf MeV\n", Ei_1);
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
