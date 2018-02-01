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

    std::cout << "Minimum energy of bremsstrahlung photons (fraction of CM energy, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempfracmin = atof(mychar);

    if (tempfracmin <= 0 || tempfracmin >= 1.0)
        E_g_frac_min = 1e-3;
    else
        E_g_frac_min = tempfracmin;

    std::cout << "Maximum energy of bremsstrahlung photons (fraction of CM energy, -1 to use default value): " << std::flush;
    std::cin.getline(mychar, 64);
    double tempfracmax = atof(mychar);

    if (tempfracmax <= 0 || tempfracmax >= 1.0)
        E_g_frac_max = 0.95;
    else
        E_g_frac_max = tempfracmax;

    std::cout << "Number of events to generate: " << std::flush;
    std::cin.getline(mychar, 64);
    int N = atoi(mychar);

    PseRan = new TRandom2();
    PseRan->SetSeed((int)(time(NULL)));

    omega = 2.0 * pi * (Cos(theta_min) - Cos(theta_max));

    vi_1.SetPxPyPzE(0.0, 0.0, Sqrt(Pow2(Ei_1) - m2), Ei_1);
    vi_2.SetPxPyPzE(0.0, 0.0, 0.0, m);

    cm = vi_1 + vi_2;

    vi_1_cm = vi_1;
    vi_1_cm.Boost(-cm.BoostVector());

    vi_2_cm = vi_2;
    vi_2_cm.Boost(-cm.BoostVector());

    E_cm = vi_1_cm.E() + vi_2_cm.E();
    E_cmp = E_cm / 2.0;
    P_cmp = Sqrt(Pow2(E_cmp) - m2);

    TLorentzVector temp, temp_cm;
    double etemp = ElasticEnergy(theta_min);
    temp.SetPxPyPzE(Sqrt(Pow2(etemp) - m2) * Sin(theta_min), 0.0, Sqrt(Pow2(etemp) - m2) * Cos(theta_min), etemp);
    temp_cm = temp;
    temp_cm.Boost(-cm.BoostVector());
    theta_min_cm = temp_cm.Theta() - 15.0 * deg;

    etemp = ElasticEnergy(theta_max);
    temp.SetPxPyPzE(Sqrt(Pow2(etemp) - m2) * Sin(theta_max), 0.0, Sqrt(Pow2(etemp) - m2) * Cos(theta_max), etemp);
    temp_cm = temp;
    temp_cm.Boost(-cm.BoostVector());
    theta_max_cm = temp_cm.Theta() + 15.0 * deg;

    E_g_min = E_g_frac_min * E_cm;
    E_g_max = E_cmp - m2 / E_cmp - m;

    if (E_g_max > (E_g_frac_max * E_cm)) E_g_max = E_g_frac_max * E_cm;

    for (int i = 0; i < InterpolPoints; i++) {
        theta_1_cm = theta_min_cm + i * (theta_max_cm - theta_min_cm) / (InterpolPoints - 1);
        theta_cm[i] = theta_1_cm;

        double sinth = Sin(theta_1_cm);

        xs_elastic_sin[i] = ElasticXS_Sin(theta_1_cm);

        double sigma_born = BornXS_Sin(theta_1_cm) / sinth;
        double sigma_elastic = xs_elastic_sin[i] / sinth;

        vf_1.SetPxPyPzE(P_cmp * Sin(theta_1_cm), 0.0, P_cmp * Cos(theta_1_cm), E_cmp);
        double q2 = -(vi_1 - vf_1) * (vi_1 - vf_1);

        printf("%8.6lf %10.4le %10.4le %10.4le %8.6lf\n", theta_cm[i] * 180.0 / pi, q2, sigma_born, sigma_elastic, sigma_elastic / sigma_born);
    }

    Interpolator_ElasticXS_Sin.SetData(InterpolPoints, theta_cm, xs_elastic_sin);

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

    TFoam *FoamBremsM = new TFoam("FoamBremsM");
    TFoamIntegrand *pFoamBremsM = new BremsIntegrandM();
    FoamBremsM->SetkDim(4);
    FoamBremsM->SetnCells(30000); // Set number of cells
    FoamBremsM->SetnSampl(1500); // Set number of samples
    FoamBremsM->SetOptRej(1); // Unweighted events in MC generation
    FoamBremsM->SetRho(pFoamBremsM); // Set distribution function
    FoamBremsM->SetPseRan(PseRan); // Set random number generator
    //FoamBremsM->SetChat(1); // Set "chat level" in the standard output
    FoamBremsM->Initialize();

    TFoam *FoamBremsP = new TFoam("FoamBremsP");
    TFoamIntegrand *pFoamBremsP = new BremsIntegrandP();
    FoamBremsP->SetkDim(4);
    FoamBremsP->SetnCells(30000); // Set number of cells
    FoamBremsP->SetnSampl(1500); // Set number of samples
    FoamBremsP->SetOptRej(1); // Unweighted events in MC generation
    FoamBremsP->SetRho(pFoamBremsP); // Set distribution function
    FoamBremsP->SetPseRan(PseRan); // Set random number generator
    //FoamBremsP->SetChat(1); // Set "chat level" in the standard output
    FoamBremsP->Initialize();

    double xsint_elastic = 2 * pi * Interpolator_ElasticXS_Sin.Integ(theta_min_cm, theta_max_cm);
    double xsint_brems_m, xsint_brems_m_error;
    double xsint_brems_p, xsint_brems_p_error;

    xsint_brems_p = FoamBremsP->GetPrimary();

    if (xsint_brems_p == 0) std::cout << "FoamBremP: no need to use." << std::endl;

    for (int i = 0; i < 10000000; i++) {
        if (i % 100000 == 0 && i != 0) std::cout << " Initializing (-):" << i << "\r" << std::flush;

        FoamBremsM->MakeEvent();
    }

    std::cout << " Initializing (-): 10000000" << std::endl;
    FoamBremsM->GetIntegMC(xsint_brems_m, xsint_brems_m_error);

    if (xsint_brems_p > 0) {
        for (int i = 0; i < 500000; i++) {
            if (i % 100000 == 0 && i != 0) std::cout << " Initializing (+):" << i << "\r" << std::flush;

            FoamBremsP->MakeEvent();
        }

        std::cout << " Initializing (+): 500000" << std::endl;
        FoamBremsP->GetIntegMC(xsint_brems_p, xsint_brems_p_error);
    }

    double xsint_brems;

    xsint_brems = xsint_brems_m + xsint_brems_p;

    double xsint = xsint_elastic + xsint_brems;

    int n_elastic = int(N * (xsint_elastic / xsint));
    int n_brems_p = int(N * (xsint_brems_p / xsint));
    int n_brems_m = N - n_elastic - n_brems_p;

    //std::cout << xsint_elastic << " " << 0.0 << " " << xsint_brems << " " << xsint_brems_error << " " << xsint << " " << xsint_error << std::endl;

    FILE *fp = fopen(filename, "w");

    int count_elastic = 0, count_brems_m = 0, count_brems_p = 0;

    for (int i = 0; i < N; ++i) {
        if (i % 10000 == 0 && i != 0) std::cout << i << std::endl;

        if ((n_elastic - count_elastic) > 0 && PseRan->Rndm() < 1.0 * (n_elastic - count_elastic) / (N - i)) {
            FoamElastic->MakeEvent();

            vf_1_cm.SetPxPyPzE(P_cmp * Sin(theta_1_cm), 0.0, P_cmp * Cos(theta_1_cm), E_cmp);
            vf_2_cm.SetPxPyPzE(-P_cmp * Sin(theta_1_cm), 0.0, -P_cmp * Cos(theta_1_cm), E_cmp);
            vf_1 = vf_1_cm;
            vf_2 = vf_2_cm;
            vf_1.Boost(cm.BoostVector());
            vf_2.Boost(cm.BoostVector());

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
            if ((n_brems_p - count_brems_p) > 0 && (count_brems_p / n_brems_p < count_brems_m / n_brems_m)) {
                FoamBremsP->MakeEvent();

                vf_1 = vf_1_cm;
                vf_2 = vf_2_cm;
                v_g = v_g_cm;
                vf_1.Boost(cm.BoostVector());
                vf_2.Boost(cm.BoostVector());
                v_g.Boost(cm.BoostVector());

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

                count_brems_p++;
            } else {
                FoamBremsM->MakeEvent();

                vf_1 = vf_1_cm;
                vf_2 = vf_2_cm;
                v_g = v_g_cm;
                vf_1.Boost(cm.BoostVector());
                vf_2.Boost(cm.BoostVector());
                v_g.Boost(cm.BoostVector());

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

                count_brems_m++;
            }
        }

        fprintf(fp, "%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_1, theta_1, phi_1, Ef_2, theta_2, phi_2, E_g, theta_g, phi_g);
        //printf("%11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf %11.5lf %10.8lf %11.8lf\n", Ef_1, theta_1, phi_1, Ef_2, theta_2, phi_2, E_g, theta_g, phi_g);
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
