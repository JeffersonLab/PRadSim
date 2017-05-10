//
// PRadRec.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PRadDataHandler.h"
#include "PRadHyCalSystem.h"
#include "PRadDetMatch.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void usage(int, char **argv)
{
    printf("usage: %s [options] FILE_NAME\n", argv[0]);
    printf("  -e, --energy=1100          Set beam energy (MeV)\n");
    printf("  -h, --help                 Print usage\n");
    printf("  -g, --gem_match=1          Do GEM matching\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double GetNonlinCorr(Double_t reconE)
{
    // reconE in MeV
    return exp(-1.0 * reconE * 1.53438e-04) + 1.11330e-04 * reconE + 7.17932e-02;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double ScaleEnergy(double e)
{
    if (e < 1600) {
        double p3 = 0.539865;
        double p2 = -0.271654;
        double p1 =  6.43518e-05;
        double p0 =  0.546174;

        return (p3 * TMath::Exp(p2 * TMath::Sqrt(e / 1100.0)) + p1 * e + p0) * e;
    } else {
        double p3 = 0.557136;
        double p2 = -0.351049;
        double p1 = 3.5904e-05;
        double p0 = 0.549588;

        return (p3 * TMath::Exp(p2 * TMath::Sqrt(e / 2141.0)) + p1 * e + p0) * e;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    std::string filename;
    bool gem_match = false;
    double ei = 1100.0;

    while (1) {
        static struct option long_options[] = {
            {"help",  no_argument, 0, 'h'},
            {"energy",  required_argument, 0, 'e'},
            {"gem_match",  no_argument, 0, 'g'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "e:gh", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'e':
            ei = atof(optarg);

            if (ei < 1.0) {
                usage(argc, argv);
                exit(1);
            }

            break;

        case 'h':
            usage(argc, argv);
            exit(0);
            break;

        case 'g':
            gem_match = true;
            break;

        case '?':
            // getopt_long already printed an error message
            break;

        default:
            usage(argc, argv);
            break;
        }
    }

    if (optind + 1 == argc)
        filename = argv[optind++];
    else {
        usage(argc, argv);
        exit(1);
    }

    // simulation data is more like raw evio data with HyCal information only,
    // so we only need hycal system to connected to the handler
    PRadDataHandler *handler = new PRadDataHandler();
    PRadHyCalSystem *hycal = new PRadHyCalSystem("config/hycal.conf");
    PRadDetMatch *det_match = new PRadDetMatch("config/det_match.conf");

    handler->SetHyCalSystem(hycal);
    handler->ReadFromEvio(filename);

    std::string tf = filename;
    size_t ppos = tf.rfind('.', tf.length());

    if (ppos != std::string::npos)
        tf = tf.substr(0, ppos);

    std::string outf = tf + "_rec.root";

    std::string gemf = tf + ".root";
    TFile *fgem = new TFile(gemf.c_str());
    TTree *tgem = (TTree *)fgem->Get("T");

    TFile *f = new TFile(outf.c_str(), "RECREATE");
    TTree *t = new TTree("T", "Reconstructed Sim Results");

    int N_HC, N_GEM, CID[100], DID_GEM[100];
    double E[100], X_HC[100], Y_HC[100], Z_HC[100], X_GEM[100], Y_GEM[100], Z_GEM[100]; // maximum number of clusters, 100 is enough
    // retrieve part of the cluster information
    t->Branch("HC.N", &N_HC, "HC.N/I");
    t->Branch("HC.X", X_HC, "HC.X[HC.N]/D");
    t->Branch("HC.Y", Y_HC, "HC.Y[HC.N]/D");
    t->Branch("HC.Z", Z_HC, "HC.Z[HC.N]/D");
    t->Branch("HC.P", E, "HC.P[HC.N]/D");
    t->Branch("HC.CID", CID, "HC.CID[HC.N]/I");
    t->Branch("GEM.N", &N_GEM, "GEM.N/I");
    t->Branch("GEM.X", X_GEM, "GEM.X[GEM.N]/D");
    t->Branch("GEM.Y", Y_GEM, "GEM.Y[GEM.N]/D");
    t->Branch("GEM.Z", Z_GEM, "GEM.Z[GEM.N]/D");
    t->Branch("GEM.DID", DID_GEM, "GEM.DID[GEM.N]/I");

    tgem->SetBranchAddress("GEM.N", &N_GEM);
    tgem->SetBranchAddress("GEM.X", X_GEM);
    tgem->SetBranchAddress("GEM.Y", Y_GEM);
    tgem->SetBranchAddress("GEM.Z", Z_GEM);
    tgem->SetBranchAddress("GEM.DID", DID_GEM);

    int i = 0;

    for (auto &event : handler->GetEventData()) {
        if (i % 1000 == 0 && i != 0)
            std::cout << i << " events processed" << std::endl;

        hycal->Reconstruct(event);
        auto &hits = hycal->GetDetector()->GetHits();

        tgem->GetEntry(i);

        if (N_GEM > 100) N_GEM = 100;

        if (gem_match) {
            std::vector<GEMHit> gem1_hits, gem2_hits;

            for (int j = 0; j < N_GEM; j++) {
                GEMHit h;
                h.x = - X_GEM[j]; // Orientation mismatch
                h.y = Y_GEM[j];
                h.z = Z_GEM[j] + 3000.0 - 89.0;

                if (DID_GEM[j]) gem2_hits.push_back(h);
                else gem1_hits.push_back(h);
            }

            for (int j = 0; j < (int)hits.size(); ++j) hits[j].z += 5640.0;

            auto matched = det_match->Match(hits, gem1_hits, gem2_hits);

            N_HC = (int)matched.size();
            N_GEM = (int)matched.size();

            for (int j = 0; j < N_HC; ++j) {

                X_HC[j] = matched[j].hycal.x;
                Y_HC[j] = matched[j].hycal.y;
                Z_HC[j] = matched[j].hycal.z;
                //E[j] = matched[j].hycal.E * GetNonlinCorr(matched[j].hycal.E);

                E[j] = ScaleEnergy(matched[j].hycal.E);

                CID[j] = matched[j].hycal.cid;

                if (matched[j].gem1.empty() && matched[j].gem2.empty()) {
                    X_GEM[j] = -10000;
                    Y_GEM[j] = -10000;
                    Z_GEM[j] = -10000;
                } else if (matched[j].gem1.empty()) {
                    X_GEM[j] = matched[j].gem2[0].x;
                    Y_GEM[j] = matched[j].gem2[0].y;
                    Z_GEM[j] = matched[j].gem2[0].z;
                } else if (matched[j].gem2.empty()) {
                    X_GEM[j] = matched[j].gem1[0].x;
                    Y_GEM[j] = matched[j].gem1[0].y;
                    Z_GEM[j] = matched[j].gem1[0].z;
                } else {
                    X_GEM[j] = 0.5 * (matched[j].gem1[0].x + matched[j].gem2[0].x);
                    Y_GEM[j] = 0.5 * (matched[j].gem1[0].y + matched[j].gem2[0].y);
                    Z_GEM[j] = 0.5 * (matched[j].gem1[0].z + matched[j].gem2[0].z);
                }
            }
        } else {
            N_HC = (int)hits.size();

            for (int j = 0; j < (int)hits.size(); ++j) {
                X_HC[j] = hits[j].x;
                Y_HC[j] = hits[j].y;
                Z_HC[j] = 5640.0 - 3000.0 + 89.0 + hits[j].z;
                //E[j] = hits[j].E * GetNonlinCorr(hits[j].E);
                E[j] = ScaleEnergy(hits[j].E);
                CID[j] = hits[j].cid;
            }
        }

        t->Fill();

        i++;
    }

    f->cd();
    t->Write();
    f->Close();

    fgem->Close();

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
