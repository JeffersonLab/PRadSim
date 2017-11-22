//
// PRadRec.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ConfigParser.h"
#include "PRadDataHandler.h"
#include "PRadHyCalSystem.h"
#include "PRadDetMatch.h"
#include "PRadCoordSystem.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TRandom2.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

#define T_BLOCKS 2156

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static TRandom2 *RandGen = new TRandom2();

double ECali[T_BLOCKS];
double nonlinConst[T_BLOCKS];

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LoadConst(double beam_energy);
Double_t EnergyCorrect(Double_t energy, Short_t cid);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void usage(int, char **argv)
{
    printf("usage: %s [options] FILE_NAME\n", argv[0]);
    printf("  -e, --energy=1100          Set beam energy (MeV)\n");
    printf("  -h, --help                 Print usage\n");
    printf("  -g, --gem_match=1          Do GEM matching\n");
    printf("  -t, --trg_eff=1            Do HyCal trigger efficiency\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    std::string filename;
    bool gem_match = false;
    bool trg_eff = false;
    double ei = 1100.0;

    while (1) {
        static struct option long_options[] = {
            {"help",  no_argument, 0, 'h'},
            {"energy",  required_argument, 0, 'e'},
            {"gem_match",  no_argument, 0, 'g'},
            {"trg_eff", no_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "e:ght", long_options, &option_index);

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

        case 't':
            trg_eff = true;
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

    // initiate random generator seed
    RandGen->SetSeed((UInt_t)time(NULL));

    // initiate calibration constants
    //LoadConst(ei);

    // simulation data is more like raw evio data with HyCal information only,
    // so we only need hycal system to connected to the handler
    PRadDataHandler *handler = new PRadDataHandler();
    PRadHyCalSystem *hycal = new PRadHyCalSystem("config/hycal.conf");
    PRadDetMatch *det_match = new PRadDetMatch("config/det_match.conf");
    PRadCoordSystem *coord_sys = new PRadCoordSystem("database/coordinates.dat", 2000);

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
        if (i % 1000 == 0 && i != 0) std::cout << i << " events processed" << std::endl;

        hycal->Reconstruct(event);
        auto &hits = hycal->GetDetector()->GetHits();
        coord_sys->Transform(PRadDetector::HyCal, hits.begin(), hits.end());

        tgem->GetEntry(i);

        if (N_GEM > 100) N_GEM = 100;

        if (gem_match) {
            std::vector<GEMHit> gem1_hits, gem2_hits;

            for (int j = 0; j < N_GEM; j++) {
                GEMHit h;
                h.x = - X_GEM[j]; // Orientation mismatch
                h.y = Y_GEM[j];
                h.z = 0;

                if (DID_GEM[j]) gem2_hits.push_back(h);
                else gem1_hits.push_back(h);
            }

            coord_sys->Transform(PRadDetector::PRadGEM1, gem1_hits.begin(), gem1_hits.end());
            coord_sys->Transform(PRadDetector::PRadGEM2, gem2_hits.begin(), gem2_hits.end());
            auto matched = det_match->Match(hits, gem1_hits, gem2_hits);

            N_HC = (int)matched.size();
            N_GEM = (int)matched.size();

            for (int j = 0; j < N_HC; ++j) {
                E[j] = EnergyCorrect(matched[j].hycal.E, matched[j].hycal.cid);

                CID[j] = matched[j].hycal.cid;

                if (trg_eff && (RandGen->Uniform() > hycal->GetModule(CID[j])->GetTriggerEfficiency(E[j]))) continue;

                X_HC[j] = matched[j].hycal.x;
                Y_HC[j] = matched[j].hycal.y;
                Z_HC[j] = matched[j].hycal.z;

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
                E[j] = hits[j].E; //EnergyCorrect(hits[j].E, hits[j].cid);
                CID[j] = hits[j].cid;

                if (trg_eff && (RandGen->Uniform() > hycal->GetModule(CID[j])->GetTriggerEfficiency(E[j]))) continue;

                X_HC[j] = hits[j].x;
                Y_HC[j] = hits[j].y;
                Z_HC[j] = hits[j].z;

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LoadConst(double beam_energy)
{
    ConfigParser parser;

    std::string path;

    if (beam_energy < 2000.) path = "./database/calibration/1GeV_mc_cali_const.dat";
    else path = "./database/calibration/2GeV_mc_cali_const.dat";

    if (!parser.OpenFile(path)) {
        std::cout << "cannot find mc calibration file" << std::endl;
        exit(0);
    }

    int count = 0;

    while (parser.ParseLine()) {
        double input[4];

        for (int i = 0; i < 4; i++) input[i] = parser.TakeFirst().Double();

        ECali[count] = input[2];
        nonlinConst[count] = input[3];
        count++;
    }

    parser.CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double EnergyCorrect(Double_t energy, Short_t cid)
{
    if (cid <= 0) return energy;

    float ecorr = 1. + nonlinConst[cid - 1] * (energy - ECali[cid - 1]) / 1000.;

    if (fabs(ecorr - 1.) < 0.6) energy /= ecorr;

    return energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
