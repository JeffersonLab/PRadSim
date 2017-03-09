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

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void usage(int, char **argv)
{
    printf("usage: %s [options] FILE_NAME\n", argv[0]);
    printf("  -h, --help                 Print usage\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    std::string filename;

    while (1) {
        static struct option long_options[] = {
            {"help",  no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'h':
            usage(argc, argv);
            exit(0);
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

    handler->SetHyCalSystem(hycal);
    handler->ReadFromEvio(filename);

    std::string tf = filename;
    size_t ppos = tf.rfind('.', tf.length());

    if (ppos != std::string::npos)
        tf = tf.substr(0, ppos);

    std::string outf = tf + "_rec.root";

    TFile *f = new TFile(outf.c_str(), "RECREATE");
    TTree *t = new TTree("T", "Reconstructed Sim Results");

    int N;
    double E[100], X[100], Y[100], Z[100]; // maximum number of clusters, 100 is enough
    // retrieve part of the cluster information
    t->Branch("HC.N", &N, "HC.N/I");
    t->Branch("HC.In.X", X, "HC.In.X[HC.N]/D");
    t->Branch("HC.In.Y", Y, "HC.In.Y[HC.N]/D");
    t->Branch("HC.In.Z", Z, "HC.In.Z[HC.N]/D");
    t->Branch("HC.In.P", E, "HC.In.P[HC.N]/D");

    int i = 1;

    for (auto &event : handler->GetEventData()) {
        if (i % 1000 == 0)
            std::cout << i << " events processed" << std::endl;

        hycal->Reconstruct(event);
        auto &hits = hycal->GetDetector()->GetHits();
        N = (int)hits.size();

        for (size_t i = 0; i < hits.size(); ++i) {
            X[i] = hits[i].x;
            Y[i] = hits[i].y;
            Z[i] = 5640.0 - 3000.0 + 88.9;
            E[i] = hits[i].E;
        }

        t->Fill();

        i++;
    }

    f->cd();
    t->Write();
    f->Close();

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
