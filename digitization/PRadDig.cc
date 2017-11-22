//
// PRadDig.cc
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PRadDigitization.hh"
#include "HyCalDigitization.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TChain.h"

#include <string>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void usage(int, char **argv)
{
    printf("usage: %s [options] FILE_NAME\n", argv[0]);
    printf("  -l, --list=run.list        Set run list\n");
    printf("  -r, --run=1288             Set run number\n");
    printf("  -h, --help                 Print usage\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    int run = 1288;
    int    sMode = 0;
    std::string list_name;
    bool use_file_list = false;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"run",  required_argument, 0, 'r'},
            {"smear", required_argument, 0, 's'},
            {"list", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "e:hl:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'r':
            run = atoi(optarg);
            break;

        case 'h':
            usage(argc, argv);
            exit(0);
            break;

        case 'l':
            list_name = optarg;
            use_file_list = true;
            break;

        case 's':
            sMode = atoi(optarg);
            break;

        case '?':
            // getopt_long already printed an error message
            break;

        default:
            usage(argc, argv);
            break;
        }
    }

    std::string file_list[50];
    int n_file_list = 0;

    if (use_file_list) {
        std::ifstream ifs(list_name);

        if (ifs.is_open()) {
            std::string line;

            while (std::getline(ifs, line)) {
                if (line.at(0) == '#')
                    continue;

                file_list[n_file_list++] = line;
            }
        }

        ifs.close();
    } else if (optind + 1 == argc)
        file_list[n_file_list++] = argv[optind++];
    else {
        usage(argc, argv);
        exit(1);
    }

    TChain *t = new TChain("T");

    for (int i = 0; i < n_file_list; i++) {
        std::cout << "Add file " << file_list[i] << std::endl;
        t->Add(file_list[i].c_str());
    }

    std::string tf = file_list[0];

    if (use_file_list)
        tf = list_name;

    size_t ipos = tf.rfind('/', tf.length());

    if (ipos != std::string::npos)
        tf = tf.substr(ipos + 1, tf.length() - ipos);

    size_t ppos = tf.rfind('.', tf.length());

    if (ppos != std::string::npos)
        tf = tf.substr(0, ppos);

    std::string outf = tf + ".evio";

    PRadDigitization *prad_digi = new PRadDigitization(t, outf);
    prad_digi->RegisterDet(new HyCalDigitization("HC", "config/hycal.conf", run, sMode));

    int N = t->GetEntries();

    for (int i = 0; i < N; i++) {
        if ((i + 1) % 1000 == 0)
            std::cout << i + 1 << " events processed" << std::endl;

        t->GetEntry(i);
        //prad_digi->Print();
        prad_digi->ProcessEvent();
    }

    delete prad_digi;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
