//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// PRadSim.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootTree.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>

#ifdef G4VIS_USE
    #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
    #include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree *gRootTree = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void usage(int, char **argv)
{
    printf("usage: %s [options] [MACRO_NAME]\n", argv[0]);
    printf("  -c, --conf=prad          Set configuration\n");
    printf("  -h, --help                 Print usage\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc, char **argv)
{
    std::string conf = "prad";
    std::string macro;
    macro.clear();

    while (1) {
        static struct option long_options[] = {
            {"help",  no_argument, 0, 'h'},
            {"conf", required_argument, 0, 'c'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "c:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'c':
            conf = optarg;
            break;

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
        macro = argv[optind++];

    // Choose the Random engine
    //
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    // Set the Random seeds
    //
    long seeds[2];
    int index;
    time_t systime = time(NULL);
    index = (int)systime;
    seeds[0] = (long)systime;
    seeds[1] = (long)(systime * G4UniformRand());
    CLHEP::HepRandom::setTheSeeds(seeds, index);

    // User Verbose output class
    //
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);

    // Initialize output root tree
    //
    std::ifstream data_file("output/file.output");
    std::string file_name;
    int run_number = 1;

    if (data_file.is_open()) {
        std::string line;

        while (std::getline(data_file, line)) {
            if (line.at(0) == '#')
                continue;

            std::stringstream iss(line);
            iss >> file_name >> run_number;
            break;
        }

        data_file.close();
    }

    if (file_name.empty()) file_name = "simrun";

    std::ofstream data_file_o("output/file.output");
    data_file_o << file_name << "  " << run_number + 1;
    data_file_o.close();

    std::string path = "output/" + file_name + "_" + std::to_string(run_number) + ".root";
    char outf[256];

    strcpy(outf, path.c_str());

    gRootTree = new RootTree(outf);

    // Construct the default run manager
    //
    G4RunManager *runManager = new G4RunManager;

    // Set mandatory initialization classes
    //
    G4String confstr = conf;
    DetectorConstruction *detector = new DetectorConstruction(confstr);
    runManager->SetUserInitialization(detector);
    //
    PhysicsList *physics = new PhysicsList;
    runManager->SetUserInitialization(physics);

    // Set user action classes
    //
    PrimaryGeneratorAction *gen_action = new PrimaryGeneratorAction(detector);
    runManager->SetUserAction(gen_action);
    //
    RunAction *run_action = new RunAction;
    runManager->SetUserAction(run_action);
    //
    EventAction *event_action = new EventAction(run_action);
    runManager->SetUserAction(event_action);
    //
    SteppingAction *stepping_action = new SteppingAction(detector, event_action);
    runManager->SetUserAction(stepping_action);

    // Initialize G4 kernel
    //
    runManager->Initialize();

#ifdef G4VIS_USE
    // Initialize visualization
    G4VisManager *visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
#endif

    // Get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (!macro.empty()) { // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macro);
    } else {
        // interactive mode : define UI session
#ifdef G4UI_USE
        G4UIExecutive *ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");
#endif

        if (ui->IsGUI())
            UImanager->ApplyCommand("/control/execute gui.mac");

        ui->SessionStart();
        delete ui;
#endif
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
#ifdef G4VIS_USE
    delete visManager;
#endif
    delete runManager;

    delete gRootTree;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
