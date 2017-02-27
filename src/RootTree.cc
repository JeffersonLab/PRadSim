//
// RootTree.cc
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RootTree.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

#include <cstring>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::RootTree()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::RootTree(const char *filename) : SD_N(0)
{
    Initialize(filename);

    Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::Initialize(const char *filename)
{
    file = new TFile(filename, "RECREATE");
    tree = new TTree("T", "Simulation Results");

    tree->Branch("N", &SD_N, "N/I");
    tree->Branch("PID", SD_PID, "PID[N]/I");
    tree->Branch("TID", SD_TID, "TID[N]/I");
    tree->Branch("PTID", SD_PTID, "PTID[N]/I");
    tree->Branch("X", SD_X, "X[N]/D");
    tree->Branch("Y", SD_Y, "Y[N]/D");
    tree->Branch("Z", SD_Z, "Z[N]/D");
    tree->Branch("P", SD_P, "P[N]/D");
    tree->Branch("Theta", SD_Theta, "Theta[N]/D");
    tree->Branch("Phi", SD_Phi, "Phi[N]/D");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::~RootTree()
{
    file->Write("", TObject::kOverwrite);
    file->Close();
    file->Delete();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::FillTree()
{
    tree->Fill();

    Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::UpdateValue(int pid, int tid, int ptid, double x, double y, double z, double p, double theta, double phi)
{
    if (SD_N < MaxSDHit) {
        SD_PID[SD_N] = pid;
        SD_TID[SD_N] = tid;
        SD_PTID[SD_N] = ptid;
        SD_X[SD_N] = x;
        SD_Y[SD_N] = y;
        SD_Z[SD_N] = z;
        SD_P[SD_N] = p;
        SD_Theta[SD_N] = theta;
        SD_Phi[SD_N] = phi;
        //SD_Edep[SD_N] = 1e+38;
        //SD_NonIonEdep[SD_N] = 1e+38;

        SD_N++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::Reset()
{
    for (int i = 0; i < SD_N; i++) {
        SD_PID[i] = -999;
        SD_TID[i] = -999;
        SD_PTID[i] = -999;
        SD_X[i] = 1e+38;
        SD_Y[i] = 1e+38;
        SD_Z[i] = 1e+38;
        SD_P[i] = 1e+38;
        SD_Theta[i] = 1e+38;
        SD_Phi[i] = 1e+38;
        //SD_Edep[i] = 1e+38;
        //SD_NonIonEdep[i] = 1e+38;
    }

    SD_N = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
