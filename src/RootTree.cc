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

#include <iostream>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::RootTree(const char *filename)
{
    Initialize(filename);

    SDMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::Initialize(const char *filename)
{
    file = new TFile(filename, "RECREATE");
    tree = new TTree("T", "Simulation Results");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::RegisterSD(const char *sdname)
{
    if (SDMap.count(sdname) > 0) {
        std::cerr << "RootTree: trying to register a sensitive detector \"" << sdname << "\", but the sensitive detector list already has it." << std::endl;
        exit(1);
    }

    SDData *newsd = new SDData();
    SDMap[sdname] = newsd;

    tree->Branch(Form("%s.N", sdname), &newsd->N, Form("%s.N/I", sdname));
    tree->Branch(Form("%s.PID", sdname), newsd->PID, Form("%s.PID[%s.N]/I", sdname, sdname));
    tree->Branch(Form("%s.TID", sdname), newsd->TID, Form("%s.TID[%s.N]/I", sdname, sdname));
    tree->Branch(Form("%s.PTID", sdname), newsd->PTID, Form("%s.PTID[%s.N]/I", sdname, sdname));
    tree->Branch(Form("%s.In.X", sdname), newsd->InPosX, Form("%s.In.X[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.In.Y", sdname), newsd->InPosY, Form("%s.In.Y[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.In.Z", sdname), newsd->InPosZ, Form("%s.In.Z[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.In.P", sdname), newsd->InMom, Form("%s.In.P[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Out.X", sdname), newsd->OutPosX, Form("%s.Out.X[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Out.Y", sdname), newsd->OutPosY, Form("%s.Out.Y[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Out.Z", sdname), newsd->OutPosZ, Form("%s.Out.Z[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Out.P", sdname), newsd->OutMom, Form("%s.Out.P[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Edep", sdname), newsd->Edep, Form("%s.Edep[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.Time", sdname), newsd->Time, Form("%s.Time[%s.N]/D", sdname, sdname));
    tree->Branch(Form("%s.DID", sdname), newsd->CopyNo, Form("%s.DID[%s.N]/I", sdname, sdname));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::~RootTree()
{
    for (std::map<const char *, SDData *>::iterator it = SDMap.begin(); it != SDMap.end(); ++it) {
        SDData *sd = it->second;

        delete sd;
    }

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

void RootTree::UpdateValue(const char *sdname, int pid, int tid, int ptid, double inx, double iny, double inz, double inp, double outx, double outy, double outz, double outp, double edep, double time, int copyno)
{
    SDData *sd = SDMap[sdname];

    if (sd->N < MaxSDHits) {
        sd->PID[sd->N] = pid;
        sd->TID[sd->N] = tid;
        sd->PTID[sd->N] = ptid;
        sd->InPosX[sd->N] = inx;
        sd->InPosY[sd->N] = iny;
        sd->InPosZ[sd->N] = inz;
        sd->InMom[sd->N] = inp;
        sd->OutPosX[sd->N] = outx;
        sd->OutPosY[sd->N] = outy;
        sd->OutPosZ[sd->N] = outz;
        sd->OutMom[sd->N] = outp;
        sd->Edep[sd->N] = edep;
        sd->Time[sd->N] = time;
        sd->CopyNo[sd->N] = copyno;

        sd->N++;
    } else
        std::cout << "RootTree: too many hits in sensitive detector \"" << sdname << "\"." << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::Reset()
{
    for (std::map<const char *, SDData *>::iterator it = SDMap.begin(); it != SDMap.end(); ++it) {
        SDData *sd = it->second;

        for (int i = 0; i < sd->N; i++) {
            sd->PID[i] = -999;
            sd->TID[i] = -999;
            sd->PTID[i] = -999;
            sd->InPosX[i] = 1e+38;
            sd->InPosY[i] = 1e+38;
            sd->InPosZ[i] = 1e+38;
            sd->InMom[i] = 1e+38;
            sd->OutPosX[i] = 1e+38;
            sd->OutPosY[i] = 1e+38;
            sd->OutPosZ[i] = 1e+38;
            sd->OutMom[i] = 1e+38;
            sd->Edep[i] = 1e+38;
            sd->Time[i] = 1e+38;
            sd->CopyNo[i] = -999;
        }

        sd->N = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
