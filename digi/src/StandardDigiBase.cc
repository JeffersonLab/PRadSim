//
// StandardDigiBase.cc
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StandardDigiBase.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TChain.h"

#include <iomanip>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDigiBase::StandardDigiBase(const char *name) : fAbbrev(name)
{
    Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDigiBase::~StandardDigiBase()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::RegisterData(TChain *t)
{
    t->SetBranchAddress(Form("%s.N", fAbbrev), &fData.N);
    t->SetBranchAddress(Form("%s.PID", fAbbrev), fData.PID);
    t->SetBranchAddress(Form("%s.TID", fAbbrev), fData.TID);
    t->SetBranchAddress(Form("%s.PTID", fAbbrev), fData.PTID);
    t->SetBranchAddress(Form("%s.In.X", fAbbrev), fData.InPosX);
    t->SetBranchAddress(Form("%s.In.Y", fAbbrev), fData.InPosY);
    t->SetBranchAddress(Form("%s.In.Z", fAbbrev), fData.InPosZ);
    t->SetBranchAddress(Form("%s.In.P", fAbbrev), fData.InMom);
    t->SetBranchAddress(Form("%s.Out.X", fAbbrev), fData.OutPosX);
    t->SetBranchAddress(Form("%s.Out.Y", fAbbrev), fData.OutPosY);
    t->SetBranchAddress(Form("%s.Out.Z", fAbbrev), fData.OutPosZ);
    t->SetBranchAddress(Form("%s.Out.P", fAbbrev), fData.OutMom);
    t->SetBranchAddress(Form("%s.Edep", fAbbrev), fData.Edep);
    t->SetBranchAddress(Form("%s.Time", fAbbrev), fData.Time);
    t->SetBranchAddress(Form("%s.DID", fAbbrev), fData.CopyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::Clear()
{
    fData.N = 0;

    for (int i = 0; i < MaxSDHits; i++) {
        fData.PID[i] = -999;
        fData.TID[i] = -999;
        fData.PTID[i] = -999;
        fData.InPosX[i] = 1e+38;
        fData.InPosY[i] = 1e+38;
        fData.InPosZ[i] = 1e+38;
        fData.InMom[i] = 1e+38;
        fData.OutPosX[i] = 1e+38;
        fData.OutPosY[i] = 1e+38;
        fData.OutPosZ[i] = 1e+38;
        fData.OutMom[i] = 1e+38;
        fData.Edep[i] = 1e+38;
        fData.Time[i] = 1e+38;
        fData.CopyNo[i] = -999;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::Print()
{
    std::cout << fAbbrev << " : " << fData.N << std::endl;

    for (int i = 0; i < fData.N; i++) {
        int prec = std::cout.precision(3);
        std::cout << std::setw(5) << fData.PID[i] << " " << std::setw(5) << fData.TID[i] << " " << std::setw(5) << fData.PTID[i] << " ";
        std::cout.setf(std::ios::fixed);
        std::cout.precision(1);
        std::cout << std::setw(6) << fData.InPosX[i] << " " << std::setw(6) << fData.InPosY[i] << " " << std::setw(7) << fData.InPosZ[i] << " ";
        std::cout << std::setw(6) << fData.OutPosX[i] << " " << std::setw(6) << fData.OutPosY[i] << " " << std::setw(7) << fData.OutPosZ[i] << " ";
        std::cout.unsetf(std::ios::fixed);
        std::cout.precision(3);
        std::cout << std::setw(8) << fData.InMom[i] << " " << std::setw(8) << fData.OutMom[i] << " " << std::setw(8) << fData.Edep[i] << " ";
        std::cout << std::setw(5) << fData.Time[i] << " " << std::setw(5) << fData.CopyNo[i] << std::endl;
        std::cout.precision(prec);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
