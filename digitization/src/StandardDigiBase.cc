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

#include "evio.h"

#include <iomanip>
#include <iostream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDigiBase::StandardDigiBase(const std::string &abbrev) : fAbbrev(abbrev.c_str())
{
    fN = 0;

    for (int i = 0; i < MaxNHits; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fDID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
        fTime[i] = 1e+38;
        fEdep[i] = 1e+38;
        fTrackL[i] = 1e+38;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StandardDigiBase::~StandardDigiBase()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::RegisterData(TChain *t)
{
    t->SetBranchAddress(Form("%s.N", fAbbrev), &fN);
    t->SetBranchAddress(Form("%s.PID", fAbbrev), fPID);
    t->SetBranchAddress(Form("%s.TID", fAbbrev), fTID);
    t->SetBranchAddress(Form("%s.PTID", fAbbrev), fPTID);
    t->SetBranchAddress(Form("%s.X", fAbbrev), fX);
    t->SetBranchAddress(Form("%s.Y", fAbbrev), fY);
    t->SetBranchAddress(Form("%s.Z", fAbbrev), fZ);
    t->SetBranchAddress(Form("%s.P", fAbbrev), fMomentum);
    t->SetBranchAddress(Form("%s.Theta", fAbbrev), fTheta);
    t->SetBranchAddress(Form("%s.Phi", fAbbrev), fPhi);
    t->SetBranchAddress(Form("%s.Time", fAbbrev), fTime);
    t->SetBranchAddress(Form("%s.Edep", fAbbrev), fEdep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::Clear()
{
    for (int i = 0; i < fN; i++) {
        fPID[i] = -9999;
        fTID[i] = -9999;
        fPTID[i] = -9999;
        fDID[i] = -9999;
        fX[i] = 1e+38;
        fY[i] = 1e+38;
        fZ[i] = 1e+38;
        fMomentum[i] = 1e+38;
        fTheta[i] = 1e+38;
        fPhi[i] = 1e+38;
        fTime[i] = 1e+38;
        fEdep[i] = 1e+38;
        fTrackL[i] = 1e+38;
    }

    fN = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandardDigiBase::Print() const
{
    std::cout << fAbbrev << " : " << fN << std::endl;

    int prec = std::cout.precision(3);

    for (int i = 0; i < fN; i++) {
        std::cout << std::setw(5) << fPID[i] << " " << std::setw(5) << fTID[i] << " " << std::setw(5) << fPTID[i] << " ";
        std::cout.setf(std::ios::fixed);
        std::cout.precision(1);
        std::cout << std::setw(6) << fX[i] << " " << std::setw(6) << fY[i] << " " << std::setw(7) << fZ[i] << " ";
        std::cout.unsetf(std::ios::fixed);
        std::cout.precision(3);
        std::cout << std::setw(8) << fMomentum[i] << " " << std::setw(8) << fEdep[i] << " ";
        std::cout << std::setw(5) << fTime[i] << " ";
        std::cout << std::endl;
    }

    std::cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
