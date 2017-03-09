//
// PRadDigitization.cc
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PRadDigitization.hh"

#include "StandardDigiBase.hh"

#include "evio.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TChain.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadDigitization::PRadDigitization(TChain *t, const std::string &filename) : fPreStart(false), fEventNumber(0), fPRadOut(-1), fData(t)
{
    OpenFile(filename);

    fDigitizer.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PRadDigitization::~PRadDigitization()
{
    for (auto &digi : fDigitizer)
        delete digi;

    uint32_t now = time(NULL);
    uint32_t end[5] = {0x00000004, 0x002001cc, now, fEventNumber, 0x00000000};
    evWrite(fPRadOut, end);
    evClose(fPRadOut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadDigitization::RegisterDet(StandardDigiBase *digi)
{
    digi->RegisterData(fData);

    fDigitizer.push_back(digi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadDigitization::PreStart()
{
    int index = 0;

    // event header
    fPRadBuffer[index++] = 0x00000000;
    fPRadBuffer[index++] = 0x008110cc;
    fEventNumberIndex = index + 2;
    index += AddEventInfoBank(&fPRadBuffer[index]);

    for (auto &digi : fDigitizer)
        index += digi->PreStart(fPRadBuffer, index);

    fPRadBuffer[0] = index - 1;

    fPreStart = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadDigitization::ProcessEvent()
{
    if (!fPreStart)
        PreStart();

    bool trigger = true;

    for (auto &digi : fDigitizer)
        trigger = trigger && digi->ProcessEvent(fPRadBuffer);

    fPRadBuffer[fEventNumberIndex] = ++fEventNumber;

    if (evWrite(fPRadOut, fPRadBuffer) != S_SUCCESS)
        std::cerr << "ERROR: cannot write event to output file" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadDigitization::Print() const
{
    for (auto &digi : fDigitizer)
        digi->Print();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PRadDigitization::OpenFile(const std::string &filename)
{
    std::string path = "output/" + filename;
    char mode[] = "w";
    char outf[256];

    strcpy(outf, path.c_str());

    int status = evOpen(outf, mode, &fPRadOut);

    if (status != S_SUCCESS)
        std::cerr << "ERROR: cannot open output file \"" << outf << "\"" << " (" << "0x" << std::hex << std::setw(8) << std::setfill('0') << status << ")" << std::endl;

    uint32_t now = time(NULL);
    uint32_t prestart[5] = {0x00000004, 0x001101cc, now, 0x00000260, 0x00000000};
    uint32_t go[5] = {0x00000004, 0x001201cc, now + 1, 0x00000000, 0x00000000};

    evWrite(fPRadOut, prestart);
    evWrite(fPRadOut, go);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int PRadDigitization::AddEventInfoBank(uint32_t *buffer)
{
    int index = 0;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0xc0000100;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x00000000;


    buffer[0] = index - 1;
    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
