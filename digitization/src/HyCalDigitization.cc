//
// HyCalDigitization.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HyCalDigitization.hh"

#include "PRadClusterProfile.h"
#include "PRadEventStruct.h"
#include "PRadHyCalSystem.h"
#include "PRadHyCalModule.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TChain.h"
#include "TRandom2.h"

#include <ctime>
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static TRandom2 *RandGen = new TRandom2();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalDigitization::HyCalDigitization(const std::string &abbrev, const std::string &path, const std::string &method) : StandardDigiBase(abbrev), fDMethod(0)
{
    if (method == "cluster") fDMethod = 1;

    RandGen->SetSeed((UInt_t)time(NULL));

    fHyCal = new PRadHyCalSystem(path);

    fModuleList = fHyCal->GetModuleList();

    if (fModuleList.size() != NModules)
        std::cout << "ERROR: number of modules do not match" << std::endl;

    fModuleHitList.clear();

    for (auto &module : fModuleList)
        fModuleHitList.push_back(ModuleHit(module->GetID(), module->GetGeometry(), module->GetLayout(), 0, false));

    fProfile = &PRadClusterProfile::Instance();

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 0;
        fModuleTrackL[i] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalDigitization::~HyCalDigitization()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::RegisterData(TChain *t)
{
    StandardDigiBase::RegisterData(t);

    t->SetBranchAddress(Form("%s.TotalEdep", fAbbrev), &fTotalEdep);

    if (fDMethod != 1) {
        t->SetBranchAddress(Form("%s.ModuleEdep", fAbbrev), fModuleEdep);
        t->SetBranchAddress(Form("%s.ModuleTrackL", fAbbrev), fModuleTrackL);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::PreStart(uint32_t *buffer, int base_index)
{
    int index = 0;

    for (int roc_id = 6; roc_id >= 4; --roc_id)
        index += addRocData(&buffer[index + base_index], roc_id, index + base_index);

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool HyCalDigitization::ProcessEvent(uint32_t *buffer)
{
    if (fTotalEdep < TRIGGER_THRESHOLD) {
        Clear();
        return false;
    }

    if (fDMethod == 1) UpdateEnergy();

    for (int i = 0; i < NModules; i++)
        FillBuffer(buffer, *(fModuleList[i]), fModuleEdep[i]);

    Clear();
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::Clear()
{
    StandardDigiBase::Clear();

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 0;
        fModuleTrackL[i] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HyCalDigitization::UpdateEnergy()
{
    for (int i = 0; i < fN; i++) {
        for (int j = 0; j < NModules; j++) {
            float fx = fX[i];
            float fy = fY[i];
            float frac = fProfile->GetProfile(fx, fy, fModuleHitList[j]).frac;
            fModuleEdep[j] += fMomentum[i] * frac;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::addRocData(uint32_t *buffer, int roc_id, int base_index)
{
    int index = 0;
    int nslot, slot[25];

    switch (roc_id) {
    default:
        return 0;

    case 4:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 22 - 2 * i;

        break;

    case 5:
    case 6:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 23 - 2 * i;

        break;
    }

    // add roc header
    buffer[index++] = 0x00000000;
    buffer[index++] = (roc_id << 16) | 0x1001;

    // add TI bank 11 words
    buffer[index++] = 0x0000000a; // 10 + 1 words in total
    buffer[index++] = 0xe10a0100; // TI bank header
    buffer[index + 2] = 2 << 24; // only 2nd word matters, it defines trigger type, here it is total sum
    index += 9; // TI bank expects 9 words

    buffer[index++] = 0x00000000;
    buffer[index++] = 0xe1200100; // Fastbus bank header
    // roc id and board number
    buffer[index++] = 0xdc0adc00 | ((roc_id & 0xf) << 20) | (nslot & 0xff);

    for (int i = 0; i < nslot; ++i) {
        buffer[index++] = (slot[i] << 27) | 65;
        data_index[(6 - roc_id) * 10 + i] = index + base_index;

        for (int ch = 0; ch < 64; ++ch)
            buffer[index++] = (slot[i] << 27) | (ch << 17);
    }

    buffer[index++] = 0xfabc0005; // end word
    buffer[0] = index - 1; // roc bank size
    buffer[13] = index - 14; // fastbus bank size

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::FillBuffer(uint32_t *buffer, const PRadHyCalModule &module, double edep)
{
    int crate = module.GetChannel()->GetAddress().crate;
    int slot = module.GetChannel()->GetAddress().slot;
    int channel = module.GetChannel()->GetAddress().channel;

    int pos = (6 - crate) * 10 + ((23 - slot) / 2);

    int index = data_index[pos] + channel;

    double ped = RandGen->Gaus(module.GetChannel()->GetPedestal().mean, module.GetChannel()->GetPedestal().sigma);
    unsigned short val = 0;

    if (!module.GetChannel()->IsDead())
        val = ped + edep / module.GetCalibrationFactor();
    else
        val = ped;

    buffer[index] = (slot << 27) | (channel << 17) | val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
