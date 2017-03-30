//
// HyCalDigitization.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HyCalDigitization_h
#define HyCalDigitization_h 1

#include "StandardDigiBase.hh"

#include "evio.h"

#include <string>
#include <vector>

#define NModules 1728
#define TRIGGER_THRESHOLD 500 // MeV

class PRadHyCalModule;
class PRadHyCalSystem;

class TChain;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HyCalDigitization : public StandardDigiBase
{
public:
    HyCalDigitization(const std::string &name, const std::string &path);
    virtual ~HyCalDigitization();

    void RegisterData(TChain *t);

    int PreStart(uint32_t *buffer, int base_index);
    bool ProcessEvent(uint32_t *buffer);

    void Clear();

private:
    int addRocData(uint32_t *buffer, int roc_id, int base_index);
    void FillBuffer(uint32_t *buffer, const PRadHyCalModule &module, double edep);

    double fTotalEdep;
    double fModuleEdep[NModules];
    double fModuleTrackL[NModules];

    int data_index[30];

    std::vector<PRadHyCalModule *> fModuleList;

    PRadHyCalSystem *fHyCal;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
