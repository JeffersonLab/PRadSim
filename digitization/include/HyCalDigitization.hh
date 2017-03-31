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

#include "PRadEventStruct.h"

#include <string>
#include <vector>

#define NModules 1728
#define TRIGGER_THRESHOLD 500 // MeV

class PRadClusterProfile;
class PRadHyCalModule;
class PRadHyCalSystem;

class TChain;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HyCalDigitization : public StandardDigiBase
{
public:
    HyCalDigitization(const std::string &abbrev, const std::string &path, const std::string &method);
    virtual ~HyCalDigitization();

    void RegisterData(TChain *t);

    int PreStart(uint32_t *buffer, int base_index);
    bool ProcessEvent(uint32_t *buffer);

    void Clear();

private:
    void UpdateEnergy();

    int addRocData(uint32_t *buffer, int roc_id, int base_index);
    void FillBuffer(uint32_t *buffer, const PRadHyCalModule &module, double edep);

    int fDMethod;

    double fTotalEdep;
    double fModuleEdep[NModules];
    double fModuleTrackL[NModules];

    int data_index[30];

    std::vector<ModuleHit> fModuleHitList;
    std::vector<PRadHyCalModule *> fModuleList;

    PRadHyCalSystem *fHyCal;
    PRadClusterProfile *fProfile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
