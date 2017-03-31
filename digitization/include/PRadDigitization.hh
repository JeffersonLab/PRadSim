//
// PRadDigitization.hh
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PRadDigitization_h
#define PRadDigitization_h 1

#include "StandardDigiBase.hh"

#include "evio.h"

#include <string>
#include <vector>

#define MAX_PRAD_BUFFER 3000

class TChain;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PRadDigitization
{
public:
    PRadDigitization(TChain *t, const std::string &filename);
    virtual ~PRadDigitization();

    void RegisterDet(StandardDigiBase *digi);

    void PreStart();
    void ProcessEvent();

    void Print() const;

private:
    void OpenFile(const std::string &filename);

    int AddEventInfoBank(uint32_t *buffer);

    bool fPreStart;

    uint32_t fEventNumber;
    int fEventNumberIndex;

    int fPRadOut;
    uint32_t fPRadBuffer[MAX_PRAD_BUFFER];

private:
    TChain *fData;

    std::vector<StandardDigiBase *> fDigitizer;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
