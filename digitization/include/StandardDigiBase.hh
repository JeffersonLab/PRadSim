//
// StandardDigiBase.hh
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StandardDigiBase_h
#define StandardDigiBase_h 1

#include "evio.h"

#include <string>

//maximum number of hits in a SD
#define MaxSDHits 5000

class TChain;

struct SDData {
    int N;
    int PID[MaxSDHits]; // Particle ID
    int TID[MaxSDHits]; // Track ID
    int PTID[MaxSDHits]; // Parent Track ID
    double InPosX[MaxSDHits];
    double InPosY[MaxSDHits];
    double InPosZ[MaxSDHits];
    double OutPosX[MaxSDHits];
    double OutPosY[MaxSDHits];
    double OutPosZ[MaxSDHits];
    double InMom[MaxSDHits];
    double OutMom[MaxSDHits];
    double Edep[MaxSDHits];
    double Time[MaxSDHits];
    int CopyNo[MaxSDHits];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StandardDigiBase
{
public:
    StandardDigiBase(const std::string &name);
    virtual ~StandardDigiBase();

    virtual void RegisterData(TChain *t);

    virtual int PreStart(uint32_t *buffer, int base_index) = 0;
    virtual bool ProcessEvent(uint32_t *buffer) = 0;

    virtual void Clear();
    virtual void Print() const;

protected:
    const char *fAbbrev;

    SDData fData;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
