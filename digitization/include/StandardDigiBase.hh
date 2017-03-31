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

#define MaxNHits 300 //maximum number of hits in a SD

class TChain;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StandardDigiBase
{
public:
    StandardDigiBase(const std::string &abbrev);
    virtual ~StandardDigiBase();

    virtual void RegisterData(TChain *t);

    virtual int PreStart(uint32_t *buffer, int base_index) = 0;
    virtual bool ProcessEvent(uint32_t *buffer) = 0;

    virtual void Clear();
    virtual void Print() const;

protected:
    const char *fAbbrev;

    int fN;
    int fPID[MaxNHits]; // Particle ID
    int fTID[MaxNHits]; // Track ID
    int fPTID[MaxNHits]; // Parent Track ID
    int fDID[MaxNHits];
    double fX[MaxNHits];
    double fY[MaxNHits];
    double fZ[MaxNHits];
    double fMomentum[MaxNHits];
    double fTheta[MaxNHits];
    double fPhi[MaxNHits];
    double fTime[MaxNHits];
    double fEdep[MaxNHits];
    double fTrackL[MaxNHits];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
