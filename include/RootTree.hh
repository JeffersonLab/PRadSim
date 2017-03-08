//
// RootTree.hh
// Developer : Chao Gu
// History:
//   Jan 2017, C. Gu, Add for ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RootTree_h
#define RootTree_h 1

#include <map>

//maximum number of hits in a SD
#define MaxSDHits 5000

class TFile;
class TTree;

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

class RootTree
{
public:
    RootTree(const char *filename);
    virtual ~RootTree();

    void Initialize(const char *filename);
    void RegisterSD(const char *sdname);
    void UpdateValue(const char *sdname, int pid, int tid, int ptid, double inx, double iny, double inz, double inp, double outx, double outy, double outz, double outp, double edep, double time, int copyno);
    void FillTree(); // fill tree

private:
    void Reset();

    TFile *file;
    TTree *tree; // hits info, event

    std::map <const char *, SDData *> SDMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
