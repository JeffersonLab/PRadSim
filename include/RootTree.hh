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

//maximum number of hits in a SD
#define MaxSDHit 1024

class TFile;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RootTree
{
public:
    RootTree(const char *filename);
    virtual ~RootTree();

    void Initialize(const char *filename);
    void UpdateValue(int pid, int tid, int ptid, double x, double y, double z, double p, double theta, double phi);

    void FillTree(); // fill tree

private:
    RootTree();

    void Reset();

    //for general sensitive detectors
    int SD_N; //must be initialized before using
    int SD_PID[MaxSDHit]; // Particle ID
    int SD_TID[MaxSDHit]; // Track ID
    int SD_PTID[MaxSDHit]; // Parent Track ID
    double SD_X[MaxSDHit];
    double SD_Y[MaxSDHit];
    double SD_Z[MaxSDHit];
    double SD_P[MaxSDHit];
    double SD_Theta[MaxSDHit];
    double SD_Phi[MaxSDHit];

private:
    TFile *file;
    TTree *tree; // hits info, event
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
