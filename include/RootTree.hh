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

class TFile;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RootTree
{
public:
    RootTree(const char *filename);
    virtual ~RootTree();

    inline TTree *GetTree() const;
    void FillTree(); // fill tree

private:
    TFile *fFile;
    TTree *fTree;
};

inline TTree *RootTree::GetTree() const
{
    return fTree;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
