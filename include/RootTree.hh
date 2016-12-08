//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: RootTree.cc,v 1.3, 2013/02/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
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
    RootTree();
    ~RootTree();

    void Initialize();
    void UpdateValue(int pid, int tid, int ptid, double x, double y, double z, double p, double theta, double phi);

    void FillTree(); // fill tree

private:
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
    //double SD_Edep[MaxSDHit]; // MeV
    //double SD_NonIonEdep[MaxSDHit]; // MeV

private:
    TFile *file;
    TTree *tree; // hits info, event
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
