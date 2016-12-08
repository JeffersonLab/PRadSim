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

#include "RootTree.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

#include <cstring>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::RootTree()
{
    Initialize();
    
    Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::Initialize()
{
    file = new TFile("test.root", "RECREATE");
    tree = new TTree("T", "PRad Simulation");

    tree->Branch("N", &SD_N, "N/I");
    tree->Branch("PID", SD_PID, "PID[N]/I");
    tree->Branch("TID", SD_TID, "TID[N]/I");
    tree->Branch("PTID", SD_PTID, "PTID[N]/I");
    tree->Branch("X", SD_X, "X[N]/D");
    tree->Branch("Y", SD_Y, "Y[N]/D");
    tree->Branch("Z", SD_Z, "Z[N]/D");
    tree->Branch("P", SD_P, "P[N]/D");
    tree->Branch("Theta", SD_Theta, "Theta[N]/D");
    tree->Branch("Phi", SD_Phi, "Phi[N]/D");
    //tree->Branch("Edep", SD_Edep, "Edep[N]/D");
    //tree->Branch("NonIoEdep", SD_NonIonEdep, "NonIonEdep[N]/D");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree::~RootTree()
{
    file->Write("", TObject::kOverwrite);
    file->Close();
    file->Delete();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::FillTree()
{
    tree->Fill();
    
    Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootTree::UpdateValue(int pid, int tid, int ptid, double x, double y, double z, double p, double theta, double phi)
{
    if (SD_N < MaxSDHit) {
        SD_PID[SD_N] = pid;
        SD_TID[SD_N] = tid;
        SD_PTID[SD_N] = ptid;
        SD_X[SD_N] = x;
        SD_Y[SD_N] = y;
        SD_Z[SD_N] = z;
        SD_P[SD_N] = p;
        SD_Theta[SD_N] = theta;
        SD_Phi[SD_N] = phi;
        //SD_Edep[SD_N] = 1e+38;
        //SD_NonIonEdep[SD_N] = 1e+38;

        SD_N++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RootTree::Reset()
{
    for (int i = 0; i < SD_N; i++) {
        SD_PID[i] = -999;
        SD_TID[i] = -999;
        SD_PTID[i] = -999;
        SD_X[i] = 1e+38;
        SD_Y[i] = 1e+38;
        SD_Z[i] = 1e+38;
        SD_P[i] = 1e+38;
        SD_Theta[i] = 1e+38;
        SD_Phi[i] = 1e+38;
        //SD_Edep[i] = 1e+38;
        //SD_NonIonEdep[i] = 1e+38;
    }

    SD_N = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
