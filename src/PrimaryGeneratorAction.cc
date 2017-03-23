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
// PrimaryGeneratorAction.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Rewrite with class PrimaryGenerator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGenerator.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String conf) : G4VUserPrimaryGeneratorAction(), fConfig(conf), fGunType("ring"), fE(0), fThetaLo(0), fThetaHi(0), fPrimaryGenerator(NULL)
{
    if (fConfig != "prad" && fConfig != "drad")
        fConfig = "prad";

    fRecoil.clear();
    fEventFile.clear();

    // create a messenger for this class
    gunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    if (fPrimaryGenerator)
        delete fPrimaryGenerator;

    delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    bool recoilon = false;
    if (!fRecoil.empty() && fRecoil != "none") recoilon = true;
    
    if (!fPrimaryGenerator) {
        if (fConfig == "prad") {
            if (fGunType == "ring")
                fPrimaryGenerator = new PrimaryGenerator(fE, fThetaLo, fThetaHi, false, "proton");
            else
                fPrimaryGenerator = new PRadPrimaryGenerator(fGunType, false, "proton", fEventFile);
        } else {
            if (fGunType == "ring")
                fPrimaryGenerator = new PrimaryGenerator(fE, fThetaLo, fThetaHi, recoilon, fRecoil);
            else
                fPrimaryGenerator = new DRadPrimaryGenerator(fGunType, recoilon, fRecoil, fEventFile);
        }
    }

    fPrimaryGenerator->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
