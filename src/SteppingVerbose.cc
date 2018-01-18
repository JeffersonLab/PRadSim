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
// SteppingVerbose.cc
// Developer : Geant4 Developers
// History:
//   Aug 2012, Copy from ExampleN02.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingVerbose.hh"

#include "G4SteppingVerbose.hh"

#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingVerbose::SteppingVerbose() : G4SteppingVerbose()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingVerbose::~SteppingVerbose()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingVerbose::StepInfo()
{
    CopyState();

    G4int prec = G4cout.precision(3);

    if (verboseLevel >= 1) {
        if (verboseLevel >= 4) VerboseTrack();

        if (verboseLevel >= 3) {
            G4cout << G4endl;
            G4cout << std::setw(5) << "Step#" << " " << std::setw(5) << "X" << "     " << std::setw(5) << "Y" << "     " << std::setw(5) << "Z" << "     " << std::setw(8) << "KineE" << "   " << std::setw(10) << "dEStep" << "  " << std::setw(8) << "StepLeng" << "  " << std::setw(8) << "TrakLeng" << "  " << std::setw(23) << "Volume" << " " << std::setw(8) << "Process" << G4endl;
        }

        G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " ";
        G4cout << std::setw(5) << G4BestUnit(fTrack->GetPosition().x(), "Length") << " " << std::setw(5) << G4BestUnit(fTrack->GetPosition().y(), "Length") << " " << std::setw(5) << G4BestUnit(fTrack->GetPosition().z(), "Length") << " ";
        G4cout << std::setw(5) << G4BestUnit(fTrack->GetKineticEnergy(), "Energy") << " " << std::setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy")  << " ";
        G4cout << std::setw(5) << G4BestUnit(fStep->GetStepLength(), "Length")  << " " << std::setw(5) << G4BestUnit(fTrack->GetTrackLength(), "Length") << " ";
        G4cout << std::setw(18) << fTrack->GetVolume()->GetName() << " ";

        if (fTrack->GetVolume()->GetName().contains("Absorber")) {
            G4int DetectorID = 0;
            G4TouchableHandle theTouchable = fStep->GetPreStepPoint()->GetTouchableHandle();

            for (G4int i = 0; i < theTouchable->GetHistoryDepth(); i++)
                DetectorID += theTouchable->GetCopyNumber(i);

            G4cout << std::setw(5) << DetectorID << " ";
        } else
            G4cout << std::setw(5) << fStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() << " ";

        const G4VProcess *process = fStep->GetPostStepPoint()->GetProcessDefinedStep();
        G4String procName = "limit";

        if (process) procName = process->GetProcessName();

        if (fStepStatus == fWorldBoundary) procName = "out";

        if (procName == "Transportation") procName = "trans";

        G4cout << std::setw(8) << procName;
        G4cout << G4endl;
    }

    G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingVerbose::TrackingStarted()
{
    CopyState();

    G4int prec = G4cout.precision(3);

    if (verboseLevel >= 1) {
        G4cout << std::setw(5) << "Step#" << " " << std::setw(5) << "X" << "     " << std::setw(5) << "Y" << "     " << std::setw(5) << "Z" << "     " << std::setw(8) << "KineE" << "   " << std::setw(10) << "dEStep" << "  " << std::setw(8) << "StepLeng" << "  " << std::setw(8) << "TrakLeng" << "  " << std::setw(23) << "Volume" << " " << std::setw(8) << "Process" << G4endl;

        G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " ";
        G4cout << std::setw(5) << G4BestUnit(fTrack->GetPosition().x(), "Length") << " " << std::setw(5) << G4BestUnit(fTrack->GetPosition().y(), "Length") << " " << std::setw(5) << G4BestUnit(fTrack->GetPosition().z(), "Length") << " ";
        G4cout << std::setw(5) << G4BestUnit(fTrack->GetKineticEnergy(), "Energy") << " " << std::setw(8) << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy")  << " ";
        G4cout << std::setw(5) << G4BestUnit(fStep->GetStepLength(), "Length")  << " " << std::setw(5) << G4BestUnit(fTrack->GetTrackLength(), "Length") << " ";
        G4cout << std::setw(18) << fTrack->GetVolume()->GetName() << " ";

        if (fTrack->GetVolume()->GetName().contains("Absorber")) {
            G4int DetectorID = 0;
            G4TouchableHandle theTouchable = fStep->GetPreStepPoint()->GetTouchableHandle();

            for (G4int i = 0; i < theTouchable->GetHistoryDepth(); i++)
                DetectorID += theTouchable->GetCopyNumber(i);

            G4cout << std::setw(5) << DetectorID << " ";
        } else
            G4cout << std::setw(5) << fStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() << " ";

        G4cout << std::setw(8) << "init";
        G4cout << G4endl;
    }

    G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
