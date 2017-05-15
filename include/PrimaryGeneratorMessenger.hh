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
// PrimaryGeneratorMessenger.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"

#include "G4String.hh"

class PrimaryGeneratorAction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger: public G4UImessenger
{
public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction *);
    virtual ~PrimaryGeneratorMessenger();

    void SetNewValue(G4UIcommand *, G4String);

private:
    PrimaryGeneratorAction    *Action;

    G4UIdirectory             *GunDir;
    G4UIcmdWithAString        *GunTypeCmd;
    G4UIcmdWithAString        *EventTypeCmd;
    G4UIcmdWithAString        *RecoilCmd;
    G4UIcmdWithADoubleAndUnit *EBeamCmd;
    G4UIcmdWith3VectorAndUnit *PosCmd;
    G4UIcmdWithADoubleAndUnit *ThetaCmd;
    G4UIcmdWithADoubleAndUnit *PhiCmd;
    G4UIcmdWithADoubleAndUnit *ThetaLowCmd;
    G4UIcmdWithADoubleAndUnit *ThetaHighCmd;
    G4UIcmdWithADoubleAndUnit *EnpLowCmd;
    G4UIcmdWithADoubleAndUnit *EnpHighCmd;
    G4UIcmdWithAString        *EventFileCmd;
    G4UIcmdWithAString        *TargetProfileCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
