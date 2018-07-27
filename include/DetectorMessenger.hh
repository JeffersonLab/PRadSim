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
// DetectorMessenger.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"

#include "G4String.hh"

class DetectorConstruction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
public:
    DetectorMessenger(DetectorConstruction *);
    ~DetectorMessenger();

    void SetNewValue(G4UIcommand *, G4String);

private:
    DetectorConstruction      *Detector;

    G4UIdirectory             *PRadSimDir;
    G4UIdirectory             *DetDir;
    G4UIdirectory             *ZDir;
    G4UIcmdWithADoubleAndUnit *TargetZCmd;
    G4UIcmdWithADoubleAndUnit *RecoilDetZCmd;
    G4UIcmdWithADoubleAndUnit *GEM1ZCmd;
    G4UIcmdWithADoubleAndUnit *GEM2ZCmd;
    G4UIcmdWithADoubleAndUnit *SciPlaneZCmd;
    G4UIcmdWithADoubleAndUnit *HyCalZCmd;
    G4UIdirectory             *TargetDir;
    G4UIcmdWithADoubleAndUnit *TargetRCmd;
    G4UIcmdWithADoubleAndUnit *TargetHalfLCmd;
    G4UIcmdWithAString        *TargetMatCmd;
    G4UIcmdWithADouble        *TargetDensityRatioCmd;
    G4UIdirectory             *RecoilDetDir;
    G4UIcmdWithAnInteger      *RecoilDetNSegCmd;
    G4UIcmdWithADoubleAndUnit *RecoilDetRCmd;
    G4UIcmdWithADoubleAndUnit *RecoilDetHalfLCmd;
    G4UIcmdWithADoubleAndUnit *RecoilDetL1ThicknessCmd;
    G4UIcmdWithADoubleAndUnit *RecoilDetL2ThicknessCmd;
    G4UIcmdWithADouble        *ExtDensityRatioCmd;
    G4UIdirectory             *SDDir;
    G4UIcmdWithABool          *TargetSDCmd;
    G4UIcmdWithABool          *RecoilDetSDCmd;
    G4UIcmdWithABool          *GEMSDCmd;
    G4UIcmdWithABool          *SciPlaneSDCmd;
    G4UIcmdWithABool          *HyCalSDCmd;
    G4UIcmdWithABool          *VirtualSDCmd;
    G4UIdirectory             *CalorimeterDir;
    G4UIcmdWithADoubleAndUnit *AttenuationLGCmd;
    G4UIcmdWithADouble        *ReflectanceLGCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
