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
// EventMessenger.cc
// Developer : Chao Peng
// History:
//   Aug 2012, C. Peng, Original version.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventMessenger.hh"

#include "EventAction.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventMessenger::EventMessenger(EventAction *act) : G4UImessenger(), Action(act)
{
    eventDir = new G4UIdirectory("/pradsim/event/");
    eventDir->SetGuidance("event control");

    PrintCmd = new G4UIcmdWithAnInteger("/pradsim/event/printmodulo", this);
    PrintCmd->SetGuidance("Print events modulo n");
    PrintCmd->SetParameterName("EventNb", false);
    PrintCmd->SetRange("EventNb>0");

    OnlyRecordHitsCmd = new G4UIcmdWithABool("/pradsim/event/onlyrecordhits", this);
    OnlyRecordHitsCmd->SetGuidance("Only write the rootfile if there is a hit on HyCal");
    OnlyRecordHitsCmd->SetParameterName("onlyrecordhits", false);
    OnlyRecordHitsCmd->SetDefaultValue(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventMessenger::~EventMessenger()
{
    delete PrintCmd;
    delete OnlyRecordHitsCmd;
    delete eventDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == PrintCmd)
        Action->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));

    if (command == OnlyRecordHitsCmd)
        Action->SetOnlyRecordHits(OnlyRecordHitsCmd->GetNewBoolValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
