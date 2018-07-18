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
// PhysicsListPureEm.cc
// Developer : Chao Gu
// History:
//   Jul 2018, C. Gu, Original version.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListPureEm.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListPureEm::PhysListPureEm(G4String type, G4int ver) : G4VModularPhysicsList()
{
    this->defaultCutValue = 0.7 * CLHEP::mm;
    this->SetVerboseLevel(ver);

    if (type == "EM")
        this->RegisterPhysics(new G4EmStandardPhysics(ver));
    else if (type == "EM_EMV")
        this->RegisterPhysics(new G4EmStandardPhysics_option1(ver));
    else if (type == "EM_EMX")
        this->RegisterPhysics(new G4EmStandardPhysics_option2(ver));
    else if (type == "EM_EMY")
        this->RegisterPhysics(new G4EmStandardPhysics_option3(ver));
    else if (type == "EM_EMZ")
        this->RegisterPhysics(new G4EmStandardPhysics_option4(ver));
    else if (type == "EM_LIV")
        this->RegisterPhysics(new G4EmLivermorePhysics(ver));
    else if (type == "EM_PEN")
        this->RegisterPhysics(new G4EmPenelopePhysics(ver));
    else if (type == "EM__GS")
        this->RegisterPhysics(new G4EmStandardPhysicsGS(ver));
    else if (type == "EM__SS")
        this->RegisterPhysics(new G4EmStandardPhysicsSS(ver));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListPureEm::~PhysListPureEm()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListPureEm::SetCuts()
{
    this->SetCutsWithDefault();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
