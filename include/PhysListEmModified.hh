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
// PhysListEmModified.hh
// Developer : Geant4
// Modified by Chao Gu from G4EmStandardPhysics_option4.hh
// History:
//   May 2018, C. Gu, For Brems test.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysListEmModified_h
#define PhysListEmModified_h 1

#include "G4VPhysicsConstructor.hh"

#include "G4String.hh"

class PhysicsListMessenger;

class G4VEmModel;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysListEmModified : public G4VPhysicsConstructor
{
public:
    PhysListEmModified(G4int ver = 1, const G4String &name = "");
    virtual ~PhysListEmModified();

    void ConstructParticle();
    void ConstructProcess();

    inline void SetBremsstrahlungAngularGenerator(G4String val);

private:
    void SetBremsstrahlungAngularGenerator(G4VEmModel *, const G4String &);

    G4String fBremsAngularGeneratorType;

    G4int verbose;

    PhysicsListMessenger *physlistMessenger;
};

inline void PhysListEmModified::SetBremsstrahlungAngularGenerator(G4String val)
{
    fBremsAngularGeneratorType = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








