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
// HyCalParameterisation.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HyCalParameterisation.hh"

#include "G4Box.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <string>
#include <unordered_map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::HyCalParameterisation(const std::string &path, G4int conf) : G4VPVParameterisation(), ConfigObject(), fConfig(conf)
{
    if (!path.empty())
        Configure(path);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::~HyCalParameterisation()
{
    //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  HyCalParameterisation::Configure(const std::string &path)
{
    ConfigObject::Configure(path);

    ReadModuleList(GetConfig<std::string>("Module List"));
    SortModuleList();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
    if ((size_t)copyNo >= module_list.size()) {
        G4cout << "ERROR: trying to load module no." << copyNo << ", but the loaded module list only has " << module_list.size() << " modules" << G4endl;
        exit(1);
    }

    auto module = module_list[copyNo];

    G4ThreeVector origin(module->GetX(), module->GetY(), module->GetZ());
    physVol->SetTranslation(G4ThreeVector(origin));
    physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::ComputeDimensions(G4Box &CalBlock, const G4int copyNo, const G4VPhysicalVolume *) const
{
    if ((size_t)copyNo >= module_list.size()) {
        G4cout << "ERROR: trying to load module no." << copyNo << ", but the loaded module list only has " << module_list.size() << " modules" << G4endl;
        exit(1);
    }

    auto module = module_list[copyNo];

    CalBlock.SetXHalfLength(module->GetSizeX() / 2.0);
    CalBlock.SetYHalfLength(module->GetSizeY() / 2.0);
    CalBlock.SetZHalfLength(module->GetSizeZ() / 2.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material *HyCalParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * /*currVol*/, const G4VTouchable * /*parentTouch*/)
{
    if ((size_t)copyNo >= module_list.size()) {
        G4cout << "ERROR: trying to load module no." << copyNo << ", but the loaded module list only has " << module_list.size() << " modules" << G4endl;
        exit(1);
    }

    switch (module_list[copyNo]->GetType()) {
    default:
        return G4Material::GetMaterial("Galaxy");

    case PRadHyCalModule::PbGlass:
        return G4Material::GetMaterial("PbGlass");

    case PRadHyCalModule::PbWO4:
        return G4Material::GetMaterial("PbWO4");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
