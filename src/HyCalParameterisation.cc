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
// $Id: HyCalParameterisation.cc,2016-03-29 $
// GEANT4 tag $Name: geant-4.10.02-p01 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HyCalParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::HyCalParameterisation(const char* filepath)
{
    std::ifstream module_list(filepath);
    if(module_list.is_open()) {
        std::string line, id;
        double size, x, y;
        int crate, slot, channel, tdc_group;
        double mean, sigma;
        while(std::getline(module_list, line))
        {
            if(line.at(0) == '#')
                continue;

            std::stringstream iss(line);
            iss >> id >> size >> x >> y
                >> crate >> slot >> channel >> tdc_group
                >> mean >> sigma;

            // default setting for lead tungstate
            HyCal_Module_Type t = Lead_Tungstate;
            double z = 0., l = 180.;

            // change if the type is lead glass
            if(id.at(0) == 'G') {
                t = Lead_Glass;
                l = 450.;
                z = 7.;
            }
            Module_DAQ daq(crate, slot, channel, tdc_group, mean, sigma);
            moduleList.push_back(HyCal_Module(id, t, size, l, x, y, z, daq));
        }
        module_list.close();
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::~HyCalParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol)
const
{
    G4ThreeVector origin(moduleList[copyNo].x*mm, moduleList[copyNo].y*mm, moduleList[copyNo].z*mm);
    physVol->SetTranslation(G4ThreeVector(origin));
    physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::ComputeDimensions(G4Box& CalBlock,
                                              const G4int copyNo,
                                              const G4VPhysicalVolume*)
const
{
    CalBlock.SetXHalfLength(moduleList[copyNo].sizeXY/2.*mm);
    CalBlock.SetYHalfLength(moduleList[copyNo].sizeXY/2.*mm);
    CalBlock.SetZHalfLength(moduleList[copyNo].length/2.*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material *HyCalParameterisation::ComputeMaterial(const G4int copyNo,
                                                   G4VPhysicalVolume * /*currVol*/,
                                                   const G4VTouchable * /*parentTouch*/)
{
    switch(moduleList[copyNo].type)
    {
    default: return G4Material::GetMaterial("Galactic");
    case Lead_Glass: return G4Material::GetMaterial("Lead Glass");
    case Lead_Tungstate: return G4Material::GetMaterial("PbWO4");
    }
}

