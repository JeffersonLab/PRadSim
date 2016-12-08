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
#include "ConfigParser.hh"
#include <fstream>
#include <iostream>
#include <unordered_map>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::HyCalParameterisation(const std::string &mod_file,
                                             const std::string &ped_file,
                                             const std::string &cal_file)
{
    if(!mod_file.empty())
        LoadModuleList(mod_file);
    if(!ped_file.empty())
        LoadPedestal(ped_file);
    if(!cal_file.empty())
        LoadCalibrationFactor(cal_file);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalParameterisation::~HyCalParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::LoadModuleList(const std::string &path)
{
    ConfigParser c_parser;
    if(!c_parser.OpenFile(path)) {
        std::cerr << "ERROR: Missing configuration file \""
                  << path << "\""
                  << std::endl;
        exit(1);
    }

    std::string name, tdc_group;
    double size_x, size_y, x, y, z, l;
    int crate, slot, channel, type;
    HyCal_Module_Type t;

    while(c_parser.ParseLine())
    {
        if(c_parser.NbofElements() < 10)
            continue;

        c_parser >> name  // module name
                 >> crate >> slot >> channel // daq setting
                 >> tdc_group // tdc group name
                 >> type // module type
                 >> size_x >> size_y >> x >> y; // geometry

        if(type == 0) {
            t = Lead_Glass;
            // length of Pb-glass module
            l = 450.; // mm
            // depth of Pb-glass module
            //z = (450. - 180.)/2. - 101.2; // mm, PWO and Pb-glass diff is 101.2 mm
            z = 0.;
        } else if(type == 1) {
            // crystal module
            t = Lead_Tungstate;
            l = 180.; // mm
            //z = 0.; // mm
            z = 101.2 - (450. - 180.)/2.;
        } else {
            continue;
        }

        //std::cout << tdc_group.substr(1) << std::endl;
        int tdc_id = std::stoi(tdc_group.substr(1));
        Module_DAQ daq(crate, slot, channel, tdc_id);
        moduleList.push_back(HyCal_Module(name, t, size_x, size_y, l, x, y, z, daq));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::LoadPedestal(const std::string &path)
{
    ConfigParser c_parser;
    if(!moduleList.size()) {
        std::cout << "WARNING: No module loaded, thus cannot load pedestal data."
                  << std::endl;
        return;
    }

    // this map is for loading the pedestal
    std::unordered_map<int, size_t> daq_map;
    for(size_t i = 0; i < moduleList.size(); ++i)
    {
        int daq_addr = (moduleList.at(i).daq_config.crate << 24) |
                       (moduleList.at(i).daq_config.slot << 16)  |
                       (moduleList.at(i).daq_config.channel);
        auto it = daq_map.find(daq_addr);
        if(it != daq_map.end()) {
            std::cout << "WARNING: DAQ map collision between module "
                      << moduleList.at(i).name
                      << " and "
                      << moduleList.at(it->second).name
                      << std::endl;
        }
        daq_map[daq_addr] = i;
    }

    // read pedestal
    if(!c_parser.OpenFile(path)) {
        std::cerr << "WARNING: Missing pedestal file \""
                  << path << "\""
                  << ", no pedestal data loaded."
                  << std::endl;
        return;
    }

    int crate, slot, channel;
    double sigma, mean;
    while(c_parser.ParseLine())
    {
        if(c_parser.NbofElements() < 5)
            continue;

        c_parser >> crate >> slot >> channel >> mean >> sigma;

        int daq_addr = (crate << 24) | (slot << 16) | channel;
        auto it = daq_map.find(daq_addr);
        if(it == daq_map.end())
            continue;
        moduleList[it->second].daq_config.ped_mean = mean;
        moduleList[it->second].daq_config.ped_sigma = sigma;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::LoadCalibrationFactor(const std::string &path)
{
    ConfigParser c_parser;
    if(!moduleList.size()) {
        std::cout << "WARNING: No module loaded, thus cannot load calibration data."
                  << std::endl;
        return;
    }

    // this map is for loading the pedestal
    std::unordered_map<std::string, size_t> name_map;
    for(size_t i = 0; i < moduleList.size(); ++i)
    {
        auto it = name_map.find(moduleList.at(i).name);
        if(it != name_map.end()) {
            std::cout << "WARNING: Name map collision between module "
                      << moduleList.at(i).name
                      << " and "
                      << moduleList.at(it->second).name
                      << std::endl;
        }
        name_map[moduleList.at(i).name] = i;
    }

    // read pedestal
    if(!c_parser.OpenFile(path)) {
        std::cerr << "WARNING: Missing calibration file \""
                  << path << "\""
                  << ", no gain factor loaded."
                  << std::endl;
        return;
    }

    std::string name;
    double factor;
    while(c_parser.ParseLine())
    {
        if(c_parser.NbofElements() < 2)
            continue;

        c_parser >> name >> factor;

        auto it = name_map.find(name);
        if(it == name_map.end())
            continue;
        moduleList[it->second].daq_config.gain_factor = factor;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol)
const
{
    if((size_t)copyNo >= moduleList.size())
    {
        std::cerr << "HyCalParameterisation: trying to load module no." << copyNo
                  << ", but the loaded module list only has " << moduleList.size()
                  << " modules, please check the setup."
                  << std::endl;
        exit(1);
    }

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
    if((size_t)copyNo >= moduleList.size())
    {
        std::cerr << "HyCalParameterisation: trying to load module no." << copyNo
                  << ", but the loaded module list only has " << moduleList.size()
                  << " modules, please check the setup."
                  << std::endl;
        exit(1);
    }

    CalBlock.SetXHalfLength(moduleList[copyNo].sizeX/2.*mm);
    CalBlock.SetYHalfLength(moduleList[copyNo].sizeY/2.*mm);
    CalBlock.SetZHalfLength(moduleList[copyNo].length/2.*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material *HyCalParameterisation::ComputeMaterial(const G4int copyNo,
                                                   G4VPhysicalVolume * /*currVol*/,
                                                   const G4VTouchable * /*parentTouch*/)
{
    if((size_t)copyNo >= moduleList.size())
    {
        std::cerr << "HyCalParameterisation: trying to load module no." << copyNo
                  << ", but the loaded module list only has " << moduleList.size()
                  << " modules, please check the setup."
                  << std::endl;
        exit(1);
    }

    switch(moduleList[copyNo].type)
    {
    default: return G4Material::GetMaterial("Galaxy");
    case Lead_Glass: return G4Material::GetMaterial("Lead Glass");
    case Lead_Tungstate: return G4Material::GetMaterial("PbWO4");
    }
}

