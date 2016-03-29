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
// $Id: HyCalParameterisation.hh, 2016-03-29$
// GEANT4 tag $Name: geant4.10.02.p01 $
// Developer: Chao Peng
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HyCalParameterisation_H
#define HyCalParameterisation_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"
#include <string>
#include <vector>

class G4VPhysicalVolume;
class G4Box;
class G4Material;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;
class G4Ellipsoid;

enum HyCal_Module_Type
{
    Lead_Glass,
    Lead_Tungstate,
};

struct Module_DAQ
{
    int crate;
    int slot;
    int channel;
    int tdc_group;
    double ped_mean;
    double ped_sigma;
    double energy;
    Module_DAQ() {};
    Module_DAQ(int c, int s, int ch, int t, double m, double sig)
    : crate(c), slot(s), channel(ch), tdc_group(t), ped_mean(m), ped_sigma{sig}, energy(0)
    {};
};

struct HyCal_Module
{
    std::string name;
    HyCal_Module_Type type;
    double sizeXY;
    double length;
    double x;
    double y;
    double z;
    Module_DAQ daq_config;
    HyCal_Module() {};
    HyCal_Module(std::string n, HyCal_Module_Type t,
                 double s, double l, double xx, double yy, double zz,
                 Module_DAQ daq)
    : name(n), type(t), sizeXY(s), length(l), x(xx), y(yy), z(zz), daq_config(daq)
    {};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HyCalParameterisation : public G4VPVParameterisation
{
public:
    HyCalParameterisation(const char* filepath);
    virtual ~HyCalParameterisation();

    void ComputeTransformation(const G4int copyNo,
                               G4VPhysicalVolume* physVol) const;
    void ComputeDimensions(G4Box & CalBlock,
                           const G4int copyNo,
                           const G4VPhysicalVolume* physVol) const;
    G4Material* ComputeMaterial(const G4int copyNo,
                                G4VPhysicalVolume *currVol,
                                const G4VTouchable *pTouch = NULL);
    const std::vector<HyCal_Module> &GetModuleList() {return moduleList;};
    const HyCal_Module &GetModule(size_t i) {return moduleList[i];};
    size_t GetNumber() {return moduleList.size();};

private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,const G4VPhysicalVolume*) const {}

private:
    std::vector<HyCal_Module> moduleList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
