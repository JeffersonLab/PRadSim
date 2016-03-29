#ifndef DIGITIZATION_H
#define DIGITIZATION_H

#include "HyCalParameterisation.hh"
#include <string>
#include <vector>
#include <fstream>

#define MAX_HYCAL_BUFFER 2000
#define TRIGGER_THRESHOLD 500 //MeV

class Digitization
{
public:
    struct HyCal_Hit
    {
        int copyNo;
        double energy;
        HyCal_Hit() : copyNo(0), energy(0) {};
        HyCal_Hit(int n, double e) : copyNo(n), energy(e) {};
    };

    struct GEM_Hit
    {
        double x;
        double y;
        double plane_z;
        GEM_Hit() : x(0), y(0), plane_z(0) {};
        GEM_Hit(double hx, double hy, double hz)
        : x(hx), y(hy), plane_z(hz) {};
    };

    enum DAQ_SubSystem
    {
        HyCal,
        GEM,
        MAX_SUBSYSTEM,
    };

public:
    Digitization();
    virtual ~Digitization();
    void EndofEvent();
    void Event(DAQ_SubSystem sub);
    void GEMHits(const double &x, const double &y, const double &z);
    void InitializeHyCalBuffer(uint32_t *buffer);
    void FillBuffer(uint32_t *buffer, const Module_DAQ &module);
    void RegisterModules(HyCalParameterisation *param);
    void UpdateEnergy(const int &copyNo, const double &energy);
    void Clear();
    unsigned short Digitize(const Module_DAQ &module);

private:
    int addEventInfoBank(uint32_t *buffer);
    int addRocData(uint32_t *buffer, int roc_id, int base_index);
    void readLeadGlassMap();
    void readModuleList();
    void initDataFile();

    uint32_t hycal_buffer[MAX_HYCAL_BUFFER];
    uint32_t event_number;
    int data_index[30];
    int event_number_index;
    int hycal_out;
    double totalEnergy;
    unsigned int system_ready;
    unsigned int system_mask;
    std::ofstream gem_out;
    std::vector<Module_DAQ> modules;
    std::vector<GEM_Hit> gem_hits;
};

#endif
