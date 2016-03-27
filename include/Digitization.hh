#ifndef DIGITIZATION_H
#define DIGITIZATION_H

#include <unordered_map>
#include <string>
#include <vector>

#define MAX_LEAD_TUNGSTATE 1152
#define MAX_LEAD_GLASS 576
#define MAX_MODULE MAX_LEAD_TUNGSTATE+MAX_LEAD_GLASS

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

    struct daq_info
    {
        std::string name;
        int crate;
        int slot;
        int channel;
        int tdc_group;
        double pedestal_mean;
        double pedestal_sigma;
        daq_info() {};
        daq_info(std::string n, int c, int s, int ch, int tdc, double mean, double sigma)
        : name(n), crate(c), slot(s), channel(ch), tdc_group(tdc),
          pedestal_mean(mean), pedestal_sigma(sigma) {};
    };

public:
    Digitization();
    virtual ~Digitization();
    void Event(double *hycal_energy, std::vector<GEM_Hit> &gem_hits);
    void Digitize();
    void InitializeHyCalBuffer(uint32_t *buffer);
    void FillBuffer(uint32_t *buffer, const daq_info &module, const double &energy);
    int IdToCopyNo(const std::string &id);
    int ReverseCalibration(const double &energy);

private:
    daq_info *modules;
    uint32_t *hycal_buffer;
    uint32_t event_number;
    std::unordered_map<int, int> leadglass_map;
    int fHandle;
    int data_index[30];
    int event_number_index;
};

#endif
