//
// HyCalDigitization.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HyCalDigitization_h
#define HyCalDigitization_h 1

#include "ConfigObject.h"
#include "StandardDigiBase.hh"

#include "evio.h"

#include <string>
#include <vector>

#define TRIGGER_THRESHOLD 500 //MeV

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HyCalDigitization : public StandardDigiBase, public ConfigObject
{
public:
    enum HyCal_Module_Type {
        Lead_Glass,
        Lead_Tungstate,
    };

    struct Module_DAQ {
        int crate;
        int slot;
        int channel;
        int tdc_group;
        double ped_mean;
        double ped_sigma;
        double energy;
        double gain_factor;
        Module_DAQ() {};
        Module_DAQ(int c, int s, int ch, int t) : crate(c), slot(s), channel(ch), tdc_group(t), ped_mean(0), ped_sigma{0}, energy(0), gain_factor(0) {};
    };

    struct HyCal_Module {
        std::string name;
        int id;
        HyCal_Module_Type type;
        double sizex;
        double sizey;
        double length;
        double x;
        double y;
        double z;
        Module_DAQ daq_config;
        HyCal_Module() {};
        HyCal_Module(std::string n, int ii, HyCal_Module_Type t, double sx, double sy, double l, double xx, double yy, double zz, Module_DAQ daq) : name(n), id(ii), type(t), sizex(sx), sizey(sy), length(l), x(xx), y(yy), z(zz), daq_config(daq) {};
    };

    struct CalPeriod {
        int begin;
        int end;
        int main;
        int sub;
        CalPeriod() : begin(0), end(0), main(0), sub(0) {};
        CalPeriod(int b, int e, int m, int s) : begin(b), end(e), main(m), sub(s) {};
    };

public:
    HyCalDigitization(const char *name, const char *path);
    virtual ~HyCalDigitization();

    void Configure(const std::string &path = "");

    int PreStart(uint32_t *buffer, int base_index);
    bool ProcessEvent(uint32_t *buffer);

    void Clear();

private:
    void LoadModuleList(const std::string &path);
    void LoadChannelList(const std::string &path);
    void LoadPedestal(const std::string &path);
    void ReadCalPeriodFile(const std::string &path);
    void ChooseRun(int run);
    void LoadCalibrationFile(const std::string &path);

    int addRocData(uint32_t *buffer, int roc_id, int base_index);
    void FillBuffer(uint32_t *buffer, const HyCal_Module &module);
    unsigned short Digitize(const HyCal_Module &module);

    int data_index[30];

    double fTotalE;

    std::vector<CalPeriod> cal_period;
    std::vector<HyCal_Module> moduleList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
