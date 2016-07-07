#include "Digitization.hh"
#include "Randomize.hh"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "sys/time.h"

Digitization::Digitization()
: event_number(0), hycal_out(-1), totalEnergy(0), system_ready(0), system_mask(0)
{
    gem_out.open("output/gem_pos.dat");
    initDataFile();
    for(int i = 0; i < int(MAX_SUBSYSTEM); ++i)
    {
        system_mask |= 1 << i;
    }
}

Digitization::~Digitization()
{
    uint32_t now = time(NULL);
    uint32_t end[5] = {0x00000004, 0x002001cc, now, event_number, 0x00000000};
    evWrite(hycal_out, end);
    evClose(hycal_out);
}

void Digitization::initDataFile()
{
    InitializeHyCalBuffer(hycal_buffer);

    std::ifstream data_file("config/data_file.txt");
    std::string file_name;
    int run_number = 1;

    if(data_file.is_open()) {
        std::string line;
        while(std::getline(data_file, line))
        {
            if(line.at(0) == '#')
                continue;

            std::stringstream iss(line);
            iss >> file_name >> run_number;
            break;
        }
        data_file.close();
    }

    if(file_name.empty()) file_name = "simrun";

    std::ofstream data_file_o("config/data_file.txt");
    data_file_o << file_name << "  " << run_number+1;
    data_file_o.close();

    std::string path = "output/" + file_name + "_" + std::to_string(run_number) + ".evio";
    char mode[] = "w";
    char outf[256];

    strcpy(outf, path.c_str());

    int status = evOpen(outf, mode, &hycal_out);
    if(status != S_SUCCESS) {
        std::cerr << "ERROR CODE "
                  << "0x" << std::hex << std::setw(8) << std::setfill('0') << status
                  << ": cannot open output file \"" << outf << "\"" << std::endl;
    }

    uint32_t now = time(NULL);
    uint32_t prestart[5] = {0x00000004, 0x001101cc, now, 0x00000260, 0x00000000};
    uint32_t go[5] = {0x00000004, 0x001201cc, now+1, 0x00000000, 0x00000000};

    evWrite(hycal_out, prestart);
    evWrite(hycal_out, go);
}

void Digitization::RegisterModules(HyCalParameterisation *param)
{
    std::vector<HyCal_Module> moduleList = param->GetModuleList();
    modules.clear();
    for(auto &module : moduleList)
    {
        modules.push_back(module.daq_config);
    }
}

void Digitization::UpdateEnergy(const G4int &copyNo, const double &edep)
{
    modules[copyNo].energy += edep;
    totalEnergy += edep;
}

void Digitization::GEMHits(const double &x, const double &y, const double &z)
{
    gem_hits.push_back(GEM_Hit(x, y, z));
}

void Digitization::Clear()
{
    for(auto &module : modules)
    {
        module.energy = 0;
    }
    gem_hits.clear();

    totalEnergy = 0;
    system_ready = 0;
}

void Digitization::Event(DAQ_SubSystem sub)
{
    system_ready |= 1<<int(sub);

    if((system_ready & system_mask) == system_mask)
        EndofEvent();
}

void Digitization::EndofEvent()
{
    if(totalEnergy < TRIGGER_THRESHOLD) {
        Clear();
        return;
    }

    hycal_buffer[event_number_index] = ++event_number;

    for(const auto &module : modules)
    {
        FillBuffer(hycal_buffer, module);
    }

    if(evWrite(hycal_out, hycal_buffer) != S_SUCCESS) {
        std::cerr << "ERROR: cannot write event to output file!" << std::endl;
    }


    for(auto &gem_hit : gem_hits)
    {
        gem_out << std::setw(12) << event_number << "  "
                << std::setw(12) << gem_hit.x
                << std::setw(12) << gem_hit.y
                << std::setw(12) << gem_hit.plane_z
                << std::endl;
    }

    Clear();
    return;
}

unsigned short Digitization::Digitize(const Module_DAQ &module)
{
    double ped = G4RandGauss::shoot(module.ped_mean, module.ped_sigma);

    return ped + module.energy*3.;
}

void Digitization::InitializeHyCalBuffer(uint32_t *buffer)
{
    int index = 0;

    // event header;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x008110cc;
    event_number_index = index + 2;
    index += addEventInfoBank(&buffer[index]);
    for(int roc_id = 6; roc_id >= 4; --roc_id)
    {
        index += addRocData(&buffer[index], roc_id, index);
    }

    buffer[0] = index - 1;
}

int Digitization::addEventInfoBank(uint32_t *buffer)
{
    int index = 0;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0xc0000100;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x00000000;


    buffer[0] = index - 1;
    return index;
}

int Digitization::addRocData(uint32_t *buffer, int roc_id, int base_index)
{
    int index = 0;
    int nslot, slot[25];

    switch(roc_id)
    {
    default: return 0;
    case 4:
        nslot = 10;
        for(int i = 0; i < nslot; ++i)
            slot[i] = 22 - 2*i;
        break;
    case 5:
    case 6:
        nslot = 10;
        for(int i = 0; i < nslot; ++i)
            slot[i] = 23 - 2*i;
        break;
    }

    buffer[index++] = 0x00000000;
    buffer[index++] = (roc_id << 16) | 0x1001;

    buffer[index++] = 0x00000000;
    buffer[index++] = 0xe1200100;
    buffer[index++] = 0xdc0adc00 | ((roc_id&0xf) << 20) | (nslot & 0xff);

    for(int i = 0; i < nslot; ++i)
    {
        buffer[index++] = (slot[i] << 27) | 65;
        data_index[(6-roc_id)*10+i] = index + base_index;
        for(int ch = 0; ch < 64; ++ch)
            buffer[index++] = (slot[i] << 27) | (ch << 17);
    }

    buffer[index++] = 0xfabb0005;
    buffer[0] = index - 1;
    buffer[2] = index - 3;
    return index;
}

void Digitization::FillBuffer(uint32_t *buffer, const Module_DAQ &module)
{
    int crate = module.crate;
    int slot = module.slot;
    int channel = module.channel;

    int pos = (6-crate)*10 + ((23-slot)/2);

    int index = data_index[pos] + channel;
    unsigned short val = Digitize(module);
    buffer[index] = (slot << 27) | (channel << 17) | val;
}
