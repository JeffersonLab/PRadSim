#include "Digitization.hh"
#include "Randomize.hh"
#include "evio.h"
#include <iostream>
#include <fstream>
#include <sstream>

Digitization::Digitization()
: event_number(0)
{
    hycal_buffer = new uint32_t[2000];
    InitializeHyCalBuffer(hycal_buffer);

    std::ifstream leadglass_id("config/leadglass_map.txt");
    if(leadglass_id.is_open()) {
        std::string line;
        int copyNo, id;
        while(std::getline(leadglass_id, line))
        {
            if(line.at(0) == '#')
                continue;
            std::stringstream iss(line);
            iss >> copyNo >> id;
            leadglass_map[id] = MAX_LEAD_TUNGSTATE + copyNo - 1;
        }
    }
    leadglass_id.close();

    modules = new daq_info[MAX_MODULE];
    std::ifstream module_list("config/module_list.txt");

    if(module_list.is_open()) {
        std::string line, id;
        int crate, slot, channel, tdc_group, type;
        double mean, sigma;
        while(std::getline(module_list, line))
        {
            if(line.at(0) == '#')
                continue;

            std::stringstream iss(line);
            iss >> id >> crate >> slot >> channel >> tdc_group
                >> mean >> sigma;
            modules[IdToCopyNo(id)] = daq_info(id, crate, slot, channel, tdc_group, mean, sigma);
        }
    }

    module_list.close();

    char outf[] = "output/simrun.dat";
    char mode[] = "w";

    int status = evOpen(outf, mode, &fHandle);
    if(status != S_SUCCESS) {
        std::cerr << "ERROR: cannot open output file \"" << outf << "\"" << std::endl;
    } else {
        std::cout << "output file is \"" << outf << "\"" << std::endl;
    }
}

Digitization::~Digitization()
{
    evClose(fHandle);
    delete modules;
    delete hycal_buffer;
}

void Digitization::Event(double *hycal_energy, std::vector<GEM_Hit> &gem_hits)
{
    ++event_number;
    for(int i = 0; i < MAX_MODULE; ++i)
    {
        FillBuffer(hycal_buffer, modules[i], hycal_energy[i]);
    }

    int status = evWrite(fHandle, hycal_buffer);
    if(status != S_SUCCESS) {
        std::cerr << "ERROR: cannot write event to output file!" << std::endl;
    }
}

int Digitization::IdToCopyNo(const std::string &id)
{
    int index = std::stoi(id.substr(1));

    if(id.at(0) == 'W') {
        if(index > 594) index -= 4;
        else if(index > 560) index -= 2;
        return --index;
    }
    else
        return leadglass_map[index];
}

int Digitization::ReverseCalibration(const double &energy)
{
    return int(energy+0.5);
}

void Digitization::Digitize()
{

}

void Digitization::InitializeHyCalBuffer(uint32_t *buffer)
{
    // event header
    buffer[0] = 0x000007cc;
    buffer[1] = 0x000110cc;

    // event info header
    buffer[2] = 0x00000004;
    buffer[3] = 0xc0000100;

    event_number_index = 4;

    // ROC 6 header
    buffer[7] = 0x00000296;
    buffer[8] = 0x00061001;

    // Fastbus bank
    buffer[9] = 0x00000294;
    buffer[10] = 0x00070100;
    buffer[11] = 0xdc6adc0a;
    for(int i = 0; i < 10; ++i)
    {
        data_index[i] = 13 + i*65;
        buffer[12+i*65] = ((23-2*i)<<27) | 65;
    }
    // End word
    buffer[661] = 0xfabc0005;

    // ROC 5 header
    buffer[669] = 0x00000296;
    buffer[670] = 0x00051001;

    // Fastbus bank
    buffer[671] = 0x00000294;
    buffer[672] = 0x00070100;
    buffer[673] = 0xdc5adc0a;
    for(int i = 0; i < 10; ++i)
    {
        data_index[10+i] = 675 + i*65;
        buffer[674 + i*65] = ((23-2*i)<<27) | 65;
    }

    // End word
    buffer[1323] = 0xfabc0005;

    // ROC 4 header
    buffer[1331] = 0x00000296;
    buffer[1332] = 0x00041001;

    // Fastbus bank
    buffer[1333] = 0x00000294;
    buffer[1334] = 0x00070100;
    buffer[1335] = 0xdc4adc0a;
    for(int i = 0; i < 10; ++i)
    {
        data_index[20+i] = 1337 + i*65;
        buffer[1336 + i*65] = ((22-2*i)<<27) | 65;
    }

    // End word
    buffer[1985] = 0xfabc0005;
}

void Digitization::FillBuffer(uint32_t *buffer, const daq_info &module, const double &energy)
{
    int crate = module.crate;
    int slot = module.slot;
    int channel = module.channel;

    int pos = (6-crate)*10 + ((23-slot)/2);

    int index = data_index[pos] + channel;
    unsigned short val = (unsigned short)(G4RandGauss::shoot(module.pedestal_mean, module.pedestal_sigma) + ReverseCalibration(energy));
    buffer[index] = (slot << 27) | (channel << 17) | val;
}
