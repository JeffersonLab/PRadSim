#include "Digitization.hh"
#include "Randomize.hh"
#include "evio.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "sys/time.h"

Digitization::Digitization()
: event_number(0), hycal_out(-1)
{
    gem_out.open("output/gem_pos.dat");
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

    std::ifstream module_list("config/module_list.txt");

    if(module_list.is_open()) {
        std::string line, id;
        int crate, slot, channel, tdc_group;
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

    char outf[] = "output/simrun.evio";
    char mode[] = "w";

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

Digitization::~Digitization()
{
    uint32_t now = time(NULL);
    uint32_t end[5] = {0x00000004, 0x002001cc, now, event_number, 0x00000000};
    evWrite(hycal_out, end);
    evClose(hycal_out);
}

void Digitization::Event(double *hycal_energy, std::vector<GEM_Hit> &gem_hits)
{
    hycal_buffer[event_number_index] = ++event_number;

    for(int i = 0; i < MAX_MODULE; ++i)
    {
        FillBuffer(hycal_buffer, modules[i], hycal_energy[i]);
    }

    if(evWrite(hycal_out, hycal_buffer) != S_SUCCESS) {
        std::cerr << "ERROR: cannot write event to output file!" << std::endl;
    }

/* text file format
    for(unsigned int i = 0; i <= hycal_buffer[0]; ++i)
    {
        out << "0x" << std::hex << std::setw(8) << std::setfill('0') << hycal_buffer[i] << std::endl;
    }
*/
    for(auto &gem_hit : gem_hits)
    {
        gem_out << std::setw(12) << event_number << "  "
                << std::setw(12) << gem_hit.x
                << std::setw(12) << gem_hit.y
                << std::setw(12) << gem_hit.plane_z
                << std::endl;
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

unsigned short Digitization::Digitize(const daq_info &module, const double &energy)
{
    double ped = G4RandGauss::shoot(module.pedestal_mean, module.pedestal_sigma);

    return ped + energy*3.;
}

void Digitization::InitializeHyCalBuffer(uint32_t *buffer)
{
    int index = 0;

    // event header;
    buffer[index++] = 0x00000000;
    buffer[index++] = 0x000110cc;
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
    buffer[index++] = 0x00070100;
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

void Digitization::FillBuffer(uint32_t *buffer, const daq_info &module, const double &energy)
{
    int crate = module.crate;
    int slot = module.slot;
    int channel = module.channel;

    int pos = (6-crate)*10 + ((23-slot)/2);

    int index = data_index[pos] + channel;
    unsigned short val = Digitize(module, energy);
    buffer[index] = (slot << 27) | (channel << 17) | val;
}
