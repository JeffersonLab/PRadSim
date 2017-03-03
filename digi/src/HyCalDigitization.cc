//
// HyCalDigitization.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HyCalDigitization.hh"

#include "ConfigObject.h"
#include "ConfigParser.h"
#include "ConfigValue.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TRandom2.h"

#include <ctime>
#include <unordered_map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static TRandom2 *RandGen = new TRandom2();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalDigitization::HyCalDigitization(const char *name, const char *path) : StandardDigiBase(name), fTotalE(0)
{
    RandGen->SetSeed((UInt_t)time(NULL));

    Configure(std::string(path));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalDigitization::~HyCalDigitization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::Configure(const std::string &path)
{
    ConfigObject::Configure(path);

    LoadModuleList(GetConfig<std::string>("Module List"));
    LoadChannelList(GetConfig<std::string>("DAQ Channel List"));
    LoadPedestal(GetConfig<std::string>("Pedestal File"));
    std::string file_path = ConfigParser::form_path(GetConfig<std::string>("Calibration Folder"), GetConfig<std::string>("Calibration Period File"));
    ReadCalPeriodFile(file_path);;
    ChooseRun(GetConfig<int>("Run Number"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::PreStart(uint32_t *buffer, int base_index)
{
    int index = 0;

    for (int roc_id = 6; roc_id >= 4; --roc_id)
        index += addRocData(&buffer[index + base_index], roc_id, index + base_index);

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool HyCalDigitization::ProcessEvent(uint32_t *buffer)
{
    for (int i = 0; i < fData.N; i++) {
        fTotalE += fData.Edep[i];
        moduleList[fData.CopyNo[i]].daq_config.energy += fData.Edep[i];
    }

    if (fTotalE < TRIGGER_THRESHOLD) {
        Clear();
        return false;
    }

    for (auto &module : moduleList) {
        if (module.daq_config.energy > 0)
            std::cout << module.name << " " << module.id << " " << module.daq_config.energy << std::endl;

        FillBuffer(buffer, module);
    }

    Clear();
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::Clear()
{
    for (auto &module : moduleList)
        module.daq_config.energy = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::LoadModuleList(const std::string &path)
{
    if (path.empty())
        return;

    ConfigParser c_parser;

    if (!c_parser.ReadFile(path)) {
        std::cerr << "Error: Failed to read module list file " << "\"" << path << "\"." << std::endl;
        return;
    }

    std::string name, type;
    double sizex, sizey, length, x, y, z;
    HyCal_Module_Type t;
    int id = 0;

    // some info that is not read from list
    while (c_parser.ParseLine()) {
        if (!c_parser.CheckElements(8))
            continue;

        c_parser >> name >> type >> sizex >> sizey >> length >> x >> y >> z;

        if (type.compare("PbGlass") == 0) {
            t = Lead_Glass;
            z = 0.;
        } else if (type.compare("PbWO4") == 0) {
            t = Lead_Tungstate;
            z = 101.2 - (450. - 180.) / 2.;
        } else
            continue;

        Module_DAQ daq(0, 0, 0, 0);
        moduleList.push_back(HyCal_Module(name, id++, t, sizex, sizey, length, x, y, z, daq));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::LoadChannelList(const std::string &path)
{
    if (path.empty())
        return;

    ConfigParser c_parser;
    // set special splitter
    c_parser.SetSplitters(",: \t");

    if (!c_parser.ReadFile(path)) {
        std::cerr << "Error: Failed to read channel list file " << "\"" << path << "\"." << std::endl;
        return;
    }

    // we accept 2 types of channels
    // tdc, adc
    std::vector<std::string> types = {"TDC", "ADC"};
    // tdc args: name crate slot channel
    // adc args: name crate slot channel tdc
    std::vector<int> expect_args = {4, 4};
    std::vector<int> option_args = {0, 1};

    // this vector is to store all the following arguments
    std::vector<std::vector<std::vector<ConfigValue>>> ch_args(types.size());

    // read all the elements in
    while (c_parser.ParseLine()) {
        std::string type = c_parser.TakeFirst();
        size_t i = 0;

        for (; i < types.size(); ++i) {
            if (ConfigParser::strcmp_case_insensitive(type, types.at(i))) {
                // only save elements from expected format
                if (c_parser.CheckElements(expect_args.at(i), option_args.at(i)))
                    ch_args[i].push_back(c_parser.TakeAll<std::vector>());

                break;
            }
        }

        if (i >= types.size())  // did not find any type
            std::cout << "Warning: Undefined channel type " << type << " in channel list file " << "\"" << path << "\"" << std::endl;
    }

    // create ADC channels, and add them to TDC groups
    for (auto &args : ch_args[1]) {
        std::string name(args[0]);

        for (auto &module : moduleList) {
            if (module.name == name) {
                module.daq_config.crate = args[1].Int();
                module.daq_config.slot = args[2].Int();
                module.daq_config.channel = args[3].Int();
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::LoadPedestal(const std::string &path)
{
    if (path.empty())
        return;

    ConfigParser c_parser;

    // read pedestal
    if (!c_parser.OpenFile(path)) {
        std::cerr << "WARNING: Missing pedestal file \"" << path << "\"" << ", no pedestal data loaded." << std::endl;
        return;
    }

    if (!moduleList.size()) {
        std::cout << "WARNING: No module loaded, thus cannot load pedestal data." << std::endl;
        return;
    }

    // this map is for loading the pedestal
    std::unordered_map<int, size_t> daq_map;

    for (size_t i = 0; i < moduleList.size(); ++i) {
        int daq_addr = (moduleList.at(i).daq_config.crate << 24) | (moduleList.at(i).daq_config.slot << 16) | (moduleList.at(i).daq_config.channel);
        auto it = daq_map.find(daq_addr);

        if (it != daq_map.end())
            std::cout << "WARNING: DAQ map collision between module " << moduleList.at(i).name << " and " << moduleList.at(it->second).name << std::endl;

        daq_map[daq_addr] = i;
    }

    int crate, slot, channel;
    double sigma, mean;

    while (c_parser.ParseLine()) {
        if (!c_parser.CheckElements(5))
            continue;

        c_parser >> crate >> slot >> channel >> mean >> sigma;

        int daq_addr = (crate << 24) | (slot << 16) | channel;
        auto it = daq_map.find(daq_addr);

        if (it == daq_map.end())
            continue;

        moduleList[it->second].daq_config.ped_mean = mean;
        moduleList[it->second].daq_config.ped_sigma = sigma;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::ReadCalPeriodFile(const std::string &path)
{
    if (path.empty())
        return;

    ConfigParser c_parser;

    if (!c_parser.ReadFile(path)) {
        std::cerr << "Error: Failed to read calibration period file " << "\"" << path << "\"" << std::endl;
        return;
    }

    cal_period.clear();

    int period, sub_period, begin, end;

    while (c_parser.ParseLine()) {
        if (!c_parser.CheckElements(4))
            continue;

        c_parser >> period >> sub_period >> begin >> end;

        // no need to update
        cal_period.emplace_back(begin, end, period, sub_period);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::ChooseRun(int run)
{
    int period = 5, subperiod = 1;

    for (auto &p : cal_period) {
        if (run <= p.end && run >= p.begin) {
            period = p.main;
            subperiod = p.sub;
        }
    }

    SetConfigValue("Period", period);
    SetConfigValue("Sub-period", subperiod);

    // run info file
    std::string file_path;

    // calibration file
    file_path = ConfigParser::form_path(GetConfig<std::string>("Calibration Folder"), GetConfig<std::string>("Calibration File"));
    LoadCalibrationFile(file_path);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::LoadCalibrationFile(const std::string &path)
{
    if (path.empty())
        return;

    ConfigParser c_parser;

    if (!c_parser.ReadFile(path)) {
        std::cerr << "Error: Failed to read calibration file " << " \"" << path << "\"" << std::endl;
        return;
    }

    // this map is for loading the pedestal
    std::unordered_map<std::string, size_t> name_map;

    for (size_t i = 0; i < moduleList.size(); ++i) {
        auto it = name_map.find(moduleList.at(i).name);

        if (it != name_map.end())
            std::cout << "WARNING: Name map collision between module " << moduleList.at(i).name << " and " << moduleList.at(it->second).name << std::endl;

        name_map[moduleList.at(i).name] = i;
    }

    std::string name;
    double factor;

    while (c_parser.ParseLine()) {
        if (!c_parser.CheckElements(4, -1)) // more than 4 elements
            continue;

        c_parser >> name >> factor;

        auto it = name_map.find(name);

        if (it == name_map.end())
            continue;

        moduleList[it->second].daq_config.gain_factor = factor;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::addRocData(uint32_t *buffer, int roc_id, int base_index)
{
    int index = 0;
    int nslot, slot[25];

    switch (roc_id) {
    default:
        return 0;

    case 4:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 22 - 2 * i;

        break;

    case 5:
    case 6:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 23 - 2 * i;

        break;
    }

    // add roc header
    buffer[index++] = 0x00000000;
    buffer[index++] = (roc_id << 16) | 0x1001;

    // add TI bank 11 words
    buffer[index++] = 0x0000000a; // 10 + 1 words in total
    buffer[index++] = 0xe10a0100; // TI bank header
    buffer[index + 2] = 2 << 24; // only 2nd word matters, it defines trigger type, here it is total sum
    index += 9; // TI bank expects 9 words

    buffer[index++] = 0x00000000;
    buffer[index++] = 0xe1200100; // Fastbus bank header
    // roc id and board number
    buffer[index++] = 0xdc0adc00 | ((roc_id & 0xf) << 20) | (nslot & 0xff);

    for (int i = 0; i < nslot; ++i) {
        buffer[index++] = (slot[i] << 27) | 65;
        data_index[(6 - roc_id) * 10 + i] = index + base_index;

        for (int ch = 0; ch < 64; ++ch)
            buffer[index++] = (slot[i] << 27) | (ch << 17);
    }

    buffer[index++] = 0xfabc0005; // end word
    buffer[0] = index - 1; // roc bank size
    buffer[13] = index - 14; // fastbus bank size

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::FillBuffer(uint32_t *buffer, const HyCal_Module &module)
{
    int crate = module.daq_config.crate;
    int slot = module.daq_config.slot;
    int channel = module.daq_config.channel;

    int pos = (6 - crate) * 10 + ((23 - slot) / 2);

    int index = data_index[pos] + channel;
    unsigned short val = Digitize(module);
    std::cout << crate << " " << slot << " " << channel << " " << val << " " << module.daq_config.ped_mean + module.daq_config.energy / module.daq_config.gain_factor << std::endl;
    buffer[index] = (slot << 27) | (channel << 17) | val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned short HyCalDigitization::Digitize(const HyCal_Module &module)
{
    //double ped = RandGen->Gaus(module.daq_config.ped_mean, module.daq_config.ped_sigma);
    double ped = module.daq_config.ped_mean;
    return ped + (module.daq_config.energy) / module.daq_config.gain_factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
