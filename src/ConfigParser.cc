//============================================================================//
// A simple parser class to read text file                                    //
//                                                                            //
// Chao Peng                                                                  //
// 06/07/2016                                                                 //
//============================================================================//

#include "ConfigParser.hh"
#include <cstring>
#include <climits>
#include <algorithm>

using namespace std;

// Config Value
ostream &operator << (ostream &os, ConfigValue &b)
{
    return  os << b._value;
}


ConfigValue::ConfigValue(const string &value)
: _value(value)
{}

ConfigValue::ConfigValue(const int &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const float &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const double &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long double &value)
: _value(to_string(value))
{}

char ConfigValue::Char()
{
    try {
       int value = stoi(_value);
       if(value > CHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to char." << endl;
       return (char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to char. 0 returned." << endl;
             return 0;
    }
}

unsigned char ConfigValue::UChar()
{
    try {
       unsigned long value = stoul(_value);
       if(value > UCHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned char." << endl;
       return (unsigned char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned char. 0 returned." << endl;
             return 0;
    }
}

short ConfigValue::Short()
{
    try {
       int value = stoi(_value);
       if(value > SHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to short." << endl;
       return (short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to short. 0 returned." << endl;
             return 0;
    }
}

unsigned short ConfigValue::UShort()
{
    try {
       unsigned long value = stoul(_value);
       if(value > USHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned short." << endl;
       return (unsigned short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned short. 0 returned." << endl;
             return 0;
    }
}

int ConfigValue::Int()
{
    try {
       return stoi(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to int. 0 returned." << endl;
             return 0;
    }
}

unsigned int ConfigValue::UInt()
{
    try {
        return (unsigned int)stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned int. 0 returned." << endl;
             return 0;
    }
}

long ConfigValue::Long()
{
    try {
        return stol(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long. 0 returned." << endl;
             return 0;
    }
}

long long ConfigValue::LongLong()
{
    try {
        return stoll(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long long. 0 returned." << endl;
             return 0;
    }
}

unsigned long ConfigValue::ULong()
{
    try {
        return stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long. 0 returned." << endl;
             return 0;
    }
}

unsigned long long ConfigValue::ULongLong()
{
    try {
        return stoull(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long long. 0 returned." << endl;
             return 0;
    }
}

float ConfigValue::Float()
{
    try {
        return stof(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to float. 0 returned." << endl;
             return 0;
    }
}

double ConfigValue::Double()
{
    try {
        return stod(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to double. 0 returned." << endl;
             return 0;
    }
}

long double ConfigValue::LongDouble()
{
    try {
        return stold(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long double. 0 returned." << endl;
             return 0;
    }
}

const char *ConfigValue::c_str()
{
    return _value.c_str();
}

// Config Parser
ConfigParser::ConfigParser(const string &s,
                           const string &w,
                           const vector<string> &c)
: splitters(s), white_space(w), comment_marks(c)
{
}

ConfigParser::~ConfigParser()
{
    CloseFile();
}

void ConfigParser::AddCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it == comment_marks.end())
        comment_marks.push_back(c);
}

void ConfigParser::RemoveCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it != comment_marks.end())
        comment_marks.erase(it);
}

void ConfigParser::EraseCommentMarks()
{
    comment_marks.clear();
}

bool ConfigParser::OpenFile(const string &path)
{
    infile.open(path);
    return infile.is_open();
}

void ConfigParser::CloseFile()
{
    infile.close();
}

void ConfigParser::OpenBuffer(char *buf)
{
    string buffer = buf;

    string line;
    for(auto c : buffer)
    {
        if(c != '\n') {
            line += c;
        } else {
            lines.push(line);
            line = "";
        }
    }
}

void ConfigParser::ClearBuffer()
{
    queue<string>().swap(lines);
}

string ConfigParser::GetLine()
{
    if(lines.size()) {
        string out = lines.front();
        lines.pop();
        return out;
    }

    return "";
}

bool ConfigParser::ParseLine()
{
    queue<string>().swap(elements);

    if(infile.is_open()) {
        string line;
        while(elements.empty())
        {
            if(!getline(infile, line))
                return false; // end of file
            ParseLine(line);
        }
    } else {
       while(elements.empty())
        {
            if(lines.empty())
                return false; // end of buffer
            ParseLine(lines.front());
            lines.pop();
        }
    }
    return true; // parsed a line
}

void ConfigParser::ParseLine(const string &line)
{
    string trim_line = trim(comment_out(line), white_space);
    queue<string> eles = split(trim_line, splitters);

    while(eles.size())
    {
        string ele = trim(eles.front(), white_space);
        if(ele.size())
            elements.push(ele);
        eles.pop();
    }
}

ConfigValue ConfigParser::TakeFirst()
{
    if(elements.empty()) {
        cout << "Config Parser: WARNING, trying to take elements while there is nothing, 0 value returned." << endl;
        return ConfigValue("0");
    }

    string output = elements.front();
    elements.pop();

    return ConfigValue(output);
}

vector<ConfigValue> ConfigParser::TakeAll()
{
    vector<ConfigValue> output;

    while(elements.size())
    {
        output.push_back(ConfigValue(elements.front()));
        elements.pop();
    }

    return output;
}

string ConfigParser::comment_out(const string &str, size_t index)
{
    if(index >= comment_marks.size())
        return str;
    else {
        string str_co = comment_out(str, comment_marks.at(index));
        return comment_out(str_co, ++index);
    }
}

string ConfigParser::comment_out(const string &str, const string &c)
{
    if(c.empty()) {
        return str;
    }

    const auto commentBegin = str.find(c);
    return str.substr(0, commentBegin);
}


string ConfigParser::trim(const string &str, const string &w)
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

queue<string> ConfigParser::split(const string &str, const string &s)
{
    queue<string> eles;

    char *cstr = new char[str.length() + 1];

    strcpy(cstr, str.c_str());

    char *pch = strtok(&cstr[0], s.c_str());
    string ele;

    while(pch != nullptr)
    {
        ele = pch;
        eles.push(pch);
        pch = strtok(nullptr, s.c_str());
    }

    delete cstr;

    return eles;
}

int ConfigParser::find_integer(const string &str, const size_t &pos)
{
    vector<int> integers = find_integers(str);
    if(pos >= integers.size())
    {
        cerr << "Config Parser: Cannot find " << pos + 1 << " integers from "
             << "\"" << str << "\"."
             << endl;
        return 0;
    }

    return integers.at(pos);
}

vector<int> ConfigParser::find_integers(const string &str)
{
    vector<int> result;

    find_integer_helper(str, result);

    return result;
}

void ConfigParser::find_integer_helper(const string &str, vector<int> &result)
{
   if(str.empty())
       return;

   int negative = 1;
   auto numBeg = str.find_first_of("-0123456789");
   if(numBeg == string::npos)
       return;

   // check negative sign
   string str2 = str.substr(numBeg);

   if(str2.at(0) == '-')
   {
       negative = -1;
       int num_check;

       do {
           str2.erase(0, 1);

           if(str2.empty())
               return;

           num_check = str2.at(0) - '0';
       } while (num_check > 9 || num_check < 0);
   }

   auto numEnd = str2.find_first_not_of("0123456789");
   if(numEnd == string::npos)
       numEnd = str2.size();

   int num = 0;
   size_t i = 0;

   for(; i < numEnd; ++i)
   {
       if( (num > INT_MAX/10) ||
           (num == INT_MAX/10 && ((str2.at(i) - '0') > (INT_MAX - num*10))) )
       {
           ++i;
           break;
       }

       num = num*10 + str2.at(i) - '0';
   }

   result.push_back(negative*num);
   find_integer_helper(str2.substr(i), result);
}

ConfigParser &operator >> (ConfigParser &c, std::string &v)
{
    v = c.TakeFirst().String();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, char &v)
{
    v = c.TakeFirst().Char();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned char &v)
{
    v = c.TakeFirst().UChar();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, short &v)
{
    v = c.TakeFirst().Short();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned short &v)
{
    v = c.TakeFirst().UShort();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, int &v)
{
    v = c.TakeFirst().Int();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned int &v)
{
    v = c.TakeFirst().UInt();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long &v)
{
    v = c.TakeFirst().Long();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned long &v)
{
    v = c.TakeFirst().ULong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long long &v)
{
    v = c.TakeFirst().LongLong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned long long &v)
{
    v = c.TakeFirst().ULongLong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, float &v)
{
    v = c.TakeFirst().Float();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, double &v)
{
    v = c.TakeFirst().Double();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long double &v)
{
    v = c.TakeFirst().LongDouble();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, const char *&v)
{
    v = c.TakeFirst().c_str();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, ConfigValue &v)
{
    v = c.TakeFirst();
    return c;
}

