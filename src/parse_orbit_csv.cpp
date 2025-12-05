#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <boost/tokenizer.hpp>

enum OrbitType{
    UNKNOWN = 0,
    RESON = 1,
    // liberation point Orbits
    LYAPONOV = 16,
    HALO = 17,
    VERTICAL = 18,
    AXIAL = 19,
    SHORT = 20,
    LONG = 21,
    BUTTERFLY = 22,
    DRAGONFLY = 23,
    // P2-Centered Orbits
    DRO = 24,
    DPO = 25,
    LPO = 26
};

struct OrbitRecord
{
    int id;
    std::array<double,6> state;
    double jacobi;
    double period_TU;
    double period_days;
    double stability_index;
    double mass_ratio;
    OrbitType orbit_type;
};

bool parseCSVRow(const std::string& line, OrbitRecord& record, OrbitType type)
{    
    boost::escaped_list_separator<char> sep('\\', ',', '"');
    boost::tokenizer<boost::escaped_list_separator<char>> tok(line, sep);

    std::vector<std::string> fields(tok.begin(), tok.end());

    if (fields.size() != 12)
        return false;

    try
    {
        record.id              = std::stoi(fields[0]);
        record.state[0]        = std::stod(fields[1]);
        record.state[1]        = std::stod(fields[2]);
        record.state[2]        = std::stod(fields[3]);
        record.state[3]        = std::stod(fields[4]);
        record.state[4]        = std::stod(fields[5]);
        record.state[5]        = std::stod(fields[6]);
        record.jacobi          = std::stod(fields[7]);
        record.period_TU       = std::stod(fields[8]);
        record.period_days     = std::stod(fields[9]);
        record.stability_index = std::stod(fields[10]);
        record.mass_ratio      = std::stod(fields[11]);
        record.orbit_type      = type;
    }
    catch (...)
    {
        return false;
    }

    return true;
}

std::vector<OrbitRecord> parse_orbit(std::string path, OrbitType type = OrbitType::UNKNOWN){
    std::ifstream file(path);
    
    if (!file)
    {
        std::cerr << "ERROR: Could not open file.\n";
        exit(-1);
    }
    
    std::vector<OrbitRecord> orbits;
    std::string line;

    std::getline(file, line); //first line is the header
    
    while (std::getline(file, line))
    {
        OrbitRecord record;
        if (parseCSVRow(line, record, type))
        {
            orbits.push_back(record);
        }
        else
        {
            std::cerr << "WARNING: Failed to parse line:\n" << line << "\n";
        }
    }
    file.close();

    return orbits;
}

void print_max_instability(std::vector<OrbitRecord> orbits){
    double max = 0; //min = 1 for real orbits
    OrbitRecord most_unstable;
    for(OrbitRecord orbit:orbits){
        if(orbit.stability_index > max){
            max = orbit.stability_index;
            most_unstable = orbit;
        }
    }
    std::cout << "id:" << most_unstable.id << ",stability_index=" << most_unstable.stability_index <<  '\n';
}

int main()
{
    auto orbits = parse_orbit("orbits/periodic_orbits_halo.csv");
    print_max_instability(orbits);
    orbits = parse_orbit("orbits/periodic_orbits_dro.csv");
    print_max_instability(orbits);
    orbits = parse_orbit("orbits/periodic_orbits_lyapynov.csv");
    print_max_instability(orbits);
    orbits = parse_orbit("orbits/periodic_orbits_vertical.csv");
    print_max_instability(orbits);
    
    return 0;
}
