#ifndef LAMMPS_H
#define LAMMPS_H

#include <vector>
#include <string>
#include "define.h"

namespace Lammps {

    void writeData(std::string, vector<dvec_t>);

    vector<dvec_t> readData(std::string);
}

#endif
