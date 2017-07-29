#ifndef PV3D_H
#define PV3D_H

#include <vector>
#include "define.h"

namespace Pv3d {

    bool inRegion(dvec_t, vector<dvec_t> int);

    vector<dvec_t> genCenters(int, dvec_t);

    vector<dvec_t> genImages(vector<dvec_t>, dvec_t);
}

#endif
