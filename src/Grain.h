#ifndef GRAIN_H
#define GRAIN_H

#include "define.h"

namespace Grain {

    dvec_t getGrainCenter(vector<dvec_t>);

    vector<dvec_t> genGrain(dvec_t, vector<dvec_t>, double, double);

    void shiftGrain(vector<dvec_t>&, dvec_t);
}

#endif
