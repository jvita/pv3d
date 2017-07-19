#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include "define.h"

using namespace std;

namespace Tools {

    vector<dvec_t> joinArrays(vector<dvec_t>, vector<dvec_t>);

    void addToVector(dvec_t&, double);

    void multiplyVector(dvec_t&, double);

    vector<dvec_t> dot(vector<dvec_t>, vector<dvec_t>);

    vector<dvec_t> dot(vector<dvec_t>, dvec_t);

    void printArr(vector<dvec_t>);

    vector<dvec_t> rotate(vector<dvec_t>, double, dvec_t);

}
#endif
