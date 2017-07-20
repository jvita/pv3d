/* Generates polycrystalline simulation cells for use in molecular dynamics and
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "define.h"
#include "Tools.h"
#include "Grain.h"
#include "Lammps.h"

using namespace std;

namespace Pv3d {

    bool inRegion(dvec_t p, vector<dvec_t> centers, int regionId) {
        /* Checks to see if point 'p' falls into the Voronoi tile specified by
         * regionId.
         *
         * Args:
         *  p           -   xyz coordinates of checked point
         *  centers     -   collection of all tile centers
         *  regionId    -   tile id to be checked
         *
         * Returns:
         *  'true' if point is in tile
         */

        dvec_t testCenter = centers[regionId];
        dvec_t diffToTest;

        for (dvec_t::size_type i=0; i<testCenter.size(); i++) {
            diffToTest.push_back(abs(testCenter[i]-p[i]));
        }

        double testDist = sqrt(diffToTest[0]*diffToTest[0] + diffToTest[1]*
                               diffToTest[1] + diffToTest[2]*diffToTest[2]);

        for (vector<dvec_t>::size_type i=0; i<centers.size(); i++) {
            dvec_t temp = centers[i];
            dvec_t diff;

            for (dvec_t::size_type j=0; j<temp.size(); j++) {
                diff.push_back(abs(temp[j]-p[j]));
            }

            double dist = sqrt(diff[0]*diff[0] + diff[1]*diff[1]
                               + diff[2]*diff[2] );

            if (dist < testDist)
                return false;
        }

        return true;
    }

    vector<dvec_t> genCenters(int nCenters, dvec_t boxDims) {
        vector<dvec_t> centers;

        srand(time(NULL));

        for (int i=0; i<nCenters; i++) {
            dvec_t temp;
            temp.reserve(3);

            for (int j=0; j<3; j++) {
                double rnd;
                
                rnd = static_cast<double>(rand()) / RAND_MAX;
                rnd *= boxDims[j];
                temp.push_back(rnd);
            }

            centers.push_back(temp);
        }

        return centers;
    }
}

int main() {

    //dvec_t axis;

    //vector<dvec_t> test_v;

    //// TEST 1: z-axis rotation
    //axis = {0,0,1};
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {0,-1,0});
    //test_v.push_back(dvec_t {-1,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {1,0,0});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //vector<dvec_t> output;
    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 2: y-axis rotation
    //axis = {0,1,0};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {-1,0,0});
    //test_v.push_back(dvec_t {1,0,0});
    //test_v.push_back(dvec_t {0,0,-1});
    //test_v.push_back(dvec_t {0,0,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 3: x-axis rotation
    //axis = {1,0,0};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {0,-1,0});
    //test_v.push_back(dvec_t {0,0,-1});
    //test_v.push_back(dvec_t {0,0,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 3: arbitrary axis rotation
    //axis = {1,1,1};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {1,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {0,0,1});
    //test_v.push_back(dvec_t {1,1,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 4: test scaleVector
    //dvec_t v = {1,2,3};

    //cout << "Old: ";
    //for (int i=0; i<3; i++) {
    //    cout << v[i] << " ";
    //}
    //cout << endl;

    //Tools::scaleVector(v,2);

    //cout << "New: ";
    //for (int i=0; i<3; i++) {
    //    cout << v[i] << " ";
    //}
    //cout << endl;
    //
    //basis.push_back(dvec_t {0.5,0.5,0.5});

    //basis2.push_back(dvec_t {0,0,0});
 
    //basis.push_back(dvec_t {1,0,0});
    //basis.push_back(dvec_t {0,1,0});
    //basis.push_back(dvec_t {0,0,1});
    //basis.push_back(dvec_t {1,1,0});
    //basis.push_back(dvec_t {0,1,1});
    //basis.push_back(dvec_t {1,0,1});
    //basis.push_back(dvec_t {1,1,1});
    
    vector<dvec_t> basis;
    vector<dvec_t> basis2;
    basis.reserve(8);
    basis2.reserve(8);
    
    basis.push_back(dvec_t {0,0,0});
    basis2.push_back(dvec_t {0.5,0.5,0.5});

    dvec_t center = {0,0,0};
    dvec_t dimensions = {10,10,10};
    double latConst = 5.0;

    // TODO: trim off atom type data for point comparison? inRegion()

    //grain1 = Grain::genGrain(center, dimensions, basis, latConst, 1.0);
    //grain2 = Grain::genGrain(center, dimensions, basis2, latConst, 2.0);
    //grain1 = Tools::joinArrays(grain1, grain2);

    //cout << "Before rotation" << endl;
    //Tools::printArr(grain1);

    //dvec_t axis = {0,0,1};
    //grain = Tools::rotate(grain, M_PI/2, axis);

    //cout << "After rotation" << endl;
    //Tools::printArr(grain);
    
    //Lammps::writeData("data.test", grain1);

    //int nCenters = 5;
    //vector<dvec_t> centers = Pv3d::genCenters(nCenters, dimensions);

    //Tools::printArr(centers);
}
