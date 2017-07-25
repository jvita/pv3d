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
        /* Randomly generates 'nCenters' number of points within 'boxDims'.
         *
         * Args:
         *  nCenters    -   the number of points to generate
         *  boxDims     -   xyz bounds of box (assumes origin as lower bound)
         */

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

    vector<dvec_t> genImages(vector<dvec_t> originals) {
        /* Produce the 26 additional images (in 3D) of a set of
         * points
         *
         * Args:
         *  originals   -   the data points to be duplicated; assumes fractional
         *                  coordinates
         *
         * Returns:
         *  images      -   the full copy of all 26 duplicated images and the
         *                  original 1 set of points
         */

        vector<dvec_t> images;
        dvec_t toAdd;

        
        for (vector<dvec_t>::size_type a=0; a<originals.size(); a++) {
            toAdd = originals[0];

            for (double i=-1; i<2; i++) {
                for (double j=-1; j<2; j++) {
                    for (double k=-1; k<2; k++) {
                        dvec_t shift = {i,j,k};
                        Tools::addVectors(toAdd, shift);

                        images.push_back(toAdd);
                    }
                }
            }
        }
        return images;
    }
}

int main() {
   
    dvec_t boxDims = {20,20,20};
    double latConst = 5.0;
    vector<dvec_t> centers = Pv3d::genCenters(4, boxDims);
    
    vector<dvec_t> basis;
    vector<dvec_t> basis2;
    basis.reserve(8);
    basis2.reserve(8);

    basis.push_back(dvec_t {0,0,0});
    basis.push_back(dvec_t {0.5,0.5,0});
    basis.push_back(dvec_t {0,0.5,0.5});
    basis.push_back(dvec_t {0.5,0,0.5});

    basis2.push_back(dvec_t {0.5,0.5,0.5});
    basis2.push_back(dvec_t {0,0,0.5});
    basis2.push_back(dvec_t {0,0.5,0});
    basis2.push_back(dvec_t {0.5,0,0});

    vector<dvec_t> fullCrystal;
    dvec_t center;

    // Creates all grains
    for (vector<dvec_t>::size_type i=0; i<centers.size(); i++) {
        center = centers[i];

        Tools::printArr(center);

        vector<dvec_t> grain1;
        vector<dvec_t> grain2;

        grain1 = Grain::genGrain(center, boxDims, basis, latConst, 1.0);
        grain2 = Grain::genGrain(center, boxDims, basis2, latConst, 2.0);
        grain1 = Tools::joinArrays(grain1, grain2);

        double theta = rand()*2*M_PI / RAND_MAX;
        double x = static_cast<double>(rand()) / RAND_MAX;
        double y = static_cast<double>(rand()) / RAND_MAX;
        double z = static_cast<double>(rand()) / RAND_MAX;
        dvec_t axis = {x,y,z};

        // TODO: rotate doesn't expect atom types
        //Tools::rotate(grain1, theta, axis);

        if (fullCrystal.size()>0) {
            fullCrystal = Tools::joinArrays(fullCrystal, grain1);
        } else {
            fullCrystal = grain1;
        }
    }
   
    Lammps::writeData("data.test", fullCrystal);
}
