/* Grain production and manipulation for use with molecular dynamics runs and
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

// TODO: create special Grain class (position, distance(), atoms)

#include <vector>
#include <cmath>
#include "define.h"
#include <iostream>
#include "Tools.h"

using namespace std;

namespace Grain {

    dvec_t getGrainCenter(vector<dvec_t> arr) {
        /* Searches the array of atom coordinates to get the true center of the
         * grain by calculating the box dimension.
         * 
         * Args:
         *  arr    -    the array of atom coordinates WITH atom types
         */

        // TODO: testing

        double xlo=0,xhi=0;
        double ylo=0,yhi=0;
        double zlo=0,zhi=0;

        for (vector<dvec_t>::size_type i=0; i<arr.size(); i++) {
            dvec_t temp = arr[i];

            if (temp[1]<xlo){
                xlo = temp[1];
            }else if (temp[1]>xhi){
                xhi = temp[1];
            }
            
            if (temp[2]<ylo){
                ylo = temp[2];
            }else if (temp[2]>yhi){
                yhi = temp[2];
            }
 
            if (temp[3]<zlo){
                zlo = temp[3];
            }else if (temp[3]>zhi){
                zhi = temp[3];
            }
        }

        double midx = abs(xhi-xlo)/2.0;
        double midy = abs(yhi-ylo)/2.0;
        double midz = abs(zhi-zlo)/2.0;

        dvec_t trueCenter = {midx,midy,midz};

        return trueCenter;
    }

    void shiftGrain(vector<dvec_t> &arr, dvec_t center) {
        /* Shifts an array of points to a new center
         *
         * Args:
         *  arr    -    array to be shifted
         *  center -    new center
         */

        for (vector<dvec_t>::size_type i=0; i<arr.size(); i++) {
            // j = direction
            for (vector<dvec_t>::size_type j=0; j<arr[i].size()-1; j++) {
                arr[i][j+1] += center[j];
            }
        }
    }

    vector<dvec_t> genGrain(dvec_t dimensions, vector<dvec_t> basis,
                                double latConst, double type) {
        /* Generates a grain of the given size, using the given basis. Note that
         * 'basis' is intended to correspond to only one atom type, to be
         * combined
         * later.
         *
         * Args:
         *  dimensions  -   xyz size of the grain (same units as latConst)
         *  basis       -   the basis set in the format [type x y z]
         *                  (fractional coordinates).
         *  latConst    -   the lattice constant of the unit cell
         *  type        -   atom type
         *
         * Returns:
         *  grain   -   an origin-centered grain
         */

        int numCells;
        dvec_t temp;

        // Initialize basis types and scale basis
        for (vector<dvec_t>::size_type i=0; i<basis.size(); i++) {
            Tools::scaleVector(basis[i], latConst);
            basis[i].insert(basis[i].begin(), type);
        }

        vector<dvec_t> grain = basis;

        // For all directions
        for (vector<dvec_t>::size_type i=0; i<dimensions.size(); i++) {
            numCells = static_cast<int>(ceil(dimensions[i]/latConst));

            // For each cell
            for (int j=1; j<numCells; j++) {
                // For each atom in basis
                for (vector<dvec_t>::size_type k=0; k<basis.size(); k++) {
                    temp = basis[k];
                    //temp[i+1] = latConst*j;
                    temp[i+1] += j*latConst;
                    grain.push_back(temp);
                }
            }
            basis.resize(grain.size());
            basis = grain;
        }

        dvec_t trueCenter = getGrainCenter(grain);
        Tools::scaleVector(trueCenter, -1);
        shiftGrain(grain, trueCenter);

        return grain;
    }
}
