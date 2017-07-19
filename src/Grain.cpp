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
    vector<dvec_t> genGrain(dvec_t center, dvec_t dimensions,
            vector<dvec_t> basis, double latConst, double type) {
        /* Generates a grain of the given size, using the given basis. Note that
         * 'basis' is intended to correspond to only one atom type, to be
         * combined
         * later.
         *
         * Args:
         *  center      -   the center of the grain
         *  dimensions  -   xyz size of the grain (same units as latConst)
         *  basis       -   the basis set in the format [type x y z]
         *  latConst    -   the lattice constant of the unit cell
         *  type        -   atom type
         *
         * Returns:
         *  grain   -   the coordinates of all atoms
         */

        int numCells;
        dvec_t temp;

        // Initialize basis types
        for (vector<dvec_t>::size_type i=0; i<basis.size(); i++) {
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
                    temp[i+1] = latConst*j;
                    grain.push_back(temp);
                }
            }
            basis.resize(grain.size());
            basis = grain;
        }

        // Shift each atom by 'center'
        for (vector<dvec_t>::size_type i=0; i<grain.size()-1; i++) {
            // j = direction
            for (vector<dvec_t>::size_type j=0; j<grain[i].size()-1; j++) {
                grain[i][j+1] += center[j];
            }
        }

        return grain;
    }
}
