/* Grain production and manipulation
:q
:q
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
    vector<doublev_t> genGrain(doublev_t center, doublev_t dimensions,
            vector<doublev_t> basis, double latConst) {
        /* Generates a grain of the given size, using the given basis. Note that
         * 'basis' is intended to correspond to only one atom type, to be
         * combined
         * later.
         *
         * Args:
         *  center      -   the center of the grain
         *  dimensions  -   xyz size of the grain (same units as latConst)
         *  basis       -   the basis set
         *  latConst    -   the lattice constant of the unit cell
         *
         * Returns:
         *  grain   -   the coordinates of all atoms
         */

        int numCells;
        doublev_t temp;
        vector<doublev_t> grain = basis;

        // For all directions
        for (vector<doublev_t>::size_type i=0; i<dimensions.size(); ++i) {
            numCells = static_cast<int>(ceil(dimensions[i]/latConst));
            cout << "Basis size: " << basis.size() << endl;
            cout << "Num cells: " << numCells << endl;

            // For each cell
            for (int i=1; i<numCells; i++) {
                // For each atom in basis
                for (vector<doublev_t>::size_type j=0; j<basis.size(); ++j) {
                    temp = grain[j];
                    Tools::addToVector(temp,latConst*i);
                    grain.push_back(temp);
                }
            }
            basis = grain;
        }

        return grain;
    }
}
