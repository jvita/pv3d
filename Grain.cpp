/* Grain production and manipulation
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

// TODO: create special Grain class (position, distance(), atoms)

#include <vector>
#include "define.h"

using namespace std;
typedef vector<double> doublev_t;

vector<doublev_t> genGrain(doublev_t center, vector<doublev_t> basis, 
                                double latConst) {
    /* Generates a grain of the given size, using the given basis. Note that
     * 'basis' is intended to correspond to only one atom type, to be combined
     * later.
     *
     * Args:
     *  center      -   the center of the grain
     *  basis       -   the basis set
     *  latConst    -   the lattice constant of the unit cell
     *
     * Returns:
     *  grain   -   the coordinates of all atoms
     */

    continue;
}
