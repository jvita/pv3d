/* Generates polycrystalline simulation cells for use in molecular dynamics and
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#include <time.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <string>
#include "define.h"
#include "Tools.h"
#include "Grain.h"
#include "Lammps.h"


//TODO: genGrain, use the cube that encapsulates the sphere... that encapsulates
//the smaller cube
using namespace std;

namespace Pv3d {

    bool inBox(dvec_t p, vector<dvec_t> boxDims) {
        /* Returns 'true' if the point falls within the rectangular box
         *
         * Args:
         *  p           -   the point to check
         *  boxDims     -   box dimensions; organized as {(xlo,xhi), (ylo,yhi),
         *                  (zlo,zhi)}
         *
         * Returns:
         *  true if point is in box
         */

        bool in = ( (p[1] >= boxDims[0][0] && p[1] <= boxDims[0][1]) && 
                    (p[2] >= boxDims[1][0] && p[2] <= boxDims[1][1]) &&
                    (p[3] >= boxDims[2][0] && p[3] <= boxDims[2][1]));

        return in;
    }

    bool inRegion(dvec_t p, vector<dvec_t> centers, int regionId) {
        /* Checks to see if point 'p' falls into the Voronoi tile specified by
         * regionId.
         *
         * Args:
         *  p           -   xyz coordinates of checked point (with atom type
         *                  info)
         *  centers     -   collection of all tile centers
         *  regionId    -   tile id to be checked
         *
         * Returns:
         *  'true' if point is in tile
         */

        // Instantiate to first distance
        dvec_t t1 = centers[0];
        dvec_t diff = {t1[0]-p[1], t1[1]-p[2], t1[2]-p[3]};     // p has atom
                                                                // type info
        double smallestDist = sqrt(diff[0]*diff[0] + diff[1]*diff[1] +
                                        diff[2]*diff[2]);
        int closestCenter = 0;

        // Compare against every other distance
        for (vector<dvec_t>::size_type i=1; i<centers.size(); i++) {
            t1 = centers[i];
            dvec_t diff = {t1[0]-p[1], t1[1]-p[2], t1[2]-p[3]};

            double checkDist = sqrt(diff[0]*diff[0] + diff[1]*diff[1] +
                                        diff[2]*diff[2]);

            // Update distance and ID if closer
            if (checkDist < smallestDist) {
                smallestDist = checkDist;

                // Each block of 27 corresponds to one 'family' of regions
                closestCenter = static_cast<int>(floor(i/27));
            }
        }

        // Return results
        if (closestCenter == regionId)
            return true;
        else {
            return false;
        }
    }

    vector<dvec_t> genCenters(int nCenters, dvec_t boxDims) {
        /* Randomly generates 'nCenters' number of points within 'boxDims'.
         *
         * Args:
         *  nCenters    -   the number of points to generate
         *  boxDims     -   xyz bounds of box (assumes origin as lower bound)
         *
         * Returns:
         *  centers     -   a set of xyz coordinates (no atom info)
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

    vector<dvec_t> genImages(vector<dvec_t> originals, dvec_t boxDims) {
        /* Produce the 26 additional images (in 3D) of a set of
         * points
         *
         * Args:
         *  originals   -   the data points to be duplicated; assumes fractional
         *                  coordinates
         *  boxDims     -   dimensions of box in xyz directions
         *
         * Returns:
         *  images      -   the full copy of all 26 duplicated images and the
         *                  original 1 set of points
         */

        vector<dvec_t> images;
        dvec_t toAdd;

        for (vector<dvec_t>::size_type a=0; a<originals.size(); a++) {
            toAdd = originals[a];

            for (double i=-1; i<2; i++) {
                for (double j=-1; j<2; j++) {
                    for (double k=-1; k<2; k++) {
                        dvec_t shift = {i*boxDims[0],j*boxDims[1],k*boxDims[2]};
                        dvec_t temp = toAdd;
                        Tools::addVectors(temp, shift);

                        images.push_back(temp);
                    }
                }
            }
        }
        return images;
    }
}

int main() {

    clock_t t1 = clock();
   
    // TODO: genImages needs to shift by boxDims; make adaptable to rectangles
    double sideLength;

    cout << "Box side length: ";
    cin >> sideLength;

    dvec_t boxDims = {sideLength,sideLength,sideLength};
    vector<dvec_t> boxMinMax;

    boxMinMax.push_back(dvec_t {0, sideLength});
    boxMinMax.push_back(dvec_t {0, sideLength});
    boxMinMax.push_back(dvec_t {0, sideLength});

    double latConst;
    int numGrains;
    string fname;

    cout << "Lattice constant: ";
    cin >> latConst;
    cout << "Number of grains: ";
    cin >> numGrains;
    cout << "Output file name: ";
    cin >> fname;

    vector<dvec_t> centers = Pv3d::genCenters(numGrains, boxDims);
    
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

    vector<dvec_t> images = Pv3d::genImages(centers, boxDims);

    // Iterate over each family of regions
    for (int j=0; j<numGrains; j++) {

        double theta = rand()*2*M_PI / RAND_MAX;
        double x = static_cast<double>(rand()) / RAND_MAX;
        double y = static_cast<double>(rand()) / RAND_MAX;
        double z = static_cast<double>(rand()) / RAND_MAX;
        dvec_t axis = {x,y,z};

        // Iterate over all 27 images
        for (int i=j*27; i<(j+1)*27; i++) {

            center = images[i];
            j = static_cast<int>(floor(i/27)); // images are in blocks of 27


            // Here grain1/2 represent the 'grains' of the 1st and 2nd basis
            // For more/less bases, delete grain2 or add grain3, grain4, ...
            vector<dvec_t> grain1;
            vector<dvec_t> grain2;

            // 1.0 and 2.0 should correspond to grain1 grain2
            grain1 = Grain::genGrain(boxDims, basis, latConst, 1.0);
            grain2 = Grain::genGrain(boxDims, basis2, latConst, 2.0);
            grain1 = Tools::joinArrays(grain1, grain2);

            Tools::rotate(grain1,theta,axis);

            Grain::shiftGrain(grain1, center);

            for (vector<dvec_t>::size_type a=0; a<grain1.size(); a++) {
                if (Pv3d::inBox(grain1[a], boxMinMax) &&
                        Pv3d::inRegion(grain1[a], images, j)) {
                    fullCrystal.push_back(grain1[a]);
                }
            }
        }
    }

    Lammps::writeData(fname, fullCrystal);

    clock_t t2 = clock();
    float diff = static_cast<float>(t2)-static_cast<float>(t1);

    cout << "Runtime: " << diff/CLOCKS_PER_SEC << " seconds"<< endl;
}
