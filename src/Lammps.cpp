/* Tools for reading/writing/editing LAMMPS style data files.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#include <vector>
#include <string>
#include <cstdio>
#include "define.h"

using namespace std;

namespace Lammps {

    void writeData(string filename, vector<dvec_t> arr) {
        /* Writes an array of atom information to a LAMMPS style data file.
         * Default output atom style is 'atomic'.
         *
         * Args:
         *  filename    -   name of output file
         *  arr         -   output array with format [type,x,y,z]
         */

        FILE * outfile;
        outfile = fopen(filename.c_str(), "w");

        // Comment lines in data file
        fprintf(outfile, "# Data file written by Lammps::writeData()\n");
        fprintf(outfile, "\n");

        // System info
        int nAtoms=0, nTypes=0;
        double xlo=0, xhi=0;
        double ylo=0, yhi=0;
        double zlo=0, zhi=0;

        dvec_t temp;

        nAtoms = static_cast<int>(arr.size());

        for (vector<dvec_t>::size_type i=0; i<arr.size(); ++i) {
            temp = arr[i];

            if (temp[0]>nTypes)
                nTypes = temp[0];

            if (temp[1]<xlo)
                xlo = temp[1];

            if (temp[1]>xhi)
                xhi = temp[1];

            if (temp[2]<ylo)
                ylo = temp[2];

            if (temp[2]>yhi)
                yhi = temp[2];

            if (temp[3]<zlo)
                zlo = temp[3];

            if (temp[3]>zhi)
                zhi = temp[3];
        }

        fprintf(outfile, "%d atoms\n", nAtoms);
        fprintf(outfile, "%d atom types\n", nTypes);
        fprintf(outfile, "\n");
        fprintf(outfile, "%f %f xlo xhi\n", xlo, xhi);
        fprintf(outfile, "%f %f ylo yhi\n", ylo, yhi);
        fprintf(outfile, "%f %f zlo zhi\n", zlo, zhi);
        fprintf(outfile, "\n");

        // Atoms
        fprintf(outfile, "Atoms # 'atomic'\n");
        fprintf(outfile, "\n");

        for (vector<dvec_t>::size_type i=0; i<arr.size(); ++i) {
            temp = arr[i];
            fprintf(outfile, "%d %d %f %f %f\n", static_cast<int>(i)+1,
                    static_cast<int>(temp[0]), temp[1], temp[2], temp[3]);
        }
    }

    vector<dvec_t> readData(string filename) {
        vector<dvec_t> temp;
        return temp;
    }
}
