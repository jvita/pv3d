/* Functions for manipulating arrays; intended for use with periodic voronoi
 * tesselation code.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <cmath>
#include "define.h"

using namespace std;

namespace Tools {

    void addToVector(dvec_t &v, double x) {
        for (vector<dvec_t>::size_type i=0; i<v.size(); ++i) {
            v[i] = v[i]+x;
        }
    }

    void multiplyVector(dvec_t &v, double x) {
        for (vector<dvec_t>::size_type i=0; i<v.size(); ++i) {
            v[i] = v[i]*x;
        }
    }

    vector<dvec_t> dot(vector<dvec_t> A, vector<dvec_t> B) {
        /* Calculates the dot product of two matrices.
         *
         * Args:
         *  A   -   left matrix of size (nxm)
         *  B   -   right matrix of size (mxp)
         *
         * Returns: 
         *  result  -   dot product matrix (A*B) of size (nxp) 
         */
        cout << "POOPSICLES" << endl;

        int rows = A.size();
        int cols = B[0].size();

        vector<dvec_t> result;
        result.reserve(rows);

        double sum;

        for (int i=0; i<rows; i++) {
            dvec_t temp;
            temp.reserve(cols);

            for (int j=0; j<cols; j++) {
                sum = 0;
                for (int k=0; k<rows; k++) {
                    sum += A[i][k]*B[k][j];
                }
                temp.push_back(sum);
            }
            result.push_back(temp);
        }

        return result;
    }

    dvec_t dot(vector<dvec_t> A, dvec_t B) {
        /* Calculates the dot product of two matrices. Function overloading for use
         * with 1D right matrix.
         *
         * Args:
         *  A   -   left matrix of size (nxm)
         *  B   -   right matrix of size (mxp)
         *
         * Returns: 
         *  result  -   dot product matrix (A*B) of size (nxp) 
         */

        int rows = A.size();
        int cols = B.size();

        dvec_t result;
        result.reserve(rows);

        double sum;

        for (int i=0; i<rows; i++) {
            for (int j=0; j<cols; j++) {
                sum = 0;
                for (int k=0; k<rows; k++) {
                    sum += A[i][k]*B[k];
                }
            }
            result.push_back(sum);
        }

        return result;
    }

    void printArr(vector<dvec_t> arr) {
        /* Prints a 2D array to the console. */

        vector<dvec_t>::size_type len = arr.size();
        vector<dvec_t>::size_type wid = arr[0].size();

        for (vector<dvec_t>::size_type i=0; i<len; ++i) {
            for (vector<dvec_t>::size_type j=0; j<wid; ++j){
                cout << arr[i][j] << " ";
            }
            cout << endl;
        }
    }

    vector<dvec_t> rotate(vector<dvec_t> arr, double theta, dvec_t axis) {
        /* Rotates a multidimensional array about a given axis by a given theta.
         * Uses the Rodrigues rotation formula.
         *
         * Args:
         *  arr     - the array to be rotated
         *  theta   - angle in radians
         *  axis    - axis of rotation
         */

        int numPoints = arr.size();

        // Initialize the K matrix
        vector<dvec_t> matrixK = {{0,-axis[2],axis[1]},
                                            {axis[2],0,-axis[0]},
                                            {-axis[1],axis[0],0}};

        vector<dvec_t> KK;
        KK.reserve(3);

        KK = dot(matrixK, matrixK);

        vector<dvec_t> rotMat;
        rotMat.reserve(3);

        double val;

        // Initialize the rotation matrix
        for (int i=0; i<3; i++) {
            dvec_t temp_v;
            temp_v.reserve(3);

            for (int j=0; j<3; j++) {
                val = matrixK[i][j]*sin(theta) + KK[i][j]*(1-cos(theta));

                if (i==j)
                    val++;

                temp_v.push_back(val);
            }
            rotMat.push_back(temp_v);
        }

        vector<dvec_t> output;
        output.reserve(numPoints);

        for (vector<dvec_t>::size_type i=0; i< static_cast<vector<dvec_t>::
                size_type>(numPoints); ++i) {
            output.push_back(dot(rotMat, arr[i]));
        }

        return output;
    }
}