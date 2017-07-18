/* Initial software design for creating 3D periodic polycrystalline materials
 * using voronoi tesselation.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last Edited: 2017/7
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
typedef vector<double> doublev_t;

vector<doublev_t> dot(vector<doublev_t> A, vector<doublev_t> B) {
    /* Calculates the dot product of two matrices.
     *
     * Args:
     *  A   -   left matrix of size (nxm)
     *  B   -   right matrix of size (mxp)
     *
     * Returns: 
     *  result  -   dot product matrix (A*B) of size (nxp) 
     */

    int rows = A.size();
    int cols = B[0].size();

    vector<doublev_t> result;
    result.reserve(rows);

    double sum;

    for (int i=0; i<rows; i++) {
        doublev_t temp;
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

doublev_t dot(vector<doublev_t> A, doublev_t B) {
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

    doublev_t result;
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

void printArr(vector<doublev_t> arr) {
    /* Prints a 2D array to the console. */

    int len = arr.size();
    int wid = arr[0].size();

    for (vector<vector <double>>::size_type i=0; i<len; ++i) {
        for (vector<doublev_t>::size_type j=0; j<wid; ++j){
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }
}

vector<doublev_t> rotate(vector<doublev_t> arr, double theta, doublev_t axis) {
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
    vector<doublev_t> matrixK = {{0,-axis[2],axis[1]},
                                        {axis[2],0,-axis[0]},
                                        {-axis[1],axis[0],0}};

    vector<doublev_t> KK;
    KK.reserve(3);

    KK = dot(matrixK, matrixK);

    vector<doublev_t> rotMat;
    rotMat.reserve(3);

    double val;

    // Initialize the rotation matrix
    for (int i=0; i<3; i++) {
        doublev_t temp_v;
        temp_v.reserve(3);

        for (int j=0; j<3; j++) {
            val = matrixK[i][j]*sin(theta) + KK[i][j]*(1-cos(theta));

            if (i==j)
                val++;

            temp_v.push_back(val);
        }
        rotMat.push_back(temp_v);
    }

    vector<doublev_t> output;
    output.reserve(numPoints);

    for (vector<doublev_t>::size_type i=0; i<numPoints; ++i) {
        output.push_back(dot(rotMat, arr[i]));
    }

    return output;
}

int main() {

    doublev_t axis;

    vector<doublev_t> test_v;

    // TEST 1: z-axis rotation
    axis = {0,0,1};
    test_v.push_back(doublev_t {0,0,0});
    test_v.push_back(doublev_t {0,-1,0});
    test_v.push_back(doublev_t {-1,0,0});
    test_v.push_back(doublev_t {0,1,0});
    test_v.push_back(doublev_t {1,0,0});

    cout << "Before:" << endl;
    printArr(test_v);

    vector<doublev_t> output;
    output = rotate(test_v, M_PI/4, axis);

    cout << "After:" << endl;
    printArr(output);

    // TEST 2: y-axis rotation
    axis = {0,1,0};
    test_v.clear();
    test_v.push_back(doublev_t {0,0,0});
    test_v.push_back(doublev_t {-1,0,0});
    test_v.push_back(doublev_t {1,0,0});
    test_v.push_back(doublev_t {0,0,-1});
    test_v.push_back(doublev_t {0,0,1});

    cout << "Before:" << endl;
    printArr(test_v);

    output = rotate(test_v, M_PI/4, axis);

    cout << "After:" << endl;
    printArr(output);

    // TEST 3: x-axis rotation
    axis = {1,0,0};
    test_v.clear();
    test_v.push_back(doublev_t {0,0,0});
    test_v.push_back(doublev_t {0,1,0});
    test_v.push_back(doublev_t {0,-1,0});
    test_v.push_back(doublev_t {0,0,-1});
    test_v.push_back(doublev_t {0,0,1});

    cout << "Before:" << endl;
    printArr(test_v);

    output = rotate(test_v, M_PI/4, axis);

    cout << "After:" << endl;
    printArr(output);

    // TEST 3: arbitrary axis rotation
    axis = {1,1,1};
    test_v.clear();
    test_v.push_back(doublev_t {0,0,0});
    test_v.push_back(doublev_t {1,0,0});
    test_v.push_back(doublev_t {0,1,0});
    test_v.push_back(doublev_t {0,0,1});
    test_v.push_back(doublev_t {1,1,1});

    cout << "Before:" << endl;
    printArr(test_v);

    output = rotate(test_v, M_PI/4, axis);

    cout << "After:" << endl;
    printArr(output);


}
