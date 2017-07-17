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

vector<vector<double> > dot(vector<vector<double> > A, vector<vector<double> > B) {
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

    vector<vector<double> > result;
    result.reserve(rows);

    double sum;

    for (int i=0; i<rows; i++) {
        vector<double> temp;
        temp.reserve(cols);

        for (int j=0; j<cols; j++) {
            sum = 0;
            for (int k=0; k<rows; k++) {
                sum += A[i][k]*B[k][j];
            }
            temp.push_back(sum);
        }
        result.push_back(temp);
        //cout << "adding temp of size: " << temp.size() << endl;
        //cout << "new size: " << result.size() << endl;
    }

    return result;
}

vector<double> dot(vector<vector<double> > A, vector<double> B) {
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

    vector<double> result;
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
        //cout << "adding temp of size: " << temp.size() << endl;
        //cout << "new size: " << result.size() << endl;
    }

    return result;
}
void printArr(vector<vector<double> > arr) {
    /* Prints a 2D array to the console. */

    int len = arr.size();
    int wid = arr[0].size();

    //cout << len << endl;
    //cout << wid << endl;
    for (vector<vector <double> >::size_type i=0; i<len; ++i) {
        for (vector<vector<double> >::size_type j=0; j<wid; ++j){
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<double> > rotate(vector<vector<double> > arr, double theta, vector<double> axis) {
    /* Rotates a multidimensional array about a given axis by a given theta.
     * Uses the Rodrigues rotation formula.
     *
     * Args:
     *  arr     - the array to be rotated
     *  theta   - angle in radians
     *  axis    - axis of rotation
     */

    // Initialize the K matrix
    //double matrixK[3][3] = {{0,-axis[0],axis[1]},
    //                        {axis[2],0,-axis[0]},
    //                        {-axis[1],axis[0],0}};
    //

    int numPoints = arr.size();

    vector<vector<double> > matrixK = {{0,-axis[0],axis[1]},
                                        {axis[2],0,-axis[0]},
                                        {-axis[1],axis[0],0}};

    cout << "matrixK:" << endl;
    printArr(matrixK);

    vector<vector<double> > KK;
    KK.reserve(3);

    KK = dot(matrixK, matrixK);

    cout << "KK:" << endl;
    printArr(KK);

    vector<vector<double> > rotMat;
    rotMat.reserve(3);

    double tmp;
    double val;

    // Initialize the rotation matrix
    for (int i=0; i<3; i++) {
        vector<double> temp_v;
        temp_v.reserve(3);

        for (int j=0; j<3; j++) {
            tmp = matrixK[i][j];
            
            val = tmp*sin(theta) + KK[i][j]*(1-cos(theta));

            if (i==j)
                val++;

            temp_v.push_back(val);
        }
        //cout << "adding temp_v of size: " << temp_v.size() << endl;
        rotMat.push_back(temp_v);
        //cout << "new size: " << rotMat.size() << endl;
    }

    cout << "rotMat:" << endl;
    printArr(rotMat);

    vector<vector<double> > output;
    output.reserve(numPoints);


    // Rotate the original array
    //for (int i=0; i<3; i++) {
    //    for (int j=0; j<3; j++) {
    //        for (int k=0; k<arrLen; k++) {
    //            temp_v[k] = arr[i][k]*rotMat[k][j];
    //        }
    //        output.push_back(temp_v);
    //    }
    //}
    
    for (vector<vector<double> >::size_type i=0; i<numPoints; ++i) {
        output.push_back(dot(rotMat, arr[i]));
    }

    return output;
}

int main() {

    //double *test[5];

    //test[0] = new double[3] {1,0,0};
    //test[1] = new double[3] {0,1,0};
    //test[2] = new double[3] {-1,0,0};
    //test[3] = new double[3] {0,-1,0};
    //test[4] = new double[3] {0,0,0};

    vector<double> axis = {0,0,1};

    vector<vector<double> > test_v;

    test_v.push_back(vector<double> {0,0,0});
    test_v.push_back(vector<double> {0,-1,0});
    test_v.push_back(vector<double> {-1,0,0});
    test_v.push_back(vector<double> {0,1,0});
    test_v.push_back(vector<double> {1,0,0});

    cout << "Before:" << endl;
    printArr(test_v);

    vector<vector<double> > output;
    output = rotate(test_v, M_PI/4, axis);

    cout << "After:" << endl;
    printArr(output);

    //for (int i=0; i<5; i++) {
    //    delete test[i];
    //}
}
