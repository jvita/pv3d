/* Generates polycrystalline simulation cells for use in molecular dynamics and
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#include <iostream>
#include <cmath>
#include "Tools.h"
#include "define.h"

using namespace std;

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
