/* Generates polycrystalline simulation cells for use in molecular dynamics and
 * first principles calculations.
 *
 * Author: Josh Vita
 * Created: 2017/7
 * Last edited: 2017/7
 */

#include <iostream>
#include <cmath>
#include "define.h"
#include "Tools.h"
#include "Grain.h"

using namespace std;

int main() {

    //dvec_t axis;

    //vector<dvec_t> test_v;

    //// TEST 1: z-axis rotation
    //axis = {0,0,1};
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {0,-1,0});
    //test_v.push_back(dvec_t {-1,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {1,0,0});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //vector<dvec_t> output;
    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 2: y-axis rotation
    //axis = {0,1,0};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {-1,0,0});
    //test_v.push_back(dvec_t {1,0,0});
    //test_v.push_back(dvec_t {0,0,-1});
    //test_v.push_back(dvec_t {0,0,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 3: x-axis rotation
    //axis = {1,0,0};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {0,-1,0});
    //test_v.push_back(dvec_t {0,0,-1});
    //test_v.push_back(dvec_t {0,0,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 3: arbitrary axis rotation
    //axis = {1,1,1};
    //test_v.clear();
    //test_v.push_back(dvec_t {0,0,0});
    //test_v.push_back(dvec_t {1,0,0});
    //test_v.push_back(dvec_t {0,1,0});
    //test_v.push_back(dvec_t {0,0,1});
    //test_v.push_back(dvec_t {1,1,1});

    //cout << "Before:" << endl;
    //Tools::printArr(test_v);

    //output = Tools::rotate(test_v, M_PI/4, axis);

    //cout << "After:" << endl;
    //Tools::printArr(output);

    //// TEST 4: test multiplyVector
    dvec_t v = {1,2,3};

    cout << "Old: ";
    for (int i=0; i<3; i++) {
        cout << v[i] << " ";
    }
    cout << endl;

    Tools::multiplyVector(v,2);

    cout << "New: ";
    for (int i=0; i<3; i++) {
        cout << v[i] << " ";
    }
    cout << endl;
    
    vector<dvec_t> basis;
    basis.reserve(8);

    basis.push_back(dvec_t {0,0,0});
    //basis.push_back(dvec_t {1,0,0});
    //basis.push_back(dvec_t {0,1,0});
    //basis.push_back(dvec_t {0,0,1});
    //basis.push_back(dvec_t {1,1,0});
    //basis.push_back(dvec_t {0,1,1});
    //basis.push_back(dvec_t {1,0,1});
    //basis.push_back(dvec_t {1,1,1});

    dvec_t center = {0,0,0};
    dvec_t dimensions = {10,10,10};
    double latConst = 5.0;

    vector<dvec_t> grain;

    grain = Grain::genGrain(center, dimensions, basis, latConst);

    Tools::printArr(grain);
    cout << grain.size() << endl;
}
