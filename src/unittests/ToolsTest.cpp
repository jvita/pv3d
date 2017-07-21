#include "UnitTest++/UnitTest++.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "define.h"
#include "Tools.h"

const double tolerance = 1e-10;

using namespace std;

SUITE(joinArrays) {
    TEST(joinEqualSized) {
        vector<dvec_t> v1;
        vector<dvec_t> v2;

        v1.push_back(dvec_t {1,1,1});
        v1.push_back(dvec_t {2,2,2});
        v1.push_back(dvec_t {3,3,3});

        v2.push_back(dvec_t {4,4,4});
        v2.push_back(dvec_t {5,5,5});
        v2.push_back(dvec_t {6,6,6});

        vector<dvec_t> ans_v;

        ans_v.push_back(dvec_t {1,1,1});
        ans_v.push_back(dvec_t {2,2,2});
        ans_v.push_back(dvec_t {3,3,3});
        ans_v.push_back(dvec_t {4,4,4});
        ans_v.push_back(dvec_t {5,5,5});
        ans_v.push_back(dvec_t {6,6,6});

        vector<dvec_t> joined = Tools::joinArrays(v1,v2);

        CHECK_ARRAY2D_CLOSE(ans_v,joined,6,3,tolerance);
    }

    TEST(joinDifferentSized) {
        vector<dvec_t> v1;
        vector<dvec_t> v2;

        v1.push_back(dvec_t {1,1,1});
        v1.push_back(dvec_t {2,2,2});
        v1.push_back(dvec_t {3,3,3});
        v1.push_back(dvec_t {4,4,4});

        v2.push_back(dvec_t {5,5,5});
        v2.push_back(dvec_t {6,6,6});

        vector<dvec_t> ans_v;

        ans_v.push_back(dvec_t {1,1,1});
        ans_v.push_back(dvec_t {2,2,2});
        ans_v.push_back(dvec_t {3,3,3});
        ans_v.push_back(dvec_t {4,4,4});
        ans_v.push_back(dvec_t {5,5,5});
        ans_v.push_back(dvec_t {6,6,6});

        vector<dvec_t> joined = Tools::joinArrays(v1,v2);

        CHECK_ARRAY2D_CLOSE(ans_v,joined,6,3,tolerance);
    }
}

SUITE(addVectors) {
    TEST(testEquals) {
        dvec_t ans;
        dvec_t v1 = {0,0,0};

        ans = {0,0,0};
        CHECK_ARRAY_EQUAL(ans,v1,3);
    }

    TEST(addPositiveInt) {
        dvec_t ans;
        dvec_t add;
        dvec_t v1 = {0,0,0};

        add = {1,1,1};
        ans = {1,1,1};
        Tools::addVectors(v1, add);
        CHECK_ARRAY_EQUAL(ans,v1,3);
    }

    TEST(addNegativeInt) {

        dvec_t ans;
        dvec_t add;
        dvec_t v1 = {0,0,0};

        add = {-2,-2,-2};
        ans = {-2,-2,-2};
        Tools::addVectors(v1,add);
        CHECK_ARRAY_EQUAL(ans,v1,3);
    }

    TEST(addPositiveFloat) {

        dvec_t ans;
        dvec_t add;
        dvec_t v1 = {0,0,0};

        add = {1.5,1.5,1.5};
        ans = {1.5,1.5,1.5};
        Tools::addVectors(v1,add);
        CHECK_ARRAY_EQUAL(ans,v1,3);
    }
}

SUITE(scaleVector) {
    TEST(scalePositive1D) {
        dvec_t vec = {1,2,3};
        double scale = 2.5;

        dvec_t ans = {2.5,5.0,7.5};
        
        Tools::scaleVector(vec, scale);
        CHECK_ARRAY_CLOSE(ans,vec,3,tolerance);
    }

    TEST(scaleNegative1D) {
        dvec_t vec = {1,2,3};
        double scale = -2.5;

        dvec_t ans = {-2.5,-5.0,-7.5};
        
        Tools::scaleVector(vec, scale);
        CHECK_ARRAY_CLOSE(ans,vec,3,tolerance);
    }
}

SUITE(dotWithArrays) {
    TEST(basicAdotB) {
        vector<dvec_t> A;
        vector<dvec_t> B;

        A.push_back(dvec_t {1,1,1});
        A.push_back(dvec_t {1,1,1});
        A.push_back(dvec_t {1,1,1});

        B.push_back(dvec_t {2,0,0});
        B.push_back(dvec_t {0,3,0});
        B.push_back(dvec_t {0,0,4});

        vector<dvec_t> dot = Tools::dot(A,B);

        vector<dvec_t> ans;

        ans.push_back(dvec_t {2,3,4});
        ans.push_back(dvec_t {2,3,4});
        ans.push_back(dvec_t {2,3,4});

        CHECK_ARRAY2D_CLOSE(ans,dot,3,3,tolerance);
    }
}

SUITE(dotArrayAndVector) {
    TEST(basicAdotBvec) {
        vector<dvec_t> A;
        dvec_t B;

        A.push_back(dvec_t {1,1,1});
        A.push_back(dvec_t {1,1,1});
        A.push_back(dvec_t {1,1,1});

        B = {2,3,4};

        dvec_t dot = Tools::dot(A,B);

        dvec_t ans;

        ans = {9,9,9};

        CHECK_ARRAY_CLOSE(ans,dot,3,tolerance);
    }
}

SUITE(rotate) {
    class BasisFixture {
        public:
            vector<dvec_t> test_v;
            dvec_t axis = {0,0,1};
            double rt2 = sqrt(2)/2.0;

            BasisFixture() {
                test_v.push_back(dvec_t {0,0,0});
                test_v.push_back(dvec_t {0,-1,0});
                test_v.push_back(dvec_t {-1,0,0});
                test_v.push_back(dvec_t {0,1,0});
                test_v.push_back(dvec_t {1,0,0});
            }
    };

    TEST_FIXTURE(BasisFixture, pi4Rotate) {
        Tools::rotate(test_v, M_PI/4, axis);

        vector<dvec_t> ans_v;
        ans_v.push_back(dvec_t {0,0,0});
        ans_v.push_back(dvec_t {rt2,-rt2,0});
        ans_v.push_back(dvec_t {-rt2,-rt2,0});
        ans_v.push_back(dvec_t {-rt2,rt2,0});
        ans_v.push_back(dvec_t {rt2,rt2,0});

        CHECK_ARRAY2D_CLOSE(ans_v, test_v, 5, 3, tolerance);
    }

    TEST_FIXTURE(BasisFixture, pi2Rotate) {
        Tools::rotate(test_v, M_PI/2, axis);

        vector<dvec_t> ans_v;
        ans_v.push_back(dvec_t {0,0,0});
        ans_v.push_back(dvec_t {1,0,0});
        ans_v.push_back(dvec_t {0,-1,0});
        ans_v.push_back(dvec_t {-1,0,0});
        ans_v.push_back(dvec_t {0,1,0});

        CHECK_ARRAY2D_CLOSE(ans_v, test_v, 5, 3, tolerance);
    }

    TEST_FIXTURE(BasisFixture, piRotate) {
        Tools::rotate(test_v, M_PI, axis);

        vector<dvec_t> ans_v;
        ans_v.push_back(dvec_t {0,0,0});
        ans_v.push_back(dvec_t {0,1,0});
        ans_v.push_back(dvec_t {1,0,0});
        ans_v.push_back(dvec_t {0,-1,0});
        ans_v.push_back(dvec_t {-1,0,0});

        CHECK_ARRAY2D_CLOSE(ans_v, test_v, 5, 3, tolerance);
    }

    TEST_FIXTURE(BasisFixture, negPi2Rotate) {
        Tools::rotate(test_v, -M_PI/2, axis);

        vector<dvec_t> ans_v;
        ans_v.push_back(dvec_t {0,0,0});
        ans_v.push_back(dvec_t {-1,0,0});
        ans_v.push_back(dvec_t {0,1,0});
        ans_v.push_back(dvec_t {1,0,0});
        ans_v.push_back(dvec_t {0,-1,0});

        CHECK_ARRAY2D_CLOSE(ans_v, test_v, 5, 3, tolerance);
    }
}

int main(int, const char *[]) {
   return UnitTest::RunAllTests();
}
