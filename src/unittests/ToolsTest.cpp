#include "UnitTest++/UnitTest++.h"
#include <vector>
#include <iostream>
#include "define.h"
#include "Tools.h"

using namespace std;

SUITE(MathTests) {
    TEST(addToVector) {

        cout << "Testing addToVector()" << endl;

        dvec_t test;
        dvec_t temp;
        dvec_t v1 = {0,0,0};

        test = {0,0,0};
        CHECK_ARRAY_EQUAL(test,v1,3);

        temp = {1,1,1};
        test = {1,1,1};
        Tools::addVectors(v1, temp);
        CHECK_ARRAY_EQUAL(test,v1,3);
    }
}

int main(int, const char *[])
{
   return UnitTest::RunAllTests();
}
