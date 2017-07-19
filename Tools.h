#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include "define.h"

using namespace std;

vector<doublev_t> dot(vector<doublev_t>, vector<doublev_t>);

vector<doublev_t> dot(vector<doublev_t>, doublev_t);

void printArr(vector<doublev_t>);

vector<doublev_t> rotate(vector<doublev_t>, double, doublev_t);

#endif
