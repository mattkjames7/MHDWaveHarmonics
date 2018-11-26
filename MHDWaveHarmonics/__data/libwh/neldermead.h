#ifndef __neldermead_h__
#define __neldermead_h__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ArgSort.h"
#include <functional>
#endif

using namespace std;
bool neldermead(float *x0, int n, std::function<float(float*,int)> F, int MaxIter, float fatol, float xatol, bool Adapt, int *nIter, float *xout);

