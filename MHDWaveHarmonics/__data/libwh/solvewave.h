#ifndef __solvewave_h__
#define __solvewave_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libwh.h"
#include "rkg.h"
#include "spline.h"
using namespace std;
void SolveWave(float f, float *x, int n, float *Va, float *dlndx, float *y); //keep this one
void SolveWaveComplex(float f, float *x, int n, float *Va, float *dlndx, float *yr, float *yi, float *phase, float *mxr, float *mxi);

#endif 
