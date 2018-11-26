#ifndef __findharmonics_h__
#define __findharmonics_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libwh.h"
#include "spline.h"
#include "solvewave.h"
#include "neldermead.h"
using namespace std;

void FindHarmonicsComplex(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, int *HarmInds, int nh, float *x0, bool *Success, int *nIter, float *freqs);
void FindHarmonics(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float df, int *HarmInds, int nh, float x0, int *nIter, float *freqs);

#endif
