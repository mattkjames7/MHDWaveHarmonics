#ifndef __rkg_h__
#define __rkg_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;
void RKG(float x, float *y, float h, void (*f)(float*,float*), float *yout);
#endif
