#ifndef __linterp_h__
#define __linterp_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

/* A simple linear interpolation object, where the original data is
 * input on creation, then the Interp function calculates a y value
 * when given an x value*/

class Linterp {
	public:
		Linterp(float*,float*,int);
		~Linterp();
		float Interp(float);
	private:
		float *yg;
		float *x, *y;
		int n;
};

#endif
