#ifndef __spline_h__
#define __spline_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

class Spline {
	public:
		Spline(float*,float*,int);
		~Spline();
		float Interp(float);
	private:
		float *y2;
		float *x, *y;
		int n;
};

#endif
