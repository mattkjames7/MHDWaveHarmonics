#include "spline.h"

Spline::Spline(float *X, float *Y, int N) {

	n = N;
	x = (float*) malloc(n*sizeof(float));
	y = (float*) malloc(n*sizeof(float));
	int i;
	float p, pu, psig, *u;
	//check everything it strictly ascending monotonically
	int sgn = 1;
	for (i=0;i<n-1;i++) {
		if (x[i] > x[i+1]) {
			sgn = -1;
		}
	}
	if (sgn == 1) {
		for (i=0;i<n;i++) {
			x[i] = X[i];
			y[i] = Y[i];
		}
	} else {
		for (i=0;i<n;i++) {
			x[i] = X[n-1-i];
			y[i] = Y[n-1-i];
		}
	}
	y2 = (float*) malloc((n)*sizeof(float));
	u = (float*) malloc((n)*sizeof(float));
	
	/* assuming that the first derivative is 0 */
	y2[0] = 0.0;
	u[0] = 0.0;
	
	/* calculate middle second derivatives */
	for (i=1;i<n-1;i++) {
		psig = (x[i] - x[i-1])/(x[i+1]-x[i-1]);
		pu = (((y[i+1]-y[i])/(x[i+1]-x[i])) - ((y[i]-y[i-1])/(x[i]-x[i-1])))/(x[i+1]-x[i-1]);
		p = psig*y2[i-1] + 2;
		y2[i] = (psig - 1)/p;
		u[i] = (6*pu - psig*u[i-1])/p;
	}

	/* also assume the last derivative is 0 */
	y2[n-1] = 0.0;
	u[n-1] = 0.0;
	
	for (i=n-2;i>=0;i--) {
		y2[i] = y2[i]*y2[i+1] + u[i];
	}

	free(u);
}

Spline::~Spline() {
	free(y2);
	free(x);
	free(y);
}

float Spline::Interp(float xi) { 
	int p, i;
	if (xi < x[0]) {
		p = 0;
	} else if (xi >= x[n-1]) { 
		p = n-2;
	} else {
		for (i=0;i<n-1;i++) {
			if ((xi >= x[i]) && (xi < x[i+1])) {
				p = i;
				break;
			}
		}
	}
	float h, a, b;
	h = x[p+1] - x[p];
	a = (x[p+1] - xi)/h;
	b = (xi - x[p])/h;
	
	return a*y[p] + b*y[p+1] + ((pow(a,3) - a)*y2[p] + (pow(b,3) - b)*y2[p+1])*pow(h,2)/6.0;
}
