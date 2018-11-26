#include "Linterp.h"

Linterp::Linterp(float *X, float *Y, int N) {

	n = N;
	x = (float*) malloc(n*sizeof(float));
	y = (float*) malloc(n*sizeof(float));
	int i;
	float p, pu, psig, *u;
	for (i=0;i<n;i++) {
		x[i] = X[i];
		y[i] = Y[i];
	}
	yg = (float*) malloc((n-1)*sizeof(float));
	

	
	/* calculate middle derivatives*/
	for (i=1;i<n-1;i++) {
		yg[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
	}


}

Linterp::~Linterp() {
	free(yg);
	free(x);
	free(y);
}

float Linterp::Interp(float xi) { 
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

	return yg[p]*(xi-x[p]) + y[p];

}
