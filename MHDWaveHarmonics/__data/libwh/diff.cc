#include "diff.h"
/* In order for this to reduce to the simple wave equation, set dln = 0.0
 * 
 * dln is the gradient of the logarithmic term of the singer equation
 * 
 * also: dln and k are global values
 * */
void Diff(float *y, float *dydx) {
	dydx[0] = y[1];
	dydx[1] = -y[1]*dln - pow(k,2.0)*y[0];
}

