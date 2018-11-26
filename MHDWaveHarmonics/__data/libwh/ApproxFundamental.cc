#include "ApproxFundamental.h"

float ApproxFundamental(float p_eq, float Smax) {
	/* This function is a very crude way of estimating the initial 
	 * value close to the fundamental frequency before shooting waves.
	 * Inputs:
	 * 		p_eq: equatorial plasma mass density
	 * 		Smax: Length of field line in km
	 * 
	 * Output:
	 * 		starting frequency (this may not be reliable though!)
	 */ 
	
	
	float L,m,a,b;
	
	L = (Smax + 9675.0)/10492.0;
	a = 0.01446;
	b = -0.0153;
	m = powf((a + b*L),2);
	return sqrt(1.0/(m*p_eq));
}
	
	
