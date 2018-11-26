#include "rkg.h"

/* the constants here have been declared in the hopoe that they will help optimize the code*/
const float root2 = sqrt(2.0);
const float _c1 = 0.5*(root2-1.0);
const float _c2 = (1.0 - 0.5*root2);
const float _c3 = 0.5*root2;
const float _c4 = (1.0 + 0.5*root2);
const float _c5 = (2.0 - root2);
const float _c6 = (2.0 + root2);
const float _c7 = 1.0/6.0;

void RKG(float x, float *y, float h, void (*f)(float*,float*), float *yout) {
	//Runge-Kutta-Gill method, here f is a function
	float tmp[2], k1[2], k2[2], k3[2], k4[2];
	
	(*f)(y,k1);
	k1[0] = h*k1[0];
	k1[1] = h*k1[1];
	tmp[0] = y[0] + 0.5*k1[0];
	tmp[1] = y[1] + 0.5*k1[1];
	(*f)(tmp,k2);
	k2[0] = h*k2[0];
	k2[1] = h*k2[1];	
	tmp[0] = y[0] + _c1*k1[0] + _c2*k2[0];
	tmp[1] = y[1] + _c1*k1[1] + _c2*k2[1];
	(*f)(tmp,k3);
	k3[0] = h*k3[0];
	k3[1] = h*k3[1];	
	tmp[0] = y[0] - _c3*k2[0] + _c4*k3[0];
	tmp[1] = y[1] - _c3*k2[1] + _c4*k3[1];
	(*f)(tmp,k4);
	k4[0] = h*k4[0];
	k4[1] = h*k4[1];
	yout[0] = y[0] + _c7*(k1[0] + _c5*k2[0] + _c6*k3[0] + k4[0]);
	yout[1] = y[1] + _c7*(k1[1] + _c5*k2[1] + _c6*k3[1] + k4[1]);
}
