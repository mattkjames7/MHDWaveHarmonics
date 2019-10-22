#include "solvewave.h"


void SolveWave(float f, float *x, int n, float *Va, float *dlndx, float *y) {
	/*
	 * Use the Runge-Kutta mehthod to solve the wave equation along
	 * an arbitrary field line.
	 * 
	 * Inputs:
	 * 		f: frequency to solve for.
	 * 		x: Distance along the field line in km.
	 * 		n: number of elements in x.
	 * 		Va: Array containing the Alfven speed at midpoints along x.
	 * 		dlndx: array containing the 2nd term of the Singer et al equation.
	 * 
	 * Outputs:
	 * 		y: field displacement array.
	 * 
	 * */
	
	float w = 2.0*M_PI*f/1000.0;
	float yt[] = {0.0,1.0}, yo[2];
	int i;
	for (i=0;i<n-1;i++) {
		k = w/Va[i];
		dln = dlndx[i];
		RKG(x[i],yt,x[i+1]-x[i],Diff,yo);
		y[i+1] = yo[0];
		yt[0] = yo[0];
		yt[1] = yo[1];
	}
}



void SolveWaveComplex(float f, float *x, int n, float *Va, float *dlndx, float *yr, float *yi, float *phase, float *mxr, float *mxi) {
	/*
	 * Use the Runge-Kutta mehthod to solve the wave equation along
	 * an arbitrary field line, with a pseudo complex result.
	 * 
	 * Inputs:
	 * 		f: frequency to solve for.
	 * 		x: Distance along the field line in km.
	 * 		n: number of elements in x.
	 * 		Va: Array containing the Alfven speed at midpoints along x.
	 * 		dlndx: array containing the 2nd term of the Singer et al equation.
	 * 
	 * Outputs:
	 * 		yr: field displacement array (real).
	 * 		yi: field displacement array (pseudo imaginary).
	 * 		phase: wave "phase" along the field line.
	 * 		mxr: maximum displacement for the real part.
	 * 		mxi: maximum displacement for the imaginary part.
	 * 
	 * */
	 
	float w = 2.0*M_PI*f/1000.0;
	float yt[] = {0.0,1.0}, yo[2];
	int i;
	yr[0] = yt[0];
	yi[0] = yt[1];
	
	for (i=0;i<n-1;i++) {
		k = w/Va[i];
		dln = dlndx[i];
		RKG(x[i],yt,x[i+1]-x[i],Diff,yo);
		yr[i+1] = yo[0];
		yi[i+1] = yo[1];
		yt[0] = yo[0];
		yt[1] = yo[1];
	}
	
	/* rescale both arrays  and calculate "phase"*/
	mxi[0] = 0.0;
	mxr[0] = 0.0;
	int np = 0;
	for (i=0;i<n;i++) {
		if (fabsf(yr[i]) > mxr[0]) {
			mxr[0] = fabsf(yr[i]);
		}
		if (fabsf(yi[i]) > mxi[0]) {
			mxi[0] = fabsf(yi[i]);
		}
	}
	for (i=0;i<n;i++) {
		yr[i] /= mxr[0];
		yi[i] /= mxi[0];
		phase[i] = atan2(yr[i],yi[i])*180.0/M_PI + np*360.0;
		if (i > 0) {
			if (phase[i] < phase[i-1] - 5) { //(what the) FUDGE!!!!!!!!!!!!!!!!
				phase[i] += 360.0;
				np++;
			}
		}
	}
	
}


void SolveWaveCC(float f, float *x, float *B0, int n, float *Va, float *dlndx, float *yr, float *yi, float *xi, float *phase, float *mxr, float *mxi, float *b, float *E) {
	/*******************************************************************
	 * Try solving a wave for the situation where the wave is closed at
	 * both ends (node at each end)
	 * 
	 * In this case:
	 * 		xi = 0 at each end
	 * 		E = 0 at each end
	 * ****************************************************************/
	
	float w = 2.0*M_PI*f/1000.0;
	float yt[] = {0.0,1.0}, yo[2];
	int i;
	b[0] = 0.0;
	bt[0] = 0.0;
	bt[1] = 1.0;
	yi[0] = b[i]/(h[i]*B0[i]);
	
	
	yr[0] = yt[0];
	yi[0] = yt[1];
	
	for (i=0;i<n-1;i++) {
		k = w/Va[i];
		dln = dlndx[i];
		RKG(x[i],yt,x[i+1]-x[i],Diff,yo);
		yr[i+1] = yo[0];
		yi[i+1] = yo[1];
		yt[0] = yo[0];
		yt[1] = yo[1];
	}

	/*calculate xi,E,b*/
	for (i=0;i<n;i++) {
		xi[i] = yr[i]*h[i]
		E[i] = -w*xi[i]*B0[i];
		b[i] = h[i]*yi[i]*B0[i];
	}

	
	/* rescale both arrays  and calculate "phase"*/
	mxi[0] = 0.0;
	mxr[0] = 0.0;
	int np = 0;
	for (i=0;i<n;i++) {
		if (fabsf(yr[i]) > mxr[0]) {
			mxr[0] = fabsf(yr[i]);
		}
		if (fabsf(yi[i]) > mxi[0]) {
			mxi[0] = fabsf(yi[i]);
		}
	}
	for (i=0;i<n;i++) {
		yr[i] /= mxr[0];
		yi[i] /= mxi[0];
		phase[i] = atan2(yr[i],yi[i])*180.0/M_PI + np*360.0;
		if (i > 0) {
			if (phase[i] < phase[i-1] - 5) { //(what the) FUDGE!!!!!!!!!!!!!!!!
				phase[i] += 360.0;
				np++;
			}
		}
	}
	
}


void SolveWaveOO(float f, float *x, float *B0, int n, float *Va, float *dlndx, float *yr, float *yi, float *xi, float *phase, float *mxr, float *mxi, float *b, float *E) {
	/*******************************************************************
	 * Try solving a wave for the situation where the wave is open at
	 * both ends (antinode at each end)
	 * 
	 * In this case:
	 * 		b = 0 at each end
	 * ****************************************************************/
	
	float w = 2.0*M_PI*f/1000.0;
	float yt[] = {0.0,1.0}, yo[2];
	int i;
	yr[0] = yt[0];
	yi[0] = yt[1];
	
	for (i=0;i<n-1;i++) {
		k = w/Va[i];
		dln = dlndx[i];
		RKG(x[i],yt,x[i+1]-x[i],Diff,yo);
		yr[i+1] = yo[0];
		yi[i+1] = yo[1];
		yt[0] = yo[0];
		yt[1] = yo[1];
	}

	/*calculate xi,E,b*/
	for (i=0;i<n;i++) {
		xi[i] = yr[i]*h[i]
		E[i] = -w*xi[i]*B0[i];
		b[i] = h[i]*yi[i]*B0[i];
	}

	
	/* rescale both arrays  and calculate "phase"*/
	mxi[0] = 0.0;
	mxr[0] = 0.0;
	int np = 0;
	for (i=0;i<n;i++) {
		if (fabsf(yr[i]) > mxr[0]) {
			mxr[0] = fabsf(yr[i]);
		}
		if (fabsf(yi[i]) > mxi[0]) {
			mxi[0] = fabsf(yi[i]);
		}
	}
	for (i=0;i<n;i++) {
		yr[i] /= mxr[0];
		yi[i] /= mxi[0];
		phase[i] = atan2(yr[i],yi[i])*180.0/M_PI + np*360.0;
		if (i > 0) {
			if (phase[i] < phase[i-1] - 5) { //(what the) FUDGE!!!!!!!!!!!!!!!!
				phase[i] += 360.0;
				np++;
			}
		}
	}
	
}
