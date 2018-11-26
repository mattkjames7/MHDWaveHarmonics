#include "misfitobject.h"


MisfitObject::MisfitObject(float *_B, float *_R, float *_s, float *_halpha, float *_InPlanet, int _n, int _ntrace, int *_nsteps, float *_maxR, float *_freqs, int *_harms, float _df) {
	/*
	 * The constructor of the MisfitObject, passes pointers to the original data tot he object.
	 * 
	 * Inputs:
	 * 		_B: Flattened array of magnetic field strength, length n*ntrace.
	 * 		_R: Flattened array of radial distance from the centre of the
	 * 			planetary dipole, length n*ntrace.
	 * 		_s: Flattened array of distance along each field line, length
	 * 			n*ntrace.
	 * 		_halpha: Flattened array of halpha values, length n*ntrace.
	 * 		_InPlanet: a flattened floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition, length n*ntrace.
	 * 		_n: maximum number of steps for the traces.
	 * 		_ntrace: number of field line traces.
	 * 		_nsteps: array containing the number of steps for each trace.
	 * 		_maxR: array containing the maximum radial distances for each trace.
	 * 		_freqs: array of frequencies to fit to for each trace.
	 * 		_harms: array containing the harmonic numbers of the freqs array.
	 * 		_df: frequency increment to make each time a wave is shot (mHz) - this
	 * 			is only used if the Complex variable is set to false.
	 */
	n = _n; //total length of arrays
	ntrace = _ntrace; //number of field line traces
	int i;
	
	R = _R;
	B = _B;
	s = _s;
	halpha = _halpha;
	InPlanet = _InPlanet;
	
	maxR = _maxR;
	nsteps = _nsteps;
	

	df = _df;
	freqs = _freqs;
	harms = _harms;
}


float MisfitObject::GetMisfit(float *Params, int nP, bool Complex) {
	/*
	 * Class member function of the MisfitObject to calculate the 
	 * misfit given a set of input parameters.
	 * 
	 * Inputs:
	 * 		Params: array of parameters (2 or 5 for each trace).
	 * 		nP: number of parameters (2 or 5).	
	 * 		Complex: Boolean value, if true the complex method of shooting waves 
	 * 			is used, otherwise Sam's method is used. 
	 * Output:
	 * 		measure of frequency misfit.
	 * */
	
	float misfit = 0.0, ftmp[1], f0c[1] = {1.0f};
	int i,p,nIter[1];
	bool Success[1];
	
	//printf("Params: %f %f\n",Params[0],Params[1]);
	p = 0;
	for (i=0;i<ntrace;i++) {
		if (Complex) {
			FindHarmonicsComplex(&B[p],&R[p],&s[p],&halpha[p],&InPlanet[p],nsteps[i],Params,nP,maxR[i],&harms[i],1,f0c,Success,nIter,ftmp);
		} else {
			FindHarmonics(&B[p],&R[p],&s[p],&halpha[p],&InPlanet[p],nsteps[i],Params,nP,maxR[i],df,&harms[i],1,df,nIter,ftmp);
			//printf("%d %f %f %f\n",nIter[0],freqs[i],ftmp[0],df);
		}
		misfit += pow(freqs[i]-ftmp[0],2.0);
		p+=(n/ntrace);
	}		
	return sqrt(misfit/ntrace);
}
