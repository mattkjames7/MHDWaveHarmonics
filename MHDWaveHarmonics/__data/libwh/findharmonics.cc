#include "findharmonics.h"


void FindHarmonics(float *B, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float df, int *HarmInds, int nh, float x0, int *nIter, float *freqs) {
	if (df <= 0.0) {
		df = 0.1;
	}
	/*	This Function should find the harmonic frequencies of a given field line using
	 * 	Sam's method.
	 *	
	 * Inputs:
	 * 		B: array containing field magnitude in nT
	 * 		R: array containing the radial distance from centre of the planet (in Rp).
	 * 		s: array which contains the distance along the field line in km
	 * 		halpha: See Singer et al 1981 - based on separation of two field lines
	 * 		InPlanet: a floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition.
	 *		n: length of the above arrays.
	 * 		Params: input parameter array - for the simple power law model, this array contains
	 * 			2 parameters: [p_eq,power]
	 * 			where mass density is given by: rho = p_eq*(Rmax/R)**power
	 * 			if provided with a 5 element array containing [n0,alpha,a,beta,mav0], 
	 * 			then a variant of the Sandhu model is used:
	 * 			where n0, alpha and a are used to calculate number density along the field line
	 * 			R/Rmax <= 0.8 : n = n0*(R/Rmax)**-alpha
	 * 			R/Rmax  > 0.8 : n = a*exp(-0.5*(((R/Rmax)-1.0)/0.1)**2) + n0
	 * 			(this is simplified to form a single equation to make the function continuuous) 
	 * 			and beta and mav0 are used to define the power law of the average mass density 
	 * 			along the field line: mav = mav0*(R/Rmax)**-beta
	 * 		nP: Integer number of parameters given in Params
	 * 		maxR: Maximum radial distance along the field line (Rp).
	 * 		df: frequency increment to make each time a wave is shot (mHz).
	 * 		HarmInds: integer array of harmonic numbers, where 1 is the fundamental
	 * 		nh: the nuber of elements of the HarmInds array
	 * 		x0: starting frequency in mHz.
	 * 		
	 * Outputs:
	 * 		nIter: integer (passed as a pointer) which will return the number of iterations
	 * 			undergone in total.
	 * 		freqs: float array (of size nh) which will contain the output frequencies in mHz.
	 * 
	 */
	int i, h = 0,CurrH = 0, nmax = 1000;
	float Va[n-1], dlndx[n-1], yarr[nmax], farr[nmax], fc = df, y[n], f0 = x0;
	
	Spline *tmpspline;
	bool CrossZero;
	
	for (i=0;i<nh;i++) {
		freqs[i] = NAN;
	//	printf("%d\n",HarmInds[i]);
	}
	
	/* populate Va array */
	CalcFieldLineVaMid(B,R,s,halpha,InPlanet,RhoBG,n,Params,nP,maxR,Va);
	
	
	/* work out dlndx array */
	Calcdlndx(B,s,halpha,InPlanet,n,dlndx);
	
		
	*nIter = 0;
	/*shoot some waves*/
	while (h < nh) {
		fc = f0;
		for (i=0;i<nmax;i++) {
			/*solve wave*/
			SolveWave(fc,s,n,Va,dlndx,y);	
			yarr[i] = y[n-1];
			farr[i] = fc;
			
			if (i > 0) {
				CrossZero = (((yarr[i] > 0) && (yarr[i-1] <= 0)) ||  ((yarr[i] < 0) && (yarr[i-1] > 0)));
			} else {
				CrossZero = false;
			}
			if (CrossZero) {
				CurrH++;
				if (CurrH == HarmInds[h]) {
					
					if (i > 2) {
						/* we can use a spline to interpolate here*/
						tmpspline = new Spline(&yarr[i-3],&farr[i-3],4);
						freqs[h] = tmpspline->Interp(0.0);
						h++;
						delete tmpspline;
					} else if (CrossZero) {
						/* use linear interpolation here */
						freqs[h] = farr[i] - yarr[i]*((farr[i]-farr[i-1])/(yarr[i]-yarr[i-1]));
						h++;
					}
				}
				if (h == nh) {
					break;
				}
			}
			fc += df;
			(*nIter)++;
			if (*nIter > 2000) {
				freqs[h] = NAN;
				break;
			}			
			if (i == nmax-1) {
				//if we have reached nmax, then all of the frequencies output after this are likely to be wrong
				return;
			}	
		}
		f0 += farr[nmax-2];

	}
	
	
}



void FindHarmonicsComplex(float *B, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, int *HarmInds, int nh, float *x0, bool *Success, int *nIter, float *freqs) {
	/*	This Function should find the harmonic frequencies of a given field line using
	 * 	the "complex" wave solution, by minimimzing the difference between the output wave
	 * 	"phase" with n*pi, where n is the requested harmonic number. It's worth noting
	 * 	that the "phase" is not a true phase, just some magic based on the gradients
	 * 	of the wave solution - the important thing it that it is continuous and it would
	 * 	be exactly equal to the phase every pi/2. The actual optimization is done using the 
	 * 	Nelder-Mead algorithm (downhill simplex).
	 *	
	 * Inputs:
	 * 		B: array containing field magnitude in nT
	 * 		R: array containing the radial distance from centre of the planet (in Rp).
	 * 		s: array which contains the distance along the field line in km
	 * 		halpha: See Singer et al 1981 - based on separation of two field lines
	 * 		InPlanet: a floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition.
	 *		n: length of the above arrays.
	 * 		Params: input parameter array - for the simple power law model, this array contains
	 * 			2 parameters: [p_eq,power]
	 * 			where mass density is given by: rho = p_eq*(Rmax/R)**power
	 * 			if provided with a 5 element array containing [n0,alpha,a,beta,mav0], 
	 * 			then a variant of the Sandhu model is used:
	 * 			where n0, alpha and a are used to calculate number density along the field line
	 * 			R/Rmax <= 0.8 : n = n0*(R/Rmax)**-alpha
	 * 			R/Rmax  > 0.8 : n = a*exp(-0.5*(((R/Rmax)-1.0)/0.1)**2) + n0
	 * 			(this is simplified to form a single equation to make the function continuuous) 
	 * 			and beta and mav0 are used to define the power law of the average mass density 
	 * 			along the field line: mav = mav0*(R/Rmax)**-beta
	 * 		nP: Integer number of parameters given in Params
	 * 		maxR: Maximum radial distance along the field line (Rp).
	 * 		HarmInds: integer array of harmonic numbers, where 1 is the fundamental
	 * 		nh: the nuber of elements of the HarmInds array
	 * 		x0: starting frequency in mHz.
	 * 		
	 * Outputs:
	 * 		Success: boolean array denoting whether the optimization was actually successful.
	 * 		nIter: integer array which will return the number of iterations
	 * 			undergone for each frequency.
	 * 		freqs: float array (of size nh) which will contain the output frequencies in mHz.
	 * 
	 */
	float Va[n-1], dlndx[n-1], yr[n], yi[n], phase[n], Rnorm, mav, ne, HarmPhase[nh], CurrTargetPhase;
	int i;

	for (i=0;i<nh;i++) {
		freqs[i] = NAN;
		HarmPhase[i] = 180.0*HarmInds[i];
	}

	/* populate Va array */
	CalcFieldLineVaMid(B,R,s,halpha,InPlanet,RhoBG,n,Params,nP,maxR,Va);
	

	/* work out dlndx array */
	Calcdlndx(B,s,halpha,InPlanet,n,dlndx);

	/* Create Lambda function which outputs the absolute difference in phase
	 * between the target harmonic (harm*180) and the output phase of the model */
	auto DiffLambda = [&CurrTargetPhase,s,n,&Va,&dlndx,&yr,&yi,&phase](float* f, int I) mutable {
		float mr, mi;
		SolveWaveComplex(f[0],s,n,Va,dlndx,yr,yi,phase,&mr,&mi);
		return fabsf(phase[n-1] - CurrTargetPhase);
	};

	float fout[1];
	float mxr,mxi;

	
	/*shoot some waves*/
	for (i=0;i<nh;i++) {
		
		CurrTargetPhase = HarmPhase[i];
		Success[i] = neldermead(&x0[i],1,DiffLambda,2000,1e-6,1e-6,true,&nIter[i],fout);
		freqs[i] = fout[0];
		
		SolveWaveComplex(fout[0],s,n,Va,dlndx,yr,yi,phase,&mxr,&mxi);

	}
	
}

