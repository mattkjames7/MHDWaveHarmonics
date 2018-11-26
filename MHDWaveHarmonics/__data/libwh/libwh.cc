#include "libwh.h"
float k = 0.0;
float dln = 0.0;
float RefIndex = 1.55; //Refractive index for SiO2 (approx) for use with field lines which penetrate the surface
vector<MisfitObject*> vMisfitObjects;
vector<int> MOInsts;

int NextInstance() {
	if (MOInsts.size() == 0) {
		return 0;
	} else {
		return MaxInstance()+1;
	}
}

int MaxInstance() { 
	int i, n, mx = 0;
	n = MOInsts.size();
	for (i=0;i<n;i++) {
		if (MOInsts[i] > mx) {
			mx = MOInsts[i];
		}
	}
	return mx;
}

int InstanceIndex(int ins) {
	int i, n;
	n = MOInsts.size();
	for (i=0;i<n;i++) {
		if (MOInsts[i] == ins) {
			return i;
		}
	}
	return -1;	
}


void SetRefractiveIndex(float NewRefIndex) {
	/* Call this function to alter the value of the refractive index used
	 * if part of the field line exists in a plasma-free region*/
	RefIndex = NewRefIndex;
}

float GetRefractiveIndex() {
	/*This timply returns the refractive index currently in use*/
	return RefIndex;
}

void CalcFieldLineVa(float *B, float *R,  float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float *Va) {
	/*
	 *	This Function will calculate the Alfven speed along the field line. 
	 * 
	 * Inputs:
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		R: radial distance of every point along the field line in Rp
	 * 		s: distance in km aloong the field line
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
	 * 
	 * Output:
	 * 		Va: Alfven speed along the field line in km/s
	 */
	
	
	int i, j;
	float Rnorm, mav, ne;
	bool HasTransition = false; // transition from in planet to out planet

	/* populate Va array */
	if (nP == 5) {
		for (i=0;i<n;i++) {
			if (InPlanet[i] == 1.0) {
				Va[i] = 3e5/RefIndex;
			} else if (InPlanet[i] == 0.0) {
				Rnorm = (R[i]/maxR);
				mav = Params[4]*pow(Rnorm,-Params[3])*amu; //this is in kg now
				ne = Params[2]*exp(-0.5*pow((Rnorm-1.0)/0.1,2.0)) + Params[0]*pow(Rnorm,-Params[1]);
				ne *= 1e6; //convert to m-3
				Va[i] = ((B[i]*1.0e-9)/sqrt(mu0*ne*mav))/1000.0;	
			} else {
				HasTransition = true;
			}		
		}
	} else if (nP == 2) {
		for (i=0;i<n;i++) {
			if (InPlanet[i] == 1.0) {
				Va[i] = 3e5/RefIndex;
			} else if (InPlanet[i] == 0.0) {
				Rnorm = (maxR/R[i]);
				Va[i] = ((B[i]*1.0e-9)/sqrt(mu0*amu*1.0e6*Params[0]*pow(Rnorm,Params[1])))/1000.0;
			} else {
				HasTransition = true;
			}
		}
	} else {
		printf("Warning! nP = %d - there should be either 2 (for power law) or 5 (for Sandhu model)!!!\n",nP);
		return;
	}
	
	if (HasTransition) {
		/* There may well be parts of the Va array currently uncalculated
		 * these are where 0 < InPlanet < 1 and we must now calculate them */
		int nJ = 0;
		//count them first
		for (i=0;i<n-1;i++) {
			if (((InPlanet[i] == 1.0) | (InPlanet[i] == 0.0)) & ((InPlanet[i+1] < 1.0) & (InPlanet[i+1] > 0.0))) {
				nJ++;
			}
		}
		int *Ji0 = new int[nJ];
		int *Ji1 = new int[nJ];
		int p0 = 0, p1 = 0;
		
		for (i=0;i<n-1;i++) {
			if (((InPlanet[i] == 1.0) | (InPlanet[i] == 0.0)) & ((InPlanet[i+1] < 1.0) & (InPlanet[i+1] > 0.0))) {
				Ji0[p0] = i+1; //inclusive index for start
				p0++;
			}
			if (((InPlanet[i+1] == 1.0) | (InPlanet[i+1] == 0.0)) & ((InPlanet[i] < 1.0) & (InPlanet[i] > 0.0))) {
				Ji1[p1] = i+1; //inclusive index for stop
				p1++;
			}
		}	
		
		//now to calculate Va between them using the surrounding values and a smoothstep function
		float dv,v0;
		for (i=0;i<nJ;i++) {
			dv = fabsf(Va[Ji1[i]+1] - Va[Ji0[i]-1]);
			if (InPlanet[Ji0[i]-1] == 1.0) {
				v0 = Va[Ji1[i]+1];
			} else {
				v0 = Va[Ji0[i]-1];
			}
			printf("j0 %d j1 %d\n",Ji0[i],Ji1[i]);
			for (j=Ji0[i];j<=Ji1[i];j++) {
				Va[j] = v0 + dv*InPlanet[j];
			}
		}
		delete Ji0;
		delete Ji1;
	}
}



void CalcFieldLineVaMid(float *B, float *R,  float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float *Vamid) {
	/*
	 *	This Function will calculate the Alfven speed along the field line at 
	 * 	the midpoints between field line steps. 
	 * 
	 * Inputs:
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		R: radial distance of every point along the field line in Rp
	 * 		s: distance in km aloong the field line
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
	 * 
	 * Output:
	 * 		Vamid: Alfven speed along the field line in km/s
	 */

	float Va[n];
	int i;
	
	CalcFieldLineVa(B,R,s,halpha,InPlanet,n,Params,nP,maxR,Va);
	
	/* create Linterp object*/
	Linterp VaInt(s,Va,n);

	for (i=0;i<n-1;i++) {	
		//Vamid[i] = VaSpline.Interp(0.5*(s[i]+s[i+1]));
		Vamid[i] = VaInt.Interp(0.5*(s[i]+s[i+1]));
	}
}

void Calcdlndx(float *B, float *s, float *halpha, float *InPlanet, int n, float *dlndx) {
	/*
	 *	This Function will calculate the Alfven speed along the field line at 
	 * 	the midpoints between field line steps. 
	 * 
	 * Inputs:
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		s: distance in km aloong the field line
	 * 		halpha: See Singer et al 1981 - based on separation of two field lines
	 * 		InPlanet: a floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition.
	 *		n: length of the above arrays.
	 * 
	 * 	Output:
	 * 		dlndx: The value of the term d/ds(ln(h^2 B))d/ds(xi/h) in the singer et al 
	 * 		equation along the field line.
	 */
	 
	int i;
	if (halpha == NULL) {
		for (i=0;i<n-1;i++) {
			dlndx[i] = 0.0;
		}
	} else {
		for (i=0;i<n-1;i++) {
			dlndx[i] = (1.0-0.5*(InPlanet[i]+InPlanet[i+1]))*(log((pow(halpha[i+1],2)*B[i+1])/(pow(halpha[i],2)*B[i])))/(s[i+1]-s[i]);
		}
	}
}


void SolveWaveWrapper(float f, float *B, float *R,  float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float *y) {
	/*
	 *	This Function will shoot a wave with a given frequency along a field line.
	 * 
	 * Inputs:
	 * 		f: frequency of the wave in mHz.
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		R: radial distance of every point along the field line in Rp
	 * 		s: distance in km aloong the field line
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
	 * 
	 * Output:
	 * 		y: The displacement of the magnetic field line.
	 */

	float Va[n-1];
	int i;
	float dlndx[n-1];
	
	/* populate Va array */
	CalcFieldLineVaMid(B,R,s,halpha,InPlanet,n,Params,nP,maxR,Va);
	
	
	/* work out dlndx array */
	Calcdlndx(B,s,halpha,InPlanet,n,dlndx);
	
	/*solve wave*/
	SolveWave(f,s,n,Va,dlndx,y);
}


void SolveWaveComplexWrapper(float f, float *B, float *R,  float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float *yr, float *yi, float *phase, float *mxr, float *mxi) {
	/*
	 *	This Function will shoot a wave with a given frequency along a field line.
	 * 
	 * Inputs:
	 * 		f: frequency of the wave in mHz.
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		R: radial distance of every point along the field line in Rp
	 * 		s: distance in km aloong the field line
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
	 * 
	 * Output:
	 * 		yr: The displacement of the magnetic field line (real).
	 * 		yi: The gradient along the field line (pseudo-imaginary).
	 * 		phase: The waves phase along the field line.
	 * 		mxr: maximum in the real part.
	 * 		mxi: maximum of imaginary part.
	 */

	float Va[n-1];
	int i;
	float dlndx[n-1];

	/* populate Va array */
	CalcFieldLineVaMid(B,R,s,halpha,InPlanet,n,Params,nP,maxR,Va);
	
	
	/* work out dlndx array */
	Calcdlndx(B,s,halpha,InPlanet,n,dlndx);
	
	/*solve wave*/
	SolveWaveComplex(f,s,n,Va,dlndx,yr,yi,phase,mxr,mxi);
	
}


void SolveWaveWrapperVa(float f, float *B, float *Va,  float *s, float *halpha, float *InPlanet, int n, float *y) {
	/*
	 *	This Function will shoot a wave with a given frequency along a field line given an alfven speed.
	 * 
	 * Inputs:
	 * 		f: frequency of the wave in mHz.
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		Va: Alfven speed.
	 * 		s: distance in km aloong the field line
	 * 		halpha: See Singer et al 1981 - based on separation of two field lines
	 * 		InPlanet: a floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition.
	 *		n: length of the above arrays.
	 * 
	 * Output:
	 * 		y: The displacement of the magnetic field line.
	 */
	int i;
	float dlndx[n-1], Vamid[n-1];
	
	/* create Linterp object*/
	Linterp VaInt(s,Va,n);

	for (i=0;i<n-1;i++) {	
		Vamid[i] = VaInt.Interp(0.5*(s[i]+s[i+1]));
	}
	
	/* work out dlndx array */
	Calcdlndx(B,s,halpha,InPlanet,n,dlndx);
	
	/*solve wave*/
	SolveWave(f,s,n,Vamid,dlndx,y);
}

void FindHarmonicsComplexWrapper(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, int *HarmInds, int nh, float *x0, bool *Success, int *nIter, float *freqs) {
	/*
	 *	This function is a wrapper for the "FindHarmonicsComplex" function which
	 * 	should find the harmonic frequencies of a given field line using
	 * 	the "complex" wave solution, by minimimzing the difference between the output wave
	 * 	"phase" with n*pi, where n is the requested harmonic number. It's worth noting
	 * 	that the "phase" is not a true phase, just some magic based on the gradients
	 * 	of the wave solution - the important thing it that it is continuous and it would
	 * 	be exactly equal to the phase every pi/2. The actual optimization is done using the 
	 * 	Nelder-Mead algorithm (downhill simplex).
	 * 
	 * Inputs:
	 * 		B: array containing the magnetic field magnitude along the field line in nT
	 * 		R: radial distance of every point along the field line in Rp
	 * 		s: distance in km aloong the field line
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
	 */
	FindHarmonicsComplex(B,R,s,halpha,InPlanet,n,Params,nP,maxR,HarmInds,nh,x0,Success,nIter,freqs);
}
void FindHarmonicsWrapper(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, float *Params, int nP, float maxR, float df, int *HarmInds, int nh, float x0, int *nIter, float *freqs) {
	/*	This a wrapper function which should find the harmonic frequencies 
	 *  of a given field line using Sam's method.
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

	FindHarmonics(B,R,s,halpha,InPlanet,n,Params,nP,maxR,df,HarmInds,nh,x0,nIter,freqs);
}

void GridMisfitWrapper(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, int ntrace, int *nsteps, float *maxR, float *freqs, int *harms, float *Params, int nP, float df, bool Complex, int nG, float *grid) {
	/*
	 * This function will calculate the frequency misfit for a grid of 
	 * different parameters.
	 * 
	 * Inputs:
	 * 		B: Flattened array of magnetic field strength, length n*ntrace.
	 * 		R: Flattened array of radial distance from the centre of the
	 * 			planetary dipole, length n*ntrace.
	 * 		s: Flattened array of distance along each field line, length
	 * 			n*ntrace.
	 * 		halpha: Flattened array of halpha values, length n*ntrace.
	 * 		InPlanet: a flattened floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition, length n*ntrace.
	 * 		n: maximum number of steps for the traces.
	 * 		ntrace: number of field line traces.
	 * 		nsteps: array containing the number of steps for each trace.
	 * 		maxR: array containing the maximum radial distances for each trace.
	 * 		freqs: array of frequencies to fit to for each trace.
	 * 		harms: array containing the harmonic numbers of the freqs array.
	 * 		Params: flattened grid of parameters (2 or 5 for each trace).
	 * 		nP: number of parameters (2 or 5).
	 * 		df: frequency increment to make each time a wave is shot (mHz) - this
	 * 			is only used if the Complex variable is set to false.
	 * 		Complex: Boolean value, if true the complex method of shooting waves 
	 * 			is used, otherwise Sam's method is used.
	 * 		nG: number of elements in the output grid.
	 * Output:
	 * 		grid: flattened output grid.
	 * 
	 */
	
	
	
	//Params should be a flattened array of nP*nG (params for each test case in the grid - gridding will be done in Python
	//B,R,s,halpha,InPlanet should each be flattened arrays of n*ntrace
	//tlens, maxR, freqs, harms should be a ntrace length array 
	
	MisfitObject *misfit = new MisfitObject(B,R,s,halpha,InPlanet,n,ntrace,nsteps,maxR,freqs,harms,df);
		
	int i;
	for (i=0;i<nG;i++) {
		grid[i] = misfit->GetMisfit(&Params[i*nP],nP,Complex);
	}
	delete misfit;
}

int InitMisfitObject(float *B, float *R, float *s, float *halpha, float *InPlanet, int n, int ntrace, int *nsteps, float *maxR, float *freqs, int *harms, float df) {
	/*
	 * This function will initialize an object to calculate the frequency 
	 * misfit for different plasma parameters.
	 * 
	 * Inputs:
	 * 		B: Flattened array of magnetic field strength, length n*ntrace.
	 * 		R: Flattened array of radial distance from the centre of the
	 * 			planetary dipole, length n*ntrace.
	 * 		s: Flattened array of distance along each field line, length
	 * 			n*ntrace.
	 * 		halpha: Flattened array of halpha values, length n*ntrace.
	 * 		InPlanet: a flattened floating point arrray containing 1.0 when the field line
	 * 			is within the surface of the planet, or 0.0 outside, values between 0.0
	 * 			and 1.0 may be used to smooth the transition, length n*ntrace.
	 * 		n: maximum number of steps for the traces.
	 * 		ntrace: number of field line traces.
	 * 		nsteps: array containing the number of steps for each trace.
	 * 		maxR: array containing the maximum radial distances for each trace.
	 * 		freqs: array of frequencies to fit to for each trace.
	 * 		harms: array containing the harmonic numbers of the freqs array.
	 * 		df: frequency increment to make each time a wave is shot (mHz) - this
	 * 			is only used if the Complex variable is set to false.
	 * Output:
	 * 		returns an integer corresponding to the index of a global vector 
	 * 		variable where a pointer to this particular object resides.
	 * 
	 */
	
	
	
	//Params should be a flattened array of nP*nG (params for each test case in the grid - gridding will be done in Python
	//B,R,s,halpha,InPlanet should each be flattened arrays of n*ntrace
	//tlens, maxR, freqs, harms should be a ntrace length array 


	int instance, ind;
	instance = NextInstance();
	MOInsts.push_back(instance);
	ind = InstanceIndex(instance);
	
	vMisfitObjects.push_back(new MisfitObject(B,R,s,halpha,InPlanet,n,ntrace,nsteps,maxR,freqs,harms,df));
	return instance;
}

float GetMisfit(int instance, float *Params, int nP, bool Complex) {
	int ind;
	float out;
	ind = InstanceIndex(instance);	
	out = vMisfitObjects[ind]->GetMisfit(Params,nP,Complex);
	return out;
}

void DestroyMisfitObject(int instance) {
	/* this will free up the memory taken up by the misfit object after it is done with*/
	
	int i, ind;
	ind = InstanceIndex(instance);	
	delete vMisfitObjects[ind];
	vMisfitObjects.erase(vMisfitObjects.begin()+ind);
	MOInsts.erase(MOInsts.begin()+ind);

}
