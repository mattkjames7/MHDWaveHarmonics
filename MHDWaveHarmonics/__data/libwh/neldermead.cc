#include "neldermead.h"


float _fError(float *fvals, int n, int *inds) {
	/* Get error for f*/
	float dfmax = 0.0, dftmp;
	int i;
	for (i=1;i<=n;i++) {
		dftmp = fabsf(fvals[inds[0]] - fvals[inds[i]]);
		if (dftmp > dfmax) {
			dfmax = dftmp;
		}
	}
	return dfmax;
}

float _xError(float **sim, int n, int *inds) {
	/* Get error for simplex in parameter space*/
	float dxmax = 0.0, dxtmp;
	int i,j;
	for (i=1;i<=n;i++) {
		for (j=0;j<n;j++) {
			dxtmp = fabsf(sim[inds[0]][j] - sim[inds[i]][j]);
			if (dxtmp > dxmax) {
				dxmax = dxtmp;
			}
		}
	}
	return dxmax;
}

void _CalcCentroid(float **sim, int n, int *inds, float *simc) {
	/* calculate the centre of all but the last simplex in parameter space */
	int i,j;
	float tmp;
	for (i=0;i<n;i++) {
		tmp = 0.0;
		for (j=0;j<n;j++) {
			tmp+= sim[inds[j]][i];
		}
		simc[i] = tmp/n;
	}
}

// can't really use float (*f)(float*,int) as the third argument here

bool neldermead(float *x0, int n, std::function<float(float*,int)> F, int MaxIter, float fatol, float xatol, bool Adapt, int *nIter, float *xout) {
	/*	If all goes to plan, this routine will take the initial guess x0
	 * and use the downhill simplex method to minimize the output of (*f).
	 *  
	 * Inputs:
	 * 		*x0: vector of initial parameter guesses, n elements long
	 * 		n: number of parameters in x0
	 * 		F: Lambda function which accepts a float array and an integer, returns float
	 * 		MaxIter: maximum number of iterations to attempt
	 * 		fatol: Absolute error in f(x) 
	 * 		xatol: Absolute error in x
	 * 		*nIter: total number of iterations performed
	 * 		*xout: return array 
	 * 
	 * Returns:
	 * 		bool indicating successful minimization or not
	 */
	
	/* define scaling parameters*/
	float alpha, gamma, rho, sigma; //reflection, expansion, contraction and shrink coeffs
	if (Adapt) {
		alpha = 1.0;
		gamma = 1.0 + 2.0/n;
		rho = 0.75 - 1.0/(2.0*n);
		sigma = 1.0 - 1.0/n;
	} else {
		alpha = 1.0;
		gamma = 2.0;
		rho = 0.5;
		sigma = 0.5;
	}		

	//Create initial simplex array (n+1 simplexes)
	float **sim = new float*[n+1];
	float sim0[n], simr[n], sime[n], simc[n];
	int i, j;
	for (i=0;i<=n;i++) {
		sim[i] = new float[n];
	}

	
	for (i=0;i<=n;i++) {
		for (j=0;j<n;j++) {
			if (j+1 == i) {
				if (x0[j] == 0.0) {
					sim[i][j] = 0.00025;
				} else {
					sim[i][j] = 1.05*x0[j];
				}
			} else {
				sim[i][j] = x0[j];
			}	 
		}
	}
	

	 
	//create array to store the results from *f
	float fvals[n+1], fr, fe, fc;
	int inds[n+1];
	for (i=0;i<n+1;i++) {
		inds[i] = i;
		fvals[i] = F(sim[i],n);
	}
	ArgSortFlt(n+1,fvals,inds);

	float fE = _fError(fvals,n,inds);

	float xE = _xError(sim,n,inds);
	bool success;

	 
	
	/* Start the loop!
	 * Here we will stop if one of two things happen:
	 * 1. We run out of iterations (i.e. nIter == MaxIter)
	 * 2. The maximum distance in parameter space between each vertex of the simplex is <= xatol
	 * and The maximum difference between the fvals for each vertex <= fatol
	 * */
	*nIter = 0;

	while ((*nIter < MaxIter) && ((fE > fatol) || (xE > xatol))) {
		/* Sort simplex vertices by their value in f */
		ArgSortFlt(n+1,fvals,inds);
		
		/* calculate centroid of all simplexes except the final one*/
		_CalcCentroid(sim,n,inds,sim0);

		/* do reflection*/
		for (i=0;i<n;i++) {
			simr[i] = sim0[i] + alpha*(sim0[i] - sim[inds[n]][i]);
		}
		fr = F(simr,n);

		if ((fr >= fvals[inds[0]]) & (fr < fvals[inds[n-1]])) {
			/* if reflected point is better than the 2nd worst point,
			 * but not better than the best point - then replace 
			 * and continue*/

			for (i=0;i<n;i++) {
				sim[inds[n]][i] = simr[i];
			}
			fvals[inds[n]] = fr;
		} else if (fr < fvals[inds[0]]) {
			/* if reflected point is the best so far then we can expand 
			 * the simplex*/

			for (i=0;i<n;i++) {
				sime[i] = sim0[i] + gamma*(simr[i] - sim0[i]);
			}
			fe = F(sime,n);
			if (fe < fr) {
				/* if expanded point is better than reflected - then use it*/

				for (i=0;i<n;i++) {
					sim[inds[n]][i] = sime[i];
				}
				fvals[inds[n]] = fe;
			} else {
				/* otherwise continue with reflected point */

				for (i=0;i<n;i++) {
					sim[inds[n]][i] = simr[i];
				}
				fvals[inds[n]] = fr;		
			}		
		} else if (fr < fvals[inds[n]]) {
			/* the reflected point is worse than the 2nd last then we try outside contraction */

			for (i=0;i<n;i++) {
				simc[i] = sim0[i] + rho*(simr[i] - sim0[i]);
			}
			fc = F(simc,n);
			if (fc <= fr) {
				/* replace worst point with this contracted point */

				for (i=0;i<n;i++) {
					sim[inds[n]][i] = simc[i];
				}
				fvals[inds[n]] = fc;
			} else {
				/* shrink if all else fails */

				for (i=0;i<=n;i++) {
					for (j=0;j<n;j++) {
						sim[i][j] = sim[inds[0]][j] + sigma*(sim[i][j] - sim[inds[0]][j]);
					}
				}
				for (i=0;i<=n;i++) {
					fvals[i] = F(&sim[i][0],n);
				}
			}				
			
		} else {
			/* try inside contraction when fr is the worst score */
			
			
			for (i=0;i<n;i++) {
				simc[i] = sim0[i] + rho*(sim[inds[n]][i] - sim0[i]);
			}
			fc = F(simc,n);
			if (fc < fvals[inds[n]]) {
				/* replace worst point with this contracted point */

				for (i=0;i<n;i++) {
					sim[inds[n]][i] = simc[i];
				}
				fvals[inds[n]] = fc;
			} else {
				/* shrink if all else fails */

				for (i=0;i<=n;i++) {
					for (j=0;j<n;j++) {
						sim[i][j] = sim[inds[0]][j] + sigma*(sim[i][j] - sim[inds[0]][j]);
					}
				}
				for (i=0;i<=n;i++) {
					fvals[i] = F(&sim[i][0],n);
				}
			}						
		}

		fE = _fError(fvals,n,inds);
		xE = _xError(sim,n,inds);
		(*nIter)++;
	}

	
	/*once optimization is finished - return the simplex with the best fval*/
	ArgSortFlt(n+1,fvals,inds);

	for (i=0;i<n;i++) {
		xout[i] = sim[inds[0]][i];
	}

	
	/* free memory*/
	for(i=0;i<=n;i++) {
		delete sim[i];
	}
	delete[] sim;
	
	
	
	if ((fE < fatol) || (xE < xatol)) {
		success = true;
	} else {
		success = false;
	}
	return success;
}
