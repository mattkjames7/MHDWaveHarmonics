#ifndef __findharmonics_h__
#define __findharmonics_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libwh.h"
#include "spline.h"
#include "solvewave.h"
#include "neldermead.h"
using namespace std;

void ShootWaves(	int n, float *Va, float *dlndx, float *s, 
					float df, float x0,
					int nh, int *HarmInds,
					int *nIter, float *freqs);

void ShootWavesComplex(	int n, float *Va, float *dlndx, float *s,
						int nh, float *HarmPhase, int *HarmInds, float *x0, 
						bool *Success, int *nIter, float *freqs);
	
extern "C" {			

	void FindHarmonics(	float *B, float *R, float *s, float *halpha, 
						float *InPlanet, float *RhoBG, int n, 
						float *Params, int nP, float maxR, float df, 
						int *HarmInds, int nh, float x0, 
						int *nIter, float *freqs);

	void FindHarmonicsComplex(	float *B, float *R, float *s, float *halpha, 
								float *InPlanet, float *RhoBG, int n, 
								float *Params, int nP, float maxR, 
								int *HarmInds, int nh, float *x0, 
								bool *Success, int *nIter, float *freqs);
								
	void FindHarmonicsPMD(	float *B, float *s, float *halpha, float *pmd,
							float *InPlanet, float *RhoBG, int n, 
							float df, 
							int *HarmInds, int nh, float x0, 
							int *nIter, float *freqs);
							
	void FindHarmonicsPMDComplex(	float *B, float *s, float *halpha, float *pmd,
									float *InPlanet, float *RhoBG, int n, 
									int *HarmInds, int nh, float *x0, 
									bool *Success, int *nIter, float *freqs);
}
#endif
