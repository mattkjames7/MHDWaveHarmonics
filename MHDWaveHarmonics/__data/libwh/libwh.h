#ifndef __libwh_h__
#define __libwh_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
extern float k;
extern float dln;
extern float RefIndex;
const float mu0 = 4*M_PI*1.0e-7;
const float amu = 1.67377e-27;
#include "diff.h"
#include "rkg.h"
#include "solvewave.h"
#include "spline.h"
#include "Linterp.h"
#include "findharmonics.h"
#include "misfitobject.h"
using namespace std;



extern "C" {
	void SolveWaveWrapper(float f, float *B, float *R,  float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float *y);
	void SolveWaveWrapperVa(float f, float *B, float *Va,  float *s, float *halpha, float *InPlanet, int n, float *y);
	void SolveWaveComplexWrapper(float f, float *B, float *R,  float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float *yr, float *yi, float *phase, float *mxr, float *mxi);
	void FindHarmonicsComplexWrapper(float *Bm, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, int *HarmInds, int nh, float *x0, bool *Success, int *nIter, float *freqs);
	void FindHarmonicsWrapper(float *Bm, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float df, int *HarmInds, int nh, float x0, int *nIter, float *freqs);
	void GridMisfitWrapper(float *B, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, int ntrace, int *nsteps, float *maxR, float *freqs, int *harms, float *Params, int nP, float df, bool Complex, int nG, float *grid);
	void SetRefractiveIndex(float NewRefIndex);
	float GetRefractiveIndex();
	void CalcFieldLineVa(float *B, float *R,  float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float *Va);
	void CalcFieldLineVaMid(float *B, float *R,  float *s, float *halpha, float *InPlanet, float *RhoBG, int n, float *Params, int nP, float maxR, float *Vamid);
	void Calcdlndx(float *B, float *s, float *halpha, float *InPlanet, int n, float *dlndx);
	int NextInstance();
	int MaxInstance();
	int InstanceIndex(int ins);
	int InitMisfitObject(float *B, float *R, float *s, float *halpha, float *InPlanet, float *RhoBG, int n, int ntrace, int *nsteps, float *maxR, float *freqs, int *harms, float df);
	float GetMisfit(int instance, float *Params, int nP, bool Complex);
	void DestroyMisfitObject(int instance);
}

#endif
