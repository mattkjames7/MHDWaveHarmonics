#ifndef __misfitobject_h__
#define __misfitobject_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libwh.h"
#include "findharmonics.h"

using namespace std;

class MisfitObject {
	public: 
		MisfitObject(float *_B, float *_R, float *_s, float *_halpha, float *_InPlanet, float *_RhoBG, int _n, int _ntrace, int *_nsteps, float *_maxR, float *_freqs, int *_harms, float _df);
		float GetMisfit(float*,int,bool);
	private:
		float *B, *R, *s, *halpha, *InPlanet, *RhoBG, *maxR, df, *freqs;
		int n, *harms, ntrace, *nsteps;
};
#endif
