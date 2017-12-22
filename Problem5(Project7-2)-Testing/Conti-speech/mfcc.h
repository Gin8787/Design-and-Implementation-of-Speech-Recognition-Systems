#ifndef _MFCC_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "fftw3.h"

using namespace std;

#define FRAME_LENGTH (200)
#define OVERLAP (100)
#define FILTERS (40)
#define MAXFREQ (4000)
//#define MAXFREQ (4000)
#define MFCC_COEF (13)

typedef struct FilterMelFreq{
	double lower_bound;
	double center;
	double higher_bound;
}FilterMelFreq;

void feature_extraction(short* waveData, int numSamples, int sampleRate, string normalFileName);

//void computeMel(float *mel, int sampleRate, const float *energySpectrum);
//void DCT(const float *mel, float *c);

#endif