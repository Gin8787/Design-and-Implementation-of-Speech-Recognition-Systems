#ifndef _MFCC_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "fftw3.h"

using namespace std;

#define FRAME_LENGTH (400)
#define OVERLAP (200)
#define FILTERS (40)
#define MAXFREQ (8000)
#define MFCC_COEF (13)

typedef struct FilterMelFreq{
	double lower_bound;
	double center;
	double higher_bound;
}FilterMelFreq;

void feature_extraction(short* waveData, int numSamples, int sampleRate);

//void computeMel(float *mel, int sampleRate, const float *energySpectrum);
//void DCT(const float *mel, float *c);

#endif