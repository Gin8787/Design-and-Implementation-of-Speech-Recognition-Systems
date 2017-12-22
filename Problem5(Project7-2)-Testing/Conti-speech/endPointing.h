#ifndef _ENDPOINTING_H_

#include <stdlib.h>
#include <cmath>
#include "readwave.h"
#include<vector>
#include<iostream>
//#include "mfcc.h"

using namespace std;

#define SAMPLE_RATE  (16000)
#define FRAMES_PER_BUFFER (400)
#define NUM_SECONDS     (60)
#define NUM_CHANNELS    (1) 

//typedef short SAMPLE;

double EnergyPerSampleDecibel(const short * frameStart, int startIndex);
bool classifyFrame(const short * frameStart, int startIndex);
void backgroundInitial(const short * firstFrame);
void backgroundInitBackward(const short *frameStart, int startIndex);
short* cutSil(string audioFile, string newAudioFile, int *numSampless, int *sampleRate);
double EnergyPerSampleDecibelBack(const short * frameEnd, int startIndex);


#endif