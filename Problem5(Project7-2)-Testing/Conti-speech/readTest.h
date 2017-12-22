#ifndef _READTEST_H

#define _CRT_SECURE_NO_WARNINGS

#include"mfcc.h"
#include"readwave.h"
#include<string>
#include<iostream>

using namespace std;

short* readAudioData(char *fileName, int *numSamples, int *sampleRate);
string extractFeatureOfAudio(short *waveData, int nbSamples, int sampleRate);

#endif