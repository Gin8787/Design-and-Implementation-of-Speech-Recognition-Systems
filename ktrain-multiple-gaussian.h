#ifndef KTRAIN_H
#define KTRAIN_H

#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>

using namespace std;

#define SEGNUM (5)
#define SAMPLENUM (5)
#define INF (100000)
#define EPSL (1)
#define INIT (1)
#define ITER (2)
#define GAUSNUM (4)

void readSampleDCT(const char* filename, int flag);
void initSeg();
void testPrintSample();
void calSegmentMeansCovar();
void calSampDistance (int num);
void kmeansDistance();
void kBackTrack(int num);
void segKMeansProcess();
void MultiGaussianProcess();
#endif