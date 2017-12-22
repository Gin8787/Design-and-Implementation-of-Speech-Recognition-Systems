#ifndef KTRAIN_H
#define KTRAIN_H

#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>
#include<map>

using namespace std;

#define STATE (5)
#define INF (10000000)
#define MEAN (1)
#define COVAR (2)

struct iso_sample{
	vector<vector<vector<double> > > s0;
	vector<vector<vector<double> > > s1;
	vector<vector<vector<double> > > s2;
	vector<vector<vector<double> > > s3;
	vector<vector<vector<double> > > s4;
	vector<vector<vector<double> > > s5;
	vector<vector<vector<double> > > s6;
	vector<vector<vector<double> > > s7;
	vector<vector<vector<double> > > s8;
	vector<vector<vector<double> > > s9;
	vector<vector<vector<double> > > so;
};

struct seg_sample{
	vector<vector<int> > s0;
	vector<vector<int> > s1;
	vector<vector<int> > s2;
    vector<vector<int> > s3;
	vector<vector<int> > s4;
	vector<vector<int> > s5;
	vector<vector<int> > s6;
	vector<vector<int> > s7;
	vector<vector<int> > s8;
	vector<vector<int> > s9;
	vector<vector<int> > so;
};

void readSampleDCT(const char* filename, int flag);
void constructSamples(const char* filename, int flag);
void initSeg();
void testPrintSample();
void calSegmentMeansCovar();
void calSampDistance (int num);
void kmeansDistance();
void kBackTrack(int num);
void segKMeansProcess();
void assigndigits();
#endif