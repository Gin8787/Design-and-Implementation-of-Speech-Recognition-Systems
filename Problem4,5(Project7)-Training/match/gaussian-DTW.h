#ifndef GAUSSIANDTW_H
#define GAUSSIANDTW_H

#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<cctype>
#include<vector>
#include<sstream>
#include<Windows.h>
#include<algorithm>

using namespace std;

#define SEGNUM (5)
#define INF (10000000)

//double min3(double superdiag, double horizon, double diag);
//double min2(double horizon, double diag);
//void readDCT(const char* filename);
//void calDistance(const char* inputfilename, const char* samplefilename);
//double DTWeditDistance();
//void clearDist();
void kProcess();
//void outTrans();
//void readTrans(const char* filename);
void singleProcess();

#endif