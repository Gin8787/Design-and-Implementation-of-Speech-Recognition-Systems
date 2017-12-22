#ifndef SINGLEDTW_H
#define SINGLEDTW_H

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

#define INF (100000)

double min3(double superdiag, double horizon, double diag);
double min2(double horizon, double diag);
void readDCT(const char* filename);
void calDistance(const char* inputfilename, const char* samplefilename);
//double DTWeditDistance(const char* filename);
double DTWeditDistance();
void clearDist();
void singleDTW();

#endif