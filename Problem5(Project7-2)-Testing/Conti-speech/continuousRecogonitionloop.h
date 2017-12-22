#ifndef _CONTINUOUSRECOGONITION_H

#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>
#include"readTest.h"

#define FRAMELENGTH 39
#define STATENUM 5
#define TOTALSTATES 55
#define TOTALNUM 11
#define Pi 3.14159
#define INF 999999

using namespace std;

void readTemplate(int flag, string fileName);
void getInputFrames(string inputFile);
void calDistMatrix();
double kmin(double diagonal, double horizon, int *operation);
double kminBackward(double backward, double horizon, int *operation);
void initialDTWMatrix();
void initialBacktrackMatrix();
int DTW();
void backtrack(vector<int > *phoneNum, int endState);
void phoneNumRecogonition();
void getFileNameList(string filelist);
void phoneNumTest();
void writeResult();
int CalEditDistance(string Template, string Input);
void readResult();
void getNumFromFileNum();
void rightSentenceCal();
void initial(int tempLen, int inputLen, int dist[][100]);
int mini(int vertical, int horizon, int diagonal);

#endif