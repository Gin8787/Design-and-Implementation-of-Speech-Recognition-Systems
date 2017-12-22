#ifndef _FIXEDDIGITRECOGNITION_H

#include "continuousRecogonition.h"
#define _CRT_SECURE_NO_WARNINGS

using namespace std;

void fixedDigitRec();
void initializeColumn();
void maximizeCost(double miniCost[]);
void initialBacktrackMatrix(int inputLen);
void initialBackwardPtr(int inputLen);
int DTWFixedRecgonition();
void backtrackFixed(vector<int > &phoneNum, int endState);
void initialbackPointer();
void backtrackWithPointer(vector<int > &phoneNum, int endState);
void initialbackPointer();
void initialStartFrame(int inputLen);

#endif