#include"FixedDigitRecognition.h"

double oddColumn7[7][50], evenColumn7[7][50];
int oddBackPointer[7][50], evenBackPointer[7][50];
vector<vector<vector<int > > > startFrame;
vector<vector<int > > backwardPtr;
vector<vector<vector<int > > > backtrackPtr;
vector<vector<int> > backtrackTable;
extern vector<vector < double> > distMatrix;
extern vector<vector < double> > inputFrames;

extern int meanFlag;
extern int corFlag;
extern int diagonalOp;
extern int horizonOp;
extern int backwardOp;
extern double backwardPenalty;

void fixedDigitRec() {
	string meanFile[10] =
	{ "0.txt", "1.txt", "2.txt",
	"3.txt", "4.txt", "5.txt",
	"6.txt", "7.txt", "8.txt",
	"9.txt"};
	string corFile[10] =
	{ "c-0.txt", "c-1.txt", "c-2.txt",
	"c-3.txt", "c-4.txt", "c-5.txt",
	"c-6.txt", "c-7.txt", "c-8.txt",
	"c-9.txt"};
	int i, endState;
	vector<int > phoneNum;
	int meanFlag = 0;
	int corFlag = 1;
	vector<int > phoneNumFixed;
	for (i = 0; i < 10; i++) {
		readTemplate(meanFlag, meanFile[i]);
		readTemplate(corFlag, corFile[i]);
	}
	getInputFrames("input.txt");
	calDistMatrix();
	endState = DTWFixedRecgonition();
	backtrackFixed(phoneNumFixed, endState);
	cout << "The phone number you said is: \n";
	for (i = phoneNumFixed.size() - 1; i >= 0; i--) {
		cout << phoneNumFixed[i] << " ";
	}
}

void initializeColumn() {
	cout << "In initialize Column\n";
	int i, stateCnt, temp;
	stateCnt = 0;
	for (i = 0; i < 7; i++) {
		for (stateCnt = 0; stateCnt < 50; stateCnt++) {
			if ((i != 0) && (i != 3)) {
				oddColumn7[i][stateCnt] = INF;
			}
			else if(i == 0) {
				if ((stateCnt >= 10) && (stateCnt % 5 == 0)) {
					oddColumn7[i][stateCnt] = distMatrix[0][stateCnt];
				}
				else {
					oddColumn7[i][stateCnt] = INF;
				}
			}
			else if ((i == 3) && (stateCnt % 5 == 0)) {
				oddColumn7[i][stateCnt] = distMatrix[0][stateCnt];
			}
			evenColumn7[i][stateCnt] = INF;
		}
	}
	cout << "initialize column finished\n";
}

void maximizeCost(double miniCost[]) {
	int i;
	for (i = 0; i < 7; i++) {
		miniCost[i] = INF;
	}
}

void initialBacktrackMatrix(int inputLen) {
	cout << "in initialize backtrack matrix\n";
	int i;
	vector<vector<int > > oneNum;
	vector<int> oneDigit;
	for (i = 0; i < 50; i++) {
		oneDigit.push_back(-1);
	}
	for (i = 0; i < 7; i++) {
		oneNum.push_back(oneDigit);
	}
	for (i = 0; i < inputLen; i++) {
		backtrackPtr.push_back(oneNum);
	}
	cout << "initialize backtrack matrix finished\n";
}

void initialStartFrame(int inputLen) {
	vector<int> stateStart;
	vector<vector<int> > digitStart;
	int i;
	for (i = 0; i < 10; i++){
		stateStart.push_back(-1);
	}
	for (i = 0; i < 7; i++) {
		digitStart.push_back(stateStart);
	}
	for (i = 0; i < inputLen; i++){
		startFrame.push_back(digitStart);
	}
}

void initialBackwardPtr(int inputLen) {
	cout << "in initialize backward matrix\n";
	int i, j;
	vector<int> digitBackward;
	for (j = 0; j < inputLen; j++) {
		digitBackward.push_back(-1);
	}
	for (i = 0; i < 7; i++) {
		backwardPtr.push_back(digitBackward);
	}
	cout << "initialize backward matrix finished\n";
}



void initialbackPointer(){
	cout << "in initialbackPointer\n";
	int i, j;
	for (i = 0; i < 7; i++){
		for (j = 0; j < 50; j++) {
			if ((j % 5 == 0) && (i == 0)){
				oddBackPointer[i][j] = 0;
			}
			else {
				oddBackPointer[i][j] = INF;
			}
			evenBackPointer[i][j] = INF;
		}
	}
	cout << "in initialbackPointer finished\n";
}



int DTWFixedRecgonition() {
	cout << "Doing DTW fixed\n";
	int frameCnt, stateCnt, digitCnt, operation, i, miniCostIndex[7];
	int inputLen = inputFrames.size();
	double miniCost[7], miniCostPre[7];
	double horizon, diagonal, backward, nodeCost, minRet, nullCost;
	vector<int> tableEntry;
	
	initialBacktrackMatrix(inputLen);
	initialBackwardPtr(inputLen);
	initializeColumn();
	initialStartFrame(inputLen);
	initialbackPointer();
	maximizeCost(miniCost);
	maximizeCost(miniCostPre);
	tableEntry.push_back(-1);
	tableEntry.push_back(-1);
	backwardPenalty = 20;
	for (frameCnt = 1; frameCnt < inputLen; frameCnt++) {
		maximizeCost(miniCost);
		for (digitCnt = 0; digitCnt < 7; digitCnt++) {
			for (stateCnt = 0; stateCnt < 50; stateCnt++) {
				nodeCost = distMatrix[frameCnt][stateCnt];
				if ((digitCnt == 0) && (stateCnt < 10))
					continue;

				if (frameCnt % 2 != 0) {
					if (stateCnt % 5 != 0) {
						diagonal = oddColumn7[digitCnt][stateCnt - 1] + nodeCost;
					}
					horizon = oddColumn7[digitCnt][stateCnt] + nodeCost;
				}
				else {
					if (stateCnt % 5 != 0) {
						diagonal = evenColumn7[digitCnt][stateCnt - 1] + nodeCost;
					}
					horizon = evenColumn7[digitCnt][stateCnt] + nodeCost;
				}

				if ((digitCnt != 0) && (stateCnt % 5 == 0)) {
					if ((frameCnt >= digitCnt * 5) || ((frameCnt >= (digitCnt - 3) * 5) && (digitCnt >= 3))){
						if (digitCnt != 1) {
							nullCost = 0;
						}
						else {
							nullCost = 5;
						}
						backward = miniCostPre[digitCnt - 1] + nodeCost + nullCost;
						minRet = kminBackward(backward, horizon, &operation);
					
					}
					else {
						minRet = INF;
						operation = horizonOp;

					}
					if (operation == backwardOp) {
						backwardPtr[digitCnt][frameCnt] = miniCostIndex[digitCnt - 1];
					}
				}
				else if ((digitCnt == 0) && (stateCnt % 5 == 0)){
					minRet = horizon;
					operation = horizonOp;
				}
				else {
					minRet = kmin(diagonal, horizon, &operation);
					if ((stateCnt + 1) % 5 == 0) {
						if (minRet < miniCost[digitCnt]) {
							miniCost[digitCnt] = minRet;
							miniCostIndex[digitCnt] = stateCnt;
						}

						 
					}
				}

				if (frameCnt % 2 != 0) {
					evenColumn7[digitCnt][stateCnt] = minRet;
				}
				else {
					oddColumn7[digitCnt][stateCnt] = minRet;
				}

				backtrackPtr[frameCnt][digitCnt][stateCnt] = operation;


				if (frameCnt % 2 != 0) {
					if (operation == backwardOp) {
						evenBackPointer[digitCnt][stateCnt] = frameCnt;
	
					}
					else {
						evenBackPointer[digitCnt][stateCnt] = oddBackPointer[digitCnt][stateCnt];
					}
					if (stateCnt == 54) {
						cout << "";
					}
					if (digitCnt == 6) {
						cout << "";
					}
					if ((stateCnt + 1) % 5 == 0) {
//						startFrame[frameCnt][digitCnt][(stateCnt - 4) / 5] = evenBackPointer[digitCnt][stateCnt];
					}
				}
				else {
					if (operation == backwardOp) {
						oddBackPointer[digitCnt][stateCnt] = frameCnt;
					} 
					else {
						oddBackPointer[digitCnt][stateCnt] = evenBackPointer[digitCnt][stateCnt];
					}
					if ((stateCnt + 1) % 5 == 0) {
//						startFrame[frameCnt][digitCnt][(stateCnt + 1) / 5] = oddBackPointer[digitCnt][stateCnt];
					}
				}
				
			}
			for (i = 0; i < 7; i++){
				miniCostPre[i] = miniCost[i];
			}
			
		}
	}
	cout << "Doing DTW fixed finished\n";
	return miniCostIndex[6];
}

void backtrackFixed(vector<int > &phoneNum, int endState) {
	int row, col, operation, number, silence, digitCnt;
	int inputLen = inputFrames.size();
	row = inputLen - 1;
	col = endState;
	silence = -1;
	digitCnt = 6;
	while (row >= 0) {
		operation = backtrackPtr[row][digitCnt][col];
		if (((col % 5 == 0) && (backwardPtr[digitCnt][row] != -1)) || (row == 0)) {
			if (col != 50) {
				number = col / 5;
				phoneNum.push_back(number);
			}
			else{
				number = silence;
			}
			operation = backwardOp;
		}
		if (operation == diagonalOp) {
			row--;
			col--;
		}
		else if (operation == horizonOp) {
			row--;
		}
		else {
			col = backwardPtr[digitCnt][row];
			row--;
			digitCnt--;
		}
	}
}

void backtrackWithPointer(vector<int > &phoneNum, int endState) {
	int frameCnt, digitCnt, stateCnt;
	int num;
	int inputLen = inputFrames.size();
	stateCnt = endState;
	frameCnt = inputLen - 1;
	digitCnt = 6;
	while (frameCnt >= 0) {
		if (frameCnt == 0) {
			break;
		}

		num = (stateCnt - 5) / 5;
		phoneNum.push_back(num);
		stateCnt = backwardPtr[digitCnt][frameCnt];
		frameCnt = startFrame[frameCnt][digitCnt][stateCnt];
		digitCnt--;
	}
}