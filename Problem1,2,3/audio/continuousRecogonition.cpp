#include"continuousRecogonition.h"

vector<vector < double> > templateMean;
vector<vector < double> > templateCor;
vector<vector < double> > inputFrames;
vector<vector < double> > distMatrix;
vector<vector < double> > DTWMatrix;
vector<vector < int> > backtrackMatrix;
vector<string > fileNameList;
vector<vector<int> > phoneNumList;
vector<int> isBackward;

int meanFlag = 0;
int corFlag = 1;
int diagonalOp = 0;
int horizonOp = 1;
int backwardOp = 2;
double backwardPenalty = 150;

void readTemplate(int flag, string fileName) {
	cout << "---------In readTemplate----------\n";
	int i, j;
	ifstream templateIn(fileName.c_str());
	if (!templateIn.is_open()) {
		cout << "Open file " << fileName << " failed!\n";
	}
	else {
		vector<double > stateTemp;
		double temp;
		for (i = 0; i < STATENUM; i++) {
			for (j = 0; j < FRAMELENGTH; j++) {
				templateIn >> temp;
				stateTemp.push_back(temp);
			}
			if (flag == meanFlag) {
				templateMean.push_back(stateTemp);
			}
			else {
				templateCor.push_back(stateTemp);
			}
			stateTemp.clear();
		}
		
	}
	cout << "---------template reading finished----------\n";
}

void getInputFrames(string inputFile) {
	int numSamples, sampleRate;
	//delete rawData at the end of the program, because it is newed in the function
	//readWavFile of readwave.cpp
//	short *rawData = ReadWavFile(inputFile.c_str(), &numSamples, &sampleRate);
//	string DCTFileName = extractFeatureOfAudio(rawData, numSamples, sampleRate);
	cout << "Reading normalization of input audio............\n";
	//used to debug
//	string DCTFileName = "4679.txt";
	ifstream inputNorm(inputFile.c_str());
	int DCTindex = 0;
	vector<double > frameTemp;
	double temp;
	int count = 0;
	if (!inputNorm.is_open()) {
		cout << "Open file " << inputFile << " failed!\n";
	}
	else {
		DCTindex = 0;
		while (!inputNorm.eof()) {
			if (DCTindex == 39) {
				inputFrames.push_back(frameTemp);
				DCTindex = 0;
				frameTemp.clear();
				count++;
			}
			inputNorm >> temp;
			frameTemp.push_back(temp);
			DCTindex++;
		}
		inputFrames.push_back(frameTemp);
		cout << "count is " << count << endl;
		cout << "inputFrames size is " << inputFrames.size() << endl;
	}
	cout << "Reading normalization of input audio finished\n";
}



void calDistMatrix() {
	cout << "Calculating every node cost now.............\n";
	int row, col, k;
	int inputLen = inputFrames.size();
	double tmpDistance = 0;
	//notice that the column here in the distanceMatrix, in fact is a row
	vector<double> column;
	for (row = 0; row < inputLen; row++) {
		for (col = 0; col < TOTALSTATES; col++) {
			tmpDistance = 0;
			for (k = 0; k < 39; k++) {
//				tmpDistance += log(2 * Pi * templateCor[col][k]) + (inputFrames[row][k] - templateMean[col][k]) * (inputFrames[row][k] - templateMean[col][k]) / templateCor[col][k];
				tmpDistance += (inputFrames[row][k] - templateMean[col][k]) * (inputFrames[row][k] - templateMean[col][k]);
			}
			tmpDistance = sqrt(tmpDistance);
//			tmpDistance = 0.5 * tmpDistance;
			column.push_back(tmpDistance);
		}
		distMatrix.push_back(column);
		column.clear();
	}
	cout << "Calculating every node cost now finished\n";
}

double kmin(double diagonal, double horizon, int *operation) {
	double result;
	if (diagonal < horizon) {
		result = diagonal;
		*operation = diagonalOp;
	}
	else {
		result = horizon;
		*operation = horizonOp;
	}
	if (result >= INF) {
		result = INF;
		*operation = -1;
	}
	return result;
}

double kminBackward(double backward, double horizon, int *operation) {
	double result;
	if (backward < horizon) {
		result = backward;
		*operation = backwardOp;
	}
	else {
		result = horizon;
		*operation = horizonOp;
	}
	if (result >= INF) {
		result = INF;
		*operation = -1;
	}
	return result;
}

/*void initialDTWMatrix() {
	int inputLen = distMatrix.size();
	vector<double> tmpCost;
	int i, j;
	for (j = 0; j < TOTALSTATES; j++) {
		tmpCost.push_back(INF);
	}
	for (i = 0; i < inputLen; i++) {
		DTWMatrix.push_back(tmpCost);
	}
}*/

void initialBacktrackMatrix() {
	int inputLen = distMatrix.size();
	vector<int > tmpOperation;
	int i;
	for (i = 0; i < inputLen; i++) {	
		isBackward.push_back(-1);    //here  -1 means no backward
	}
}

int DTW() {
	int row, col, operation;
	int inputLen = inputFrames.size();
	vector<double > tmpCost;
	vector<int > tmpOperation;
	int digitIndex;
	double miniCost, miniCostPre, diagonal, horizon, backward, curNode;
	miniCost = miniCostPre = INF;
	digitIndex = -1;
//	initialDTWMatrix();
	initialBacktrackMatrix();
	for (col = 0; col < TOTALSTATES; col++) {
		tmpCost.push_back(INF);
		tmpOperation.push_back(-1);
	}
	cout << "Calculating DTW matrix now.............\n";
	for (row = 0; row < inputLen; row++) {
//		cout << "row is " << row << endl;
		miniCost = INF;
		operation = -1;
		for (col = 0; col < TOTALSTATES; col++) {
			curNode = distMatrix[row][col];
			if (row == 0) {
				if (col % 5 == 0){
					tmpCost[col] = distMatrix[row][col];
				}
				else {
					tmpCost[col] = INF;
				}
			}
			else if ((row != 0) && (col % 5 == 0)) {
				backward = miniCostPre + curNode + backwardPenalty;
				horizon = DTWMatrix[row - 1][col] + curNode;
				tmpCost[col] = kminBackward(backward, horizon, &operation);
				if (operation == backwardOp) {
					isBackward[row] = digitIndex;
				}
			}
			else {
				diagonal = DTWMatrix[row - 1][col - 1] + curNode;
				horizon = DTWMatrix[row - 1][col] + curNode;
				tmpCost[col] = kmin(diagonal, horizon, &operation);
				if ((col + 1) % 5 == 0) {
					if (miniCost > tmpCost[col]) {
						miniCost = tmpCost[col];
						digitIndex = col;
					}
				}
			}
			tmpOperation[col] = operation;
		}
		DTWMatrix.push_back(tmpCost);
		backtrackMatrix.push_back(tmpOperation);
		miniCostPre = miniCost;
	}
	cout << "Calculating DTW matrix now finished\n";
	return digitIndex;
}

void backtrack(vector<int > &phoneNum, int endState) {
	cout << "Doing backtrack now.............\n";
	int row, col, number, operation;
	int silence = -1;
	row = DTWMatrix.size() - 1;
	col = endState;
	number = silence;
//	cout << "row is " << row << endl;
	while (row >= 0) {
//		cout << "col is " << col << endl;
//		cout << "row is " << row << endl;
		operation = backtrackMatrix[row][col];
//		if (((col % 5 == 0) &&(isBackward[row] == backwardOp)&&(row != 0)) || (row == 0)) {
		if (((col % 5 == 0) && (operation == backwardOp) && (row != 0)) || (row == 0)) {
		    if (col < 50) {
				number = col / 5;
				phoneNum.push_back(number);
			}
/*			else if((col < 55) && (col >= 50)){
				number = 0;
				phoneNum.push_back(number);
			}*/
			else {
				number = silence;
			}
//			operation = backwardOp;
		}		
		if (operation == diagonalOp) {
			row--;
			col--;
		}
		else if (operation == horizonOp) {
			row--;
		}
		else {
			col = isBackward[row];
			row--;
		}
	}
	cout << "Doing backtrack finished\n";
}

void phoneNumRecogonition() {
	string meanFile[11] =
	{ "m-0.txt", "m-1.txt", "m-2.txt",
	"m-3.txt", "m-4.txt", "m-5.txt",
	"m-6.txt", "m-7.txt", "m-8.txt" ,
	"m-9.txt", "silence-mean.txt"};
	string corFile[11] =
	{ "c-0.txt", "c-1.txt", "c-2.txt",
	"c-3.txt", "c-4.txt", "c-5.txt",
	"c-6.txt", "c-7.txt", "c-8.txt",
	"c-9.txt", "silence-cov.txt" };
	int i, endState;
	vector<int > phoneNum;
	for (i = 0; i < 11; i++) {
		readTemplate(meanFlag, meanFile[i]);
		readTemplate(corFlag, corFile[i]);
	}
	getInputFrames("next/09-5.txt");
	calDistMatrix();
	endState = DTW();
	backtrack(phoneNum, endState);
	cout << "The phone number you said is: \n";
	for (i = 0; i < phoneNum.size(); i++) {
		cout << phoneNum[i] << " ";
	}
}

//调用这个函数来完成多次循环
void phoneNumTest() {
	string path = "hwdata/";
	string dataFile, wavNameNew;
	string filelist;
/*	string meanFile[12] =
	{ "n-0.txt", "n-1.txt", "n-2.txt",
	"n-3.txt", "n-4.txt", "n-5.txt",
	"n-6.txt", "n-7.txt", "n-8.txt",
	"n-9.txt","silence-mean.txt" };*/
	string meanFile[12] =
	{ "0.txt", "1.txt", "2.txt",
	"3.txt", "4.txt", "5.txt",
	"6.txt", "7.txt", "8.txt",
	"9.txt", "o.txt", "silence-mean.txt" };
/*	string corFile[12] =
	{ "nc-0.txt", "nc-1.txt", "nc-2.txt",
	"nc-3.txt", "nc-4.txt", "nc-5.txt",
	"nc-6.txt", "nc-7.txt", "nc-8.txt",
	"nc-9.txt", "silence-cov.txt" };*/
	string corFile[12] =
	{ "c-0.txt", "c-1.txt", "c-2.txt",
	"c-3.txt", "c-4.txt", "c-5.txt",
	"c-6.txt", "c-7.txt", "c-8.txt",
	"c-9.txt", "c-o.txt","silence-cov.txt" };
	int i, endState, fileNum, nameLen;
	vector<int > phoneNum;
	dataFile = "test_mfcc/";
	filelist = "TEST.filelist";
	filelist = path + filelist;
	getFileNameList(filelist);
	fileNum = fileNameList.size();
	for (i = 0; i < 12; i++) {
		readTemplate(meanFlag, meanFile[i]);
		readTemplate(corFlag, corFile[i]);
	}
	for (i = 0; i < 1; i++) {                     //循环的次数如果是所有test文件，就把fileNum放入for条件，写i<fileNum
		nameLen = fileNameList[i].size();
		if (nameLen < 9) {
			backwardPenalty = 200;
		}
		else if (nameLen < 12) {
			backwardPenalty = 100;
		}
		else {
			backwardPenalty = 50;
		}
		backwardPenalty = 50;
//		wavNameNew = path + dataFile + fileNameList[i] + "_new.txt";
		getInputFrames("hwdata/test_mfcc/MCP_ZB_new.txt");       //WAVNameNew是每一个input文件的名字,多次循环时，传入wavNameNew做参数
		getInputFrames(wavNameNew);
		calDistMatrix();
		endState = DTW();
		backtrack(phoneNum, endState);
		inputFrames.clear();
		distMatrix.clear();
		DTWMatrix.clear();
		backtrackMatrix.clear();
		isBackward.clear();
		phoneNumList.push_back(phoneNum);
		phoneNum.clear();
	}
	writeReuslt();
}

void getFileNameList(string filelist){
	ifstream filelist_in(filelist.c_str());
	if (!filelist_in.is_open()) {
		cout << "Open file " << filelist << "failed\n";
	}
	else {
		string wavName;
		while (!filelist_in.eof()) {
			filelist_in >> wavName;
			fileNameList.push_back(wavName);
		}
	}
	filelist_in.close();
}

void writeReuslt(){
	ofstream out("recgonition_list.txt");
	if (!out.is_open()) {
		cout << "open file recgonition_list.txt failed\n";
	}
	else {
		int testNum = phoneNumList.size();
		int i, j, numLen;
		for (i = 0; i < testNum; i++) {
			numLen = phoneNumList[i].size();
			for (j = 0; j < numLen; j++) {
				out << phoneNumList[i][j] << " ";
			}
			out << endl;
		}
	}
}