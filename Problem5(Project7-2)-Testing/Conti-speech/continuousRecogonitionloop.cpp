#include"continuousRecogonitionloop.h"

vector<vector < double> > templateMean;
vector<vector < double> > templateCor;
vector<vector < double> > inputFrames;
vector<vector < double> > distMatrix;
vector<vector < double> > DTWMatrix;
vector<vector < int> > backtrackMatrix;
vector<string > fileNameList;
vector<vector<int> > phoneNumList;
vector<int> isBackward;
vector<string > phoneNumString;
vector<string > rightNum;

int meanFlag = 0;
int minCurrent;
int corFlag = 1;
int diagonalOp = 0;
int horizonOp = 1;
int backwardOp = 2;
int rightSetence = 0;
int wrongWords = 0;
int totalWords = 0;
bool singleword = false;
double backwardPenalty = 50;

void readTemplate(int flag, string fileName) {
//	cout << "---------In readTemplate----------\n";
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
//	cout << "---------template reading finished----------\n";
}

void getInputFrames(string inputFile) {
	int numSamples, sampleRate;
	//delete rawData at the end of the program, because it is newed in the function
	//readWavFile of readwave.cpp
//	short *rawData = ReadWavFile(inputFile.c_str(), &numSamples, &sampleRate);
//	string DCTFileName = extractFeatureOfAudio(rawData, numSamples, sampleRate);
//	cout << "Reading normalization of input audio............\n";
	//used to debug
//	string DCTFileName = "4679.txt";
//	ifstream inputNorm("hwdata/test_mfcc/FCS_4B_new.txt");
	
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
		//Optional Optimization for X mfcc!!!!!!!
		inputFrames.pop_back();
//		cout << "count is " << count << endl;
//		cout << "inputFrames size is " << inputFrames.size() << endl;
	}
//	cout << "Reading normalization of input audio finished\n";
}



void calDistMatrix() {
//	cout << "Calculating every node cost now.............\n";
	int row, col, k;
	int inputLen = inputFrames.size();
	double tmpDistance = 0;
	//notice that the column here in the distanceMatrix, in fact is a row
	vector<double> column;
	for (row = 0; row < inputLen; row++) {
//		cout << "row is " << row << endl;
		for (col = 0; col < TOTALSTATES; col++) {
			tmpDistance = 0;
//			cout << "col is " << col << endl;
			for (k = 0; k < 39; k++) {
				tmpDistance += log(2 * Pi * templateCor[col][k]) + (inputFrames[row][k] - templateMean[col][k]) * (inputFrames[row][k] - templateMean[col][k]) / templateCor[col][k];
				//tmpDistance += (inputFrames[row][k] - templateMean[col][k]) * (inputFrames[row][k] - templateMean[col][k]);
			}
			//tmpDistance = sqrt(tmpDistance);
			tmpDistance = 0.5 * tmpDistance;
			column.push_back(tmpDistance);
		}
		distMatrix.push_back(column);
		column.clear();
	}
//	cout << "Calculating every node cost now finished\n";
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
	if(result >= INF) {
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
	if(result >= INF) {
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
//	cout << "Calculating DTW matrix now.............\n";
	for (row = 0; row < inputLen; row++) {
//		cout << "row is " << row << endl;
		miniCost = INF;
		operation = -1;
		for (col = 0; col < TOTALSTATES; col++) {
			curNode = distMatrix[row][col];
			if (row == 0) {
				//cout << "row == 0 broken\n";
				//cout << "inputLen is " << inputLen << ", row is " << row << ", col is " << col << endl;
				tmpCost[col] = distMatrix[row][col];
			}
			else if ((row != 0) && (col % 5 == 0)) {
				//cout << "row == 0 && (col % 5 == 0) broken\n";
				//cout << "inputLen is " << inputLen << ", row is " << row << ", col is " << col << endl;
				backward = miniCostPre + curNode + backwardPenalty;
				horizon = DTWMatrix[row - 1][col] + curNode;
				tmpCost[col] = kminBackward(backward, horizon, &operation);
				if (operation == backwardOp) {
					isBackward[row] = digitIndex;
				}
			}
			else {
				//cout << "else broken\n";
				//cout << "inputLen is " << inputLen << ", row is " << row << ", col is " << col << endl;
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
		//cout << "push data brocken\n";
		//cout << "inputLen is " << inputLen << ", row is " << row << ", col is " << col << endl;
		DTWMatrix.push_back(tmpCost);
		backtrackMatrix.push_back(tmpOperation);
		miniCostPre = miniCost;
	}
//	cout << "Calculating DTW matrix now finished\n";
	return digitIndex;
}

void backtrack(vector<int > &phoneNum, int endState) {
//	cout << "Doing backtrack now.............\n";
	int row, col, number, operation;
	int silence = -1;
	row = DTWMatrix.size() - 1;
	col = endState;
//    cout << "row is " << row << endl;
	number = silence;
	while (row >= 0) {
		operation = backtrackMatrix[row][col];
//		cout << "col is " << col << endl;
		if (((col % 5 == 0) &&(operation == backwardOp)&&(row != 0)) || (row == 0)) {
		    if (col < 50) {
				number = col / 5;
				cout << number << " ";
				phoneNum.push_back(number);
			}
			else if((col < 55) && (col >= 50)){
				number = 0;
				cout << number << " ";
				phoneNum.push_back(number);
			}
			else {
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
			col = isBackward[row];
			row--;
		}
	}
//	cout << "\nDoing backtrack finished\n";
}

//调用这个函数来完成多次循环
void phoneNumTest() {
	ofstream out_correspond("rrrr.txt");
	if (!out_correspond.is_open()) {
		cout << "open file recgonition_list.txt failed\n";
	}
	string path = "hwdata/";
	string dataFile, wavNameNew;
	string filelist;
	/*string meanFile[TOTALNUM] =
	{ "n-0.txt", "n-1.txt", "n-2.txt",
	"n-3.txt", "n-4.txt", "n-5.txt",
	"n-6.txt", "n-7.txt", "n-8.txt",
	"n-9.txt", "n-o.txt"};*/
	string meanFile[TOTALNUM] =
	{ "0.txt", "1.txt", "2.txt",
	"3.txt", "4.txt", "5.txt",
	"6.txt", "7.txt", "8.txt",
	"9.txt", "o.txt" };
	/*string corFile[TOTALNUM] =
	{ "nc-0.txt", "nc-1.txt", "nc-2.txt",
	"nc-3.txt", "nc-4.txt", "nc-5.txt",
	"nc-6.txt", "nc-7.txt", "nc-8.txt",
	"nc-9.txt", "nc-o.txt"};*/
	string corFile[TOTALNUM] =
	{ "c-0.txt", "c-1.txt", "c-2.txt",
	"c-3.txt", "c-4.txt", "c-5.txt",
	"c-6.txt", "c-7.txt", "c-8.txt",
	"c-9.txt", "c-o.txt"};
	int i, endState, fileNum, nameLen;
	vector<int > phoneNum;
	dataFile = "test/";
	filelist = "TEST.filelist";
	filelist = path + filelist;
	getFileNameList(filelist);
	fileNum = fileNameList.size();
	for (i = 0; i < TOTALNUM; i++) {
		cout << "mean file is " << meanFile[i] << endl;
		cout << "cov file is "<< corFile[i] << endl;
		readTemplate(meanFlag, meanFile[i]);
		readTemplate(corFlag, corFile[i]);
	}
	//for (i = 0; i < fileNum; i++) {                     //循环的次数如果是所有test文件，就把fileNum放入for条件，写i<fileNum
	for (i = 0; i < fileNum; i++){
		nameLen = fileNameList[i].size();
		if(singleword){
			if (nameLen > 6) {
			//if (nameLen <=6 || nameLen >=8){
				continue;
			}
		}
		out_correspond << fileNameList[i] << endl;

		//wavNameNew = path + dataFile + fileNameList[i] + "_new.txt";
		wavNameNew = path + dataFile + fileNameList[i];
		cout << wavNameNew << endl;
//		getInputFrames("hwdata/test_mfcc/MBE_OB_new.txt");       //WAVNameNew是每一个input文件的名字,多次循环时，传入wavNameNew做参数
		getInputFrames(wavNameNew);
		calDistMatrix();
		endState = DTW();
		backtrack(phoneNum, endState);
		inputFrames.clear();
		distMatrix.clear();
		DTWMatrix.clear();
		backtrackMatrix.clear();
		isBackward.clear();
		if(singleword){
			phoneNum.push_back((int)fileNameList[i][4] - '0');
			//phoneNum.push_back((int)fileNameList[i][5] - '0');

		}
		phoneNumList.push_back(phoneNum);
		phoneNum.clear();
	}
	writeResult();
	if (!singleword){
		readResult();
		getNumFromFileNum();
		rightSentenceCal();
		cout << "number of right sectences is " << rightSetence << endl;
		cout << "total number of words is " << totalWords << endl;
		cout << "number of wrong words is " << wrongWords << endl;
	}

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

void writeResult(){
	int correct = 0;
	ofstream out("recgonition_list.txt");
	if (!out.is_open()) {
		cout << "open file recgonition_list.txt failed\n";
	}
	else {
		int testNum = phoneNumList.size();
		int i, j, numLen;
		for (i = 0; i < testNum; i++) {
			numLen = phoneNumList[i].size();
			if (singleword){
				if (	
					(phoneNumList[i][numLen-1] == phoneNumList[i][0]) ||
					(phoneNumList[i][numLen-1] == 31 && phoneNumList[i][0] == 0) ||
					(phoneNumList[i][numLen-1] == 42 && phoneNumList[i][0] == 0) //|| 
					/*(phoneNumList[i][numLen-2] == phoneNumList[i][1]) ||
					(phoneNumList[i][numLen-2] == 31 && phoneNumList[i][1] == 0) ||
					(phoneNumList[i][numLen-2] == 42 && phoneNumList[i][1] == 0)  */
				)
				{
					correct++;
				}
			}

			for (j = numLen-1; j >= 0; j--) {
				
				out << phoneNumList[i][j];
			}
			if (i != testNum - 1){
				out << endl;
			}

		}
		if(singleword){
			out << correct << endl;
		}
	}
	out.close();
}

void initial(int tempLen, int inputLen, int dist[][100]) {
	int i, j;
	minCurrent = INF;
	for (i = 0; i <= tempLen; i++) {
		for (j = 0; j <= inputLen; j++) {
			dist[i][j] = INF;
//			prune[i][j] = false;
		}
	}
	dist[0][0] = 0;
}

int mini(int vertical, int horizon, int diagonal) {
	int temp = (vertical < horizon) ? vertical : horizon;
	return (temp < diagonal) ? temp : diagonal;
}

int CalEditDistance(string Template, string Input){
	int row, col, insertCost, deleteCost, subCost;
	int tempLen = Template.size();
	int inputLen = Input.size();
	int cost = 0;
	int dist[100][100];

	initial(tempLen, inputLen, dist);

	for (col = 0; col < inputLen; col++) {
		minCurrent = INF;
		for (row = 0; row < tempLen; row++) {
			if (Template[row] == Input[col])
				cost = 0;
			else
				cost = 1;
			if ((col != 0) && (row != 0)) {
				insertCost = dist[row - 1][col] + 1;
				deleteCost = dist[row][col - 1] + 1;
				subCost = dist[row - 1][col - 1] + cost;
				dist[row][col] = mini(insertCost, deleteCost, subCost);
			}
			else if ((col == 0) && (row != 0)) {      //Calculate the first column
				dist[row][col] = dist[row - 1][col] + 1;
			}
			else if ((row == 0) && (col != 0)) {     //Calculate the first row
				dist[row][col] = dist[row][col - 1] + 1;
			}		
			if (dist[row][col] < minCurrent) minCurrent = dist[row][col];		//update the minimum edit distance 
		}
		//prune
/*		for (row = 0; row < tempLen; row++) {
			if (pruneOption == fix) {
				if (dist[row][col] > threshold) {
					prune[row][col] = true;
				}
			}
			else {
				if (dist[row][col] > (minCurrent + threshold))
					prune[row][col] = true;
			}
		}*/
	}
//	printTrellis(tempLen, inputLen);
	return dist[row - 1][col - 1];
}

void readResult() {
	cout << "readResult\n";
    ifstream in("recgonition_list.txt");
	string temp;
	if (!in.is_open()) {
		cout << "open file recgonition_list.txt failed\n";
	} else {
		while(!in.eof()) {
		    in >> temp;
			temp = "#" + temp;
			phoneNumString.push_back(temp);
			temp.clear();
		}
	}
	cout << "readResult finished\n";
}

void getNumFromFileNum() {
	cout << "getFileNum" << endl;
	int i, fileNum, phoneLen, j, nameLen;
	string temp;
	fileNum = fileNameList.size();
//	fileNum = 100;
	for(i = 0; i < fileNum; i++) {		
	    nameLen = fileNameList[i].size();
		if (singleword){
			if (nameLen > 6) {
				continue;
			}
		}
		j = nameLen - 1;
		while(fileNameList[i][j] != '_') {
			char ch = fileNameList[i][j];
			if((( ch >= '0') && (ch <= '9'))) {
			    temp.push_back(ch);
			}
			if(ch == 'Z') {
			    temp.push_back('0');
			}
			if(ch == 'O') {
			    temp.push_back('0');
			}
			j--;
		}		
		totalWords += temp.size();
		reverse(temp.begin(), temp.end());
		temp = "#" + temp;
		rightNum.push_back(temp);
		temp.clear();
	}
	cout << "getFileNum  finished\n" << endl;
}

void rightSentenceCal() {
	cout << "rightSentenceCal\n";
	int cost, i;
	int len = phoneNumString.size();
	rightSetence = 0;
    wrongWords = 0;
	cout << "len is " << len << endl;
	for(i = 0; i < len; i++) {
//		cout << "i is " << i << endl; 
	    string inputS = phoneNumString[i];
		string templateS = rightNum[i];
		cost = CalEditDistance(templateS, inputS);
		if(cost != 0) {
		    wrongWords += cost;
		} else {
		    rightSetence++;
		}
	}
	cout << "rightSentenceCal finished\n";
}