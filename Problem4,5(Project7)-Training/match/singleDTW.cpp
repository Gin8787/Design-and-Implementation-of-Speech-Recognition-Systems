#include "singleDTW.h"
vector<vector<double>> input(39);
vector<vector<double>> sample(39);
vector<vector<double>> dist;

int templateStart[51] = {};

int SPEECH = 1;
int SAMP = 2;

void readDCT(const char* filename, int flag){
	//Read MFCC Spectrum from file.
	ifstream in_DCT(filename);
	string temp; 
	double DCT_value;
	int count = 0;
	if (!in_DCT.is_open()) {
		cout << "Open file failed! \n";
	}
	else {
		while (!in_DCT.eof()) {
			if (count >= 39){
				count = 0;
			}
		    in_DCT >> temp;
			DCT_value = std::stod(temp);
			if (flag == SPEECH){
				input[count].push_back(DCT_value);
			}
			else if (flag == SAMP){
				sample[count].push_back(DCT_value);
			}
			else{
				break;
			}
			count ++;
		}
	}
	// input[i][j]: iMAX = 38 (0-38, 39 Dimensions); jMAX = Number of Frames.
	
}


//Read template in time sync way.
void readAllTemplate(string samplefilename[10][5]) {
	int number, index;
	string filename;
	string temp;
	double DCT_val;
	int frameTotalCnt, templateCnt, featureCnt;
	templateCnt = 0;
	frameTotalCnt = 0;
	featureCnt = 0;
	templateStart[templateCnt] = frameTotalCnt;
	for (number = 0; number < 10; number++) {
		for (index = 0; index < 5; index++) {
			filename = samplefilename[number][index];
			ifstream in_template(filename);
			if (!in_template.is_open()) {
				cout << "Open file failed! \n";
			}
			else {
				while (!in_template.eof()) {
					if (featureCnt >= 39) {
						featureCnt = 0;
					}
					in_template >> temp;
					DCT_val = std::stod(temp);
					sample[featureCnt].push_back(DCT_val);
					featureCnt++;
					frameTotalCnt++;
				}
			}
			in_template.close();
			templateCnt++;
			templateStart[templateCnt] = frameTotalCnt + 1;
		}
	}
}

void calDistance(const char* inputfilename, const char* samplefilename){
	//Using global variable to calculate the distance between samples and input.
	ofstream disfile("distance.txt");
	readDCT(inputfilename, SPEECH);
	readDCT(samplefilename, SAMP);
	//read ALL template, time sync way, samplefilename is a string array
	//readAllTemplate(samplefilename);
	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = input[1].size();
	int samplesize = sample[1].size();
	//Calculate the Euclidean Distance between all frames and store it in vector.
	dist.resize(inputsize);
	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < samplesize; j++)
		{
			//Calculate the Euclidean Distance between given frames. 
			for ( k = 0; k < 39; k++)
			{
				tmpDistance += (input[k][i] - sample[k][j]) * (input[k][i] - sample[k][j]); 
			}
			tmpDistance = sqrt(tmpDistance);
			dist[i].push_back(tmpDistance);
			tmpDistance = 0;
		}
	}
	//Output distance information to file
	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < inputsize; i++)
		{
			for (j = 0; j < samplesize; j++)
			{
				disfile << dist[i][j] << ' ';
			}
			disfile << endl;
		}
	}
}

double min3(double superdiag, double horizon, double diag){
	double temp = (superdiag < horizon) ? superdiag : horizon;
	return (temp < diag) ? temp : diag;
}

double min2(double horizon, double diag){
	return (horizon < diag) ? horizon : diag;
}

double DTWeditDistance(){
	int insize = dist.size();
	int tmpsize = dist[0].size();
	int i, j;
	//For Pruning
	double mindistCur = 10000.0;
	double mindistPre = 10000.0;
	//Threshold of pruning
	double prCo = 10;
	for (i = 0; i < insize; i++){
		for(j = 0; j < tmpsize; j++){
			//If i!= 0 and j!= 0;
			if (!(i == 0 && j == 0)) {
				//If (i,j) is in the range of calculate
				if ( j > i * 2 ){
					dist[i][j] = INF;
				}
				else{
					if (j == 0){
						dist[i][j] = dist[i-1][j] + dist[i][j];
					}
					else if (j == 1){
						dist[i][j] = min2(dist[i-1][j], dist[i-1][j-1]) + dist[i][j];
					}
					else{
						dist[i][j] = min3(dist[i-1][j], dist[i-1][j-1], dist[i-1][j-2]) + dist[i][j];
					}

					if( dist[i][j] > mindistPre * (1 + prCo) ){
						dist[i][j] = INF;
					}

					if (mindistCur > dist[i][j]){
						mindistCur = dist[i][j];
					}
				}
			}
		}
		mindistPre = mindistCur;
		mindistCur = INF;
	}
	for(i = 0; i < 39; i++){
		input[i].clear();
        sample[i].clear();
	}
	return dist[insize-1][tmpsize-1];
}


double DTWeditDistanceSyn(){
	//Use syncronized placeholder to indicate the time sync template.
	int insize = dist.size();
	int tmpsize = dist[0].size();
	int i, j;
	int frameCnt, templateCnt;
	frameCnt = 0;
	templateCnt = 0;
	for (i = 0; i < insize; i++){
		for (j = 0; j < tmpsize; j++){
			//If i!= 0 and j!= 0;
			if (!(i == 0 && j == templateStart[templateCnt++])) {
				//If (i,j) is in the range of calculate
				if (j > i * 2){
					dist[i][j] = INF;
				}
				else{
					if (j == templateStart[templateCnt++]){
						dist[i][j] = dist[i - 1][j] + dist[i][j];
					}
					else if (j == templateStart[templateCnt++] + 1){
						dist[i][j] = min2(dist[i - 1][j], dist[i - 1][j - 1]) + dist[i][j];
					}
					else{
						dist[i][j] = min3(dist[i - 1][j], dist[i - 1][j - 1], dist[i - 1][j - 2]) + dist[i][j];
					}
				}
			}
		}
	}
	for (i = 0; i < 39; i++){
		input[i].clear();
		sample[i].clear();
	}
	return dist[insize - 1][tmpsize - 1];
}


void clearDist(){
	int i;
	for (i = 0; i < dist.size(); i++){
		dist[i].clear();
	}
}

void process(){
	int i, j;
	int a, b;
	int countright = 0;
	int count = 0;
     string sample[10][5] =
	{"s-0-1.txt","s-0-2.txt","s-0-3.txt","s-0-4.txt","s-0-5.txt",
	"s-1-1.txt","s-1-2.txt","s-1-3.txt","s-1-4.txt","s-1-5.txt",
	"s-2-1.txt","s-2-2.txt","s-2-3.txt","s-2-4.txt","s-2-5.txt",
    "s-3-1.txt","s-3-2.txt","s-3-3.txt","s-3-4.txt","s-3-5.txt",
	"s-4-1.txt","s-4-2.txt","s-4-3.txt","s-4-4.txt","s-4-5.txt",
	"s-5-1.txt","s-5-2.txt","s-5-3.txt","s-5-4.txt","s-5-5.txt",
	"s-6-1.txt","s-6-2.txt","s-6-3.txt","s-6-4.txt","s-6-5.txt",
	"s-7-1.txt","s-7-2.txt","s-7-3.txt","s-7-4.txt","s-7-5.txt",
	"s-8-1.txt","s-8-2.txt","s-8-3.txt","s-8-4.txt","s-8-5.txt",
	"s-9-1.txt","s-9-2.txt","s-9-3.txt","s-9-4.txt","s-9-5.txt"};

	 string input[10][5] =
	{"i-0-1.txt","i-0-2.txt","i-0-3.txt","i-0-4.txt","i-0-5.txt",
	"i-1-1.txt","i-1-2.txt","i-1-3.txt","i-1-4.txt","i-1-5.txt",
	"i-2-1.txt","i-2-2.txt","i-2-3.txt","i-2-4.txt","i-2-5.txt",
    "i-3-1.txt","i-3-2.txt","i-3-3.txt","i-3-4.txt","i-3-5.txt",
	"i-4-1.txt","i-4-2.txt","i-4-3.txt","i-4-4.txt","i-4-5.txt",
	"i-5-1.txt","i-5-2.txt","i-5-3.txt","i-5-4.txt","i-5-5.txt",
	"i-6-1.txt","i-6-2.txt","i-6-3.txt","i-6-4.txt","i-6-5.txt",
	"i-7-1.txt","i-7-2.txt","i-7-3.txt","i-7-4.txt","i-7-5.txt",
	"i-8-1.txt","i-8-2.txt","i-8-3.txt","i-8-4.txt","i-8-5.txt",
	"i-9-1.txt","i-9-2.txt","i-9-3.txt","i-9-4.txt","i-9-5.txt"};

	double min[10][5] = {};
	double tmpMin = INF;
	int number = 0;
	for (a = 0; a < 10; a++){
		for (b = 0; b < 5; b++){
			number = 0;
			tmpMin = INF;
			for (i = 0; i < 10; i++)
			{
				for (j = 0; j < 4; j++)
				{
					calDistance(input[a][b].c_str(), sample[i][j].c_str());
					min[i][j] = DTWeditDistance();
					//cout << min[i][j] << " ";
					clearDist();
					if (tmpMin > min[i][j]) 
					{
							tmpMin = min[i][j];
							number = i;
					}	
				}
			}
			cout << "Guess:"<< number << " Real:"<< a << "|";
			if (number == a) countright++;
			count ++;
		}
		cout << endl;
	}
	cout << "Correct Ratio is: "<< countright << " out of " << count << endl;
}

void kProcess(){
	int i;
	int a, b;
	int countright = 0;
	int count = 0;
	string ksample[10] = {"0.txt","1.txt","2.txt","3.txt","4.txt","5.txt","6.txt","7.txt","8.txt","9.txt"};

	string input[10][5] =
	{"i-0-1.txt","i-0-2.txt","i-0-3.txt","i-0-4.txt","i-0-5.txt",
	"i-1-1.txt","i-1-2.txt","i-1-3.txt","i-1-4.txt","i-1-5.txt",
	"i-2-1.txt","i-2-2.txt","i-2-3.txt","i-2-4.txt","i-2-5.txt",
    "i-3-1.txt","i-3-2.txt","i-3-3.txt","i-3-4.txt","i-3-5.txt",
	"i-4-1.txt","i-4-2.txt","i-4-3.txt","i-4-4.txt","i-4-5.txt",
	"i-5-1.txt","i-5-2.txt","i-5-3.txt","i-5-4.txt","i-5-5.txt",
	"i-6-1.txt","i-6-2.txt","i-6-3.txt","i-6-4.txt","i-6-5.txt",
	"i-7-1.txt","i-7-2.txt","i-7-3.txt","i-7-4.txt","i-7-5.txt",
	"i-8-1.txt","i-8-2.txt","i-8-3.txt","i-8-4.txt","i-8-5.txt",
	"i-9-1.txt","i-9-2.txt","i-9-3.txt","i-9-4.txt","i-9-5.txt"};

	double min[10] = {};
	double tmpMin = INF;
	int number = 0;
	for (a = 0; a < 10; a++){
		for (b = 0; b < 5; b++){
			number = 0;
			tmpMin = INF;
			for (i = 0; i < 10; i++)
			{
				calDistance(input[a][b].c_str(), ksample[i].c_str());
				min[i] = DTWeditDistance();
				clearDist();
				if (tmpMin > min[i]) 
				{
						tmpMin = min[i];
						number = i;
				}	
			}
			cout << "Guess:"<< number << " Real:"<< a << "|";
			if (number == a) countright++;
			count ++;
		}
		cout << endl;
	}
	cout << "Correct Ratio is: "<< countright << " out of " << count << endl;
}


void singleProcess(){
	int i;
	//int a, b;
	int countright = 0;
	int count = 0;
	string ksample[10] = {"0.txt","1.txt","2.txt","3.txt","4.txt","5.txt","6.txt","7.txt","8.txt","9.txt"};

	/*string input[10][5] =
	{"i-0-1.txt","i-0-2.txt","i-0-3.txt","i-0-4.txt","i-0-5.txt",
	"i-1-1.txt","i-1-2.txt","i-1-3.txt","i-1-4.txt","i-1-5.txt",
	"i-2-1.txt","i-2-2.txt","i-2-3.txt","i-2-4.txt","i-2-5.txt",
    "i-3-1.txt","i-3-2.txt","i-3-3.txt","i-3-4.txt","i-3-5.txt",
	"i-4-1.txt","i-4-2.txt","i-4-3.txt","i-4-4.txt","i-4-5.txt",
	"i-5-1.txt","i-5-2.txt","i-5-3.txt","i-5-4.txt","i-5-5.txt",
	"i-6-1.txt","i-6-2.txt","i-6-3.txt","i-6-4.txt","i-6-5.txt",
	"i-7-1.txt","i-7-2.txt","i-7-3.txt","i-7-4.txt","i-7-5.txt",
	"i-8-1.txt","i-8-2.txt","i-8-3.txt","i-8-4.txt","i-8-5.txt",
	"i-9-1.txt","i-9-2.txt","i-9-3.txt","i-9-4.txt","i-9-5.txt"};*/
	string input = "i-7-1.txt";
	double min[10] = {};
	double tmpMin = INF;
	int number = 0;
	//for (a = 0; a < 10; a++){
		//for (b = 0; b < 5; b++){
			number = 0;
			tmpMin = INF;
			for (i = 0; i < 10; i++)
			{
				calDistance(input.c_str(), ksample[i].c_str());
				min[i] = DTWeditDistance();
				clearDist();
				if (tmpMin > min[i]) 
				{
						tmpMin = min[i];
						number = i;
				}	
			}
			cout << "Guess:"<< number << " Real: 7 ";
			//if (number == a) countright++;
			//count ++;
		//}
	//	cout << endl;
	//}
	//cout << "Correct Ratio is: "<< countright << " out of " << count << endl;
}