#include "singleDTW.h"

vector<vector<double>> input(39);
vector<vector<double>> sample(39);
vector<vector<double>> dist;

int SPEECH = 1;
int SAMP = 2;

//void readDCT(string dct, vector<vector<double>> dctvector) {
void readDCT(const char* filename, int flag){
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

void calDistance(const char* inputfilename, const char* samplefilename){
	//readDCT("1.txt",input);
	//const char* inputfilename = "sample-6.txt";
	//const char* samplefilename = "input-6-1.txt";
	//ofstream disfile("distance.txt");
	readDCT(inputfilename, SPEECH);
	//cout << "Speech success!" << endl;
	readDCT(samplefilename, SAMP);
	//cout << "Sample success!" << endl;
	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = input[1].size();
	int samplesize = sample[1].size();
	//Calculate the Euclidean Distance between all frames and store it in vector.
    //cout << input[1].size();
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
/*	if(!disfile) {
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
	//cout << input[38][76];
	//cout << sample[38][92];*/
}

double min3(double superdiag, double horizon, double diag){
	double temp = (superdiag < horizon) ? superdiag : horizon;
	return (temp < diag) ? temp : diag;
}

double min2(double horizon, double diag){
	return (horizon < diag) ? horizon : diag;
}

//double DTWeditDistance(const char* filename){
double DTWeditDistance(){
	//ofstream disfile(filename);
	int insize = dist.size();
	int tmpsize = dist[0].size();
	int i, j;
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
				}
			}
		}
	}
	for(i = 0; i < 39; i++){
		input[i].clear();
        sample[i].clear();
	}

/*	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < insize; i++)
		{
			for (j = 0; j < tmpsize; j++)
			{
				disfile << dist[i][j] << ' ';
			}
			disfile << endl;
		}
	}*/

	return dist[insize-1][tmpsize-1];
}

void clearDist(){
	int i;
	for (i = 0; i < dist.size(); i++){
		dist[i].clear();
	}
}

void singleDTW(){
	int i, j;
     string sample[10] =
	{"0.txt",
	"1.txt",
	"2.txt",
    "3.txt",
	"4.txt",
	"5.txt",
	"6.txt",
	"7.txt",
	"8.txt",
	"9.txt"};
	double min[10] = {};
	double tmpMin = INF;
	int number = 0;
	for (i = 0; i < 10; i++)
	{
		calDistance("input.txt", sample[i].c_str());
	    min[i] = DTWeditDistance();
	    cout << min[i] << " ";
		//cout << sample[i][j].c_str() << " ";
	    clearDist();
	    if (tmpMin > min[i]) 
		{ 
			tmpMin = min[i];
			number = i;
		}	
	    cout << endl;
	}
	
	cout << "You said: "<< number << endl;
}