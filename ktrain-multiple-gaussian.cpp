/*
*
* This project is aiming for multiple(4) gaussian mixture
* recognition and its sample training.
*
*/

#include "ktrain.h"

vector<vector<vector<double> > > samples;

vector<vector<vector<double> > >state1;//Plus
vector<vector<vector<double> > >state2;//Sub
vector<vector<vector<double> > >state3;
vector<vector<vector<double> > >state4;
vector<vector<vector<double> > >tempstate;
vector<vector<vector<double> > >tempstate2;

vector<vector<double> > mean_st1;
vector<vector<double> > mean_st2;
vector<vector<double> > mean_st3;
vector<vector<double> > mean_st4;
vector<vector<double> > covar_st1;
vector<vector<double> > covar_st2;
vector<vector<double> > covar_st3;
vector<vector<double> > covar_st4;

vector<vector<double> > means;
vector<vector<double> > covars;
vector<vector<double> > kdist;

double trans [SEGNUM+1][SEGNUM];
int segdiv [SAMPLENUM][SEGNUM+1];
//0 for state1-num; 1 for state2-num; 2 for overall-num;
//39 for dimensions;
double weight [39][SEGNUM][GAUSNUM];

void readSampleDCT(const char* filename, int flag){
	samples.resize(SAMPLENUM);
	int i;
	for (i = 0; i < SAMPLENUM; i++){
		samples[i].resize(39);
	}
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
			samples[flag][count].push_back(DCT_value);
			count ++;
		}
	}
}

void initSeg()
{
	//Read all samples from file
	readSampleDCT("s-0-1.txt", 0);
    readSampleDCT("s-0-2.txt", 1);
	readSampleDCT("s-0-3.txt", 2);
    readSampleDCT("s-0-4.txt", 3);
    readSampleDCT("s-0-5.txt", 4);
	int i, j;
	int frameLen;
	//Initialize all of the segments in even distribution
	for (i = 0; i <SAMPLENUM; i++){
		frameLen = samples[i][0].size();
		for (j = 0; j <SEGNUM; j++){
			segdiv[i][j] = (frameLen / 5) * (j);
			//cout << segdiv[i][j] << " ";
	    }
		segdiv[i][j] = frameLen;
		//cout << segdiv[i][j] << " ";
	}
	cout << endl;
}

void initTrans(){
	//Transition cost initialization
    int i, j, k, l;
	double count = 0.0;
	for (i = 0; i < SEGNUM + 1; i++){
		for(j = 0; j < SEGNUM; j++){
			if(i == 0){
				if (j == 0){
				    trans[i][j] = 0;
				}
				else{
					trans[i][j] = INF;
				}
			}
			else{
				if ((i - j) == 1 || (i == j)){
					
					for (k = 0; k < SAMPLENUM; k++){
						for (l = segdiv[k][i-1]; l < segdiv[k][i]; l++)
						{count++;}
					}
					//Initialization of transition cost based on the number of samples.
					//e.g. T(11)
					if (( i - j ) == 1){
						trans[i][j] = -log((count-SAMPLENUM)/count);
					}
					//e.g. T(12)
					else{
						trans[i][j] = -log(SAMPLENUM / count);
					}

					count = 0;

				}
				else{
					trans[i][j] = INF;
				}
			}
		}

	}
}

void calSegmentMeansCovar(){
	means.resize(39);
	covars.resize(39);
	double tmpMean = 0.0;
	int count = 0;
	int samp, dimens, frame, seg;
	//Calculate 5 means of each segment.
	for(seg = 0; seg < SEGNUM; seg ++){
	    //Calculate 39 means of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment means of all samples.
	        for(samp = 0; samp < SAMPLENUM; samp++){
			//Get the segment range.
			    for (frame = segdiv[samp][seg]; frame < segdiv[samp][seg+1]; frame++){
				    tmpMean += samples[samp][dimens][frame];
				    count++;
			    }
		    }
		    means[dimens].push_back(tmpMean/count);
		    tmpMean = 0;
		    count = 0;
	    }
	}

    double tmpCovar = 0.0;

    //Calculate 5 covars of each segment.
	for(seg = 0; seg < SEGNUM; seg ++){
	    //Calculate 39 covars of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment covars of all samples.
	        for(samp = 0; samp < SAMPLENUM; samp++){
			//Get the segment range.
			    for (frame = segdiv[samp][seg]; frame < segdiv[samp][seg+1]; frame++){
				    tmpCovar += (samples[samp][dimens][frame] - means[dimens][seg]) * (samples[samp][dimens][frame] - means[dimens][seg]);
				    count++;
			    }
		    }
		    covars[dimens].push_back(tmpCovar/count);
		    tmpCovar = 0;
		    count = 0;
	    }
	}
}

void calSampDistance (int num){
	//Warning! Call this function only after means and covar is fully initialized!
	ofstream disfile("distance.txt");

	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = samples[num][1].size();//Frame Length
	int meansize = means[1].size();
	kdist.resize(inputsize);

	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < meansize; j++)
		{
			for ( k = 0; k < 39; k++)
			{
				//The single Guassian distance distribution formula.
				tmpDistance += log( 2 * 3.14159 * covars[k][j] ) + (samples[num][k][i] - means[k][j]) * (samples[num][k][i] - means[k][j]) / covars[k][j];
			}
			tmpDistance = 0.5 * tmpDistance;
			kdist[i].push_back(tmpDistance);
			tmpDistance = 0;
		}
	}
	// Distance file test output.
	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < inputsize; i++)
		{
			for (j = 0; j < meansize; j++)
			{
				disfile << kdist[i][j] << ' ';
			}
			disfile << endl;
		}
	}
}

double kmin3(double superdiag, double horizon, double diag){
	double temp = (superdiag < horizon) ? superdiag : horizon;
	return (temp < diag) ? temp : diag;
}

double kmin2(double horizon, double diag){
	return (horizon < diag) ? horizon : diag;
}

void kmeansDistance(){
	//Warning! Call this function only after kdist vector is fully initialized!
	ofstream disfile("keditdistance.txt");
	int insize = kdist.size();
	int tmpsize = kdist[0].size();
	int i, j;
	for (i = 0; i < insize; i++){
		for(j = 0; j < tmpsize; j++){
			//If i!= 0 and j!= 0;
			if (!(i == 0 && j == 0)) {
				//If (i,j) is in the range of calculate
				if ( j > i * 2 ){
					kdist[i][j] = INF;
				}
				else{
					//Add (transition cost)/(edge cost) based on DTW.
					if (j == 0){
						kdist[i][j] = kdist[i-1][j] + kdist[i][j];
					}
					else{

						kdist[i][j] = kmin2(kdist[i-1][j], kdist[i-1][j-1]) + kdist[i][j];
					}
				}
			}
		}
	}
}

void kBackTrack(int num){
	//Warning! Only after kmeansdistance() invokes, the function can be called.
	//Will update segdiv[][].
	//num is the (num)th sample.
	int insize = kdist.size();
	int tmpsize = kdist[0].size();
	int inp = insize - 1;
	int tmp = tmpsize - 1;
	while(inp >=0){
		//(tmp)th segdiv element will be updated
		//The boundary between segment will be updated.
		segdiv[num][tmp] = inp;
		if (tmp == 0){
			inp--;
		}
		else{
			//If the state has not change
			if (kmin2(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1]) == kdist[inp - 1][tmp]){
				inp--;
			}
			//If the state has change for 1.
			else{
				inp--;
				tmp--;
			}
		}
		/*else{
			//If the state has not change
			if (kmin3(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1], kdist[inp - 1][tmp - 2]) == kdist[inp - 1][tmp]){
				inp--;
			}
			//If the state is down for 1 state.
			else if (kmin3(kdist[inp - 1][tmp], kdist[inp - 1][tmp - 1], kdist[inp - 1][tmp - 2]) == kdist[inp - 1][tmp - 1]){
				inp--;
				tmp--;
			}
			//If the state is down for 2 states.
			else{
				inp--;
				tmp=tmp-2;
			}
		}*/
	}
}

void clearforLoop(){
	//Clear global vectors for further usage.
	means.clear();
    covars.clear();
    kdist.clear();
}

void outputMean(const char* filename){
	//Output mean data to file.
	int i, j;
	ofstream meanfile(filename);
	if(!meanfile) {
			cout << "Fail to open distance inputfile! \n";
		}
		else{
			for (j = 0; j < means[0].size(); j++){
			    for (i = 0; i <39; i++){
					meanfile << means[i][j] << " ";
				}
				meanfile << endl;
			}
		}
}

void outputCovar(const char* filename){
	//Output covar data to file.
	int i, j;
	ofstream covarfile(filename);
	if(!covarfile) {
			cout << "Fail to open distance inputfile! \n";
		}
		else{
			for (j = 0; j < covars[0].size(); j++){
			    for (i = 0; i <39; i++){
					covarfile << covars[i][j] << " ";
				}
				covarfile << endl;
			}
		}
}

void testPrintSample(){
    int i, j;
	//cout << samples.size() << endl; //==2
	//cout << samples[0].size() << endl; //==39
	//cout << samples[1].size() << endl; //==39
	//cout << samples[0][0].size() << endl; //==Frame Length
	//这样输出来的每一个元素都是第i维的第j帧的元素
	//The (1)st element belongs to the (i)th dimension and (j)th frame.
    /*for (i = 0; i <39; i++){
		for (j = 0; j < samples[0][0].size(); j++){
			cout << samples[1][i][j] << " ";
	    }
		cout << endl;
	}*/

  for(i = 0; i < SAMPLENUM; i++){
	for (j = 0; j <SEGNUM+1; j++){
		cout << segdiv[i][j] << " ";
	}
	cout << endl;
  }
}

void segKMeansProcess(){
	//The training process of seg-kmeans.
	//Segment evenly divided.
	initSeg();
	int i, j;
	/*int sum1 = 0;
	int sum2 = 0;
	for(i = 0; i < SAMPLENUM; i++){
		for (j = 0; j <SEGNUM+1; j++){
				sum2 += segdiv[i][j];
			}
	}*/
	//Testing iteration process
	for (i = 0; i < 10; i++){
	/*while (1){*/
		//Clear vectors
		clearforLoop();
		//Init transition costs.
	    initTrans();
		//Calculate means an covariances
        calSegmentMeansCovar();
		for (j = 0; j < SAMPLENUM; j++){
		    //In each sample, doing sample distance calculation.
	        calSampDistance(j);
			//Update the kmeans distance.
            kmeansDistance();
			//Update segment division array.
	        kBackTrack(j);
			kdist.clear();
		}
		//testPrintSample();
		//Calculate the sum of segments to see if it is converged.
		/*for(i = 0; i < SAMPLENUM; i++){
			for (j = 0; j <SEGNUM+1; j++){
				sum1 += segdiv[i][j];
			}
		}
	    testPrintSample();
		if (sum2 == sum1) break;
		sum2 = sum1;
		sum1 = 0;*/
	}
	outputMean("means.txt");
	outputCovar("covars.txt");
}

void initState(){

	int i;
	state1.resize(SEGNUM);
	state2.resize(SEGNUM);
	for (i = 0; i < SEGNUM; i++){
		state1[i].resize(39);
		state2[i].resize(39);
	}
}

void initMultiMeanCovar(){

	mean_st1.resize(39);
	mean_st2.resize(39);
	covar_st1.resize(39);
	covar_st2.resize(39);

}

void clearWeight(){
	int i, j, k;
	for(i = 0; i < 39; i++){
		for(j = 0; j < SEGNUM; j++){
			for(k = 0; k < GAUSNUM/2; k++){
				weight[i][j][k] = 0;
			}
		}
	}
}

void clearStates(){
	state1.clear();
	state2.clear();
}

void clearMultiMeanCovar(){
	mean_st1.clear();
	mean_st2.clear();
	covar_st1.clear();
	covar_st2.clear();
}

void outputWeight(const char* weightfile){
	int i, j, k;
	ofstream file(weightfile);
	for(i = 0; i < 39; i++){
		for(j = 0; j < SEGNUM; j++){
			for(k = 0; k < GAUSNUM; k++){
				file << weight[i][j][k] << " ";
			}
			if(!(i == 38 && j == SEGNUM-1)){
				file << endl;
			}

		}
	}
}

void splitGaussian(int flag){
	ofstream file1("sp1.txt");
	ofstream file2("sp2.txt");
	int i, j, k;
	//Use this function only after finishing SINGLE GAUSSIAN segmental k-means operation!
	int dimens, seg, frame, samp;
	double meanplus = 0.0, meansub = 0.0;
	double covarplus = 0.0, covarsub = 0.0;
	double clusterdist_plus = 0.0, clusterdist_sub = 0.0;
	double current = 0.0;

	
	//Each segment has 39 means value with respect to each dimension.
	for (dimens = 0; dimens < 39; dimens++){
		//Operation in every dimension
		for (seg = 0; seg < SEGNUM; seg++){
			//add a small epsilon in every means, split the gaussian into half.
			if (flag == INIT){
				meanplus = means[dimens][seg] + EPSL;
				meansub = means[dimens][seg] - EPSL;
				covarplus = covars[dimens][seg];
				covarsub = covars[dimens][seg];
			}
			if (flag == ITER){
				meanplus = mean_st1[dimens][seg];
				meansub = mean_st2[dimens][seg];
				covarplus = covar_st1[dimens][seg];
				covarsub = covar_st2[dimens][seg];
			}
			//The number of samples is SAMPLENUM(e.g. 5)
			for(samp = 0; samp < SAMPLENUM; samp++){
				//In each sample, frames have already segmented by the first step of Single Guassian.
				for (frame = segdiv[samp][seg]; frame < segdiv[samp][seg+1]; frame++){
					current = samples[samp][dimens][frame];
					clusterdist_plus = log( 2 * 3.14159 * covarplus ) + 
							(current - meanplus) *  (current - meanplus) / covarplus;
					clusterdist_sub = log( 2 * 3.14159 * covarsub ) + 
							(current - meansub) * (current - meansub) / covarsub;
					//Split into two states.
					if(clusterdist_plus < clusterdist_sub){
						state1[seg][dimens].push_back(current);
						weight[dimens][seg][0] += 1;
					}
					else{
						state2[seg][dimens].push_back(current);
						weight[dimens][seg][1] += 1;
					}
				}
			}
		}
	}

	file1 << "STATE-1 :" << endl;
	for (i = 0; i < state1.size(); i++){
		for (j = 0; j < state1[i].size(); j++){
			for( k = 0; k < state1[i][j].size(); k++ ){

				file1 << state1[i][j][k] << " ";

			}
			file1 << " | ";
		}
		
	}

	file2 << "STATE-2 :" << endl;

	for (i = 0; i < state2.size(); i++){
		for (j = 0; j < state2[i].size(); j++){
			for( k = 0; k < state2[i][j].size(); k++ ){

				file2 << state2[i][j][k] << " ";
			}
			file2 << " | ";
		}
	
	}
}

void splitfromStates(int flag){
	ofstream file1("sp1.txt");
	ofstream file2("sp2.txt");
	int i, j, k;
	//Use this function only after finishing SINGLE GAUSSIAN segmental k-means operation!
	int dimens, seg, frame, samp;
	double meanplus = 0.0, meansub = 0.0;
	double covarplus = 0.0, covarsub = 0.0;
	double clusterdist_plus = 0.0, clusterdist_sub = 0.0;
	double current = 0.0;

	
	//Each segment has 39 means value with respect to each dimension.
	for (dimens = 0; dimens < 39; dimens++){
		//Operation in every dimension
		for (seg = 0; seg < SEGNUM; seg++){
			//add a small epsilon in every means, split the gaussian into half.
			if (flag == INIT){
				meanplus = means[dimens][seg] + EPSL;
				meansub = means[dimens][seg] - EPSL;
				covarplus = covars[dimens][seg];
				covarsub = covars[dimens][seg];
			}
			if (flag == ITER){
				meanplus = mean_st1[dimens][seg];
				meansub = mean_st2[dimens][seg];
				covarplus = covar_st1[dimens][seg];
				covarsub = covar_st2[dimens][seg];
			}
			for(frame = 0; frame < tempstate[seg][dimens].size(); frame++){
				current = tempstate[seg][dimens][frame];
				clusterdist_plus = log( 2 * 3.14159 * covarplus ) + 
						(current - meanplus) *  (current - meanplus) / covarplus;
				clusterdist_sub = log( 2 * 3.14159 * covarsub ) + 
						(current - meansub) * (current - meansub) / covarsub;
				//Split into two states.
				if(clusterdist_plus < clusterdist_sub){
					state1[seg][dimens].push_back(current);
					weight[dimens][seg][0] += 1;
				}
				else{
					state2[seg][dimens].push_back(current);
					weight[dimens][seg][1] += 1;
				}
			}
		}
	}

	file1 << "STATE-1 :" << endl;
	for (i = 0; i < state1.size(); i++){
		for (j = 0; j < state1[i].size(); j++){
			for( k = 0; k < state1[i][j].size(); k++ ){

				file1 << state1[i][j][k] << " ";

			}
			file1 << " | ";
		}
		
	}

	file2 << "STATE-2 :" << endl;

	for (i = 0; i < state2.size(); i++){
		for (j = 0; j < state2[i].size(); j++){
			for( k = 0; k < state2[i][j].size(); k++ ){

				file2 << state2[i][j][k] << " ";
			}
			file2 << " | ";
		}
	
	}
}

void calState1_MeansCovar(){
	//Under STATE 1!!!
	double tmpMean = 0.0;
	int size = 0;
	int count = 0;
	int i, j;
	int seg, frames, dimens;
	
	//Calculate 5 means of each segment.
	for(seg = 0; seg < state1.size(); seg++){
	    //Calculate 39 means of 39 dimensions.
	    for (dimens = 0; dimens < state1[seg].size(); dimens++){
		    //Calculate segment means of all samples.
			size = state1[seg][dimens].size();
	        for(frames = 0; frames < size; frames++){
				tmpMean += state1[seg][dimens][frames];
			}
		    mean_st1[dimens].push_back(tmpMean/size);
		    tmpMean = 0;
	    }
	}

    double tmpCovar = 0.0;

    //Calculate 5 covars of each segment.
	for(seg = 0; seg < state1.size(); seg ++){
	    //Calculate 39 covars of 39 dimensions.
	    for (dimens = 0; dimens < state1[seg].size(); dimens++){
		    //Calculate segment covars of all samples.
			size = state1[seg][dimens].size();
	        for(frames = 0; frames < size; frames++){
				tmpCovar += (state1[seg][dimens][frames] - mean_st1[dimens][seg]) * (state1[seg][dimens][frames] - mean_st1[dimens][seg]);
			}
		    covar_st1[dimens].push_back(tmpCovar/size);
		    tmpCovar = 0;
	    }
	}
}

void calState2_MeansCovar(){
	//UNDER state 2!!!
	double tmpMean = 0.0;
	int size = 0;
	int count = 0;
	int i, j;
	int seg, frames, dimens;
	
	//Calculate 5 means of each segment.
	for(seg = 0; seg < state2.size(); seg++){
	    //Calculate 39 means of 39 dimensions.
	    for (dimens = 0; dimens < state2[seg].size(); dimens++){
		    //Calculate segment means of all samples.
			size = state2[seg][dimens].size();
	        for(frames = 0; frames < size; frames++){
				tmpMean += state2[seg][dimens][frames];
			}
		    mean_st2[dimens].push_back(tmpMean/size);
		    tmpMean = 0;
	    }
	}

    double tmpCovar = 0.0;

    //Calculate 5 covars of each segment.
	for(seg = 0; seg < state2.size(); seg ++){
	    //Calculate 39 covars of 39 dimensions.
	    for (dimens = 0; dimens < state2[seg].size(); dimens++){
		    //Calculate segment covars of all samples.
			size = state2[seg][dimens].size();
	        for(frames = 0; frames < size; frames++){
				tmpCovar += (state2[seg][dimens][frames] - mean_st2[dimens][seg]) * (state2[seg][dimens][frames] - mean_st2[dimens][seg]);
			}
		    covar_st2[dimens].push_back(tmpCovar/size);
		    tmpCovar = 0;
	    }
	}

	/*ofstream file1("s2m.txt");
	ofstream file2("s2c.txt");

	file1 << "STATE-2 mean:" << endl;
	for (j = 0; j < mean_st2[0].size(); j++){
		for (i = 0; i < mean_st2.size(); i++){
			file1 << mean_st2[i][j] << " ";
		}
		file1 << endl;
	}

	file2 << "STATE-2 covar:" << endl;
	for (j = 0; j < covar_st2[0].size(); j++){
		for (i = 0; i < covar_st2.size(); i++){
				file2 << covar_st2[i][j] << " ";
		}
		file2 << endl;
	}*/
}

void outputMulti_MeansCovar(const char* means, const char* covars, int flag){
	int i, j;
	ofstream file1(means);
	ofstream file2(covars);

	switch (flag)
	{
		case 1:
		{
			for (j = 0; j < mean_st1[0].size(); j++){
				for (i = 0; i < mean_st1.size(); i++){
					file1 << mean_st1[i][j] << " ";
				}
				if(j != mean_st1[0].size() - 1){
					file1 << endl;
				}
			}

			for (j = 0; j < covar_st1[0].size(); j++){
				for (i = 0; i < covar_st1.size(); i++){
						file2 << covar_st1[i][j] << " ";
				}
				if(j != covar_st1[0].size() - 1){
					file2 << endl;
				}
			}
			break;
		}

		case 2:
		{
			for (j = 0; j < mean_st2[0].size(); j++){
				for (i = 0; i < mean_st2.size(); i++){
					file1 << mean_st2[i][j] << " ";
				}
				if(j != mean_st2[0].size() - 1){
					file1 << endl;
				}
			}

			for (j = 0; j < covar_st2[0].size(); j++){
				for (i = 0; i < covar_st2.size(); i++){
						file2 << covar_st2[i][j] << " ";
				}
				if(j != covar_st2[0].size() - 1){
					file2 << endl;
				}
			}
			break;
		}

		case 3:
		{
			for (j = 0; j < mean_st3[0].size(); j++){
				for (i = 0; i < mean_st3.size(); i++){
					file1 << mean_st3[i][j] << " ";
				}
				if(j != mean_st3[0].size() - 1){
					file1 << endl;
				}
			}

			for (j = 0; j < covar_st3[0].size(); j++){
				for (i = 0; i < covar_st3.size(); i++){
						file2 << covar_st3[i][j] << " ";
				}
				if(j != covar_st3[0].size() - 1){
					file2 << endl;
				}
			}
			break;
		}

		case 4:
		{
			for (j = 0; j < mean_st4[0].size(); j++){
				for (i = 0; i < mean_st4.size(); i++){
					file1 << mean_st4[i][j] << " ";
				}
				if(j != mean_st4[0].size() - 1){
					file1 << endl;
				}
			}

			for (j = 0; j < covar_st4[0].size(); j++){
				for (i = 0; i < covar_st4.size(); i++){
						file2 << covar_st4[i][j] << " ";
				}
				if(j != covar_st4[0].size() - 1){
					file2 << endl;
				}
			}
			break;
		}

	default:
		break;
	}
}

void calDoubleDistance(int num){

	ofstream disfile("distance_double.txt");

	int i, j, k;

	double tmpDistance = 0.0;
	int inputsize = samples[num][1].size();//Frame Length
	int meansize = means[1].size();
	kdist.resize(inputsize);

	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < meansize; j++)
		{
			for ( k = 0; k < 39; k++)
			{
				//The single Guassian distance distribution formula.
				tmpDistance += ( ( log( 2 * 3.14159 * covar_st1[k][j] ) + 
					(samples[num][k][i] - mean_st1[k][j]) * (samples[num][k][i] - mean_st1[k][j]) / covar_st1[k][j] ) *
					 (  weight[k][j][0] / (weight[k][j][0] + weight[k][j][1]) ) ) +

					(  ( log( 2 * 3.14159 * covar_st2[k][j] ) + 
					(samples[num][k][i] - mean_st2[k][j]) * (samples[num][k][i] - mean_st2[k][j]) / covar_st2[k][j] ) *
					( weight[k][j][1] / (weight[k][j][0] + weight[k][j][1]) ) );
			}
			tmpDistance = 0.5 * tmpDistance;
			kdist[i].push_back(tmpDistance);
			tmpDistance = 0;
		}
	}
	// Distance file test output.
	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < inputsize; i++)
		{
			for (j = 0; j < meansize; j++)
			{
				disfile << kdist[i][j] << ' ';
			}
			disfile << endl;
		}
	}
	/*int dimens, seg;
	means.clear();
	covars.clear();
	means.resize(39);
	covars.resize(39);
	for (dimens = 0; dimens < 39; dimens++){
		means[dimens].resize(SEGNUM);
		covars[dimens].resize(SEGNUM);
	}
	//Call this function only after means and covars of 2 states are fully updated!
	for (dimens = 0; dimens < means.size(); dimens++){
		for (seg = 0; seg < means[dimens].size(); seg++){
			means[dimens][seg] =  * mean_st1[dimens][seg] 
				+ weight[dimens][seg][1] / (weight[dimens][seg][0]+weight[dimens][seg][1]) * mean_st2[dimens][seg];
			covars[dimens][seg] = weight[dimens][seg][0] / (weight[dimens][seg][0]+weight[dimens][seg][1]) * covar_st1[dimens][seg] 
				+ weight[dimens][seg][1] / (weight[dimens][seg][0]+weight[dimens][seg][1]) * covar_st2[dimens][seg];
		}
	}*/
}

void calQuadDistance(int num){

	ofstream disfile("distance_quad.txt");

	int i, j, k;

	double tmpDistance = 0.0;
	double overallWeight = 0.0;
	int inputsize = samples[num][1].size();//Frame Length
	int meansize = means[1].size();
	kdist.resize(inputsize);

	for (i = 0; i < inputsize; i++)
	{
		for (j = 0; j < meansize; j++)
		{
			for ( k = 0; k < 39; k++)
			{
				overallWeight = weight[k][j][0] + weight[k][j][1] + weight[k][j][2] + weight[k][j][3];
				//The single Guassian distance distribution formula.
				tmpDistance += ( ( log( 2 * 3.14159 * covar_st1[k][j] ) + 
					(samples[num][k][i] - mean_st1[k][j]) * (samples[num][k][i] - mean_st1[k][j]) / covar_st1[k][j] ) *
					 (  weight[k][j][0] / overallWeight ) ) +

					(  ( log( 2 * 3.14159 * covar_st2[k][j] ) + 
					(samples[num][k][i] - mean_st2[k][j]) * (samples[num][k][i] - mean_st2[k][j]) / covar_st2[k][j] ) *
					( weight[k][j][1] / overallWeight ) ) +

					(  ( log( 2 * 3.14159 * covar_st3[k][j] ) + 
					(samples[num][k][i] - mean_st3[k][j]) * (samples[num][k][i] - mean_st3[k][j]) / covar_st3[k][j] ) *
					( weight[k][j][3] / overallWeight ) ) +

					(  ( log( 2 * 3.14159 * covar_st4[k][j] ) + 
					(samples[num][k][i] - mean_st4[k][j]) * (samples[num][k][i] - mean_st4[k][j]) / covar_st4[k][j] ) *
					( weight[k][j][4] / overallWeight ) );
			}
			tmpDistance = 0.5 * tmpDistance;
			kdist[i].push_back(tmpDistance);
			tmpDistance = 0;
		}
	}
	// Distance file test output.
	if(!disfile) {
		cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < inputsize; i++)
		{
			for (j = 0; j < meansize; j++)
			{
				disfile << kdist[i][j] << ' ';
			}
			disfile << endl;
		}
	}

}

void copyState(vector<vector<vector<double> > > &v1, vector<vector<vector<double> > > &v2){
	v2.assign(v1.begin(), v1.end());
}

void copy2d(vector<vector<double> > &v1, vector<vector<double> > &v2){
	v2.assign(v1.begin(), v1.end());
}

void copyWeight(){
	int i, j, k;
	for(i = 0; i < 39; i++){
		for(j = 0; j < SEGNUM; j++){
			for(k = 0; k < GAUSNUM / 2; k++){
				weight[i][j][k+GAUSNUM/2] = weight[i][j][k];
			}
		}
	}
}

void splitInitProcess(int flag){
	clearStates();
	initState();
	clearMultiMeanCovar();
	initMultiMeanCovar();
	clearWeight();
	if (flag == INIT){
		splitGaussian(INIT);
	}
	if (flag == ITER){
		splitfromStates(INIT);
	}
	calState1_MeansCovar();
	calState2_MeansCovar();
}

void splitIterProcess(int flag){

	clearStates();
	initState();
	clearWeight();
	if (flag == INIT){
		splitGaussian(ITER);
	}
	if (flag == ITER){
		splitfromStates(ITER);
	}
	clearMultiMeanCovar();
	initMultiMeanCovar();
	calState1_MeansCovar();
	calState2_MeansCovar();
	//outputWeight();
}

void clearall(){
	state1.clear();//Plus
	state2.clear();//Sub
	state3.clear();
	state4.clear();
	tempstate.clear();
	tempstate2.clear();

	mean_st1.clear();
	mean_st2.clear();
	mean_st3.clear();
	mean_st4.clear();
	covar_st1.clear();
	covar_st2.clear();
	covar_st3.clear();
	covar_st4.clear();
}


void outputall(){
	outputMulti_MeansCovar("multiple\\m-0-1.txt", "multiple\\c-0-1.txt", 1);
	outputMulti_MeansCovar("multiple\\m-0-2.txt", "multiple\\c-0-2.txt", 2);
	outputMulti_MeansCovar("multiple\\m-0-3.txt", "multiple\\c-0-3.txt", 3);
	outputMulti_MeansCovar("multiple\\m-0-4.txt", "multiple\\c-0-4.txt", 4);
	outputWeight("multiple\\w-0.txt");
}

/*void MultiGaussianProcess(){

	int a, b;
	int i;
	splitInitProcess(INIT);

	//Iterate inside every state.
	for (i = 0; i < 7; i++){
		splitIterProcess(INIT);
		//outputWeight();
	}
	testPrintSample();
	initTrans();

	//Split two gaussians into four

	//State 2 will split into State 3 and State 4.
	clearforLoop();

	copyState(state2, tempstate);
	copyState(state2, state4); //Preserve state2
	copyState(state1, tempstate2);//temp store
	copy2d(mean_st2, means);
	copy2d(mean_st1, mean_st4);//temp store
	copy2d(covar_st2, covars);
	copy2d(covar_st1, covar_st4);//temp store
	
	splitInitProcess(ITER);

	for (i = 0; i < 7; i++){
		splitIterProcess(ITER);
		//outputWeight();
	}

	//clearforLoop();
	//State1 split into State 1 and State 2.
	copyState(tempstate2, tempstate);
	copyState(state4, tempstate2);
	copyState(state1, state3);
	copyState(state2, state4);

	copy2d(mean_st4, means);
	copy2d(mean_st1, mean_st3);
	copy2d(mean_st2, mean_st4);
	copy2d(covar_st4, covars);
	copy2d(covar_st1, covar_st3);
	copy2d(covar_st2, covar_st4);

	copyWeight();
	splitInitProcess(ITER);
	for (i = 0; i < 7; i++){
		splitIterProcess(ITER);
	}
	initTrans();
	testPrintSample();
	cout << endl;
	outputall();
}*/

void MultiGaussianProcess(){

	int a, b;
	int i;
	//for(b = 0; b < 2; b++){
		clearWeight();
		for(a = 0; a < 4; a++){
			splitInitProcess(INIT);

			//Iterate inside every state.
			for (i = 0; i < 4; i++){
				splitIterProcess(INIT);
//				outputWeight();
				//cout << endl;
			}
			//testPrintSample();
			//outputMean("meansplit.txt");
			//cout << endl;
			initTrans();

			for (i = 0; i < SAMPLENUM; i++){
				//In each sample, doing sample distance calculation.
				calDoubleDistance(i);
				//Update the kmeans distance.
				kmeansDistance();
				//Update segment division array.
				kBackTrack(i);
				kdist.clear();
			}
			testPrintSample();
			cout << "------------------------" << endl;
		}

		//Split two gaussians into four

		//State 2 will split into State 3 and State 4.
		clearforLoop();

		copyState(state2, tempstate);
		copyState(state2, state4); //Preserve state2
		copyState(state1, tempstate2);//temp store
		copy2d(mean_st2, means);
		copy2d(mean_st1, mean_st4);//temp store
		copy2d(covar_st2, covars);
		copy2d(covar_st1, covar_st4);//temp store
	
		splitInitProcess(ITER);

		for (i = 0; i < 4; i++){
			splitIterProcess(ITER);
			//outputWeight();
		}

		//clearforLoop();
		//State1 split into State 1 and State 2.
		copyState(tempstate2, tempstate);
		copyState(state4, tempstate2);
		copyState(state1, state3);
		copyState(state2, state4);

		copy2d(mean_st4, means);
		copy2d(mean_st1, mean_st3);
		copy2d(mean_st2, mean_st4);
		copy2d(covar_st4, covars);
		copy2d(covar_st1, covar_st3);
		copy2d(covar_st2, covar_st4);

		copyWeight();
		//outputWeight();
		splitInitProcess(ITER);
		for (i = 0; i < 4; i++){
			splitIterProcess(ITER);
		}
		//initTrans();
		testPrintSample();
		//outputMean("meansplit.txt");
		cout << endl;

		//if( b == 0 ){
			for (i = 0; i < SAMPLENUM; i++){
				//In each sample, doing sample distance calculation.
				calQuadDistance(i);
				//Update the kmeans distance.
				kmeansDistance();
				//Update segment division array.
				kBackTrack(i);
				kdist.clear();
			}
		//}
		//testPrintSample();
		cout << "------------------------" << endl;
		//Clear vectors
		clearforLoop();
		//Init transition costs.
	    //initTrans();
		//Calculate means an covariances
        calSegmentMeansCovar();
	//}
	outputall();
}