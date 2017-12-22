#include "ktrain.h"

vector<vector<vector<double> > > samples;

iso_sample sample_all;
seg_sample seg_all;

vector<vector<double> > means(39);
vector<vector<double> > covars(39);
vector<vector<double> > kdist;

vector<vector<int> > segdiv;
vector<vector<double> > trans;

vector<vector<string > > transcripts;

void readSampleDCT(const char* filename){
	vector<vector<double> > sample_tmp(39);
	sample_tmp.resize(39);

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
			if(!in_DCT.eof()){
				in_DCT >> temp;
				DCT_value = std::stod(temp);
				sample_tmp[count].push_back(DCT_value);
				count ++;
			}
		}
		//If MFCC file has extra space/enter at the end of line
		//ADD THIS LINE!!!
		sample_tmp[count-1].pop_back();
	}
	samples.push_back(sample_tmp);
}

void constructSamples(const char* filename, int flag){
	//Construct means and covar file according to input.
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
			if(!in_DCT.eof()){
				if (flag == MEAN){
					means[count].push_back(DCT_value);
				}
				else if (flag == COVAR){
					covars[count].push_back(DCT_value);
				}
				else{
					break;
				}
			}
			count ++;
		}
	}
	// input[i][j]: iMAX = 38 (0-38, 39 Dimensions); jMAX = Number of Frames.
}

void continue2isolate(int samp, int push_num, int position){
	//Warning! Position cannot exceed the segment max length! 
	vector<int> seg_tmp;
	seg_tmp.resize(STATE+1);
	vector<vector<double> > sample_tmp(39);
	sample_tmp.resize(39);
	int dimens, frame;
	int i, tmp = 0;
	tmp = segdiv[samp][ position * STATE ];
	for (i = 0; i <= STATE; i++){
		seg_tmp[i] = ( segdiv[samp][position * STATE + i] - tmp );
	}

	int size = segdiv[samp].size() + 1;
	//int size = segdiv[samp].size() - 1;
	//Calculate 39 means of 39 dimensions.
	for (dimens = 0; dimens < 39; dimens++){
		//Calculate segment means of all samples.
		//The state is set to 5.
		for (frame = segdiv[samp][ position * STATE ]; frame < segdiv[samp][ (position + 1) * STATE ]; frame++){
			sample_tmp[dimens].push_back(samples[samp][dimens][frame]);
		}
	}
	switch(push_num){
	case 0:
		{
			sample_all.s0.push_back(sample_tmp);
			seg_all.s0.push_back(seg_tmp);
			break;
		}
	case 1:
		{
			sample_all.s1.push_back(sample_tmp);
			seg_all.s1.push_back(seg_tmp);
			break;
		}
	case 2:
		{
			sample_all.s2.push_back(sample_tmp);
			seg_all.s2.push_back(seg_tmp);
			break;
		}
	case 3:
		{
			sample_all.s3.push_back(sample_tmp);
			seg_all.s3.push_back(seg_tmp);
			break;
		}
	case 4:
		{
			sample_all.s4.push_back(sample_tmp);
			seg_all.s4.push_back(seg_tmp);
			break;
		}
	case 5:
		{
			sample_all.s5.push_back(sample_tmp);
			seg_all.s5.push_back(seg_tmp);
			break;
		}
	case 6:
		{
			sample_all.s6.push_back(sample_tmp);
			seg_all.s6.push_back(seg_tmp);
			break;
		}
	case 7:
		{
			sample_all.s7.push_back(sample_tmp);
			seg_all.s7.push_back(seg_tmp);
			break;
		}
	case 8:
		{
			sample_all.s8.push_back(sample_tmp);
			seg_all.s8.push_back(seg_tmp);
			break;
		}
	case 9:
		{
			sample_all.s9.push_back(sample_tmp);
			seg_all.s9.push_back(seg_tmp);
			break;
		}
	case 10:
		{
			sample_all.so.push_back(sample_tmp);
			seg_all.so.push_back(seg_tmp);
			break;
		}
	default:
		break;
	}
}

void initSeg(int samplenum, int segnum)
{
	int i, j;
	int frameLen;
	//Initialize all of the segments in even distribution
	for (i = 0; i <samplenum; i++){
		frameLen = samples[i][0].size();
		for (j = 0; j <segnum; j++){
			segdiv[i][j] = (frameLen / segnum) * (j);
	    }
		segdiv[i][j] = frameLen;
	}
}

void initTrans(int segnum, int samplenum){
	//Transition cost initialization
    int i, j, k, l;
	double count = 0.0;
	for (i = 0; i < segnum + 1; i++){
		for(j = 0; j < segnum; j++){
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
					
					for (k = 0; k < samplenum; k++){
						for (l = segdiv[k][i-1]; l < segdiv[k][i]; l++)
						{count++;}
					}
					//Initialization of transition cost based on the number of samples.
					//e.g. T(11)
					if (( i - j ) == 1){
						trans[i][j] = -log((count-samplenum)/count);
					}
					//e.g. T(12)
					else{
						trans[i][j] = -log(samplenum / count);
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

void calSegmentMeansCovar(int samplenum, int segnum){
	double tmpMean = 0.0;
	int count = 0;
	int samp, dimens, frame, seg;
	//Calculate 5 means of each segment.
	for(seg = 0; seg < segnum; seg ++){
	    //Calculate 39 means of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment means of all samples.
	        for(samp = 0; samp < samplenum; samp++){
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
	for(seg = 0; seg < segnum; seg ++){
	    //Calculate 39 covars of 39 dimensions.
	    for (dimens = 0; dimens < 39; dimens++){
		    //Calculate segment covars of all samples.
	        for(samp = 0; samp < samplenum; samp++){
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
	int inputsize = samples[num][0].size();//Frame Length
	int meansize = means[0].size();
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
						//kdist[i][j] = kdist[i-1][j] + trans[j+1][j] + kdist[i][j];
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
	segdiv[num][tmpsize] = insize;
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
	}
}

void clearforLoop(){
	//Clear global vectors for further usage.
	means.clear();
	means.resize(39);
    covars.clear();
	covars.resize(39);
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
				if( j != means[0].size() - 1 ){
					meanfile << endl;
				}

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
				if( j != covars[0].size() - 1 ){
					covarfile << endl;
				}
			}
		}
}

void outputTrans(const char* filename, int segnum){
	int i, j;
	ofstream transfile(filename);
	if(!transfile) {
			cout << "Fail to open distance inputfile! \n";
	}
	else{
		for (i = 0; i < segnum + 1; i++){
			for(j = 0; j < segnum; j++){
				transfile << trans[i][j] << " ";
			}
			if ( i != segnum ){
				transfile << endl;
			}
		}
	}
}

void testPrintSample(int samplenum, int segnum){
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
	
	for(i = 0; i < samplenum; i++){
		for (j = 0; j <segnum+1; j++){
			cout << segdiv[i][j] << " ";
		}
		cout << endl;
	}
}

void outputSample(int samplenum, int segnum, const char* filename){
    int i, j;
	ofstream file(filename);
	if(!file) {
			cout << "Fail to open sample output file! \n";
	}
	else{
		for(i = 0; i < samplenum; i++){
			for (j = 0; j <segnum+1; j++){
				file << segdiv[i][j] << " ";
			}
			file << endl;
		}
	}
}


void splitWords(){
	int i, j;
	int samplenum = samples.size();
	int segnum = means[0].size();
	segdiv.resize(samplenum);
	trans.resize(segnum+1);
	for(i = 0; i < segdiv.size(); i++){
		segdiv[i].resize(segnum+1);
	}
	for(i = 0; i < trans.size(); i++){
		trans[i].resize(segnum);
	}

	//initSeg(samplenum, segnum);

	for (j = 0; j < samplenum; j++){
		    //In each sample, doing sample distance calculation.
	        calSampDistance(j);
			//Update the kmeans distance.
            kmeansDistance();
			//Update segment division array.
	        kBackTrack(j);
			kdist.clear();
	}
	//testPrintSample(samplenum, segnum);
}

int transcript2int(string transcript){
	map<string, int> m;
	m[string("one")] = 1;
	m[string("two")] = 2;
	m[string("three")] = 3;
	m[string("four")] = 4;
	m[string("five")] = 5;
	m[string("six")] = 6;
	m[string("seven")] = 7;
	m[string("eight")] = 8;
	m[string("nine")] = 9;
	m[string("zero")] = 0;
	m[string("oh")] = 10;
	m[string("sil")] = -1;
	return (m.find(transcript)->second);
}

void getFileNameList(string filelist){
	vector<string> transcriptline;
	char charsToTrim[] = { '(', ')'};
	ifstream filelist_in(filelist.c_str());
	if (!filelist_in.is_open()) {
		cout << "Open file " << filelist << "failed\n";
	}
	else {
		string input;
		while (!filelist_in.eof()) {
			filelist_in >> input;
			if( input[0] == '(' ){
				input.erase(0,1);
				input.erase(input.find(')'));
				input = "train/" + input;
				transcriptline.push_back(input);
				transcripts.push_back(transcriptline);
				transcriptline.clear();
			}
			else{
				if (input != "sil"){
					transcriptline.push_back(input);
				}

			}
		}
	}
	filelist_in.close();
}

void segKMeansProcess(){
	//The training process of seg-kmeans.
	int i, j, k;
	int tsize, ssize, curnum;
	string msample[11] = {"0.txt","1.txt","2.txt","3.txt","4.txt","5.txt","6.txt","7.txt","8.txt","9.txt","o.txt"};
	string csample[11] = {"c-0.txt","c-1.txt","c-2.txt","c-3.txt","c-4.txt","c-5.txt","c-6.txt","c-7.txt","c-8.txt","c-9.txt","c-o.txt"};
	
	string filelist = "TRAIN.transcripts";
	getFileNameList(filelist);

	tsize = transcripts.size();

	for (i = 0; i < tsize; i++){
		ssize = transcripts[i].size();
		/*  Only for training isolated words.
			if (ssize >2) {
			cout << "(SKIP) |"; 
			continue;
			
		}*/
		for(j = 0; j < ssize - 1; j++){
			curnum = transcript2int(transcripts[i][j]);
			constructSamples(msample[curnum].c_str(), MEAN);
			constructSamples(csample[curnum].c_str(), COVAR);
		}
		readSampleDCT(transcripts[i][ssize-1].c_str());
		splitWords();
		for (j = 0; j < ssize - 1; j++){
			curnum = transcript2int(transcripts[i][j]);
			continue2isolate(0, curnum, j);
		}
		cout << i << ".. |";
		clearforLoop();
		samples.clear();
		segdiv.clear();
	}
}

void eachSeg(const char* meanfile, const char* covarfile, int segnum){
	int i, j;
	int samplenum = samples.size();
	segdiv.resize(samplenum);
	trans.resize(segnum+1);
	for(i = 0; i < segdiv.size(); i++){
		segdiv[i].resize(segnum+1);
	}
	for(i = 0; i < trans.size(); i++){
		trans[i].resize(segnum);
	}
	initSeg(samplenum, segnum);
	//Testing iteration process
	for (i = 0; i < 4; i++){
		//Clear vectors
		clearforLoop();
		//Init transition costs.
	    initTrans(segnum, samplenum);
		//Calculate means an covariances
        calSegmentMeansCovar(samplenum, segnum);
		for (j = 0; j < samplenum; j++){
			if(!samples[j].empty()){
				//In each sample, doing sample distance calculation.
				calSampDistance(j);
				//Update the kmeans distance.
				kmeansDistance();
				//Update segment division array.
				kBackTrack(j);
				kdist.clear();
			}
		}
		//testPrintSample(samplenum, segnum);
		if(i % 10 == 0){
			cout << i << ".. ";
		}
	}
	outputMean(meanfile);
	outputCovar(covarfile);
	string logname = meanfile;
	logname = "log-" + logname;
	outputSample(samplenum, segnum, logname.c_str());
	cout << meanfile << endl;
}


void assigndigits(){
	segdiv.assign(seg_all.s0.begin(), seg_all.s0.end());
	samples.assign(sample_all.s0.begin(), sample_all.s0.end());
	eachSeg("0.txt", "c-0.txt", STATE);

	segdiv.assign(seg_all.s1.begin(), seg_all.s1.end());
	samples.assign(sample_all.s1.begin(), sample_all.s1.end());
	eachSeg("1.txt", "c-1.txt", STATE);

	segdiv.assign(seg_all.s2.begin(), seg_all.s2.end());
	samples.assign(sample_all.s2.begin(), sample_all.s2.end());
	eachSeg("2.txt", "c-2.txt", STATE);

	segdiv.assign(seg_all.s3.begin(), seg_all.s3.end());
	samples.assign(sample_all.s3.begin(), sample_all.s3.end());
	eachSeg("3.txt", "c-3.txt", STATE);

	segdiv.assign(seg_all.s4.begin(), seg_all.s4.end());
	samples.assign(sample_all.s4.begin(), sample_all.s4.end());
	eachSeg("4.txt", "c-4.txt", STATE);

	segdiv.assign(seg_all.s5.begin(), seg_all.s5.end());
	samples.assign(sample_all.s5.begin(), sample_all.s5.end());
	eachSeg("5.txt", "c-5.txt", STATE);

	segdiv.assign(seg_all.s6.begin(), seg_all.s6.end());
	samples.assign(sample_all.s6.begin(), sample_all.s6.end());
	eachSeg("6.txt", "c-6.txt", STATE);

	segdiv.assign(seg_all.s7.begin(), seg_all.s7.end());
	samples.assign(sample_all.s7.begin(), sample_all.s7.end());
	eachSeg("7.txt", "c-7.txt", STATE);

	segdiv.assign(seg_all.s8.begin(), seg_all.s8.end());
	samples.assign(sample_all.s8.begin(), sample_all.s8.end());
	eachSeg("8.txt", "c-8.txt", STATE);

	segdiv.assign(seg_all.s9.begin(), seg_all.s9.end());
	samples.assign(sample_all.s9.begin(), sample_all.s9.end());
	eachSeg("9.txt", "c-9.txt", STATE);

	segdiv.assign(seg_all.so.begin(), seg_all.so.end());
	samples.assign(sample_all.so.begin(), sample_all.so.end());
	eachSeg("o.txt", "c-o.txt", STATE);
	cout << "end" ;
}


void clearall(){
	samples.clear();
	means.clear();
	covars.clear();
	kdist.clear();
	segdiv.clear();
	trans.clear();
	transcripts.clear();

	sample_all.s0.clear();
	seg_all.s0.clear();

	sample_all.s1.clear();
	seg_all.s1.clear();

	sample_all.s2.clear();
	seg_all.s2.clear();

	sample_all.s3.clear();
	seg_all.s3.clear();

	sample_all.s4.clear();
	seg_all.s4.clear();

	sample_all.s5.clear();
	seg_all.s5.clear();

	sample_all.s6.clear();
	seg_all.s6.clear();

	sample_all.s7.clear();
	seg_all.s7.clear();

	sample_all.s8.clear();
	seg_all.s8.clear();

	sample_all.s9.clear();
	seg_all.s9.clear();

	sample_all.so.clear();
	seg_all.so.clear();
}
