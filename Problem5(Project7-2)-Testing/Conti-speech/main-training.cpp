#include "mfcc.h"
//#include "readwave.h"
#include "endPointing.h"
#include <iomanip>

vector<string > fileNameList;

void mfccBatch(bool isTrain);
void getFileNameList(string filelist);

int main(){
	mfccBatch(0);
//	cutSil("word/zero/MPB_ZB.wav", "word/zero/MPB_ZB_new.wav");
/*	string wavName = "FAC_OA_new.wav";
	string normalizeName = "FAC_OA_new.txt";
	int numSamples, sampleRate;
	short* rawData = NULL;
	rawData = ReadWavFile(wavName.c_str(), &numSamples, &sampleRate);
	feature_extraction(rawData, numSamples, sampleRate, normalizeName);*/
	cout << "main going to end" << endl;
	return 0;
}

void mfccBatch(bool isTrain) {
	string path = "hwdata/";
	string storePath;
	string dataFile, wavName, normalizeName, wavNameNew;
	string filelist;
	int i, fileNum, numSamples, sampleRate;
	if (isTrain) {
		dataFile = "train/";
		filelist = "TRAIN.filelist";
		storePath = path + "train_mfcc/";
	}
	else {
		dataFile = "test/";
		filelist = "TEST.filelist";
		storePath = path + "test_mfcc/";
	}
	filelist = path + filelist;
	getFileNameList(filelist);
	fileNum = fileNameList.size();
	for (i = 0; i < fileNum; i++) {
		wavName = path + dataFile + fileNameList[i] + ".wav";
		wavNameNew = path + dataFile + fileNameList[i] + "_new.wav";
		normalizeName = storePath + fileNameList[i] + "_new.txt";
		cout << "normal file name is " << normalizeName << endl;
		short* rawData = NULL;
		rawData = cutSil(wavName.c_str(), wavNameNew.c_str(), &numSamples, &sampleRate);
		feature_extraction(rawData, numSamples, sampleRate, normalizeName);
		delete rawData;
	}
	cout << "mfccBatch end" << endl;
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