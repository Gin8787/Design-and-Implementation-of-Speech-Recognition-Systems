#include"readTest.h"

short* readAudioData(char *fileName, int *numSamples, int *sampleRate) {
	short *audioData = ReadWavFile(fileName, numSamples, sampleRate);
	return audioData;
}

string extractFeatureOfAudio(short *waveData, int nbSamples, int sampleRate) {
	string normalizedFileName;
	cout << "Please input the file name where you want to store the normalized feature data!" << endl;
	cin >> normalizedFileName;
	feature_extraction(waveData, nbSamples, sampleRate, normalizedFileName);
	return normalizedFileName;
}