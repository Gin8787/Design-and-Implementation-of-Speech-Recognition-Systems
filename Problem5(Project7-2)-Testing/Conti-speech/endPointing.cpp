#include"endPointing.h"

double forgetfactor = 1.1;
double level = 0;
double threshold = 12.5;
double background = 0;
double adjustment = 0.05;
unsigned long silenceFrames = 0;
bool isSpeechOver = false;
double silenceAverage = 0.0;
unsigned long sampleIndex = 0;
vector<short > audioNoSil;
//vector<short> forwardSil;
//vector<short> backwardSil;

double EnergyPerSampleDecibel(const short * frameStart,  int startIndex) {
	int i;
	double totalEnergy = 0;
	short temp;
	for (i = 0; i < FRAMES_PER_BUFFER; i++) {
		//		cout << "sample index is " << sampleIndex++ << endl;
		temp = frameStart[startIndex + i];
		//	if(i < 10) cout << "temp is" << temp << endl;
		totalEnergy += (double)temp * (double)temp;
	}
	totalEnergy = 10 * log(totalEnergy);
	return totalEnergy;
}

double EnergyPerSampleDecibelBack(const short * frameEnd, int startIndex) {
	int i;
	double totalEnergy = 0;
	short temp;
	for (i = 0; i < FRAMES_PER_BUFFER; i++) {
		//		cout << "sample index is " << sampleIndex++ << endl;
		temp = frameEnd[startIndex - i];
		//	if(i < 10) cout << "temp is" << temp << endl;
		totalEnergy += (double)temp * (double)temp;
	}
	totalEnergy = 10 * log(totalEnergy);
	return totalEnergy;
}

bool classifyFrame(const short * frameStart, bool isStart, int startIndex) {
	double current;
	if (isStart)
		current = EnergyPerSampleDecibel(frameStart, startIndex);
	else
		current = EnergyPerSampleDecibelBack(frameStart, startIndex);
	bool isSpeech = false;
	level = ((level * forgetfactor) + current) / (forgetfactor + 1);
	if (current < background)
		background = current;
	else
		background += (current - background) * adjustment;

	if (level < background) level = background;
	if ((level - background) > threshold) isSpeech = true;
	//	cout << "level - background is " << level - background << endl;
	silenceAverage += (level - background);
	return isSpeech;
}

/* Use the firt 16000 samples to initial background*/
void backgroundInitial(const short * firstFrame) {
	int i;
	short temp;
	for (i = 0; i < (FRAMES_PER_BUFFER * 2); i++) {
		//cout << "i is " << i << endl;
		temp = *firstFrame++;
		background += (double)temp * (double)temp;
	}
	background = 10 * log(background);
}

void backgroundInitBackward(const short *frameStart, int startIndex) {
	int i;
	short temp;
	for (i = 0; i < (FRAMES_PER_BUFFER * 2); i++) {
		//cout << "i is " << i << endl;
		temp = frameStart[startIndex - i];
		background += (double)temp * (double)temp;
	}
	background = 10 * log(background);
}

short* cutSil(string audioFile, string newAudioFile, int *numSamples, int *sampleRate) {
	int i, speechStartIndex, speechEndIndex;
	short *wavData = ReadWavFile(audioFile.c_str(), numSamples, sampleRate);
	bool isStart = true;
	backgroundInitial(wavData);
	for (i = 0; i < *numSamples; i+=400) {
		if (classifyFrame(wavData, isStart, i)) {
			speechStartIndex = i;
			break;
		}
	}
/*	for (i = 0; i < numSamples; i++) {
		cout << "i is " << i << " sample is " << hex << wavData[i] << endl;
		if (i == 7167) {
			cout << "hehe\n";
		}
	}*/
	isStart = false;
//	backgroundInitial(wavData);
	backgroundInitBackward(wavData, *numSamples - 1);
	level = 0;
//	short *wavDataEnd = wavData + numSamples - 1;
	for (i = *numSamples - 1; i >= 0; i-=400) {
		if (classifyFrame(wavData, isStart, i)){
			speechEndIndex = i;
			break;
		}
	}
/*	for (i = 0; i < speechStartIndex; i++) {
		forwardSil.push_back(wavData[i]);
	}*/
	for (i = speechStartIndex; i < speechEndIndex; i++) {
		audioNoSil.push_back(wavData[i]);
	}
/*	for (i = speechEndIndex; i < *numSamples; i++) {
		backwardSil.push_back(wavData[i]);
	}*/
	int audioLen = audioNoSil.size();
	short *audio = new short[audioLen];
	for (i = 0; i < audioLen; i++){
		audio[i] = audioNoSil[i];
	}
/*	int silForLen = forwardSil.size();
	short *forSil = new short[silForLen];
	for (i = 0; i < speechStartIndex; i++) {
		forSil[i] = forwardSil[i];
	}*/
//	int forwardLen = 
	if (WriteWave(newAudioFile.c_str(), audio, speechEndIndex - speechStartIndex, *sampleRate))
		printf("File is saved!");
	else
		printf("Fail to save the file");
//	delete audio;
/*	if (WriteWave(newAudioFile.c_str(), forSil, speechStartIndex - 1, sampleRate))
		printf("File is saved!");
	else
		printf("Fail to save the file");*/
//	delete forSil;
	audioNoSil.clear();
	return audio;
}