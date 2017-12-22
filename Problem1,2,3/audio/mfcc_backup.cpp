#include "mfcc.h"

const double PREEMP_COEF = 0.95;
const double PI = 3.1416;
const double epsilon = 0.00001;
const int FFT_LENGTH = 512;

vector< vector<double> > normalizedMel(13);

void feature_extraction(short *waveData, int nbSamples, int sampleRate)
{
	printf("-------------------In MFCC Now------------------\n");
	int i, j;
	int nbFrames = 0;
	int countFrame = 0;

	int countLength = FRAME_LENGTH;
	double frameData[FFT_LENGTH] = {0};

	double hamming_win[FRAME_LENGTH] = {0};

	double energy[FFT_LENGTH] = {0};

	double mel_spec[FILTERS] = {0};
	double mel_cept[MFCC_COEF] = {0};

	double DCTSum[13] = {0.0};              //Sum for each feature point from every frame ater DCT transiformation
	vector<double> tempNormalize(1);
	ofstream normalization("normalize.txt");

	printf("-------------------Open Files Now------------------\n");

	FILE *original = fopen("original.txt", "w");
//	FILE *preemp = fopen("preemp.txt", "w");
//	FILE *hamming = fopen("hamming.txt", "w");
//	FILE *fftin = fopen("fftin.txt", "w");
//	FILE *fftout = fopen("fftout.txt", "w");
	FILE *energyfile = fopen("energy.txt", "w");
//	FILE *melbound = fopen("melb.txt", "w");
	FILE *mel = fopen("mel.txt", "w");
	FILE *melspecfile = fopen("melspecfile.txt", "w");
//	FILE *normalization = fopen("normalization.txt", "w");
	

	printf("-------------------Open File Succeed------------------\n");

	/*Test Output for original data*/
	for (i = 0; i < FRAME_LENGTH; i++)
	{
		fprintf(original,"%d ",waveData[i]);
	}
	/*End Test Output for original data*/

	/**Pre-emphasize**/
	double *pre = new double[nbSamples];
	pre[0] = waveData[0];
	for (i = 1; i < nbSamples; i++)
	{
		pre[i] = waveData[i] - PREEMP_COEF * waveData[i-1];
		//printf("%d ",pre[i]);
	}
	/*Test Output for preemp data*
	for (i = 0; i < FRAME_LENGTH; i++)
	{
		fprintf(preemp,"%f ",pre[i]);
	}
	/*End Test Output for preemp data*/

	/**Calculate Frames**/
	nbFrames = ceil((double)((nbSamples - FRAME_LENGTH) / (FRAME_LENGTH - OVERLAP))) + 1; //Quest?

	/**Define Hamming Window**/
	for (i = 0; i < FRAME_LENGTH; i++)
	{
		hamming_win[i] = 0.54 - 0.46 * cos(2 * PI * i / (FRAME_LENGTH - 1));
	}

	/*Test Output for hamming data*
	for (i = 0; i < FRAME_LENGTH; i++)
	{
		fprintf(hamming,"%f ",hamming_win[i]);
	}
	/*End Test Output for hamming data*/

	/** Define and Calculate Mel-Filter Boundaries **/
	//The mel scale frequency can be choose wisely by ourselves, set as 8000Hz in our exp. 
	//The requirement is 6855.76Hz. 
	double maxMelfreq = 1127 * log( (double) (1 + ( MAXFREQ / 700) ));
	double perMelfreq = maxMelfreq / (FILTERS + 1);
	//Define the boundaries of Melfilters in Melscale.
	//Transfer MelScale frequency to normal scale frequency.
	//Transfer frequency to dots on Frequency Domain.(QUESTION: 
	FilterMelFreq filter[FILTERS];
	for (i = 0; i < FILTERS; i++)
	{
		//Get rid of floor option
		filter[i].lower_bound = perMelfreq * i;
		filter[i].lower_bound = 700 * (exp(filter[i].lower_bound / 1127) - 1);
		filter[i].lower_bound = ((FFT_LENGTH / 2 + 1) * filter[i].lower_bound / sampleRate);

		filter[i].center = perMelfreq * (i+1);
		filter[i].center = 700 * (exp(filter[i].center / 1127) - 1);
		filter[i].center = ((FFT_LENGTH / 2 + 1) * filter[i].center / sampleRate);

		filter[i].higher_bound = perMelfreq * (i+2);
		filter[i].higher_bound = 700 * (exp(filter[i].higher_bound / 1127) - 1);
		filter[i].higher_bound = ((FFT_LENGTH / 2 + 1) * filter[i].higher_bound / sampleRate);
		//fprintf(melbound,"%f\n",filter[i].higher_bound);
	}
	/** End of Defining and Calculating Mel-Filter Boundaries**/

	printf("-------------------Initial every fileter succeed------------------\n");

	printf("-------------------handle every frame------------------\n");
	while (countFrame < nbFrames){
		/**Operation Inside Each Frames**/

		if (countFrame == (nbFrames - 1)) 
			countLength = (nbSamples - (countFrame * (FRAME_LENGTH - OVERLAP)));
		else 
			countLength = FRAME_LENGTH;

		for (i = 0; i < countLength; i++)
		{
			frameData[i] = pre[countFrame * (FRAME_LENGTH - OVERLAP) + i] * hamming_win[i];
		}

		//printf("%d, %d, %d",countFrame, countLength, countFrame * (FRAME_LENGTH - OVERLAP));

		/**Zero Padding**/
		for (i = countLength; i < FFT_LENGTH; i++)
		{
			frameData[i] = 0;
		}

		/**FFT Transform**/
		fftw_complex *in,*out;
		fftw_plan p;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_LENGTH);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_LENGTH);

		if((in == NULL)||(out == NULL)){
			printf("Error:insufficient available memory\n");
		}
		else{
			for(i = 0; i < FFT_LENGTH; i++)/* Testing data */
			{
				in[i][0] = frameData[i];
				in[i][1] = 0;
			}
		}
		p = fftw_plan_dft_1d(FFT_LENGTH, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p); /* repeat as needed */
		fftw_destroy_plan(p);
		fftw_cleanup();

		/**Test Input and output of FFT transform**
		for(i = 0; i < FFT_LENGTH; i++) //Test INPUT
		{
			fprintf(fftin,"%f ",in[i][0]);
		}
		for(i = 0; i < FFT_LENGTH; i++)//Test OUTPUT
		{
			fprintf(fftout,"%f,%fi ",out[i][0],out[i][1]);
		}
		/**End Test Input and output of FFT transform**/

		/**Calculate the energy distribution**/
		for (i = 0; i < (FFT_LENGTH / 2); i++)
		{
			energy[i] = (out[i][0] * out[i][0]) + (out[i][1] * out[i][1]);
		}

		for(i = 0; i < ( (FFT_LENGTH / 2) - 1); i++)//Energy OUTPUT to file
		{
			fprintf(energyfile,"%f,",log10(energy[i]));
		}
		fprintf(energyfile,"%f\n",log10(energy[i]));

		if(in != NULL) fftw_free(in);
		if(out != NULL) fftw_free(out);
		/**End of FFT Transform**/

		/** Melfilter and energy add-up **/

		for (i = 0; i < FILTERS; i++)
		{
			mel_spec[i] = 0;
			for (j = 0; j < (FFT_LENGTH / 2); j++)  //Question: Only half of the dots on Freq Domain was calculated?
			{
				if(j <= filter[i].center && j>= filter[i].lower_bound)
				{
					mel_spec[i] += energy[j] * ((j - filter[i].lower_bound) / (filter[i].center - filter[i].lower_bound));
				}
				else if(j > filter[i].center && j<= filter[i].higher_bound)
				{
					mel_spec[i] += energy[j] * ((filter[i].higher_bound - j) / (filter[i].higher_bound - filter[i].center));
				}
				else
				{
					mel_spec[i] += 0;
				}
			}
		}

		/** End of melfilter and energy addup**/

		for(i = 0; i < (FILTERS - 1); i++)//Test OUTPUT
		{
			fprintf(melspecfile,"%f,",log(mel_spec[i]));
		}
		fprintf(melspecfile,"%f\n",log(mel_spec[i]));

		/** DCT Transform **/

		for (i = 0; i < MFCC_COEF; i++)
		{
			mel_cept[i] = 0;
			for(int j = 0; j < FILTERS; j++)
			{
				if(mel_spec[j] <= epsilon || mel_spec[j] >= epsilon)

					mel_cept[i] += cos(i * PI / FILTERS * (j + 0.5)) * log(mel_spec[j]); 

			}
		}
        
		for(i = 0; i < (MFCC_COEF - 1); i++)//Test OUTPUT
		{
			normalizedMel[i].push_back(mel_cept[i]);
			DCTSum[i] += mel_cept[i];
			fprintf(mel,"%f,",mel_cept[i]);
		}
		normalizedMel[i].push_back(mel_cept[i]);  //Store very feature point to normalizedMel
		DCTSum[i] += mel_cept[i];
		fprintf(mel,"%f\n",mel_cept[i]);

		/**End of DCT Transform**/
		countFrame ++;
		//break;
	}

	printf("-------------------Every frame handle finished------------------\n");
	printf("countFrame is %d\n", countFrame);
	printf("Normalize column is %lu\n", normalizedMel[0].size());
	printf("Normalize row is %lu\n", normalizedMel.size());
	printf("-------------------Normalization now------------------\n");
	/*Normalization*/
	for(i = 0; i < MFCC_COEF; i++) {
//		printf("i is %d \n", i);
		for(j = 0; j < nbFrames; j++) {
			normalizedMel[i][j] = normalizedMel[i][j] - (DCTSum[i] / (double)countFrame);
		}
	}

	/** End of Nomalization **/

	printf("-------------------normalization finished------------------\n");

	/*Write normalizationMel to file*/

	if(!normalization) {
		cout << "fail to open normalization.txt \n";
		exit(0);
	}
	printf("--------------Open normalization succeed---------------\n");
	for(i = 0; i < nbFrames; i++) {
		for(j = 0; j < (MFCC_COEF - 1); j++) {
			normalization << normalizedMel[j][i] << ',';
		}
		normalization << normalizedMel[j][i] << '\n';
	}
	/*end of writing*/

	/*delete normalizedMel vector*/
	for(i = 0; i < MFCC_COEF; i++)
		normalizedMel[i].clear();

	delete[] pre;

//	fclose(fftin);
//	fclose(fftout);
	fclose(original);
//	fclose(preemp);
	fclose(energyfile);
//	fclose(hamming);
	fclose(mel);
	fclose(melspecfile);
//	fclose(melbound);

	printf("-------------------\nEnd of mfcc\n\n---------------");
}