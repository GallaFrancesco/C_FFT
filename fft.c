#include "fft.h"

/* if flag is EVEN (0), it takes only the even elements
 * otherwise if flag is ODD (1) it takes only the odd ones
 */
cplx *split_array(cplx *a, int len, int flag)
{
	int i, cnt = 0;
	cplx *ret = malloc((len/2)*sizeof(cplx));

	for(i=0+flag; i<len; i=i+2){
		ret[cnt] = a[i];
		cnt++;
	}
	return ret;
}

/* recursively compute the fft on an array of complex numbers
 * splitting the array in two parts each recursion
 */
cplx *_fast_ft(cplx *compArray, int len)
{
	cplx omegaN, omega;
	cplx *evenA, *oddA, *transformedA;
	int i;

	/*termination*/
	if(len == 1){
		return compArray;
	}

	omega = 1;
	omegaN = cexp(2*PI*I/len); //the fourier coefficient
    evenA = _fast_ft(split_array(compArray, len, EVEN), len/2);
	oddA = _fast_ft(split_array(compArray, len, ODD), len/2);

	/*the final array*/
	transformedA = malloc(len*sizeof(cplx));

	for(i=0; i<(len/2); i++){
		transformedA[i] = evenA[i] + omega*oddA[i];
		transformedA[i+(len/2)] = evenA[i] - omega*oddA[i];
		omega = omegaN*omega;
	}
	free(evenA);
	free(oddA);
	free(compArray);
	return transformedA;
}

void
print_components(cplx *a, int len)
{
	int i;
	for(i=0; i<len; i++){
		/*creal and cimag extract the real and imaginary parts of a[i]*/
		fprintf(stdout, "%g, %g\n", creal(a[i]), cimag(a[i])); 
	}
	fprintf(stdout, "\n");
}

unsigned int
amplitude(cplx c, unsigned int n)
{
    // compute a normalized amplitude (actually power spectrum, since it's not squared
    // normalize on n (N_SAMPLES)
	double sq = 0;
	unsigned int res = 0;

	/*compute amplitude*/
    sq = sqrt(pow(creal(c)/n, 2) + pow(cimag(c)/n, 2));
	res = round(20*log10(sq)); // dB scale
    /*fprintf(stderr, "%f - %d\n", sq, res);*/
	return res;
}

void
normalize_fft(int inLen, unsigned int* fftSig)
{
	// normalize by dividing for the amplitude
	int i;
	for(i=0; i<inLen; i++){
		fftSig[i] = (int)(10*fftSig[i]/(20*log10(fftSig[512])));
		fprintf(stderr, "%d\n", fftSig[i]);
	}
}

void
fast_fft(int inLen, uint16_t *sig, unsigned int *fftSig)
{
	int i;
	cplx *inputComponents;
	cplx *outputComponents;

	if(inLen % 2 != 0){
		fprintf(stderr, "Note that the length of the array MUST be a power of 2.");
		exit(EXIT_FAILURE);
	}
	inputComponents = (cplx*)malloc((inLen)*sizeof(cplx));

	for(i=1; i<inLen; i++){
		inputComponents[i] = sig[i];
	}

	outputComponents = _fast_ft(inputComponents, inLen);
	for(i=0; i<inLen; i++){
		fftSig[i] = amplitude(outputComponents[i], inLen);
	}
	/*normalize_fft(inLen, fftSig);*/
	free(outputComponents);
}

void
average_signal(unsigned int *fftBuf, int inLen, int max, unsigned int* fftAvg)
{
	int i, j, step, k=0;
	unsigned int avg;
	int maxfreq = 0;

	int NADJ=4;
	int boost=8;
	step=1;

	for(i=0; i<inLen/2; i=i+2*NADJ) {
		if(i - NADJ >= 0 && i+NADJ < k && \
				fftBuf[i] >= fftBuf[i-NADJ] + 2 && \
				fftBuf[i] >= fftBuf[i+NADJ] + 2) {
			fftBuf[i] += boost;
			for(j=1; j<NADJ; ++j) {
				fftBuf[i-j] += boost/(boost-step);
				fftBuf[i+j] = fftBuf[i-j];
				step += boost/NADJ;
				fftBuf[i+inLen/2+j] = fftBuf[i+j];
				fftBuf[i+inLen/2-j] = fftBuf[i-j];
			}
		}
	}

    // N_SAMPLES / maximum number of columns
    step = inLen*FOCUS/max;
	for(i=0; i<inLen; i=i+step){
		maxfreq = 0;
		for(j=0; j<step; j++){
			if (fftBuf[i+j] > maxfreq) {
				maxfreq = fftBuf[i+j];
			}
		}
		fftAvg[k] = maxfreq-step/FOCUS; //the 80 is a correction for the display

		for(j=1; j<FOCUS; ++j) {
			fftAvg[k+j] = fftAvg[k]; //the 80 is a correction for the display
		}
		k += FOCUS;
	}

}

void
energy_sub_signal(unsigned int *fftBuf, int inLen, int outLen, unsigned int* fftSEnergy)
{
    // given sound amplitudes in a buffer, divide it into `outLen` subbands,
    // compute the energy of the subband
    // store it in fftEnergy
    int i, k;
    double val;

    for(i=0; i<outLen; i++)
    {
        fftSEnergy[i] = 0;
        for(k=outLen*i; k<outLen*(i+1); k++) {
            val = outLen*fftBuf[k]/inLen;
            fftSEnergy[i] += (int)(val);
        }
    }
}
