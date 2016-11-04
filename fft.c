#include <stdio.h>
#include <stdlib.h> // malloc
#include <math.h>
#include <complex.h>

#define PI acos(-1.0) // even if it is defined as M_PI in math.h, this is more precise computationally
#define EVEN 0
#define ODD 1

typedef double complex cplx;

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
 * this algorithm involves splitting the array in two parts each recursion
 * to be more efficient 
 */
cplx *fast_ft(cplx *compArray, int len)
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
	
	evenA = fast_ft(split_array(compArray, len, EVEN), len/2);
	oddA = fast_ft(split_array(compArray, len, ODD), len/2);
	
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

int 
main(int argc, char *argv[])
{
	int i, inLen;
	cplx *inputComponents;
	cplx *outputComponents;

	inLen = argc - 1; // inLen must be a power of 2
	if(inLen % 2 != 0 || argc < 2){
		fprintf(stderr, "Usage: fft [array of input components]\nNote that the length of the array MUST be a power of 2.\nEXAMPLE:\n$ ./fft 1 1 1 1 0 0 0 0\n");
		exit(EXIT_FAILURE);
	}
	inputComponents = (cplx*)malloc((inLen)*sizeof(cplx));
	
	for(i=1; i<argc; i++){
		inputComponents[i-1] = atoi(argv[i]);
	}
	
	fprintf(stdout, "in:\n");
	print_components(inputComponents, inLen);
	
	outputComponents = fast_ft(inputComponents, inLen);

	fprintf(stdout, "out:\n");
	print_components(outputComponents, inLen);
	
	free(outputComponents);
	
	return 0;	
}
