#pragma once
#include <stdio.h>
#include <stdlib.h> // malloc
#include <math.h>
#include <complex.h>
#include <stdint.h>		// uint16_t

#define PI acos(-1.0) // even if it is defined as M_PI in math.h, this is more precise computationally
#define EVEN 0
#define ODD 1
#define FOCUS 2

typedef double complex cplx;

cplx *split_array(cplx *a, int len, int flag);
cplx *_fast_ft(cplx *compArray, int len);
void print_components(cplx *a, int len);
unsigned int amplitude(cplx c, unsigned int n);
void fast_fft(int inLen, uint16_t *sig, unsigned int* fftSig);
void average_signal(unsigned int *fftBuf, int inLen, int maxC, unsigned int *fftAvg);
void energy_sub_signal(unsigned int *fftBuf, int inLen, int outLen, unsigned int* fftSEnergy);
