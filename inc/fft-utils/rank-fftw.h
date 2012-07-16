#ifndef RANKFFTW_H_
#define RANKFFTW_H_

#include <stdio.h>
#include <stdlib.h>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include "fftw3.h"
#include "fftwPrecision.h"
#include "../utils/utils.h"

#define NO_ERROR 0
#define INVALID_PARAMETER -1
#define MEMORY_ALLOCATION_FAILED -2
#define FFTW_FORWARD_FAILED -3
#define FFTW_BACKWARD_FAILED -4
#define FFTW_FORWARD_WISDOM_FAILED -5
#define FFTW_BACKWARD_WISDOM_FAILED -6
#define WISDOM_EXPORT_FILE_OPEN_FAILED -7
#define WISDOM_IMPORT_FILE_OPEN_FAILED -8

#ifndef MIN_SIZE
   #define MIN_SIZE 2
#endif

#ifndef MAX_SIZE
   #define MAX_SIZE 256
#endif

#ifndef IN_PLACE
   #define IN_PLACE true
#endif

#ifndef MAX_ITER
   #define MAX_ITER 3
#endif

#ifndef WISDOM_FILE
   #define WISDOM_FILE "wisdom.txt"
#endif

int computeEffGrid( int minSize, int maxSize );

#endif /* RANKFFTW_H_ */
