#ifndef CCV_NDFT_H
#define CCV_NDFT_H

#include "fftw3.h"
#include "fftwPrecision.h"
#include "math/SmoothingFunction.h"
#include "sparsefft3.h"


void gridding( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f, 
               double blobbiness, int n, int m, 
               bool smoothSkin, SmoothingFunction* smoothingFunction, 
               FFTW_complex* gHat, bool spreadSkin );

void griddingElec( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f, 
                   double blobbiness, int n, double elecRadiusInGrids, FFTW_DATA_TYPE* gHat, bool forInPlaceFFT, bool movingMol );

void griddingHbond( int M, double* x, double* y, double* z, float *r, double rExp, FFTW_complex* f, 
                    double blobbiness, int n, FFTW_complex* gHat, bool movingMol );

void griddingHydrophobicity( int M, double* x, double* y, double* z, float *r, FFTW_complex* f, 
                             double blobbiness, int n, FFTW_complex* gHat, double hydroRadExt, bool pairWise );

void griddingSimpleComplementarity( int M, double* x, double* y, double* z, float *r, FFTW_complex* f, 
                             double blobbiness, int n, FFTW_complex* gHat, double simpleRadExt );

bool getCenterFrequencies( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f,
                           double blobbiness, double alpha, int N, int m, 
			   FFTW_complex* hHat, bool smoothSkin, SmoothingFunction* smoothingFunction,
			   FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			   bool griddingDone, bool spreadSkin );

bool getCenterElecFrequencies( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f,
                               double blobbiness, double alpha, int N, double elecRadiusInGrids, 
			       FFTW_complex* hHat, FFTW_DATA_TYPE* gHat, FFTW_plan gHatPlan, bool griddingDone, bool movingMol );
			       
bool getCenterHbondFrequencies( int M, double* x, double* y, double* z, float *r, double rExp, FFTW_complex* f,
                                double blobbiness, double alpha, int N, 
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone, bool movingMol );

bool getCenterHydrophobicityFrequencies( int M, double* x, double* y, double* z, float *r, FFTW_complex* f,
                                double blobbiness, double alpha, int N, double hydroRadExt, bool pairWise,
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone );

bool getCenterSimpleComplementarityFrequencies( int M, double* x, double* y, double* z, float *r, FFTW_complex* f,
                                double blobbiness, double alpha, int N, double simpleRadExt,
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone );

#endif
