/*
   Copyright 2011 The University of Texas at Austin

Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

This file is part of F2Dock.

F2Dock is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License version 2.1 as published by the Free Software Foundation.

F2Dock is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */


#include "Docking.h"

#include "TopValues.h"

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

#include <time.h>
#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#ifndef M_PI
#define 	M_PI   3.14159265358979323846
#endif

// variables for the rotation server for pthreads
pthread_mutex_t rotLock;
int curRotation = 0, maxRotation = 0, processedRotations = 0;

// priority queue locks for filtering
pthread_mutex_t globalInLock, globalOutLock;

//int *curRotation, *maxRotation, processedRotations = 0;

extern "C"
{
	void test1(int i)
	{
		printf(" GOT: %d\n", i);
	}

	void test2(float *d, int nb)
	{
		int i;
		for (i=0; i<nb; i++)
			printf(" GOT: %d %f\n", i, d[i]);
	}

	typedef struct
	{
		int i;
		float f;
		float *vect;
	} T1;

	void test3(T1 *t1)
	{
		printf(" GOT: %d %f\n", t1->i, t1->f);
		int i;
		printf(" pointer %p\n", t1->vect);
		if (t1->vect!=NULL)
			for (i=0; i<10; i++)
				printf(" Vect: %d %f\n", i, t1->vect[i]);
	}

	typedef struct
	{
		int breakDownScores;
		int numberOfPositions;
		int gridSize;
		int numFreq;
		int interpFuncExtent;
	} IN_PARAMS;

	int test4(IN_PARAMS *t1)
	{
		printf(" GOT: %d %d %d %d %d\n", t1->breakDownScores,
				t1->numberOfPositions, t1->gridSize, t1->numFreq,
				t1->interpFuncExtent);
		return 1;
	}

}

/**************************************************************************/
/*                                                                        */
/*  Get the current time in seconds as a double value                     */
/*                                                                        */
/**************************************************************************/
/*
   double getTime()
   {
#ifdef _WIN32
time_t ltime;
_timeb tstruct;
time( &ltime );
_ftime( &tstruct );
return (double) (ltime + 1e-3*(tstruct.millitm));
#else
struct timeval t;
gettimeofday( &t, NULL );
return (double)(t.tv_sec + 1e-6*t.tv_usec);
#endif
}
 */

// we need to add neg of the dot product ( since + + should repel etc )
void updateSCwithElec( FFTW_complex* sparseProfile, FFTW_complex* sparseElecProfile, int numFreq, double scale )
{

	int c;
	for ( c = 0; c < numFreq * numFreq * numFreq; c++ )
	{
		sparseProfile[ c ][ 0 ] -= sparseElecProfile[ c ][ 0 ] * scale;
		/*		if( scale < -0.01 )
				{
				if( sparseProfile[c][0] > 0 ) sparseProfile[c][0] = -sparseElecProfile[c][0];
				}
				else
				{
				if( sparseProfile[c][0] > 0 ) sparseProfile[c][0] -= sparseElecProfile[c][0]*scale;
				}
		 */
	}

}

bool compute3DFFT( FFTW_complex* dest, FFTW_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags )
{
	if( !dest || !src || dim1<1 || dim2<1 || dim3<1 ) return false;

	FFTW_plan p;
	p = FFTW_plan_dft_3d(dim1, dim2, dim3, src, dest, type, flags);
	FFTW_execute(p);

	FFTW_destroy_plan(p);
	return true;
}


bool fftshift3D( FFTW_complex* data, int length )
{
	if( !data || (length<1)) return false;
	FFTW_complex* temp = (FFTW_complex*)malloc( sizeof(FFTW_complex) * length * length * length );

	int c;
	for( c=0; c< length * length * length; c++ )
	{
		temp[c][0] = data[c][0];
		temp[c][1] = data[c][1];
	}

	int i, j, k;
	for( i=0; i<length; i++ )
	{
		for( j=0; j<length; j++ )
		{
			for( k=0; k<length; k++ )
			{
				int shift_i = (i+length/2)%length;
				int shift_j = (j+length/2)%length;
				int shift_k = (k+length/2)%length;

				data[shift_i*length*length + shift_j*length + shift_k][0] = temp[i*length*length + j*length + k][0];
				data[shift_i*length*length + shift_j*length + shift_k][1] = temp[i*length*length + j*length + k][1];
			}
		}
	}
	free( temp);
	return true;
}



/*******************************************************************/
/*                                                                 */
/*  The fft of a truncated gaussian is returned in freqHat         */
/*                                                                 */
/*******************************************************************/
void fftTruncHat( FFTW_complex* freqHat, int N, double blobbyness, double scale )
{
	int i;
	FFTW_plan p;

	int width = N/2;
	//double sigma = ((double)width) / 3.0;
	double sigma = 1.0 / (sqrt(2*-blobbyness));

	for (i=0;i<N;i++)
	{
		freqHat[i][0] = 0;
		freqHat[i][1] = 0;
	}

	for (i=0;i<(width)+1;i++)
	{
		freqHat[i][0] = exp(-i*i/(2.0*sigma*sigma*scale*scale));
	}

	for (i=N-1;i>(N-1)-(width); i-- )
	{
		freqHat[i][0] = exp(-(i-N)*(i-N)/(2.0*sigma*sigma*scale*scale));
	}
	//for (i=0;i<N;i++)
	//printf("[%d] = %f\n", i, freqHat[i][0] );
	p = FFTW_plan_dft_1d( N, freqHat, freqHat, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );
	FFTW_execute( p );
	FFTW_destroy_plan(p);
}


void fftTruncHat( FFTW_complex* freqHat, FFTW_plan freqHatPlan, int N, double blobbyness, double scale )
{
	int i;
	//	FFTW_plan p;

	int width = N/2;
	//double sigma = ((double)width) / 3.0;
	double sigma = 1.0 / (sqrt(2*-blobbyness));

	for (i=0;i<N;i++)
	{
		freqHat[i][0] = 0;
		freqHat[i][1] = 0;
	}

	for (i=0;i<(width)+1;i++)
	{
		freqHat[i][0] = exp(-i*i/(2.0*sigma*sigma*scale*scale));
	}

	for (i=N-1;i>(N-1)-(width); i-- )
	{
		freqHat[i][0] = exp(-(i-N)*(i-N)/(2.0*sigma*sigma*scale*scale));
	}
	//for (i=0;i<N;i++)
	//printf("[%d] = %f\n", i, freqHat[i][0] );
	//	p = FFTW_plan_dft_1d( N, freqHat, freqHat, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );
	FFTW_execute( freqHatPlan );
	//	FFTW_destroy_plan(p);
}


double getDistFromCenter( int gridPoint, int size, double scale )
{
	return ( size / 2.0 - gridPoint ) * scale;
}

/*************************************************************************/
/*                                                                       */
/*  The fft of the elec kernel is returned in 'kernel'                   */
/*                                                                       */
/*************************************************************************/
bool computeElecKernel_bak( FFTW_complex* kernel, int size, double scale )
{
	double D = 2.0;

	for ( int i = 0, c = 0; i < size; i++ )
	{
		double iDist = getDistFromCenter( i, size - 1, 1.0 / ( size - 1 ) );
		for ( int j = 0; j < size; j++ )
		{
			double jDist = getDistFromCenter( j, size - 1, 1.0 / ( size - 1 ) );
			for ( int k = 0; k < size; k++, c++ )
			{
				double kDist = getDistFromCenter( k, size - 1, 1.0 / ( size - 1 ) );
				double dist2 = iDist * iDist + jDist * jDist + kDist * kDist;
				double dist = sqrt( dist2 ) / scale;
				double distA = sqrt( dist2 + D * D * exp( - dist2 / ( 4 * D * D ) ) ) / scale;

				if ( distA < 1.0 ) kernel[ c ][ 0 ] = 1;
				else kernel[ c ][ 0 ] = 1.0 / distA;

				if ( dist <= 6.0 ) kernel[ c ][ 0 ] /= 4.0;
				else if ( dist <= 8.0 ) kernel[ c ][ 0 ] /= ( 38 * dist - 224 );
				else kernel[ c ][ 0 ] /= 80.0;

				kernel[ c ][ 1 ] = 0;
			}
		}
	}

	if ( !fftshift3D( kernel, size ) ) return false;
	if ( !compute3DFFT( kernel, kernel, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE ) ) return false;

	return true;
}




bool computeElecKernel_orig( FFTW_complex* kernel, int size, double scale )
{
	int i, j, k, c=0;

	for( i=0; i<size; i++ )
	{
		double iDist = getDistFromCenter( i, size - 1, 1.0 / (size - 1) /*scale*/ );
		for( j=0; j<size; j++ )
		{
			double jDist = getDistFromCenter( j, size - 1, 1.0 / (size - 1) /*scale*/ );
			for( k=0; k<size; k++ )
			{
				double kDist = getDistFromCenter( k, size - 1, 1.0 / (size - 1) /*scale*/ );
				double dist = sqrt( iDist*iDist + jDist*jDist + kDist*kDist );

				double radiusContribution = 0;
				double dielectricContribution = 0;

				{
					double distA = dist/scale;
					if( distA<1.0 )
						radiusContribution = 1;
					else
						radiusContribution = 1.0/distA;

					if( distA<=6.0 )
						dielectricContribution = 4;
					else if( distA < 8.0 )
						dielectricContribution = 38*distA-224;
					else
						dielectricContribution = 80;
				}

				kernel[c][0] = radiusContribution * ( 1.0/dielectricContribution );
				kernel[c][1] = 0;
				c++;
			}
		}
	}


	if( !fftshift3D( kernel, size ) ) return false;
	//printDataChars( "elec kernel after shift", kernel, size, false );
	if( !compute3DFFT( kernel, kernel, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE ) ) return false;
	//printDataChars( "elec kernel after fft", kernel, size, false );
	return true;
}



bool computeElecKernel( FFTW_complex* kernel, int size, double scale, double distVoid, double ldist, double lval, double hdist, double hval )
{
	double c1 = ( hval - lval ) / ( hdist - ldist );
	double c2 = lval - ldist * c1;

	for ( int i = 0, c = 0; i < size; i++ )
	{
		double iDist = getDistFromCenter( i, size - 1, 1.0 / ( size - 1 ) );

		for ( int j = 0; j < size; j++ )
		{
			double jDist = getDistFromCenter( j, size - 1, 1.0 / ( size - 1 ) );

			for ( int k = 0; k < size; k++, c++ )
			{
				double kDist = getDistFromCenter( k, size - 1, 1.0 / ( size - 1 ) );
				double dist = sqrt( iDist * iDist + jDist * jDist + kDist * kDist );

				double distA = dist / scale;

				if ( distA >= distVoid )
				{
					double radiusContribution, dielectricContribution;

					if ( distA < 1.0 ) radiusContribution = 1;
					else radiusContribution = 1.0 / distA;

					if ( distA <= ldist ) dielectricContribution = lval;
					else if ( distA < hdist ) dielectricContribution = c1 * distA + c2;
					else dielectricContribution = hval;

					kernel[ c ][ 0 ] = radiusContribution * ( 1.0 / dielectricContribution );
				}
				else kernel[ c ][ 0 ] = 0;

				kernel[ c ][ 1 ] = 0;
			}
		}
	}

	if ( !fftshift3D( kernel, size ) ) return false;

	if ( !compute3DFFT( kernel, kernel, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE ) ) return false;

	return true;
}



bool multiplyFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, int gridSize, double blobbiness, bool smoothSkin, FFTW_complex* frequenciesProduct, double scale)
{
	if ( !frequenciesA || !frequenciesB || !frequenciesProduct || ( numFreq < 1 ) || ( gridSize < 1 ) || ( numFreq > gridSize ) ) return false;

	if ( smoothSkin )
	{
		FFTW_complex* freqHat;
		freqHat = ( FFTW_complex* ) malloc( sizeof( FFTW_complex ) * gridSize );
		fftTruncHat( freqHat, gridSize, blobbiness, scale );

		for ( int k = 0, c = 0; k < numFreq; k++ )
		{
			int kIndex = k;
			if ( kIndex >= numFreq / 2 ) kIndex += ( gridSize - numFreq );

			for ( int j = 0; j < numFreq; j++ )
			{
				int jIndex = j;
				if ( jIndex >= numFreq / 2 ) jIndex += ( gridSize - numFreq );

				for ( int i = 0; i < numFreq; i++, c++ )
				{
					int iIndex = i;
					if ( iIndex >= numFreq / 2 ) iIndex += ( gridSize - numFreq );

					double val = freqHat[ iIndex ][ 0 ] * freqHat[ jIndex ][ 0 ] * freqHat[ kIndex ][ 0 ];

					frequenciesProduct[ c ][ 0 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ] ) * val * val;
					frequenciesProduct[ c ][ 1 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ] ) * val * val;
				}
			}
		}

		free( freqHat );
	}
	else
	{
		for ( int k = 0, c = 0; k < numFreq; k++ )
			for ( int j = 0; j < numFreq; j++ )
				for ( int i = 0; i < numFreq; i++, c++ )
				{
					frequenciesProduct[ c ][ 0 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ];
					frequenciesProduct[ c ][ 1 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ];
				}
	}

	return true;
}


bool multiplyFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, int gridSize, double blobbiness, bool smoothSkin, FFTW_complex* frequenciesProduct, FFTW_complex* freqHat, FFTW_plan freqHatPlan, double scale )
{
	if ( !frequenciesA || !frequenciesB || !frequenciesProduct || ( numFreq < 1 ) || ( gridSize < 1 ) || ( numFreq > gridSize ) ) return false;

	if ( smoothSkin )
	{
		fftTruncHat( freqHat, freqHatPlan, gridSize, blobbiness, scale );

		for ( int k = 0, c = 0; k < numFreq; k++ )
		{
			int kIndex = k;
			if ( kIndex >= numFreq / 2 ) kIndex += ( gridSize - numFreq );

			for ( int j = 0; j < numFreq; j++ )
			{
				int jIndex = j;
				if ( jIndex >= numFreq / 2 ) jIndex += ( gridSize - numFreq );

				for ( int i = 0; i < numFreq; i++, c++ )
				{
					int iIndex = i;
					if ( iIndex >= numFreq / 2 ) iIndex += ( gridSize - numFreq );

					double val = freqHat[ iIndex ][ 0 ] * freqHat[ jIndex ][ 0 ] * freqHat[ kIndex ][ 0 ];

					frequenciesProduct[ c ][ 0 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ] ) * val * val;
					frequenciesProduct[ c ][ 1 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ] ) * val * val;
				}
			}
		}
	}
	else
	{
		for ( int k = 0, c = 0; k < numFreq; k++ )
			for ( int j = 0; j < numFreq; j++ )
				for ( int i = 0; i < numFreq; i++, c++ )
				{
					frequenciesProduct[ c ][ 0 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ];
					frequenciesProduct[ c ][ 1 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ];
				}
	}

	return true;
}


bool multiplyElecFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, int kernelSize, FFTW_complex* frequenciesElecKernel, FFTW_complex* frequenciesProduct, double scale )
{
	if ( !frequenciesA || !frequenciesB || !frequenciesElecKernel || !frequenciesProduct || ( numFreq < 1 ) || ( kernelSize < 1 ) || ( numFreq > kernelSize ) ) return false;

	for ( int k = 0, c = 0; k < numFreq; k++ )
	{
		int kIndex = k;
		if ( kIndex >= numFreq / 2 ) kIndex += ( kernelSize - numFreq );

		for ( int j = 0; j < numFreq; j++ )
		{
			int jIndex = j;
			if ( jIndex >= numFreq / 2 ) jIndex += ( kernelSize - numFreq );

			for ( int i = 0; i < ( numFreq / 2 ) + 1; i++, c++ )
			{
				int iIndex = i;
				if ( iIndex >= numFreq / 2 ) iIndex += ( kernelSize - numFreq );

				double val = frequenciesElecKernel[ kIndex * numFreq * numFreq + jIndex * numFreq + iIndex ][ 0 ];

				frequenciesProduct[ c ][ 0 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ] ) * val;
				frequenciesProduct[ c ][ 1 ] = ( frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ] ) * val;
			}
		}
	}

	return true;
}


//bool multiplyHbondFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, FFTW_complex* frequenciesProduct )
//{
//	if ( !frequenciesA || !frequenciesB || !frequenciesProduct || ( numFreq < 1 ) ) return false;
//
//        for ( int k = 0, c = 0; k < numFreq; k++ )
//      	   for ( int j = 0; j < numFreq; j++ )
//      	      for ( int i = 0; i < numFreq; i++, c++ )
//      		 {
//      		   frequenciesProduct[ c ][ 0 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ];
//      		   frequenciesProduct[ c ][ 1 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ];
//      		 }
//
//	return true;
//}
//
//
//bool multiplyHydrophobicityFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, FFTW_complex* frequenciesProduct )
//{
//	if ( !frequenciesA || !frequenciesB || !frequenciesProduct || ( numFreq < 1 ) ) return false;
//
//        for ( int k = 0, c = 0; k < numFreq; k++ )
//      	   for ( int j = 0; j < numFreq; j++ )
//      	      for ( int i = 0; i < numFreq; i++, c++ )
//      		 {
//      		   frequenciesProduct[ c ][ 0 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ];
//      		   frequenciesProduct[ c ][ 1 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ];
//      		 }
//
//	return true;
//}


bool multiplyFrequencyMaps( FFTW_complex* frequenciesA, FFTW_complex* frequenciesB, int numFreq, FFTW_complex* frequenciesProduct )
{
	if ( !frequenciesA || !frequenciesB || !frequenciesProduct || ( numFreq < 1 ) ) return false;

	for ( int k = 0, c = 0; k < numFreq; k++ )
		for ( int j = 0; j < numFreq; j++ )
			for ( int i = 0; i < numFreq; i++, c++ )
			{
				frequenciesProduct[ c ][ 0 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 0 ] - frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 1 ];
				frequenciesProduct[ c ][ 1 ] = frequenciesA[ c ][ 0 ] * frequenciesB[ c ][ 1 ] + frequenciesA[ c ][ 1 ] * frequenciesB[ c ][ 0 ];
			}

	return true;
}


// FIXME .. what is that for ???
bool initializeFFTW(int size)
{
	if( size<1) return false;

	FFTW_complex* data = (FFTW_complex*)FFTW_malloc(sizeof(FFTW_complex)*size*size*size);
	{
		int i;
		for( i=0; i<size*size*size; i++ )
		{
			data[i][0] = (double)rand()/((double)(RAND_MAX));
			data[i][1] = (double)rand()/((double)(RAND_MAX));
		}
	}

	//double t = getTime();
	if( !compute3DFFT( data,data, size, size, size, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE ) ) return false;
	if( !compute3DFFT( data,data, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE ) ) return false;
	//printf("FFT Took: [%dx%dx%d] %lf\n", size, size, size, getTime()-t);

	FFTW_free( data );

	return true;
}


void initializeCurvatureGrid( PARAMS_IN *pr, int pw, double maxRatio, double upScale )
{
	printf( "initializing curvature grid " );
	fflush( stdout );

	double *x = pr->xkAOrig, *y = pr->ykAOrig, *z = pr->zkAOrig;
	float *r = pr->radiiA;
	char *type = pr->typeA;
	int M = pr->numCentersA;
	CURVATURE_GRID *cGrid = &( pr->cGrid );

	double minX, minY, minZ;
	double maxX, maxY, maxZ;
	double dimX, dimY, dimZ;

	minX = minY = minZ =  1000000000.0;
	maxX = maxY = maxZ = -1000000000.0;

	for ( int i = 0; i < M; i++ )
		if ( type[ i ] == 'I' )
		{
			if ( x[ i ] - r[ i ] < minX ) minX = x[ i ] - r[ i ];
			if ( y[ i ] - r[ i ] < minY ) minY = y[ i ] - r[ i ];
			if ( z[ i ] - r[ i ] < minZ ) minZ = z[ i ] - r[ i ];

			if ( x[ i ] + r[ i ] > maxX ) maxX = x[ i ] + r[ i ];
			if ( y[ i ] + r[ i ] > maxY ) maxY = y[ i ] + r[ i ];
			if ( z[ i ] + r[ i ] > maxZ ) maxZ = z[ i ] + r[ i ];
		}

	minX -= 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );
	minY -= 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );
	minZ -= 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );

	maxX += 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );
	maxY += 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );
	maxZ += 2 * ( pr->pseudoAtomDistance + pr->curvatureWeightingRadius );

	dimX = ( int ) ceil( ( maxX - minX ) / pr->gridSpacing );
	dimY = ( int ) ceil( ( maxY - minY ) / pr->gridSpacing );
	dimZ = ( int ) ceil( ( maxZ - minZ ) / pr->gridSpacing );

	cGrid->minX = minX;
	cGrid->minY = minY;
	cGrid->minZ = minZ;

	cGrid->dimX = dimX;
	cGrid->dimY = dimY;
	cGrid->dimZ = dimZ;

	cGrid->gridSpacing = pr->gridSpacing;

	int gSize = dimX * dimY * dimZ;

	cGrid->grid = ( int * ) malloc( gSize * sizeof( int ) );

	for ( int c = 0; c < gSize; c++ )
		cGrid->grid[ c ] = 0;

	cGrid->cRad = pr->curvatureWeightingRadius;

	if ( pw < 1 ) pw = 1;
	cGrid->pw = pw;

	if ( maxRatio < 2 ) maxRatio = 2;
	cGrid->maxRatio = maxRatio;

	if ( upScale < 2 ) upScale = 2;
	cGrid->upScale = upScale;

	for ( int i = 0; i < M; i++ )
		if ( type[ i ] == 'I' )
		{
			double xc = x[ i ] - minX, yc = y[ i ] - minY, zc = z[ i ] - minZ;
			double rc = r[ i ] + pr->pseudoAtomDistance;
			double rc2 = rc * rc;

			int xL, xU, yL, yU, zL, zU;

			xL = floor( ( xc - rc ) / pr->gridSpacing ) - 1, xU = ceil( ( xc + rc ) / pr->gridSpacing ) + 1;
			yL = floor( ( yc - rc ) / pr->gridSpacing ) - 1, yU = ceil( ( yc + rc ) / pr->gridSpacing ) + 1;
			zL = floor( ( zc - rc ) / pr->gridSpacing ) - 1, zU = ceil( ( zc + rc ) / pr->gridSpacing ) + 1;

			// apply atom s's influence
			for ( int zt = zL; zt <= zU; zt++ )
				for ( int yt = yL; yt <= yU; yt++ )
					for ( int xt = xL; xt <= xU; xt++ )
					{
						double d2 = ( xt * pr->gridSpacing - xc ) * ( xt * pr->gridSpacing - xc )
							+ ( yt * pr->gridSpacing - yc ) * ( yt * pr->gridSpacing - yc )
							+ ( zt * pr->gridSpacing - zc ) * ( zt * pr->gridSpacing - zc );

						if ( d2 > rc2 ) continue;

						int index = ( zt * dimY + yt ) * dimX + xt;

						cGrid->grid[ index ] = 1;
					}

			printf( "." );
			fflush( stdout );
		}

	printf( " done\n\n" );
	fflush( stdout );
}


void destroyCurvatureGrid( PARAMS_IN *pr )
{
	CURVATURE_GRID *cGrid = &( pr->cGrid );

	free( cGrid->grid );
}


double getCurvature( double x, double y, double z, PARAMS_IN *pr, int *lCount, int *hCount )
{
	CURVATURE_GRID *cGrid = &( pr->cGrid );
	double xc = x - cGrid->minX, yc = y - cGrid->minY, zc = z - cGrid->minZ;
	double rc = cGrid->cRad;
	double rc2 = rc * rc;
	double maxRatio = cGrid->maxRatio;
	double upScale = cGrid->upScale;
	int pw = cGrid->pw;

	int xL, xU, yL, yU, zL, zU;

	xL = floor( ( xc - rc ) / pr->gridSpacing ) - 1, xU = ceil( ( xc + rc ) / pr->gridSpacing ) + 1;
	yL = floor( ( yc - rc ) / pr->gridSpacing ) - 1, yU = ceil( ( yc + rc ) / pr->gridSpacing ) + 1;
	zL = floor( ( zc - rc ) / pr->gridSpacing ) - 1, zU = ceil( ( zc + rc ) / pr->gridSpacing ) + 1;

	*lCount = *hCount = 0;

	for ( int zt = zL; zt <= zU; zt++ )
		for ( int yt = yL; yt <= yU; yt++ )
			for ( int xt = xL; xt <= xU; xt++ )
			{
				double d2 = ( xt * cGrid->gridSpacing - xc ) * ( xt * cGrid->gridSpacing - xc )
					+ ( yt * cGrid->gridSpacing - yc ) * ( yt * cGrid->gridSpacing - yc )
					+ ( zt * cGrid->gridSpacing - zc ) * ( zt * cGrid->gridSpacing - zc );

				if ( d2 > rc2 ) continue;

				int index = ( zt * cGrid->dimY + yt ) * cGrid->dimX + xt;

				( *lCount ) += ( 1 - cGrid->grid[ index ] );
				( *hCount ) += cGrid->grid[ index ];
			}

	double u, v, w;

	if ( *hCount > *lCount ) v = ( ( ( *hCount ) > maxRatio * ( *lCount ) ) ? maxRatio : ( ( ( *hCount ) * 1.0 ) / ( *lCount ) ) );
	else v = ( ( ( *lCount ) > maxRatio * ( *hCount ) ) ? maxRatio : ( ( ( *lCount ) * 1.0 ) / ( *hCount ) ) );

	u = 1;

	for ( int i = 0; i < pw; i++ )
		u *= v;

	v = maxRatio;
	w = 1;

	for ( int i = 0; i < pw; i++ )
		w *= v;

	u = 1 + ( ( u - 1 ) * ( upScale - 1 ) ) / ( w - 1 );

	return u;
}


double getCurvature( double x, double y, double z, PARAMS_IN *pr )
{
	int lCount, hCount;

	return getCurvature( x, y, z, pr, &lCount, &hCount );
}


double getCurvature( int i, PARAMS_IN *pr, int *lCount, int *hCount )
{
	return getCurvature( pr->xkAOrig[ i ], pr->ykAOrig[ i ], pr->zkAOrig[ i ], pr, lCount, hCount );
}


double getCurvature( int i, PARAMS_IN *pr )
{
	int lCount, hCount;

	return getCurvature( i, pr, &lCount, &hCount );
}


bool build_fks_bak( char *type, int numCenters, double *xk, double *yk, double *zk, float *rk, float *charges, char *hbondType,
		FFTW_complex** fk, double elecWeight, FFTW_complex** fkElec, double hbondWeight, FFTW_complex** fkHbond,
		double real_magnitude, double imag_magnitude, bool staticMol, bool singleLayerLigandSkin, bool curvatureWeighting,
		double bandwidth, double gradFactor )
{
	if ( !type ) return false;

	if ( staticMol ) printf( "building fks for static molecule...\n" );
	else  printf( "building fks for moving molecule...\n" );

	fflush( stdout );

	( *fk ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( elecWeight != 0 ) ( *fkElec ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( hbondWeight != 0 ) ( *fkHbond ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );

	char *typ = ( char * ) FFTW_malloc( sizeof( char ) * numCenters );
	double *fkE2 = ( double * ) FFTW_malloc( sizeof( double ) * numCenters );

	memcpy( typ, type, numCenters * sizeof( char ) );

	if ( staticMol && curvatureWeighting )
	{
		int *idxE = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );
		double *fkE1 = ( double * ) FFTW_malloc( sizeof( double ) * numCenters );

		int num_E = 0;

		for ( int i = 0; i < numCenters; i++ )
			if ( type[ i ] == 'E' )
			{
				idxE[ num_E ] = i;
				num_E++;
			}

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE1[ ii ] = 0;

			for ( int j = 0; j < num_E; j++ )
			{
				int jj = idxE[ j ];
				double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
					+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
					+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

				if ( d2 < ( 4 * rk[ ii ] ) * ( 4 * rk[ ii ] ) ) fkE1[ ii ]++;
			}
		}

		double min_fkE2 = num_E * num_E;

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE2[ ii ] = 0;

			for ( int j = 0; j < num_E; j++ )
			{
				int jj = idxE[ j ];
				double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
					+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
					+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

				if ( d2 < ( 3 * rk[ ii ] ) * ( 3 * rk[ ii ] ) ) fkE2[ ii ] += fkE1[ jj ];
			}

			fkE2[ ii ] = sqrt( fkE2[ ii ] );

			fkE2[ ii ] = fkE1[ ii ];

			if ( fkE2[ ii ] < min_fkE2 ) min_fkE2 = fkE2[ ii ];
		}


		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE2[ ii ] /= min_fkE2;
			printf( "| %lf ", fkE2[ ii ] );
			fflush( stdout );
		}

		printf( "|\n" );
		fflush( stdout );

		FFTW_free( idxE );
		FFTW_free( fkE1 );
	}
	else
	{
		for ( int i = 0; i < numCenters; i++ )
			fkE2[ i ] = 1;
	}


	int num_I = 0;
	double imag_mag = ( gradFactor == 1.0 ) ? imag_magnitude : 0.0;

	for ( int i = 0; i < numCenters; i++ )
	{
		if ( elecWeight != 0 )
		{
			if ( staticMol && ( type[ i ] == 'E' ) ) ( *fkElec )[ i ][ 0 ] = 0;
			else ( *fkElec )[ i ][ 0 ] = ( ( FFTW_DATA_TYPE ) charges[ i ] );
		}

		if ( hbondWeight != 0 )
		{
			( *fkHbond )[ i ][ 0 ] = ( *fkHbond )[ i ][ 1 ] = 0;

			if ( staticMol )
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 0 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 1 ] = -1;
			}
			else
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 1 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 0 ] = 1;
			}
		}

		if ( type[ i ] == 'I' )
		{
			if ( staticMol || singleLayerLigandSkin )
			{
				( *fk )[ i ][ 0 ] = 0;
				( *fk )[ i ][ 1 ] = imag_mag;
			}
			else
			{
				int j;
				for ( j = 0; j < numCenters; j++ )
				{
					if ( type[ j ] != 'E' ) continue;

					double d2 = ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] )
						+ ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] )
						+ ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] );

					if ( sqrt( d2 ) <= rk[ i ] * rk[ i ] + rk[ j ] * rk[ j ] ) break;
				}

				if ( j < numCenters )
				{
					( *fk )[ i ][ 0 ] = real_magnitude;
					( *fk )[ i ][ 1 ] = 0;
				}
				else
				{
					( *fk )[ i ][ 0 ] = 0;
					( *fk )[ i ][ 1 ] = imag_mag;
				}
			}
			num_I++;
		}
		else if ( type[ i ] == 'E' )
		{
			( *fk )[ i ][ 0 ] = fkE2[ i ] * real_magnitude;
			( *fk )[ i ][ 1 ] = 0;
		}
		else
		{
			return false;
		}
	}

	FFTW_free( fkE2 );

	double bw2 = bandwidth * bandwidth;
	int num_I_prev;

	if ( gradFactor == 1.0 ) num_I = 0;

	while ( num_I > 0 )
	{
		num_I_prev = num_I;

		for ( int i = 0; i < numCenters; i++ )
		{
			if ( ( type[ i ] == 'I' ) && ( ( *fk )[ i ][ 1 ] == 0 ) )
			{
				for ( int j = 0; j < numCenters; j++ )
				{
					if ( ( j != i ) && ( ( type[ j ] == 'E' ) || ( ( *fk )[ j ][ 1 ] != 0 ) ) && ( ( *fk )[ j ][ 1 ] != ( FFTW_DATA_TYPE ) imag_magnitude ) )
					{
						if ( ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] ) + ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] ) + ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] ) <= bw2 )
						{
							( *fk )[ i ][ 1 ] = imag_magnitude;
							num_I--;
							break;
						}
					}
				}
			}
		}

		if ( num_I == num_I_prev ) bw2 *= 2;
		else
		{
			printf( "bandwidth = %lf, imag_magnitude = %lf, #atoms = %d\n", sqrt( bw2 ), imag_magnitude, num_I_prev - num_I, num_I );
			fflush( stdout );
			bw2 = bandwidth * bandwidth;
			imag_magnitude *= gradFactor;
		}
	}

	return true;
}


bool build_fks( char *type, int numCenters, double *xk, double *yk, double *zk, float *rk, float *charges, char *hbondType, float *hydrophobicity,
		FFTW_complex** fk, double elecWeight, FFTW_complex** fkElec, double hbondWeight, FFTW_complex** fkHbond,
		double hydrophobicityWeight, double hydroPhobicPhobicWeight, double hydroPhobicPhilicWeight, double hydroPhilicPhilicWeight,
		double staticMolHydroDistCutoff, FFTW_complex** fkHydrophobicity, bool twoWayHydrophobicity, FFTW_complex** fkHydrophobicityTwo,
		double simpleShapeWeight, double simpleChargeWeight, FFTW_complex** fkSimpleComplementarity,
		double real_magnitude, double imag_magnitude, bool staticMol, bool singleLayerLigandSkin,
		bool curvatureWeighting, double curvatureWeightingRadius, double bandwidth, double gradFactor, PARAMS_IN *pr )
{
	if ( !type ) return false;

	if ( staticMol ) printf( "building fks for static molecule...\n" );
	else  printf( "building fks for moving molecule...\n" );

	fflush( stdout );

	( *fk ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( elecWeight != 0 ) ( *fkElec ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( hbondWeight != 0 ) ( *fkHbond ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );

	double hydrophobicity_real_magnitude = sqrt( hydroPhobicPhobicWeight ),
	       hydrophobicity_imag_magnitude = sqrt( hydroPhilicPhilicWeight );

	if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) )
	{
		( *fkHydrophobicity ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
		if ( twoWayHydrophobicity ) ( *fkHydrophobicityTwo ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	}

	double simpleComplementarity_real_magnitude = sqrt( simpleShapeWeight ),
	       simpleComplementarity_imag_magnitude = sqrt( simpleChargeWeight );

	if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		( *fkSimpleComplementarity ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );

	double *fkE = ( double * ) FFTW_malloc( sizeof( double ) * numCenters );

	int *idxE = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );
	int *idxI = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );

	char *typ = ( char * ) FFTW_malloc( sizeof( char ) * numCenters );


	memcpy( typ, type, numCenters * sizeof( char ) );

	int num_E = 0, num_I = 0;

	for ( int i = 0; i < numCenters; i++ )
		if ( typ[ i ] == 'E' ) idxE[ num_E++ ] = i;
		else if ( typ[ i ] == 'I' ) idxI[ num_I++ ] = i;

	if ( staticMol )
	{
		double trans = 0;

		for ( int i = 0; i < numCenters; i++ )
		{
			if ( xk[ i ] < trans ) trans = xk[ i ];
			if ( yk[ i ] < trans ) trans = yk[ i ];
			if ( zk[ i ] < trans ) trans = zk[ i ];
		}

		PG *skinPG = new PG( 10.0, -trans, 5.0 );
		Point pt;

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];

			pt.x = xk[ ii ];
			pt.y = yk[ ii ];
			pt.z = zk[ ii ];

			skinPG->addPoint( &pt );
		}

		//             double distCutoff = 6.0;

		for ( int i = 0; i < num_I; i++ )
		{
			int ii = idxI[ i ];

			pt.x = xk[ ii ];
			pt.y = yk[ ii ];
			pt.z = zk[ ii ];

			double distCutoff = staticMolHydroDistCutoff + rk[ ii ];

			if ( skinPG->pointsWithinRange( &pt, distCutoff ) )
				typ[ ii ] = 'i';
		}

		delete skinPG;
	}


	if ( !staticMol && !singleLayerLigandSkin )
	{
		for ( int i = 0; i < num_I; i++ )
		{
			int ii = idxI[ i ];
			for ( int j = 0; j < num_E; j++ )
			{
				int jj = idxE[ j ];
				double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
					+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
					+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

				if ( d2 <= ( rk[ ii ] + rk[ jj ] ) * ( rk[ ii ] + rk[ jj ] ) )
				{
					typ[ ii ] = 'E';
					break;
				}
			}
		}
	}

	for ( int i = 0; i < numCenters; i++ )
		fkE[ i ] = 1;

	if ( /*staticMol &&*/ curvatureWeighting )
	{
		//            if ( staticMol )
		//              {
		//                initializeCurvatureGrid( pr, 3, 3, 9 );
		//
		//                for ( int i = 0; i < num_E; i++ )
		//                     {
		//                       int ii = idxE[ i ];
		//
		//                       int l, h;
		//
		//                       fkE[ ii ] = getCurvature( ii, pr, &l, &h );
		//
		//                       printf( "| %lf ( %d, %d ) ", fkE[ ii ], l, h );
		//                       fflush( stdout );
		//                     }
		//
		//                printf( "|\n" );
		//                fflush( stdout );
		//
		//                destroyCurvatureGrid( pr );
		//              }
		//            else
		{
			double min_fkE = num_E;

			for ( int i = 0; i < num_E; i++ )
			{
				int ii = idxE[ i ];
				fkE[ ii ] = 0;

				for ( int j = 0; j < num_E; j++ )
				{
					int jj = idxE[ j ];
					double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
						+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
						+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

					if ( d2 < curvatureWeightingRadius * curvatureWeightingRadius /*( 4 * rk[ ii ] ) * ( 4 * rk[ ii ] )*/ ) fkE[ ii ]++;
				}

				if ( fkE[ ii ] < min_fkE ) min_fkE = fkE[ ii ];
			}

			for ( int i = 0; i < num_E; i++ )
			{
				int ii = idxE[ i ];
				fkE[ ii ] /= min_fkE;
				printf( "| %lf ", fkE[ ii ] );
				fflush( stdout );
			}

			printf( "|\n" );
			fflush( stdout );
		}
	}

	std::cout<<"A "<<std::endl;

	double imag_mag = ( gradFactor == 1.0 ) ? imag_magnitude : 0.0;

	for ( int i = 0; i < numCenters; i++ )
	{
		if ( elecWeight != 0 )
		{
			if ( staticMol && ( typ[ i ] == 'E' ) ) ( *fkElec )[ i ][ 0 ] = 0;
			else ( *fkElec )[ i ][ 0 ] = ( ( FFTW_DATA_TYPE ) charges[ i ] );
		}

		if ( hbondWeight != 0 )
		{
			( *fkHbond )[ i ][ 0 ] = ( *fkHbond )[ i ][ 1 ] = 0;

			if ( staticMol )
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 0 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 1 ] = -1;
			}
			else
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 1 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 0 ] = 1;
			}
		}

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) )
		{
			( *fkHydrophobicity )[ i ][ 0 ] = ( *fkHydrophobicity )[ i ][ 1 ] = 0;

			if ( hydrophobicityWeight != 0 )
			{
				if ( staticMol )
				{
					if ( typ[ i ] == 'i' )
					{
						if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicity )[ i ][ 0 ] = hydrophobicity[ i ];
						else if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicity )[ i ][ 1 ] = -hydrophobicity[ i ];
					}
				}
				else ( *fkHydrophobicity )[ i ][ 0 ] = 1.0;

				if ( twoWayHydrophobicity )
				{
					( *fkHydrophobicityTwo )[ i ][ 0 ] = ( *fkHydrophobicityTwo )[ i ][ 1 ] = 0;

					if ( staticMol )
					{
						if ( typ[ i ] == 'i' ) ( *fkHydrophobicityTwo )[ i ][ 0 ] = 1.0;
					}
					else
					{
						if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicityTwo )[ i ][ 0 ] = hydrophobicity[ i ];
						else if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicityTwo )[ i ][ 1 ] = -hydrophobicity[ i ];
					}
				}
			}
			else
			{
				if ( ( staticMol && ( typ[ i ] == 'i' ) ) || ( !staticMol && ( type[ i ] == 'E' ) ) )
				{
					if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicity )[ i ][ 0 ] = - hydrophobicity_real_magnitude * hydrophobicity[ i ];
					else if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicity )[ i ][ 1 ] = hydrophobicity_imag_magnitude * hydrophobicity[ i ];
				}
			}
		}

		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		{
			( *fkSimpleComplementarity )[ i ][ 0 ] = ( *fkSimpleComplementarity )[ i ][ 1 ] = 0;

			if ( staticMol )
			{
				if ( ( typ[ i ] == 'I' ) || ( typ[ i ] == 'i' ) )
				{
					if ( typ[ i ] == 'i' ) ( *fkSimpleComplementarity )[ i ][ 0 ] = simpleComplementarity_real_magnitude;
					( *fkSimpleComplementarity )[ i ][ 1 ] = simpleComplementarity_imag_magnitude * charges[ i ];
				}
			}
			else
			{
				if ( typ[ i ] == 'E' ) ( *fkSimpleComplementarity )[ i ][ 0 ] = simpleComplementarity_real_magnitude;
				( *fkSimpleComplementarity )[ i ][ 1 ] = simpleComplementarity_imag_magnitude * charges[ i ];
			}
		}

		if ( ( typ[ i ] == 'I' ) || ( typ[ i ] == 'i' ) )
		{
			( *fk )[ i ][ 0 ] = 0;
			( *fk )[ i ][ 1 ] = imag_mag;
		}
		else if ( typ[ i ] == 'E' )
		{
			( *fk )[ i ][ 0 ] = fkE[ i ] * real_magnitude;
			( *fk )[ i ][ 1 ] = 0;
		}
		else
		{
			return false;
		}
	}
	std::cout<<"A1 "<<std::endl;

	FFTW_free( fkE );

	std::cout<<"B "<<std::endl;

	if ( gradFactor != 1.0 )
	{
		//            num_E = num_I = 0;
		//
		//            for ( int i = 0; i < numCenters; i++ )
		//              if ( typ[ i ] == 'E' ) idxE[ num_E++ ] = i;
		//              else if ( typ[ i ] == 'I' ) idxI[ num_I++ ] = i;
		//
		//            printf( "\nnum_E = %d, num_I = %d\n", num_E, num_I );
		//            fflush( stdout );

		num_E = num_I = 0;

		for ( int i = 0; i < numCenters; i++ )
			if ( type[ i ] == 'E' ) idxE[ num_E++ ] = i;
			else if ( type[ i ] == 'I' ) idxI[ num_I++ ] = i;

		printf( "\nnum_E = %d, num_I = %d\n\n", num_E, num_I );
		fflush( stdout );

		double bw2 = bandwidth * bandwidth;
		int j1 = 0, j2;

		while ( num_I > 0 )
		{
			j2 = num_E;

			for ( int i = 0; i < num_I; i++ )
			{
				int ii = idxI[ i ];

				for ( int j = j1; j < j2; j++ )
				{
					int jj = idxE[ j ];
					double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
						+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
						+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

					if ( d2 <= bw2 )
					{
						( *fk )[ ii ][ 1 ] = imag_magnitude;
						idxE[ num_E++ ] = ii;
						idxI[ i-- ] = idxI[ --num_I ];
						break;
					}
				}
			}

			if ( j2 == num_E ) bw2 *= 2;
			else
			{
				printf( "bandwidth = %lf, imag_magnitude = %lf, #atoms = %d\n", sqrt( bw2 ), imag_magnitude, num_E - j2 );
				fflush( stdout );
				bw2 = bandwidth * bandwidth;
				imag_magnitude *= gradFactor;
				j1 = j2;
			}
		}
	}

	FFTW_free( typ );
	FFTW_free( idxE );
	FFTW_free( idxI );

	return true;
}



bool build_fks( char *type, int numCenters, double *xk, double *yk, double *zk, float *rk, float *charges, char *hbondType, float *hydrophobicity,
		FFTW_complex** fk, FFTW_complex** fk_01, FFTW_complex** fk_10, FFTW_complex** fk_11,
		double elecWeight, FFTW_complex** fkElec, double hbondWeight, FFTW_complex** fkHbond,
		double hydrophobicityWeight, double hydroPhobicPhobicWeight, double hydroPhobicPhilicWeight, double hydroPhilicPhilicWeight,
		double staticMolHydroDistCutoff, FFTW_complex** fkHydrophobicity, bool twoWayHydrophobicity, FFTW_complex** fkHydrophobicityTwo,
		double simpleShapeWeight, double simpleChargeWeight, FFTW_complex** fkSimpleComplementarity,
		double real_magnitude, double imag_magnitude,
		bool staticMol, bool singleLayerLigandSkin, bool curvatureWeighting, double curvatureWeightingRadius,
		double bandwidth, double gradFactor, PARAMS_IN *pr )
{
	if ( !type ) return false;

	if ( staticMol ) printf( "building fks for static molecule...\n" );
	else  printf( "building fks for moving molecule...\n" );

	fflush( stdout );

	( *fk ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	( *fk_01 ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	( *fk_10 ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	( *fk_11 ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( elecWeight != 0 ) ( *fkElec ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	if ( hbondWeight != 0 ) ( *fkHbond ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );

	double hydrophobicity_real_magnitude = sqrt( hydroPhobicPhobicWeight ),
	       hydrophobicity_imag_magnitude = sqrt( hydroPhilicPhilicWeight );

	if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) )
	{
		( *fkHydrophobicity ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
		if ( twoWayHydrophobicity ) ( *fkHydrophobicityTwo ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );
	}

	double simpleComplementarity_real_magnitude = sqrt( simpleShapeWeight ),
	       simpleComplementarity_imag_magnitude = sqrt( simpleChargeWeight );

	if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		( *fkSimpleComplementarity ) = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numCenters );

	double *fkE2 = ( double * ) FFTW_malloc( sizeof( double ) * numCenters );

	int *idxE = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );
	int *idxI = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );

	char *typ = ( char * ) FFTW_malloc( sizeof( char ) * numCenters );

	memcpy( typ, type, numCenters * sizeof( char ) );

	int num_E = 0, num_I = 0;

	for ( int i = 0; i < numCenters; i++ )
		if ( typ[ i ] == 'E' ) idxE[ num_E++ ] = i;
		else if ( typ[ i ] == 'I' ) idxI[ num_I++ ] = i;

	if ( staticMol )
	{
		double trans = 0;

		for ( int i = 0; i < numCenters; i++ )
		{
			if ( xk[ i ] < trans ) trans = xk[ i ];
			if ( yk[ i ] < trans ) trans = yk[ i ];
			if ( zk[ i ] < trans ) trans = zk[ i ];
		}

		PG *skinPG = new PG( 10.0, -trans, 5.0 );
		Point pt;

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];

			pt.x = xk[ ii ];
			pt.y = yk[ ii ];
			pt.z = zk[ ii ];

			skinPG->addPoint( &pt );
		}

		//             double distCutoff = 6.0;

		for ( int i = 0; i < num_I; i++ )
		{
			int ii = idxI[ i ];

			pt.x = xk[ ii ];
			pt.y = yk[ ii ];
			pt.z = zk[ ii ];

			double distCutoff = staticMolHydroDistCutoff + rk[ ii ];

			if ( skinPG->pointsWithinRange( &pt, distCutoff ) )
				typ[ ii ] = 'i';
		}

		delete skinPG;
	}

	if ( staticMol && curvatureWeighting )
	{
		//            int *idxE = ( int * ) FFTW_malloc( sizeof( int ) * numCenters );
		double *fkE1 = ( double * ) FFTW_malloc( sizeof( double ) * numCenters );

		//            int num_E = 0;
		//
		//            for ( int i = 0; i < numCenters; i++ )
		//              if ( type[ i ] == 'E' )
		//                {
		//                  idxE[ num_E ] = i;
		//                  num_E++;
		//                }

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE1[ ii ] = 0;

			for ( int j = 0; j < num_E; j++ )
			{
				int jj = idxE[ j ];
				double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
					+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
					+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

				if ( d2 < curvatureWeightingRadius * curvatureWeightingRadius /*( 4 * rk[ ii ] ) * ( 4 * rk[ ii ] )*/ ) fkE1[ ii ]++;
			}
		}

		double min_fkE2 = num_E * num_E;

		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE2[ ii ] = 0;

			for ( int j = 0; j < num_E; j++ )
			{
				int jj = idxE[ j ];
				double d2 = ( xk[ ii ] - xk[ jj ] ) * ( xk[ ii ] - xk[ jj ] )
					+ ( yk[ ii ] - yk[ jj ] ) * ( yk[ ii ] - yk[ jj ] )
					+ ( zk[ ii ] - zk[ jj ] ) * ( zk[ ii ] - zk[ jj ] );

				if ( d2 < ( 3 * rk[ ii ] ) * ( 3 * rk[ ii ] ) ) fkE2[ ii ] += fkE1[ jj ];
			}

			fkE2[ ii ] = sqrt( fkE2[ ii ] );

			if ( fkE2[ ii ] < min_fkE2 ) min_fkE2 = fkE2[ ii ];
		}


		for ( int i = 0; i < num_E; i++ )
		{
			int ii = idxE[ i ];
			fkE2[ ii ] /= min_fkE2;
			printf( "| %lf ", fkE2[ ii ] );
			fflush( stdout );
		}

		printf( "|\n" );
		fflush( stdout );

		FFTW_free( idxE );
		FFTW_free( fkE1 );
	}
	else
	{
		for ( int i = 0; i < numCenters; i++ )
			fkE2[ i ] = 1;
	}

	int i, j;
	//        int num_I = 0;
	double imag_mag = ( gradFactor == 1.0 ) ? imag_magnitude : 0.0;
	double imag_mag2 = ( gradFactor == 1.0 ) ? 1.0 : 0.0;

	for ( i = 0; i < numCenters; i++ )
	{
		if ( elecWeight != 0 )
		{
			if ( staticMol && ( type[ i ] == 'E' ) ) ( *fkElec )[ i ][ 0 ] = 0;
			else ( *fkElec )[ i ][ 0 ] = ( ( FFTW_DATA_TYPE ) charges[ i ] );
		}

		if ( hbondWeight != 0 )
		{
			( *fkHbond )[ i ][ 0 ] = ( *fkHbond )[ i ][ 1 ] = 0;

			if ( staticMol )
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 0 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 1 ] = -1;
			}
			else
			{
				if ( hbondType[ i ] == 'D' ) ( *fkHbond )[ i ][ 1 ] = 1;
				else if ( hbondType[ i ] == 'A' ) ( *fkHbond )[ i ][ 0 ] = 1;
			}
		}

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) )
		{
			( *fkHydrophobicity )[ i ][ 0 ] = ( *fkHydrophobicity )[ i ][ 1 ] = 0;

			if ( hydrophobicityWeight != 0 )
			{
				if ( staticMol )
				{
					if ( typ[ i ] == 'i' )
					{
						if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicity )[ i ][ 0 ] = hydrophobicity[ i ];
						else if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicity )[ i ][ 1 ] = -hydrophobicity[ i ];
					}
				}
				else ( *fkHydrophobicity )[ i ][ 0 ] = 1.0;

				if ( twoWayHydrophobicity )
				{
					( *fkHydrophobicityTwo )[ i ][ 0 ] = ( *fkHydrophobicityTwo )[ i ][ 1 ] = 0;

					if ( staticMol )
					{
						if ( typ[ i ] == 'i' ) ( *fkHydrophobicityTwo )[ i ][ 0 ] = 1.0;
					}
					else
					{
						if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicityTwo )[ i ][ 0 ] = hydrophobicity[ i ];
						else if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicityTwo )[ i ][ 1 ] = -hydrophobicity[ i ];
					}
				}
			}
			else
			{
				if ( ( staticMol && ( typ[ i ] == 'i' ) ) || ( !staticMol && ( type[ i ] == 'E' ) ) )
				{
					if ( hydrophobicity[ i ] < 0 ) ( *fkHydrophobicity )[ i ][ 0 ] = - hydrophobicity_real_magnitude * hydrophobicity[ i ];
					else if ( hydrophobicity[ i ] > 0 ) ( *fkHydrophobicity )[ i ][ 1 ] = hydrophobicity_imag_magnitude * hydrophobicity[ i ];
				}
			}
		}


		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		{
			( *fkSimpleComplementarity )[ i ][ 0 ] = ( *fkSimpleComplementarity )[ i ][ 1 ] = 0;

			if ( staticMol )
			{
				if ( type[ i ] == 'I' )
				{
					( *fkSimpleComplementarity )[ i ][ 0 ] = simpleComplementarity_real_magnitude;
					( *fkSimpleComplementarity )[ i ][ 1 ] = simpleComplementarity_imag_magnitude * charges[ i ];
				}
			}
			else
			{
				if ( type[ i ] == 'E' ) ( *fkSimpleComplementarity )[ i ][ 0 ] = simpleComplementarity_real_magnitude;
				( *fkSimpleComplementarity )[ i ][ 1 ] = simpleComplementarity_imag_magnitude * charges[ i ];
			}
		}

		if ( type[ i ] == 'I' )
		{
			if ( staticMol || singleLayerLigandSkin )
			{
				( *fk )[ i ][ 0 ] = 0;
				( *fk )[ i ][ 1 ] = imag_mag;

				( *fk_01 )[ i ][ 0 ] = 0;
				( *fk_01 )[ i ][ 1 ] = imag_mag2;

				( *fk_10 )[ i ][ 0 ] = 0;
				( *fk_10 )[ i ][ 1 ] = 0;

				( *fk_11 )[ i ][ 0 ] = 0;
				( *fk_11 )[ i ][ 1 ] = imag_mag2;
			}
			else
			{
				int j;
				for ( j = 0; j < numCenters; j++ )
				{
					if ( type[ j ] != 'E' ) continue;

					double d2 = ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] )
						+ ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] )
						+ ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] );

					if ( sqrt( d2 ) <= rk[ i ] * rk[ i ] + rk[ j ] * rk[ j ] ) break;
				}

				if ( j < numCenters )
				{
					( *fk )[ i ][ 0 ] = real_magnitude;
					( *fk )[ i ][ 1 ] = 0;

					( *fk_01 )[ i ][ 0 ] = 0;
					( *fk_01 )[ i ][ 1 ] = 0;

					( *fk_10 )[ i ][ 0 ] = 1;
					( *fk_10 )[ i ][ 1 ] = 0;

					( *fk_11 )[ i ][ 0 ] = 1;
					( *fk_11 )[ i ][ 1 ] = 0;
				}
				else
				{
					( *fk )[ i ][ 0 ] = 0;
					( *fk )[ i ][ 1 ] = imag_mag;

					( *fk_01 )[ i ][ 0 ] = 0;
					( *fk_01 )[ i ][ 1 ] = imag_mag2;

					( *fk_10 )[ i ][ 0 ] = 0;
					( *fk_10 )[ i ][ 1 ] = 0;

					( *fk_11 )[ i ][ 0 ] = 0;
					( *fk_11 )[ i ][ 1 ] = imag_mag2;
				}
			}

			//                num_I++;
		}
		else if ( type[ i ] == 'E' )
		{
			( *fk )[ i ][ 0 ] = fkE2[ i ] * real_magnitude;
			( *fk )[ i ][ 1 ] = 0;

			( *fk_01 )[ i ][ 0 ] = 0;
			( *fk_01 )[ i ][ 1 ] = 0;

			( *fk_10 )[ i ][ 0 ] = 1;
			( *fk_10 )[ i ][ 1 ] = 0;

			( *fk_11 )[ i ][ 0 ] = 1;
			( *fk_11 )[ i ][ 1 ] = 0;
		}
		else
		{
			return false;
		}
	}

	FFTW_free( fkE2 );

	double bw2 = bandwidth * bandwidth;
	double imag_magnitude2 = 1.0;
	int num_I_prev;

	if ( gradFactor == 1.0 ) num_I = 0;

	while ( num_I )
	{
		num_I_prev = num_I;

		for ( i = 0; i < numCenters; i++ )
		{
			if ( ( type[ i ] == 'I' ) && ( ( *fk )[ i ][ 1 ] == 0 ) )
			{
				for ( j = 0; j < numCenters; j++ )
				{
					if ( ( j != i ) && ( ( type[ j ] == 'E' ) || ( ( *fk )[ j ][ 1 ] != 0 ) ) && ( ( *fk )[ j ][ 1 ] != ( FFTW_DATA_TYPE ) imag_magnitude ) )
					{
						if ( ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] ) + ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] ) + ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] ) <= bw2 )
						{
							( *fk )[ i ][ 1 ] = imag_magnitude;
							( *fk_01 )[ i ][ 1 ] = imag_magnitude2;
							( *fk_11 )[ i ][ 1 ] = imag_magnitude2;

							num_I--;
							break;
						}
					}
				}
			}
		}

		if ( num_I == num_I_prev ) bw2 *= 2;
		else
		{
			printf( "bandwidth = %lf, imag_magnitude = %lf, #atoms = %d\n", sqrt( bw2 ), imag_magnitude, num_I_prev - num_I );
			fflush( stdout );
			bw2 = bandwidth * bandwidth;
			imag_magnitude *= gradFactor;
			imag_magnitude2 *= gradFactor;
		}
	}

	return true;
}




/*******************************************************************/
/*                                                                 */
/* Rotates a given point about the origin by the given angles      */
/*                                                                 */
/*******************************************************************/
void rotatePointAboutOrigin( double* x, double* y, double* z, double theta1, double theta2, double theta3 )
{
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;

	// multiply by theta1
	x1 =  cos( theta1 ) * (*x)  + sin( theta1 ) * (*y)  + 0         * (*z);
	y1 = -sin( theta1 ) * (*x)  + cos( theta1 ) * (*y)  + 0         * (*z);
	z1 = 0              * (*x)  + 0             * (*y)  + 1         * (*z);

	// multiply by theta2
	x2 =  cos( theta2 ) * x1 + 0              * y1 + -sin( theta2 ) * z1;
	y2 = 0              * x1 + 1              * y1 + 0              * z1;
	z2 =  sin( theta2 ) * x1 + 0              * y1 +  cos( theta2 ) * z1;

	// multiply by theta3
	x3 = 1              * x2 + 0              * y2 + 0              * z2;
	y3 = 0              * x2 +  cos( theta3 ) * y2 +  sin( theta3 ) * z2;
	z3 = 0              * x2 + -sin( theta3 ) * y2 +  cos( theta3 ) * z2;

	*x = x3; *y = y3; *z = z3;
}
// here we are given the rotation  matrix itself.
void rotatePointAboutOrigin( double* x, double* y, double* z,
		float r1, float r2, float r3,
		float r4, float r5, float r6,
		float r7, float r8, float r9 )
{
	double x1, y1, z1;

	x1 = r1*(*x) + r2*(*y) + r3*(*z);
	y1 = r4*(*x) + r5*(*y) + r6*(*z);
	z1 = r7*(*x) + r8*(*y) + r9*(*z);

	(*x) = x1;
	(*y) = y1;
	(*z) = z1;
}


bool closeRotations( float *rot, int i, int j, double cosTheta )
{
	// Julie C. Mitchell ( Rotation Samples )
	// Lattman Formula?
	double d = 0;

	for ( int k = 0; k < 9; k++ )
		d += rot[ i * 9 + k ] * rot[ j * 9 + k ];

	d = ( d - 1.0 ) / 2.0 ;

	if ( d >  1.0 ) d =  1.0;
	if ( d < -1.0 ) d = -1.0;

	if ( d >= cosTheta ) return true;
	else return false;
}


//bool closeRotations( float *rot, int i, int j, double cosTheta )
//{
//   double dotP = 0;
//
//   for ( int k = 0; k < 3; k++ )
//     {
//      double di = 0, dj = 0;
//
//      for ( int l = 0; l < 3; l++ )
//        {
//         di += rot[ i * 9 + k * 3 + l ];
//         dj += rot[ j * 9 + k * 3 + l ];
//        }
//
//      dotP += di * dj;
//     }
//
//   if ( dotP >= 3 * cosTheta ) return true;
//   else return false;
//}


bool closePeaksExist( int x, int y, int z, int r, float *rotations, int *grid, int *peakList, int numFreq, double gridSpacing, double clusterTransRad, double cosTheta )
{
	double dist = clusterTransRad / gridSpacing;
	double dist2 = dist * dist;
	int dCells = ( int ) ceil( dist );

	if ( x > ( numFreq >> 1 ) ) x -= numFreq;
	if ( y > ( numFreq >> 1 ) ) y -= numFreq;
	if ( z > ( numFreq >> 1 ) ) z -= numFreq;

	int lx = x - dCells, hx = x + dCells;
	int ly = y - dCells, hy = y + dCells;
	int lz = z - dCells, hz = z + dCells;

	if ( lx <= - ( numFreq >> 1 ) ) lx = - ( numFreq >> 1 ) + 1;
	if ( ly <= - ( numFreq >> 1 ) ) ly = - ( numFreq >> 1 ) + 1;
	if ( lz <= - ( numFreq >> 1 ) ) lz = - ( numFreq >> 1 ) + 1;

	if ( hx > ( numFreq >> 1 ) ) hx = ( numFreq >> 1 );
	if ( hy > ( numFreq >> 1 ) ) hy = ( numFreq >> 1 );
	if ( hz > ( numFreq >> 1 ) ) hz = ( numFreq >> 1 );

	for ( int xx = lx; xx <= hx; xx++ )
		for ( int yy = ly; yy <= hy; yy++ )
			for ( int zz = lz; zz <= hz; zz++ )
			{
				double d = ( x - xx ) * ( x - xx )
					+ ( y - yy ) * ( y - yy )
					+ ( z - zz ) * ( z - zz );

				if ( d > dist2 ) continue;

				int nx = xx, ny = yy, nz = zz;

				if ( nx < 0 ) nx += numFreq;
				if ( ny < 0 ) ny += numFreq;
				if ( nz < 0 ) nz += numFreq;

				int c = ( nx * numFreq + ny ) * numFreq + nz;

				int l = grid[ c ];

				while ( l >= 0 )
				{
					int rr = peakList[ l ];

					if ( closeRotations( rotations, r, rr, cosTheta ) ) return true;

					l = peakList[ l + 1 ];
				}
			}

	return false;
}



// return the min max of the array
void getMinMax( double* x, double* minVal, double* maxVal, int n)
{
	(*minVal) = (*maxVal) = x[0];

	int i;
	for( i=1; i<n; i++ )
	{
		if( (*minVal) > x[i] ) (*minVal) = x[i];
		if( (*maxVal) < x[i] ) (*maxVal) = x[i];
	}
}


bool findCenter( double* x, double* y, double* z, int n, double *xC, double *yC, double *zC )
{
	if ( !x || !y || !z || ( n < 1 ) ) return false;

	double minXVal, maxXVal;
	double minYVal, maxYVal;
	double minZVal, maxZVal;

	getMinMax( x, &minXVal, &maxXVal, n );
	getMinMax( y, &minYVal, &maxYVal, n );
	getMinMax( z, &minZVal, &maxZVal, n );

	*xC = ( minXVal + maxXVal ) / 2.0;
	*yC = ( minYVal + maxYVal ) / 2.0;
	*zC = ( minZVal + maxZVal ) / 2.0;

	//	*xC = - minXVal - ( maxXVal - minXVal ) / 2.0;
	//	*yC = - minYVal - ( maxYVal - minYVal ) / 2.0;
	//	*zC = - minZVal - ( maxZVal - minZVal ) / 2.0;

	return true;
}


bool center( double* x, double* y, double* z, int n, double *xTrans, double *yTrans, double *zTrans )
{
	if ( !x || !y || !z || ( n < 1 ) ) return false;

	findCenter( x, y, z, n, xTrans, yTrans, zTrans );

	*xTrans = - ( *xTrans );
	*yTrans = - ( *yTrans );
	*zTrans = - ( *zTrans );

	for ( int i = 0; i < n; i++ )
	{
		x[ i ] = x[ i ] + ( *xTrans );
		y[ i ] = y[ i ] + ( *yTrans );
		z[ i ] = z[ i ] + ( *zTrans );
	}

	return true;
}




bool rotateAboutOrigin( double* x, double* y, double* z, int n, Matrix &rotMat )
{
	if ( !x || !y || !z || ( n < 1 ) ) return false;

	for ( int i = 0; i < n; i++ )
	{
		Vector oldPos( x[ i ], y[ i ], z[ i ], 1 );
		Vector newPos = rotMat * oldPos;

		x[ i ] = newPos[ 0 ];
		y[ i ] = newPos[ 1 ];
		z[ i ] = newPos[ 2 ];
	}

	return true;
}


bool invert( double* x, double* y, double* z, int n )
{
	if( !x || !y || !z || (n<1) ) return false;

	int i;
	for( i=0; i<n; i++ )
	{
		x[i] = -x[i];
		y[i] = -y[i];
		z[i] = -z[i];
	}
	return true;
}

bool getLargestSize( double* x, double* y, double* z, int n, double* maxLength )
{
	double minVal, maxVal;

	getMinMax( x, &minVal, &maxVal, n);
	double xLen = maxVal - minVal;

	getMinMax( y, &minVal, &maxVal, n);
	double yLen = maxVal - minVal;

	getMinMax( z, &minVal, &maxVal, n);
	double zLen = maxVal - minVal;

	if( xLen >= yLen && xLen >= zLen )
		(*maxLength) = xLen;
	else if( yLen >= xLen && yLen >= zLen )
		(*maxLength) = yLen;
	else
		(*maxLength) = zLen;
	return true;
}


/* Given n atom centers with the center of the i-th atom given by ( x[ i ], y[ i ], z[ i ] ),
   returns the distance between the two atoms that lie the furthest from each other. */

bool getLargestPairwiseDistance( double* x, double* y, double* z, int n, double* maxLength )
{
	*maxLength = 0;

	for ( int i = 0; i < n; i++ )
		for ( int j = 0; j < n; j++ )
		{
			double d = ( x[ i ] - x[ j ] ) * ( x[ i ] - x[ j ] )
				+ ( y[ i ] - y[ j ] ) * ( y[ i ] - y[ j ] )
				+ ( z[ i ] - z[ j ] ) * ( z[ i ] - z[ j ] );

			if ( d > *maxLength ) *maxLength = d;
		}

	*maxLength = sqrt( *maxLength );

	return true;
}


bool getMaxDistanceFromCenter( double* x, double* y, double* z, int n, double* maxDist )
{
	double xC, yC, zC;

	if ( !findCenter( x, y, z, n, &xC, &yC, &zC ) ) return false;

	*maxDist = 0;

	for ( int i = 0; i < n; i++ )
	{
		double d = ( x[ i ] - xC ) * ( x[ i ] - xC )
			+ ( y[ i ] - yC ) * ( y[ i ] - yC )
			+ ( z[ i ] - zC ) * ( z[ i ] - zC );

		if ( d > *maxDist ) *maxDist = d;
	}

	*maxDist = sqrt( *maxDist );

	return true;
}



/* Given two arrays of values (i.e., radii) returns the largest value (i.e., radius) among them. */
/* Array r1 contains n1 floating point values (i.e., doubles), and array r2 contains n2 of them. */

bool getMaxRadius( float* r1, int n1, float* r2, int n2, float *maxRadius )
{
	if ( !r1 || !r2 || ( n1 < 1 ) || ( n2 < 1 ) || !maxRadius ) return false;

	*maxRadius = r1[ 0 ];

	for ( int i = 1; i < n1; i++ )
		if ( *maxRadius < r1[ i ] ) *maxRadius = r1[ i ];

	for ( int i = 0; i < n2; i++ )
		if ( *maxRadius < r2[ i ] ) *maxRadius = r2[ i ];

	return true;
}


/* The inputs are two molecules, i.e., sets of atoms. The center of the i-th atom of molecule 1 is given by ( x1[ i ], y1[ i ], z1[ i ] )
   and its radius is r1[ i ]. The center and radius of the i-th atom of molecule 2 are given by ( x2[ i ], y2[ i ], z2[ i ] ) and r2[ i ],
   respectively. This function returns the dimension of the smallest axis-aligned cube that can contain the first molecule in the given
   orientation and the second molecule in any orientation side by side assuming that each atom has radius r_max, where r_max is the
   radius of the largest atom in the two molecules. */

bool getLargestEdge( double* x1, double* y1, double* z1, float *r1, int n1, double* x2, double* y2, double* z2, float *r2, int n2, double* largestEdge )
{
	if ( !x1 || !y1 || !z1 || !r1 || !x2 || !y2 || !z2 || !r2 || ( n1 < 1 ) || ( n2 < 1 ) || !largestEdge ) return false;

	double l1 = 0;
	// compute the dimension l1 of the smallest cube that can contain the atom centers of the first molecule in the given orientation
	if ( !getLargestSize( x1, y1, z1, n1, &l1 ) ) return false;

	double l2 = 0;
	// compute the dimension l2 of the smallest cube that can contain the atom centers of the second molecule in any orientation
	//	if ( !getLargestSize( x2, y2, z2, n2, &l2 ) ) return false;
	if ( !getLargestPairwiseDistance( x2, y2, z2, n2, &l2 ) ) return false;

	float r_max = 0;
	// compute the radius of the largest atom in the two molecules
	if ( !getMaxRadius( r1, n1, r2, n2, &r_max ) ) return false;

	*largestEdge = l1 + 2 * r_max + 2 * l2;

	return true;
}


// Computes the value of the interpFuncExtent which is basically the radius of the largest atom in terms of grid points.
// This value is used by the smoothing function.

bool computeInterpFuncExtent( double* x1, double* y1, double* z1, float *r1, int n1, double* x2, double* y2, double* z2, float *r2, int n2, int numFreq, int* interpFuncExtent )
{
	double largestEdge = 0;
	// compute the dimension of the smallest axis-aligned cube that can contain either molecule in any orientation
	// assuming that each atom has radius r_max, where r_max is the radius of the largest atom in the two molecules
	if ( !getLargestEdge( x1, y1, z1, r1, n1, x2, y2, z2, r2, n2, &largestEdge ) ) return false;

	float r_max = 0;
	// compute the radius of the largest atom in the two molecules
	if ( !getMaxRadius( r1, n1, r2, n2, &r_max ) ) return false;

	*interpFuncExtent = ceil( r_max * ( numFreq / largestEdge ) );
	printf("interpFuncExtent %d\n", *interpFuncExtent);

	return true;
}


bool computeGridParametersFromVolume( PARAMS_IN *pr, int *interpFuncExtent, double *scale )
{
	int xDim, yDim, zDim;
	int numFreq;
	double xCenter, yCenter, zCenter;
	double xc, yc, zc, sc, EPS=0.0001;

	// read header of static skin
	if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xCenter, &yCenter, &zCenter, scale, pr->staticMoleculeSCReRaw ) ) return false;
	if ( ( xDim != yDim ) || ( yDim != zDim ) || ( zDim != xDim ) )
	{
		printf( "\nError: The grid in %s is not cubic!\n", pr->staticMoleculeSCReRaw );
		return false;
	}
	numFreq = xDim;
	if ( numFreq & 1 )
	{
		printf( "\nError: The cubic grid in %s does not even dimensions!\n", pr->staticMoleculeSCReRaw );
		return false;
	}

	// read header of static core
	if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, pr->staticMoleculeSCImRaw ) ) return false;
	if ( ( xDim != numFreq ) || ( yDim != numFreq ) || ( zDim != numFreq ) )
	{
		printf( "\nError1: The dimensions of the grids in %s and %s are different!\n", pr->staticMoleculeSCReRaw, pr->staticMoleculeSCImRaw );
		return false;
	}
	if ( ( fabs( xCenter -xc ) > EPS) || ( fabs( yCenter - yc) > EPS ) || ( fabs( zCenter - zc ) > EPS ) || ( fabs( *scale - sc ) > EPS ) )
	{
		printf( "\nError2: The grids in %s and %s have different centers and/or scales!\n", pr->staticMoleculeSCReRaw, pr->staticMoleculeSCImRaw );
		return false;
	}

	if ( pr->elecScale != 0 )
	{
		// read header of static elect grid
		if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, pr->staticMoleculeElecReRaw ) ) return false;
		if ( ( xDim != numFreq ) || ( yDim != numFreq ) || ( zDim != numFreq ) )
		{
			printf( "\nError3: The dimensions of the grids in %s and %s are different!\n", pr->staticMoleculeSCReRaw, pr->staticMoleculeElecReRaw );
			return false;
		}
		if ( ( fabs( xCenter -xc ) > EPS) || ( fabs( yCenter - yc) > EPS) || ( fabs( zCenter - zc) > EPS ) || ( fabs( *scale - sc ) > EPS ) )
		{
			printf( "\nError4: The grids in %s and %s have different centers and/or scales!\n", pr->staticMoleculeSCReRaw, pr->staticMoleculeElecReRaw );
			return false;
		}
	}


	// read header of moving skin
	if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xCenter, &yCenter, &zCenter, &sc, pr->movingMoleculeSCReRaw ) ) return false;
	if ( ( xDim != numFreq ) || ( yDim != numFreq ) || ( zDim != numFreq ) )
	{
		printf( "\nError5: The dimensions of the grids in %s and %s are different!\n", pr->staticMoleculeSCReRaw, pr->movingMoleculeSCReRaw );
		return false;
	}
	if ( fabs( *scale - sc ) > EPS )
	{
		printf( "\nError: The grids in %s and %s have different scales %f %f %f!\n", pr->staticMoleculeSCReRaw, pr->movingMoleculeSCReRaw, *scale, sc, *scale-sc );
		return false;
	}

	// read header of moving core
	if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, pr->movingMoleculeSCImRaw ) ) return false;
	if ( ( xDim != numFreq ) || ( yDim != numFreq ) || ( zDim != numFreq ) )
	{
		printf( "\nError6: The dimensions of the grids in %s and %s are different!\n", pr->movingMoleculeSCReRaw, pr->movingMoleculeSCImRaw );
		return false;
	}
	if ( ( fabs( xCenter -xc ) > EPS) || ( fabs( yCenter - yc) > EPS) || ( fabs( zCenter - zc) > EPS ) || ( fabs( *scale - sc ) > EPS ) )
	{
		printf( "\nError7: The grids in %s and %s have different centers and/or scales %d %d %d %d, %f %f %f!\n", pr->movingMoleculeSCReRaw, pr->movingMoleculeSCImRaw,
				( xCenter != xc ),( yCenter != yc ), ( zCenter != zc ), ( fabs( *scale - sc ) > EPS ), yCenter, yc, yCenter-yc);
		return false;
	}

	if ( pr->elecScale != 0 )
	{
		if ( !readRAWIVHeader( &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, pr->movingMoleculeElecReRaw ) ) return false;
		if ( ( xDim != numFreq ) || ( yDim != numFreq ) || ( zDim != numFreq ) )
		{
			printf( "\nError8: The dimensions of the grids in %s and %s are different!\n", pr->movingMoleculeSCReRaw, pr->movingMoleculeElecReRaw );
			return false;
		}
		if ( ( fabs( xCenter -xc ) > EPS) || ( fabs( yCenter - yc) > EPS) || ( fabs( zCenter - zc) > EPS ) || ( fabs( *scale - sc ) > EPS ) )
		{
			printf( "\nError9: The grids in %s and %s have different centers and/or scales!\n", pr->movingMoleculeSCReRaw, pr->movingMoleculeElecReRaw );
			return false;
		}
	}

	pr->numFreq = numFreq;
	double gridLength1D = 1.0 / ( *scale );
	pr->gridSpacing = gridLength1D / ( pr->numFreq - 1 );

	*interpFuncExtent = ceil( ( pr->numFreq - 1 ) * ( pr->interpFuncExtentInAngstroms / gridLength1D ) );

	if ( !pr->gridSizeSpecified || ( pr->gridSize < pr->numFreq ) ) pr->gridSize = 4 * pr->numFreq;

	return true;
}


/* The inputs are two molecules, i.e., sets of atoms. The center of the i-th atom of molecule 1 is given by ( x1[ i ], y1[ i ], z1[ i ] )
   and its radius is r1[ i ]. The center and radius of the i-th atom of molecule 2 are given by ( x2[ i ], y2[ i ], z2[ i ] ) and r2[ i ],
   respectively. This function first computes the dimension of the smallest axis-aligned cube that can contain the first molecule in the given
   orientation and the second molecule in any orientation side by side assuming that each atom has radius r_max, where r_max is the
   radius of the largest atom in the two molecules. Then it determines the size of the spatial FFTW grid and grid spacing based
   on what the user has requested and the grid sizes for which FFTW is efficient. It also computes the interpFuncExtent parameter
   and scale which is the inverse of the length of the spatial grid. */

bool computeGridParameters( PARAMS_IN *pr, int *interpFuncExtent, double *scale )
{
	double r_maxFactor = 1.5, BFactor = 2.0;

	if ( pr->dockVolume )
	{
		bool retVal = computeGridParametersFromVolume( pr, interpFuncExtent, scale );
		if ( !pr->smoothSkin || !retVal || ( *interpFuncExtent > 0 ) ) return retVal;
	}

	double l1 = 0;
	// compute the dimension l1 of the smallest cube that can contain the atom centers of the first molecule in the given orientation
	if ( !getLargestSize( pr->xkAOrig, pr->ykAOrig, pr->zkAOrig, pr->numCentersA, &l1 ) ) return false;

	double l2 = 0;
	// compute the dimension l2 of the smallest cube that can contain the atom centers of the second molecule in any orientation
	// if ( !getLargestPairwiseDistance( pr->xkBOrig, pr->ykBOrig, pr->zkBOrig, pr->numCentersB, &l2 ) ) return false;

	// compute the distance of the atom farthest from the center of the second molecule
	if ( !getMaxDistanceFromCenter( pr->xkBOrig, pr->ykBOrig, pr->zkBOrig, pr->numCentersB, &l2 ) ) return false;


	float r_max = 0;
	// compute the radius of the largest atom in the two molecules
	if ( !getMaxRadius( pr->radiiA, pr->numCentersA, pr->radiiB, pr->numCentersB, &r_max ) ) return false;

	if ( pr->dockVolume )
	{
		double gridLength1DVol = 1.0 / ( *scale );
		double interpFuncExtentInAngstroms = ( gridLength1DVol - ( l1 + 2 * r_max + BFactor * l2 ) ) / 2.0;
		if ( interpFuncExtentInAngstroms <= 0 )
		{
			printf( "\nError: interpFuncExtentInAngstroms = %lf!\n", interpFuncExtentInAngstroms );
			return false;
		}
		*interpFuncExtent = ceil( ( pr->numFreq - 1 ) * ( interpFuncExtentInAngstroms / gridLength1DVol ) );
	}
	else
	{
		double gridLength1D = l1 + 2 * r_max + BFactor * l2;
		double interpFuncExtentInAngstroms = r_maxFactor * r_max;
		gridLength1D += 2 * interpFuncExtentInAngstroms;   // for interpFuncExtent

		if ( ( pr->gridSpacingSpecified ) || ( !pr->gridSpacingSpecified && !pr->numFreqSpecified ) )
		{
			pr->numFreq = ceil( gridLength1D / pr->gridSpacing ) + 1;
			pr->numFreq += ( pr->numFreq & 1 );
			for ( int i = 0; i < pr->numEfficientGridSizes; i++ )
				if ( pr->efficientGridSizes[ i ] >= pr->numFreq )
				{
					pr->numFreq = pr->efficientGridSizes[ i ];
					break;
				}

			if ( pr->enforceExactGridSpacing ) gridLength1D = ( pr->numFreq - 1 ) * pr->gridSpacing;
		}

		pr->gridSpacing = gridLength1D / ( pr->numFreq - 1 );

		*interpFuncExtent = ceil( ( pr->numFreq - 1 ) * ( interpFuncExtentInAngstroms / gridLength1D ) );
		*scale = 1 / gridLength1D;

		if ( !pr->gridSizeSpecified || ( pr->gridSize < pr->numFreq ) ) pr->gridSize = 4 * pr->numFreq;
	}

	if ( pr->peaksPerRotation <= 0 ) pr->peaksPerRotation = pr->numFreq * pr->numFreq * pr->numFreq;
	if ( pr->filterDepth <= 0 ) pr->filterDepth = pr->numFreq * pr->numFreq * pr->numFreq;

	return true;
}




bool transformAndNormalize( double *xkOrig, double* ykOrig, double* zkOrig, float *rkOrig,
		double* xk, double* yk, double* zk, float *rk,
		int numCenters,
		float r1, float r2, float r3, float r4, float r5, float r6, float r7, float r8, float r9,
		int gridSize, int numFreq, int interpFuncExtent, double largestSize, int protein,
		//						   double* xTrans, double* yTrans, double* zTrans,
		double* scaleVal )
{
	if( !xkOrig || !ykOrig || !zkOrig || !xk || !yk || !zk || (numCenters<1) || (gridSize<1) ) return false;

	// copy
	{
		int i;
		for( i=0; i<numCenters; i++ )
		{
			xk[ i ] = xkOrig[ i ];
			yk[ i ] = ykOrig[ i ];
			zk[ i ] = zkOrig[ i ];
			rk[ i ] = rkOrig[ i ];
		}
	}


	// center
	//	if( !center( xk, yk, zk, numCenters, xTrans, yTrans, zTrans ) ) return false;

	if( protein == 2 )
	{
		// if we are looking at the second molecule, we flip it.
		// This is because the convolution flips one molecule.
		if( !invert( xk, yk, zk, numCenters ) ) return false;
		//		printf("Flipped about the origin\n");

		// rotate the second protein as well
		int i;
		for( i=0; i<numCenters; i++ )
		{
			rotatePointAboutOrigin( &(xk[i]), &(yk[i]), &(zk[i]), r1, r2, r3, r4, r5, r6, r7, r8, r9 );
		}
	}
	// normalize
	{
		double scale = 1.0/largestSize;
		/*		// [-0.5 .. 0.5)
		//(for all rotations, for both molecules )
		//		scale *= 0.5;                                                       // [-0.25 .. 0.25) */
		scale *= numFreq / ( ( double ) ( numFreq + 2 * interpFuncExtent ) );     // add space to perform gridding

		*scaleVal = scale;

		//		printf("Scaled by : [%15.10lf]\n", scale );

		int i;
		for( i=0; i<numCenters; i++ )
		{
			xk[ i ] *= scale;
			yk[ i ] *= scale;
			zk[ i ] *= scale;
			rk[ i ] *= scale;
		}
	}

	return true;
}



bool transformAndNormalize( double *xkOrig, double* ykOrig, double* zkOrig, float *rkOrig,
		double* xk, double* yk, double* zk, float *rk,
		int numCenters,
		float r1, float r2, float r3, float r4, float r5, float r6, float r7, float r8, float r9,
		int protein, double scale )
{
	if( !xkOrig || !ykOrig || !zkOrig || !xk || !yk || !zk || (numCenters<1) ) return false;

	// copy
	{
		for( int i=0; i<numCenters; i++ )
		{
			xk[ i ] = xkOrig[ i ];
			yk[ i ] = ykOrig[ i ];
			zk[ i ] = zkOrig[ i ];
			rk[ i ] = rkOrig[ i ];
		}
	}


	if( protein == 2 )
	{
		// if we are looking at the second molecule, we flip it.
		// This is because the convolution flips one molecule.
		if( !invert( xk, yk, zk, numCenters ) ) return false;

		// rotate the second protein as well
		for( int i=0; i<numCenters; i++ )
		{
			rotatePointAboutOrigin( &(xk[i]), &(yk[i]), &(zk[i]), r1, r2, r3, r4, r5, r6, r7, r8, r9 );
		}
	}
	// normalize
	{
		for( int i=0; i<numCenters; i++ )
		{
			xk[ i ] *= scale;
			yk[ i ] *= scale;
			zk[ i ] *= scale;
			rk[ i ] *= scale;
		}
	}

	return true;
}


float getRMSD(PARAMS_IN *pr, Matrix transformation )
{
	float *oldposx = pr->xRef;
	float *oldposy = pr->yRef;
	float *oldposz = pr->zRef;
	double rmsd=0.0;
	int *atIndex = pr->atNums;
	int i;

	if  (pr->nbRMSDAtoms == 0) return -1.;

	for( i=0; i < pr->nbRMSDAtoms; i++ )
	{
		Vector oldPos( *oldposx, *oldposy, *oldposz, 1.);
		Vector newPos = transformation * oldPos;
		rmsd += ( oldPos[ 0 ] - newPos[ 0 ] ) * ( oldPos[ 0 ] - newPos[ 0 ] )
			+	( oldPos[ 1 ] - newPos[ 1 ] ) * ( oldPos[ 1 ] - newPos[ 1 ] )
			+	( oldPos[ 2 ] - newPos[ 2 ] ) * ( oldPos[ 2 ] - newPos[ 2 ] );
		oldposx++;
		oldposy++;
		oldposz++;
	}

	return floor( 10 * sqrt( rmsd / ( ( double ) pr->nbRMSDAtoms ) ) ) / 10.0;
}


void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		int n, float scale, float *translate_A, float *translate_B, float *rotations,
		double functionScaleFactor, Matrix randRot,
		double &rx, double &ry, double &rz,
		Matrix &trans )
{
	c = 0;

	trans = Matrix::translation( translate_B[ 3 * c + 0 ], translate_B[ 3 * c + 1 ], translate_B[ 3 * c + 2 ] );

	trans = trans.preMultiplication( Matrix( rotations[ r * 9 + 0 ], rotations[ r * 9 + 1 ], rotations[ r * 9 + 2 ], 0,
				rotations[ r * 9 + 3 ], rotations[ r * 9 + 4 ], rotations[ r * 9 + 5 ], 0,
				rotations[ r * 9 + 6 ], rotations[ r * 9 + 7 ], rotations[ r * 9 + 8 ], 0,
				0,                      0,                      0, 1 ) );

	trans = trans.preMultiplication( randRot );

	//       transf = trans.preMultiplication( Matrix( fineRotations[ f * 9 + 0 ], fineRotations[ f * 9 + 1 ], fineRotations[ f * 9 + 2 ], 0,
	// 		 			           fineRotations[ f * 9 + 3 ], fineRotations[ f * 9 + 4 ], fineRotations[ f * 9 + 5 ], 0,
	// 						   fineRotations[ f * 9 + 6 ], fineRotations[ f * 9 + 7 ], fineRotations[ f * 9 + 8 ], 0,
	// 						                            0,                      0,                      0, 1 ) );

	trans = trans.preMultiplication( Matrix::translation( -translate_B[ 3 * c + 0 ], -translate_B[ 3 * c + 1 ], -translate_B[ 3 * c + 2 ] ) );

	// now get the translation in real space
	// get the real pos
	double realx = translate_A[ 0 ] - translate_B[ 3 * c + 0 ];
	double realy = translate_A[ 1 ] - translate_B[ 3 * c + 1 ];
	double realz = translate_A[ 2 ] - translate_B[ 3 * c + 2 ];

	double x = ix, y = iy, z = iz;

	// normalize it if needed
	if ( x > ( n - 1 ) / 2.0 ) x -= ( n - 1 );
	if ( y > ( n - 1 ) / 2.0 ) y -= ( n - 1 );
	if ( z > ( n - 1 ) / 2.0 ) z -= ( n - 1 );

	//    int nn = ( n >> 1 ), odd = ( n & 1 );
	//
	//    if ( x >= nn + odd ) x -= n;
	//    if ( y >= nn + odd ) y -= n;
	//    if ( z >= nn + odd ) z -= n;

	x /= ( ( n - 1 ) * scale );
	y /= ( ( n - 1 ) * scale );
	z /= ( ( n - 1 ) * scale );

	// now get real
	rx = x - realx;
	ry = y - realy;
	rz = z - realz;

	// find final transformation by adding the translation
	trans = trans.preMultiplication( Matrix::translation( rx, ry, rz ) );
}


void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		PARAMS *pr,
		double &rx, double &ry, double &rz,
		Matrix &trans )
{
	retrieveTransformation( ix, iy, iz, c, r, f,
			pr->numFreq, pr->scaleB, pr->translate_A, pr->translate_B, pr->rotations,
			pr->functionScaleFactor, pr->randRot,
			rx, ry, rz,
			trans );
}


void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		FILTER_PARAMS *pr,
		double &rx, double &ry, double &rz,
		Matrix &trans )
{
	retrieveTransformation( ix, iy, iz, c, r, f,
			pr->numFreq, pr->scaleB, pr->translate_A, pr->translate_B, pr->rotations,
			pr->functionScaleFactor, pr->randRot,
			rx, ry, rz,
			trans );
}


void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		PARAMS *pr,
		double &rx, double &ry, double &rz,
		Matrix &trans, double *dtrans )
{
	retrieveTransformation( ix, iy, iz, c, r, f,
			pr,
			rx, ry, rz,
			trans );

	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			dtrans[ i * 4 + j ] = trans.get( i, j );
}

void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		PARAMS *pr,
		Matrix &trans )
{
	double rx, ry, rz;

	retrieveTransformation( ix, iy, iz, c, r, f,
			pr,
			rx, ry, rz,
			trans );
}


void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
		PARAMS *pr,
		Matrix &trans, double *dtrans )
{
	double rx, ry, rz;

	retrieveTransformation( ix, iy, iz, c, r, f,
			pr,
			rx, ry, rz,
			trans, dtrans );
}


void retrieveTransformation( ValuePosition3D &sol,
		FILTER_PARAMS *pr,
		Matrix &trans, double *dtrans )
{
	double rx, ry, rz;

	retrieveTransformation( ( int ) sol.m_Translation[ 0 ], ( int ) sol.m_Translation[ 1 ], ( int ) sol.m_Translation[ 2 ],
			( int ) sol.m_ConformationIndex, ( int ) sol.m_RotationIndex, ( int ) sol.m_FineRotationIndex,
			pr,
			rx, ry, rz,
			trans );

	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			dtrans[ i * 4 + j ] = trans.get( i, j );
}


void retrieveTransformation( ValuePosition3D &sol,
		PARAMS *pr,
		Matrix &trans, double *dtrans )
{
	double rx, ry, rz;

	retrieveTransformation( ( int ) sol.m_Translation[ 0 ], ( int ) sol.m_Translation[ 1 ], ( int ) sol.m_Translation[ 2 ],
			( int ) sol.m_ConformationIndex, ( int ) sol.m_RotationIndex, ( int ) sol.m_FineRotationIndex,
			pr,
			rx, ry, rz,
			trans );

	for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			dtrans[ i * 4 + j ] = trans.get( i, j );
}

/*
   void f_printf( FILE *fp, char *format, ... )
   {
   va_list args, args2;

   va_start( args, format );
   va_copy( args2, args );

   vprintf( format, args );
   vfprintf( fp, format, args2 );

   va_end( args );
   }
 */

void printIntermediateStats( FILE *fp,
		int rmsdToReport,
		int rmsdGood,
		PARAMS *pr )
{
	bool reDocking = (pr->pri->nbRMSDAtoms > 0) ? true:false;
	bool breakDownScores = (pr->pri->breakDownScores==1) ? true:false;
	double skinSkinWeight = pr->skinSkinWeight;
	double skinCoreWeight = pr->skinCoreWeight;
	double coreCoreWeight = pr->coreCoreWeight;

	double unrealWeight = 0;

	if ( !breakDownScores )
		unrealWeight = ( skinSkinWeight * coreCoreWeight > 0 ) ? ( skinCoreWeight / sqrt( skinSkinWeight * coreCoreWeight ) ) : 0.0;

	float scale_B = pr->scaleB;
	double functionScaleFactor = pr->functionScaleFactor;
	TopValues* curTopValues = pr->localTopValues;

	int n = curTopValues->getCurrentNumberOfPositions( );

	double v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc, elec, hbond, hydrop, vdw, scomp, pgsol, pgsolh, dispe, x, y, z, rmsd;
	double rx, ry, rz;
	int nclashes, r, f, c;
	Matrix transformation;

	double mv, mrv, mrv_ss, mrv_cc, mrv_sc, miv, miv_ss, miv_cc, miv_sc, melec, mhbond, mhydrop, mvdw, mscomp, mpgsol, mpgsolh, mdispe, mx, my, mz, mrmsd = 100000000;
	int mnclashes, mr, mf, mc, rank;
	Matrix mtransformation;

#if defined(__APPLE__)
	double min_rv_ss = 1000000000, min_rv_cc = 1000000000, min_rv_sc = 1000000000;
	double max_rv_ss = -1000000000, max_rv_cc = -1000000000, max_rv_sc = -1000000000;
	double min_iv_ss = 1000000000, min_iv_cc = 1000000000, min_iv_sc = 1000000000;
	double max_iv_ss = -1000000000, max_iv_cc = -1000000000, max_iv_sc = -1000000000;
#else
	double min_rv_ss = 10000000000, min_rv_cc = 10000000000, min_rv_sc = 10000000000;
	double max_rv_ss = -10000000000, max_rv_cc = -10000000000, max_rv_sc = -10000000000;
	double min_iv_ss = 10000000000, min_iv_cc = 10000000000, min_iv_sc = 10000000000;
	double max_iv_ss = -10000000000, max_iv_cc = -10000000000, max_iv_sc = -10000000000;
#endif

	int counter[ rmsdToReport ]; //keep until 20A
	int hitsInRange[ 6 ]; // 0, 0-9, 0-99, 0-999, 0-9999, 0-99999
	int highestPos[ rmsdToReport ];

	int maxRank = n + 1;
	for ( int i = 0; i < rmsdToReport; i++ )
	{
		counter[ i ] = 0;
		highestPos[ i ] = maxRank;
	}

	for ( int i = 0; i < 6; i++ )
		hitsInRange[ i ] = 0;

	f_printf( fp, (char *)"\n# grid spacing = %f angstrom\n# ", 1.0 / ( ( double ) scale_B * curTopValues->getGridSize( ) ) );
	f_printf( fp, (char *)"\n# number of peaks = %d\n# ", n );
	f_printf( fp, (char *)"\n# score scale down factor = %lf\n# ", functionScaleFactor );
	f_printf( fp, (char *)"\n# number of rotation matrices processed = %d ( %0.2lf\% )\n# ", processedRotations, ( processedRotations * 100.0 ) / pr->pri->numberOfRotations );

	fprintf( fp, (char *)"\n# START PEAKS" );

	int numberOfGoodPeaks = 0;
	int highestRank = maxRank;

	while ( n-- )
	{
		curTopValues->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss,
				&iv_cc, &iv_sc, &elec, &hbond, &hydrop, &vdw, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
				&x, &y, &z, &r, &f, &c );

		if ( rv_ss < min_rv_ss ) min_rv_ss = rv_ss;
		if ( rv_cc < min_rv_cc ) min_rv_cc = rv_cc;
		if ( rv_sc < min_rv_sc ) min_rv_sc = rv_sc;

		if ( rv_ss > max_rv_ss ) max_rv_ss = rv_ss;
		if ( rv_cc > max_rv_cc ) max_rv_cc = rv_cc;
		if ( rv_sc > max_rv_sc ) max_rv_sc = rv_sc;

		if ( iv_ss < min_iv_ss ) min_iv_ss = iv_ss;
		if ( iv_cc < min_iv_cc ) min_iv_cc = iv_cc;
		if ( iv_sc < min_iv_sc ) min_iv_sc = iv_sc;

		if ( iv_ss > max_iv_ss ) max_iv_ss = iv_ss;
		if ( iv_cc > max_iv_cc ) max_iv_cc = iv_cc;
		if ( iv_sc > max_iv_sc ) max_iv_sc = iv_sc;

		retrieveTransformation( x, y, z, c, r, f,
				pr,
				rx, ry, rz,
				transformation );

		if ( reDocking )
			rmsd = getRMSD( pr->pri, transformation );
		else
			rmsd = 0.0;
		// MS .. WHAT IS THAT ?? FIXME
		//          rmsd = baseComplex->getRMSD( transformation, unboundLigandAtomList[ c ], unboundLigandInterfaceAtomIndex[ c ] );

		if ( !breakDownScores ) {
			double riv = rv + iv * unrealWeight;
			fprintf( fp, (char *)"\n%6d %16.5lf %16.5lf %16.5lf %16.5lf %16.5lf %20.5lf %20.5lf %20.5lf %20.5lf %16.5lf %10d %20.5lf %20.5lf %20.5lf ",
					( n + 1 ), ( double ) v / functionScaleFactor, ( double ) riv / functionScaleFactor, ( double ) rv_ss / functionScaleFactor,
					( double ) rv_cc / functionScaleFactor, ( double ) iv / functionScaleFactor, ( double ) elec / functionScaleFactor,
					( double ) hbond / functionScaleFactor, ( double ) hydrop / functionScaleFactor, ( double ) scomp / functionScaleFactor,
					vdw, nclashes, pgsol, pgsolh, dispe );

			//	fprintf( fp, (char *)"\n%6d %16.5lf %16.5lf %16.5lf %16.5lf %20.5lf %20.5lf %16.5lf %10d ", ( n + 1 ), ( double ) v / functionScaleFactor,
			//	       ( double ) riv / functionScaleFactor, ( double ) rv / functionScaleFactor, ( double ) iv / functionScaleFactor,
			//	       ( double ) elec / functionScaleFactor, ( double ) hbond / functionScaleFactor, vdw, nclashes );
		} else {
			double riv = skinSkinWeight * rv_ss + coreCoreWeight * rv_cc + skinCoreWeight * iv_sc;
			fprintf( fp, (char *)"\n%6d %16.5lf %16.5lf %16.5lf %16.5lf %16.5lf %16.5lf %16.5lf %16.5lf %20.5lf %20.5lf %20.5lf %20.5lf %16.5lf %10d %20.5lf %20.5lf %20.5lf ",
					( n + 1 ), ( double ) v / functionScaleFactor,
					( double ) riv / functionScaleFactor,
					( double ) rv_ss / functionScaleFactor, ( double ) rv_cc / functionScaleFactor,
					( double ) rv_sc / functionScaleFactor,
					( double ) iv_ss / functionScaleFactor, ( double ) iv_cc / functionScaleFactor,
					( double ) iv_sc / functionScaleFactor,
					( double ) elec / functionScaleFactor,
					( double ) hbond / functionScaleFactor,
					( double ) hydrop / functionScaleFactor,
					( double ) scomp / functionScaleFactor,
					vdw, nclashes, pgsol, pgsolh, dispe );
		}

		for ( int i = 0; i < 3; i++ )
			for ( int j = 0; j < 4; j++ )
				fprintf( fp, (char *)"%9.3f ", transformation.get( i, j ) );
		fprintf( fp, (char *)"%2d %9.2f", c, rmsd );

		if ( mrmsd > rmsd )
		{
			mrmsd = rmsd;
			mv = v;
			mrv = rv; mrv_ss = rv_ss; mrv_cc = rv_cc; mrv_sc = rv_sc;
			miv = iv; miv_ss = iv_ss; miv_cc = iv_cc; miv_sc = iv_sc;
			melec = elec;
			mhbond = hbond;
			mhydrop = hydrop;
			mvdw = vdw;
			mnclashes = nclashes;
			mscomp = scomp;
			mpgsol = pgsol;
			mpgsolh = pgsolh;
			mdispe = dispe;
			mx = x; my = y; mz = z;
			mx = rx; my = ry; mz = rz;
			mr = r; mf = f; mc = c;
			mtransformation = transformation;
			rank = n + 1;
		}

		int rmsdIntF = ( int ) floor( rmsd );
		int rmsdIntC = ( int ) ceil( rmsd );

		//      if ( rmsdInt > 0 ) rmsdInt--;

		if ( rmsdIntF < rmsdToReport )
		{
			counter[ rmsdIntF ]++;
			if ( highestPos[ rmsdIntF ] > n ) highestPos[ rmsdIntF ] = n;
			if ( rmsdIntC <= rmsdGood )
			{
				numberOfGoodPeaks++;
				if ( n < highestRank ) highestRank = n;
			}
		}

		if ( rmsdIntC <= rmsdGood )
		{
			if ( n <= 0 ) hitsInRange[ 0 ]++;
			if ( n <= 9 ) hitsInRange[ 1 ]++;
			if ( n <= 99 ) hitsInRange[ 2 ]++;
			if ( n <= 999 ) hitsInRange[ 3 ]++;
			if ( n <= 9999 ) hitsInRange[ 4 ]++;
			if ( n <= 99999 ) hitsInRange[ 5 ]++;
		}
	}
	fprintf( fp, (char *)"\n# END PEAKS" );


	//   int numberOfGoodPeaks = 0;
	//   int highestRank = maxRank;
	//
	//   for ( int i = 0; i < rmsdGood; i++ )
	//     {
	//      numberOfGoodPeaks += counter[ i ];
	//      if ( highestPos[ i ] < highestRank ) highestRank = highestPos[ i ];
	//     }

	f_printf( fp, (char *)"\n#\n#" );
	for ( int i = 0; i < rmsdToReport; i++ )
		f_printf( fp, (char *)" %d --> %d %d\n# ", i, counter[ i ], ( highestPos[ i ] == maxRank ) ? -1 : highestPos[ i ] );

	f_printf( fp, (char *)"\n#\n#Hits in Range:\n#" );
	for ( int i = 0, j = 1; i < 6; i++ )
	{
		f_printf( fp, (char *)" [%d, %d] --> %d\n#", 1, j, hitsInRange[ i ] );
		j *= 10;
	}

	f_printf( fp, (char *)"\n# good peaks under %d A: count = %d highest rank = %d min RMSD = %f\n# ", rmsdGood, numberOfGoodPeaks, highestRank + 1, mrmsd );

	f_printf( fp, (char *)"\n# best peak: rmsd = %f rank = %d score = %lf ", mrmsd, rank, ( double ) mv / functionScaleFactor );

	if ( breakDownScores )
		fprintf( fp, (char *)"realScore = < skin-skin = %lf core-core = %lf skin-core = %lf >, unrealScore = < skin-skin = %lf core-core = %lf skin-core = %lf >, elecScore = %lf, hbondScore = %lf, hydrophobicityScore = %lf, simpleComplementarityScore = %lf, vdWPotential = %lf, nClashes = %d, pGsol = %lf, pGsolH = %lf, delDispE = %lf ",
				( double ) mrv_ss / functionScaleFactor, ( double ) mrv_cc / functionScaleFactor, ( double ) mrv_sc / functionScaleFactor,
				( double ) miv_ss / functionScaleFactor, ( double ) miv_cc / functionScaleFactor, ( double ) miv_sc / functionScaleFactor,
				( double ) melec / functionScaleFactor,
				( double ) mhbond / functionScaleFactor,
				( double ) mhydrop / functionScaleFactor,
				( double ) mscomp / functionScaleFactor,
				( double ) mvdw, ( int ) nclashes, mpgsol, mpgsolh, mdispe );
	else
		fprintf( fp, (char *)"realScore = %lf unrealScore = %lf elecScore = %lf hbondScore = %lf  hydrophobicityScore = %lf simpleComplementarityScore = %lf vdWPotential = %lf nClashes = %d pGsol = %lf pGsolH = %lf delDispE = %lf ",
				( double ) mrv / functionScaleFactor, ( double ) miv / functionScaleFactor,
				( double ) melec / functionScaleFactor, ( double ) mhbond / functionScaleFactor,
				( double ) mhydrop / functionScaleFactor, ( double ) mscomp / functionScaleFactor,
				( double ) mvdw, ( int ) nclashes, mpgsol, mpgsolh, mdispe );


	fprintf( fp, (char *)"transformation matrix = [ " );

	for ( int i = 0; i < 4; i++ )
	{
		fprintf( fp, (char *)"< " );

		for ( int j = 0; j < 4; j++ )
			fprintf( fp, (char *)"%f ", mtransformation.get( i, j ) );

		fprintf( fp, (char *)"> " );
	}

	fprintf( fp, (char *)"], conformationIndex = %d\n# ",  mc );

	if ( breakDownScores )
	{
		min_rv_ss /= functionScaleFactor;
		min_rv_cc /= functionScaleFactor;
		min_rv_sc /= functionScaleFactor;

		max_rv_ss /= functionScaleFactor;
		max_rv_cc /= functionScaleFactor;
		max_rv_sc /= functionScaleFactor;

		min_iv_ss /= functionScaleFactor;
		min_iv_cc /= functionScaleFactor;
		min_iv_sc /= functionScaleFactor;

		max_iv_ss /= functionScaleFactor;
		max_iv_cc /= functionScaleFactor;
		max_iv_sc /= functionScaleFactor;

		double max_abs_rv_ss;

		if ( fabs( max_rv_ss ) > fabs( min_rv_ss ) ) max_abs_rv_ss = fabs( max_rv_ss ) / 100;
		else max_abs_rv_ss = fabs( min_rv_ss ) / 100;

		fprintf( fp, (char *)"\n#\n# max_rv_ss = %lf ( %lf ), min_rv_ss = %lf ( %lf )",
				( double ) max_rv_ss, ( double ) max_rv_ss / max_abs_rv_ss, ( double ) min_rv_ss, ( double ) min_rv_ss / max_abs_rv_ss  );
		fprintf( fp, (char *)"\n# max_rv_cc = %lf ( %lf ), min_rv_cc = %lf ( %lf )",
				( double ) max_rv_cc, ( double ) max_rv_cc / max_abs_rv_ss, ( double ) min_rv_cc, ( double ) min_rv_cc / max_abs_rv_ss  );
		fprintf( fp, (char *)"\n# max_iv_sc = %lf ( %lf ), min_iv_sc = %lf ( %lf )\n",
				( double ) max_iv_sc, ( double ) max_iv_sc / max_abs_rv_ss, ( double ) min_iv_sc, ( double ) min_iv_sc / max_abs_rv_ss  );

		fprintf( fp, (char *)"\n# max_rv_sc = %lf ( %lf ), min_rv_sc = %lf ( %lf )",
				( double ) max_rv_sc, ( double ) max_rv_sc / max_abs_rv_ss, ( double ) min_rv_sc, ( double ) min_rv_sc / max_abs_rv_ss  );
		fprintf( fp, (char *)"\n# max_iv_ss = %lf ( %lf ), min_iv_ss = %lf ( %lf )",
				( double ) max_iv_ss, ( double ) max_iv_ss / max_abs_rv_ss, ( double ) min_iv_ss, ( double ) min_iv_ss / max_abs_rv_ss  );
		fprintf( fp, (char *)"\n# max_iv_cc = %lf ( %lf ), min_iv_cc = %lf ( %lf )\n\n",
				( double ) max_iv_cc, ( double ) max_iv_cc / max_abs_rv_ss, ( double ) min_iv_cc, ( double ) min_iv_cc / max_abs_rv_ss  );
	}

	//  printf( "m_CurMin = %d\n", curTopValues->getCurMin( ) );
}



void printUntransformedScore( FILE *fp, PARAMS *pr )
{
	bool reDocking = ( pr->pri->nbRMSDAtoms > 0 ) ? true : false;
	bool breakDownScores = ( pr->pri->breakDownScores == 1 ) ? true : false;

	float scale_B = pr->scaleB;
	double functionScaleFactor = pr->functionScaleFactor;
	TopValues* curTopValues = pr->localTopValues;

	int n = curTopValues->getCurrentNumberOfPositions( );

	double v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc, elec, hbond, hydrop, vdw, scomp, pgsol, pgsolh, dispe, x, y, z, rmsd;
	double mv, mrv, mrv_ss, mrv_cc, mrv_sc, miv, miv_ss, miv_cc, miv_sc, melec, mhbond, mhydrop, mvdw, mscomp, mpgsol, mpgsolh, mdispe, mrmsd;
	double md;
	double rx, ry, rz;
	double mrx, mry, mrz;
	int nclashes, mnclashes, r, f, c;
	Matrix transformation;

	f_printf( fp, (char *)"\n# grid spacing = %f angstrom\n# ", 1.0 / ( ( double ) scale_B * curTopValues->getGridSize( ) ) );
	f_printf( fp, (char *)"\n# score scale down factor = %lf\n# ", functionScaleFactor );

	fprintf( fp, (char *)"\n# START UNTRANSFORMED SCORE" );

	md = 10000000000000.0;

	while ( n-- )
	{
		curTopValues->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss,
				&iv_cc, &iv_sc, &elec, &hbond, &hydrop, &vdw, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
				&x, &y, &z, &r, &f, &c );

		retrieveTransformation( x, y, z, c, r, f,
				pr,
				rx, ry, rz,
				transformation );

		double d = rx * rx + ry * ry + rz * rz;

		if ( d > md ) continue;

		md = d;

		if ( reDocking )
			rmsd = getRMSD( pr->pri, transformation );
		else
			rmsd = 0.0;

		mrmsd = rmsd;
		mv = v;
		mrv = rv; mrv_ss = rv_ss; mrv_cc = rv_cc; mrv_sc = rv_sc;
		miv = iv; miv_ss = iv_ss; miv_cc = iv_cc; miv_sc = iv_sc;
		melec = elec;
		mhbond = hbond;
		mhydrop = hydrop;
		mvdw = vdw;
		mnclashes = nclashes;
		mscomp = scomp;
		mpgsol = pgsol;
		mpgsolh = pgsolh;
		mdispe = dispe;
		mrx = rx; mry = ry; mrz = rz;
	}

	if ( breakDownScores )
	{
		fprintf( fp, (char *)"\nShape Complementarity Score:\t total = %lf,\n\t\t\t\t real = < skin-skin = %lf core-core = %lf skin-core = %lf >,\n\t\t\t\t imaginary = < skin-skin = %lf core-core = %lf skin-core = %lf >",
				( double ) mv / functionScaleFactor,
				( double ) mrv_ss / functionScaleFactor, ( double ) mrv_cc / functionScaleFactor, ( double ) mrv_sc / functionScaleFactor,
				( double ) miv_ss / functionScaleFactor, ( double ) miv_cc / functionScaleFactor, ( double ) miv_sc / functionScaleFactor );

		//      fprintf( fp, (char *)"\nRMSD = %lf ( rx = %lf, ry = %lf, rz = %lf )", mrmsd, mrx, mry, mrz );
	}
	else
	{
		fprintf( fp, (char *)"\nShape Complementarity Score:\t total = %lf, real = %lf, imaginary = %lf",
				( double ) mv / functionScaleFactor,
				( double ) mrv / functionScaleFactor,
				( double ) miv / functionScaleFactor );
	}

	fprintf( fp, (char *)"\nElectrostatics Score:\t\t total ( real ) = %lf",
			( double ) melec / functionScaleFactor );

	fprintf( fp, (char *)"\nHydrogen Bonding Score:\t\t total ( real ) = %lf",
			( double ) mhbond / functionScaleFactor );

	fprintf( fp, (char *)"\nHydrophobicity Score:\t\t total ( real ) = %lf",
			( double ) mhydrop / functionScaleFactor );

	fprintf( fp, (char *)"\nSimple Complementarity Score:\t\t total ( real ) = %lf",
			( double ) mscomp / functionScaleFactor );

	fprintf( fp, (char *)"\nClash Filter Score:\t\t total ( real ) = %d",
			mnclashes);

	fprintf( fp, (char *)"\nLennard-Jones Filter Score:\t\t total ( real ) = %lf",
			( double ) mvdw / functionScaleFactor );

	fprintf( fp, (char *)"\nPseudo Gsol Filter Score:\t\t total ( real ) = %lf",
			( double )  mpgsol / functionScaleFactor );

	fprintf( fp, (char *)"\n# END UNTRANSFORMED SCORE" );
}


int nextRotation( PARAMS *pr )
{
	int nextRot = -1;

	pthread_mutex_lock( &rotLock );
	if ( pr->pri->pruneAngle > 0  )
	{
		while ( ( curRotation < maxRotation ) && ( pr->pri->rotGraph[ curRotation ].badNode ) ) curRotation++;
	}
	if ( curRotation < maxRotation )
	{
		nextRot = curRotation++;
		processedRotations++;
	}
	pthread_mutex_unlock( &rotLock );

	return nextRot;
}


//int nextRotation( PARAMS *pr )
//{
//   int nextRot = -1, k = pr->threadID;
//
//   pthread_mutex_lock( &rotLock );
//
//   if ( pr->pri->pruneAngle > 0  )
//     {
//       while ( ( curRotation[ k ] < maxRotation[ k ] ) && ( pr->pri->rotGraph[ curRotation[ k ] ].badNode ) ) curRotation[ k ]++;
//     }
//   if ( curRotation[ k ] < maxRotation[ k ] )
//     {
//       nextRot = curRotation[ k ]++;
//       processedRotations++;
//     }
//   else
//     {
//       for ( int i = 0; i < pr->pri->numThreads; i++ )
//         {
//           if ( pr->pri->pruneAngle > 0  )
//             {
//               while ( ( curRotation[ i ] < maxRotation[ i ] ) && ( pr->pri->rotGraph[ curRotation[ i ] ].badNode ) ) curRotation[ i ]++;
//             }
//
//          if ( curRotation[ i ] < maxRotation[ i ] )
//             {
//               nextRot = curRotation[ i ]++;
//               processedRotations++;
//               break;
//             }
//         }
//     }
//
//   pthread_mutex_unlock( &rotLock );
//
//   return nextRot;
//}


void markBadRotGraphNodes( int k, PARAMS *pr )
{
	if ( pr->pri->pruneAngle <= 0  ) return;

	pthread_mutex_lock( &rotLock );
	for ( int i = pr->pri->rotGraph[ k ].startNeighbors; i <= pr->pri->rotGraph[ k ].endNeighbors; i++ )
		pr->pri->rotGraph[ pr->pri->rotNeighbors[ i ] ].badNode = true;
	pthread_mutex_unlock( &rotLock );
}


void initRotationServer( int numRot, int numThreads )
{
	pthread_mutex_init( &rotLock, NULL );
	curRotation = 0;
	processedRotations = 0;
	maxRotation = numRot;
}


//void initRotationServer( int numRot, int numThreads )
//{
//  pthread_mutex_init( &rotLock, NULL );
//
//  curRotation = ( int * ) malloc( numThreads * sizeof( int ) );
//  maxRotation = ( int * ) malloc( numThreads * sizeof( int ) );
//
//  int len = ( int ) floor( ( 1.0 * numRot ) / numThreads );
//  int res = numRot - ( len * numThreads );
//
//  processedRotations = 0;
//
//  for ( int i = 0, j = 0; i < numThreads; i++ )
//    {
//      curRotation[ i ] = j;
//
//      j += len;
//      if ( res-- > 0 ) j++;
//
//      maxRotation[ i ] = j;
//    }
//}

/*
 *   Rotates rotGrid to grid, based on the angles theta
 */
void rotateGrid( FFTW_complex *rotGrid, FFTW_complex *grid,  int n, double theta1, double theta2, double theta3 )
{
	int nn = n >> 1, odd = n & 1;
	double eps = 0.001;
	double ofs;

	if ( odd ) ofs = 0;
	else ofs = 0.5;

	// 3D iterations, assuming that 0 is the cartesian origin
	// rotates each point of rotGrid?
	for ( int z = -nn; z < nn + odd; z++ )
		for ( int y = -nn; y < nn + odd; y++ )
			for ( int x = -nn; x < nn + odd; x++ )
			{
				double xf = x + ofs, yf = y + ofs, zf = z + ofs;

				rotatePointAboutOrigin( &xf, &yf, &zf, theta1, theta2, theta3 );

				xf -= ofs; yf -= ofs; zf -= ofs;

				int xl = floor( xf ), xh = ceil( xf );
				int yl = floor( yf ), yh = ceil( yf );
				int zl = floor( zf ), zh = ceil( zf );

				int i = ( ( z + nn ) * n + ( y + nn ) ) * n + ( x + nn );

				rotGrid[ i ][ 0 ] = rotGrid[ i ][ 1 ] = 0;

				if ( ( xl < -nn ) || ( xh >= nn + odd )
						|| ( yl < -nn ) || ( yh >= nn + odd )
						|| ( zl < -nn ) || ( zh >= nn + odd ) )
					continue;

				double d2[ 8 ], ds = 0;

				// these 2 sets of loops could be merged, and d2 can be only a double variable
				for ( int zz = zl, c = 0; zz <= zh; zz++ )
					for ( int yy = yl; yy <= yh; yy++ )
						for ( int xx = xl; xx <= xh; xx++, c++ )
						{
							d2[ c ] = ( xx - xf ) * ( xx - xf ) + ( yy - yf ) * ( yy - yf ) + ( zz - zf ) * ( zz - zf );
							d2[ c ] = pow( d2[ c ], 6 );
							if ( d2[ c ] >= eps ) ds += 1 / d2[ c ];
						}

				//         double dd = 10000000000000000.0;

				for ( int zz = zl, c = 0; zz <= zh; zz++ )
					for ( int yy = yl; yy <= yh; yy++ )
						for ( int xx = xl; xx <= xh; xx++, c++ )
						{
							int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );

							//                if ( d2[ c ] < dd )
							//                  {
							//                   rotGrid[ i ][ 0 ] = grid[ j ][ 0 ];
							//                   rotGrid[ i ][ 1 ] = grid[ j ][ 1 ];
							//                   dd = d2[ c ];
							//                  }

							if ( d2[ c ] < eps )
							{
								rotGrid[ i ][ 0 ] = grid[ j ][ 0 ];
								rotGrid[ i ][ 1 ] = grid[ j ][ 1 ];
								xx = xh + 1; yy = yh + 1; zz = zh + 1;
							}
							else
							{
								rotGrid[ i ][ 0 ] += ( ( 1 / d2[ c ] ) / ds ) * grid[ j ][ 0 ];
								rotGrid[ i ][ 1 ] += ( ( 1 / d2[ c ] ) / ds ) * grid[ j ][ 1 ];
							}
						}
			}
}


int collectNonzeroGridCells( NONZERO_GRIDCELLS **gridCells, FFTW_complex *grid, int n )
{
	int numNonzeroCells = 0;
	double eps = 0.0;

	for ( int c = 0; c < n * n * n; c++ )
		if ( ( fabs( ( double ) grid[ c ][ 0 ] ) > eps ) || ( fabs( ( double ) grid[ c ][ 1 ] ) > eps ) ) numNonzeroCells++;

	*gridCells = ( NONZERO_GRIDCELLS * ) malloc( numNonzeroCells * sizeof( NONZERO_GRIDCELLS ) );

	if ( *gridCells == NULL ) return 0;

	int nn = n >> 1, odd = n & 1;

	for ( int z = -nn, c = 0, k = 0; z < nn + odd; z++ )
		for ( int y = -nn; y < nn + odd; y++ )
			for ( int x = -nn; x < nn + odd; x++, c++ )
				if ( ( fabs( ( double ) grid[ c ][ 0 ] ) > eps ) || ( fabs( ( double ) grid[ c ][ 1 ] ) > eps ) )
				{
					( *gridCells )[ k ].x = x;
					( *gridCells )[ k ].y = y;
					( *gridCells )[ k ].z = z;

					( *gridCells )[ k ].v[ 0 ] = grid[ c ][ 0 ];
					( *gridCells )[ k ].v[ 1 ] = grid[ c ][ 1 ];

					k++;
				}

	return numNonzeroCells;
}


int collectNonzeroElecGridCells( NONZERO_GRIDCELLS **gridCells, FFTW_DATA_TYPE *grid, int n )
{
	int numNonzeroCells = 0;
	double eps = 0.0;

	for ( int c = 0; c < n * n * n; c++ )
		if ( fabs( ( double ) grid[ c ] ) > eps ) numNonzeroCells++;

	*gridCells = ( NONZERO_GRIDCELLS * ) malloc( numNonzeroCells * sizeof( NONZERO_GRIDCELLS ) );

	if ( *gridCells == NULL ) return 0;

	int nn = n >> 1, odd = n & 1;

	for ( int z = -nn, c = 0, k = 0; z < nn + odd; z++ )
		for ( int y = -nn; y < nn + odd; y++ )
			for ( int x = -nn; x < nn + odd; x++, c++ )
				if ( fabs( ( double ) grid[ c ] ) > eps )
				{
					( *gridCells )[ k ].x = x;
					( *gridCells )[ k ].y = y;
					( *gridCells )[ k ].z = z;

					( *gridCells )[ k ].v[ 0 ] = grid[ c ];

					k++;
				}

	return numNonzeroCells;
}


int collectNonzeroHbondGridCells( NONZERO_GRIDCELLS **gridCells, FFTW_complex *grid, int n )
{
	return collectNonzeroGridCells( gridCells, grid, n );
}


int collectNonzeroHydrophobicityGridCells( NONZERO_GRIDCELLS **gridCells, FFTW_complex *grid, int n )
{
	return collectNonzeroGridCells( gridCells, grid, n );
}


int collectNonzeroSimpleComplementarityGridCells( NONZERO_GRIDCELLS **gridCells, FFTW_complex *grid, int n )
{
	return collectNonzeroGridCells( gridCells, grid, n );
}


double findNonzeroGridCellsRadius( NONZERO_GRIDCELLS *gridCells, int numNonzeroCells, int n )
{
	double maxDist2 = 0;
	double nc = ( n - 1 ) / 2.0;
	int nn = n >> 1;

	for ( int k = 0; k < numNonzeroCells; k++ )
	{
		int x = gridCells[ k ].x, y = gridCells[ k ].y, z = gridCells[ k ].z;
		double d2 = ( x + nn - nc ) * ( x + nn - nc ) + ( y + nn - nc ) * ( y + nn - nc ) + ( z + nn - nc ) * ( z + nn - nc );

		if ( d2 > maxDist2 ) maxDist2 = d2;
	}

	return sqrt( maxDist2 );
}


void findMovingMolMinMaxRadius( NONZERO_GRIDCELLS *gridCells, int numNonzeroCells, int n , double *minRad, double *maxRad )
{
	double maxDist2 = 0, minDist2 = ( n << 2 ) * ( n << 2 );
	//   double eps = 0.00001;
	double nc = ( n - 1 ) / 2.0;
	int nn = n >> 1;
	double maxRe = 0, maxIm = 0;

	for ( int k = 0; k < numNonzeroCells; k++ )
	{
		if ( fabs( gridCells[ k ].v[ 0 ] ) > maxRe ) maxRe = fabs( gridCells[ k ].v[ 0 ] );
		if ( fabs( gridCells[ k ].v[ 1 ] ) > maxIm ) maxIm = fabs( gridCells[ k ].v[ 1 ] );
	}

	for ( int k = 0; k < numNonzeroCells; k++ )
	{
		int x = gridCells[ k ].x, y = gridCells[ k ].y, z = gridCells[ k ].z;
		double d2 = ( x + nn - nc ) * ( x + nn - nc ) + ( y + nn - nc ) * ( y + nn - nc ) + ( z + nn - nc ) * ( z + nn - nc );

		if ( d2 > maxDist2 ) maxDist2 = d2;
		if ( ( fabs( gridCells[ k ].v[ 0 ] ) > maxRe / 10 ) && ( fabs( gridCells[ k ].v[ 1 ] ) > maxIm / 100 ) && ( d2 < minDist2 ) )
		{
			minDist2 = d2;
			//          printf( "gridCells[ k ].v[ 0 ] = %lf, gridCells[ k ].v[ 1 ] = %lf, d2 = %lf\n", gridCells[ k ].v[ 0 ], gridCells[ k ].v[ 1 ], d2 );
		}
	}

	//   exit( 1 );

	*minRad = floor( sqrt( minDist2 ) );
	*maxRad = ceil( sqrt( maxDist2 ) );
}


void findMovingMolMinMaxRadius( int M, double* x, double* y, double* z, float *r, char *type, int n, double scale, double *minRad, double *maxRad )
{
	*minRad = -1.0;
	*maxRad = 0.0;

	for ( int s = 0; s < M; s++ )
		if ( type[ s ] == 'E' )
		{
			double d = sqrt( x[ s ] * x[ s ] + y[ s ] * y[ s ] + z[ s ] * z[ s ] ) + r[ s ];

			if ( d > *maxRad ) *maxRad = d;
			if ( ( *minRad == -1.0 ) || ( d < *minRad ) ) *minRad = d;
		}

	*minRad = floor( ( *minRad ) * scale * ( n - 1 ) );
	*maxRad = ceil( ( *maxRad ) * scale * ( n - 1 ) );
}


void flipGrid( NONZERO_GRIDCELLS *gridCells, int numNonzeroCells )
{
	for ( int k = 0; k < numNonzeroCells; k++ )
	{
		gridCells[ k ].x = -gridCells[ k ].x;
		gridCells[ k ].y = -gridCells[ k ].y;
		gridCells[ k ].z = -gridCells[ k ].z;
	}
}

/*
 *   Rotates the rotGrid to grid, considering only the nonzero cells (optimization)
 *   and a rotation matrix, instead of specifying the angles.
 *   CleanupOnly only cleanups the rotateGrid instead of rotating the grid.
 */

void rotateGrid( FFTW_complex *rotGrid, int n, NONZERO_GRIDCELLS *grid, int numNonzeroCells, Matrix &rotMat, bool cleanupRotGridFirst, bool cleanupOnly )
{
	if ( cleanupRotGridFirst )
	{
		for ( int i = 0; i < n * n * n; i++ )
			rotGrid[ i ][ 0 ] = rotGrid[ i ][ 1 ] = 0;
	}

	int nn = n >> 1, odd = n & 1;
	double eps = 0.000001, pw = 3.0;
	double ofs;

	if ( odd ) ofs = 0;
	else ofs = 0.5;

	if ( cleanupOnly )
	{
		for ( int i = 0; i < numNonzeroCells; i++ )
		{
			Vector oldPos( ( float ) grid[ i ].x + ofs, ( float ) grid[ i ].y + ofs, ( float ) grid[ i ].z + ofs, 1.0 );
			Vector newPos = rotMat * oldPos;

			double xf = newPos[ 0 ] - ofs,
			       yf = newPos[ 1 ] - ofs,
			       zf = newPos[ 2 ] - ofs;

			int xl = floor( xf ), xh = ceil( xf );
			int yl = floor( yf ), yh = ceil( yf );
			int zl = floor( zf ), zh = ceil( zf );

			if ( ( xl < -nn ) || ( xh >= nn + odd )
					|| ( yl < -nn ) || ( yh >= nn + odd )
					|| ( zl < -nn ) || ( zh >= nn + odd ) )
				continue;

			for ( int zz = zl; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++ )
					{
						int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );
						rotGrid[ j ][ 0 ] = rotGrid[ j ][ 1 ] = 0;
					}
		}
	}
	else
	{
		for ( int i = 0; i < numNonzeroCells; i++ )
		{
			//         printf( "\n< %d, %d, %d > ==> < %f, %f >", grid[ i ].x, grid[ i ].y, grid[ i ].z,
			//                                                    ( float ) grid[ i ].v[ 0 ], ( float ) grid[ i ].v[ 1 ] );
			Vector oldPos( grid[ i ].x + ofs, grid[ i ].y + ofs, grid[ i ].z + ofs, 1.0 );
			Vector newPos = rotMat * oldPos;

			double xf = newPos[ 0 ] - ofs,
			       yf = newPos[ 1 ] - ofs,
			       zf = newPos[ 2 ] - ofs;

			int xl = floor( xf ), xh = ceil( xf );
			int yl = floor( yf ), yh = ceil( yf );
			int zl = floor( zf ), zh = ceil( zf );

			if ( ( xl < -nn ) || ( xh >= nn + odd )
					|| ( yl < -nn ) || ( yh >= nn + odd )
					|| ( zl < -nn ) || ( zh >= nn + odd ) )
				continue;

			//         printf( " :: < %f, %f, %f >", xf, yf, zf );

			double d2[ 8 ], ds = 0;

			// Can be merged
			for ( int zz = zl, c = 0; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++, c++ )
					{
						d2[ c ] = ( xx - xf ) * ( xx - xf ) + ( yy - yf ) * ( yy - yf ) + ( zz - zf ) * ( zz - zf );
						d2[ c ] = pow( d2[ c ], pw );
						if ( d2[ c ] < eps ) d2[ c ] = eps;
						ds += 1 / d2[ c ];
					}

			for ( int zz = zl, c = 0; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++, c++ )
					{
						int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );

						rotGrid[ j ][ 0 ] += ( ( 1 / d2[ c ] ) / ds ) * grid[ i ].v[ 0 ];
						rotGrid[ j ][ 1 ] += ( ( 1 / d2[ c ] ) / ds ) * grid[ i ].v[ 1 ];
					}
		}
	}
}

void rotateElecGrid( FFTW_DATA_TYPE *rotGrid, int n, NONZERO_GRIDCELLS *grid, int numNonzeroCells, Matrix &rotMat, bool cleanupRotGridFirst, bool cleanupOnly )
{
	if ( cleanupRotGridFirst )
	{
		for ( int i = 0; i < n * n * n; i++ )
			rotGrid[ i ] = 0;
	}

	int nn = n >> 1, odd = n & 1;
	double eps = 0.000001, pw = 3.0;
	double ofs;

	if ( odd ) ofs = 0;
	else ofs = 0.5;

	if ( cleanupOnly )
	{
		for ( int i = 0; i < numNonzeroCells; i++ )
		{
			Vector oldPos( ( float ) grid[ i ].x + ofs, ( float ) grid[ i ].y + ofs, ( float ) grid[ i ].z + ofs, 1.0 );
			Vector newPos = rotMat * oldPos;

			double xf = newPos[ 0 ] - ofs,
			       yf = newPos[ 1 ] - ofs,
			       zf = newPos[ 2 ] - ofs;

			int xl = floor( xf ), xh = ceil( xf );
			int yl = floor( yf ), yh = ceil( yf );
			int zl = floor( zf ), zh = ceil( zf );

			if ( ( xl < -nn ) || ( xh >= nn + odd )
					|| ( yl < -nn ) || ( yh >= nn + odd )
					|| ( zl < -nn ) || ( zh >= nn + odd ) )
				continue;

			for ( int zz = zl; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++ )
					{
						int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );
						rotGrid[ j ] = 0;
					}
		}
	}
	else
	{
		for ( int i = 0; i < numNonzeroCells; i++ )
		{
			Vector oldPos( grid[ i ].x + ofs, grid[ i ].y + ofs, grid[ i ].z + ofs, 1.0 );
			Vector newPos = rotMat * oldPos;

			double xf = newPos[ 0 ] - ofs,
			       yf = newPos[ 1 ] - ofs,
			       zf = newPos[ 2 ] - ofs;

			int xl = floor( xf ), xh = ceil( xf );
			int yl = floor( yf ), yh = ceil( yf );
			int zl = floor( zf ), zh = ceil( zf );

			if ( ( xl < -nn ) || ( xh >= nn + odd )
					|| ( yl < -nn ) || ( yh >= nn + odd )
					|| ( zl < -nn ) || ( zh >= nn + odd ) )
				continue;

			double d2[ 8 ], ds = 0;

			for ( int zz = zl, c = 0; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++, c++ )
					{
						d2[ c ] = ( xx - xf ) * ( xx - xf ) + ( yy - yf ) * ( yy - yf ) + ( zz - zf ) * ( zz - zf );
						d2[ c ] = pow( d2[ c ], pw );
						if ( d2[ c ] < eps ) d2[ c ] = eps;
						ds += 1 / d2[ c ];
					}

			for ( int zz = zl, c = 0; zz <= zh; zz++ )
				for ( int yy = yl; yy <= yh; yy++ )
					for ( int xx = xl; xx <= xh; xx++, c++ )
					{
						int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );

						rotGrid[ j ] += ( ( 1 / d2[ c ] ) / ds ) * grid[ i ].v[ 0 ];
					}
		}
	}
}



void rotateHbondGrid( FFTW_complex *rotGrid, int n, NONZERO_GRIDCELLS *grid, int numNonzeroCells, Matrix &rotMat, bool cleanupRotGridFirst, bool cleanupOnly )
{
	rotateGrid( rotGrid, n, grid, numNonzeroCells, rotMat, cleanupRotGridFirst, cleanupOnly );
}



void rotateHydrophobicityGrid( FFTW_complex *rotGrid, int n, NONZERO_GRIDCELLS *grid, int numNonzeroCells, Matrix &rotMat, bool cleanupRotGridFirst, bool cleanupOnly )
{
	rotateGrid( rotGrid, n, grid, numNonzeroCells, rotMat, cleanupRotGridFirst, cleanupOnly );
}


void rotateSimpleComplementarityGrid( FFTW_complex *rotGrid, int n, NONZERO_GRIDCELLS *grid, int numNonzeroCells, Matrix &rotMat, bool cleanupRotGridFirst, bool cleanupOnly )
{
	rotateGrid( rotGrid, n, grid, numNonzeroCells, rotMat, cleanupRotGridFirst, cleanupOnly );
}


void shapeScoreSingleConfiguration( FFTW_complex *gridA, FFTW_complex *gridB, int n, int xt, int yt, int zt,
		double *scoreSS, double *scoreCC, double *scoreSC )
{
	int nn = n >> 1, odd = n & 1;
	*scoreSS = *scoreCC = *scoreSC = 0;

	if ( xt >= nn + odd ) xt -= n;
	if ( yt >= nn + odd ) yt -= n;
	if ( zt >= nn + odd ) zt -= n;

	for ( int z = -nn, i = 0; z < nn + odd; z++ )
		for ( int y = -nn; y < nn + odd; y++ )
			for ( int x = -nn; x < nn + odd; x++, i++ )
			{
				int xx = x - xt, yy = y - yt, zz = z - zt;

				xx = - xx;
				yy = - yy;
				zz = - zz;

				if ( ( xx < -nn ) || ( xx >= nn + odd )
						|| ( yy < -nn ) || ( yy >= nn + odd )
						|| ( zz < -nn ) || ( zz >= nn + odd ) )
					continue;

				int j = ( ( zz + nn ) * n + ( yy + nn ) ) * n + ( xx + nn );

				( *scoreSS ) += gridA[ i ][ 0 ] * gridB[ j ][ 0 ];
				( *scoreSC ) += gridA[ i ][ 0 ] * gridB[ j ][ 1 ] + gridA[ i ][ 1 ] * gridB[ j ][ 0 ];
				( *scoreCC ) += gridA[ i ][ 1 ] * gridB[ j ][ 1 ];
			}
}


typedef union
{
	unsigned long ix;
	float fx;
} UL_F_union;


inline float invSqrt( float x )
{
	if ( x < 0 ) return INFINITY;

	volatile UL_F_union xx;

	xx.ix = 0;
	xx.fx = x;
	xx.ix = ( 0xBE6EB50CUL - xx.ix ) >> 1;

	return ( 0.5f * xx.fx * ( 3.0f - x * xx.fx * xx.fx ) );
}

inline float fastExp( float x )
{
	if ( x > 85 ) return INFINITY;
	if ( x < -85 ) return 0.0;

	volatile UL_F_union xx;

	// 12102203.16156148555068672305845f = ( 2^23 ) / ln(2);
	unsigned long i = ( unsigned long ) ( x * 12102203.16156f );

	// 361007 = ( 0.08607133 / 2 ) * ( 2^23 )
	xx.ix = i + ( 127L << 23 ) - 361007;

	return xx.fx;
}



void elecScoreSingleConfiguration( int nA, double *xA, double *yA, double *zA, float *rA, float *qA, char *typeA,
		int nB, double *xB, double *yB, double *zB, float *rB, float *qB, char *typeB, Matrix transMat,
		double *elecPos, double *elecNeg )
{
	*elecPos = *elecNeg = 0;

	//  int numStaticAtoms = 0;
	//
	//  for ( int i = 0; i < nA; i++ )
	//     if ( typeA[ i ] == 'I' ) numStaticAtoms++;
	//
	//  printf( "\nnumStaticAtoms = %d\n", numStaticAtoms );
	//
	//  for ( int i = 0; i < nA; i++ )
	//    if ( typeA[ i ] == 'I' )
	//       printf( "\n< %lf, %lf, %lf, %lf, %lf >\n", xA[ i ], yA[ i ], zA[ i ], ( double ) qA[ i ], ( double ) rA[ i ] );
	//
	//  printf( "\nnumCentersB = %d\n", nB );
	//
	//  for ( int i = 0; i < nB; i++ )
	//       printf( "\n< %lf, %lf, %lf, %lf, %lf >\n", xB[ i ], yB[ i ], zB[ i ], ( double ) qB[ i ], ( double ) rB[ i ] );
	//
	//  printf( "\n\n" );


	for ( int i = 0; i < nB; i++ )
	{
		Vector oldPos( xB[ i ], yB[ i ], zB[ i ], 1.0 );
		Vector newPos = transMat * oldPos;
		double xxB = newPos[ 0 ], yyB = newPos[ 1 ], zzB = newPos[ 2 ];

		for ( int j = 0; j < nA; j++ )
			if ( typeA[ j ] == 'I' )
			{
				float dist2 = ( xA[ j ] - xxB ) * ( xA[ j ] - xxB )
					+ ( yA[ j ] - yyB ) * ( yA[ j ] - yyB )
					+ ( zA[ j ] - zzB ) * ( zA[ j ] - zzB );

				//            double dist = sqrt( dist2 ); //1.0 / invSqrt( dist2 );

				if ( dist2 < ( 0.8 * ( rA[ j ] + rB[ i ] ) ) * ( 0.8 * ( rA[ j ] + rB[ i ] ) ) )
					dist2 = ( 0.8 * ( rA[ j ] + rB[ i ] ) ) * ( 0.8 * ( rA[ j ] + rB[ i ] ) );
				//            if ( dist < 1.0 ) dist = 1.0;

				//            double cons;
				//
				//            if ( dist <= 6 ) cons = 4;
				//            else if ( dist < 8 ) cons = 38 * dist - 224;
				//                 else cons = 80;
				//
				//            double val = ( qA[ j ] * qB[ i ] ) / ( cons * dist );

				//            double val = ( qA[ j ] * qB[ i ] ) / dist2;
				//
				//            if ( val >= 0 ) ( *elecPos ) += val;
				//            else ( *elecNeg ) += val;

				( *elecPos ) += ( 1.0 / dist2 );
			}
	}
}


bool applyRotations( int threadID, int confID, float *rotations,
		int numberOfRotations, int numberOfPositions,
		double *xkBOrig, double *ykBOrig, double *zkBOrig, float *rkBOrig,
		double *xkB, double *ykB, double *zkB, float *rkB, int numCentersB,
		bool rotateVolume,
		NONZERO_GRIDCELLS *gridBCells, int numNonzeroGridBCells,
		NONZERO_GRIDCELLS *gridBCells_01, int numNonzeroGridBCells_01,
		NONZERO_GRIDCELLS *gridBCells_10, int numNonzeroGridBCells_10,
		NONZERO_GRIDCELLS *gridBCells_11, int numNonzeroGridBCells_11,
		NONZERO_GRIDCELLS *elecGridBCells, int numNonzeroElecGridBCells,
		NONZERO_GRIDCELLS *hbondGridBCells, int numNonzeroHbondGridBCells,
		NONZERO_GRIDCELLS *hydrophobicityGridBCells, int numNonzeroHydrophobicityGridBCells,
		NONZERO_GRIDCELLS *hydrophobicityTwoGridBCells, int numNonzeroHydrophobicityTwoGridBCells,
		NONZERO_GRIDCELLS *simpleComplementarityGridBCells, int numNonzeroSimpleComplementarityGridBCells,
		int gridSize, int numFreq, int interpFuncExtent,
		double alpha, double blobbiness,
		double skinSkinWeight, double coreCoreWeight, double skinCoreWeight,
		double realSCWeight, double imaginarySCWeight,
		double elecScale, double elecRadiusInGrids,
		double hbondWeight, double hbondDistanceCutoff,
		double hydrophobicityWeight, double hydroPhobicPhobicWeight, double hydroPhobicPhilicWeight, double hydroPhilicPhilicWeight,
		double simpleShapeWeight, double simpleChargeWeight,
		double scaleA,
		FFTW_complex* centerFrequenciesA, FFTW_complex* centerElecFrequenciesA, FFTW_complex* centerHbondFrequenciesA,
		FFTW_complex* centerHydrophobicityFrequenciesA, FFTW_complex* centerHydrophobicityTwoFrequenciesA, FFTW_complex* centerSimpleComplementarityFrequenciesA,
		FFTW_complex* centerFrequenciesB, FFTW_complex* centerElecFrequenciesB,
		FFTW_complex* centerFrequenciesProduct, FFTW_complex* centerFrequenciesElecProduct,
		FFTW_complex* sparseProfile, FFTW_complex* sparseShapeProfile, FFTW_complex* sparseElecProfile, FFTW_complex* sparseHbondProfile,
		FFTW_complex* sparseHydrophobicityProfile, FFTW_complex* sparseHydrophobicityTwoProfile, FFTW_complex* sparseSimpleComplementarityProfile,
		int *validOutputMap,
		FFTW_plan freqPlan, FFTW_plan elecFreqPlan,
		sparse3DFFT_plan sparseFreqPlanBackward,
		FFTW_complex *fkB, FFTW_complex *fkBElec, FFTW_complex *fkBHbond, FFTW_complex *fkBHydrophobicity, FFTW_complex *fkBHydrophobicityTwo,
		FFTW_complex *fkBSimpleComplementarity,
		bool breakDownScores,
		FFTW_complex *centerFrequenciesA_01, FFTW_complex *centerFrequenciesA_10, FFTW_complex *centerFrequenciesA_11,
		FFTW_complex *fkB_01, FFTW_complex *fkB_10, FFTW_complex *fkB_11,
		FFTW_complex *sparseProfile_01, FFTW_complex *sparseProfile_10, FFTW_complex *sparseProfile_11,
		FFTW_complex* freqHat, FFTW_plan freqHatPlan,
		FFTW_complex* ourMoreFrequencies, FFTW_complex* ourMoreFrequenciesOut, FFTW_plan moreFreqPlan, sparse3DFFT_plan sparseFreqPlanForward,
		FFTW_complex* elecGridB, FFTW_plan moreElecFreqPlan,
		FFTW_complex* smallElectrostaticsKernel,
		bool smoothSkin, SmoothingFunction* smoothingFunction, TopValues *localTopValues, double functionScaleFactor, PARAMS *pr )
{
	double scaleB = scaleA;
	Matrix rotMat;

	if ( rotateVolume )
	{
		for ( int c = 0; c < numFreq * numFreq * numFreq; c++ )
			ourMoreFrequencies[ c ][ 0 ] = ourMoreFrequencies[ c ][ 1 ] = 0;

		if ( elecScale != 0 )
			for ( int c = 0; c < numFreq * numFreq * numFreq; c++ )
				elecGridB[ c ][ 0 ] = elecGridB[ c ][ 1 ] = 0;
	}

	while ( 1 )
	{
		int r = nextRotation( pr );

		if ( r < 0 ) break;

		printf("# \n# THREAD = %d, ROTATION = %d\n# \n", threadID, r + 1 );
		fflush( stdout );

		if ( rotateVolume )
		{
			rotMat = Matrix( rotations[ r * 9 + 0 ], rotations[ r * 9 + 1 ], rotations[ r * 9 + 2 ], 0.0,
					rotations[ r * 9 + 3 ], rotations[ r * 9 + 4 ], rotations[ r * 9 + 5 ], 0.0,
					rotations[ r * 9 + 6 ], rotations[ r * 9 + 7 ], rotations[ r * 9 + 8 ], 0.0,                                                        0.0,                    0.0,                    0.0, 1.0 );
		}
		else
		{
			// transform B
			if( !transformAndNormalize( xkBOrig, ykBOrig, zkBOrig, rkBOrig, xkB, ykB, zkB, rkB, numCentersB,
						rotations[ r * 9 + 0 ],	rotations[ r * 9 + 1 ], rotations[ r * 9 + 2 ],
						rotations[ r * 9 + 3 ], rotations[ r * 9 + 4 ], rotations[ r * 9 + 5 ],
						rotations[ r * 9 + 6 ], rotations[ r * 9 + 7 ], rotations[ r * 9 + 8 ],
						2, scaleB ) ) return false;
		}


		if ( elecScale != 0 )
		{
			if ( rotateVolume ) rotateElecGrid( ( FFTW_DATA_TYPE * ) elecGridB, numFreq, elecGridBCells, numNonzeroElecGridBCells, rotMat, false, false );

			if ( !getCenterElecFrequencies( numCentersB, xkB, ykB, zkB, rkB, 0, fkBElec, blobbiness, alpha, numFreq, elecRadiusInGrids,
						centerElecFrequenciesB, ( FFTW_DATA_TYPE * ) elecGridB, moreElecFreqPlan, rotateVolume, true ) ) return false;

			if ( !multiplyElecFrequencyMaps( centerElecFrequenciesA, centerElecFrequenciesB, numFreq, numFreq, smallElectrostaticsKernel, centerFrequenciesElecProduct, scaleA) ) return false;

			FFTW_execute( elecFreqPlan );

			if ( rotateVolume ) rotateElecGrid( ( FFTW_DATA_TYPE * ) elecGridB, numFreq, elecGridBCells, numNonzeroElecGridBCells, rotMat, false, true );
		}

		if ( hbondWeight != 0 )
		{
			if ( rotateVolume ) rotateHbondGrid( ourMoreFrequencies, numFreq, hbondGridBCells, numNonzeroHbondGridBCells, rotMat, false, false );

			if ( !getCenterHbondFrequencies( numCentersB, xkB, ykB, zkB, rkB, hbondDistanceCutoff * scaleA * ( numFreq - 1 ), fkBHbond,
						blobbiness, alpha, numFreq,
						centerFrequenciesB, ourMoreFrequencies, ourMoreFrequenciesOut, moreFreqPlan, sparseFreqPlanForward, rotateVolume, true ) ) return false;

			if ( !multiplyFrequencyMaps( centerHbondFrequenciesA, centerFrequenciesB, numFreq, centerFrequenciesProduct ) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseHbondProfile, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseHbondProfile );

			if ( rotateVolume ) rotateHbondGrid( ourMoreFrequencies, numFreq, hbondGridBCells, numNonzeroHbondGridBCells, rotMat, false, true );
		}

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
		{
			if ( rotateVolume ) rotateHydrophobicityGrid( ourMoreFrequencies, numFreq, hydrophobicityGridBCells, numNonzeroHydrophobicityGridBCells, rotMat, false, false );

			if ( !getCenterHydrophobicityFrequencies( numCentersB, xkB, ykB, zkB, rkB, fkBHydrophobicity,
						blobbiness, alpha, numFreq, pr->pri->hydroRadExt, ( bool ) ( hydrophobicityWeight == 0 ),
						centerFrequenciesB, ourMoreFrequencies, ourMoreFrequenciesOut,
						moreFreqPlan, sparseFreqPlanForward, rotateVolume ) ) return false;

			if ( !multiplyFrequencyMaps( centerHydrophobicityFrequenciesA, centerFrequenciesB, numFreq, centerFrequenciesProduct ) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseHydrophobicityProfile, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseHydrophobicityProfile );

			if ( rotateVolume ) rotateHydrophobicityGrid( ourMoreFrequencies, numFreq, hydrophobicityGridBCells, numNonzeroHydrophobicityGridBCells, rotMat, false, true );

			if ( pr->pri->twoWayHydrophobicity )
			{
				if ( rotateVolume ) rotateHydrophobicityGrid( ourMoreFrequencies, numFreq, hydrophobicityTwoGridBCells, numNonzeroHydrophobicityTwoGridBCells, rotMat, false, false );

				if ( !getCenterHydrophobicityFrequencies( numCentersB, xkB, ykB, zkB, rkB, fkBHydrophobicityTwo,
							blobbiness, alpha, numFreq, 0, false,
							centerFrequenciesB, ourMoreFrequencies,
							ourMoreFrequenciesOut,
							moreFreqPlan, sparseFreqPlanForward, rotateVolume ) ) return false;

				if ( !multiplyFrequencyMaps( centerHydrophobicityTwoFrequenciesA, centerFrequenciesB, numFreq, centerFrequenciesProduct ) ) return false;

				if ( sparseFreqPlanBackward == NULL )
				{
					FFTW_execute( freqPlan );
					memcpy( sparseHydrophobicityTwoProfile, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
				}
				else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseHydrophobicityTwoProfile );

				for ( int c = 0; c < numFreq * numFreq * numFreq; c++ )
				{
					sparseHydrophobicityProfile[ c ][ 0 ] += sparseHydrophobicityTwoProfile[ c ][ 0 ];
					sparseHydrophobicityProfile[ c ][ 1 ] += sparseHydrophobicityTwoProfile[ c ][ 1 ];
				}

				if ( rotateVolume ) rotateHydrophobicityGrid( ourMoreFrequencies, numFreq, hydrophobicityTwoGridBCells, numNonzeroHydrophobicityTwoGridBCells, rotMat, false, true );
			}
		}

		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		{
			if ( rotateVolume ) rotateSimpleComplementarityGrid( ourMoreFrequencies, numFreq, simpleComplementarityGridBCells, numNonzeroSimpleComplementarityGridBCells, rotMat, false, false );

			if ( !getCenterSimpleComplementarityFrequencies( numCentersB, xkB, ykB, zkB, rkB, fkBSimpleComplementarity,
						blobbiness, alpha, numFreq, pr->pri->simpleRadExt,
						centerFrequenciesB, ourMoreFrequencies, ourMoreFrequenciesOut,
						moreFreqPlan, sparseFreqPlanForward, rotateVolume ) ) return false;

			if ( !multiplyFrequencyMaps( centerSimpleComplementarityFrequenciesA, centerFrequenciesB, numFreq, centerFrequenciesProduct ) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseSimpleComplementarityProfile, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseSimpleComplementarityProfile );

			if ( rotateVolume ) rotateSimpleComplementarityGrid( ourMoreFrequencies, numFreq, simpleComplementarityGridBCells, numNonzeroSimpleComplementarityGridBCells, rotMat, false, true );
		}

		if ( breakDownScores )
		{
			if ( rotateVolume ) rotateGrid( ourMoreFrequencies, numFreq, gridBCells_01, numNonzeroGridBCells_01, rotMat, false, false );

			// get frequencies of B
			if ( !getCenterFrequencies( numCentersB, xkB, ykB, zkB, rkB, 0, fkB_01, blobbiness, alpha, numFreq, interpFuncExtent, centerFrequenciesB, smoothSkin, smoothingFunction, ourMoreFrequencies, ourMoreFrequenciesOut, moreFreqPlan, sparseFreqPlanForward, rotateVolume, false ) ) return false;

			// get product with gaussian squared and two frequency maps
			if ( !multiplyFrequencyMaps( centerFrequenciesA_01, centerFrequenciesB, numFreq, gridSize*4, blobbiness, smoothSkin, centerFrequenciesProduct, freqHat, freqHatPlan, scaleA) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseProfile_01, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseProfile_01 );

			if ( rotateVolume )
			{
				rotateGrid( ourMoreFrequencies, numFreq, gridBCells_01, numNonzeroGridBCells_01, rotMat, false, true );
				rotateGrid( ourMoreFrequencies, numFreq, gridBCells_10, numNonzeroGridBCells_10, rotMat, false, false );
			}

			// get frequencies of B
			if ( !getCenterFrequencies( numCentersB, xkB, ykB, zkB, rkB, 0, fkB_10, blobbiness, alpha, numFreq, interpFuncExtent, centerFrequenciesB, smoothSkin, smoothingFunction, ourMoreFrequencies, ourMoreFrequenciesOut, moreFreqPlan, sparseFreqPlanForward, rotateVolume, false ) ) return false;

			// get product with gaussian squared and two frequency maps
			if ( !multiplyFrequencyMaps( centerFrequenciesA_10, centerFrequenciesB, numFreq, gridSize * 4, blobbiness, smoothSkin, centerFrequenciesProduct, freqHat, freqHatPlan, scaleA ) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseProfile_10, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseProfile_10 );

			if ( rotateVolume )
			{
				rotateGrid( ourMoreFrequencies, numFreq, gridBCells_10, numNonzeroGridBCells_10, rotMat, false, true );
				rotateGrid( ourMoreFrequencies, numFreq, gridBCells_11, numNonzeroGridBCells_11, rotMat, false, false );
			}

			// get frequencies of B
			if ( !getCenterFrequencies( numCentersB, xkB, ykB, zkB, rkB, 0, fkB_11, blobbiness, alpha, numFreq, interpFuncExtent, centerFrequenciesB, smoothSkin, smoothingFunction, ourMoreFrequencies, ourMoreFrequenciesOut, moreFreqPlan, sparseFreqPlanForward, rotateVolume, false ) ) return false;

			// get product with gaussian squared and two frequency maps
			if ( !multiplyFrequencyMaps( centerFrequenciesA_11, centerFrequenciesB, numFreq, gridSize * 4, blobbiness, smoothSkin, centerFrequenciesProduct, freqHat, freqHatPlan, scaleA) ) return false;

			if ( sparseFreqPlanBackward == NULL )
			{
				FFTW_execute( freqPlan );
				memcpy( sparseProfile_11, sparseProfile, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseProfile_11 );

			for ( int c = 0; c < numFreq * numFreq * numFreq ; c++ )
			{
				sparseProfile_11[ c ][ 0 ] -= ( sparseProfile_10[ c ][ 0 ] + sparseProfile_01[ c ][ 0 ] );
				sparseProfile_11[ c ][ 1 ] -= ( sparseProfile_10[ c ][ 1 ] + sparseProfile_01[ c ][ 1 ] );

				sparseShapeProfile[ c ][ 0 ] = skinSkinWeight * sparseProfile_10[ c ][ 0 ] + coreCoreWeight * sparseProfile_01[ c ][ 0 ]
					+ skinCoreWeight * sparseProfile_11[ c ][ 1 ];
			}

			if ( !localTopValues->updateTopPositions( validOutputMap, sparseShapeProfile, sparseProfile_10, sparseProfile_01, sparseProfile_11,
						sparseElecProfile, elecScale, sparseHbondProfile, hbondWeight,
						sparseHydrophobicityProfile, ( ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) ) ? 1.0 : 0.0,
						( hydroPhobicPhobicWeight * hydroPhilicPhilicWeight > 0 ) ? ( hydroPhobicPhilicWeight / sqrt( hydroPhobicPhobicWeight * hydroPhilicPhilicWeight ) ) : 0.0,
						sparseSimpleComplementarityProfile,
						r, 0, confID, pr ) )
				markBadRotGraphNodes( r, pr );

			if ( rotateVolume ) rotateGrid( ourMoreFrequencies, numFreq, gridBCells_11, numNonzeroGridBCells_11, rotMat, false, true );
		}
		else
		{
			if ( rotateVolume ) rotateGrid( ourMoreFrequencies, numFreq, gridBCells, numNonzeroGridBCells, rotMat, false, false );

			// get frequencies of B
			if ( !getCenterFrequencies( numCentersB, xkB, ykB, zkB, rkB, 0, fkB, blobbiness, alpha, numFreq, interpFuncExtent, centerFrequenciesB, smoothSkin, smoothingFunction, ourMoreFrequencies, ourMoreFrequenciesOut, moreFreqPlan, sparseFreqPlanForward, rotateVolume, false ) ) return false;

			// get product with gaussian squared and two frequency maps
			if ( !multiplyFrequencyMaps( centerFrequenciesA, centerFrequenciesB, numFreq, gridSize * 4, blobbiness, smoothSkin, centerFrequenciesProduct, freqHat, freqHatPlan, scaleA) ) return false;

			if ( sparseFreqPlanBackward == NULL ) FFTW_execute( freqPlan );
			else sparse3DFFT( sparseFreqPlanBackward, centerFrequenciesProduct, sparseProfile );

			if ( !localTopValues->updateTopPositions( validOutputMap,
						sparseProfile, 1.0, ( skinSkinWeight * coreCoreWeight > 0 ) ? ( skinCoreWeight / sqrt( skinSkinWeight * coreCoreWeight ) ) : 0.0,
						sparseElecProfile, elecScale, sparseHbondProfile, hbondWeight,
						sparseHydrophobicityProfile, ( ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) ) ? 1.0 : 0.0,
						( hydroPhobicPhobicWeight * hydroPhilicPhilicWeight > 0 ) ? ( hydroPhobicPhilicWeight / sqrt( hydroPhobicPhobicWeight * hydroPhilicPhilicWeight ) ) : 0.0,
						sparseSimpleComplementarityProfile,
						r, 0, confID, pr ) )
				markBadRotGraphNodes( r, pr );

			if ( rotateVolume ) rotateGrid( ourMoreFrequencies, numFreq, gridBCells, numNonzeroGridBCells, rotMat, false, true );
		}
	}

	return true;
}



static void *startApplyRotationsThread( void *v )
{
	PARAMS *pr = ( PARAMS * ) v;

	applyRotations( pr->threadID, pr->confID, pr->rotations,
			pr->numberOfRotations, pr->numberOfPositions,
			pr->xkBOrig, pr->ykBOrig, pr->zkBOrig, pr->radiiB,
			pr->xkB, pr->ykB, pr->zkB, pr->rkB, pr->numCentersB,
			pr->rotateVolume,
			pr->gridBCells, pr->numNonzeroGridBCells,
			pr->gridBCells_01, pr->numNonzeroGridBCells_01,
			pr->gridBCells_10, pr->numNonzeroGridBCells_10,
			pr->gridBCells_11, pr->numNonzeroGridBCells_11,
			pr->elecGridBCells, pr->numNonzeroElecGridBCells,
			pr->hbondGridBCells, pr->numNonzeroHbondGridBCells,
			pr->hydrophobicityGridBCells, pr->numNonzeroHydrophobicityGridBCells,
			pr->hydrophobicityTwoGridBCells, pr->numNonzeroHydrophobicityTwoGridBCells,
			pr->simpleComplementarityGridBCells, pr->numNonzeroSimpleComplementarityGridBCells,
			pr->gridSize, pr->numFreq, pr->interpFuncExtent,
			pr->alpha, pr->blobbiness,
			pr->skinSkinWeight, pr->coreCoreWeight, pr->skinCoreWeight,
			pr->realSCWeight, pr->imaginarySCWeight,
			pr->elecScale, pr->elecRadiusInGrids,
			pr->hbondWeight, pr->hbondDistanceCutoff,
			pr->hydrophobicityWeight, pr->hydroPhobicPhobicWeight, pr->hydroPhobicPhilicWeight, pr->hydroPhilicPhilicWeight,
			pr->simpleShapeWeight, pr->simpleChargeWeight,
			pr->scaleA,
			pr->centerFrequenciesA, pr->centerElecFrequenciesA, pr->centerHbondFrequenciesA, pr->centerHydrophobicityFrequenciesA,
			pr->centerHydrophobicityTwoFrequenciesA, pr->centerSimpleComplementarityFrequenciesA,
			pr->centerFrequenciesB, pr->centerElecFrequenciesB,
			pr->centerFrequenciesProduct, pr->centerFrequenciesElecProduct,
			pr->sparseProfile, pr->sparseShapeProfile, pr->sparseElecProfile, pr->sparseHbondProfile,
			pr->sparseHydrophobicityProfile, pr->sparseHydrophobicityTwoProfile,
			pr->sparseSimpleComplementarityProfile,
			pr->validOutputMap,
			pr->freqPlan, pr->elecFreqPlan,
			pr->sparseFreqPlanBackward,
			pr->fkB, pr->fkBElec, pr->fkBHbond, pr->fkBHydrophobicity, pr->fkBHydrophobicityTwo, pr->fkBSimpleComplementarity,
			pr->breakDownScores,
			pr->centerFrequenciesA_01, pr->centerFrequenciesA_10, pr->centerFrequenciesA_11,
			pr->fkB_01, pr->fkB_10, pr->fkB_11,
			pr->sparseProfile_01, pr->sparseProfile_10, pr->sparseProfile_11,
			pr->freqHat, pr->freqHatPlan,
			pr->ourMoreFrequencies, pr->ourMoreFrequenciesOut, pr->moreFreqPlan, pr->sparseFreqPlanForward,
			pr->elecGridB, pr->moreElecFreqPlan,
			pr->smallElectrostaticsKernel,
			pr->smoothSkin, pr->smoothingFunction, pr->localTopValues,
			pr->functionScaleFactor, pr );
}

void printInputParamters( PARAMS_IN *pr, FILE* fp )
	// print the values in PARAMS_IN to a given file
	// print description of columns in a line describing a docking result
{

	fprintf( fp, (char *)"# INPUT PARAMETERS:\n#\n" );
	//	fprintf( fp, (char *)"# \t performDocking = %d\n", pr->performDocking );
	fprintf( fp, (char *)"# \t numThreads = %d\n", pr->numThreads );
	fprintf( fp, (char *)"# \t breakDownScores = %d\n", pr->breakDownScores );
	fprintf( fp, (char *)"# \t numberOfPositions = %d\n", pr->numberOfPositions );
	fprintf( fp, (char *)"# \t gridSize = %d\n", pr->gridSize );
	fprintf( fp, (char *)"# \t gridSizeSpecified = %d\n", pr->gridSizeSpecified );
	fprintf( fp, (char *)"# \t enforceExactGridSpacing = %d\n", pr->enforceExactGridSpacing );
	fprintf( fp, (char *)"# \t gridSpacingSpecified = %d\n", pr->gridSpacingSpecified );
	fprintf( fp, (char *)"# \t numFreq = %d\n", pr->numFreq );
	//	fprintf( fp, (char *)"# \t numFreqSpecified = %d\n", pr->numFreqSpecified );
	fprintf( fp, (char *)"# \t smoothSkin = %d\n", pr->smoothSkin );
	fprintf( fp, (char *)"# \t singleLayerLigandSkin = %d\n", pr->singleLayerLigandSkin );
	fprintf( fp, (char *)"# \t pseudoAtomRadius = %lf\n", pr->pseudoAtomRadius );
	fprintf( fp, (char *)"# \t rotateVolume = %d\n", pr->rotateVolume );
	fprintf( fp, (char *)"# \t dockVolume = %d\n", pr->dockVolume );
	fprintf( fp, (char *)"# \t numberOfRotations = %d\n", pr->numberOfRotations );
	fprintf( fp, (char *)"# \t pruneAngle = %lf\n", pr->pruneAngle );
	fprintf( fp, (char *)"# \t curvatureWeightedStaticMol = %d\n", pr->curvatureWeightedStaticMol );
	fprintf( fp, (char *)"# \t curvatureWeightedMovingMol = %d\n", pr->curvatureWeightedMovingMol );
	fprintf( fp, (char *)"# \t curvatureWeightingRadius = %lf\n", pr->curvatureWeightingRadius );
	fprintf( fp, (char *)"# \t spreadReceptorSkin = %d\n", pr->spreadReceptorSkin );
	fprintf( fp, (char *)"# \t randomRotate = %d\n", pr->randomRotate );
	fprintf( fp, (char *)"# \t initialRotation = < %f %f %f %f %f %f %f %f %f >\n",
			pr->initRot.get( 0, 0 ), pr->initRot.get( 0, 1 ), pr->initRot.get( 0, 2 ),
			pr->initRot.get( 1, 0 ), pr->initRot.get( 1, 1 ), pr->initRot.get( 1, 2 ),
			pr->initRot.get( 2, 0 ), pr->initRot.get( 2, 1 ), pr->initRot.get( 2, 2 ) );
	fprintf( fp, (char *)"# \t useSparseFFT = %d\n", pr->useSparseFFT );
	fprintf( fp, (char *)"# \t narrowBand = %d\n", pr->narrowBand );
	//fprintf( fp, (char *)"# \t  = %d\n", pr-> );

	fprintf( fp, (char *)"# \t gridSpacing = %lf\n", pr->gridSpacing );
	fprintf( fp, (char *)"# \t interpFuncExtentInAngstroms = %lf\n", pr->interpFuncExtentInAngstroms );
	//	fprintf( fp, (char *)"# \t distanceCutoff = %lf\n", pr->distanceCutoff );
	fprintf( fp, (char *)"# \t alpha = %lf\n", pr->alpha );
	fprintf( fp, (char *)"# \t blobbiness = %lf\n", pr->blobbiness );
	fprintf( fp, (char *)"# \t blobbinessSpecified = %d\n", pr->blobbinessSpecified );
	fprintf( fp, (char *)"# \t skinSkinWeight = %lf\n", pr->skinSkinWeight );
	fprintf( fp, (char *)"# \t coreCoreWeight = %lf\n", pr->coreCoreWeight );
	fprintf( fp, (char *)"# \t skinCoreWeight = %lf\n", pr->skinCoreWeight );
	//	fprintf( fp, (char *)"# \t realSCWeight = %lf\n", pr->realSCWeight );
	//	fprintf( fp, (char *)"# \t imaginarySCWeight = %lf\n", pr->imaginarySCWeight );
	fprintf( fp, (char *)"# \t elecScale = %lf\n", pr->elecScale );
	fprintf( fp, (char *)"# \t elecRadiusInGrids = %lf\n", pr->elecRadiusInGrids );
	fprintf( fp, (char *)"# \t elecKernelVoidRad = %lf\n", pr->elecKernelVoidRad );
	fprintf( fp, (char *)"# \t elecKernelDistLow = %lf\n", pr->elecKernelDistLow );
	fprintf( fp, (char *)"# \t elecKernelDistHigh = %lf\n", pr->elecKernelDistHigh );
	fprintf( fp, (char *)"# \t elecKernelValLow = %lf\n", pr->elecKernelValLow );
	fprintf( fp, (char *)"# \t elecKernelValHigh = %lf\n", pr->elecKernelValHigh );
	fprintf( fp, (char *)"# \t hbondWeight = %lf\n", pr->hbondWeight );
	fprintf( fp, (char *)"# \t hbondDistanceCutoff = %lf\n", pr->hbondDistanceCutoff );
	fprintf( fp, (char *)"# \t hydrophobicityWeight = %lf\n", pr->hydrophobicityWeight );
	fprintf( fp, (char *)"# \t hydrophobicityProductWeight = %lf\n", pr->hydrophobicityProductWeight );
	fprintf( fp, (char *)"# \t hydroRatioTolerance = %lf\n", pr->hydroRatioTolerance );
	fprintf( fp, (char *)"# \t twoWayHydrophobicity = %d\n", pr->twoWayHydrophobicity );
	fprintf( fp, (char *)"# \t hydroPhobicPhobicWeight = %lf\n", pr->hydroPhobicPhobicWeight );
	fprintf( fp, (char *)"# \t hydroPhobicPhilicWeight = %lf\n", pr->hydroPhobicPhilicWeight );
	fprintf( fp, (char *)"# \t hydroPhilicPhilicWeight = %lf\n", pr->hydroPhilicPhilicWeight );
	fprintf( fp, (char *)"# \t hydroRadExt = %lf\n", pr->hydroRadExt );
	fprintf( fp, (char *)"# \t useInterfacePropensity = %d\n", pr->useInterfacePropensity );
	fprintf( fp, (char *)"# \t perResidueHydrophobicity = %d\n", pr->perResidueHydrophobicity );
	fprintf( fp, (char *)"# \t numTopHydrophobicResidues = %d\n", pr->numTopHydrophobicResidues );
	fprintf( fp, (char *)"# \t staticMolHydroDistCutoff = %lf\n", pr->staticMolHydroDistCutoff );
	fprintf( fp, (char *)"# \t simpleShapeWeight = %lf\n", pr->simpleShapeWeight );
	fprintf( fp, (char *)"# \t simpleChargeWeight = %lf\n", pr->simpleChargeWeight );
	fprintf( fp, (char *)"# \t simpleRadExt = %lf\n", pr->simpleRadExt );
	fprintf( fp, (char *)"# \t scoreScaleUpFactor = %lf\n", pr->scoreScaleUpFactor );
	fprintf( fp, (char *)"# \t bandwidth = %lf\n", pr->bandwidth );
	fprintf( fp, (char *)"# \t gradFactor = %lf\n", pr->gradFactor );
	fprintf( fp, (char *)"# \t applyVdWFilter = %d\n", pr->applyVdWFilter );
	fprintf( fp, (char *)"# \t vdWCutoffLow = %lf\n", pr->vdWCutoffLow );
	fprintf( fp, (char *)"# \t vdWCutoffHigh = %lf\n", pr->vdWCutoffHigh );
	fprintf( fp, (char *)"# \t vdWWellWidth = %lf\n", pr->vdWWellWidth );
	fprintf( fp, (char *)"# \t vdWTolerance = %lf\n", pr->vdWTolerance );
	fprintf( fp, (char *)"# \t compQuadVdW = %lf\n", pr->compQuadVdW );
	fprintf( fp, (char *)"# \t vdWGridSize = %d\n", pr->vdWGridSize );
	fprintf( fp, (char *)"# \t surfaceBasedVdW = %d\n", pr->surfaceBasedVdW );
	fprintf( fp, (char *)"# \t applyClashFilter = %d\n", pr->applyClashFilter );
	fprintf( fp, (char *)"# \t eqmDistFrac = %lf\n", pr->eqmDistFrac );
	fprintf( fp, (char *)"# \t clashTolerance = %d\n", pr->clashTolerance );
	fprintf( fp, (char *)"# \t filterDepth = %d\n", pr->filterDepth );
	fprintf( fp, (char *)"# \t clashWeight = %lf\n", pr->clashWeight );
	fprintf( fp, (char *)"# \t clusterTransRad = %lf\n", pr->clusterTransRad );
	fprintf( fp, (char *)"# \t clusterRotRad = %lf\n", pr->clusterRotRad );
	fprintf( fp, (char *)"# \t applyPseudoGsolFilter = %d\n", pr->applyPseudoGsolFilter );
	fprintf( fp, (char *)"# \t pseudoGsolCutoff = %lf\n", pr->pseudoGsolCutoff );
	fprintf( fp, (char *)"# \t pseudoGsolWeight = %lf\n", pr->pseudoGsolWeight );
	fprintf( fp, (char *)"# \t pseudoGsolFilterLowestRank = %d\n", pr->pseudoGsolFilterLowestRank );
	fprintf( fp, (char *)"# \t applyDispersionFilter = %d\n", pr->applyDispersionFilter );
	fprintf( fp, (char *)"# \t dispersionCutoff = %lf\n", pr->dispersionCutoff );
	fprintf( fp, (char *)"# \t dispersionWeight = %lf\n", pr->dispersionWeight );
	fprintf( fp, (char *)"# \t filterScaleDownFactor = %lf\n", pr->filterScaleDownFactor );
	// fprintf( fp, (char *)"# \t  = %lf\n", pr-> );

	fprintf( fp, (char *)"# \t peaksPerRotation = %d\n", pr->peaksPerRotation );
	fprintf( fp, (char *)"# \t spectrum = %s\n", pr->spectrum );
	// fprintf( fp, (char *)"# \t  = %d\n", pr-> );

	fprintf( fp, (char *)"# \t staticMoleculeF2d = %s\n", pr->staticMoleculeF2d );
	fprintf( fp, (char *)"# \t movingMoleculeF2d = %s\n", pr->movingMoleculeF2d );

	int numStaticMolSkinAtoms = 0, numStaticMolCoreAtoms = 0;

	for ( int i = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' ) numStaticMolCoreAtoms++;
		else if ( pr->typeA[ i ] == 'E' ) numStaticMolSkinAtoms++;

	int numMovingMolSkinAtoms = 0, numMovingMolCoreAtoms = 0;

	for ( int i = 0; i < pr->numCentersB; i++ )
		if ( pr->typeB[ i ] == 'I' ) numMovingMolCoreAtoms++;
		else if ( pr->typeB[ i ] == 'E' ) numMovingMolSkinAtoms++;

	fprintf( fp, (char *)"# \t numStaticMolSkinAtoms = %d\n", numStaticMolSkinAtoms );
	fprintf( fp, (char *)"# \t numStaticMolCoreAtoms = %d\n", numStaticMolCoreAtoms );
	fprintf( fp, (char *)"# \t numMovingMolSkinAtoms = %d\n", numMovingMolSkinAtoms );
	fprintf( fp, (char *)"# \t numMovingMolCoreAtoms = %d\n", numMovingMolCoreAtoms );

	double sumStaticMolPosCharges = 0, sumStaticMolNegCharges = 0;

	for ( int i = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' )
		{
			if ( pr->chargesA[ i ] > 0 ) sumStaticMolPosCharges += pr->chargesA[ i ];
			else sumStaticMolNegCharges += pr->chargesA[ i ];
		}

	double sumMovingMolPosCharges = 0, sumMovingMolNegCharges = 0;

	for ( int i = 0; i < pr->numCentersB; i++ )
	{
		if ( pr->chargesB[ i ] > 0 ) sumMovingMolPosCharges += pr->chargesB[ i ];
		else sumMovingMolNegCharges += pr->chargesB[ i ];
	}

	fprintf( fp, (char *)"# \t sumStaticMolCharges = %lf ( %lf, %lf )\n", sumStaticMolPosCharges + sumStaticMolNegCharges,
			sumStaticMolPosCharges, sumStaticMolNegCharges );
	fprintf( fp, (char *)"# \t sumMovingMolCharges = %lf ( %lf, %lf )\n", sumMovingMolPosCharges + sumMovingMolNegCharges,
			sumMovingMolPosCharges, sumMovingMolNegCharges );

	fprintf( fp, (char *)"# \t staticMoleculeSCReRaw = %s\n", pr->staticMoleculeSCReRaw );
	fprintf( fp, (char *)"# \t staticMoleculeSCImRaw = %s\n", pr->staticMoleculeSCImRaw );
	fprintf( fp, (char *)"# \t staticMoleculeElecReRaw = %s\n", pr->staticMoleculeElecReRaw );
	fprintf( fp, (char *)"# \t movingMoleculeSCReRaw = %s\n", pr->movingMoleculeSCReRaw );
	fprintf( fp, (char *)"# \t movingMoleculeSCImRaw = %s\n", pr->movingMoleculeSCImRaw );
	fprintf( fp, (char *)"# \t movingMoleculeElecReRaw = %s\n", pr->movingMoleculeElecReRaw );
	fprintf( fp, (char *)"# \t outputFilename = %s\n", pr->outputFilename );

	fprintf( fp, (char *)"#\n" );
	fprintf( fp, (char *)"#\n" );
	if (pr->breakDownScores)
		fprintf( fp, (char *)"# OUTPUT FORMAT: 32\n" );
	else
		fprintf( fp, (char *)"# OUTPUT FORMAT: 29\n" );

	fprintf( fp, (char *)"#\t COLNAME rank int rank of this docking result according to the F2Dock scoring function\n");
	fprintf( fp, (char *)"#\t COLNAME score float overall F2Dock score, the higher the better\n");
	fprintf( fp, (char *)"#\t COLNAME shape float shape complementarity score\n");
	if ( pr->breakDownScores )
	{
		fprintf( fp, (char *)"#\t COLNAME ssr float \n");
		fprintf( fp, (char *)"#\t COLNAME ccr float \n");
		fprintf( fp, (char *)"#\t COLNAME scr float \n");
		fprintf( fp, (char *)"#\t COLNAME ssi float \n");
		fprintf( fp, (char *)"#\t COLNAME cci float \n");
		fprintf( fp, (char *)"#\t COLNAME sci float \n");
	}
	else
	{
		//	   fprintf( fp, (char *)"#\t COLNAME realshape float real part of shape complementarity score\n");
		//	   fprintf( fp, (char *)"#\t COLNAME unrealshape float imaginary part of shape complementarity score\n");
		fprintf( fp, (char *)"#\t COLNAME ssr float skin-skin overlap score\n");
		fprintf( fp, (char *)"#\t COLNAME ccr float core-core overlap score\n");
		fprintf( fp, (char *)"#\t COLNAME scr float skin-core overlap score\n");
	}
	fprintf( fp, (char *)"#\t COLNAME elec float \n");
	fprintf( fp, (char *)"#\t COLNAME hbond float \n");
	fprintf( fp, (char *)"#\t COLNAME hydrophobicity float \n");
	fprintf( fp, (char *)"#\t COLNAME smplcomp float \n");
	fprintf( fp, (char *)"#\t COLNAME vdw float \n");
	fprintf( fp, (char *)"#\t COLNAME clashes int \n");
	fprintf( fp, (char *)"#\t COLNAME pgsol float \n");
	fprintf( fp, (char *)"#\t COLNAME pgsolh float \n");
	fprintf( fp, (char *)"#\t COLNAME deldispe float \n");
	fprintf( fp, (char *)"#\t COLNAME mat1 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat2 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat3 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat4 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat5 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat6 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat7 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat8 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat9 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat10 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat11 float \n");
	fprintf( fp, (char *)"#\t COLNAME mat12 float \n");
	fprintf( fp, (char *)"#\t COLNAME conf int \n");
	fprintf( fp, (char *)"#\t COLNAME rmsd float \n");

	fprintf( fp, (char *)"#\n" );
	fprintf( fp, (char *)"#\n" );
}


int getFromGlobalIn( TopValues *globalIn, ValuePosition3D &sol, bool setVars, int maxCount )
{
	static int count = 0;
	static int cmax = globalIn->getCurrentNumberOfPositions( );

	int retVal = 0;

	if ( setVars )
	{
		count = 0;
		cmax = globalIn->getCurrentNumberOfPositions( );

		if ( maxCount < cmax ) cmax = maxCount;

		pthread_mutex_init( &globalInLock, NULL );
	}
	else
	{
		pthread_mutex_lock( &globalInLock );

		if ( count == cmax ) retVal = 0;
		else
		{
			if ( globalIn->extractMin( sol ) )
			{
				retVal = globalIn->getCurrentNumberOfPositions( ) + 1;
				count++;
			}
		}

		pthread_mutex_unlock( &globalInLock );
	}

	return retVal;
}


bool insertIntoGlobalOut( TopValues *globalOut, ValuePosition3D &sol, bool setVars )
{
	bool retVal = true;

	if ( setVars )
	{
		pthread_mutex_init( &globalOutLock, NULL );
	}
	else
	{
		pthread_mutex_lock( &globalOutLock );

		retVal = globalOut->updateTopValues( sol );

		pthread_mutex_unlock( &globalOutLock );
	}

	return retVal;
}


void applyFilters( FILTER_PARAMS *pr )
{
	ValuePosition3D sol;
	int retVal;

	while ( ( retVal =  getFromGlobalIn( pr->TopValuesIn, sol, false, 0 ) ) > 0 )
	{
		printf("# \n# THREAD = %d, SOL = %d\n# \n", pr->threadID, retVal );
		fflush( stdout );

		Matrix transM;
		double transD[ 16 ];

		retrieveTransformation( sol, pr, transM, transD );

		bool filtered = false;
		bool fFiltered = false;

		if ( pr->pri->applyForbiddenVolumeFilter )
		{

			if(pr->pri->forbiddenVolumeFileType%2==0){
				double intVal;
				int nclashes = 0, nclose = 0;

				pr->pri->fFilter->computeInteractions( transM, &nclashes, &nclose, &intVal );

				sol.m_nClashes = nclashes;

				if ( nclashes > pr->pri->forbiddenVolClashTolerance ) fFiltered = true;
			}

			else {

				Vector oldPos( pr->pri->centroidxB,pr->pri->centroidyB, pr->pri->centroidzB,  1.0 );
				Vector newPos = transM * oldPos;
				double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];


				if ( xM< pr->pri->forbiddenBBoxXMax && yM < pr->pri->forbiddenBBoxYMax && zM < pr->pri->forbiddenBBoxZMax && 
						xM > pr->pri->forbiddenBBoxXMin && yM > pr->pri->forbiddenBBoxYMin && zM > pr->pri->forbiddenBBoxZMin)
					fFiltered = true;
			}
		}

		if ( pr->pri->applyClashFilter )
		{
			double intVal;
			int nclashes = 0, nclose = 0;

			pr->pri->cFilter->computeInteractions( transM, &nclashes, &nclose, &intVal );

			sol.m_nClashes = nclashes;

			if ( nclashes > pr->pri->clashTolerance ) filtered = true;
		}

		if ( pr->pri->applyVdWFilter )
		{
			double vdwp = 0, vdwpQ = 0;
			int nclashesT = 0, ncloseT = 0;

			//computeSingleVdWPotential( pr->pri->vdWParams, transM, &vdwp, &vdwpQ, &nclashesT, &ncloseT );
			pr->pri->ljFilter->computePotential( pr->threadID - 1, transD, &vdwp );

			sol.m_vdWPotential = vdwp;

			//           if ( ( sol.m_ConformationIndex < 3 ) && ( !pr->pri->applyClashFilter || ( sol.m_nClashes < ( pr->pri->clashTolerance >> 1 ) ) ) )
			//              {
			//                if ( vdwp > pr->pri->vdWCutoff ) filtered = true;
			//              }

			if ( !pr->pri->applyClashFilter || ( sol.m_nClashes < ( pr->pri->clashTolerance >> 1 ) ) ) 
			{
				if ( vdwp > pr->pri->vdWCutoffLow ) filtered = true;
			}
			else   
			{
				if ( vdwp > pr->pri->vdWCutoffHigh ) filtered = true;
			}
		}

		if ( pr->pri->applyPseudoGsolFilter )
		{
			double pseudoGsol, pGsolHStaticPos, pGsolHStaticNeg, pGsolHMovingPos, pGsolHMovingNeg;

			pr->pri->pGsol->getPseudoGsol( pr->threadID - 1, transD, &pseudoGsol,
					&pGsolHStaticPos, &pGsolHStaticNeg, &pGsolHMovingPos, &pGsolHMovingNeg );

			if ( pGsolHStaticPos + pGsolHMovingPos > 0 )
			{
				sol.m_pGsol = - ( pGsolHStaticNeg + pGsolHMovingNeg ) / ( pGsolHStaticPos + pGsolHMovingPos );
				sol.m_pGsolH = - ( pGsolHStaticNeg + pGsolHMovingNeg ) * sol.m_pGsol;
			}
			else
			{
				sol.m_pGsol = 0;//-1000.0;
				sol.m_pGsolH = -1000.0;
			}

			pr->pgsolSum += sol.m_pGsol;

			//           filtered = updatePseudoGsolCutoff( pr->pri, pseudoGsol );
		}

#ifdef RERANK_DEBUG
		if ( pr->pri->applyDispersionFilter )
			sol.m_delDispE = pr->pri->dispF->getDelDispE( pr->threadID - 1, transD );
#endif

		if ( sol.m_simpComp > 0 ) filtered = true;

		if ( filtered ) sol.m_Value *= pr->pri->filterScaleDownFactor;
		if(fFiltered) {
			std::cout<<"Filtering solution that enters the forbidden volume.\n";
			sol.m_Value = 0;
		}

		//       double delDispE = sol.m_delDispE;
		//
		//       if ( delDispE < - pr->pri->dispersionEnergyLimit ) delDispE = - pr->pri->dispersionEnergyLimit;
		//       else if ( delDispE > pr->pri->dispersionEnergyLimit ) delDispE = pr->pri->dispersionEnergyLimit;

		sol.m_Value -= ( ( pr->pri->clashWeight * sol.m_nClashes + pr->pri->pseudoGsolWeight * sol.m_pGsolH
					/*+ pr->pri->dispersionWeight * delDispE*/ ) * pr->functionScaleFactor );

		insertIntoGlobalOut( pr->TopValuesOut, sol, false );
	}
}


static void *startApplyFiltersThread( void *v )
{
	FILTER_PARAMS *pr = ( FILTER_PARAMS * ) v;

	applyFilters( pr );
}


void filterPoses( PARAMS *params, TopValues *globalIn, TopValues *globalOut,
		int numFreq, float scale, float *translate_A, float *translate_B, float *rotations,
		double functionScaleFactor, Matrix randRot )
{
	PARAMS_IN *pr = ( PARAMS_IN * ) params->pri;
	int numThreads = pr->numThreads;
	FILTER_PARAMS prT[ numThreads ];
	pthread_t p[ numThreads ];

	ValuePosition3D sol;

	getFromGlobalIn( globalIn, sol, true, pr->numberOfPositions );

	insertIntoGlobalOut( globalOut, sol, true );

	for ( int i = 0; i < numThreads; i++ )
	{
		prT[ i ].threadID = i + 1;
		prT[ i ].TopValuesIn = globalIn;
		prT[ i ].TopValuesOut = globalOut;
		prT[ i ].rotations = rotations;
		prT[ i ].numFreq = numFreq;
		prT[ i ].translate_A = translate_A;
		prT[ i ].scaleB = scale;
		prT[ i ].translate_B = translate_B;
		prT[ i ].functionScaleFactor = functionScaleFactor;
		prT[ i ].randRot = randRot;
		prT[ i ].pri = pr;
		prT[ i ].pgsolSum = 0;

		pthread_create( &p[ i ], NULL, startApplyFiltersThread, ( void * ) &prT[ i ] );
	}

	for ( int i = 0; i < numThreads; i++ )
		pthread_join( p[ i ], NULL );

	double pgsolSum = 0;

	for ( int i = 0; i < numThreads; i++ )
		pgsolSum += prT[ i ].pgsolSum;

	int n = globalIn->getCurrentNumberOfPositions( );

	while ( n-- ) globalIn->extractMin( sol );

	n = globalOut->getCurrentNumberOfPositions( );
	int m = n - pr->pseudoGsolFilterLowestRank;

	double pgsolAvg = pgsolSum / n;
	double factor = 0.55;

	params->pGsolAvg = pgsolAvg;

	while ( n-- )
	{
		globalOut->extractMin( sol );

		//       if ( ( n >= m ) && pr->applyPseudoGsolFilter && ( sol.m_pGsol < factor * pgsolAvg ) )
		//          {
		//            double pseudoGsolPenalty = ( ( ( factor * pgsolAvg ) - sol.m_pGsol ) / ( factor * pgsolAvg ) ) * 0.5 * ( - sol.m_Value );
		//            sol.m_Value += pseudoGsolPenalty;
		//          }

		globalIn->updateTopValues( sol );
	}
}


void rerankPoses( PARAMS *pr, TopValues *globalTV )
{
	int n = globalTV->getCurrentNumberOfPositions( );
	int nr = pr->pri->numRerank;

	if ( nr > n ) nr = n;

	ValuePosition3D *sols = new ValuePosition3D[ n ];
	ValuePosition3D *sols2 = new ValuePosition3D[ nr ];
	ValuePosition3D sol;

	clashFilter *cFilter = NULL, *cFilterL3 = NULL, *cFilterH3 = NULL, *cFilterH2 = NULL, *cFilterL1H1 = NULL;

	if ( ( pr->pri->complexType == 'A' ) && pr->pri->applyAntibodyFilter )
	{
		initAntibodyClashFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, true, false, 3, 3, &cFilterL3 );
		initAntibodyClashFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, false, true, 3, 3, &cFilterH3 );
		initAntibodyClashFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, false, true, 2, 2, &cFilterH2 );
		initAntibodyClashFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, true, true, 1, 1, &cFilterL1H1 );
	}
	else if ( ( pr->pri->complexType == 'E' ) && pr->pri->applyEnzymeFilter )
		initEnzymeClashFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, &cFilter );

	resContFilter *rcFilter = NULL;

	if ( pr->pri->applyResidueContactFilter ) initResContFilter( pr->pri->staticMoleculePQR, pr->pri->movingMoleculePQR, NULL, &rcFilter );

	int m = n;

	TopValues *tempTV = new TopValues( nr, pr->pri->numFreq );

	while ( m-- )
	{
		globalTV->extractMin( sol );
		sol.m_origRank = m + 1;
		sol.m_origScore = sol.m_Value;

		if ( m >= nr )
		{
			sol.m_rerankerScore = 0;
			sols[ m ] = sol;
		}
		else
		{
			printf("# \n# SOL = %d\n# \n", m + 1 );
			fflush( stdout );

			Matrix transM;
			double transD[ 16 ];

			retrieveTransformation( sol, pr, transM, transD );

			double origScore = sol.m_origScore / pr->functionScaleFactor;

			sol.m_rerankerScore = pr->pri->rerankerPseudoGsolWeight * sol.m_pGsolH;

			double lowDiv = 0.0;
			int lowDivID = 0;

			if ( pr->pri->applyPseudoGsolFilter )
			{
				double pgsolFactor = 0.2;

				if ( pr->pri->complexType == 'A' ) pgsolFactor = 0.55;
				else if ( pr->pri->complexType == 'E' ) pgsolFactor = 0.9;

				double pgsolAvgF = pgsolFactor * pr->pGsolAvg;

				if ( sol.m_pGsol < pgsolAvgF )
				{
					double pseudoGsolPenalty = ( ( pgsolAvgF - sol.m_pGsol ) / pgsolAvgF ) * 0.5 * ( - origScore );
					origScore += pseudoGsolPenalty;
					lowDivID = 1;
				}             

				if ( pr->pri->complexType == 'A' )
				{
					if ( ( sol.m_pGsol < 1.0 ) || ( sol.m_pGsol > 8.0 ) || ( sol.m_pGsolH < 450 ) )
					{ 
						lowDiv += 1.0, lowDivID = 2;
					}
				}
				else if ( pr->pri->complexType == 'E' )
				{
					if ( ( sol.m_pGsol < 1.0 ) || ( sol.m_pGsol > 7.0 ) || ( sol.m_pGsolH < 450 ) )
					{ 
						lowDiv += 1.0, lowDivID = 2;
					}
				}
				else
				{
					if ( ( sol.m_pGsol < 0.4 ) || ( sol.m_pGsol > 7.0 ) || ( sol.m_pGsolH < 90 ) )
					{ 
						lowDiv += 1.0, lowDivID = 2;
					}
				}                  
			}

			if ( sol.m_clusterPenalty > 2.0 ) lowDiv += 2.0, lowDivID = 3;
			else if ( sol.m_clusterPenalty > 1.0 ) lowDiv += 1.0, lowDivID = 4;

			if ( sol.m_simpComp > 0 ) lowDiv += 1.0, lowDivID = 5;

			if ( sol.m_nClashes > pr->pri->clashTolerance ) lowDiv += 2.0, lowDivID = 6;

			if ( ( ( pr->pri->complexType == 'E' ) && ( sol.m_elecValue > 20.0 * pr->functionScaleFactor ) ) 
					|| ( ( pr->pri->complexType != 'E' ) && ( sol.m_elecValue > 2.0 * pr->functionScaleFactor ) )  )
				lowDiv += 2.0, lowDivID = 7;

			if ( ( pr->pri->complexType == 'A' ) && pr->pri->applyAntibodyFilter )
			{
				double intVal;
				int nInterfaceL3 = 0, nInterfaceH3 = 0, nInterfaceH2 = 0, nInterfaceL1H1 = 0;
				int nCloseL3 = 0, nCloseH3 = 0, nCloseH2 = 0, nCloseL1H1 = 0;

				if ( cFilterL3 != NULL ) cFilterL3->computeInteractions( transM, &nInterfaceL3, &nCloseL3, &intVal );
				if ( cFilterH3 != NULL ) cFilterH3->computeInteractions( transM, &nInterfaceH3, &nCloseH3, &intVal );
				if ( cFilterH2 != NULL ) cFilterH2->computeInteractions( transM, &nInterfaceH2, &nCloseH2, &intVal );
				if ( cFilterL1H1 != NULL ) cFilterL1H1->computeInteractions( transM, &nInterfaceL1H1, &nCloseL1H1, &intVal );

				if ( ( ( cFilterL3 != NULL ) && ( ( nInterfaceL3 < 150 ) || ( nCloseL3 < 5 ) ) ) 
						|| ( ( cFilterH3 != NULL ) && ( ( nInterfaceH3 < 150 ) || ( nCloseH3 < 10 ) ) )
						|| ( ( cFilterH2 != NULL ) && ( nInterfaceH2 <  50 ) ) 
						|| ( ( cFilterL1H1 != NULL ) && ( nInterfaceL1H1 < 100 ) ) 
						|| ( ( cFilterL3 != NULL ) && ( cFilterH3 != NULL ) && ( ( nInterfaceL3 + nInterfaceH3 < 500 ) || ( nCloseL3 + nCloseH3 < 50 ) ) ) 
						|| ( ( cFilterL3 != NULL ) && ( cFilterH3 != NULL ) && ( cFilterL1H1 != NULL ) && ( ( nInterfaceL3 + nInterfaceH3 + nInterfaceL1H1 < 1200 ) || ( nCloseL3 + nCloseH3 + nCloseL1H1 < 90 ) ) ) 
						|| ( ( cFilterL3 != NULL ) && ( cFilterH3 != NULL ) && ( cFilterH2 != NULL ) && ( cFilterL1H1 != NULL ) && ( nInterfaceL3 + nInterfaceH3 + nInterfaceH2 + nInterfaceL1H1 < 1500 ) ) )
					lowDiv += 2.0, lowDivID = 8;
				else sol.m_rerankerScore += 10.0 * ( nInterfaceL3 + nInterfaceH3 + nInterfaceL1H1 );    
			}

			if ( ( pr->pri->complexType == 'E' ) && pr->pri->applyEnzymeFilter )
			{
				double intVal;
				int nInterface = 0, nClose = 0;

				cFilter->computeInteractions( transM, &nInterface, &nClose, &intVal );

				if ( nInterface < 40 ) lowDiv += 1.0, lowDivID = 9;
			}

			if ( pr->pri->applyResidueContactFilter )
			{
				double resContValPos = 0, resContValNeg = 0;

				rcFilter->computeInteractions( transM, &resContValPos, &resContValNeg );

				if ( ( pr->pri->complexType == 'A' ) && ( ( resContValPos < 25.0 ) || ( resContValNeg > 15.0 ) || ( resContValPos - resContValNeg < 10 )
							|| ( resContValPos < 3 * resContValNeg ) ) )
					lowDiv += 2.0, lowDivID = 10;

				if ( ( pr->pri->complexType == 'E' ) && ( ( resContValPos < 45.0 ) || ( resContValNeg > 25.0 ) || ( resContValPos - resContValNeg < 30.0 )
							|| ( resContValPos < 3 * resContValNeg ) ) )
					lowDiv += 2.0, lowDivID = 10;

				if ( ( pr->pri->complexType == 'G' ) && ( ( resContValPos < 15.0 ) || ( resContValNeg > 25.0 ) || ( resContValPos - resContValNeg < 15.0 )
							|| ( resContValPos < 3 * resContValNeg ) ) )
					lowDiv += 2.0, lowDivID = 10;
			}

			if ( ( lowDiv == 0 ) && pr->pri->applyDispersionFilter )
			{
				sol.m_delDispE = pr->pri->dispF->getDelDispE( 0, transD );

				double delDispE = sol.m_delDispE;

				//                if ( delDispE < - pr->pri->dispersionEnergyLimit ) delDispE = - pr->pri->dispersionEnergyLimit;
				//                else if ( delDispE > pr->pri->dispersionEnergyLimit ) delDispE = pr->pri->dispersionEnergyLimit;
				//
				//                if ( sol.m_origRank < pr->pri->rerankerMinF2DockRank ) sol.m_rerankerScore += pr->pri->rerankerDispersionWeightLow * delDispE;
				//                else sol.m_rerankerScore += pr->pri->rerankerDispersionWeightHigh * delDispE;
				//
				//                if ( ( pr->pri->complexType == 'A' ) && ( delDispE > 100.0 ) )
				//                  lowDiv += 2.0, lowDivID = 11;
			}
			else sol.m_delDispE = 0;

			if ( pr->pri->applyVdWFilter )
			{
				if ( pr->pri->complexType == 'A' )
				{
					if ( sol.m_vdWPotential > 0.0 ) lowDiv += 1.0, lowDivID = 12;
				}
				else if ( pr->pri->complexType == 'E' )
				{
					if ( ( ( sol.m_nClashes < 4 ) && ( sol.m_vdWPotential > 0.0 ) ) 
							|| ( ( sol.m_nClashes >= 4 ) && ( sol.m_vdWPotential > 20.0 ) ) )
						lowDiv += 1.0, lowDivID = 12;
				}
				else
				{
					if ( ( ( sol.m_nClashes < 4 ) && ( sol.m_vdWPotential > 0.0 ) ) 
							|| ( ( sol.m_nClashes >= 4 ) && ( sol.m_vdWPotential > 5.0 ) ) )
						lowDiv += 1.0, lowDivID = 12;
				}                                 
			}

#ifdef LIBMOL_FOUND
			if( pr->pri->applyHbondFilter )
			{
				double deltaHbondEnergy;
				pr->pri->hbFilter->getEnergy(transD,&deltaHbondEnergy);

				/* Must add code to update the score here */
			}
#else
			printf("LibMol was not found. Could not apply hygrogen bond filter\n");
#endif


			sol.m_rerankerScore += pr->pri->rerankerF2DockScoreWeight * origScore;

			double lowScoreFactor = 0.05;

			if ( lowDiv > 0 ) sol.m_rerankerScore *= ( lowScoreFactor / lowDiv );

			sol.m_Value = sol.m_rerankerScore;

			//           sol.m_hbondValue = ( lowDivID * 100 + sol.m_clusterPenalty ) * pr->functionScaleFactor;           

			tempTV->updateTopValues( sol );
		}
	}

	m = tempTV->getCurrentNumberOfPositions( );

	while ( m-- )
	{
		tempTV->extractMin( sol );

		sol.m_Value = 1;
		sols2[ m ] = sol;
	}

	pr->pri->bands[ 2 * pr->pri->numBands ] = pr->pri->bands[ 2 * pr->pri->numBands + 1 ] = nr;

	for ( int j = 0, k = 0; j <= pr->pri->numBands; j++ )
	{
		int maxK = ( ( pr->pri->bands[ 2 * j ] > nr ) ? nr : pr->pri->bands[ 2 * j ] );
		int l = 0;

		while ( k < maxK )
		{
			while ( ( sols2[ l ].m_Value == 0 ) || ( sols2[ l ].m_origRank > pr->pri->bands[ 2 * j + 1 ] ) ) l++;

			sols[ k++ ] = sols2[ l ];
			sols2[ l++ ].m_Value = 0;
		}
	}

	m = n;

	while ( m-- )
	{
		sol = sols[ m ];
		sol.m_Value = ( n - m ) * pr->functionScaleFactor;
		sol.m_rerankerRank = m + 1;

		globalTV->updateTopValues( sol );
	}

	if ( ( pr->pri->complexType == 'A' ) && pr->pri->applyAntibodyFilter )
	{
		if ( cFilterL3 != NULL ) delete cFilterL3;
		if ( cFilterH3 != NULL ) delete cFilterH3;
		if ( cFilterH2 != NULL ) delete cFilterH2;
		if ( cFilterL1H1 != NULL ) delete cFilterL1H1;
	}
	else if ( ( pr->pri->complexType == 'E' ) && pr->pri->applyEnzymeFilter )
		delete cFilter;

	if ( pr->pri->applyResidueContactFilter )
		delete rcFilter;

	delete[ ] sols;
	delete[ ] sols2;
	delete tempTV;
}





int saveGrid( PARAMS_IN *pr )
{
	int interpFuncExtent;
	double scale;
	pr->dockVolume = false;
	if ( !computeGridParameters( pr, &interpFuncExtent, &scale ) ) return -1;

	double blobbiness = pr->blobbiness;
	int gridSize = pr->gridSize;
	int numFreq = pr->numFreq;
	double alpha = pr->alpha;

	double skinSkinWeight = pr->skinSkinWeight;
	double coreCoreWeight = pr->coreCoreWeight;

	double elecRadiusInGrids = pr->elecRadiusInGrids;

	double bandwidth = pr->bandwidth;
	double gradFactor = pr->gradFactor;

	bool singleLayerLigandSkin = pr->singleLayerLigandSkin;
	bool smoothSkin = pr->smoothSkin;
	bool curvatureWeightedStaticMol = pr->curvatureWeightedStaticMol;
	bool curvatureWeightedMovingMol = pr->curvatureWeightedMovingMol;
	double curvatureWeightingRadius = pr->curvatureWeightingRadius;
	bool spreadReceptorSkin = pr->spreadReceptorSkin;

	double real_magnitude = sqrt( skinSkinWeight );
	double imag_magnitude = sqrt( coreCoreWeight );

	double numFreq3 = numFreq * numFreq * numFreq;

	FFTW_complex* gridA = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
	FFTW_DATA_TYPE* elecGridA = ( FFTW_DATA_TYPE * ) FFTW_malloc( sizeof( FFTW_DATA_TYPE ) * numFreq3 );

	FFTW_complex* gridB = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
	FFTW_DATA_TYPE* elecGridB = ( FFTW_DATA_TYPE * ) FFTW_malloc( sizeof( FFTW_DATA_TYPE ) * numFreq3 );

	// static molecule
	double *xkAOrig = pr->xkAOrig;
	double *ykAOrig = pr->ykAOrig;
	double *zkAOrig =pr->zkAOrig;
	int numCentersA = pr->numCentersA;
	char *typeA = pr->typeA;
	char *hbondTypeA = pr->hbondTypeA;
	float *radiiA = pr->radiiA;
	float *chargesA = pr->chargesA;
	char *fileNameA = pr->staticMoleculeF2d;

	// moving molecule
	double *xkBOrig = pr->xkBOrig;
	double *ykBOrig = pr->ykBOrig;
	double *zkBOrig = pr->zkBOrig;
	int numCentersB = pr->numCentersB;
	char *typeB = pr->typeB;
	char *hbondTypeB = pr->hbondTypeB;
	float *radiiB = pr->radiiB;
	float *chargesB = pr->chargesB;

	char *fileNameB = pr->movingMoleculeF2d;

	double *xkA = new double[ numCentersA ];
	double *ykA = new double[ numCentersA ];
	double *zkA = new double[ numCentersA ];
	float *rkA = new float[ numCentersA ];

	double *xkB = new double[ numCentersB ];
	double *ykB = new double[ numCentersB ];
	double *zkB = new double[ numCentersB ];
	float *rkB = new float[ numCentersB ];

	FFTW_complex *fkA = 0, *fkB = 0;
	FFTW_complex *fkAElec = 0, *fkBElec = 0;
	FFTW_complex *fkAHbond = 0, *fkBHbond = 0;

	SmoothingFunction *smoothingFunction;

	double n = ( int ) ( alpha * numFreq );

	smoothingFunction = new Gaussian( alpha, interpFuncExtent, n, gridSize );
	std::cout<<"typeA "<<typeA<<" typeB "<<typeB<<" numCentersA "<<numCentersA<<" "<<"numCentersB " <<numCentersB
		<<" A origin "<<*xkAOrig<<" "<<*ykAOrig<<" "<<*zkAOrig<<" B origin "<<*xkBOrig<<" "<<*ykBOrig<<" "<<*zkBOrig
		<<std::endl;

	if ( !build_fks( typeA, numCentersA, xkAOrig, ykAOrig, zkAOrig, radiiA, chargesA, hbondTypeA, NULL,
				&fkA, 1, &fkAElec, 1, &fkAHbond,
				0, 0, 0, 0, 0, NULL, false, NULL,
				0, 0, NULL,
				real_magnitude, imag_magnitude, true,
				singleLayerLigandSkin, curvatureWeightedStaticMol, curvatureWeightingRadius,
				bandwidth, gradFactor, pr ) ) return -1;

	double xTransA = 0, yTransA = 0, zTransA = 0;
	// center A
	if ( !center( xkAOrig, ykAOrig, zkAOrig, numCentersA, &xTransA, &yTransA, &zTransA ) ) return -1;

	// transform A
	if ( !transformAndNormalize( xkAOrig, ykAOrig, zkAOrig, radiiA, xkA, ykA, zkA, rkA, numCentersA, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, scale ) ) return -1;

	gridding( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA, blobbiness, numFreq, interpFuncExtent, smoothSkin, smoothingFunction, gridA, spreadReceptorSkin );

	griddingElec( numCentersA, xkA, ykA, zkA, rkA, typeA, fkAElec,
			blobbiness, numFreq, elecRadiusInGrids, elecGridA, false, false );

	if ( !writeGrid( gridA, elecGridA, n, - xTransA, - yTransA, - zTransA, scale, fileNameA,
				pr->staticMoleculeSCReRaw, pr->staticMoleculeSCImRaw, pr->staticMoleculeElecReRaw ) ) return -1;

	if ( !build_fks( typeB, numCentersB, xkBOrig, ykBOrig, zkBOrig, radiiB, chargesB, hbondTypeB, NULL,
				&fkB, 1, &fkBElec, 1, &fkBHbond,
				0, 0, 0, 0, 0, NULL, false, NULL,
				0, 0, NULL,
				real_magnitude, imag_magnitude, false,
				singleLayerLigandSkin, curvatureWeightedMovingMol, curvatureWeightingRadius,
				bandwidth, gradFactor, pr ) ) return -1;

	double xTransB = 0, yTransB = 0, zTransB = 0;
	// center B
	if ( !center( xkBOrig, ykBOrig, zkBOrig, numCentersB, &xTransB, &yTransB, &zTransB ) ) return -1;

	// transform B
	if ( !transformAndNormalize( xkBOrig, ykBOrig, zkBOrig, radiiB, xkB, ykB, zkB, rkB, numCentersB, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1 /* do not want to flip (2) */, scale ) ) return -1;

	gridding( numCentersB, xkB, ykB, zkB, rkB, typeB, fkB, blobbiness, numFreq, interpFuncExtent,
			smoothSkin, smoothingFunction, gridB, false );

	griddingElec( numCentersB, xkB, ykB, zkB, rkB, typeB, fkBElec,
			blobbiness, numFreq, elecRadiusInGrids, elecGridB, false, true );

	if ( !writeGrid( gridB, elecGridB, n, - xTransB, - yTransB, - zTransB, scale, fileNameB,
				pr->movingMoleculeSCReRaw, pr->movingMoleculeSCImRaw, pr->movingMoleculeElecReRaw ) ) return -1;


	FFTW_free( gridA );
	FFTW_free( elecGridA );
	FFTW_free( fkA );
	FFTW_free( fkAElec );

	delete [ ] xkA;
	delete [ ] ykA;
	delete [ ] zkA;

	FFTW_free( gridB );
	FFTW_free( elecGridB );
	FFTW_free( fkB );
	FFTW_free( fkBElec );

	delete [ ] xkB;
	delete [ ] ykB;
	delete [ ] zkB;

	delete smoothingFunction;

	return 0;
}


int inSphere( int x[ 3 ], void *data )
{
	INSPHERE_DATA *d = ( INSPHERE_DATA * ) data;
	double dx2, dy2, dz2;

	dx2 = ( x[ 0 ] - d->cx ) * ( x[ 0 ] - d->cx );
	dy2 = ( x[ 1 ] - d->cy ) * ( x[ 1 ] - d->cy );
	dz2 = ( x[ 2 ] - d->cz ) * ( x[ 2 ] - d->cz );

	return ( dx2 + dy2 + dz2 <= d->r2 );
}


void createValidOutputMap( FFTW_complex *staticMol, int n, int minRadMovingMol, int maxRadMovingMol, bool narrowBand, int *validOutputMap )
{
	int n3 = n * n * n;

	for ( int c = 0; c < n3; c++ )
		validOutputMap[ c ] = 1;

	if ( narrowBand )
	{
		int nn = n >> 1;
		double ofs = ( n & 1 ) ? 0.0 : 0.5;
		double eps = 0.00001;

		FFTW_complex *sMol = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * n3 );
		FFTW_complex *mMol = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * n3 );
		FFTW_complex *smProd = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * n3 );

		for ( int c = 0; c < n3; c++ )
		{
			if ( ( staticMol[ c ][ 0 ] > eps ) || ( - staticMol[ c ][ 0 ] > eps )
					|| ( staticMol[ c ][ 1 ] > eps ) || ( - staticMol[ c ][ 1 ] > eps ) )
				sMol[ c ][ 0 ] = 1.0;
			else sMol[ c ][ 0 ] = 0.0;

			sMol[ c ][ 1 ] = 0.0;

			mMol[ c ][ 0 ] = mMol[ c ][ 1 ] = 0.0;
		}

		//      maxRadMovingMol >>= 1;
		//      maxRadMovingMol += 2;
		//      if ( maxRadMovingMol > nn ) maxRadMovingMol = nn;

		int maxRad2 = maxRadMovingMol * maxRadMovingMol;
		int nl = nn, nr = nn + ( n & 1 );

		if ( maxRadMovingMol < nl )
		{
			nl = maxRadMovingMol;
			if ( maxRadMovingMol + 1 <= nr ) nr = maxRadMovingMol + 1;
			else nr = maxRadMovingMol;
		}

		for ( int iz = -nl; iz < nr; iz++ )
			for ( int iy = -nl; iy < nr; iy++ )
				for ( int ix = -nl; ix < nr; ix++ )
				{
					double d2 = ( ix + ofs ) * ( ix + ofs ) + ( iy + ofs ) * ( iy + ofs ) + ( iz + ofs ) * ( iz + ofs );

					if ( d2 <= maxRad2 )
					{
						int c = ( iz + nn ) * n * n + ( iy + nn ) * n + ( ix + nn );

						mMol[ c ][ 0 ] = 1;
					}
				}

		FFTW_plan sForwardPlan = FFTW_plan_dft_3d( n, n, n, sMol, sMol, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );
		FFTW_plan mForwardPlan = FFTW_plan_dft_3d( n, n, n, mMol, mMol, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );

		FFTW_plan smBackwardPlan = FFTW_plan_dft_3d( n, n, n, smProd, smProd, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE );

		FFTW_execute( sForwardPlan );
		FFTW_execute( mForwardPlan );

		for ( int c = 0; c < n3; c++ )
		{
			smProd[ c ][ 0 ] = sMol[ c ][ 0 ] * mMol[ c ][ 0 ] - sMol[ c ][ 1 ] * mMol[ c ][ 1 ];
			smProd[ c ][ 1 ] = sMol[ c ][ 0 ] * mMol[ c ][ 1 ] + sMol[ c ][ 1 ] * mMol[ c ][ 0 ];
		}

		FFTW_execute( smBackwardPlan );

		int t = 0;
		for ( int c = 0; c < n3; c++ )
			if ( smProd[ c ][ 0 ] < 0.5 * n3 )
			{
				validOutputMap[ c ] = 0;
				t++;
			}

		printf( "\nmaxRadMovingMol = %d", maxRadMovingMol );
		printf( "\nExcluded ( outer band ): %lf \%\n", ( t * 100.0 ) / n3 );

		//      for ( int c = 0; c < n3; c++ )
		//          mMol[ c ][ 0 ] = mMol[ c ][ 1 ] = 0.0;
		//
		//      minRadMovingMol -= 2;
		//      if ( minRadMovingMol < 0 ) minRadMovingMol = 0;
		//
		//      int minRad2 = minRadMovingMol * minRadMovingMol;
		//      nl = nn;
		//      nr = nn + ( n & 1 );
		//
		//      if ( minRadMovingMol < nl )
		//        {
		//          nl = minRadMovingMol;
		//          if ( minRadMovingMol + 1 <= nr ) nr = minRadMovingMol + 1;
		//          else nr = minRadMovingMol;
		//        }
		//
		//      for ( int iz = -nl; iz < nr; iz++ )
		//        for ( int iy = -nl; iy < nr; iy++ )
		//          for ( int ix = -nl; ix < nr; ix++ )
		//            {
		//              double d2 = ( ix + ofs ) * ( ix + ofs ) + ( iy + ofs ) * ( iy + ofs ) + ( iz + ofs ) * ( iz + ofs );
		//
		//              if ( d2 <= minRad2 )
		//                 {
		//                   int c = ( iz + nn ) * n * n + ( iy + nn ) * n + ( ix + nn );
		//
		//                   mMol[ c ][ 0 ] = 1;
		//                 }
		//            }
		//
		//      FFTW_execute( mForwardPlan );
		//
		//      for ( int c = 0; c < n3; c++ )
		//	 {
		//	   smProd[ c ][ 0 ] = sMol[ c ][ 0 ] * mMol[ c ][ 0 ] - sMol[ c ][ 1 ] * mMol[ c ][ 1 ];
		//	   smProd[ c ][ 1 ] = sMol[ c ][ 0 ] * mMol[ c ][ 1 ] + sMol[ c ][ 1 ] * mMol[ c ][ 0 ];
		//	 }
		//
		//      FFTW_execute( smBackwardPlan );
		//
		//      t = 0;
		//      for ( int c = 0; c < n3; c++ )
		//         if ( smProd[ c ][ 0 ] > 0.5 * n3 )
		//           {
		//             validOutputMap[ c ] = 0;
		//             t++;
		//           }
		//
		//      printf( "\nminRadMovingMol = %d", minRadMovingMol );
		//      printf( "\nExcluded ( inner band ): %lf \%\n", ( t * 100.0 ) / n3 );
		//
		//      t = 0;
		//      for ( int z = 0, k = 0; z < n; z++ )
		//        {
		//          int iz = z;
		//          if ( iz > n / 2 ) iz -= n;
		//
		//          for ( int y = 0; y < n; y++ )
		//            {
		//              int iy = y;
		//              if ( iy > n / 2 ) iy -= n;
		//
		//              for ( int x = 0; x < n; x++, k++ )
		//                {
		//                  int ix = x;
		//                  if ( ix > n / 2 ) ix -= n;
		//
		//                  int ik = ( ( iz + ( n / 2 ) ) * n + ( iy + ( n / 2 ) ) ) * n + ( ix + ( n / 2 ) );
		//
		//                  if ( ( staticMol[ ik ][ 1 ] > eps ) || ( - staticMol[ ik ][ 1 ] > eps ) )
		//                    {
		//                     if ( validOutputMap[ k ] )
		//                       {
		//                         validOutputMap[ k ] = 0;
		//                         t++;
		//                       }
		//                    }
		//                }
		//            }
		//        }
		//
		//      if ( t > 0 ) printf( "\nExcluded ( core ): %lf \%\n", ( t * 100.0 ) / n3 );

		FFTW_destroy_plan( sForwardPlan );
		FFTW_destroy_plan( mForwardPlan );
		FFTW_destroy_plan( smBackwardPlan );

		FFTW_free( sMol );
		FFTW_free( mMol );
		FFTW_free( smProd );
	}
}


int outputValid( int x[ 3 ], void *data )
{
	VALID_OUTPUT_DATA *d = ( VALID_OUTPUT_DATA * ) data;
	int k = ( x[ 2 ] * ( d->n ) + x[ 1 ] ) * ( d->n ) + x[ 0 ];

	return ( d->validOutputMap[ k ] );
}


//void markCluster( int x, int y, int z, int *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad, int addVal )
//{
//    double dist = clusterTransRad / gridSpacing;
//    double dist2 = dist * dist;
//    int dCells = ( int ) ceil( dist );
//
//    if ( x > ( numFreq >> 1 ) ) x -= numFreq;
//    if ( y > ( numFreq >> 1 ) ) y -= numFreq;
//    if ( z > ( numFreq >> 1 ) ) z -= numFreq;
//
//    int lx = x - dCells, hx = x + dCells;
//    int ly = y - dCells, hy = y + dCells;
//    int lz = z - dCells, hz = z + dCells;
//
//    if ( lx <= - ( numFreq >> 1 ) ) lx = - ( numFreq >> 1 ) + 1;
//    if ( ly <= - ( numFreq >> 1 ) ) ly = - ( numFreq >> 1 ) + 1;
//    if ( lz <= - ( numFreq >> 1 ) ) lz = - ( numFreq >> 1 ) + 1;
//
//    if ( hx > ( numFreq >> 1 ) ) hx = ( numFreq >> 1 );
//    if ( hy > ( numFreq >> 1 ) ) hy = ( numFreq >> 1 );
//    if ( hz > ( numFreq >> 1 ) ) hz = ( numFreq >> 1 );
//
//    for ( int xx = lx, dx = ( x - lx ) * ( x - lx ); xx <= hx; xx++ )
//      {
//       for ( int yy = ly, dxy = dx + ( y - ly ) * ( y - ly ); yy <= hy; yy++ )
//         {
//          int dxy = dx + ( y - yy ) * ( y - yy );
//
//          for ( int zz = lz, dxyz = dxy + ( z - lz ) * ( z - lz ); zz <= hz; zz++ )
//             {
//              if ( dxyz <= dist2 )
//                {
//                  int nx = xx, ny = yy, nz = zz;
//
//                  if ( nx < 0 ) nx += numFreq;
//                  if ( ny < 0 ) ny += numFreq;
//                  if ( nz < 0 ) nz += numFreq;
//
//                  int cc = ( nx * numFreq + ny ) * numFreq + nz;
//
//                  markedPeaks[ cc ] += addVal;
//                }
//
//              dxyz -= ( ( ( z - zz ) << 1 ) - 1 );
//             }
//
//           dxy -= ( ( ( y - yy ) << 1 ) - 1 );
//         }
//
//        dx -= ( ( ( x - xx ) << 1 ) - 1 );
//      }
//}
//
//
//
//int getClusterValue( int x, int y, int z, int *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad )
//{
//    double dist = clusterTransRad / gridSpacing;
//    double dist2 = dist * dist;
//    int dCells = ( int ) ceil( dist );
//
//    if ( x > ( numFreq >> 1 ) ) x -= numFreq;
//    if ( y > ( numFreq >> 1 ) ) y -= numFreq;
//    if ( z > ( numFreq >> 1 ) ) z -= numFreq;
//
//    int lx = x - dCells, hx = x + dCells;
//    int ly = y - dCells, hy = y + dCells;
//    int lz = z - dCells, hz = z + dCells;
//
//    if ( lx <= - ( numFreq >> 1 ) ) lx = - ( numFreq >> 1 ) + 1;
//    if ( ly <= - ( numFreq >> 1 ) ) ly = - ( numFreq >> 1 ) + 1;
//    if ( lz <= - ( numFreq >> 1 ) ) lz = - ( numFreq >> 1 ) + 1;
//
//    if ( hx > ( numFreq >> 1 ) ) hx = ( numFreq >> 1 );
//    if ( hy > ( numFreq >> 1 ) ) hy = ( numFreq >> 1 );
//    if ( hz > ( numFreq >> 1 ) ) hz = ( numFreq >> 1 );
//
//    int cVal = 0;
//
//    for ( int xx = lx, dx = ( x - lx ) * ( x - lx ); xx <= hx; xx++ )
//      {
//       for ( int yy = ly, dxy = dx + ( y - ly ) * ( y - ly ); yy <= hy; yy++ )
//         {
//          int dxy = dx + ( y - yy ) * ( y - yy );
//
//          for ( int zz = lz, dxyz = dxy + ( z - lz ) * ( z - lz ); zz <= hz; zz++ )
//             {
//              if ( dxyz <= dist2 )
//                {
//                  int nx = xx, ny = yy, nz = zz;
//
//                  if ( nx < 0 ) nx += numFreq;
//                  if ( ny < 0 ) ny += numFreq;
//                  if ( nz < 0 ) nz += numFreq;
//
//                  int cc = ( nx * numFreq + ny ) * numFreq + nz;
//
//                  cVal += markedPeaks[ cc ];
//                }
//
//              dxyz -= ( ( ( z - zz ) << 1 ) - 1 );
//             }
//
//           dxy -= ( ( ( y - yy ) << 1 ) - 1 );
//         }
//
//        dx -= ( ( ( x - xx ) << 1 ) - 1 );
//      }
//
//    return cVal;
//}


void markCluster( int x, int y, int z, double *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad, double val )
{
	double dist = clusterTransRad / gridSpacing;
	double dist2 = dist * dist;
	int dCells = ( int ) ceil( dist );

	if ( x > ( numFreq >> 1 ) ) x -= numFreq;
	if ( y > ( numFreq >> 1 ) ) y -= numFreq;
	if ( z > ( numFreq >> 1 ) ) z -= numFreq;

	int lx = x - dCells, hx = x + dCells;
	int ly = y - dCells, hy = y + dCells;
	int lz = z - dCells, hz = z + dCells;

	if ( lx <= - ( numFreq >> 1 ) ) lx = - ( numFreq >> 1 ) + 1;
	if ( ly <= - ( numFreq >> 1 ) ) ly = - ( numFreq >> 1 ) + 1;
	if ( lz <= - ( numFreq >> 1 ) ) lz = - ( numFreq >> 1 ) + 1;

	if ( hx > ( numFreq >> 1 ) ) hx = ( numFreq >> 1 );
	if ( hy > ( numFreq >> 1 ) ) hy = ( numFreq >> 1 );
	if ( hz > ( numFreq >> 1 ) ) hz = ( numFreq >> 1 );

	for ( int xx = lx, dx = ( x - lx ) * ( x - lx ); xx <= hx; xx++ )
	{
		for ( int yy = ly, dxy = dx + ( y - ly ) * ( y - ly ); yy <= hy; yy++ )
		{
			int dxy = dx + ( y - yy ) * ( y - yy );

			for ( int zz = lz, dxyz = dxy + ( z - lz ) * ( z - lz ); zz <= hz; zz++ )
			{
				if ( dxyz <= dist2 )
				{
					int nx = xx, ny = yy, nz = zz;

					if ( nx < 0 ) nx += numFreq;
					if ( ny < 0 ) ny += numFreq;
					if ( nz < 0 ) nz += numFreq;

					int cc = ( nx * numFreq + ny ) * numFreq + nz;

					if ( val > markedPeaks[ cc ] ) markedPeaks[ cc ] = val;
				}

				dxyz -= ( ( ( z - zz ) << 1 ) - 1 );
			}

			dxy -= ( ( ( y - yy ) << 1 ) - 1 );
		}

		dx -= ( ( ( x - xx ) << 1 ) - 1 );
	}
}



int getClusterValue( int x, int y, int z, double *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad, double val )
{
	double dist = clusterTransRad / gridSpacing;
	double dist2 = dist * dist;
	int dCells = ( int ) ceil( dist );

	if ( x > ( numFreq >> 1 ) ) x -= numFreq;
	if ( y > ( numFreq >> 1 ) ) y -= numFreq;
	if ( z > ( numFreq >> 1 ) ) z -= numFreq;

	int lx = x - dCells, hx = x + dCells;
	int ly = y - dCells, hy = y + dCells;
	int lz = z - dCells, hz = z + dCells;

	if ( lx <= - ( numFreq >> 1 ) ) lx = - ( numFreq >> 1 ) + 1;
	if ( ly <= - ( numFreq >> 1 ) ) ly = - ( numFreq >> 1 ) + 1;
	if ( lz <= - ( numFreq >> 1 ) ) lz = - ( numFreq >> 1 ) + 1;

	if ( hx > ( numFreq >> 1 ) ) hx = ( numFreq >> 1 );
	if ( hy > ( numFreq >> 1 ) ) hy = ( numFreq >> 1 );
	if ( hz > ( numFreq >> 1 ) ) hz = ( numFreq >> 1 );

	int cVal = 0;

	for ( int xx = lx, dx = ( x - lx ) * ( x - lx ); xx <= hx; xx++ )
	{
		for ( int yy = ly, dxy = dx + ( y - ly ) * ( y - ly ); yy <= hy; yy++ )
		{
			int dxy = dx + ( y - yy ) * ( y - yy );

			for ( int zz = lz, dxyz = dxy + ( z - lz ) * ( z - lz ); zz <= hz; zz++ )
			{
				if ( dxyz <= dist2 )
				{
					int nx = xx, ny = yy, nz = zz;

					if ( nx < 0 ) nx += numFreq;
					if ( ny < 0 ) ny += numFreq;
					if ( nz < 0 ) nz += numFreq;

					int cc = ( nx * numFreq + ny ) * numFreq + nz;

					if ( markedPeaks[ cc ] >= val ) cVal++;
				}

				dxyz -= ( ( ( z - zz ) << 1 ) - 1 );
			}

			dxy -= ( ( ( y - yy ) << 1 ) - 1 );
		}

		dx -= ( ( ( x - xx ) << 1 ) - 1 );
	}

	return cVal;
}



bool initClashFilter( PARAMS_IN *pr, clashFilter **cFilter )
{
	int numStaticAtoms = 0, numMovingAtoms = 0;//pr->numCentersB;

	for ( int i = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' ) numStaticAtoms++;

	for ( int i = 0; i < pr->numCentersB; i++ )
		if ( pr->typeB[ i ] == 'I' ) numMovingAtoms++;

	double *staticAtoms = ( double * ) FFTW_malloc( 5 * numStaticAtoms * sizeof( double ) );
	double *movingAtoms = ( double * ) FFTW_malloc( 5 * numMovingAtoms * sizeof( double ) );

	if ( ( staticAtoms == NULL ) || ( movingAtoms == NULL ) ) return false;

	printf( "\nnumStaticAtoms = %d ( numCentersA = %d ), numMovingAtoms = %d ( numCentersB = %d )\n", numStaticAtoms, pr->numCentersA, numMovingAtoms, pr->numCentersB );

	for ( int i = 0, j = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' )
		{
			staticAtoms[ j++ ] = pr->xkAOrig[ i ];
			staticAtoms[ j++ ] = pr->ykAOrig[ i ];
			staticAtoms[ j++ ] = pr->zkAOrig[ i ];

			staticAtoms[ j++ ] = pr->chargesA[ i ];

			staticAtoms[ j++ ] = pr->radiiA[ i ];
		}

	for ( int i = 0, j = 0; i < pr->numCentersB; i++ )
		if ( pr->typeB[ i ] == 'I' )
		{
			movingAtoms[ j++ ] = pr->xkBOrig[ i ];
			movingAtoms[ j++ ] = pr->ykBOrig[ i ];
			movingAtoms[ j++ ] = pr->zkBOrig[ i ];

			movingAtoms[ j++ ] = pr->chargesB[ i ];

			movingAtoms[ j++ ] = pr->radiiB[ i ];
		}

	( *cFilter ) = new clashFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, false );

	//  ( *cFilter )->setNumThreads( 1 );
	//  ( *cFilter )->setProximityFactors( 0.4, 0.25, 0.8 );

	FFTW_free( staticAtoms );
	FFTW_free( movingAtoms );

	return true;
}

//essentially a duplicate of initClashFilter, except that one of the inputs is a raw file enclosing the forbidden volume
bool initForbiddenVolumeFilter( PARAMS_IN *pr, clashFilter **cFilter )
{

	std::cout<<"Entering initForbiddenVolumeFilter().\n";
	int numStaticAtoms = 0, numMovingAtoms = 0;//pr->numCentersB;

	//std::ifstream rawfile(pr->forbiddenVolFileName);

	double* staticAtoms;
	double x,y,z;
	double charge = 0.0;
	double radius = 3.0;

	if(pr->forbiddenVolumeFileType==0) { //rawfile

		std::cout<<"Forbidden volume bounded by surface "<<pr->forbiddenVolFileName<<std::endl;

		FILE* rawfp = fopen(pr->forbiddenVolFileName, "r");

		//  int i = 0;
		if(rawfp) 
		{
			int dum;

			if(!(fscanf(rawfp, "%d %d", &numStaticAtoms, &dum)==2))
				std::cerr<<"Cannot read first line\n";
		}

		else { 

			std::cerr<<"Cannot read forbidden volume raw file. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;

		}

		staticAtoms = ( double * ) FFTW_malloc( 5 * numStaticAtoms * sizeof( double ) );

		if ( ( staticAtoms == NULL ) ) { 

			std::cerr<<"Cannot allocate memory. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;

		}

		for ( int i = 0, j = 0; i < numStaticAtoms; i++ )
		{
			if((fscanf(rawfp, "%lf %lf %lf", &x, &y, &z)==3)) {

				std::cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<std::endl; 
				staticAtoms[ j++ ] = x;
				staticAtoms[ j++ ] = y;
				staticAtoms[ j++ ] = z;

				staticAtoms[ j++ ] = charge;

				staticAtoms[ j++ ] = radius;

			}
		}

		fclose(rawfp);

	}

	else if(pr->forbiddenVolumeFileType==2) {//molecule
		FILE *fp;
		std::cout<<"Forbidden volume bounded by molecule  "<<pr->forbiddenVolFileName<<std::endl;
		if( (fp = fopen(pr->forbiddenVolFileName, "r")) == NULL ) {

			std::cerr<<"Cannot read forbidden volume pqr file. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;

		}

		while(true) {
			char str[40];

			if(fscanf(fp,"%s",str)!=0){
				break;
			}

			string temp(str);

			if(temp.compare("ATOM")==0 || temp.compare("HETATM")==0){
				numStaticAtoms++;
				fgets(str, 1999, fp);
			}
		}
		fclose(fp);

		staticAtoms = ( double * ) FFTW_malloc( 5 * numStaticAtoms * sizeof( double ) );

		if ( ( staticAtoms == NULL ) ) {

			std::cerr<<"Cannot allocate memory. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;

		}

		if( (fp = fopen(pr->forbiddenVolFileName, "r")) == NULL ) 
			return false;

		int j=0;

		while(true) {
			char str[40];

			if(fscanf(fp,"%s",str)!=0){
				break;
			}

			string temp(str);
			int atomNum;
			char atomName[50];
			char resName[50];

			if(temp.compare("ATOM")==0 || temp.compare("HETATM")==0){
				if(fscanf(fp,"%d %s %s %lf %lf %lf %lf %lf", &atomNum, &atomName, &resName, &x, &y, &z, &charge, &radius)==8){
					staticAtoms[ j++ ] = x;
					staticAtoms[ j++ ] = y;
					staticAtoms[ j++ ] = z;

					staticAtoms[ j++ ] = charge;

					staticAtoms[ j++ ] = radius;
				}
			}
		}

		fclose(fp);
	}

	else if(pr->forbiddenVolumeFileType==1) {//bounding box

		FILE *fp;
		std::cout<<"Forbidden volume bounded by bbox  "<<pr->forbiddenVolFileName<<std::endl;
		if( (fp = fopen(pr->forbiddenVolFileName, "r")) == NULL ) {

			std::cerr<<"Cannot read forbidden volume bbox file. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;

		}

		if(fscanf(fp,"%lf %lf %lf %lf %lf %lf", &pr->forbiddenBBoxXMin, &pr->forbiddenBBoxYMin, &pr->forbiddenBBoxZMin, &pr->forbiddenBBoxXMax, &pr->forbiddenBBoxYMax, &pr->forbiddenBBoxZMax)!=6) {

			std::cerr<<"Wrong bbox file format. Forbidden volume filter will not be applied.\n";
			pr->applyForbiddenVolumeFilter = false;
			return false;
		}

		pr->centroidxB = 0;
		pr->centroidyB = 0;
		pr->centroidzB = 0;

		for ( int i = 0; i < pr->numCentersB; i++ )
		{
			if(i==0) {

				pr->bboxXBMin = pr->xkBOrig[ i ]; pr->bboxYBMin  = pr->ykBOrig[ i ]; pr->bboxZBMin = pr->zkBOrig[ i ];
				pr->bboxXBMax = pr->bboxXBMin;  pr->bboxYBMax = pr->bboxYBMin; pr->bboxZBMax=pr->bboxZBMin;

			}

			pr->centroidxB += pr->xkBOrig[ i ];
			pr->centroidyB += pr->ykBOrig[ i ];
			pr->centroidzB += pr->zkBOrig[ i ];

			pr->bboxXBMin = std::min(pr->xkBOrig[i], pr->bboxXBMin);
			pr->bboxYBMin = std::min(pr->ykBOrig[i], pr->bboxYBMin);
			pr->bboxZBMin = std::min(pr->zkBOrig[i], pr->bboxZBMin);

			pr->bboxXBMax = std::max(pr->xkBOrig[i], pr->bboxXBMax);
			pr->bboxYBMax = std::max(pr->ykBOrig[i], pr->bboxYBMax);
			pr->bboxZBMax = std::max(pr->zkBOrig[i], pr->bboxZBMax);

		}

		pr->centroidxB /= pr->numCentersB;
		pr->centroidyB /= pr->numCentersB;
		pr->centroidzB /= pr->numCentersB;


		return true;

	}

	else {

		pr->applyForbiddenVolumeFilter = false;
		return false;

	}

	for ( int i = 0; i < pr->numCentersB; i++ )
		if ( pr->typeB[ i ] == 'I' ) numMovingAtoms++;

	double *movingAtoms = ( double * ) FFTW_malloc( 5 * numMovingAtoms * sizeof( double ) );

	if ( ( movingAtoms == NULL ) ) return false;

	for ( int i = 0, j = 0; i < pr->numCentersB; i++ )
		if ( pr->typeB[ i ] == 'I' )
		{
			movingAtoms[ j++ ] = pr->xkBOrig[ i ];
			movingAtoms[ j++ ] = pr->ykBOrig[ i ];
			movingAtoms[ j++ ] = pr->zkBOrig[ i ];

			movingAtoms[ j++ ] = pr->chargesB[ i ];

			movingAtoms[ j++ ] = pr->radiiB[ i ];
		}

	( *cFilter ) = new clashFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, false, 1.0 );

	//  ( *cFilter )->setNumThreads( 1 );
	//  ( *cFilter )->setProximityFactors( 0.4, 0.25, 0.8 );

	FFTW_free( staticAtoms );
	FFTW_free( movingAtoms );

	std::cout<<"Leaving initForbiddenVolumeFilter().\n";

	return true;
}


void initLJFilter( PARAMS_IN *pr, fastLJ **ljFilter )
{
	( *ljFilter ) = new fastLJ( pr->paramFile );

	( *ljFilter )->setPrintStatus( false );
}



//bool updatePseudoGsolCutoff( PARAMS_IN *pr, double pseudoGsol )
//{
//   double factor = 0.55, minVal = -5000;
//   bool filtered = false;
//
//   pthread_mutex_lock( &( pr->pseudoGsolLock ) );
//   if ( pseudoGsol > pr->pseudoGsolCutoff ) filtered = true;
//   else if ( ( pseudoGsol >= minVal ) && ( factor * pseudoGsol < pr->pseudoGsolCutoff ) ) pr->pseudoGsolCutoff = factor * pseudoGsol;
//   pthread_mutex_unlock( &( pr->pseudoGsolLock ) );
//
//   return filtered;
//}


bool updatePseudoGsolCutoff( PARAMS_IN *pr, double pseudoGsol )
{
	double factor = 4, maxVal = pr->hydroRatioTolerance;
	static int numAttempts = 0;
	static double sum = 0;
	static int maxAttempts = 0;
	bool filtered = false;

	pthread_mutex_lock( &( pr->pseudoGsolLock ) );

	//   numAttempts++;
	//   sum += pseudoGsol;
	//   if ( pseudoGsol < pr->pseudoGsolCutoff * ( sum / numAttempts ) ) filtered = true;

	if ( pseudoGsol < pr->pseudoGsolCutoff ) filtered = true;
	else if ( ( pseudoGsol <= maxVal ) && ( pseudoGsol > factor * pr->pseudoGsolCutoff ) )
	{
		//             pr->pseudoGsolCutoff = pseudoGsol / factor;
		numAttempts++;
		sum += pseudoGsol;

		if ( numAttempts > maxAttempts )
		{
			pr->pseudoGsolCutoff = ( sum / numAttempts ) / factor;
			numAttempts = 0;
			sum = 0;
			maxAttempts += 2;
		}
	}

	pthread_mutex_unlock( &( pr->pseudoGsolLock ) );

	return filtered;
}


void initPseudoGsolFilter( PARAMS_IN *pr )
{
	pthread_mutex_init( &( pr->pseudoGsolLock ), NULL );

	int numCentersA = 0;

	for ( int i = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' ) numCentersA++;

	double *atomsA, *atomsB;

	atomsA = ( double * ) malloc( ( 5 * numCentersA ) * sizeof( double ) );
	atomsB = ( double * ) malloc( ( 5 * pr->numCentersB ) * sizeof( double ) );

	if ( ( atomsA == NULL ) || ( atomsB == NULL ) )
	{
		printf( "\nError: Memory allocation failed in initPseudoGsolFilter!\n" );
		exit( 1 );
	}

	for ( int i = 0, j = 0; i < pr->numCentersA; i++ )
		if ( pr->typeA[ i ] == 'I' )
		{
			atomsA[ 5 * j + 0 ] = pr->xkAOrig[ i ];
			atomsA[ 5 * j + 1 ] = pr->ykAOrig[ i ];
			atomsA[ 5 * j + 2 ] = pr->zkAOrig[ i ];
			atomsA[ 5 * j + 3 ] = pr->radiiA[ i ];
			atomsA[ 5 * j + 4 ] = pr->hydrophobicityA[ i ];
			j++;
		}

	for ( int i = 0; i < pr->numCentersB; i++ )
	{
		atomsB[ 5 * i + 0 ] = pr->xkBOrig[ i ];
		atomsB[ 5 * i + 1 ] = pr->ykBOrig[ i ];
		atomsB[ 5 * i + 2 ] = pr->zkBOrig[ i ];
		atomsB[ 5 * i + 3 ] = pr->radiiB[ i ];
		atomsB[ 5 * i + 4 ] = pr->hydrophobicityB[ i ];
	}

	pr->pGsol = new pseudoGsol( pr->paramFile, numCentersA, atomsA, pr->numCentersB, atomsB, pr->numThreads );

	free( atomsA );
	free( atomsB );
}

double getRadius( char *atomName, double r )
{
	switch ( atomName[ 0 ] )
	{
		case 'N' : return 1.55;

		case 'C' : return 1.70;

		case 'O' : return 1.40;

		case 'H' : return 1.20;

		case 'S' : return 1.85;

		case 'P' : return 1.90;

		default  : return r;
	}
}


void initDispersionFilter( PARAMS_IN *pr )
{
	int numCentersA = 0, numCentersB = 0;
	double minRad = pr->dispersionMinAtomRadius;

	for ( int i = 0; i < pr->numCentersA; i++ )
		if ( ( pr->typeA[ i ] == 'I' ) && ( ( pr->radiiA[ i ] > 0 ) || ( minRad > 0 ) ) ) numCentersA++;

	for ( int i = 0; i < pr->numCentersB; i++ )
		if ( ( pr->radiiB[ i ] > 0 ) || ( minRad > 0 ) ) numCentersB++;

	double *atomsA, *atomsB;

	atomsA = ( double * ) malloc( ( 5 * numCentersA ) * sizeof( double ) );
	atomsB = ( double * ) malloc( ( 5 * numCentersB ) * sizeof( double ) );

	if ( ( atomsA == NULL ) || ( atomsB == NULL ) )
	{
		printf( "\nError: Memory allocation failed in initDispersionFilter!\n" );
		exit( 1 );
	}

	for ( int i = 0, j = 0; i < pr->numCentersA; i++ )
		if ( ( pr->typeA[ i ] == 'I' ) && ( ( pr->radiiA[ i ] > 0 ) || ( minRad > 0 ) ) )
		{
			atomsA[ 5 * j + 0 ] = pr->xkAOrig[ i ];
			atomsA[ 5 * j + 1 ] = pr->ykAOrig[ i ];
			atomsA[ 5 * j + 2 ] = pr->zkAOrig[ i ];
			atomsA[ 5 * j + 3 ] = pr->chargesA[ i ];
			if ( pr->radiiA[ i ] > 0 ) atomsA[ 5 * j + 4 ] = pr->radiiA[ i ]; //getRadius( pr->atNamesA[ i ], pr->radiiA[ i ] );
			else atomsA[ 5 * j + 4 ] = minRad;
			j++;
		}

	for ( int i = 0, j = 0; i < pr->numCentersB; i++ )
		if ( ( pr->radiiB[ i ] > 0 ) || ( minRad > 0 ) )
		{
			atomsB[ 5 * j + 0 ] = pr->xkBOrig[ i ];
			atomsB[ 5 * j + 1 ] = pr->ykBOrig[ i ];
			atomsB[ 5 * j + 2 ] = pr->zkBOrig[ i ];
			atomsB[ 5 * j + 3 ] = pr->chargesB[ i ];
			if ( pr->radiiB[ i ] > 0 ) atomsB[ 5 * j + 4 ] = pr->radiiB[ i ]; //getRadius( pr->atNamesB[ i ], pr->radiiB[ i ] );
			else atomsB[ 5 * j + 4 ] = minRad;
			j++;
		}

#ifdef RERANK_DEBUG
	pr->dispF = new fastGB::fastDispE( pr->paramFile, numCentersA, atomsA, numCentersB, atomsB, pr->numThreads, 1 );
#else
	pr->dispF = new fastGB::fastDispE( pr->paramFile, numCentersA, atomsA, numCentersB, atomsB, 1, pr->numThreads );
#endif

	free( atomsA );
	free( atomsB );
}


void initHBondFilter( PARAMS_IN *pr )
{    
#ifdef LIBMOL_FOUND
	pr->hbFilter = new hbondFilter(pr->staticMoleculePdb, pr->movingMoleculePdb, pr->staticMoleculePSF, pr->movingMoleculePSF, pr->staticMoleculeMol2, pr->movingMoleculeMol2, pr->prmFile, pr->aprmFile, pr->rtfFile);	
	pr->hbFilter->initializeFilter();
#else
	printf("LibMol was not found. Could not apply hygrogen bond filter\n");
#endif
}


int dockingMain( PARAMS_IN *pr, bool scoreUntransformed )
{
	int interpFuncExtent;
	double scale;
	if ( !computeGridParameters( pr, &interpFuncExtent, &scale ) ) return -1;

	int numFreq = pr->numFreq;
	double numFreq3 = numFreq * numFreq * numFreq;

	if ( pr->numberOfPositions == -1 ) pr->numberOfPositions = numFreq3;

	if ( scoreUntransformed )
	{
		for ( int i = 0; i < 9; i++ )
			pr->rotations[ i ] = ( i % 4 ) ? 0 : 1;

		pr->numberOfRotations = 1;
		pr->numThreads = 1;

		pr->numberOfPositions = numFreq3;
	}

	//#ifdef DEBUG
	printf("\n\ndata in PR in dockingMain\n");
	printf("performDocking: %d\n", pr->performDocking);
	printf("numThreads: %d\n", pr->numThreads);
	printf("breakDownScores: %d\n", pr->breakDownScores);
	printf("numberOfPositions: %d\n", pr->numberOfPositions);
	printf("gridSize: %d\n", pr->gridSize);
	printf("numFreq: %d\n", pr->numFreq);
	//  printf("interpFuncExtent: %d\n", pr->interpFuncExtent);
	printf("numCentersA: %i\n", pr->numCentersA);
	printf("numCentersB: %i\n", pr->numCentersB);
	printf("numberOfRotations: %d\n", pr->numberOfRotations);

	printf("distanceCutoff: %lf\n", pr->distanceCutoff);
	printf("alpha: %lf\n", pr->alpha);
	printf("blobbiness: %lf\n", pr->blobbiness);
	printf("skinSkinWeight: %lf\n", pr->skinSkinWeight);
	printf("coreCoreWeight: %lf\n", pr->coreCoreWeight);
	printf("skinCoreWeight: %lf\n", pr->skinCoreWeight);
	printf("realSCWeight: %lf\n", pr->realSCWeight);
	printf("imaginarySCWeight: %lf\n", pr->imaginarySCWeight);
	printf("elecScale: %lf\n", pr->elecScale);
	printf("scoreScaleUpFactor: %lf\n", pr->scoreScaleUpFactor);

	printf("bandwidth: %lf\n", pr->bandwidth);
	printf("gradFactor: %lf\n", pr->gradFactor);

	printf("outputFilename: %s\n", pr->outputFilename);
	printf("staticMoleculePdb: %s\n", pr->staticMoleculePdb);
	printf("staticMoleculeF2d: %s\n", pr->staticMoleculeF2d);
	printf("movingMoleculePdb: %s\n", pr->movingMoleculePdb);
	printf("movingMoleculeF2d: %s\n", pr->movingMoleculeF2d);
	printf("typeA: %s\n", pr->typeA);
	printf("typeB: %s\n", pr->typeB);

	printf("chargesA: %p %f %f\n", pr->chargesA, pr->chargesA[0],
			pr->chargesA[pr->numCentersA-1]);
	printf("radiiA: %p %f %f\n", pr->radiiA, pr->radiiA[0],
			pr->radiiA[pr->numCentersA-1]);
	printf("chargesB: %p %f %f\n", pr->chargesB, pr->chargesB[0],
			pr->chargesB[pr->numCentersB-1]);
	printf("radiiB: %p %f %f\n", pr->radiiB, pr->radiiB[0],
			pr->radiiB[pr->numCentersB-1]);
	printf("rotations: %p\n", pr->rotations);

	printf("xkAOrig: %p\n", pr->xkAOrig);
	printf("xkAOrig: %f %f\n", pr->xkAOrig[0], pr->xkAOrig[pr->numCentersA-1]);
	printf("ykAOrig: %p\n", pr->ykAOrig);
	printf("xkAOrig: %f %f\n", pr->ykAOrig[0], pr->ykAOrig[pr->numCentersA-1]);
	printf("zkAOrig: %p\n", pr->zkAOrig);
	printf("xkAOrig: %f %f\n", pr->zkAOrig[0], pr->zkAOrig[pr->numCentersA-1]);
	printf("atName: %s\n", pr->atNamesA[pr->numCentersA-1]);
	printf("resType: %s\n", pr->resTypesA[pr->numCentersA-1]);
	printf("resNum: %d\n", pr->resNumsA[pr->numCentersA-1]);

	printf("xkBOrig: %p\n", pr->xkBOrig);
	printf("xkBOrig: %f %f\n", pr->xkBOrig[0], pr->xkBOrig[pr->numCentersB-1]);
	printf("ykBOrig: %p\n", pr->ykBOrig);
	printf("xkBOrig: %f %f\n", pr->ykBOrig[0], pr->ykBOrig[pr->numCentersB-1]);
	printf("zkBOrig: %p\n", pr->zkBOrig);
	printf("xkBOrig: %f %f\n", pr->zkBOrig[0], pr->zkBOrig[pr->numCentersB-1]);
	printf("atName: %s\n", pr->atNamesB[pr->numCentersB-1]);
	printf("resType: %s\n", pr->resTypesB[pr->numCentersB-1]);
	printf("resNum: %d\n", pr->resNumsB[pr->numCentersB-1]);

	printf("nbRMSDAtoms: %d\n", pr->nbRMSDAtoms);
	if (pr->nbRMSDAtoms>0) {
		printf("atNums: %d %d\n", pr->atNums[0], pr->atNums[pr->nbRMSDAtoms-1]);
		printf("xRef: %f %f\n", pr->xRef[0], pr->xRef[pr->nbRMSDAtoms-1]);
		printf("yRef: %f %f\n", pr->yRef[0], pr->yRef[pr->nbRMSDAtoms-1]);
		printf("zRef: %f %f\n", pr->zRef[0], pr->zRef[pr->nbRMSDAtoms-1]);
	}
	//  printf("computevdw: %d\n", pr->computevdw);
	//  printf("vdwSmoothWidth: %lf\n", pr->vdwSmoothWidth);

	printf("control %d\n", pr->control);
	printf("numFreq %d\n", numFreq);
	//#endif

	int numThreads = pr->numThreads;
	char *outputFilename = pr->outputFilename;
	double blobbiness = pr->blobbiness;
	bool smoothSkin = (bool)pr->smoothSkin;
	bool rotateVolume = (bool) pr->rotateVolume;
	bool dockVolume = (bool) pr->dockVolume;
	bool singleLayerLigandSkin = ( bool ) pr->singleLayerLigandSkin;
	double distanceCutoff = pr->distanceCutoff;
	bool performDocking = (bool)pr->performDocking;
	bool breakDownScores = (bool)pr->breakDownScores;
	int numberOfPositions = pr->numberOfPositions;
	int gridSize = pr->gridSize;
	double gridSpacing = pr->gridSpacing;
	//  int numFreq = pr->numFreq;
	//  int interpFuncExtent = pr->interpFuncExtent;

	bool useSparseFFT = ( bool ) pr->useSparseFFT;
	bool narrowBand = ( bool ) pr->narrowBand;

	double alpha = pr->alpha;

	float*   = pr->rotations;
	int numberOfRotations = pr->numberOfRotations;

	double skinSkinWeight = pr->skinSkinWeight;
	double coreCoreWeight = pr->coreCoreWeight;
	double skinCoreWeight = pr->skinCoreWeight;
	double realSCWeight = pr->realSCWeight;
	double imaginarySCWeight = pr->imaginarySCWeight;

	double elecScale = pr->elecScale;
	double elecRadiusInGrids = pr->elecRadiusInGrids;

	double hydrophobicityWeight = pr->hydrophobicityWeight;
	double hydroPhobicPhobicWeight = pr->hydroPhobicPhobicWeight;
	double hydroPhobicPhilicWeight = pr->hydroPhobicPhilicWeight;
	double hydroPhilicPhilicWeight = pr->hydroPhilicPhilicWeight;

	double hydroRadExt = pr->hydroRadExt;

	bool twoWayHydrophobicity = pr->twoWayHydrophobicity;

	double staticMolHydroDistCutoff = pr->staticMolHydroDistCutoff;

	double simpleShapeWeight = pr->simpleShapeWeight;
	double simpleChargeWeight = pr->simpleChargeWeight;

	double simpleRadExt = pr->simpleRadExt;

	double hbondWeight = pr->hbondWeight;
	double hbondDistanceCutoff = pr->hbondDistanceCutoff;

	double clashWeight = pr->clashWeight;

	double pseudoGsolWeight = pr->pseudoGsolWeight;

	double dispersionWeight = pr->dispersionWeight;

	double bandwidth = pr->bandwidth;
	double gradFactor = pr->gradFactor;

	bool curvatureWeightedStaticMol = pr->curvatureWeightedStaticMol;
	bool curvatureWeightedMovingMol = pr->curvatureWeightedMovingMol;
	double curvatureWeightingRadius = pr->curvatureWeightingRadius;
	bool spreadReceptorSkin = pr->spreadReceptorSkin;
	bool randomRotate = pr->randomRotate;
	Matrix randRot = pr->initRot;

	double clusterTransRad = pr->clusterTransRad;
	int clusterTransSize = pr->clusterTransSize;
	double clusterRotRad = pr->clusterRotRad;
	int peaksPerRotation = pr->peaksPerRotation;

	char *outputFileName = outputFilename;

	double real_magnitude = sqrt( skinSkinWeight );
	double imag_magnitude = sqrt( coreCoreWeight );

	double hydrophobicity_real_magnitude = sqrt( hydroPhobicPhobicWeight );
	double hydrophobicity_imag_magnitude = sqrt( hydroPhilicPhilicWeight );

	// results
	TopValues *localTopValues[ numThreads ];
	TopValues *globalTopValues = 0;

	TopValues *funnel = 0;
	int *peakList;

	double mainStartTime = getTime();

	for ( int i = 0; i < numThreads; i++ )
		localTopValues[ i ] = 0;

	//  double numFreq3 = numFreq*numFreq*numFreq;
	FFTW_complex* = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 ) ;

	FFTW_complex* centerElecFrequenciesA = 0;
	FFTW_complex* smallElectrostaticsKernel = 0;

	if ( elecScale != 0 )
	{
		centerElecFrequenciesA = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );
		smallElectrostaticsKernel = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );
	}

	FFTW_complex* centerHbondFrequenciesA = 0;

	if ( hbondWeight != 0 )
		centerHbondFrequenciesA = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );

	FFTW_complex* centerHydrophobicityFrequenciesA = 0;
	FFTW_complex* centerHydrophobicityTwoFrequenciesA = 0;

	if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
	{
		centerHydrophobicityFrequenciesA = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );
		if ( twoWayHydrophobicity )
			centerHydrophobicityTwoFrequenciesA = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );
	}

	FFTW_complex* centerSimpleComplementarityFrequenciesA = 0;

	if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		centerSimpleComplementarityFrequenciesA = (FFTW_complex*)FFTW_malloc( sizeof(FFTW_complex)*numFreq3 );

	FFTW_complex* ourMoreFrequenciesA;

	FFTW_complex* centerFrequenciesB[ numThreads ];
	FFTW_complex* centerFrequenciesProduct[ numThreads ];
	FFTW_complex* sparseProfile[ numThreads ];

	FFTW_complex* sparseShapeProfile[ numThreads ];

	FFTW_complex* centerElecFrequenciesB[ numThreads ];
	FFTW_complex* centerFrequenciesElecProduct[ numThreads ];
	FFTW_complex* sparseElecProfile[ numThreads ];

	FFTW_complex* sparseHbondProfile[ numThreads ];
	FFTW_complex* sparseHydrophobicityProfile[ numThreads ];
	FFTW_complex* sparseHydrophobicityTwoProfile[ numThreads ];
	FFTW_complex* sparseSimpleComplementarityProfile[ numThreads ];

	FFTW_complex* freqHat[ numThreads ];

	FFTW_complex* ourMoreFrequencies[ numThreads ];
	FFTW_complex* ourMoreFrequenciesOut[ numThreads ];

	FFTW_complex* elecGridB[ numThreads ];

	FFTW_plan freqPlan[ numThreads ];
	sparse3DFFT_plan sparseFreqPlanBackward;
	FFTW_plan elecFreqPlan[ numThreads ];

	FFTW_plan freqHatPlan[ numThreads ];

	FFTW_plan moreFreqPlan[ numThreads ];
	sparse3DFFT_plan sparseFreqPlanForward;
	FFTW_plan moreFreqPlanA;

	FFTW_plan moreElecFreqPlan[ numThreads ];
	FFTW_plan moreElecFreqPlanA;

	double *sortedPeaks[ numThreads ];

	// static molecule
	double *xkAOrig = pr->xkAOrig;
	double *ykAOrig = pr->ykAOrig;
	double *zkAOrig =pr->zkAOrig;
	int numCentersA = pr->numCentersA;
	char *typeA = pr->typeA;
	char *hbondTypeA = pr->hbondTypeA;
	float *chargesA = pr->chargesA;
	float *hydrophobicityA = pr->hydrophobicityA;
	float *radiiA = pr->radiiA;
	int _numCentersA = 0;

	while ( typeA[ _numCentersA ] == 'I' ) 
		_numCentersA++;

	// moving molecule
	double *xkBOrig = pr->xkBOrig;
	double *ykBOrig = pr->ykBOrig;
	double *zkBOrig = pr->zkBOrig;
	int numCentersB = pr->numCentersB;
	double *xkB[ numThreads ], *ykB[ numThreads ], *zkB[ numThreads ];
	float *rkB[ numThreads ];
	char *typeB = pr->typeB;
	char *hbondTypeB = pr->hbondTypeB;
	float *chargesB = pr->chargesB;
	float *hydrophobicityB = pr->hydrophobicityB;
	float *radiiB = pr->radiiB;

	NONZERO_GRIDCELLS *gridBCells = 0;
	int numNonzeroGridBCells = 0;

	NONZERO_GRIDCELLS *gridBCells_01 = 0;
	int numNonzeroGridBCells_01 = 0;

	NONZERO_GRIDCELLS *gridBCells_10 = 0;
	int numNonzeroGridBCells_10 = 0;

	NONZERO_GRIDCELLS *gridBCells_11 = 0;
	int numNonzeroGridBCells_11 = 0;

	NONZERO_GRIDCELLS *elecGridBCells = 0;
	int numNonzeroElecGridBCells = 0;

	NONZERO_GRIDCELLS *hbondGridBCells = 0;
	int numNonzeroHbondGridBCells = 0;

	NONZERO_GRIDCELLS *hydrophobicityGridBCells = 0;
	int numNonzeroHydrophobicityGridBCells = 0;

	NONZERO_GRIDCELLS *hydrophobicityTwoGridBCells = 0;
	int numNonzeroHydrophobicityTwoGridBCells = 0;

	NONZERO_GRIDCELLS *simpleComplementarityGridBCells = 0;
	int numNonzeroSimpleComplementarityGridBCells = 0;

	fastLJ *ljFilter;

	if ( pr->applyVdWFilter )
	{
		initLJFilter( pr, &ljFilter );
		pr->ljFilter = ljFilter;
	}

	clashFilter *cFilter;
	clashFilter* fFilter;

	if ( pr->applyClashFilter )
	{
		initClashFilter( pr, &cFilter );
		pr->cFilter = cFilter;
	}

	if(pr->applyForbiddenVolumeFilter) 
	{
		initForbiddenVolumeFilter( pr, &fFilter );
		std::cout<<"Forbidden volume filter enabled.\n";	
		pr->fFilter = fFilter;
	}

	if ( pr->applyPseudoGsolFilter ) initPseudoGsolFilter( pr );

	if ( pr->applyDispersionFilter ) initDispersionFilter( pr );

#ifdef LIBMOL_FOUND
	if ( pr->applyHbondFilter ) initHBondFilter( pr );
#else
	printf("LibMol was not found. Could not apply hygrogen bond filter\n");
#endif

	//#ifdef DEBUG
	printf("\n\nLocal variable in main\n");
	printf("numThreads: %d\n", numThreads);
	printf("performDocking: %d\n", performDocking);
	printf("blobbiness: %lf\n", blobbiness);
	printf("distanceCutoff: %lf\n", distanceCutoff);
	printf("outputFilename: %s\n", outputFilename);
	printf("breakDownScores: %d\n", breakDownScores);
	printf("numberOfPositions: %d\n", numberOfPositions);
	printf("gridSize: %d\n", gridSize);
	printf("gridSpacing: %lf\n", gridSpacing);
	printf("numFreq: %d\n", numFreq);
	printf("interpFuncExtent: %d\n", interpFuncExtent);

	printf("bandwidth: %lf\n", bandwidth);
	printf("gradFactor: %lf\n", gradFactor);

	printf("numberOfRotations: %d\n",numberOfRotations );
	printf("skinSkinWeight: %f\n", skinSkinWeight);
	printf("coreCoreWeight: %f\n", coreCoreWeight);
	printf("skinCoreWeight: %f\n", skinCoreWeight);
	printf("realSCWeight: %f\n", realSCWeight);
	printf("imaginarySCWeight: %f\n", imaginarySCWeight);
	printf("elecScale: %f\n",elecScale );
	printf("hbondWeight: %f\n",hbondWeight );
	printf("hydrophobicityWeight = %f\n", hydrophobicityWeight );
	printf("hydroPhobicPhobicWeight: %f\n", hydroPhobicPhobicWeight);
	printf("hydroPhobicPhilicWeight: %f\n", hydroPhobicPhilicWeight);
	printf("hydroPhilicPhilicWeight: %f\n", hydroPhilicPhilicWeight);
	printf("simpleShapeWeight: %f\n", simpleShapeWeight);
	printf("simpleChargeWeight: %f\n", simpleChargeWeight);
	printf("Molecule A\n");
	printf("   numCenters:%d", numCentersA);
	printf("   coords  0: %lf %lf %lf\n", xkAOrig[0], ykAOrig[0], zkAOrig[0]);
	printf("   coords -1: %lf %lf %lf\n", xkAOrig[numCentersA-1],
			ykAOrig[numCentersA-1], zkAOrig[numCentersA-1]);
	printf("   types: %s\n", typeA);
	printf("   charges: %f %f\n", chargesA[0], chargesA[numCentersA-1]);
	printf("   radii: %f %f\n", radiiA[0], radiiA[numCentersA-1]);

	printf("Molecule B\n");
	printf("   numCenters:%d", numCentersB);
	printf("   coords  0: %lf %lf %lf\n", xkBOrig[0], ykBOrig[0], zkBOrig[0]);
	printf("   coords -1: %lf %lf %lf\n", xkBOrig[numCentersB-1],
			ykBOrig[numCentersB-1], zkBOrig[numCentersB-1]);
	printf("   types: %s\n", typeB);
	printf("   charges: %f %f\n", chargesB[0], chargesB[numCentersB-1]);
	printf("   radii: %f %f\n", radiiB[0], radiiB[numCentersB-1]);
	//#endif


	double *xkA = NULL, *ykA = NULL, *zkA = NULL;
	float *rkA = NULL;

	if ( numCentersA > 0 )
	{
		xkA = new double[ numCentersA ];
		ykA = new double[ numCentersA ];
		zkA = new double[ numCentersA ];
		rkA = new float[ numCentersA ];
	}

	for ( int i = 0; i < numThreads; i++ )
	{
		if ( numCentersB > 0  )
		{
			xkB[ i ] = new double[ numCentersB ];
			ykB[ i ] = new double[ numCentersB ];
			zkB[ i ] = new double[ numCentersB ];
			rkB[ i ] = new float[ numCentersB ];
		}
		else
		{
			xkB[ i ] = ykB[ i ] = zkB[ i ] = NULL;
			rkB[ i ] = NULL;
		}
	}


	FFTW_complex *fkA = 0, 					  *fkB = 0;
	FFTW_complex *fkAElec = 0, 				  *fkBElec = 0;
	FFTW_complex *fkAHbond = 0, 				  *fkBHbond = 0;
	FFTW_complex *fkAHydrophobicity = 0, 		  *fkBHydrophobicity = 0;
	FFTW_complex *fkAHydrophobicityTwo = 0, 	  *fkBHydrophobicityTwo = 0;
	FFTW_complex *fkASimpleComplementarity = 0, *fkBSimpleComplementarity = 0;

	FFTW_complex* gridA = 0;

	gridA = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );

	globalTopValues = new TopValues( numberOfPositions, numFreq );
	if ( clusterRotRad > 0 )
	{
		funnel = new TopValues( 2 * numThreads, numFreq );
		peakList = new int[ 2 * numberOfPositions ];
	}

	SmoothingFunction* smoothingFunction[ numThreads ];

	double n = ( int ) ( alpha * numFreq );
	int alphaM = ( int ) n;

	int *validOutputMap = ( int * ) FFTW_malloc( sizeof( int ) * numFreq3 );

	for ( int i = 0; i < numThreads; i++ )
	{
		centerFrequenciesB[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		centerFrequenciesProduct[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		sparseProfile[ i ] = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		sparseShapeProfile[ i ] = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );

		freqHat[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * gridSize * 4 );

		ourMoreFrequencies[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * alphaM * alphaM * alphaM );

		if ( !useSparseFFT ) freqPlan[ i ] = FFTW_plan_dft_3d( numFreq, numFreq, numFreq,
				centerFrequenciesProduct[ i ], sparseProfile[ i ],
				FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE );

		freqHatPlan[ i ] = FFTW_plan_dft_1d( gridSize * 4, freqHat[ i ], freqHat[ i ], FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );

		if ( rotateVolume )
			ourMoreFrequenciesOut[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * alphaM * alphaM * alphaM );
		else
			ourMoreFrequenciesOut[ i ] = ourMoreFrequencies[ i ];

		if ( !useSparseFFT ) moreFreqPlan[ i ] = FFTW_plan_dft_3d( alphaM, alphaM, alphaM,
				ourMoreFrequencies[ i ], ourMoreFrequenciesOut[ i ],
				FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );

		if ( elecScale != 0 )
		{
			centerElecFrequenciesB[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
			centerFrequenciesElecProduct[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
			sparseElecProfile[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );

			elecFreqPlan[ i ] = FFTW_plan_dft_c2r_3d( numFreq, numFreq, numFreq,
					centerFrequenciesElecProduct[ i ],
					( FFTW_DATA_TYPE * ) sparseElecProfile[ i ],
					OPT_FFTW_SEARCH_TYPE );

			if ( rotateVolume )
				elecGridB[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * alphaM * alphaM * alphaM );
			else elecGridB[ i ] = centerElecFrequenciesB[ i ];

			moreElecFreqPlan[ i ] = FFTW_plan_dft_r2c_3d( alphaM, alphaM, alphaM,
					( FFTW_DATA_TYPE * ) elecGridB[ i ], centerElecFrequenciesB[ i ],
					OPT_FFTW_SEARCH_TYPE );
		}
		else
			centerElecFrequenciesB[ i ] = centerFrequenciesElecProduct[ i ] = sparseElecProfile[ i ] = elecGridB[ i ] = NULL;

		if ( hbondWeight != 0 )
			sparseHbondProfile[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		else
			sparseHbondProfile[ i ] = NULL;

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
		{
			sparseHydrophobicityProfile[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3);

			if ( twoWayHydrophobicity )
				sparseHydrophobicityTwoProfile[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
			else
				sparseHydrophobicityTwoProfile[ i ] = NULL;
		}
		else
			sparseHydrophobicityProfile[ i ] = sparseHydrophobicityTwoProfile[ i ] = NULL;


		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
			sparseSimpleComplementarityProfile[ i ] = ( FFTW_complex * ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		else
			sparseSimpleComplementarityProfile[ i ] = NULL;

		localTopValues[ i ] = new TopValues( numberOfPositions, numFreq );

		if ( smoothSkin ) smoothingFunction[ i ] = new Gaussian( alpha, interpFuncExtent, (int)(alpha*numFreq), gridSize );
		else smoothingFunction[ i ] = NULL;

		if ( ( clusterTransRad > 0 ) || ( peaksPerRotation < numFreq3 ) || ( ( clusterRotRad > 0 ) && ( i == 0 ) ) )
		{
			sortedPeaks[ i ] = ( double * ) malloc( sizeof( double ) * numFreq3 );
			for ( int j = 0; j < numFreq3; j++ )
				sortedPeaks[ i ][ j ] = 0;
		}
		else sortedPeaks[ i ] = NULL;
	}


	if ( !useSparseFFT ) sparseFreqPlanForward = sparseFreqPlanBackward = NULL;

	moreFreqPlanA = FFTW_plan_dft_3d( alphaM, alphaM, alphaM,
			ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ],
			FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE );

	if ( elecScale != 0 ) moreElecFreqPlanA = FFTW_plan_dft_r2c_3d( alphaM, alphaM, alphaM,
			( FFTW_DATA_TYPE * ) elecGridB[ 0 ],
			centerElecFrequenciesB[ 0 ],
			OPT_FFTW_SEARCH_TYPE );

	FFTW_complex *fkA_01 = 0, *fkA_10 = 0, *fkA_11 = 0;
	FFTW_complex *fkB_01 = 0, *fkB_10 = 0, *fkB_11 = 0;

	FFTW_complex *centerFrequenciesA_01 = 0, *centerFrequenciesA_10 = 0, *centerFrequenciesA_11 = 0;

	FFTW_complex* sparseProfile_01[ numThreads ];
	FFTW_complex* sparseProfile_10[ numThreads ];
	FFTW_complex* sparseProfile_11[ numThreads ];

	if ( breakDownScores )
	{
		centerFrequenciesA_01 = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		centerFrequenciesA_10 = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		centerFrequenciesA_11 = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );

		for ( int i = 0; i < numThreads; i++ )
		{
			sparseProfile_01[ i ] = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
			sparseProfile_10[ i ] = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
			sparseProfile_11[ i ] = ( FFTW_complex* ) FFTW_malloc( sizeof( FFTW_complex ) * numFreq3 );
		}
	}
	else
	{
		centerFrequenciesA_01 = centerFrequenciesA_10 = centerFrequenciesA_11 = NULL;

		for ( int i = 0; i < numThreads; i++ )
			sparseProfile_01[ i ] = sparseProfile_10[ i ] = sparseProfile_11[ i ] = NULL;
	}

	// report peak
	// FPRINTF CUIDADO COM MPI
	FILE* fpOpt = fopen( outputFilename, "w");
	printInputParamters( pr, fpOpt );

	// read molecules
	if ( breakDownScores )
	{
		if ( !dockVolume )
		{
			if ( !build_fks( typeA, numCentersA, xkAOrig, ykAOrig, zkAOrig, radiiA, chargesA, hbondTypeA, hydrophobicityA,
						&fkA, &fkA_01, &fkA_10, &fkA_11,
						elecScale, &fkAElec, hbondWeight, &fkAHbond,
						hydrophobicityWeight, hydroPhobicPhobicWeight, hydroPhobicPhilicWeight, hydroPhilicPhilicWeight,
						staticMolHydroDistCutoff, &fkAHydrophobicity, twoWayHydrophobicity, &fkAHydrophobicityTwo,
						simpleShapeWeight, simpleChargeWeight, &fkASimpleComplementarity,
						real_magnitude, imag_magnitude,
						true, singleLayerLigandSkin, curvatureWeightedStaticMol, curvatureWeightingRadius,
						bandwidth, gradFactor, pr ) ) 
				return -1;

			if ( !build_fks( typeB, numCentersB, xkBOrig, ykBOrig, zkBOrig, radiiB, chargesB, hbondTypeB, hydrophobicityB,
						&fkB, &fkB_01, &fkB_10, &fkB_11,
						elecScale, &fkBElec, hbondWeight, &fkBHbond,
						hydrophobicityWeight, hydroPhobicPhobicWeight, hydroPhobicPhilicWeight, hydroPhilicPhilicWeight,
						staticMolHydroDistCutoff, &fkBHydrophobicity, twoWayHydrophobicity, &fkBHydrophobicityTwo,
						simpleShapeWeight, simpleChargeWeight, &fkBSimpleComplementarity,
						real_magnitude, imag_magnitude,
						false, singleLayerLigandSkin, curvatureWeightedMovingMol, curvatureWeightingRadius,
						bandwidth, gradFactor, pr ) ) 
				return -1;
		}
	}
	else
	{
		if ( !dockVolume )
		{
			if ( !build_fks( typeA, numCentersA, xkAOrig, ykAOrig, zkAOrig, radiiA, chargesA, hbondTypeA, hydrophobicityA,
						&fkA, elecScale, &fkAElec, hbondWeight, &fkAHbond,
						hydrophobicityWeight, hydroPhobicPhobicWeight, hydroPhobicPhilicWeight, hydroPhilicPhilicWeight,
						staticMolHydroDistCutoff, &fkAHydrophobicity, twoWayHydrophobicity, &fkAHydrophobicityTwo,
						simpleShapeWeight, simpleChargeWeight, &fkASimpleComplementarity,
						real_magnitude, imag_magnitude,
						true, singleLayerLigandSkin, curvatureWeightedStaticMol, curvatureWeightingRadius,
						bandwidth, gradFactor, pr ) ) 
				return -1;

			if ( !build_fks( typeB, numCentersB, xkBOrig, ykBOrig, zkBOrig, radiiB, chargesB, hbondTypeB, hydrophobicityB,
						&fkB, elecScale, &fkBElec, hbondWeight, &fkBHbond,
						hydrophobicityWeight, hydroPhobicPhobicWeight, hydroPhobicPhilicWeight, hydroPhilicPhilicWeight,
						staticMolHydroDistCutoff, &fkBHydrophobicity, twoWayHydrophobicity, &fkBHydrophobicityTwo,
						simpleShapeWeight, simpleChargeWeight, &fkBSimpleComplementarity,
						real_magnitude, imag_magnitude,
						false, singleLayerLigandSkin, curvatureWeightedMovingMol, curvatureWeightingRadius,
						bandwidth, gradFactor, pr ) ) 
				return -1;
		}
	}


	float translate_A[ 3 ];
	float translate_B[ 3 ];

	/*
	 *	A partir daqui parece fazer coisas
	 */

	if ( elecScale != 0 )
	{
		// get electrostatics kernel using the right scale and grid sizes
		// get elec kernel FFT defined by distance weighted dielectric.
		// this function does not return Fourier coefficients, but FFT of the function sampled on the grid. This is in
		// accordance with the fastsum paper by nfft folks.
		if ( !computeElecKernel( smallElectrostaticsKernel, numFreq, scale,
								pr->elecKernelVoidRad,
								pr->elecKernelDistLow, pr->elecKernelValLow,
								pr->elecKernelDistHigh, pr->elecKernelValHigh ) ) 
			return -1;
	}

	INSPHERE_DATA isd;
	sparse3DFFT_nonzero_func sparsityFuncForward = inSphere;
	void *sparsityFuncForwardData = &isd;

	isd.cx = isd.cy = isd.cz = ( alphaM - 1 ) / 2.0;

	VALID_OUTPUT_DATA voutd;
	sparse3DFFT_nonzero_func sparsityFuncBackward = outputValid;
	void *sparsityFuncBackwardData = &voutd;

	voutd.n = numFreq;
	voutd.validOutputMap = validOutputMap;

	if ( dockVolume )
	{
		double skinSkinWeightSqrt = sqrt( pr->skinSkinWeight );
		double coreCoreWeightSqrt = sqrt( pr->coreCoreWeight );

		double xCenter, yCenter, zCenter;

		if ( !readShapeCompGrid( sparseProfile[ 0 ], &xCenter, &yCenter, &zCenter, pr->movingMoleculeSCReRaw, pr->movingMoleculeSCImRaw ) ) 
			return -1;

		for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
		{
			centerFrequenciesB[ 0 ][ i ][ 0 ] = sparseProfile[ 0 ][ i ][ 0 ] * skinSkinWeightSqrt;
			centerFrequenciesB[ 0 ][ i ][ 1 ] = sparseProfile[ 0 ][ i ][ 1 ] * coreCoreWeightSqrt;
		}

		translate_B[ 0 ] = -xCenter;
		translate_B[ 1 ] = -yCenter;
		translate_B[ 2 ] = -zCenter;

		if ( breakDownScores ) memcpy( sparseProfile_01[ 0 ], sparseProfile[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

		if ( !readShapeCompGrid( sparseProfile[ 0 ], &xCenter, &yCenter, &zCenter, pr->staticMoleculeSCReRaw, pr->staticMoleculeSCImRaw ) ) 
			return -1;

		for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
		{
			gridA[ i ][ 0 ] = ourMoreFrequencies[ 0 ][ i ][ 0 ] = sparseProfile[ 0 ][ i ][ 0 ] * skinSkinWeightSqrt;
			gridA[ i ][ 1 ] = ourMoreFrequencies[ 0 ][ i ][ 1 ] = sparseProfile[ 0 ][ i ][ 1 ] * coreCoreWeightSqrt;
		}

		translate_A[ 0 ] = -xCenter;
		translate_A[ 1 ] = -yCenter;
		translate_A[ 2 ] = -zCenter;

		if ( !getCenterFrequencies(  numCentersA, xkA, ykA, zkA, rkA, typeA, fkA, blobbiness, alpha, numFreq, interpFuncExtent,
					centerFrequenciesA, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ],
					moreFreqPlanA, NULL, true, spreadReceptorSkin ) ) 
			return -1;

		if ( breakDownScores )
		{
			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				ourMoreFrequencies[ 0 ][ i ][ 0 ] = 0;
				ourMoreFrequencies[ 0 ][ i ][ 1 ] = sparseProfile[ 0 ][ i ][ 1 ];
			}

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_01, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_01, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, true, spreadReceptorSkin ) ) 
				return -1;

			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				ourMoreFrequencies[ 0 ][ i ][ 0 ] = sparseProfile[ 0 ][ i ][ 0 ];
				ourMoreFrequencies[ 0 ][ i ][ 1 ] = 0;
			}

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_10, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_10, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, true, spreadReceptorSkin ) ) 
				return -1;


			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				ourMoreFrequencies[ 0 ][ i ][ 0 ] = sparseProfile[ 0 ][ i ][ 0 ];
				ourMoreFrequencies[ 0 ][ i ][ 1 ] = sparseProfile[ 0 ][ i ][ 1 ];
			}

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_11, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_11, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, true, spreadReceptorSkin ) ) 
				return -1;
		}

		if ( elecScale != 0 )
		{
			if ( !readElecGrid( ( FFTW_DATA_TYPE * ) elecGridB[ 0 ], &xCenter, &yCenter, &zCenter, pr->staticMoleculeElecReRaw ) ) 
				return -1;

			if ( !getCenterElecFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkAElec, blobbiness, alpha, numFreq, elecRadiusInGrids,
						centerElecFrequenciesB[ 0 ], ( FFTW_DATA_TYPE * ) elecGridB[ 0 ], moreElecFreqPlanA, true, false ) ) 
				return -1;

			memcpy( centerElecFrequenciesA, centerElecFrequenciesB[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
		}

		numNonzeroGridBCells = collectNonzeroGridCells( &gridBCells, centerFrequenciesB[ 0 ], numFreq );

		double minRadB, maxRadB;

		findMovingMolMinMaxRadius( gridBCells, numNonzeroGridBCells, numFreq, &minRadB, &maxRadB );

		createValidOutputMap( sparseProfile[ 0 ], numFreq, ( int ) minRadB, ( int ) maxRadB, narrowBand, validOutputMap );

		if ( useSparseFFT )
		{
			isd.r2 = ( ceil( maxRadB ) + 2 ) * ( ceil( maxRadB ) + 2 );

			sparseFreqPlanForward = sparse3DFFT_create_plan( alphaM, alphaM, alphaM, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE,
					SPARSE3DFFT_SPARSEINPUT,
					sparsityFuncForward, sparsityFuncForwardData, ourMoreFrequencies[ 0 ], ourMoreFrequenciesOut[ 0 ] );

			sparseFreqPlanBackward = sparse3DFFT_create_plan( alphaM, alphaM, alphaM, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE,
					SPARSE3DFFT_SPARSEOUTPUT,
					sparsityFuncBackward, sparsityFuncBackwardData, ourMoreFrequencies[ 0 ], ourMoreFrequenciesOut[ 0 ] );
		}

		flipGrid( gridBCells, numNonzeroGridBCells );

		if ( breakDownScores )
		{
			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				centerFrequenciesB[ 0 ][ i ][ 0 ] = 0;
				centerFrequenciesB[ 0 ][ i ][ 1 ] = sparseProfile_01[ 0 ][ i ][ 1 ];
			}

			numNonzeroGridBCells_01 = collectNonzeroGridCells( &gridBCells_01, centerFrequenciesB[ 0 ], numFreq );

			flipGrid( gridBCells_01, numNonzeroGridBCells_01 );

			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				centerFrequenciesB[ 0 ][ i ][ 0 ] = sparseProfile_01[ 0 ][ i ][ 0 ];
				centerFrequenciesB[ 0 ][ i ][ 1 ] = 0;
			}

			numNonzeroGridBCells_10 = collectNonzeroGridCells( &gridBCells_10, centerFrequenciesB[ 0 ], numFreq );

			flipGrid( gridBCells_10, numNonzeroGridBCells_10 );

			for ( int i = 0; i < numFreq * numFreq * numFreq; i++ )
			{
				centerFrequenciesB[ 0 ][ i ][ 0 ] = sparseProfile_01[ 0 ][ i ][ 0 ];
				centerFrequenciesB[ 0 ][ i ][ 1 ] = sparseProfile_01[ 0 ][ i ][ 1 ];
			}

			numNonzeroGridBCells_11 = collectNonzeroGridCells( &gridBCells_11, centerFrequenciesB[ 0 ], numFreq );

			flipGrid( gridBCells_11, numNonzeroGridBCells_11 );
		}

		if ( randomRotate )
		{
			rotateGrid( centerFrequenciesB[ 0 ], numFreq, gridBCells, numNonzeroGridBCells, randRot, true, false );
			free( gridBCells );
			numNonzeroGridBCells = collectNonzeroGridCells( &gridBCells, centerFrequenciesB[ 0 ], numFreq );

			if ( breakDownScores )
			{
				rotateGrid( centerFrequenciesB[ 0 ], numFreq, gridBCells_01, numNonzeroGridBCells_01, randRot, true, false );
				free( gridBCells_01 );
				numNonzeroGridBCells_01 = collectNonzeroGridCells( &gridBCells_01, centerFrequenciesB[ 0 ], numFreq );

				rotateGrid( centerFrequenciesB[ 0 ], numFreq, gridBCells_10, numNonzeroGridBCells_10, randRot, true, false );
				free( gridBCells_10 );
				numNonzeroGridBCells_10 = collectNonzeroGridCells( &gridBCells_10, centerFrequenciesB[ 0 ], numFreq );

				rotateGrid( centerFrequenciesB[ 0 ], numFreq, gridBCells_11, numNonzeroGridBCells_11, randRot, true, false );
				free( gridBCells_11 );
				numNonzeroGridBCells_11 = collectNonzeroGridCells( &gridBCells_11, centerFrequenciesB[ 0 ], numFreq );
			}
		}

		printf( "\nnumNonzeroGridBCells = %d\n\n", numNonzeroGridBCells );

		if ( elecScale != 0 )
		{
			if ( !readElecGrid( ( FFTW_DATA_TYPE * ) centerElecFrequenciesB[ 0 ], &xCenter, &yCenter, &zCenter, pr->movingMoleculeElecReRaw ) ) 
				return -1;

			numNonzeroElecGridBCells = collectNonzeroElecGridCells( &elecGridBCells, ( FFTW_DATA_TYPE * ) centerElecFrequenciesB[ 0 ], numFreq );

			flipGrid( elecGridBCells, numNonzeroElecGridBCells );

			if ( randomRotate )
			{
				rotateElecGrid( ( FFTW_DATA_TYPE * ) centerElecFrequenciesB[ 0 ], numFreq, elecGridBCells, numNonzeroElecGridBCells, randRot, true, false );
				free( elecGridBCells );
				numNonzeroElecGridBCells = collectNonzeroElecGridCells( &elecGridBCells, ( FFTW_DATA_TYPE * ) centerElecFrequenciesB[ 0 ], numFreq );
			}

			printf( "numNonzeroElecGridBCells = %d\n\n", numNonzeroElecGridBCells );
		}
	}
	else
	{
		double xTrans = 0, yTrans = 0, zTrans = 0;

		// center A
		if( !center( xkAOrig, ykAOrig, zkAOrig, numCentersA, &xTrans, &yTrans, &zTrans ) ) 
			return -1;

		translate_A[ 0 ] = xTrans;
		translate_A[ 1 ] = yTrans;
		translate_A[ 2 ] = zTrans;

		// transform A
		if ( !transformAndNormalize( xkAOrig, ykAOrig, zkAOrig, radiiA, xkA, ykA, zkA, rkA, numCentersA,
					1, 0, 0, 0, 1, 0, 0, 0, 1, 1, scale ) ) 
			return -1;

		gridding( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA, blobbiness, numFreq, interpFuncExtent,
				smoothSkin, smoothingFunction[ 0 ], gridA /*sparseProfile[ 0 ]*/, spreadReceptorSkin );

		memcpy( ourMoreFrequencies[ 0 ], gridA /*sparseProfile[ 0 ]*/, numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

		// get frequencies of A.
		if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA, blobbiness, alpha, numFreq, interpFuncExtent,
					centerFrequenciesA, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, true, spreadReceptorSkin ) ) 
			return -1;

		if ( breakDownScores )
		{
			memcpy( ourMoreFrequencies[ 0 ], sparseProfile[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_01, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_01, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false, spreadReceptorSkin ) ) 
				return -1;

			memcpy( ourMoreFrequencies[ 0 ], sparseProfile[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_10, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_10, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false, spreadReceptorSkin ) ) 
				return -1;

			memcpy( ourMoreFrequencies[ 0 ], sparseProfile[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

			if ( !getCenterFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkA_11, blobbiness, alpha, numFreq, interpFuncExtent,
						centerFrequenciesA_11, smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false, spreadReceptorSkin ) ) 
				return -1;
		}

		if ( elecScale != 0 )
		{
			if ( !getCenterElecFrequencies( numCentersA, xkA, ykA, zkA, rkA, typeA, fkAElec, blobbiness, alpha, numFreq, elecRadiusInGrids,
						centerElecFrequenciesB[ 0 ], ( FFTW_DATA_TYPE * ) elecGridB[ 0 ], moreElecFreqPlanA, false, false ) ) 
				return -1;

			memcpy( centerElecFrequenciesA, centerElecFrequenciesB[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
		}

		if ( hbondWeight != 0 )
		{
			if ( !getCenterHbondFrequencies( numCentersA, xkA, ykA, zkA, rkA, hbondDistanceCutoff * scale * ( numFreq - 1 ), fkAHbond, blobbiness, alpha, numFreq,
						ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false, false ) ) 
				return -1;

			memcpy( centerHbondFrequenciesA, ourMoreFrequencies[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
		}

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
		{
			if ( !getCenterHydrophobicityFrequencies( numCentersA, xkA, ykA, zkA, rkA, fkAHydrophobicity, blobbiness, alpha, numFreq,
						( hydrophobicityWeight == 0 ) ? hydroRadExt : 0.0, ( bool ) ( hydrophobicityWeight == 0 ),
						ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false ) ) 
				return -1;

			memcpy( centerHydrophobicityFrequenciesA, ourMoreFrequencies[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );

			if ( twoWayHydrophobicity )
			{
				if ( !getCenterHydrophobicityFrequencies( numCentersA, xkA, ykA, zkA, rkA, fkAHydrophobicityTwo, blobbiness, alpha, numFreq,
							hydroRadExt, false, ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ],
							ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false ) ) 
					return -1;

				memcpy( centerHydrophobicityTwoFrequenciesA, ourMoreFrequencies[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
			}
		}

		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
		{
			if ( !getCenterSimpleComplementarityFrequencies( numCentersA, xkA, ykA, zkA, rkA, fkASimpleComplementarity, blobbiness, alpha, numFreq, simpleRadExt,
						ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], ourMoreFrequencies[ 0 ], moreFreqPlanA, NULL, false ) ) 
				return -1;

			memcpy( centerSimpleComplementarityFrequenciesA, ourMoreFrequencies[ 0 ], numFreq * numFreq * numFreq * sizeof( FFTW_complex ) );
		}

		// center B
		if ( !center( xkBOrig, ykBOrig, zkBOrig, numCentersB, &xTrans, &yTrans, &zTrans ) ) 
			return -1;

		translate_B[ 0 ] = xTrans;
		translate_B[ 1 ] = yTrans;
		translate_B[ 2 ] = zTrans;

		if ( randomRotate )
		{
			if ( !rotateAboutOrigin( xkBOrig, ykBOrig, zkBOrig, numCentersB, randRot ) ) 
				return -1;
		}

		if ( !transformAndNormalize( xkBOrig, ykBOrig, zkBOrig, radiiB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], numCentersB,
					1, 0, 0, 0, 1, 0, 0, 0, 1, 2, scale ) ) 
			return -1;

		gridding( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], typeB, fkB, blobbiness, numFreq, interpFuncExtent,
				smoothSkin, smoothingFunction[ 0 ], ourMoreFrequencies[ 0 ], false );

		numNonzeroGridBCells = collectNonzeroGridCells( &gridBCells, ourMoreFrequencies[ 0 ], numFreq );

		if ( rotateVolume )
		{
			if ( breakDownScores )
			{
				gridding( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], typeB, fkB_01, blobbiness, numFreq, interpFuncExtent,
						smoothSkin, smoothingFunction[ 0 ], centerFrequenciesB[ 0 ], false );

				numNonzeroGridBCells_01 = collectNonzeroGridCells( &gridBCells_01, centerFrequenciesB[ 0 ], numFreq );

				gridding( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], typeB, fkB_10, blobbiness, numFreq, interpFuncExtent,
						smoothSkin, smoothingFunction[ 0 ], centerFrequenciesB[ 0 ], false );

				numNonzeroGridBCells_10 = collectNonzeroGridCells( &gridBCells_10, centerFrequenciesB[ 0 ], numFreq );

				gridding( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], typeB, fkB_11, blobbiness, numFreq, interpFuncExtent,
						smoothSkin, smoothingFunction[ 0 ], centerFrequenciesB[ 0 ], false );

				numNonzeroGridBCells_11 = collectNonzeroGridCells( &gridBCells_11, centerFrequenciesB[ 0 ], numFreq );
			}

			if ( elecScale != 0 )
			{
				griddingElec( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], typeB, fkBElec,
						blobbiness, numFreq, elecRadiusInGrids, ( FFTW_DATA_TYPE * ) elecGridB[ 0 ], false, true );

				numNonzeroElecGridBCells = collectNonzeroElecGridCells( &elecGridBCells, ( FFTW_DATA_TYPE * ) elecGridB[ 0 ], numFreq );

				printf( "numNonzeroElecGridBCells = %d\n\n", numNonzeroElecGridBCells );
			}


			if ( hbondWeight != 0 )
			{
				griddingHbond( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], hbondDistanceCutoff * scale * ( numFreq - 1 ), fkBHbond,
						blobbiness, numFreq, ourMoreFrequencies[ 0 ], true );

				numNonzeroHbondGridBCells = collectNonzeroHbondGridCells( &hbondGridBCells, ourMoreFrequencies[ 0 ], numFreq );

				printf( "numNonzeroHbondGridBCells = %d\n\n", numNonzeroHbondGridBCells );
			}


			if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
			{
				griddingHydrophobicity( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], fkBHydrophobicity,
						blobbiness, numFreq, ourMoreFrequencies[ 0 ], hydroRadExt, ( bool ) ( hydrophobicityWeight == 0 ) );

				numNonzeroHydrophobicityGridBCells = collectNonzeroHydrophobicityGridCells( &hydrophobicityGridBCells, ourMoreFrequencies[ 0 ], numFreq );

				printf( "numNonzeroHydrophobicityGridBCells = %d\n\n", numNonzeroHydrophobicityGridBCells );

				if ( twoWayHydrophobicity )
				{
					griddingHydrophobicity( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], fkBHydrophobicityTwo,
							blobbiness, numFreq, ourMoreFrequencies[ 0 ], 0, false );

					numNonzeroHydrophobicityTwoGridBCells = collectNonzeroHydrophobicityGridCells( &hydrophobicityTwoGridBCells, ourMoreFrequencies[ 0 ], numFreq );

					printf( "numNonzeroHydrophobicityTwoGridBCells = %d\n\n", numNonzeroHydrophobicityTwoGridBCells );
				}
			}


			if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
			{
				griddingSimpleComplementarity( numCentersB, xkB[ 0 ], ykB[ 0 ], zkB[ 0 ], rkB[ 0 ], fkBSimpleComplementarity,
						blobbiness, numFreq, ourMoreFrequencies[ 0 ], simpleRadExt );

				numNonzeroSimpleComplementarityGridBCells = collectNonzeroSimpleComplementarityGridCells( &simpleComplementarityGridBCells, ourMoreFrequencies[ 0 ], numFreq );

				printf( "numNonzeroSimpleComplementarityGridBCells = %d\n\n", numNonzeroSimpleComplementarityGridBCells );
			}
		}


		double minRadB, maxRadB;

		findMovingMolMinMaxRadius( numCentersB, xkBOrig, zkBOrig, zkBOrig, radiiB, typeB, numFreq, scale, &minRadB, &maxRadB );

		createValidOutputMap( gridA /*sparseProfile[ 0 ]*/, numFreq, ( int ) minRadB, ( int ) maxRadB, narrowBand, validOutputMap );

		if ( rotateVolume || useSparseFFT )
		{
			if ( useSparseFFT )
			{
				isd.r2 = ( ceil( maxRadB ) + 2 ) * ( ceil( maxRadB ) + 2 );

				sparseFreqPlanForward = sparse3DFFT_create_plan( alphaM, alphaM, alphaM, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE,
						SPARSE3DFFT_SPARSEINPUT,
						sparsityFuncForward, sparsityFuncForwardData, ourMoreFrequencies[ 0 ], ourMoreFrequenciesOut[ 0 ] );

				sparseFreqPlanBackward = sparse3DFFT_create_plan( alphaM, alphaM, alphaM, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE,
						SPARSE3DFFT_SPARSEOUTPUT,
						sparsityFuncBackward, sparsityFuncBackwardData, ourMoreFrequencies[ 0 ], ourMoreFrequenciesOut[ 0 ] );
			}

			printf( "\nnumNonzeroGridBCells = %d\n\n", numNonzeroGridBCells );
		}
	}



	PARAMS prT[ numThreads ];
	pthread_t p[ numThreads ];

	/*
	 *	Faz realmente coisas
	 */

	initRotationServer( numberOfRotations, numThreads );

	double functionScaleFactor = pow( numFreq * alpha, 6 ) / pr->scoreScaleUpFactor;

	for ( int i = 0; i < numThreads; i++ )
	{
		prT[ i ].threadID = i + 1;
		prT[ i ].confID = 0;
		prT[ i ].rotations = rotations;
		prT[ i ].numberOfRotations = numberOfRotations;
		prT[ i ].numberOfPositions = numberOfPositions;
		prT[ i ].xkBOrig = xkBOrig;
		prT[ i ].ykBOrig = ykBOrig;
		prT[ i ].zkBOrig = zkBOrig;
		prT[ i ].radiiB = radiiB;
		prT[ i ].xkB = xkB[ i ];
		prT[ i ].ykB = ykB[ i ];
		prT[ i ].zkB = zkB[ i ];
		prT[ i ].rkB = rkB[ i ];
		prT[ i ].numCentersB = numCentersB;
		prT[ i ].gridSize = gridSize;
		prT[ i ].numFreq = numFreq;
		prT[ i ].interpFuncExtent = interpFuncExtent;
		prT[ i ].alpha = alpha;
		prT[ i ].blobbiness = blobbiness;
		prT[ i ].skinSkinWeight = skinSkinWeight;
		prT[ i ].coreCoreWeight = coreCoreWeight;
		prT[ i ].skinCoreWeight = skinCoreWeight;
		prT[ i ].realSCWeight = realSCWeight;
		prT[ i ].imaginarySCWeight = imaginarySCWeight;
		prT[ i ].elecScale = elecScale;
		prT[ i ].elecRadiusInGrids = elecRadiusInGrids;
		prT[ i ].hbondWeight = hbondWeight;
		prT[ i ].hbondDistanceCutoff = hbondDistanceCutoff;
		prT[ i ].hydrophobicityWeight = hydrophobicityWeight;
		prT[ i ].hydroPhobicPhobicWeight = hydroPhobicPhobicWeight;
		prT[ i ].hydroPhobicPhilicWeight = hydroPhobicPhilicWeight;
		prT[ i ].hydroPhilicPhilicWeight = hydroPhilicPhilicWeight;
		prT[ i ].simpleShapeWeight = simpleShapeWeight;
		prT[ i ].simpleChargeWeight = simpleChargeWeight;
		prT[ i ].scaleA = scale;
		prT[ i ].translate_A = translate_A;
		prT[ i ].scaleB = scale;
		prT[ i ].translate_B = translate_B;
		prT[ i ].centerFrequenciesA = centerFrequenciesA;
		prT[ i ].centerElecFrequenciesA = centerElecFrequenciesA;
		prT[ i ].centerHbondFrequenciesA = centerHbondFrequenciesA;
		prT[ i ].centerHydrophobicityFrequenciesA = centerHydrophobicityFrequenciesA;
		prT[ i ].centerHydrophobicityTwoFrequenciesA = centerHydrophobicityTwoFrequenciesA;
		prT[ i ].centerSimpleComplementarityFrequenciesA = centerSimpleComplementarityFrequenciesA;

		prT[ i ].gridA = gridA;
		prT[ i ].gridB = ourMoreFrequencies[ i ];

		prT[ i ].centerFrequenciesB = centerFrequenciesB[ i ];
		prT[ i ].centerElecFrequenciesB = centerElecFrequenciesB[ i ];

		prT[ i ].rotateVolume = rotateVolume;
		prT[ i ].numNonzeroGridBCells = numNonzeroGridBCells;
		prT[ i ].gridBCells = gridBCells;

		prT[ i ].numNonzeroGridBCells_01 = numNonzeroGridBCells_01;
		prT[ i ].gridBCells_01 = gridBCells_01;
		prT[ i ].numNonzeroGridBCells_10 = numNonzeroGridBCells_10;
		prT[ i ].gridBCells_10 = gridBCells_10;
		prT[ i ].numNonzeroGridBCells_11 = numNonzeroGridBCells_11;
		prT[ i ].gridBCells_11 = gridBCells_11;

		prT[ i ].elecGridB = elecGridB[ i ];
		prT[ i ].numNonzeroElecGridBCells = numNonzeroElecGridBCells;
		prT[ i ].elecGridBCells = elecGridBCells;

		prT[ i ].numNonzeroHbondGridBCells = numNonzeroHbondGridBCells;
		prT[ i ].hbondGridBCells = hbondGridBCells;

		prT[ i ].numNonzeroHydrophobicityGridBCells = numNonzeroHydrophobicityGridBCells;
		prT[ i ].hydrophobicityGridBCells = hydrophobicityGridBCells;

		prT[ i ].numNonzeroHydrophobicityTwoGridBCells = numNonzeroHydrophobicityTwoGridBCells;
		prT[ i ].hydrophobicityTwoGridBCells = hydrophobicityTwoGridBCells;

		prT[ i ].numNonzeroSimpleComplementarityGridBCells = numNonzeroSimpleComplementarityGridBCells;
		prT[ i ].simpleComplementarityGridBCells = simpleComplementarityGridBCells;

		prT[ i ].validOutputMap = validOutputMap;

		//            prT[ i ].cFilter = cFilter;

		prT[ i ].centerFrequenciesProduct = centerFrequenciesProduct[ i ];
		prT[ i ].centerFrequenciesElecProduct = centerFrequenciesElecProduct[ i ];
		prT[ i ].sparseProfile = sparseProfile[ i ];
		prT[ i ].sparseShapeProfile = sparseShapeProfile[ i ];
		prT[ i ].sparseElecProfile = sparseElecProfile[ i ];
		prT[ i ].sparseHbondProfile = sparseHbondProfile[ i ];
		prT[ i ].sparseHydrophobicityProfile = sparseHydrophobicityProfile[ i ];
		prT[ i ].sparseHydrophobicityTwoProfile = sparseHydrophobicityTwoProfile[ i ];
		prT[ i ].sparseSimpleComplementarityProfile = sparseSimpleComplementarityProfile[ i ];
		prT[ i ].freqHat = freqHat[ i ];
		prT[ i ].freqPlan = freqPlan[ i ];
		prT[ i ].sparseFreqPlanBackward = sparseFreqPlanBackward;
		prT[ i ].elecFreqPlan = elecFreqPlan[ i ];

		prT[ i ].freqHatPlan = freqHatPlan[ i ];

		prT[ i ].ourMoreFrequencies = ourMoreFrequencies[ i ];
		prT[ i ].ourMoreFrequenciesOut = ourMoreFrequenciesOut[ i ];

		prT[ i ].moreFreqPlan = moreFreqPlan[ i ];
		prT[ i ].sparseFreqPlanForward = sparseFreqPlanForward;
		prT[ i ].moreElecFreqPlan = moreElecFreqPlan[ i ];
		prT[ i ].fkB = fkB;
		prT[ i ].fkBElec = fkBElec;
		prT[ i ].fkBHbond = fkBHbond;
		prT[ i ].fkBHydrophobicity = fkBHydrophobicity;
		prT[ i ].fkBHydrophobicityTwo = fkBHydrophobicityTwo;
		prT[ i ].fkBSimpleComplementarity = fkBSimpleComplementarity;
		prT[ i ].smallElectrostaticsKernel = smallElectrostaticsKernel;
		prT[ i ].smoothSkin = smoothSkin;
		prT[ i ].smoothingFunction = smoothingFunction[ i ];
		prT[ i ].localTopValues = localTopValues[ i ];
		prT[ i ].sortedPeaks = sortedPeaks[ i ];
		
		if ( clusterTransRad > 0 ) 
			prT[ i ].clustPG = new PG( 10.0, - numFreq * gridSpacing, 5.0 );
		else 
			prT[ i ].clustPG = NULL;
		
		double _gridFactor = ( double ) gridSize / ( double ) localTopValues[ i ]->getGridSize();
		
		prT[ i ].functionScaleFactor = functionScaleFactor; //pow( numFreq * alpha, 6 ) / pr->scoreScaleUpFactor;
		prT[ i ].gridFactor = _gridFactor;
		prT[ i ].pri = pr;
		prT[ i ].randRot = randRot;
		prT[ i ].breakDownScores = breakDownScores;

		if ( breakDownScores )
		{
			prT[ i ].centerFrequenciesA_01 = centerFrequenciesA_01;
			prT[ i ].centerFrequenciesA_10 = centerFrequenciesA_10;
			prT[ i ].centerFrequenciesA_11 = centerFrequenciesA_11;

			prT[ i ].fkB_01 = fkB_01;
			prT[ i ].fkB_10 = fkB_10;
			prT[ i ].fkB_11 = fkB_11;

			prT[ i ].sparseProfile_01 = sparseProfile_01[ i ];
			prT[ i ].sparseProfile_10 = sparseProfile_10[ i ];
			prT[ i ].sparseProfile_11 = sparseProfile_11[ i ];
		}

		pthread_create( &p[ i ], NULL, startApplyRotationsThread, ( void * ) &prT[ i ] );
	}

	// Threads finish execution; synchronize point
	for ( int i = 0; i < numThreads; i++ )
		pthread_join( p[ i ], NULL );

	if ( clusterTransRad > 0 )
	{
		for ( int i = 0; i < numThreads; i++ )
			delete prT[ i ].clustPG;
	}

	if ( !scoreUntransformed )
	{
		if ( clusterRotRad > 0  )
		{
			double cosTheta = cos( ( clusterRotRad / 180.0 ) * M_PI );
			double v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc, elec, hbond, hydrop, vdwp, scomp, pgsol, pgsolh, dispe, x, y, z;
			int nclashes, r, f, c, k;
			double vMax;

			for ( int i = 0; i < numThreads; i++ )
			{
				int n = localTopValues[ i ]->getCurrentNumberOfPositions( );

				while ( n-- )
				{
					localTopValues[ i ]->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					globalTopValues->updateTopValues( v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, c );
				}

				n = globalTopValues->getCurrentNumberOfPositions( );

				if ( n > 0 )
				{
					if ( i == 0 ) 
						vMax = v;
					else 
						if ( v > vMax ) 
							vMax = v;
				}

				while ( n-- )
				{
					globalTopValues->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					localTopValues[ i ]->updateTopValues( v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, c );
				}
			}

			for ( int i = 0, j = 0; i < numThreads; i++ )
			{
				int n = localTopValues[ i ]->getCurrentNumberOfPositions( );

				while ( n-- )
				{
					localTopValues[ i ]->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					globalTopValues->updateTopValues( v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, c );
				}

				n = globalTopValues->getCurrentNumberOfPositions( );

				while ( n-- )
				{
					globalTopValues->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					localTopValues[ i ]->updateTopValues( vMax - v + 1, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, c );
				}
			}

			for ( int i = 0; i < numThreads; i++ )
				if ( localTopValues[ i ]->getCurrentNumberOfPositions( ) > 0 )
				{
					localTopValues[ i ]->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					funnel->updateTopValues( v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, i );
				}

			int *sPeaks = ( int * ) malloc( sizeof( int ) * numFreq3 );

			for ( int i = 0; i < numFreq3; i++ )
				sPeaks[ i ] = -1;

			for ( int i = 0; i < 2 * numberOfPositions; i++ )
				peakList[ i ] = -1;

			int numInserted = 0;
			while ( ( numInserted < numberOfPositions ) && ( funnel->getCurrentNumberOfPositions( ) > 0 ) )
			{
				funnel->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
						&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
						&x, &y, &z, &r, &f, &k );

				if ( !closePeaksExist( z, y, x, r, rotations, sPeaks, peakList, numFreq, gridSpacing, clusterTransRad, cosTheta ) )
				{
					int cc = ( z * numFreq + y ) * numFreq + x;

					peakList[ 2 * numInserted ] = r;
					peakList[ 2 * numInserted + 1 ] = sPeaks[ cc ];
					sPeaks[ cc ] = 2 * numInserted;

					globalTopValues->updateTopValues( vMax - v + 1, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, 0 );
					numInserted++;
				}

				if ( localTopValues[ k ]->getCurrentNumberOfPositions( ) > 0 )
				{
					localTopValues[ k ]->extractMin( &v, &rv, &rv_ss, &rv_cc, &rv_sc, &iv, &iv_ss, &iv_cc, &iv_sc,
							&elec, &hbond, &hydrop, &vdwp, &nclashes, &scomp, &pgsol, &pgsolh, &dispe,
							&x, &y, &z, &r, &f, &c );
					funnel->updateTopValues( v, rv, rv_ss, rv_cc, rv_sc, iv, iv_ss, iv_cc, iv_sc,
							elec, hbond, hydrop, vdwp, nclashes, scomp, pgsol, pgsolh, dispe,
							x, y, z, r, f, k );
				}

			}

			free( sPeaks );
		}
		else
		{
			TopValues *globalTopValuesT = new TopValues( numThreads * numberOfPositions, numFreq );

			for ( int i = 0; i < numThreads; i++ )
			{
				printf("# \n# PROCESSING THREAD = %d\n# \n", i + 1 );
				fflush( stdout );

				int n = localTopValues[ i ]->getCurrentNumberOfPositions( );
				ValuePosition3D sol;

				while ( n-- )
				{
					localTopValues[ i ]->extractMin( sol );

					sol.m_Value = - sol.m_Value;

					globalTopValuesT->updateTopValues( sol );
				}
			}

			filterPoses( &prT[ 0 ], globalTopValuesT, globalTopValues, numFreq, scale, translate_A, translate_B,
					rotations, functionScaleFactor, randRot );

			if ( clusterTransRad > 0 )
			{
				double *hash3d = ( double * ) malloc( sizeof( double ) * numFreq3 );

				for ( int c = 0; c < numFreq3; c++ )
					hash3d[ c ] = 0;

				int n = globalTopValuesT->getCurrentNumberOfPositions( );
				ValuePosition3D sol;

				while ( n-- )
				{
					globalTopValuesT->extractMin( sol );

					double v = sol.m_Value;
					double x = sol.m_Translation[ 0 ], y = sol.m_Translation[ 1 ], z = sol.m_Translation[ 2 ];

					int clusterPenalty = 0;

					int cVal1 = getClusterValue( x, y, z, hash3d, pr->numFreq, gridSpacing, clusterTransRad * 1, -v );
					int cVal2 = getClusterValue( x, y, z, hash3d, pr->numFreq, gridSpacing, clusterTransRad * 2, -v );
					int cVal3 = getClusterValue( x, y, z, hash3d, pr->numFreq, gridSpacing, clusterTransRad * 3, -v );

					if ( cVal1 >= clusterTransSize * 1 + 0 ) clusterPenalty = 8;
					else if ( cVal2 >= clusterTransSize * 2 + 1 ) clusterPenalty = 5;
					else if ( cVal3 >= clusterTransSize * 3 + 3 ) clusterPenalty = 1;

					sol.m_clusterPenalty = clusterPenalty;
#ifdef RERANK_DEBUG
					sol.m_ConformationIndex = clusterPenalty;  
#endif
					sol.m_Value = ( - sol.m_Value ) * ( 1 - 0.10 * clusterPenalty );

					globalTopValues->updateTopValues( sol );

					if ( clusterPenalty == 0 ) markCluster( x, y, z, hash3d, pr->numFreq, gridSpacing, 0, -v );
				}

				free( hash3d );
			}
			else
			{
				int n = globalTopValuesT->getCurrentNumberOfPositions( );
				ValuePosition3D sol;

				while ( n-- )
				{
					globalTopValuesT->extractMin( sol );

					sol.m_Value = - sol.m_Value;

					globalTopValues->updateTopValues( sol );
				}

				fflush( stdout );
			}

			delete globalTopValuesT;

			if ( pr->rerank )
			{
				rerankPoses( &prT[ 0 ], globalTopValues );
			}
		}
	}

	// Free some of the memory allocated for the FFTs
	FFTW_free( fkB );
	if ( elecScale != 0 ) 
		FFTW_free( fkBElec );
	if ( hbondWeight != 0 ) 
		FFTW_free( fkBHbond );
	if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
	{
		FFTW_free( fkBHydrophobicity );
		if ( twoWayHydrophobicity ) 
			FFTW_free( fkBHydrophobicityTwo );
	}
	if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) ) 
		FFTW_free( fkBSimpleComplementarity );

	if ( numCentersB > 0  )
		for ( int i = 0; i < numThreads; i++ )
		{
			delete [ ] xkB[ i ];
			delete [ ] ykB[ i ];
			delete [ ] zkB[ i ];
		}

	if ( breakDownScores )
	{
		FFTW_free( fkB_01 );
		FFTW_free( fkB_10 );
		FFTW_free( fkB_11 );
	}

	if ( !scoreUntransformed )
	{
		printf("# \n# Computation Time = %f sec\n# ", getTime() - mainStartTime);

		fprintf( fpOpt, (char *)"# computation time = %f sec\n# \n# ", getTime() - mainStartTime);

		fprintf( fpOpt, (char *)"# center of the 2nd protein:\n# \n# " );
		fprintf( fpOpt, (char *)"     conformation 0: %f %f %f\n# ", -translate_B[0], -translate_B[1], -translate_B[2] );
		printf( "\n# " );
	}


	prT[ 0 ].confID = 0;
	prT[ 0 ].rotations = rotations;
	prT[ 0 ].numFreq = numFreq;
	prT[ 0 ].scaleB = scale;
	prT[ 0 ].translate_A = translate_A;
	prT[ 0 ].translate_B = translate_B;
	prT[ 0 ].functionScaleFactor = functionScaleFactor;
	prT[ 0 ].localTopValues = globalTopValues;
	prT[ 0 ].gridFactor = ( double ) gridSize / ( double ) globalTopValues->getGridSize( );
	prT[ 0 ].pri = pr;
	prT[ 0 ].randRot = randRot;


	if ( !scoreUntransformed )
	{
		prT[ 0 ].localTopValues = globalTopValues;
		printIntermediateStats( fpOpt, 20, 5, &prT[ 0 ] );
	}
	else
	{
		prT[ 0 ].localTopValues = localTopValues[ 0 ];
		printUntransformedScore( fpOpt, &prT[ 0 ] );
	}

	fflush( fpOpt );

	// Free allocated memory
	FFTW_free( centerFrequenciesA );
	FFTW_destroy_plan( moreFreqPlanA );
	FFTW_free( fkA );
	FFTW_free( gridA );
	if ( rotateVolume )
	{
		free( gridBCells );
		if ( breakDownScores )
		{
			free( gridBCells_01 );
			free( gridBCells_10 );
			free( gridBCells_11 );
		}
	}

	if ( elecScale != 0 )
	{
		FFTW_free( centerElecFrequenciesA );
		FFTW_destroy_plan( moreElecFreqPlanA );
		FFTW_free( fkAElec );
		if ( rotateVolume ) free( elecGridBCells );
	}

	if ( hbondWeight != 0 )
	{
		FFTW_free( centerHbondFrequenciesA );
		FFTW_free( fkAHbond );
		if ( rotateVolume ) free( hbondGridBCells );
	}

	if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
	{
		FFTW_free( centerHydrophobicityFrequenciesA );
		FFTW_free( fkAHydrophobicity );
		if ( rotateVolume ) free( hydrophobicityGridBCells );

		if ( twoWayHydrophobicity )
		{
			FFTW_free( centerHydrophobicityTwoFrequenciesA );
			FFTW_free( fkAHydrophobicityTwo );
			if ( rotateVolume ) free( hydrophobicityTwoGridBCells );
		}
	}

	if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
	{
		FFTW_free( centerSimpleComplementarityFrequenciesA );
		FFTW_free( fkASimpleComplementarity );
		if ( rotateVolume ) free( simpleComplementarityGridBCells );
	}

	for ( int i = 0; i < numThreads; i++ )
	{
		FFTW_free( centerFrequenciesB[ i ] );
		FFTW_free( centerFrequenciesProduct[ i ] );
		FFTW_free( sparseProfile[ i ] );
		FFTW_free( sparseShapeProfile[ i ] );

		FFTW_free( freqHat[ i ] );

		FFTW_free( ourMoreFrequencies[ i ] );

		if ( rotateVolume )
		{
			FFTW_free( ourMoreFrequenciesOut[ i ] );
			if ( elecScale != 0 ) FFTW_free( elecGridB[ i ] );
		}

		FFTW_destroy_plan( freqHatPlan[ i ] );

		if ( !useSparseFFT )
		{
			FFTW_destroy_plan( freqPlan[ i ] );
			FFTW_destroy_plan( moreFreqPlan[ i ] );
		}

		if ( elecScale != 0 )
		{
			FFTW_free( centerElecFrequenciesB[ i ] );
			FFTW_free( centerFrequenciesElecProduct[ i ] );
			FFTW_free( sparseElecProfile[ i ] );

			FFTW_destroy_plan( elecFreqPlan[ i ] );
			FFTW_destroy_plan( moreElecFreqPlan[ i ] );
		}

		if ( hbondWeight != 0 )
			FFTW_free( sparseHbondProfile[ i ] );

		if ( ( hydrophobicityWeight != 0 ) || ( hydroPhobicPhobicWeight != 0 ) || ( hydroPhilicPhilicWeight != 0 ) || ( hydroPhobicPhilicWeight != 0 ) )
		{
			FFTW_free( sparseHydrophobicityProfile[ i ] );

			if ( twoWayHydrophobicity ) FFTW_free( sparseHydrophobicityTwoProfile[ i ] );
		}

		if ( ( simpleShapeWeight > 0 ) || ( simpleChargeWeight > 0 ) )
			FFTW_free( sparseSimpleComplementarityProfile[ i ] );

		delete localTopValues[ i ];
		if ( smoothSkin ) delete smoothingFunction[ i ];

		if ( ( clusterTransRad > 0 ) || ( peaksPerRotation < numFreq3 ) || ( ( clusterRotRad > 0 ) && ( i == 0 ) ) )
			free( sortedPeaks[ i ] );
	}

	if ( breakDownScores )
	{
		FFTW_free( fkA_01 );
		FFTW_free( fkA_10 );
		FFTW_free( fkA_11 );

		FFTW_free( centerFrequenciesA_01 );
		FFTW_free( centerFrequenciesA_10 );
		FFTW_free( centerFrequenciesA_11 );

		for ( int i = 0; i < numThreads; i++ )
		{
			FFTW_free( sparseProfile_01[ i ] );
			FFTW_free( sparseProfile_10[ i ] );
			FFTW_free( sparseProfile_11[ i ] );
		}
	}

	if ( useSparseFFT )
	{
		sparse3DFFT_destroy_plan( sparseFreqPlanForward );
		sparse3DFFT_destroy_plan( sparseFreqPlanBackward );
	}

	FFTW_free( validOutputMap );

	delete globalTopValues;

	if ( clusterRotRad > 0 )
	{
		delete funnel;
		delete peakList;
	}

	if ( pr->applyVdWFilter ) delete ljFilter;

	if ( pr->applyClashFilter ) delete cFilter;

	if ( pr->applyPseudoGsolFilter ) delete pr->pGsol;

	if ( pr->applyDispersionFilter ) delete pr->dispF;

	if ( !scoreUntransformed )
	{
		printf("# \n# Total Time = %f sec\n# \n# ", getTime( ) - mainStartTime);

		fprintf( fpOpt, (char *)"\n# total time = %f sec\n# ", getTime( ) - mainStartTime);
	}

	fclose( fpOpt );

	return 0;
}



int scoreUntransformed( PARAMS_IN *pr )
{
	return dockingMain( pr, true );
}


int dock( PARAMS_IN *pr )
{
	return dockingMain( pr, false );
}
