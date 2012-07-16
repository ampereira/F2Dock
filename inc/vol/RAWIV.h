/*
  Copyright 2011 The University of Texas at Austin

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


#ifndef __RAWIV_H__
#define __RAWIV_H__

#include "../fft-utils/fftw3.h"
#include "../fft-utils/fftwPrecision.h"

#define SWAP_64( a ) \
  { \
    unsigned char tmp[ 8 ]; \
    unsigned char *ch; \
    ch = ( unsigned char * )( a ); \
    tmp[ 0 ] = ch[ 0 ]; tmp[ 1 ] = ch[ 1 ]; tmp[ 2 ] = ch[ 2 ]; tmp[ 3 ] = ch[ 3 ]; \
    tmp[ 4 ] = ch[ 4 ]; tmp[ 5 ] = ch[ 5 ]; tmp[ 6 ] = ch[ 6 ]; tmp[ 7 ] = ch[ 7 ]; \
    ch[ 0 ] = tmp[ 7 ]; ch[ 1 ] = tmp[ 6 ]; ch[ 2 ] = tmp[ 5 ]; ch[ 3 ] = tmp[ 4 ]; \
    ch[ 4 ] = tmp[ 3 ]; ch[ 5 ] = tmp[ 2 ]; ch[ 6 ] = tmp[ 1 ]; ch[ 7 ] = tmp[ 0 ]; \
  }
#define SWAP_32( a ) \
  { \
    unsigned char tmp[ 4 ]; \
    unsigned char *ch; \
    ch = ( unsigned char * )( a ); \
    tmp[ 0 ] = ch[ 0 ]; tmp[ 1 ] = ch[ 1 ]; tmp[ 2 ] = ch[ 2 ]; tmp[ 3 ] = ch[ 3 ]; \
    ch[ 0 ] = tmp[ 3 ]; ch[ 1 ] = tmp[ 2 ]; ch[ 2 ] = tmp[ 1 ]; ch[ 3 ] = tmp[ 0 ]; \
  }

#ifdef FFTW_SINGLE_PRECISION
   #define SWAP_FFTW_DATA SWAP_32
#else
   #define SWAP_FFTW_DATA SWAP_64
#endif

struct RAWIVHeader
{
  float minExt[ 3 ];
  float maxExt[ 3 ];
  unsigned int numVertices;
  unsigned int numCells;
  unsigned int dim[ 3 ];
  float origin[ 3 ];
  float span[ 3 ];
}; 

static inline int bigEndian( void )
{
  long one = 1;
  return !( *( ( char * )( &one ) ) );
}

int writeGrid( FFTW_complex *scGrid, FFTW_DATA_TYPE *elecGrid, int n, double xCenter, double yCenter, double zCenter, double scale, 
               char *fileName, char *fileNameSCRe, char *fileNameSCIm, char *fileNameElecRe );
int readRAWIVHeader( int *xDim, int *yDim, int *zDim, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileName );
int readShapeCompGrid( FFTW_complex **scGrid, int *n, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileNameRe, char *fileNameIm );
int readShapeCompGrid( FFTW_complex *scGrid, double *xCenter, double *yCenter, double *zCenter, char *fileNameRe, char *fileNameIm );
int readElecGrid( FFTW_DATA_TYPE **elecGrid, int *n, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileNameRe );
int readElecGrid( FFTW_DATA_TYPE *elecGrid, double *xCenter, double *yCenter, double *zCenter, char *fileNameRe );

#endif
