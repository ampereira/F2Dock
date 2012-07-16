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


/* If you want to compile by itself do the following:
Please compile it as follows for double precision floats:
g++ rank-fftw.cpp -O2 -lfftw3 -o rank-fftw.exe
For single precision floats:
g++ rank-fftw.cpp -O2 -DFFTW_SINGLE_PRECISION -lfftw3f -o rank-fftw.exe
It will tell you about other options when you try to run the compiled
executable. The executable should be run as follows:
./rank-fftw.exe minGridSize maxGridSize
It will measure FFTW for all even grid sizes from minGridSize to
maxGridSize.
*/

#include "rank-fftw.h"
#include <time.h>

#ifdef _WIN32
    #include <sys/types.h>
    #include <sys/timeb.h>
#else
    #include <sys/time.h>
#endif

char errorMsg[ ][ 100 ] = { "No Error",
                            "Invalid Parameter",
                            "Memory Allocation Failed",
                            "Forward FFTW Failed",
                            "Backward FFTW Failed",
                            "Forward FFTW Wisdom Failed",
                            "Backward FFTW Wisdom Failed",
                            "Failed to Open File for Exporting Wisdom",
                            "Failed to Open File for Importing Wisdom" };


/*double getTime( void )
{
#ifdef _WIN32
    time_t ltime;
    _timeb tstruct;
    time( &ltime );
    _ftime( &tstruct );
    return ( double ) ( ltime + 1e-3*( tstruct.millitm ) );
#else
    struct timeval t;
    gettimeofday( &t, NULL );
    return ( double ) ( t.tv_sec + 1e-6 * t.tv_usec );
#endif
}*/


void printError( int errorID )
{
   fprintf( stderr, "\nError: %s!\n\n", errorMsg[ - errorID ] );
}


bool create3DFFTWisdom( FFTW_complex* dest, FFTW_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags, FILE *wisdomFP )
{
    if ( !dest || !src || dim1 < 1 || dim2 < 1 || dim3 < 1 ) return false;

    FFTW_plan p;

    p = FFTW_plan_dft_3d( dim1, dim2, dim3, src, dest, type, flags );

#ifdef SAVE_WISDOM_TO_FILE
    FFTW_export_wisdom_to_file( wisdomFP );
#endif

    FFTW_destroy_plan( p );

    return true;
}



bool compute3DFFT( FFTW_complex* dest, FFTW_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags, FILE *wisdomFP )

{
    if ( !dest || !src || dim1 < 1 || dim2 < 1 || dim3 < 1 ) return false;

    FFTW_plan p;

//    if ( !FFTW_import_wisdom_from_file( wisdomFP ) ) return false;

#ifdef SAVE_WISDOM_TO_FILE
    FFTW_import_wisdom_from_file( wisdomFP );
#endif

    p = FFTW_plan_dft_3d( dim1, dim2, dim3, src, dest, type, flags | FFTW_DESTROY_INPUT );

    FFTW_execute( p );

    FFTW_destroy_plan( p );

    return true;
}



void initializeData( FFTW_complex *data, int size )
{
    for ( int i = 0; i < size * size * size; i++ )
      {
       data[ i ][ 0 ] = ( FFTW_DATA_TYPE ) rand( ) / ( ( FFTW_DATA_TYPE ) ( RAND_MAX ) );
       data[ i ][ 1 ] = ( FFTW_DATA_TYPE ) rand( ) / ( ( FFTW_DATA_TYPE ) ( RAND_MAX ) );
      }
}



double timeFFTW( int size, bool inPlace, int iter, char *wisdomFile )

{
    if ( ( size < 1 ) || ( iter < 1 ) ) return INVALID_PARAMETER;

    FFTW_complex* src = ( FFTW_complex* ) FFTW_malloc( size * size * size * sizeof( FFTW_complex ) );

    if ( src == NULL ) return MEMORY_ALLOCATION_FAILED;

    FFTW_complex* dest;

    if ( inPlace ) dest = src;

    else dest = ( FFTW_complex* ) FFTW_malloc( size * size * size * sizeof( FFTW_complex ) );

    if ( dest == NULL )
      {
       FFTW_free( src );
       return MEMORY_ALLOCATION_FAILED;
      }

    initializeData( src, size );

    FILE *fp = NULL;

#ifdef SAVE_WISDOM_TO_FILE
    if ( ( fp = fopen( wisdomFile, "at" ) ) == NULL )
       {
        FFTW_free( src );
        if ( !inPlace ) FFTW_free( dest );
        return WISDOM_EXPORT_FILE_OPEN_FAILED;
       }

#endif

    if ( !create3DFFTWisdom( dest, src, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE, fp ) )
       {
        FFTW_free( src );
        if ( !inPlace ) FFTW_free( dest );
        return FFTW_FORWARD_WISDOM_FAILED;
       }

    if ( !create3DFFTWisdom( dest, src, size, size, size, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE, fp ) )
       {
        FFTW_free( src );
        if ( !inPlace ) FFTW_free( dest );
        return FFTW_BACKWARD_WISDOM_FAILED;
       }

#ifdef SAVE_WISDOM_TO_FILE
    fclose( fp );
#endif


#ifdef SAVE_WISDOM_TO_FILE
    if ( ( fp = fopen( wisdomFile, "rt" ) ) == NULL )
       {
        FFTW_free( src );
        if ( !inPlace ) FFTW_free( dest );
        return WISDOM_IMPORT_FILE_OPEN_FAILED;
       }

#endif

    double t = getTime( );

    for ( int i = 0; i < iter; i++ )
       {
        if ( !compute3DFFT( dest, src, size, size, size, FFTW_FORWARD, OPT_FFTW_SEARCH_TYPE, fp ) )
           {
            FFTW_free( src );
            if ( !inPlace ) FFTW_free( dest );
            return FFTW_FORWARD_FAILED;
           }

        if ( !compute3DFFT( dest, src, size, size, size, FFTW_BACKWARD, OPT_FFTW_SEARCH_TYPE, fp ) )
           {
            FFTW_free( src );
            if ( !inPlace ) FFTW_free( dest );
            return FFTW_BACKWARD_FAILED;
           }
       }

    t = getTime( ) - t;

#ifdef SAVE_WISDOM_TO_FILE
    fclose( fp );
#endif

    FFTW_free( src );

    if ( !inPlace ) FFTW_free( dest );

    return ( t / iter );
}



int rankFFTWSizesByTime( int minSize, int maxSize )
{
    int *size;
    double *runTime;

    size = ( int * ) malloc( maxSize * sizeof( int ) );
    runTime = ( double * ) malloc( maxSize * sizeof( double ) );

    if ( ( size == NULL ) || ( runTime == NULL ) )
      {
       if ( size == NULL ) free( size );
       if ( runTime == NULL ) free( runTime );
       return MEMORY_ALLOCATION_FAILED;
      }

    srand( time( NULL ) );

#ifdef FFTW_SINGLE_PRECISION
    printf( "# Measuring 3D FFTW performance on single precision floats\n" );
#else
    printf( "# Measuring 3D FFTW performance on double precision floats (define FFTW_SINGLE_PRECISION and link to fftw3f for single precision floats)\n" );
#endif


#if OPT_FFTW_SEARCH_TYPE == FFTW_MEASURE
    printf( "# Using flag FFTW_MEASURE (redefine OPT_FFTW_SEARCH_TYPE for other options: FFTW_ESTIMATE, FFTW_PATIENT and FFTW_EXHAUSTIVE)\n" );
#else
 #if OPT_FFTW_SEARCH_TYPE == FFTW_PATIENT
    printf( "# Using flag FFTW_PATIENT (redefine OPT_FFTW_SEARCH_TYPE for other options: FFTW_ESTIMATE, FFTW_MEASURE and FFTW_EXHAUSTIVE)\n" );
 #else
  #if OPT_FFTW_SEARCH_TYPE == FFTW_EXHAUSTIVE
    printf( "# Using flag FFTW_EXHAUSTIVE (redefine OPT_FFTW_SEARCH_TYPE for other options: FFTW_ESTIMATE, FFTW_MEASURE and FFTW_PATIENT)\n" );
  #else
    printf( "# Using flag FFTW_ESTIMATE (redefine OPT_FFTW_SEARCH_TYPE for other options: FFTW_MEASURE, FFTW_PATIENT and FFTW_EXHAUSTIVE)\n" );
  #endif
 #endif
#endif


#ifdef SAVE_WISDOM_TO_FILE
    printf( "# Saving wisdom to %s (filename can be changed by redefining WISDOM_FILE)\n#\n", ( char * ) WISDOM_FILE );
#else
    printf( "# Wisdom will not be saved (since SAVE_WISDOM_TO_FILE is not defined)\n#\n" );
#endif

    double t = getTime( );

    printf( "# size -> time ( sec )\n" );

    for ( int i = minSize; i <= maxSize; i += 2 )
      {
       int j = i >> 1;

       size[ j ] = i;
       runTime[ j ] = timeFFTW( i, IN_PLACE, MAX_ITER, ( char * ) WISDOM_FILE );

       if ( runTime[ j ] < 0 )
         {
          free( size );
          free( runTime );
          return ( ( int ) runTime[ j ] );
         }

       printf( "# %4d -> %12.6lf\n", i, runTime[ j ] );
       fflush( stdout );
      }

    t = getTime( ) - t;

    printf( "#\n# Time elapsed = %lf sec\n", t );
    fflush( stdout );

    for ( int i = ( maxSize >> 1 ); i > ( minSize >> 1 ); )
      {
       int j = i--;
       while ( ( i >= 1 ) && ( runTime[ i ] >= runTime[ j ] ) ) size[ i-- ] = -1;
      }

    printf( "#\n# List of efficient sizes...\n" );

    for ( int i = ( minSize >> 1 ); i <= ( maxSize >> 1 ); i++ )
      if ( size[ i ] != -1 ) printf( "%d\n", size[ i ] );

    fflush( stdout );

    free( size );
    free( runTime );

    return NO_ERROR;
}



//int computeEffGrid( PARAMS_IN *pr )
int computeEffGrid( int minSize, int maxSize )
{
//    int minSize = pr->minEffGridSize;
//    int maxSize = pr->maxEffGridSize;
    int errID = NO_ERROR;

    if ( errID == NO_ERROR ) errID = rankFFTWSizesByTime( minSize, maxSize );
    if ( errID != NO_ERROR ) printError( errID );

    return errID;
}

