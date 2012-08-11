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

#include <stdio.h>
#include <stdlib.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#endif

#include <string.h>

#include "RAWIV.h"

int writeRAWIV( FFTW_DATA_TYPE *vol, int n, double xCenter, double yCenter, double zCenter, double scale, char *fileName )
{
   FILE *fp;
   
   if ( ( fp = fopen( fileName, "wb" ) ) == NULL )
     {
      fprintf( stderr, "\nError: Failed to create file %s!\n\n", fileName );
      return 0;
     }
   
   RAWIVHeader header;
   float gridLength1D = 1.0 / scale;
   
   header.minExt[ 0 ] = xCenter - gridLength1D / 2.0;
   header.minExt[ 1 ] = yCenter - gridLength1D / 2.0;
   header.minExt[ 2 ] = zCenter - gridLength1D / 2.0;

   header.maxExt[ 0 ] = xCenter + gridLength1D / 2.0;
   header.maxExt[ 1 ] = yCenter + gridLength1D / 2.0;
   header.maxExt[ 2 ] = zCenter + gridLength1D / 2.0;
  
   header.numVertices = n * n * n;
   header.numCells = ( n - 1 ) * ( n - 1 ) * ( n - 1 );   
   
   for ( int i = 0; i < 3; i++ )
     {
      header.dim[ i ] = n;
      header.origin[ i ] = header.minExt[ i ]; 
      header.span[ i ] = gridLength1D / ( n - 1 );   
     } 
   
   if ( !bigEndian( ) )
     {
      for ( int i = 0; i < 3; i++ )
        {
         SWAP_32( &( header.minExt[ i ] ) );
         SWAP_32( &( header.maxExt[ i ] ) );         
         SWAP_32( &( header.dim[ i ] ) );         
         SWAP_32( &( header.origin[ i ] ) );                           
         SWAP_32( &( header.span[ i ] ) );         
        }
        
      SWAP_32( &( header.numVertices ) );  
      SWAP_32( &( header.numCells ) );        
     }
     
   if ( fwrite( &header, sizeof( header ), 1, fp ) != 1 )  
     {
      fprintf( stderr, "\nError: Failed to write header to file %s!\n\n", fileName );
      fclose( fp );      
      return 0;      
     }
   
   int rowLen = n * sizeof( FFTW_DATA_TYPE );
   unsigned char *buf = ( unsigned char * ) calloc( rowLen, 1 );
   
   if ( buf == NULL )
     {
      fprintf( stderr, "\nError: Failed to allocate temporary buffer space!\n\n" );
      fclose( fp );      
      return 0;      
     }

   FFTW_DATA_TYPE *volPtr = vol;
   
   for ( int k = 0; k < n; k++ )
     for ( int j = 0; j < n; j++ )
       {
        memcpy( buf, volPtr, rowLen );
        volPtr += n;
       
        if ( !bigEndian( ) )
          {
           for ( int i = 0; i < n; i++ )
              SWAP_FFTW_DATA( buf + i * sizeof( FFTW_DATA_TYPE ) );
          }

        if ( fwrite( buf, rowLen, 1, fp ) != 1 )  
          {
           fprintf( stderr, "\nError: Failed to write data to file %s!\n\n", fileName );
           free( buf );           
           fclose( fp );      
           return 0;      
          }
       }     
     
   free( buf );
   fclose( fp );
   
   return 1;
}

int writeGrid( FFTW_complex *scGrid, FFTW_DATA_TYPE *elecGrid, int n, double xCenter, double yCenter, double zCenter, double scale, 
               char *fileName, char *fileNameSCRe, char *fileNameSCIm, char *fileNameElecRe )
{
   if ( ( scGrid == NULL ) || ( elecGrid == NULL ) || ( n < 1 ) || ( scale <= 0 ) || ( fileName == NULL ) ) return 0;
   
   FFTW_DATA_TYPE *vol;
   int nCube = n * n * n;
   
   vol = ( FFTW_DATA_TYPE * ) malloc( 2 * nCube * sizeof( FFTW_DATA_TYPE ) );
   
   if ( vol == NULL ) return 0;
   
   for ( int i = 0; i < nCube; i++ )
     {
      vol[ i ] = scGrid[ i ][ 0 ];
      vol[ nCube + i ] = scGrid[ i ][ 1 ]; 
     } 
   
   int l = strlen( fileName );
   
   char fileNameSCRe2[ l + 15 ], fileNameSCIm2[ l + 15 ], fileNameElecRe2[ l + 15 ];
   
   strcpy( fileNameSCRe2, fileName );   
   strcpy( fileNameSCIm2, fileName );   
   strcpy( fileNameElecRe2, fileName );   
   
   if ( ( l > 3 ) && ( fileName[ l - 4 ] == '.' ) )
     fileNameSCRe2[ l - 4 ] = fileNameSCIm2[ l - 4 ] = fileNameElecRe2[ l - 4 ] = 0;
   
   strcat( fileNameSCRe2, "-SC-Re.rawiv" );
   strcat( fileNameSCIm2, "-SC-Im.rawiv" );   
   strcat( fileNameElecRe2, "-Elec-Re.rawiv" );

   if ( fileNameSCRe == NULL ) fileNameSCRe = fileNameSCRe2;
   if ( fileNameSCIm == NULL ) fileNameSCIm = fileNameSCIm2;
   if ( fileNameElecRe == NULL ) fileNameElecRe = fileNameElecRe2;   
   
   int retVal = 1;
   
   if ( !writeRAWIV( vol, n, xCenter, yCenter, zCenter, scale, fileNameSCRe ) 
     || !writeRAWIV( vol + nCube, n, xCenter, yCenter, zCenter, scale, fileNameSCIm ) 
     || !writeRAWIV( elecGrid, n, xCenter, yCenter, zCenter, scale, fileNameElecRe ) ) 
       retVal = 0;
   
   free( vol );
   
   return retVal;
}


int readRAWIVHeader( int *xDim, int *yDim, int *zDim, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileName )
{
   FILE *fp;
   
   if ( ( fp = fopen( fileName, "rb" ) ) == NULL )
     {
      fprintf( stderr, "\nError: Failed to open file %s!\n\n", fileName );
      return 0;
     }
     
   RAWIVHeader header;

   if ( fread( &header, sizeof( header ), 1, fp ) != 1 )  
     {
      fprintf( stderr, "\nError: Failed to read header from file %s!\n\n", fileName );
      fclose( fp );      
      return 0;      
     }

   if ( !bigEndian( ) )
     {
      for ( int i = 0; i < 3; i++ )
        {
         SWAP_32( &( header.minExt[ i ] ) );
         SWAP_32( &( header.maxExt[ i ] ) );         
         SWAP_32( &( header.dim[ i ] ) );         
         SWAP_32( &( header.origin[ i ] ) );                           
         SWAP_32( &( header.span[ i ] ) );         
        }
        
      SWAP_32( &( header.numVertices ) );  
      SWAP_32( &( header.numCells ) );        
     }
   
   *xCenter = ( header.minExt[ 0 ] + header.maxExt[ 0 ] ) / 2.0;
   *yCenter = ( header.minExt[ 1 ] + header.maxExt[ 1 ] ) / 2.0;
   *zCenter = ( header.minExt[ 2 ] + header.maxExt[ 2 ] ) / 2.0;      
      
   float gridLength1D = header.maxExt[ 0 ] - header.minExt[ 0 ];   
   *scale = 1.0 / gridLength1D;
   
   *xDim = header.dim[ 0 ];
   *yDim = header.dim[ 1 ];
   *zDim = header.dim[ 2 ];      
   
   fclose( fp );
   
   return 1;
}


int readRAWIV( FFTW_DATA_TYPE **vol, int *xDim, int *yDim, int *zDim, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileName )
{
   FILE *fp;
   
   if ( ( fp = fopen( fileName, "rb" ) ) == NULL )
     {
      fprintf( stderr, "\nError: Failed to open file %s!\n\n", fileName );
      return 0;
     }
     
   fseek( fp, 0, SEEK_END );  
   long fileSize = ftell( fp );
   
   fseek( fp, 0, SEEK_SET );     
   
   RAWIVHeader header;

   if ( fread( &header, sizeof( header ), 1, fp ) != 1 )  
     {
      fprintf( stderr, "\nError: Failed to read header from file %s!\n\n", fileName );
      fclose( fp );      
      return 0;      
     }

   if ( !bigEndian( ) )
     {
      for ( int i = 0; i < 3; i++ )
        {
         SWAP_32( &( header.minExt[ i ] ) );
         SWAP_32( &( header.maxExt[ i ] ) );         
         SWAP_32( &( header.dim[ i ] ) );         
         SWAP_32( &( header.origin[ i ] ) );                           
         SWAP_32( &( header.span[ i ] ) );         
        }
        
      SWAP_32( &( header.numVertices ) );  
      SWAP_32( &( header.numCells ) );        
     }
   
   *xCenter = ( header.minExt[ 0 ] + header.maxExt[ 0 ] ) / 2.0;
   *yCenter = ( header.minExt[ 1 ] + header.maxExt[ 1 ] ) / 2.0;
   *zCenter = ( header.minExt[ 2 ] + header.maxExt[ 2 ] ) / 2.0;      
      
   float gridLength1D = header.maxExt[ 0 ] - header.minExt[ 0 ];   
   *scale = 1.0 / gridLength1D;
   
   *xDim = header.dim[ 0 ];
   *yDim = header.dim[ 1 ];
   *zDim = header.dim[ 2 ];      
   
   int dataTypeSize = ( fileSize - sizeof( header ) ) / ( ( *xDim ) * ( *yDim ) * ( *zDim ) );
   
   int rowLen;
   
   if ( dataTypeSize == 4 ) rowLen = ( *xDim ) * sizeof( float );
   else rowLen = ( *xDim ) * sizeof( double );
   
   unsigned char *buf = ( unsigned char * ) calloc( rowLen, 1 );
   
   if ( buf == NULL )
     {
      fprintf( stderr, "\nError: Failed to allocate temporary buffer space!\n\n" );
      fclose( fp );      
      return 0;      
     }
     
   *vol = ( FFTW_DATA_TYPE * ) malloc( ( *xDim ) * ( *yDim ) * ( *zDim ) * sizeof( FFTW_DATA_TYPE ) );  
   
   if ( *vol == NULL )
     {
      fprintf( stderr, "\nError: Failed to allocate buffer space for reading file *s!\n\n", fileName );
      free( buf );
      fclose( fp );      
      return 0;      
     }   

   FFTW_DATA_TYPE *volPtr = ( FFTW_DATA_TYPE * ) *vol;
   
   for ( int k = 0; k < *zDim; k++ )
     for ( int j = 0; j < *yDim; j++ )
       {
        if ( fread( buf, rowLen, 1, fp ) != 1 )  
          {
           fprintf( stderr, "\nError: Failed to read data from file %s!\n\n", fileName );
           free( buf );           
           free( *vol );
           fclose( fp );      
           return 0;      
          }

        if ( dataTypeSize == 4 )
          {
           if ( !bigEndian( ) )
             {
              for ( int i = 0; i < *xDim; i++ )
                SWAP_32( buf + i * sizeof( float ) );
             }
          
           float *bufPtr = ( float * ) buf;

           for ( int i = 0; i < *xDim; i++ )
              volPtr[ i ] = ( FFTW_DATA_TYPE ) bufPtr[ i ];
          }   
        else  
          {
           if ( !bigEndian( ) )
             {
              for ( int i = 0; i < *xDim; i++ )
                SWAP_64( buf + i * sizeof( double ) );
             }
          
           double *bufPtr = ( double * ) buf;

           for ( int i = 0; i < *xDim; i++ )
              volPtr[ i ] = ( FFTW_DATA_TYPE ) bufPtr[ i ];
          }   
           
        volPtr += ( *xDim );
       }     
     
   free( buf );
   fclose( fp );
 
   printf( "\n\nREAD %s ( xDim = %d, yDim = %d, zDim = %d, xCenter = %lf, yCenter = %lf, zCenter = %lf, scale = %lf, dataTypeSize = %d )!\n\n", 
           fileName, *xDim, *yDim, *zDim, *xCenter, *yCenter, *zCenter, *scale, dataTypeSize );
   
   return 1;
}


int readShapeCompGrid( FFTW_complex **scGrid, int *n, double *xCenter, double *yCenter, double *zCenter, double *scale, 
                       char *fileNameRe, char *fileNameIm )
{
   if ( ( fileNameRe == NULL ) || ( fileNameIm == NULL ) ) return 0;

   FFTW_DATA_TYPE *vol;
   int xDim, yDim, zDim;

   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter, scale, fileNameRe ) ) return 0;
   
   if ( ( xDim != yDim ) || ( yDim != zDim ) || ( zDim != xDim ) ) 
     {
      free( vol );
      return 0;
     } 
   
   *n = xDim;
   int n3 = ( *n ) * ( *n ) * ( *n );
   *scGrid = ( FFTW_complex * ) FFTW_malloc( n3 * sizeof( FFTW_complex ) );
   
   if ( *scGrid == NULL )
     {
      free( vol );
      return 0;
     }
   
   for ( int c = 0; c < n3; c++ )
     ( *scGrid )[ c ][ 0 ] = vol[ c ];
   
   free( vol );

   double xc, yc, zc, sc;
   
   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, fileNameIm ) ) return 0;
   
   if ( ( xDim != ( *n ) ) || ( yDim != ( *n ) ) || ( zDim != ( *n ) ) 
     || ( xc != ( *xCenter ) ) || ( yc != ( *yCenter ) ) || ( zc != ( *zCenter ) ) || ( sc != ( *scale ) ) ) 
     {
      printf( "\nError: Parameter mismatch between files %s and %s!\n", fileNameRe, fileNameIm );
      FFTW_free( scGrid );
      free( vol );
      return 0;
     } 
      
   for ( int c = 0; c < n3; c++ )
     ( *scGrid )[ c ][ 1 ] = vol[ c ];

   free( vol );
      
   return 1;
}


int readShapeCompGrid( FFTW_complex *scGrid, double *xCenter, double *yCenter, double *zCenter, 
                       char *fileNameRe, char *fileNameIm )
{
   if ( ( scGrid == NULL ) || ( fileNameRe == NULL ) || ( fileNameIm == NULL ) ) return 0;

   FFTW_DATA_TYPE *vol;
   int xDim, yDim, zDim;
   double scale;

   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter, &scale, fileNameRe ) ) return 0;
   
   if ( ( xDim != yDim ) || ( yDim != zDim ) || ( zDim != xDim ) ) 
     {
      free( vol );
      return 0;
     } 
   
   int n = xDim;
   int n3 = n * n * n;
   
   for ( int c = 0; c < n3; c++ )
     scGrid[ c ][ 0 ] = vol[ c ];
   
   free( vol );

   double xc, yc, zc, sc;
   
   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, fileNameIm ) ) return 0;
   
   if ( ( xDim != n ) || ( yDim != n ) || ( zDim != n ) 
     || ( xc != ( *xCenter ) ) || ( yc != ( *yCenter ) ) || ( zc != ( *zCenter ) ) || ( sc != scale ) ) 
     {
      printf( "\nError: Parameter mismatch between files %s and %s!\n", fileNameRe, fileNameIm );
      free( vol );
      return 0;
     } 
      
   for ( int c = 0; c < n3; c++ )
     scGrid[ c ][ 1 ] = vol[ c ];

   free( vol );
      
   return 1;
}




int readElecGrid( FFTW_DATA_TYPE **elecGrid, int *n, double *xCenter, double *yCenter, double *zCenter, double *scale, char *fileNameRe )
{
   if ( fileNameRe == NULL ) return 0;

   FFTW_DATA_TYPE *vol;
   int xDim, yDim, zDim;

   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter, scale, fileNameRe ) ) return 0;
   
   if ( ( xDim != yDim ) || ( yDim != zDim ) || ( zDim != xDim ) ) 
     {
      free( vol );
      return 0;
     } 
   
   *n = xDim;
   int n3 = ( *n ) * ( *n ) * ( *n );
   *elecGrid = ( FFTW_DATA_TYPE * ) FFTW_malloc( n3 * sizeof( FFTW_complex ) );
   
   if ( *elecGrid == NULL )
     {
      free( vol );
      return 0;
     }
   
   for ( int c = 0; c < n3; c++ )
     ( *elecGrid )[ c ] = vol[ c ];
   
   free( vol );
      
   return 1;
}


int readElecGrid( FFTW_DATA_TYPE *elecGrid, double *xCenter, double *yCenter, double *zCenter, char *fileNameRe )
{
   if ( ( elecGrid == NULL ) || ( fileNameRe == NULL ) ) return 0;

   FFTW_DATA_TYPE *vol;
   int xDim, yDim, zDim;
   double scale;

   if ( !readRAWIV( &vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter, &scale, fileNameRe ) ) return 0;
   
   if ( ( xDim != yDim ) || ( yDim != zDim ) || ( zDim != xDim ) ) 
     {
      free( vol );
      return 0;
     } 
   
   int n = xDim;
   int n3 = n * n * n;
   
   for ( int c = 0; c < n3; c++ )
     elecGrid[ c ] = vol[ c ];
   
   free( vol );
      
   return 1;
}

