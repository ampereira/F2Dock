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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#endif

#include "fastfft.h"

#ifndef M_PI
#define 	M_PI   3.14159265358979323846
#endif

#ifndef min
#define min( a, b ) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#endif


/* 
   Embed a sum of Gaussians into a grid of size n^3.

   The following function performs STEP 1 (gridding) of algorithm 1.2 in section 1.2 of the following paper:

   Potts, D., Steidl G., and Tasche M.
   Fast Fourier transforms for nonequispaced data: A tutorial.
   in: Modern Sampling Theory: Mathematics and Applications, J.J. Benedetto and P. Ferreira (Eds.), 
   Chapter 12, pages 249-274, 1998.

   We are given M atoms, with the center of atom s \in [0, M - 1] in < x[ s ], y[ s ], z[ s ] >, radius
   in r[ s ], and the weight at the atom center in f[ s ]. Each center < x[ s ], y[ s ], z[ s ] >
   \in ( < -1/2 + m / n, -1/2 + m / n, -1/2 + m / n >, < 1/2 - m / n, 1/2 - m / n, 1/2 - m / n > ),
   where m is large enough so that m > n * max_s{ r[ s ] }, and n^3 is the size of the 3D output grid gHat.
   The smoothingFunction object gives access to precomputed Gaussian bell function
   Phi (which is actually truncated to [-m, m], and is referred to as Psi in the paper).

   We compute gHat( l ) = \sum_{ s \in [0, M): l - m <= n * w_s <= l + m }{ f_s * Phi( w_s - l / n ) },
   where l \in [ - n/2, n/2 )^3, w_s = < x[ s ], y[ s ], z[ s ] >, f_s = exp( - blobbiness * r[ s ] * r[ s ] ) * f[ s ], 
   and Phi( ) is the truncated [-m, m] Gaussian bell function used for smoothing.
   Actually, the gHat array computed below by the gridding function is shifted by ( n/2, n/2, n/2 ),
   that is, the origin is at ( n/2, n/2, n/2 ), and not at ( 0, 0, 0 ).
*/

void gridding( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f, double blobbiness, int n, int m, bool smoothSkin, SmoothingFunction* smoothingFunction, FFTW_complex* gHat, bool spreadSkin )
{
    double ofs = 0.5, m2 = m * m;
    
    for ( int c = 0; c < n * n * n; c++ )
        gHat[ c ][ 0 ] = gHat[ c ][ 1 ] = 0;

    //printf("gridding %d atoms: smoothSkin %d, spreadSkin %d numfreq %d m %d\n",M,  smoothSkin, spreadSkin, n, m);

    for ( int s = 0; s < M; s++ )
       {		
         // map from ( < -1/2 + m / n, -1/2 + m / n, -1/2 + m / n >, < 1/2 - m / n, 1/2 - m / n, 1/2 - m / n > )
         //       to ( < -n/2 + m, -n/2 + m, -n/2 + m >, < n/2 - m, n/2 - m, n/2 - m > )
         double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
         double rc = r[ s ] * ( n - 1 );
         double rc2 = rc * rc;

//         double xc = x[ s ] * n, yc = y[ s ] * n, zc = z[ s ] * n;
//         double rc = r[ s ];

         // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
         // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
  	 int xL, xU, yL, yU, zL, zU;

         xL = floor( xc ) - m, xU = ceil( xc ) + m;
	 yL = floor( yc ) - m, yU = ceil( yc ) + m;
	 zL = floor( zc ) - m, zU = ceil( zc ) + m;
         
         // apply atom s's influence
	 for ( int zt = zL; zt <= zU; zt++ )
	    for ( int yt = yL; yt <= yU; yt++ )
 	       for ( int xt = xL; xt <= xU; xt++ )
		  {
   	            double d2 = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc )
	  	              + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc )
	  	              + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );                                   
	  	              		  
                    // map from [ -n/2, n/2 ) to [ 0, n )
      		    int index = ( zt + n / 2 ) * n * n + ( yt + n / 2 ) * n + ( xt + n / 2 );
	 				
	            if ( smoothSkin )				
	              {
   	  	        if ( d2 > m2 ) continue;         
	              
                        // value of the smoothing function at [ xt + n / 2, yt + n / 2, zt + n / 2 ]                                        
              	        double bellVal = smoothingFunction->getPhi( ( ( zt + ofs ) - zc ) /*/ rc*/ ) 
              	                       * smoothingFunction->getPhi( ( ( yt + ofs ) - yc ) /*/ rc*/ ) 
              	                       * smoothingFunction->getPhi( ( ( xt + ofs ) - xc ) /*/ rc*/ );
	  	                   
                        // add atom s's contribution to gHat[ xt + n / 2, yt + n / 2, zt + n / 2 ]
                        gHat[ index ][ 0 ] += f[ s ][ 0 ] * bellVal;// * exp( -blobbiness * rc * rc );
                        gHat[ index ][ 1 ] += f[ s ][ 1 ] * bellVal;// * exp( -blobbiness * rc * rc ); 		        
 		      }  
 		    else
 		      {  
   	  	        if ( d2 > rc2 ) continue;         
 		      
                        double v = exp( blobbiness * ( d2 / ( rc * rc ) - 1 ) );

                        gHat[ index ][ 0 ] += v * f[ s ][ 0 ];
                        gHat[ index ][ 1 ] += v * f[ s ][ 1 ];
                      }
		  }
       }
       
      
    if ( spreadSkin )
      {
       printf( "building spread receptor skin" );
       fflush( stdout );
       int MCore = 0;
       
       for ( ; MCore < M; MCore++ )
          if ( type[ MCore ]  == 'E' ) break;
          
       double maxRe = gHat[ 0 ][ 0 ];    
             
       for ( int c = 1; c < n * n * n; c++ )
          if ( maxRe < gHat[ c ][ 0 ] ) maxRe = gHat[ c ][ 0 ];          
                                    
       double wss = 1, p = 2;
         
       for ( int s = 0; s < M; s++ )
         if ( type[ s ] == 'E' )
           {	         	
             // map from ( < -1/2 + m / n, -1/2 + m / n, -1/2 + m / n >, < 1/2 - m / n, 1/2 - m / n, 1/2 - m / n > )
             //       to ( < -n/2 + m, -n/2 + m, -n/2 + m >, < n/2 - m, n/2 - m, n/2 - m > )
             double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
             double rc = r[ s ] * ( n - 1 );
             double rc2 = rc * rc;
                          
             double md = 0.5 * m + rc; //ceil( pow( ( wss * f[ s ][ 0 ] ) / ( maxRe / ( ( M - MCore ) * 1000.0 ) ),  1 / p ) );
             double md2 = md * md;

             // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
             // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
  	     int xL = floor( xc ) - md, xU = ceil( xc ) + md;
	     int yL = floor( yc ) - md, yU = ceil( yc ) + md;
	     int zL = floor( zc ) - md, zU = ceil( zc ) + md;
             
             if ( xL < - n / 2 ) xL = - n / 2;
             if ( yL < - n / 2 ) yL = - n / 2;
             if ( zL < - n / 2 ) zL = - n / 2;
                                  
             if ( xU > n / 2 - 1 ) xU = n / 2 - 1;
             if ( yU > n / 2 - 1 ) yU = n / 2 - 1;
             if ( zU > n / 2 - 1 ) zU = n / 2 - 1;
                                       
             // apply atom s's influence
	     for ( int zt = zL; zt <= zU; zt++ )
	        for ( int yt = yL; yt <= yU; yt++ )
 	           for ( int xt = xL; xt <= xU; xt++ )
		      {
		        double d = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc ) 
		                 + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc ) 
		                 + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );
		        
		        if ( smoothSkin )
		          {
		            if ( d <= m2 ) continue;
		          }
		        else
		          {
		            if ( d <= rc2 ) continue;		          
		          }  
		                 
		        if ( d > md2 ) continue;        
		      
                        // map from [ -n/2, n/2 ) to [ 0, n - 1 )
      		        int index = ( zt + n / 2 ) * n * n + ( yt + n / 2 ) * n + ( xt + n / 2 );
      		        
      		        if ( gHat[ index ][ 1 ] > 0 ) continue;

      		        if ( gHat[ index ][ 0 ] == 0 ) continue;
				
//		        if ( p != 2 ) d = pow( d, p / 2 );
//		        
//		        if ( d < 1 ) d = 1.0;

//                        printf( "!" );     
//                        fflush( stdout ); 

                        double v = 1;//exp( blobbiness * ( d / ( rc * rc ) - 1 ) );
		      
		        // add atom s's contribution to gHat[ xt + n / 2, yt + n / 2, zt + n / 2 ]
		        gHat[ index ][ 0 ] += f[ s ][ 0 ] * v; //( wss * f[ s ][ 0 ] ) / d;
		      }
		      
             printf( "." );     
             fflush( stdout ); 
           }
           
       printf( "\n" );    
       fflush( stdout );                  
      }         
} 


void griddingElec( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f, double blobbiness, int n, double elecRadiusInGrids, FFTW_DATA_TYPE* gHat, bool forInPlaceFFT, bool movingMol )
{
    int nn;
    
    if ( forInPlaceFFT ) nn = 2 * ( n / 2 + 1 );
    else nn = n;
    
    for ( int c = 0; c < n * n * nn; c++ )
        gHat[ c ] = 0;

    for ( int s = 0; s < M; s++ )
       {		
         // map from ( < -1/2 + m / n, -1/2 + m / n, -1/2 + m / n >, < 1/2 - m / n, 1/2 - m / n, 1/2 - m / n > )
         //       to ( < -n/2 + m, -n/2 + m, -n/2 + m >, < n/2 - m, n/2 - m, n/2 - m > )
         double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
         double rc = elecRadiusInGrids;
         
         // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
         // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
         int xL = floor( xc ) - rc, xU = ceil( xc ) + rc;
         int yL = floor( yc ) - rc, yU = ceil( yc ) + rc;
         int zL = floor( zc ) - rc, zU = ceil( zc ) + rc;
     	 
         if ( xL < - ( n >> 1 ) ) xL = -( n >> 1 );
         if ( yL < - ( n >> 1 ) ) yL = -( n >> 1 );
         if ( zL < - ( n >> 1 ) ) zL = -( n >> 1 );
         	     	 
         if ( xU > ( n >> 1 ) - 1 ) xU = ( n >> 1 ) - 1;	     	 
         if ( yU > ( n >> 1 ) - 1 ) yU = ( n >> 1 ) - 1;
         if ( zU > ( n >> 1 ) - 1 ) zU = ( n >> 1 ) - 1;
         
         double ofs = 0.5;	     
         	     	     
         // apply atom s's influence
         for ( int zt = zL; zt <= zU; zt++ )
            for ( int yt = yL; yt <= yU; yt++ )
               for ( int xt = xL; xt <= xU; xt++ )
        	  { 
   	            double d2 = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc )
          	              + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc )
          	              + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );                                   
          	              
          	    if ( d2 > rc * rc ) continue;         
        	  
                    // map from [ -n/2, n/2 ) to [ 0, n )
      		    int index = ( zt + n / 2 ) * n * nn + ( yt + n / 2 ) * nn + ( xt + n / 2 );
        		        		
                    double v = exp( blobbiness * ( d2 / ( rc * rc ) - 1 ) );
  
                    gHat[ index ] += v * f[ s ][ 0 ];                      
                  }                                              
       }       
} 



void griddingHbond( int M, double* x, double* y, double* z, float *r, double rExp, FFTW_complex* f, 
                    double blobbiness, int n, FFTW_complex* gHat, bool movingMol )
{
    for ( int c = 0; c < n * n * n; c++ )
        gHat[ c ][ 0 ] = gHat[ c ][ 1 ] = 0;

    for ( int s = 0; s < M; s++ )
       {		
         if ( ( f[ s ][ 0 ] == 0 ) && ( f[ s ][ 1 ] == 0 ) ) continue;
//         if ( ( hbondType[ s ] != 'A' ) && ( hbondType[ s ] != 'D' ) ) continue;
       
         double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
         double rc = r[ s ] * ( n - 1 );
         
//         if ( hbondType[ s ] == 'A' ) rc += rExp;
         if ( ( movingMol && ( f[ s ][ 1 ] != 0 ) ) || ( !movingMol && ( f[ s ][ 0 ] != 0 ) ) ) rc += rExp;
         
         // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
         // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
         int xL = floor( xc ) - rc, xU = ceil( xc ) + rc;
         int yL = floor( yc ) - rc, yU = ceil( yc ) + rc;
         int zL = floor( zc ) - rc, zU = ceil( zc ) + rc;
     	 
         if ( xL < - ( n >> 1 ) ) xL = -( n >> 1 );
         if ( yL < - ( n >> 1 ) ) yL = -( n >> 1 );
         if ( zL < - ( n >> 1 ) ) zL = -( n >> 1 );
         	     	 
         if ( xU > ( n >> 1 ) - 1 ) xU = ( n >> 1 ) - 1;	     	 
         if ( yU > ( n >> 1 ) - 1 ) yU = ( n >> 1 ) - 1;
         if ( zU > ( n >> 1 ) - 1 ) zU = ( n >> 1 ) - 1;
         
         double ofs = 0.5;	     
         	     	     
         // apply atom s's influence
         for ( int zt = zL; zt <= zU; zt++ )
            for ( int yt = yL; yt <= yU; yt++ )
               for ( int xt = xL; xt <= xU; xt++ )
        	  { 
   	            double d2 = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc )
          	              + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc )
          	              + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );                                   
          	              
          	    if ( d2 > rc * rc ) continue; 
        	  
                    // map from [ -n/2, n/2 ) to [ 0, n )
      		    int index = ( zt + n / 2 ) * n * n + ( yt + n / 2 ) * n + ( xt + n / 2 );
        		        		
                    double v = exp( blobbiness * ( d2 / ( rc * rc ) - 1 ) );
  
                    gHat[ index ][ 0 ] += v * f[ s ][ 0 ];                      
                    gHat[ index ][ 1 ] += v * f[ s ][ 1 ];                                          
                  }                                              
       }       
} 



void griddingHydrophobicity( int M, double* x, double* y, double* z, float *r, FFTW_complex* f, 
                             double blobbiness, int n, FFTW_complex* gHat, double hydroRadExt, bool pairWise )
{
    for ( int c = 0; c < n * n * n; c++ )
        gHat[ c ][ 0 ] = gHat[ c ][ 1 ] = 0;

    for ( int s = 0; s < M; s++ )
       {		
         if ( ( f[ s ][ 0 ] == 0 ) && ( f[ s ][ 1 ] == 0 ) ) continue;
       
         double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
         double rc = r[ s ] * ( n - 1 ) + hydroRadExt;
         
         // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
         // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
         int xL = floor( xc ) - rc, xU = ceil( xc ) + rc;
         int yL = floor( yc ) - rc, yU = ceil( yc ) + rc;
         int zL = floor( zc ) - rc, zU = ceil( zc ) + rc;
     	 
         if ( xL < - ( n >> 1 ) ) xL = -( n >> 1 );
         if ( yL < - ( n >> 1 ) ) yL = -( n >> 1 );
         if ( zL < - ( n >> 1 ) ) zL = -( n >> 1 );
         	     	 
         if ( xU > ( n >> 1 ) - 1 ) xU = ( n >> 1 ) - 1;	     	 
         if ( yU > ( n >> 1 ) - 1 ) yU = ( n >> 1 ) - 1;
         if ( zU > ( n >> 1 ) - 1 ) zU = ( n >> 1 ) - 1;
         
         double ofs = 0.5;	     
         	     	     
         // apply atom s's influence
         for ( int zt = zL; zt <= zU; zt++ )
            for ( int yt = yL; yt <= yU; yt++ )
               for ( int xt = xL; xt <= xU; xt++ )
        	  { 
   	            double d2 = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc )
          	              + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc )
          	              + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );                                   
          	              
          	    if ( d2 > rc * rc ) continue; 
        	  
                    // map from [ -n/2, n/2 ) to [ 0, n )
      		    int index = ( zt + n / 2 ) * n * n + ( yt + n / 2 ) * n + ( xt + n / 2 );

                    if ( pairWise )        		        		
                      {
                        double v = exp( blobbiness * ( d2 / ( rc * rc ) - 1 ) );
  
                        gHat[ index ][ 0 ] += v * f[ s ][ 0 ];                      
                        gHat[ index ][ 1 ] += v * f[ s ][ 1 ];                                          
                      }
                    else
                      {    
                        if ( hydroRadExt > 0 ) 
                          {
                            gHat[ index ][ 0 ] = f[ s ][ 0 ];
                            gHat[ index ][ 1 ] = f[ s ][ 1 ];
                          }  
                        else
                          {
                            gHat[ index ][ 0 ] += f[ s ][ 0 ];
                            gHat[ index ][ 1 ] += f[ s ][ 1 ];                      
                          }
                      }  
                  }                                              
       }       
} 


void griddingSimpleComplementarity( int M, double* x, double* y, double* z, float *r, FFTW_complex* f, 
                                    double blobbiness, int n, FFTW_complex* gHat, double simpleRadExt )
{
    for ( int c = 0; c < n * n * n; c++ )
        gHat[ c ][ 0 ] = gHat[ c ][ 1 ] = 0;

    for ( int s = 0; s < M; s++ )
       {		
         if ( ( f[ s ][ 0 ] == 0 ) && ( f[ s ][ 1 ] == 0 ) ) continue;
       
         double xc = x[ s ] * ( n - 1 ), yc = y[ s ] * ( n - 1 ), zc = z[ s ] * ( n - 1 );
         double rc = r[ s ] * ( n - 1 ) + simpleRadExt;
         
         // compute the region [ xL, xU ] x [ yL, yU ] x [ zL, zU ] of atom s's influence;
         // observe that xL, xU, yL, yU, zL, zU \in [ -n/2, n/2 ) 
         int xL = floor( xc ) - rc, xU = ceil( xc ) + rc;
         int yL = floor( yc ) - rc, yU = ceil( yc ) + rc;
         int zL = floor( zc ) - rc, zU = ceil( zc ) + rc;
     	 
         if ( xL < - ( n >> 1 ) ) xL = -( n >> 1 );
         if ( yL < - ( n >> 1 ) ) yL = -( n >> 1 );
         if ( zL < - ( n >> 1 ) ) zL = -( n >> 1 );
         	     	 
         if ( xU > ( n >> 1 ) - 1 ) xU = ( n >> 1 ) - 1;	     	 
         if ( yU > ( n >> 1 ) - 1 ) yU = ( n >> 1 ) - 1;
         if ( zU > ( n >> 1 ) - 1 ) zU = ( n >> 1 ) - 1;
         
         double ofs = 0.5;	     
         	     	     
         // apply atom s's influence
         for ( int zt = zL; zt <= zU; zt++ )
            for ( int yt = yL; yt <= yU; yt++ )
               for ( int xt = xL; xt <= xU; xt++ )
        	  { 
   	            double d2 = ( ( xt + ofs ) - xc ) * ( ( xt + ofs ) - xc )
          	              + ( ( yt + ofs ) - yc ) * ( ( yt + ofs ) - yc )
          	              + ( ( zt + ofs ) - zc ) * ( ( zt + ofs ) - zc );                                   
          	              
          	    if ( d2 > rc * rc ) continue; 
        	  
                    // map from [ -n/2, n/2 ) to [ 0, n )
      		    int index = ( zt + n / 2 ) * n * n + ( yt + n / 2 ) * n + ( xt + n / 2 );
        		        		
                    double v = 1;//exp( blobbiness * ( d2 / ( rc * rc ) - 1 ) );
  
                    gHat[ index ][ 0 ] += v * f[ s ][ 0 ];                      
                    gHat[ index ][ 1 ] += v * f[ s ][ 1 ];                                          
                  }                                              
       }       
} 


/*
   NDFT algorithm is used to obtain the low frequencies from the input molecule.
   The frequencies are returned in the array ourFrequencies.
 
   It is an implementation of algorithm 1.2 in section 1.2 of the following paper:

   Potts, D., Steidl G., and Tasche M.
   Fast Fourier transforms for nonequispaced data: A tutorial.
   in: Modern Sampling Theory: Mathematics and Applications, J.J. Benedetto and P. Ferreira (Eds.), 
   Chapter 12, pages 249-274, 1998.

   We are given M atoms, with the center of atom s \in [0, M - 1] in < x[ s ], y[ s ], z[ s ] >, radius
   in r[ s ], and the weight at the atom center in f[ s ]. Each center < x[ s ], y[ s ], z[ s ] >
   \in ( < -1/2 + m / n, -1/2 + m / n, -1/2 + m / n >, < 1/2 - m / n, 1/2 - m / n, 1/2 - m / n > ),
   where m is large enough so that m > n * max_s{ r[ s ] }, n = alpha * N, alpha (>= 1) is the oversampling
   factor for NFFT, and N^3 is the number of frequencies we want to compute.
   The computed frequencies are returned in hHat, while gHat (of size n^3) and the 3D FFTW plan
   gHatPlan are used for intermediate computation (gridding and 3D FFT) by NFFT.
   The smoothingFunction object gives access to precomputed Gaussian bell function
   Phi (which is actually truncated to [-m, m], and is referred to as Psi in the paper),
   and its FFT PhiHat.
*/

bool getCenterFrequencies( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f,
                           double blobbiness, double alpha, int N, int m, 
			   FFTW_complex* hHat, bool smoothSkin, SmoothingFunction* smoothingFunction,
			   FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			   bool griddingDone, bool spreadSkin )
{
	if ( ( !griddingDone && ( !x || !y || !z || !f ) ) 
	  || ( smoothSkin && !smoothingFunction ) || !hHat || !gHat ) return false;

        int n = ( int ) ( alpha * N );
	n -= (n % 2);     // keep n even

        if ( !griddingDone )
          {
            // STEP 1: Compute gHat( l ) = \sum_{ s \in [0, M): l - m <= n * w_s <= l + m }{ f_s * Psi( w_s - l / n ) },
            //         where l \in [ - n/2, n/2 )^3, w_s = < x[ s ], y[ s ], z[ s ] >, f_s = exp( - blobbiness * r[ s ] * r[ s ] ) * f[ s ], 
            //         and Psi( ) is the truncated [-m, m] Gaussian bell function used for smoothing.
            //         Actually, the gHat array computed below by the gridding function is shifted by ( n/2, n/2, n/2 ),
            //         that is, the origin is at ( n/2, n/2, n/2 ), and not at ( 0, 0, 0 )
    
            gridding( M, x, y, z, r, type, f, blobbiness, n, m, smoothSkin, smoothingFunction, gHat, spreadSkin );       	    
          }  

        // STEP 2: Compute 3D FFT of gHat: for t \in [ - N/2, N/2 )^3,
        //               cHat( t ) = ( 1 / n^3 ) * \sum_{ l \in I_n }{ gHat( l ) * exp( - 2 * \Pi * i * t * l / n ) }.
        //         We actually compute cHat( t ) for t \in [ - n/2, n/2 )^3, and retain only the entries for 
        //         t \in [ - N/2, N/2 )^3 in step 3.
        //         The cHat array is computed inplace on gHat, but the origin of cHat is no longer at ( n/2, n/2, n/2 ).
        //         The function has period n in all 3 directions, and the origin is at 
        //           ( 0, 0, 0 ) = ( 0, 0, n - 1 ) = ( 0, n - 1, n - 1 ) = ( 0, n - 1, 0 ) 
        //         = ( n - 1, 0, 0 ) = ( n - 1, 0, n - 1 ) = ( n - 1, n - 1, n - 1 ) = ( n - 1, n - 1, 0 ).

       	if ( gHatSparsePlan != NULL ) sparse3DFFT( gHatSparsePlan, gHat, cHat ); 
       	else FFTW_execute( gHatPlan );	
        
        // STEP 3: Compute hHat( t ) = cHat( t ) / PhiHat( t ), for t \in [ - N/2, N/2 )^3.
        //         Observe that since the origin of cHat is at its eight corners, we only need to copy a cube of size (N/2)^3
        //         from each corner of gHat (= cHat) to the corresponding corner of hHat.

        if ( griddingDone )         
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
        		        if ( smoothSkin )
        		          {
        			    double ijkPsiBar = smoothingFunction->getPhiHat( min( i, n - i ) ) * smoothingFunction->getPhiHat( min( j, n - j ) ) * smoothingFunction->getPhiHat( min( k, n - k ) );
        
        			    hHat[  d  ][ 0 ] = cHat[ c ][ 0 ] / ijkPsiBar;
        			    hHat[ d++ ][ 1 ] = cHat[ c ][ 1 ] / ijkPsiBar;
                                  }
                                else
                                  {  			
        			    hHat[  d  ][ 0 ] = cHat[ c ][ 0 ];
        			    hHat[ d++ ][ 1 ] = cHat[ c ][ 1 ];			
        			  }
        		      }
          }		      
        else  
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
        		        if ( smoothSkin )
        		          {
        			    double ijkPsiBar = smoothingFunction->getPhiHat( min( i, n - i ) ) * smoothingFunction->getPhiHat( min( j, n - j ) ) * smoothingFunction->getPhiHat( min( k, n - k ) );
        
        			    hHat[  d  ][ 0 ] = gHat[ c ][ 0 ] / ijkPsiBar;
        			    hHat[ d++ ][ 1 ] = gHat[ c ][ 1 ] / ijkPsiBar;
                                  }
                                else
                                  {  			
        			    hHat[  d  ][ 0 ] = gHat[ c ][ 0 ];
        			    hHat[ d++ ][ 1 ] = gHat[ c ][ 1 ];			
        			  }
        		      }
          }		      
									
	return true;
}



bool getCenterElecFrequencies( int M, double* x, double* y, double* z, float *r, char *type, FFTW_complex* f,
                               double blobbiness, double alpha, int N, double elecRadiusInGrids, 
			       FFTW_complex* hHat, FFTW_DATA_TYPE* gHat, FFTW_plan gHatPlan, bool griddingDone, bool movingMol )
{
	if ( ( !griddingDone && ( !x || !y || !z || !f ) ) || !gHat ) return false;

        int n = ( int ) ( alpha * N );
	n -= (n % 2);     // keep n even

        if ( !griddingDone )
          {
            bool forInPlaceFFT;
            
            if ( gHat == ( FFTW_DATA_TYPE * ) hHat ) forInPlaceFFT = true;
            else forInPlaceFFT = false;
            
            griddingElec( M, x, y, z, r, type, f, blobbiness, n, elecRadiusInGrids, gHat, forInPlaceFFT, movingMol );       	    
          }

       	FFTW_execute( gHatPlan );	
        
	return true;
}



bool getCenterHbondFrequencies( int M, double* x, double* y, double* z, float *r, double rExp, FFTW_complex* f,
                                double blobbiness, double alpha, int N, 
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone, bool movingMol )
{
	if ( ( !griddingDone && ( !x || !y || !z || !f ) ) || !gHat ) return false;

        int n = ( int ) ( alpha * N );
	n -= (n % 2);     // keep n even

        if ( !griddingDone ) griddingHbond( M, x, y, z, r, rExp, f, blobbiness, n, gHat, movingMol );       	    

       	if ( gHatSparsePlan != NULL ) sparse3DFFT( gHatSparsePlan, gHat, cHat ); 
       	else FFTW_execute( gHatPlan );	

        if ( griddingDone )         
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
       			        hHat[  d  ][ 0 ] = cHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = cHat[ c ][ 1 ];			
        		      }
          }		      
        else  
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
        			hHat[  d  ][ 0 ] = gHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = gHat[ c ][ 1 ];			
        		      }
          }		      
									
	return true;
}



bool getCenterHydrophobicityFrequencies( int M, double* x, double* y, double* z, float *r, FFTW_complex* f,
                                double blobbiness, double alpha, int N, double hydroRadExt, bool pairWise,
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone )
{
	if ( ( !griddingDone && ( !x || !y || !z || !f ) ) || !gHat ) return false;

        int n = ( int ) ( alpha * N );
	n -= (n % 2);     // keep n even

        if ( !griddingDone ) griddingHydrophobicity( M, x, y, z, r, f, blobbiness, n, gHat, hydroRadExt, pairWise );       	    

       	if ( gHatSparsePlan != NULL ) sparse3DFFT( gHatSparsePlan, gHat, cHat ); 
       	else FFTW_execute( gHatPlan );	

        if ( griddingDone )         
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
       			        hHat[  d  ][ 0 ] = cHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = cHat[ c ][ 1 ];			
        		      }
          }		      
        else  
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
        			hHat[  d  ][ 0 ] = gHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = gHat[ c ][ 1 ];			
        		      }
          }		      
									
	return true;
}


bool getCenterSimpleComplementarityFrequencies( int M, double* x, double* y, double* z, float *r, FFTW_complex* f,
                                double blobbiness, double alpha, int N, double simpleRadExt,
			        FFTW_complex* hHat, FFTW_complex* gHat, FFTW_complex* cHat, FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, 
			        bool griddingDone )
{
	if ( ( !griddingDone && ( !x || !y || !z || !f ) ) || !gHat ) return false;

        int n = ( int ) ( alpha * N );
	n -= (n % 2);     // keep n even

        if ( !griddingDone ) griddingSimpleComplementarity( M, x, y, z, r, f, blobbiness, n, gHat, simpleRadExt );       	    

       	if ( gHatSparsePlan != NULL ) sparse3DFFT( gHatSparsePlan, gHat, cHat ); 
       	else FFTW_execute( gHatPlan );	

        if ( griddingDone )         
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
       			        hHat[  d  ][ 0 ] = cHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = cHat[ c ][ 1 ];			
        		      }
          }		      
        else  
          {
        	for ( int c = 0, d = 0, i = 0; i < n; i++ )
        	   for ( int j = 0; j < n; j++ )
        	      for ( int k = 0; k < n; k++, c++ )
                          if ( ( i < N / 2 || n - i <= N / 2 ) && ( j < N / 2 || n - j <= N / 2 ) && ( k < N / 2 || n - k <= N / 2 ) )
        		      {
        			hHat[  d  ][ 0 ] = gHat[ c ][ 0 ];
        			hHat[ d++ ][ 1 ] = gHat[ c ][ 1 ];			
        		      }
          }		      
									
	return true;
}
