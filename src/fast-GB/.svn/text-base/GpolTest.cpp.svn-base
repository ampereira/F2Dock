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
#include <string.h>

#include "fastBornRadius.h"
#include "fastGpol.h"
#include "../utils/utils.h"

using fastGB::fastBornRadius;
using fastGB::fastGpol;

typedef struct
  {
    char *pqrFile;    
    char *quadFile;

    int numThreadsBR, numThreadsGpol;
    double epsilonBR, epsilonGpol;
    bool useApproxMath;
    
    bool printStatus;

  } PARAMS_IN;


double getGpol( PARAMS_IN *params )
{
  fastBornRadius fastBR( params->quadFile, params->pqrFile, params->printStatus );
  
  fastBR.setNumThreads( params->numThreadsBR );
  fastBR.setEpsilon( params->epsilonBR );  

  fastBR.computeBornRadii( );
  
  fflush( stdout );
  
  int nAtoms;
  double *atomsPQRR = NULL;
  
  fastBR.getAtomsPQRR( &nAtoms, &atomsPQRR, false );
  
  fastGpol fastGp( nAtoms, atomsPQRR, params->printStatus );

  fastGp.setNumThreads( params->numThreadsGpol );
  fastGp.setEpsilon( params->epsilonGpol );
  fastGp.useApproxMathFunctions( params->useApproxMath );

  double Gpol;

  fastGp.computeFastGpol( &Gpol );
  
  fflush( stdout );  
  
  freeMem( atomsPQRR );  
  
  return Gpol;
}


double getNaiveGpol( PARAMS_IN *params )
{
  fastBornRadius fastBR( params->quadFile, params->pqrFile, params->printStatus );
  
  fastBR.computeBornRadiiNaively( );
  
  fflush( stdout );
  
  int nAtoms;
  double *atomsPQRR = NULL;
  
  fastBR.getAtomsPQRR( &nAtoms, &atomsPQRR, true );
  
  fastGpol fastGp( nAtoms, atomsPQRR, params->printStatus );

  double Gpol;

  fastGp.computeQuadGpol( &Gpol );

  fflush( stdout );
  
  freeMem( atomsPQRR );  
  
  return Gpol;
}


bool getParamsFromFile( PARAMS_IN *p, char *paramFile )
{
  char s[ 2000 ];
  char key[ 500 ], val[ 500 ];
  FILE *fp;

  fp = fopen( paramFile, "r" );

  if ( fp == NULL )
    {
      printError( (char *)"Failed to open parameter file %s!", paramFile );
      return false;
    }

  p->pqrFile = NULL;
  p->quadFile = NULL;
  p->numThreadsBR = p->numThreadsGpol = 4;
  p->epsilonBR = 0.5;
  p->epsilonGpol = 0.7;
  p->useApproxMath = false;
  p->printStatus = true;

  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, (char *)"%s %s", key, val ) != 2 ) continue;
    
      if ( !strcasecmp( key, (char *)"pqrFile" ) ) p->pqrFile = strdup( val );
      else if ( !strcasecmp( key, (char *)"quadFile" ) ) p->quadFile = strdup( val );
      else if ( !strcasecmp( key, (char *)"numThreadsBR" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numThreadsBR = v;  
             }
      else if ( !strcasecmp( key, (char *)"numThreadsGpol" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numThreadsGpol = v;                 
             }
      else if ( !strcasecmp( key, (char *)"epsilonBR" ) )
             {
               double v = atof( val );
               
               if ( v < 0 )
                 {
                   printError( (char *)"%s must be a non-negative float!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->epsilonBR = v;                 
             }
      else if ( !strcasecmp( key, (char *)"epsilonGpol" ) )
             {
               double v = atof( val );
               
               if ( v < 0 )
                 {
                   printError( (char *)"%s must be a non-negative float!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->epsilonGpol = v;                 
             }                    
      else if ( !strcasecmp( key, (char *)"useApproxMath" ) ) 
             {
	       if ( !strcasecmp( val, (char *)"true" ) ) p->useApproxMath = true;
	       else if ( !strcasecmp( val, (char *)"false" ) ) p->useApproxMath = false;
	            else  
	              {
		        printError( (char *)"%s must be a Boolean value!", key );
		        fclose( fp );
		        return false;
	              }	    
	     }
      else if ( !strcasecmp( key, (char *)"printStatus" ) ) 
             {
	       if ( !strcasecmp( val, (char *)"true" ) ) p->printStatus = true;
	       else if ( !strcasecmp( val, (char *)"false" ) ) p->printStatus = false;
	            else  
	              {
		        printError( (char *)"%s must be a Boolean value!", key );
		        fclose( fp );
		        return false;
	              }	    
	     }	     
    }

  fclose( fp );
    
  if ( p->pqrFile == NULL )  
    { 
      printError( (char *)"Missing PQR file!" );
      return false;
    }

  if ( p->quadFile == NULL )  
    { 
      printError( (char *)"Missing QUAD file!" );
      return false;
    }
        
  return true;
}



int main( int argc, char *argv[ ] )
{
  if ( argc < 2 )
    {
      printError( (char *)"Input text file not specified!" );
      return 1;
    }
     
  PARAMS_IN params;
   
  if ( !getParamsFromFile( &params, argv[ 1 ] ) ) return 1; 
 
  for ( int i = 9; i >= 1; i-- ) 
    for ( int j = 9; j >= 1; j-- )
      {
        params.epsilonBR = ( i * 1.0 ) / 10;
        params.epsilonGpol = ( j * 1.0 ) / 10;        
        printf( (char *)"G_pol = %lf kcal/mol\n\n", getGpol( &params ) );       
        fflush( stdout );       
      }  

  printf( (char *)"Naive G_pol = %lf kcal/mol\n\n", getNaiveGpol( &params ) );  
  fflush( stdout );
      
  return 0;
}
