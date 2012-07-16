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


#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "fastDispE.h"

using fastGB::fastDispE;


void fastDispE::freeMemory( void )
{
   freeMem( atoms );
   freeMem( qPoints );
   freeMem( staticMol.atoms );
   freeMem( staticMol.qPoints );   
   freeMem( movingMol.atoms );
   freeMem( movingMol.qPoints );   
}


fastDispE::fastDispE( char *paramFile, int numStaticAtoms, double *staticAtomsPQR, int numMovingAtoms, double *movingAtomsPQR, 
                      int numThreads, int numThreadsPerThread )
{
  if ( !getParamsFromFile( &params, paramFile ) ) exit( 1 );
  
  params.numThreads = numThreads;
  params.numThreadsPerThread = numThreadsPerThread;
  
  atoms = qPoints = NULL;  
  staticMol.atoms = staticMol.qPoints = NULL;
  movingMol.atoms = movingMol.qPoints = NULL;  
  
  getAtomsQuadsAndDispE( params.staticMoleculeQUAD, numStaticAtoms, staticAtomsPQR, &staticMol );

  getAtomsQuadsAndDispE( params.movingMoleculeQUAD, numMovingAtoms, movingAtomsPQR, &movingMol );  
  
  numAtoms = staticMol.numAtoms + movingMol.numAtoms;
  numQPoints = staticMol.numQPoints + movingMol.numQPoints;
  
  atoms = ( double * ) malloc( params.numThreads * ( staticMol.numAtoms + movingMol.numAtoms ) * 5 * sizeof( double ) );
  qPoints = ( double * ) malloc( params.numThreads * ( staticMol.numQPoints + movingMol.numQPoints ) * 7 * sizeof( double ) );  
  
  if ( ( atoms == NULL ) || ( qPoints == NULL ) )
    {
      printError( (char *)"Memory Allocation Failed!" );
      freeMemory( );
      exit( 1 );
    }
    
  copyStaticMoleculeToAllThreads( );
}


fastDispE::~fastDispE( )
{
   freeMemory( );
}



void fastDispE::getAtomsQuadsAndDispE( char *quadFile, int numAtoms, double *atomsPQR, MOLECULE_INFO *mol )
{
  fastGB::fastBornRadius fastBR( quadFile, numAtoms, atomsPQR, false );
  
  fastBR.getAtomsPQR( &( mol->numAtoms ), &( mol->atoms ) );
  fastBR.getQPoints( &( mol->numQPoints ), &( mol->qPoints ) );  
  
  fastBR.setNumThreads( params.numThreadsPerThread );
  fastBR.setEpsilon( params.epsilonBR );

  fastBR.computeBornRadii( );
  
  mol->dispE = fastBR.getDispersionEnergy( );    
}


bool fastDispE::getParamsFromFile( PARAMS_IN *p, char *paramFile )
{
  char s[ 2000 ];
  char key[ 500 ], val[ 500 ];
  FILE *fp;

  fp = fopen( paramFile, (char *)"r" );

  if ( fp == NULL )
    {
      printError( (char *)"Failed to open parameter file %s!", paramFile );
      return false;
    }

  p->staticMoleculeQUAD = p->staticMoleculeQUAD = NULL;
  p->numThreads = 4;
  p->epsilonBR = 0.3;
  
  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, (char *)"%s %s", key, val ) != 2 ) continue;
    
      if ( !strcasecmp( key, (char *)"staticMoleculeQUAD" ) ) p->staticMoleculeQUAD = strdup( val );
      else if ( !strcasecmp( key, (char *)"movingMoleculeQUAD" ) ) p->movingMoleculeQUAD = strdup( val );
//      else if ( !strcasecmp( key, (char *)"numThreads" ) )
//             {
//               int v = atoi( val );
//               
//               if ( v < 1 )
//                 {
//                   printError( (char *)"%s must be a positive integer!", key );
//                   fclose( fp );
//                   return false;
//                 }
//                 
//               p->numThreads = v;  
//             }
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
    }

  fclose( fp );
    
  if ( p->staticMoleculeQUAD == NULL )  
    { 
      printError( (char *)"Missing QUAD file name for the static molecule!" );
      return false;
    }

  if ( p->movingMoleculeQUAD == NULL )  
    { 
      printError( (char *)"Missing QUAD file name for the moving molecule!" );
      return false;
    }
        
  return true;
}


void fastDispE::copyStaticMoleculeToAllThreads( void )
{
   for ( int t = 0; t < params.numThreads; t++ )
     {   
       int j = t * numAtoms;
       
       for ( int i = 0; i < staticMol.numAtoms; i++, j++ )
         for ( int k = 0; k < 5; k++ )
           atoms[ 5 * j + k ] = staticMol.atoms[ 5 * i + k ];
           
       j = t * numQPoints;
       
       for ( int i = 0; i < staticMol.numQPoints; i++, j++ )
         for ( int k = 0; k < 7; k++ )
           qPoints[ 7 * j + k ] = staticMol.qPoints[ 7 * i + k ];           
     }          
}



void fastDispE::transformAndCopyMovingMolecule( int threadID, double *transMat )
{
   int j = threadID * numAtoms + staticMol.numAtoms;

   for ( int i = 0; i < movingMol.numAtoms; i++, j++ )
     {
       double x = movingMol.atoms[ 5 * i + 0 ],
              y = movingMol.atoms[ 5 * i + 1 ],
              z = movingMol.atoms[ 5 * i + 2 ],
              q = movingMol.atoms[ 5 * i + 3 ],
              r = movingMol.atoms[ 5 * i + 4 ];
              
       double xx = transMat[  0 ] * x + transMat[  1 ] * y + transMat[  2 ] * z + transMat[  3 ],
              yy = transMat[  4 ] * x + transMat[  5 ] * y + transMat[  6 ] * z + transMat[  7 ],
              zz = transMat[  8 ] * x + transMat[  9 ] * y + transMat[ 10 ] * z + transMat[ 11 ];       
              
       atoms[ 5 * j + 0 ] = xx;       
       atoms[ 5 * j + 1 ] = yy;
       atoms[ 5 * j + 2 ] = zz;
       atoms[ 5 * j + 3 ] = q;
       atoms[ 5 * j + 4 ] = r;                            
     }

   j = threadID * numQPoints + staticMol.numQPoints;

   for ( int i = 0; i < movingMol.numQPoints; i++, j++ )
     {
       double x  = movingMol.qPoints[ 7 * i + 0 ],
              y  = movingMol.qPoints[ 7 * i + 1 ],
              z  = movingMol.qPoints[ 7 * i + 2 ],
              nx = movingMol.qPoints[ 7 * i + 3 ],
              ny = movingMol.qPoints[ 7 * i + 4 ],
              nz = movingMol.qPoints[ 7 * i + 5 ],
              w  = movingMol.qPoints[ 7 * i + 6 ];
              
       double xx = transMat[  0 ] * x + transMat[  1 ] * y + transMat[  2 ] * z + transMat[  3 ],
              yy = transMat[  4 ] * x + transMat[  5 ] * y + transMat[  6 ] * z + transMat[  7 ],
              zz = transMat[  8 ] * x + transMat[  9 ] * y + transMat[ 10 ] * z + transMat[ 11 ];       
              
       double nxx = ( transMat[  0 ] * ( x + nx ) + transMat[  1 ] * ( y + ny ) + transMat[  2 ] * ( z + nz ) + transMat[  3 ] ) - xx,
              nyy = ( transMat[  4 ] * ( x + nx ) + transMat[  5 ] * ( y + ny ) + transMat[  6 ] * ( z + nz ) + transMat[  7 ] ) - yy,
              nzz = ( transMat[  8 ] * ( x + nx ) + transMat[  9 ] * ( y + ny ) + transMat[ 10 ] * ( z + nz ) + transMat[ 11 ] ) - zz;              
              
       qPoints[ 7 * j + 0 ] = xx;       
       qPoints[ 7 * j + 1 ] = yy;
       qPoints[ 7 * j + 2 ] = zz;
       qPoints[ 7 * j + 3 ] = nxx;
       qPoints[ 7 * j + 4 ] = nyy;
       qPoints[ 7 * j + 5 ] = nzz;
       qPoints[ 7 * j + 6 ] = w;                                   
     }
}


double fastDispE::getDispE( int threadID, double *transMat )
{
  transformAndCopyMovingMolecule( threadID, transMat );

  fastGB::fastBornRadius fastBR( numQPoints, qPoints + threadID * numQPoints * 7, numAtoms, atoms + threadID * numAtoms * 5, false );
  
  fastBR.setNumThreads( params.numThreadsPerThread );
  fastBR.setEpsilon( params.epsilonBR );

  fastBR.computeBornRadii( );
  
  double dispE = fastBR.getDispersionEnergy( );  
  
  return dispE;
}


double fastDispE::getDispE( int threadID )
{
  double trans[ ] = { 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0 };

  return getDispE( threadID, trans );
}


double fastDispE::getDelDispE( int threadID, double *transMat )
{
  transformAndCopyMovingMolecule( threadID, transMat );

  fastGB::fastBornRadius fastBR( numQPoints, qPoints + threadID * numQPoints * 7, numAtoms, atoms + threadID * numAtoms * 5, false );
  
  fastBR.setNumThreads( params.numThreadsPerThread );
  fastBR.setEpsilon( params.epsilonBR );

  fastBR.computeBornRadii( );
  
  double delDispE = fastBR.getDispersionEnergy( ) - ( staticMol.dispE + movingMol.dispE );  
  
  return delDispE;
}


double fastDispE::getDelDispE( int threadID )
{
  double trans[ ] = { 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0 };

  return getDelDispE( threadID, trans );
}
