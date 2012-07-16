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

#include "fastBornRadius.h"

using fastGB::fastBornRadius;


void fastBornRadius::freeMemory( void )
{
   freeMem( qPoints );
   freeMem( atoms );
   freeMem( qPointsOctree );
   freeMem( atomsOctree );  
   freeMem( atomsSubtreeRoots ); 
}


bool fastBornRadius::allocateAtomsSubtreeRootsArray( int nThreads )
{
   int arraySize = 1;
   
   while ( arraySize < 8 * nThreads )
      arraySize *= 8;
      
   atomsSubtreeRoots = ( int * ) malloc( arraySize * sizeof( int ) );      
   
   if ( atomsSubtreeRoots == NULL )
     {
      printError( (char *)"Failed to allocate memory!" );
      return false;
     }
     
   atomsSubtreeRootsSize = arraySize;  
        
   return true;  
}


void fastBornRadius::setDefaults( void )
{
   qPoints = NULL;
   atoms = NULL;
   qPointsOctree = NULL;
   atomsOctree = NULL;
   atomsSubtreeRoots = NULL;
   
   minRadius = 2.0;
   minRadiusUsed = -1;   
   
   maxLeafSize = 10;
   maxLeafSizeUsed = -1;   

   epsilon = 0.1;
   epsilonUsed = -1;
   
   maxBornRadius = 1000.0;
   maxBornRadiusInv = 1.0 / maxBornRadius;
   
   qPointsOctreeBuilt = false;
   atomsOctreeBuilt = false;
   
   priorComputationCleared = true;

   computedBornRadiiNaively = false;
   
   curNode = maxNode = 0;
   curSubtreeRoot = maxSubtreeRoot = 0;   
   
   qPointsOctreeRoot = -1;
   atomsOctreeRoot = -1;
   
   fastBRTime = naiveBRTime = 0;
   
   numThreads = 4;
   
   printStatus = true;

   vdwRad[ N ] = 1.55;
   vdwRad[ C ] = 1.70;
   vdwRad[ O ] = 1.40;
   vdwRad[ H ] = 1.20;
   vdwRad[ S ] = 1.85;
   vdwRad[ P ] = 1.90;   
}


void fastBornRadius::printCurrentSettings( void )
{
   printf( (char *)"\nCurrent Parameter Settings:\n" );
   printf( (char *)"\tminRadius = %lf, maxLeafSize = %d, epsilon = %lf, numThreads = %d\n\n", minRadius, maxLeafSize, epsilon, numThreads );
}


fastBornRadius::fastBornRadius( char *qPtsFile, char *atmsPQRFile, bool printStat )
{
   setDefaults( );
   
   printStatus = printStat;
      
   if ( !readQPoints( qPtsFile ) || !readAtoms( atmsPQRFile ) )
      {
       freeMemory( );
       exit( 1 );
      }
            
   if ( !allocateAtomsSubtreeRootsArray( numThreads ) )
     {
      freeMemory( );
      exit( 1 );
     }
   
   if ( printStatus ) printCurrentSettings( );
}


fastBornRadius::fastBornRadius( int numQPoints, double *qPts, int numAtoms, double *atms, bool printStat )
{
   setDefaults( );

   printStatus = printStat;
         
   if ( !copyQPointsFromArray( numQPoints, qPts ) || !copyAtomsFromArray( numAtoms, atms ) )
      {
       freeMemory( );
       exit( 1 );
      }
            
   if ( !allocateAtomsSubtreeRootsArray( numThreads ) )
     {
      freeMemory( );
      exit( 1 );
     }
   
   if ( printStatus ) printCurrentSettings( );
}



fastBornRadius::fastBornRadius( char *qPtsFile, int numAtoms, double *atms, bool printStat )
{
   setDefaults( );

   printStatus = printStat;
         
   if ( !readQPoints( qPtsFile ) || !copyAtomsFromArray( numAtoms, atms ) )
      {
       freeMemory( );
       exit( 1 );
      }
            
   if ( !allocateAtomsSubtreeRootsArray( numThreads ) )
     {
      freeMemory( );
      exit( 1 );
     }
   
   if ( printStatus ) printCurrentSettings( );
}



fastBornRadius::~fastBornRadius( )
{
   freeMemory( );
   if ( printStatus ) printf( (char *)"\n" );
}


bool fastBornRadius::setMinRadius( double minRad )
{
   if ( minRad < 0 )
     {
      printError( (char *)"minRadius must be a non-negative real number!" );
      return false;     
     }
     
   minRadius = minRad;

   if ( printStatus ) printf( (char *)"\nminRadius is set to %lf\n", minRad );
   
   return true;
}


bool fastBornRadius::setMaxLeafSize( int maxLfSize )
{
   if ( maxLfSize <= 0 )
     {
      printError( (char *)"maxLeafSize must be a positive integer!" );
      return false;     
     }
     
   maxLeafSize = maxLfSize;

   if ( printStatus ) printf( (char *)"\nmaxLeafSize is set to %d\n", maxLfSize );

   return true;
}


bool fastBornRadius::setEpsilon( double eps )
{
   if ( eps < 0 )
     {
      printError( (char *)"epsilon must be a non-negative real number!" );
      return false;     
     }
     
   epsilon = eps;

   if ( printStatus ) printf( (char *)"\nepsilon is set to %lf\n", eps );

   return true;
}


bool fastBornRadius::setMaxBornRadius( double maxBornRad )
{
   if ( maxBornRad <= 0 )
     {
      printError( (char *)"maxBornRadius must be a positive real number!" );
      return false;     
     }
     
   maxBornRadius = maxBornRad;
   maxBornRadiusInv = 1.0 / maxBornRad;   

   if ( printStatus ) printf( (char *)"\nmaxBornRadius is set to %lf\n", maxBornRad );
   
   return true;
}


bool fastBornRadius::setNumThreads( int nThreads )
{
   if ( nThreads < 1 )
     {
      printError( (char *)"numThreads must be a positive integer!" );
      return false;     
     }
     
   int *atomsSubtreeRootsT;  
   
   atomsSubtreeRootsT = atomsSubtreeRoots;
     
   if ( allocateAtomsSubtreeRootsArray( nThreads ) )
     {
      numThreads = nThreads;
      freeMem( atomsSubtreeRootsT );
      if ( atomsOctreeBuilt ) fillAtomsSubtreeRootsArray( atomsOctreeRoot, 1 );   
      if ( printStatus ) printf( (char *)"\nnumThreads is set to %d\n", numThreads );
      return true;
     }
   else
     {
      atomsSubtreeRoots = atomsSubtreeRootsT;
      if ( printStatus ) printf( (char *)"\nnumThreads remains unchanged ( %d )\n", numThreads );      
      return false;
     }            
}


void fastBornRadius::setPrintStatus( bool printStat )
{
   printStatus = printStat;
   
   if ( printStatus ) printf( (char *)"\nprintStatus is set to true\n" );
}


int fastBornRadius::nextFreeNode( void )
{
   int nextNode = -1;
   
   pthread_mutex_lock( &nodesLock );
   if ( curNode < maxNode ) nextNode = curNode++;
   pthread_mutex_unlock( &nodesLock );
   
   return nextNode;
}


void fastBornRadius::initFreeNodeServer( int numNodes )
{
   pthread_mutex_init( &nodesLock, NULL );
   curNode = 0;
   maxNode = numNodes;
}


int fastBornRadius::nextSubtreeRoot( void )
{
   int nextRoot = -1;
   
   pthread_mutex_lock( &subtreeRootsLock );
   if ( curSubtreeRoot < maxSubtreeRoot ) nextRoot = curSubtreeRoot++;
   pthread_mutex_unlock( &subtreeRootsLock );
   
   return ( nextRoot == -1 ) ? ( - 1 ) : atomsSubtreeRoots[ nextRoot ];
}


void fastBornRadius::initSubtreeRootServer( void )
{
   pthread_mutex_init( &subtreeRootsLock, NULL );
   curSubtreeRoot = 0;
}


bool fastBornRadius::readQPoints( char *qPtsFile )
{
   if ( printStatus ) printf( (char *)"\nreading quadrature points from %s... ", qPtsFile );

   double startT = getTime( );
   
   FILE *fp;
   
   fp = fopen( qPtsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open quadrature points file (%s)!", qPtsFile );
      return false;
     }
   
   nQPoints = 0;
   QPOINT qPt;  
          
   while ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf %lf %lf", &qPt.x, &qPt.y, &qPt.z, &qPt.nx, &qPt.ny, &qPt.nz, &qPt.w ) == 7 ) nQPoints++; 

   fclose( fp );  
   
   qPoints = ( QPOINT * ) malloc( nQPoints * sizeof( QPOINT ) );
   
   if ( qPoints == NULL )
     {
      printError( (char *)"Failed to allocate memory for quadrature points!" );
      return false;
     }
     
   fp = fopen( qPtsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open quadrature points file (%s)!", qPtsFile );
      return false;
     }

   for ( int i = 0; i < nQPoints; i++ )   
     {          
      if ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf %lf %lf", &qPt.x, &qPt.y, &qPt.z, &qPt.nx, &qPt.ny, &qPt.nz, &qPt.w ) != 7 )
        {
         printError( (char *)"Failed to read the quadrature points file (%s)!", qPtsFile );
         return false;
        }
         
      qPoints[ i ] = qPt;         
     }    
     
   fclose( fp );

   computedBornRadiiNaively = false;
   
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, read %d quadrature points )\n", endT - startT, nQPoints );
         
   return true;
}



bool fastBornRadius::copyQPointsFromArray( int numQPoints, double *qPts )
{
   if ( printStatus ) printf( (char *)"\ncopying quadrature points from input array... " );

   double startT = getTime( );
   
   if ( numQPoints <= 0 )
     {
      printError( (char *)"No quadrature points to copy!" );
      return false;
     }
   
   nQPoints = numQPoints;
   qPoints = ( QPOINT * ) malloc( nQPoints * sizeof( QPOINT ) );
   
   if ( qPoints == NULL )
     {
      printError( (char *)"Failed to allocate memory for quadrature points!" );
      return false;
     }

   for ( int i = 0; i < nQPoints; i++ )   
     {          
      qPoints[ i ].x = qPts[ 7 * i + 0 ];
      qPoints[ i ].y = qPts[ 7 * i + 1 ];
      qPoints[ i ].z = qPts[ 7 * i + 2 ];

      qPoints[ i ].nx = qPts[ 7 * i + 3 ];
      qPoints[ i ].ny = qPts[ 7 * i + 4 ];
      qPoints[ i ].nz = qPts[ 7 * i + 5 ];
      
      qPoints[ i ].w = qPts[ 7 * i + 6 ];
     }    
     
   computedBornRadiiNaively = false;
   
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, copied %d quadrature points )\n", endT - startT, nQPoints );
         
   return true;
}


bool fastBornRadius::readAtomsFromPQR( char *atomsFile )
{
   if ( printStatus ) printf( (char *)"\nreading atoms from %s... ", atomsFile );

   double startT = getTime( );

   FILE *fp;
   
   fp = fopen( atomsFile, (char *)"rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PQR file (%s)!", atomsFile );
      return false;
     }
   
   nAtoms = 0;
   
   char line[ 101 ];

   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0;
      while ( isspace( line[ i ] ) ) i++;
      if ( strstr( line + i, (char *)"ATOM" ) == line + i ) nAtoms++;
     }
   
   fclose( fp );  
   
   atoms = ( ATOM * ) malloc( nAtoms * sizeof( ATOM ) );
   
   if ( atoms == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }
     
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PQR file (%s)!", atomsFile );
      return false;
     }

   ATOM atm;  
   
   int k = 0;
   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j;
      double v;
      char tmp[ 100 ], anm[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      i = getInt( line, i, &j );            // get atom number, and ignore

      i = getString( line, i, anm );        // get atom name

      i = getAlphaString( line, i, tmp );   // get residue name, and ignore

      i = getInt( line, i, &j );            // get residue number, and ignore

      i = getDouble( line, i, &v );         // get X coordinate
      atoms[ k ].x = v;    
    
      i = getDouble( line, i, &v );         // get Y coordinate
      atoms[ k ].y = v;    

      i = getDouble( line, i, &v );         // get Z coordinate
      atoms[ k ].z = v;    

      i = getDouble( line, i, &v );         // get charge
      atoms[ k ].q = v;    

      i = getDouble( line, i, &v );         // get radius
//      atoms[ k ].r = v;  
      
      switch ( anm[ 0 ] )
        { 
          case 'N' : atoms[ k ].r = vdwRad[ N ];
                     break;

          case 'C' : atoms[ k ].r = vdwRad[ C ];
                     break;

          case 'O' : atoms[ k ].r = vdwRad[ O ];
                     break;

          case 'H' : atoms[ k ].r = vdwRad[ H ];
                     break;

          case 'S' : atoms[ k ].r = vdwRad[ S ];
                     break;

          case 'P' : atoms[ k ].r = vdwRad[ P ];
                     break;
                     
          default  : atoms[ k ].r = v;    
                     break;           
        }        
      
      atoms[ k ].id = k + 1;
      atoms[ k ].invR = atoms[ k ].invCR = 0;      
      atoms[ k ].R = atoms[ k ].naiveR = 0;

      k++;
     }
     
   fclose( fp );

   computedBornRadiiNaively = false;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, read %d atoms )\n", endT - startT, nAtoms );
         
   return true;
}


bool fastBornRadius::readAtoms( char *atomsFile )
{
   int l = strlen( atomsFile );
   
   if ( ( l > 3 ) && ( atomsFile[ l - 4 ] == '.' ) && ( toupper( atomsFile[ l - 3 ] ) == 'P' ) 
     && ( toupper( atomsFile[ l - 2 ] ) == 'Q' ) && ( toupper( atomsFile[ l - 1 ] ) == 'R' ) )
     return readAtomsFromPQR( atomsFile );
   else  
     {
      printError( (char *)"Unknown file type (%s)!", atomsFile );
      return false;
     }
}


bool fastBornRadius::copyAtomsFromArray( int numAtoms, double *atms )
{
   if ( printStatus ) printf( (char *)"\ncopying atoms from array... " );

   double startT = getTime( );
   
   if ( numAtoms <= 0 )
     {
      printError( (char *)"No atoms to copy!" );
      return false;
     }
   
   nAtoms = numAtoms;
   atoms = ( ATOM * ) malloc( nAtoms * sizeof( ATOM ) );
   
   if ( atoms == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }

   for ( int i = 0; i < nAtoms; i++ )   
     {          
      atoms[ i ].x = atms[ 5 * i + 0 ];
      atoms[ i ].y = atms[ 5 * i + 1 ];
      atoms[ i ].z = atms[ 5 * i + 2 ];

      atoms[ i ].q = atms[ 5 * i + 3 ];

      atoms[ i ].r = atms[ 5 * i + 4 ];

      atoms[ i ].id = i + 1;
      atoms[ i ].invR = atoms[ i ].invCR = 0;      
      atoms[ i ].R = atoms[ i ].naiveR = 0;
     }    
     
   computedBornRadiiNaively = false;
   
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, copied %d atoms )\n", endT - startT, nAtoms );
         
   return true;
}



void fastBornRadius::transformAtomsOctree( int nodeID, double *transMatrix )
{
   if ( nodeID < 0 ) return;
   
   double cx = atomsOctree[ nodeID ].cx;
   double cy = atomsOctree[ nodeID ].cy;
   double cz = atomsOctree[ nodeID ].cz;
   
   atomsOctree[ nodeID ].cx = transMatrix[  0 ] * cx + transMatrix[  1 ] * cy + transMatrix[  2 ] * cz + transMatrix[  3 ];
   atomsOctree[ nodeID ].cy = transMatrix[  4 ] * cx + transMatrix[  5 ] * cy + transMatrix[  6 ] * cz + transMatrix[  7 ];
   atomsOctree[ nodeID ].cz = transMatrix[  8 ] * cx + transMatrix[  9 ] * cy + transMatrix[ 10 ] * cz + transMatrix[ 11 ];            

   for ( int i = 0; i < 8; i++ ) 
     if ( atomsOctree[ nodeID ].cPtr[ i ] >= 0 )
        transformAtomsOctree( atomsOctree[ nodeID ].cPtr[ i ], transMatrix );
}


void fastBornRadius::transformAtoms( double *transMatrix )
{
   if ( printStatus ) 
     {
       printf( (char *)"\ntransforming atoms using the following matrix...\n" );

       for ( int i = 0; i < 16; i++ )
         {
          printf( (char *)"%20lf  ", transMatrix[ i ] );
          if ( i % 4 == 3 ) printf( (char *)"\n" );
         }
    
       printf( "\n" );
     }  

   double startT = getTime( );
   
   for ( int i = 0; i < nAtoms; i++ )
     {
      double x = atoms[ i ].x, y = atoms[ i ].y, z = atoms[ i ].z;
      
      atoms[ i ].x = transMatrix[  0 ] * x + transMatrix[  1 ] * y + transMatrix[  2 ] * z + transMatrix[  3 ];
      atoms[ i ].y = transMatrix[  4 ] * x + transMatrix[  5 ] * y + transMatrix[  6 ] * z + transMatrix[  7 ];
      atoms[ i ].z = transMatrix[  8 ] * x + transMatrix[  9 ] * y + transMatrix[ 10 ] * z + transMatrix[ 11 ];            
     }
     
   if ( atomsOctreeBuilt ) transformAtomsOctree( atomsOctreeRoot, transMatrix ); 

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
}


void fastBornRadius::countQPointsOctreeNodesAndSortQPoints( QPOINT *qPts, int qPtsStartID, int qPtsEndID, 
                                                            QPOINT *qPtsT, int *numNodes )
{
   double minX = qPts[ qPtsStartID ].x, minY = qPts[ qPtsStartID ].y, minZ = qPts[ qPtsStartID ].z;
   double maxX = qPts[ qPtsStartID ].x, maxY = qPts[ qPtsStartID ].y, maxZ = qPts[ qPtsStartID ].z;
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      if ( qPts[ i ].x < minX ) minX = qPts[ i ].x;      
      if ( qPts[ i ].x > maxX ) maxX = qPts[ i ].x;      
      
      if ( qPts[ i ].y < minY ) minY = qPts[ i ].y;      
      if ( qPts[ i ].y > maxY ) maxY = qPts[ i ].y;      

      if ( qPts[ i ].z < minZ ) minZ = qPts[ i ].z;      
      if ( qPts[ i ].z > maxZ ) maxZ = qPts[ i ].z;      
     } 
   
   double cx = ( minX + maxX ) / 2,
          cy = ( minY + maxY ) / 2,
          cz = ( minZ + maxZ ) / 2;

   double r2 = ( qPts[ qPtsStartID ].x - cx ) * ( qPts[ qPtsStartID ].x - cx )
             + ( qPts[ qPtsStartID ].y - cy ) * ( qPts[ qPtsStartID ].y - cy )
             + ( qPts[ qPtsStartID ].z - cz ) * ( qPts[ qPtsStartID ].z - cz );
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      double r2T = ( qPts[ i ].x - cx ) * ( qPts[ i ].x - cx )
                 + ( qPts[ i ].y - cy ) * ( qPts[ i ].y - cy )
                 + ( qPts[ i ].z - cz ) * ( qPts[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   *numNodes = 1;
      
   if ( ( qPtsEndID - qPtsStartID + 1 > maxLeafSize ) && ( r2 > minRadius * minRadius ) )
     {
      int qPtsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {
         qPtsT[ i ] = qPts[ i ];
         
         int j = ( zeroIfLess( qPts[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPts[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPts[ i ].x, cx ) );
         
         qPtsCount[ j ]++;
        }
      
      int qPtsStartIndex[ 8 ];
      int qPtsCurIndex[ 8 ];              
      
      qPtsCurIndex[ 0 ] = qPtsStartIndex[ 0 ] = qPtsStartID;
      for ( int i = 1; i < 8; i++ )
        qPtsCurIndex[ i ] = qPtsStartIndex[ i ] = qPtsStartIndex[ i - 1 ] + qPtsCount[ i - 1 ];

      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {        
         int j = ( zeroIfLess( qPtsT[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPtsT[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPtsT[ i ].x, cx ) );
           
         qPts[ qPtsCurIndex[ j ] ] = qPtsT[ i ];
         qPtsCurIndex[ j ]++;  
        }        
        
      for ( int i = 0; i < 8; i++ ) 
        if ( qPtsCount[ i ] > 0 )
          {
           int numNodesT = 0;
           
           countQPointsOctreeNodesAndSortQPoints( qPts, qPtsStartIndex[ i ], qPtsStartIndex[ i ] + qPtsCount[ i ] - 1, qPtsT, &numNodesT );
         
           *numNodes += numNodesT;
          }
     }  
}



int fastBornRadius::constructQPointsOctree( int qPtsStartID, int qPtsEndID, QPOINT *qPts )
{
   int nodeID = nextFreeNode( );
   
   qPointsOctree[ nodeID ].qPtsStartID = qPtsStartID;
   qPointsOctree[ nodeID ].qPtsEndID = qPtsEndID;

   qPointsOctree[ nodeID ].wnx = qPointsOctree[ nodeID ].wny = qPointsOctree[ nodeID ].wnz = 0;
   
   double minX = qPts[ qPtsStartID ].x, minY = qPts[ qPtsStartID ].y, minZ = qPts[ qPtsStartID ].z;
   double maxX = qPts[ qPtsStartID ].x, maxY = qPts[ qPtsStartID ].y, maxZ = qPts[ qPtsStartID ].z;
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      if ( qPts[ i ].x < minX ) minX = qPts[ i ].x;      
      if ( qPts[ i ].x > maxX ) maxX = qPts[ i ].x;      
      
      if ( qPts[ i ].y < minY ) minY = qPts[ i ].y;      
      if ( qPts[ i ].y > maxY ) maxY = qPts[ i ].y;      

      if ( qPts[ i ].z < minZ ) minZ = qPts[ i ].z;      
      if ( qPts[ i ].z > maxZ ) maxZ = qPts[ i ].z;      
     } 
   
   double cx = qPointsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = qPointsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = qPointsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   double r2 = ( qPts[ qPtsStartID ].x - cx ) * ( qPts[ qPtsStartID ].x - cx )
             + ( qPts[ qPtsStartID ].y - cy ) * ( qPts[ qPtsStartID ].y - cy )
             + ( qPts[ qPtsStartID ].z - cz ) * ( qPts[ qPtsStartID ].z - cz );
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      double r2T = ( qPts[ i ].x - cx ) * ( qPts[ i ].x - cx )
                 + ( qPts[ i ].y - cy ) * ( qPts[ i ].y - cy )
                 + ( qPts[ i ].z - cz ) * ( qPts[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   double cr = qPointsOctree[ nodeID ].cr = sqrt( r2 );
         
   if ( ( qPtsEndID - qPtsStartID +  1 > maxLeafSize ) && ( r2 > minRadius * minRadius ) /*( cr > minRadius )*/ )   
     {
      qPointsOctree[ nodeID ].leaf = false;     
      
      int qPtsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {
         int j = ( zeroIfLess( qPts[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPts[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPts[ i ].x, cx ) );
         
         qPtsCount[ j ]++;
        }

      int qPtsStartIndex[ 8 ];

      qPtsStartIndex[ 0 ] = qPtsStartID;
      for ( int i = 1; i < 8; i++ )
        qPtsStartIndex[ i ] = qPtsStartIndex[ i - 1 ] + qPtsCount[ i - 1 ];

      for ( int i = 0; i < 8; i++ ) 
        if ( qPtsCount[ i ] > 0 )
          {
           int j = constructQPointsOctree( qPtsStartIndex[ i ], qPtsStartIndex[ i ] + qPtsCount[ i ] - 1, qPts );

           qPointsOctree[ nodeID ].cPtr[ i ] = j; 
           
           qPointsOctree[ nodeID ].wnx += qPointsOctree[ j ].wnx;
           qPointsOctree[ nodeID ].wny += qPointsOctree[ j ].wny;           
           qPointsOctree[ nodeID ].wnz += qPointsOctree[ j ].wnz;           
          }
        else qPointsOctree[ nodeID ].cPtr[ i ] = -1;           
     }  
   else
     {
      qPointsOctree[ nodeID ].leaf = true;
              
      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
         {
          qPointsOctree[ nodeID ].wnx += qPts[ i ].w * qPts[ i ].nx;
          qPointsOctree[ nodeID ].wny += qPts[ i ].w * qPts[ i ].ny;
          qPointsOctree[ nodeID ].wnz += qPts[ i ].w * qPts[ i ].nz;                    
         }
     }
     
   return nodeID;  
}



bool fastBornRadius::buildQPointsOctree( void )
{  
   if ( printStatus ) printf( (char *)"\nbuilding quadrature points octree... " );

   double startT = getTime( );

   QPOINT *qPointsT;         
   qPointsT = ( QPOINT * ) malloc( nQPoints * sizeof( QPOINT ) );
   
   if ( qPointsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for quadrature points!" );
      if ( !qPointsOctreeBuilt ) exit( 1 );
      return false;
     }

   countQPointsOctreeNodesAndSortQPoints( qPoints, 0, nQPoints - 1, qPointsT, &numQPointsOctreeNodes );

   freeMem( qPointsT );
   
   QPOINTS_OCTREE_NODE *qPointsOctreeT;
   
   qPointsOctreeT = ( QPOINTS_OCTREE_NODE * ) malloc( numQPointsOctreeNodes * sizeof( QPOINTS_OCTREE_NODE ) );
   
   if ( qPointsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s quadrature points octree - memory allocation failed!", ( qPointsOctreeBuilt ) ? (char *)"rebuild" : (char *)"build" );
      if ( !qPointsOctreeBuilt ) exit( 1 );
      return false;
     }
 
   freeMem( qPointsOctree );

   qPointsOctree = qPointsOctreeT;  

   initFreeNodeServer( nQPoints );
   
   qPointsOctreeRoot = constructQPointsOctree( 0, nQPoints - 1, qPoints );
   
   qPointsOctreeBuilt = true;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}




void fastBornRadius::countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID, 
                                                        ATOM *atomsT, int *numNodes )
{
   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      
     } 
   
   double cx = ( minX + maxX ) / 2,
          cy = ( minY + maxY ) / 2,
          cz = ( minZ + maxZ ) / 2;

   double r2 = ( atoms[ atomsStartID ].x - cx ) * ( atoms[ atomsStartID ].x - cx )
             + ( atoms[ atomsStartID ].y - cy ) * ( atoms[ atomsStartID ].y - cy )
             + ( atoms[ atomsStartID ].z - cz ) * ( atoms[ atomsStartID ].z - cz );
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      double r2T = ( atoms[ i ].x - cx ) * ( atoms[ i ].x - cx )
                 + ( atoms[ i ].y - cy ) * ( atoms[ i ].y - cy )
                 + ( atoms[ i ].z - cz ) * ( atoms[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   *numNodes = 1;
      
   if ( ( atomsEndID - atomsStartID +  1 > maxLeafSize ) && ( r2 > minRadius * minRadius ) )
     {
      int atomsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {
         atomsT[ i ] = atoms[ i ];
         
         int j = ( zeroIfLess( atoms[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atoms[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atoms[ i ].x, cx ) );
         
         atomsCount[ j ]++;
        }
        
      int atomsStartIndex[ 8 ];
      int atomsCurIndex[ 8 ];              
      
      atomsCurIndex[ 0 ] = atomsStartIndex[ 0 ] = atomsStartID;
      for ( int i = 1; i < 8; i++ )
        atomsCurIndex[ i ] = atomsStartIndex[ i ] = atomsStartIndex[ i - 1 ] + atomsCount[ i - 1 ];

      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {        
         int j = ( zeroIfLess( atomsT[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atomsT[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atomsT[ i ].x, cx ) );
           
         atoms[ atomsCurIndex[ j ] ] = atomsT[ i ];
         atomsCurIndex[ j ]++;  
        }        
        
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCount[ i ] > 0 )
          {
           int numNodesT = 0;
           
           countAtomsOctreeNodesAndSortAtoms( atoms, atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atomsT, &numNodesT );
         
           *numNodes += numNodesT;
          }
     }  
}



int fastBornRadius::constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms )
{
   int nodeID = nextFreeNode( );
   
   atomsOctree[ nodeID ].atomsStartID = atomsStartID;
   atomsOctree[ nodeID ].atomsEndID = atomsEndID;
   
   atomsOctree[ nodeID ].invR = atomsOctree[ nodeID ].invCR = 0.0;   

   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      
     } 
   
   double cx = atomsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = atomsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = atomsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   double r2 = ( atoms[ atomsStartID ].x - cx ) * ( atoms[ atomsStartID ].x - cx )
             + ( atoms[ atomsStartID ].y - cy ) * ( atoms[ atomsStartID ].y - cy )
             + ( atoms[ atomsStartID ].z - cz ) * ( atoms[ atomsStartID ].z - cz );
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      double r2T = ( atoms[ i ].x - cx ) * ( atoms[ i ].x - cx )
                 + ( atoms[ i ].y - cy ) * ( atoms[ i ].y - cy )
                 + ( atoms[ i ].z - cz ) * ( atoms[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   double cr = atomsOctree[ nodeID ].cr = sqrt( r2 );
         
   if ( ( atomsEndID - atomsStartID +  1 <= maxLeafSize ) || ( cr <= minRadius ) )
      atomsOctree[ nodeID ].leaf = true;
   else
     {
      atomsOctree[ nodeID ].leaf = false;     
      
      int atomsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {
         int j = ( zeroIfLess( atoms[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atoms[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atoms[ i ].x, cx ) );
         
         atomsCount[ j ]++;
        }
        
      int atomsStartIndex[ 8 ];

      atomsStartIndex[ 0 ] = atomsStartID;
      for ( int i = 1; i < 8; i++ )
        atomsStartIndex[ i ] = atomsStartIndex[ i - 1 ] + atomsCount[ i - 1 ];

      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCount[ i ] > 0 )
          {
           int j = constructAtomsOctree( atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atoms );
           atomsOctree[ nodeID ].cPtr[ i ] = j; 
          }
        else atomsOctree[ nodeID ].cPtr[ i ] = -1;           
     }  
     
   return nodeID;  
}


void fastBornRadius::fillAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel )
{
   if ( nodeID < 0 ) return;
   
   if ( nodeID == atomsOctreeRoot ) maxSubtreeRoot = 0;
   
   if ( atomsOctree[ nodeID ].leaf || ( maxNodesInLevel >= atomsSubtreeRootsSize ) ) atomsSubtreeRoots[ maxSubtreeRoot++ ] = nodeID;
   else 
     {
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsOctree[ nodeID ].cPtr[ i ] >= 0 )
           fillAtomsSubtreeRootsArray( atomsOctree[ nodeID ].cPtr[ i ], 8 * maxNodesInLevel );
     }
}


bool fastBornRadius::buildAtomsOctree( void )
{  
   if ( printStatus ) printf( (char *)"\nbuilding atoms octree... " );
   
   double startT = getTime( );
   
   ATOM *atomsT;         
   atomsT = ( ATOM * ) malloc( nAtoms * sizeof( ATOM ) );
   
   if ( atomsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for atoms!" );
      if ( !atomsOctreeBuilt ) exit( 1 );
      return false;
     }

   countAtomsOctreeNodesAndSortAtoms( atoms, 0, nAtoms - 1, atomsT, &numAtomsOctreeNodes );
   freeMem( atomsT );
   
   ATOMS_OCTREE_NODE *atomsOctreeT;
   
   atomsOctreeT = ( ATOMS_OCTREE_NODE * ) malloc( numAtomsOctreeNodes * sizeof( ATOMS_OCTREE_NODE ) );

   if ( atomsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s atoms octree - memory allocation failed!", ( atomsOctreeBuilt ) ? (char *)"rebuild" : (char *)"build" );
      if ( !atomsOctreeBuilt ) exit( 1 );
      return false;
     }
 
   freeMem( atomsOctree );
   atomsOctree = atomsOctreeT;  

   initFreeNodeServer( nAtoms );
   
   atomsOctreeRoot = constructAtomsOctree( 0, nAtoms - 1, atoms );
   
   atomsOctreeBuilt = true;
   
   fillAtomsSubtreeRootsArray( atomsOctreeRoot, 1 );   

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}


bool fastBornRadius::buildOctrees( void )
{
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !qPointsOctreeBuilt ) ) buildQPointsOctree( );  
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !atomsOctreeBuilt ) ) buildAtomsOctree( );
   minRadiusUsed = minRadius;
   maxLeafSizeUsed = maxLeafSize;
}


void fastBornRadius::cleanupAtomsOctree( int nodeID )
{
   if ( nodeID < 0 ) return;
   
   atomsOctree[ nodeID ].invR = atomsOctree[ nodeID ].invCR = 0.0;

   if ( !atomsOctree[ nodeID ].leaf )
     {
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsOctree[ nodeID ].cPtr[ i ] >= 0 )
           cleanupAtomsOctree( atomsOctree[ nodeID ].cPtr[ i ] );
     }      
}


void fastBornRadius::cleanupPriorComputation( void )
{
   if ( priorComputationCleared ) return;

   if ( printStatus ) printf( (char *)"\ncleaning up prior computation... " );

   double startT = getTime( );
   
   for ( int i = 0; i < nAtoms; i++ )
     {
      atoms[ i ].invR = atoms[ i ].invCR = 0;
      atoms[ i ].R = 0;
     }
     
   if ( atomsOctreeBuilt ) cleanupAtomsOctree( atomsOctreeRoot );
   
   priorComputationCleared = true;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
}


void fastBornRadius::approximateIntegrals( int nodeA, int nodeQ )
{
   double sumRad = atomsOctree[ nodeA ].cr + qPointsOctree[ nodeQ ].cr;
   double sumRad2 = sumRad * sumRad;    
   double dx = qPointsOctree[ nodeQ ].cx - atomsOctree[ nodeA ].cx,
          dy = qPointsOctree[ nodeQ ].cy - atomsOctree[ nodeA ].cy,
          dz = qPointsOctree[ nodeQ ].cz - atomsOctree[ nodeA ].cz;
   double d2 = dx * dx + dy * dy + dz * dz;
   double wnx = qPointsOctree[ nodeQ ].wnx,
          wny = qPointsOctree[ nodeQ ].wny,
          wnz = qPointsOctree[ nodeQ ].wnz;                            
   bool farEnough = false;
   
   if ( ( sumRad2 > 0.25 ) && ( d2 > 4 * sumRad2 ) && ( d2 > ( sumRad2 / ( epsilon * epsilon ) ) ) ) farEnough = true;
   
   if ( farEnough )
      {
#ifdef BRR4      
       atomsOctree[ nodeA ].invR += ( ( wnx * dx + wny * dy + wnz * dz ) / ( d2 * d2 ) );
       atomsOctree[ nodeA ].invCR += ( ( wnx * dx + wny * dy + wnz * dz ) / ( d2 * d2 * d2 * sqrt( d2 ) ) );
#else       
       atomsOctree[ nodeA ].invR += ( ( wnx * dx + wny * dy + wnz * dz ) / ( d2 * d2 * d2 ) );
#endif       
      } 
   else
      {
       if ( atomsOctree[ nodeA ].leaf && qPointsOctree[ nodeQ ].leaf )
         {            
           for ( int i = atomsOctree[ nodeA ].atomsStartID; i <= atomsOctree[ nodeA ].atomsEndID; i++ )
              for ( int j = qPointsOctree[ nodeQ ].qPtsStartID; j <= qPointsOctree[ nodeQ ].qPtsEndID; j++ )
                {                            
                 dx = qPoints[ j ].x - atoms[ i ].x;
                 dy = qPoints[ j ].y - atoms[ i ].y;
                 dz = qPoints[ j ].z - atoms[ i ].z;                     
                        
                 double nx = qPoints[ j ].nx;
                 double ny = qPoints[ j ].ny;
                 double nz = qPoints[ j ].nz;                            
                        
                 double w = qPoints[ j ].w;
                        
                 d2 = dx * dx + dy * dy + dz * dz;       
                 
                 if ( d2 < atoms[ i ].r * atoms[ i ].r ) d2 = atoms[ i ].r * atoms[ i ].r;
                  
#ifdef BRR4                        
                 atoms[ i ].invR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 ) );            
                 atoms[ i ].invCR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 * d2 * sqrt( d2 ) ) );            
#else                 
                 atoms[ i ].invR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 * d2 ) );            
#endif                                  
                }                 
         }
       else if ( !atomsOctree[ nodeA ].leaf && !qPointsOctree[ nodeQ ].leaf )         
              {
                for ( int i = 0; i < 8; i++ )
                  if ( atomsOctree[ nodeA ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( qPointsOctree[ nodeQ ].cPtr[ j ] >= 0 ) 
                          approximateIntegrals( atomsOctree[ nodeA ].cPtr[ i ], qPointsOctree[ nodeQ ].cPtr[ j ] );                         
              }
            else if ( !atomsOctree[ nodeA ].leaf )         
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( atomsOctree[ nodeA ].cPtr[ i ] >= 0 )
                          approximateIntegrals( atomsOctree[ nodeA ].cPtr[ i ], nodeQ );                         
                   }
                 else 
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( qPointsOctree[ nodeQ ].cPtr[ j ] >= 0 )
                          approximateIntegrals( nodeA, qPointsOctree[ nodeQ ].cPtr[ j ] );                         
                   }              
      }      
}



void fastBornRadius::pushIntegralsToAtoms( int nodeID, double intVal, double intCVal )
{
   if ( nodeID < 0 ) return;

   intVal += atomsOctree[ nodeID ].invR;      
   intCVal += atomsOctree[ nodeID ].invCR;         
   atomsOctree[ nodeID ].invR = atomsOctree[ nodeID ].invCR = 0;
   
   if ( atomsOctree[ nodeID ].leaf )
     {
       for ( int i = atomsOctree[ nodeID ].atomsStartID; i <= atomsOctree[ nodeID ].atomsEndID; i++ )
         {                            
           atoms[ i ].invR += intVal;            
           atoms[ i ].invCR += intCVal;                       
           if ( atoms[ i ].invCR < 0.0 ) atoms[ i ].invCR = 0;
           double invR = COMBINE_INTEGRALS( atoms[ i ].invR, atoms[ i ].invCR );
           if ( invR > maxBornRadiusInv ) atoms[ i ].R = 1.0 / invR;
           else atoms[ i ].R = atoms[ i ].r;//maxBornRadius;
           if ( atoms[ i ].R < atoms[ i ].r ) atoms[ i ].R = atoms[ i ].r;
         }                       
     }
   else
     {
       for ( int i = 0; i < 8; i++ ) 
         if ( atomsOctree[ nodeID ].cPtr[ i ] >= 0 )
            pushIntegralsToAtoms( atomsOctree[ nodeID ].cPtr[ i ], intVal, intCVal );
     }  
}



void fastBornRadius::pushIntegralsAndFillAtomsSubtreeRootsArray( int nodeID, double intVal, double intCVal, int maxNodesInLevel )
{
   if ( nodeID < 0 ) return;
   
   if ( nodeID == atomsOctreeRoot ) maxSubtreeRoot = 0;
   
   if ( atomsOctree[ nodeID ].leaf || ( maxNodesInLevel >= atomsSubtreeRootsSize ) ) 
     {
       atomsOctree[ nodeID ].invR += intVal;
       atomsOctree[ nodeID ].invCR += intCVal;       
       atomsSubtreeRoots[ maxSubtreeRoot++ ] = nodeID;
      } 
   else 
     {
      intVal += atomsOctree[ nodeID ].invR;      
      intCVal += atomsOctree[ nodeID ].invCR;         
      atomsOctree[ nodeID ].invR = atomsOctree[ nodeID ].invCR = 0;
     
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsOctree[ nodeID ].cPtr[ i ] >= 0 )
           pushIntegralsAndFillAtomsSubtreeRootsArray( atomsOctree[ nodeID ].cPtr[ i ], intVal, intCVal, 8 * maxNodesInLevel );
     }
}


bool fastBornRadius::computeBornRadii( void )
{
   if ( !priorComputationCleared && ( minRadius == minRadiusUsed ) && ( maxLeafSize == maxLeafSizeUsed ) && ( epsilon == epsilonUsed ) ) return false;

   cleanupPriorComputation( );
   
   buildOctrees( );

   if ( printStatus ) printf( (char *)"\napproximating integrals and Born radii... " );

   double startT = getTime( );
   
   if ( numThreads == 1 ) 
      {
       approximateIntegrals( atomsOctreeRoot, qPointsOctreeRoot );
       pushIntegralsToAtoms( atomsOctreeRoot, 0, 0 );       
      } 
   else
      {
       initSubtreeRootServer( );
             
       pthread_t p[ numThreads ];           
       
       for ( int i = 0; i < numThreads; i++ )
          pthread_create( &p[ i ], NULL, approximateIntegralsThread, ( void * ) this );
          
       for ( int i = 0; i < numThreads; i++ )          
          pthread_join( p[ i ], NULL );
          
       pushIntegralsAndFillAtomsSubtreeRootsArray( atomsOctreeRoot, 0, 0, 1 );      
       
       initSubtreeRootServer( );          

       for ( int i = 0; i < numThreads; i++ )
          pthread_create( &p[ i ], NULL, pushIntegralsToAtomsThread, ( void * ) this );
          
       for ( int i = 0; i < numThreads; i++ )          
          pthread_join( p[ i ], NULL );
      }
      
   
   epsilonUsed = epsilon;
   priorComputationCleared = false;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   fastBRTime = endT - startT;   
      
   return true;
}



bool fastBornRadius::computeBornRadiiNaively( void )
{
   if ( computedBornRadiiNaively ) return false;
   
   if ( printStatus ) printf( (char *)"\ncomputing Born radii naively... " );

   double startT = getTime( );
   
   for ( int i = 0; i < nAtoms; i++ )
     {
      double invR = 0, invCR = 0;
      
      for ( int j = 0; j < nQPoints; j++ )
        {                            
         double dx = qPoints[ j ].x - atoms[ i ].x;
         double dy = qPoints[ j ].y - atoms[ i ].y;
         double dz = qPoints[ j ].z - atoms[ i ].z;                     
                
         double nx = qPoints[ j ].nx;
         double ny = qPoints[ j ].ny;
         double nz = qPoints[ j ].nz;                            
                
         double w = qPoints[ j ].w;
                
         double d2 = dx * dx + dy * dy + dz * dz;       

         if ( d2 < atoms[ i ].r * atoms[ i ].r ) d2 = atoms[ i ].r * atoms[ i ].r;
         
#ifdef BRR4                                  
         invR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 ) );            
         invCR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 * d2 * sqrt( d2 ) ) );                     
#else         
         invR += ( ( w * ( nx * dx + ny * dy + nz * dz ) ) / ( d2 * d2 * d2 ) );            
#endif                  
        }                 
        
      if ( invCR < 0.0 ) invCR = 0;
      
      invR = COMBINE_INTEGRALS( invR, invCR );
      
      if ( invR > maxBornRadiusInv ) atoms[ i ].naiveR = 1.0 / invR;
      else atoms[ i ].naiveR = maxBornRadius;                    
      
      if ( atoms[ i ].naiveR < atoms[ i ].r ) atoms[ i ].naiveR = atoms[ i ].r;      
     }   

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   naiveBRTime = endT - startT;
   
   return true;
}


bool fastBornRadius::getQPoints( int *numQPoints, double **qPts )
{
   if ( printStatus ) printf( (char *)"\ngetting quadrature points... " );
   
   double startT = getTime( );   
   
   *numQPoints = nQPoints;
   
   *qPts = ( double * ) malloc( 7 * nQPoints * sizeof( double ) );

   if ( *qPts == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }
        
   for ( int i = 0; i < nQPoints; i++ )
     {
       ( *qPts )[ 7 * i + 0 ] = qPoints[ i ].x;
       ( *qPts )[ 7 * i + 1 ] = qPoints[ i ].y;
       ( *qPts )[ 7 * i + 2 ] = qPoints[ i ].z;

       ( *qPts )[ 7 * i + 3 ] = qPoints[ i ].nx;
       ( *qPts )[ 7 * i + 4 ] = qPoints[ i ].ny;
       ( *qPts )[ 7 * i + 5 ] = qPoints[ i ].nz;

       ( *qPts )[ 7 * i + 6 ] = qPoints[ i ].w;
     } 

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, got %d quadrature points )\n", endT - startT, nQPoints );

   return true;     
}


bool fastBornRadius::getAtomsPQR( int *numAtoms, double **atomsPQR )
{
   if ( printStatus ) printf( (char *)"\ngetting atoms ( PQR )... " );
   
   double startT = getTime( );   
   
   *numAtoms = nAtoms;
   
   *atomsPQR = ( double * ) malloc( 5 * nAtoms * sizeof( double ) );

   if ( *atomsPQR == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }
        
   for ( int i = 0; i < nAtoms; i++ )
     {
       ( *atomsPQR )[ 5 * i + 0 ] = atoms[ i ].x;
       ( *atomsPQR )[ 5 * i + 1 ] = atoms[ i ].y;
       ( *atomsPQR )[ 5 * i + 2 ] = atoms[ i ].z;

       ( *atomsPQR )[ 5 * i + 3 ] = atoms[ i ].q;

       ( *atomsPQR )[ 5 * i + 4 ] = atoms[ i ].r;
     } 

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, got %d atoms )\n", endT - startT, nAtoms );

   return true;     
}



bool fastBornRadius::getAtomsPQRR( int *numAtoms, double **atomsPQRR, bool getNaiveR )
{
   if ( printStatus ) printf( (char *)"\ngetting atoms ( PQRR )... " );
   
   double startT = getTime( );   
   
   *numAtoms = nAtoms;
   
   *atomsPQRR = ( double * ) malloc( 6 * nAtoms * sizeof( double ) );

   if ( *atomsPQRR == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }
        
   for ( int i = 0; i < nAtoms; i++ )
     {
       ( *atomsPQRR )[ 6 * i + 0 ] = atoms[ i ].x;
       ( *atomsPQRR )[ 6 * i + 1 ] = atoms[ i ].y;
       ( *atomsPQRR )[ 6 * i + 2 ] = atoms[ i ].z;

       ( *atomsPQRR )[ 6 * i + 3 ] = atoms[ i ].q;

       ( *atomsPQRR )[ 6 * i + 4 ] = atoms[ i ].r;
       
       if ( getNaiveR ) ( *atomsPQRR )[ 6 * i + 5 ] = atoms[ i ].naiveR;       
       else ( *atomsPQRR )[ 6 * i + 5 ] = atoms[ i ].R;       
     } 

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, got %d atoms )\n", endT - startT, nAtoms );

   return true;     
}



bool fastBornRadius::getAtomsPQRR( int *numAtoms, double **atomsPQRR )
{
   return getAtomsPQRR( numAtoms, atomsPQRR, false );
}


double *fastBornRadius::getBornRadii( void )
{
   if ( printStatus ) printf( (char *)"\ncopying computed Born radii... " );
   
   double startT = getTime( );   
   
   if ( priorComputationCleared ) 
     {
      printError( (char *)"Need to (re)compute Born radii first!" );
      return NULL;
     } 
     
   if ( epsilon != epsilonUsed ) 
     {
      printError( (char *)"Born radii for epsilon = %lf have not yet been computed!", epsilon );
      return NULL;
     } 

   double *BornR;
   
   BornR = ( double * ) malloc( nAtoms * sizeof( double ) );

   if ( BornR == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return NULL;
     }
        
   for ( int i = 0; i < nAtoms; i++ )
      BornR[ atoms[ i ].id - 1 ] = atoms[ i ].R;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );

   return BornR;     
}



double *fastBornRadius::getNaiveBornRadii( void )
{
   if ( !computedBornRadiiNaively ) computeBornRadiiNaively( );

   if ( printStatus ) printf( (char *)"\ncopying naively computed Born radii... " );
   
   double startT = getTime( );   
   
   double *naiveBornR;
   
   naiveBornR = ( double * ) malloc( nAtoms * sizeof( double ) );

   if ( naiveBornR == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return NULL;
     }
        
   for ( int i = 0; i < nAtoms; i++ )
      naiveBornR[ atoms[ i ].id - 1 ] = atoms[ i ].naiveR;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );

   return naiveBornR;     
}



double *fastBornRadius::getIntegrals( void )
{
   if ( printStatus ) printf( (char *)"\ncopying computed integrals... " );
   
   double startT = getTime( );   

   if ( priorComputationCleared ) 
     {
      printError( (char *)"Need to (re)compute integrals first!" );
      return NULL;
     } 

   if ( epsilon != epsilonUsed ) 
     {
      printError( (char *)"Integrals for relative error %lf have not yet been computed!", epsilon );
      return NULL;
     } 
     
   double *integrals;
   
   integrals = ( double * ) malloc( 2 * nAtoms * sizeof( double ) );

   if ( integrals == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return NULL;
     }
        
   for ( int i = 0; i < nAtoms; i++ )
     {
      integrals[ 2 * atoms[ i ].id - 2 ] = atoms[ i ].invR / ( 4 * M_PI );     
      integrals[ 2 * atoms[ i ].id - 1 ] = atoms[ i ].invCR / ( 16 * M_PI );
     } 

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );

   return integrals;     
}


double fastBornRadius::getDispersionEnergy( void )
{
#ifdef BRR4
   return 0;
#else  
   double dispE = 0;
    
   for ( int i = 0; i < nAtoms; i++ )
      dispE += atoms[ i ].invR;     
     
   return dispE;  
#endif   
}


bool fastBornRadius::writeBornRadiiToFile( char *outFile )
{
   double *BornR;
   
   BornR = getBornRadii( );
   
   if ( BornR == NULL ) return false;
   
   if ( printStatus ) printf( (char *)"\nwriting computed Born radii to %s... ", outFile );   
   
   double startT = getTime( );   

   FILE *fp;
   
   fp = fopen( outFile, "wt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to create output file (%s)!", outFile );
      return false;
     }
     
   for ( int i = 0; i < nAtoms; i++ )
     fprintf( fp, (char *)"%lf\n", BornR[ i ] );
     
   fclose( fp );
     
   free( BornR );

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
        
   return true;     
}


bool fastBornRadius::writePQRRFile( char *outFile )
{
   if ( priorComputationCleared ) 
     {
      printError( (char *)"Need to (re)compute Born radii first!" );
      return false;
     } 
     
   if ( epsilon != epsilonUsed ) 
     {
      printError( (char *)"Born radii for epsilon = %lf have not yet been computed!", epsilon );
      return false;
     } 

   if ( printStatus ) printf( (char *)"\nwriting computed Born radii to %s... ", outFile );   
   
   double startT = getTime( );   

   FILE *fp;
   
   fp = fopen( outFile, (char *)"wt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to create output file (%s)!", outFile );
      return false;
     }
     
   fprintf( fp, (char *)"%d\n", nAtoms );
     
   for ( int i = 0; i < nAtoms; i++ )
     fprintf( fp, (char *)"%18.6lf %18.6lf %18.6lf %18.6lf %18.6lf %18.6lf\n", atoms[ i ].x, atoms[ i ].y, atoms[ i ].z, atoms[ i ].q, atoms[ i ].r, atoms[ i ].R );
     
   fclose( fp );
     
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
        
   return true;     
}


bool fastBornRadius::writeFastAndNaiveBornRadiiToFile( char *outFile )
{
   double *BornR;
   
   BornR = getBornRadii( );
   
   if ( BornR == NULL ) return false;

   double *naiveBornR;
   
   naiveBornR = getNaiveBornRadii( );
   
   if ( naiveBornR == NULL ) 
     {
      free( BornR );
      return false;
     } 
   
   if ( printStatus ) printf( (char *)"\nwriting both fast and naively computed Born radii to %s... ", outFile );   
   
   double startT = getTime( );   

   FILE *fp;
   
   fp = fopen( outFile, "wt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to create output file (%s)!", outFile );
      return false;
     }

   fprintf( fp, (char *)"Time for computing Born radii naively = %lf sec\n\n", naiveBRTime );

   fprintf( fp, (char *)"Time for computing Born radii fast = %lf sec ( speed-up factor = %0.2lf )\n", fastBRTime, ( fastBRTime != 0 ) ? ( naiveBRTime / fastBRTime ) : 0 );
   fprintf( fp, (char *)"( Parameters: minRadius = %lf, maxLeafSize = %d, epsilon = %lf, numThreads = %d )\n\n", minRadiusUsed, maxLeafSizeUsed, epsilonUsed, numThreads );
  
   fprintf( fp, (char *)"        FAST         NAIVE  (       ERROR )\n" );
        
   double sumErrorPercent = 0;
       
   for ( int i = 0; i < nAtoms; i++ )
     {
#if defined(__APPLE__)
      double errorPercent = ( naiveBornR[ i ] != 0 ) ? ( ( ( BornR[ i ] - naiveBornR[ i ] ) / naiveBornR[ i ] ) * 100 ) : ( ( BornR[ i ] == 0 ) ? 0 : 1000000000 );
#else
      double errorPercent = ( naiveBornR[ i ] != 0 ) ? ( ( ( BornR[ i ] - naiveBornR[ i ] ) / naiveBornR[ i ] ) * 100 ) : ( ( BornR[ i ] == 0 ) ? 0 : 10000000000 );
#endif

      fprintf( fp, (char *)"%12.6lf  %12.6lf  ( %9.3lf \% )\n", BornR[ i ], naiveBornR[ i ], errorPercent );
      sumErrorPercent += fabs( errorPercent );
     } 
     
   fprintf( fp, (char *)"\nAverage Absolute Error = %0.3lf \%\n", sumErrorPercent / nAtoms );
     
   fclose( fp );

   free( BornR );     
   free( naiveBornR );

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
        
   return true;     
}



bool fastBornRadius::writeIntegralsToFile( char *outFile )
{
   double *integrals;
   
   integrals = getIntegrals( );
   
   if ( integrals == NULL ) return false;

   if ( printStatus ) printf( (char *)"\nwriting computed integrals to %s... ", outFile );
   
   double startT = getTime( );   

   FILE *fp;
   
   fp = fopen( outFile, "wt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to create output file (%s)!", outFile );
      return false;
     }

   fprintf( fp, (char *)"    ORIGINAL    CORRECTION\n" );
     
   for ( int i = 0; i < nAtoms; i++ )   
     fprintf( fp, (char *)"%12.6lf  %12.6lf\n", integrals[ 2 * i ], integrals[ 2 * i + 1 ] );
     
   fclose( fp );
     
   free( integrals );

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
        
   return true;     
}
