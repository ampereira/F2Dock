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

#include "clashFilter.h"


void clashFilter::printError( char *format, ... )
{
   char eMsg[ 500 ];
   va_list args;
   
   va_start( args, format );
   
   vsprintf( eMsg, format, args );
   
   va_end( args );
   
   printf( (char *)"\nError: %s\n\n", eMsg );   
}


double clashFilter::getTime( void )
{
#ifdef _WIN32
   time_t ltime;
   _timeb tstruct;
   time( &ltime );
   _ftime( &tstruct );
   return ( double ) ( ltime + 1e-3 * ( tstruct.millitm ) );
#else
   struct timeval t;
   gettimeofday( &t, NULL );
   return ( double )( t.tv_sec + 1e-6 * t.tv_usec );
#endif
}


void clashFilter::freeMemory( void )
{
   freeMem( staticAtoms );
   freeMem( movingAtoms );
   freeMem( staticAtomsOctree );
   freeMem( movingAtomsOctree );  
   freeMem( movingAtomsSubtreeRoots ); 
}


bool clashFilter::allocateMovingAtomsSubtreeRootsArray( int nThreads )
{
   int arraySize = 1;
   
   while ( arraySize < 8 * nThreads )
      arraySize *= 8;
      
   movingAtomsSubtreeRoots = ( int * ) malloc( arraySize * sizeof( int ) );      
   
   if ( movingAtomsSubtreeRoots == NULL )
     {
      printError( (char *)"Failed to allocate memory!" );
      return false;
     }
     
   movingAtomsSubtreeRootsSize = arraySize;  
        
   return true;  
}


void clashFilter::setDefaults( void )
{
   staticAtoms = NULL;
   movingAtoms = NULL;
      
   staticAtomsOctree = NULL;
   movingAtomsOctree = NULL;
   
   movingAtomsSubtreeRoots = NULL;
   
   minRadius = 2.0;
   minRadiusUsed = -1;   
   
   maxLeafSize = 10;
   maxLeafSizeUsed = -1;   

   epsilon = 0.5;
   
   clashFrac = 0.4;
   severeClashFrac = 0.25;
   fuzzyFrac = 1.0;     
      
   staticAtomsOctreeBuilt = false;
   movingAtomsOctreeBuilt = false;
   
   priorComputationCleared = true;

   curNode = maxNode = 0;
   curSubtreeRoot = maxSubtreeRoot = 0;   
   
   staticAtomsOctreeRoot = -1;
   movingAtomsOctreeRoot = -1;
   
   numThreads = 1;

   transMatrix.reset( );
      
   printStatus = true;
   
   pthread_mutex_init( &nodesLock, NULL );   
}


void clashFilter::printCurrentSettings( void )
{
   printf( (char *)"\nCurrent Parameter Settings:\n" );
   printf( (char *)"\tminRadius = %lf, maxLeafSize = %d, epsilon = %lf, numThreads = %d\n", minRadius, maxLeafSize, epsilon, numThreads );
   printf( (char *)"\tclasfFrac = %lf, severeClashFrac = %lf, fuzzyFrac = %lf\n\n", clashFrac, severeClashFrac, fuzzyFrac );   
}


clashFilter::clashFilter( int numStaticAtoms, double *stAtoms, int numMovingAtoms, double *mvAtoms, bool printStat )
{
   setDefaults( );

   printStatus = printStat;
         
   if ( !copyAtomsFromArray( numStaticAtoms, stAtoms, &nStaticAtoms, &staticAtoms ) 
     || !copyAtomsFromArray( numMovingAtoms, mvAtoms, &nMovingAtoms, &movingAtoms ) )
      {
       freeMemory( );
       exit( 1 );
      }
            
   if ( !allocateMovingAtomsSubtreeRootsArray( numThreads ) )
     {
      freeMemory( );
      exit( 1 );
     }
 
   buildOctrees( );
  
   if ( printStatus ) printCurrentSettings( );
}


clashFilter::clashFilter( int numStaticAtoms, double *stAtoms, int numMovingAtoms, double *mvAtoms )
{
   clashFilter( numStaticAtoms, stAtoms, numMovingAtoms, mvAtoms, true );
}



clashFilter::~clashFilter( )
{
   freeMemory( );
   if ( printStatus ) printf( (char *)"\n" );
}


bool clashFilter::setProximityFactors( double clashFactor, double severeClashFactor, double fuzzyFactor )
{
   if ( clashFactor < 0 )
     {
      printError( (char *)"clashFactor must be a non-negative real number!" );
      return false;     
     }
     
   clashFrac = clashFactor;
   if ( printStatus ) printf( (char *)"\nclashFactor is set to %lf\n", clashFactor );

   if ( severeClashFactor < 0 )
     {
      printError( (char *)"severeClashFactor must be a non-negative real number!" );
      return false;     
     }

   if ( severeClashFactor > clashFactor )
     {
      printError( (char *)"severeClashFactor cannot be larger than clashFactor!" );
      return false;     
     }
     
   severeClashFrac = severeClashFactor;
   if ( printStatus ) printf( (char *)"\nsevereClashFactor is set to %lf\n", severeClashFactor );
     
   if ( fuzzyFactor < 0 )
     {
      printError( (char *)"fuzzyFactor must be a non-negative real number!" );
      return false;     
     }
     
   fuzzyFrac = fuzzyFactor;
   if ( printStatus ) printf( (char *)"\nfuzzyFactor is set to %lf\n", fuzzyFactor );   
   
   return true;
}


bool clashFilter::setMinRadius( double minRad )
{
   if ( minRad < 0 )
     {
      printError( (char *)"minRadius must be a non-negative real number!" );
      return false;     
     }
     
   minRadius = minRad;

   buildOctrees( );

   if ( printStatus ) printf( (char *)"\nminRadius is set to %lf\n", minRad );
   
   return true;
}


bool clashFilter::setMaxLeafSize( int maxLfSize )
{
   if ( maxLfSize <= 0 )
     {
      printError( (char *)"maxLeafSize must be a positive integer!" );
      return false;     
     }
     
   maxLeafSize = maxLfSize;

   buildOctrees( );

   if ( printStatus ) printf( (char *)"\nmaxLeafSize is set to %d\n", maxLfSize );

   return true;
}


bool clashFilter::setEpsilon( double eps )
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


bool clashFilter::setNumThreads( int nThreads )
{
   if ( nThreads < 1 )
     {
      printError( (char *)"numThreads must be a positive integer!" );
      return false;     
     }
     
   int *movingAtomsSubtreeRootsT;  
   
   movingAtomsSubtreeRootsT = movingAtomsSubtreeRoots;
     
   if ( allocateMovingAtomsSubtreeRootsArray( nThreads ) )
     {
      numThreads = nThreads;
      freeMem( movingAtomsSubtreeRootsT );
      if ( movingAtomsOctreeBuilt ) fillMovingAtomsSubtreeRootsArray( movingAtomsOctreeRoot, 1 );   
      if ( printStatus ) printf( (char *)"\nnumThreads is set to %d\n", numThreads );
      return true;
     }
   else
     {
      movingAtomsSubtreeRoots = movingAtomsSubtreeRootsT;
      if ( printStatus ) printf( (char *)"\nnumThreads remains unchanged ( %d )\n", numThreads );      
      return false;
     }            
}



bool clashFilter::setTransformationMatrix( Matrix transMat )
{
   transMatrix.set( transMat );
   
   if ( printStatus ) 
     {
       printf( (char *)"\nTransformation Matrix is set to:\n" );
       transMatrix.print( );
     }  
   
   return true;
}


void clashFilter::setPrintStatus( bool printStat )
{
   printStatus = printStat;
   
   if ( printStatus ) printf( (char *)"\nprintStatus is set to true\n" );
}


int clashFilter::nextFreeNode( void )
{
   int nextNode = -1;
   
   pthread_mutex_lock( &nodesLock );
   if ( curNode < maxNode ) nextNode = curNode++;
   pthread_mutex_unlock( &nodesLock );
   
   return nextNode;
}


void clashFilter::initFreeNodeServer( int numNodes )
{
//   pthread_mutex_init( &nodesLock, NULL );
   curNode = 0;
   maxNode = numNodes;
}


int clashFilter::nextSubtreeRoot( void )
{
   int nextRoot = -1;
   
   pthread_mutex_lock( &subtreeRootsLock );
   if ( curSubtreeRoot < maxSubtreeRoot ) nextRoot = curSubtreeRoot++;
   pthread_mutex_unlock( &subtreeRootsLock );
   
   return ( nextRoot == -1 ) ? ( - 1 ) : movingAtomsSubtreeRoots[ nextRoot ];
}


void clashFilter::initSubtreeRootServer( void )
{
   pthread_mutex_init( &subtreeRootsLock, NULL );
   curSubtreeRoot = 0;
}


bool clashFilter::copyAtomsFromArray( int numAtomsSrc, double *atmsSrc, int *numAtomsDest, ATOM **atmsDest )
{
   if ( printStatus ) printf( (char *)"\ncopying atoms from array... " );

   double startT = getTime( );
   
   if ( numAtomsSrc <= 0 )
     {
      printError( (char *)"No atoms to copy!" );
      return false;
     }
   
   *numAtomsDest = numAtomsSrc;
   ( *atmsDest ) = ( ATOM * ) malloc( numAtomsSrc * sizeof( ATOM ) );
   
   if ( ( *atmsDest ) == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }

   for ( int i = 0; i < numAtomsSrc; i++ )   
     {          
      ( *atmsDest )[ i ].x = atmsSrc[ 5 * i + 0 ];
      ( *atmsDest )[ i ].y = atmsSrc[ 5 * i + 1 ];
      ( *atmsDest )[ i ].z = atmsSrc[ 5 * i + 2 ];

      ( *atmsDest )[ i ].q = atmsSrc[ 5 * i + 3 ];

      ( *atmsDest )[ i ].r = atmsSrc[ 5 * i + 4 ];
     }    
     
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec, copied %d atoms )\n", endT - startT, numAtomsSrc );
         
   return true;
}


void clashFilter::countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID, 
                                                     ATOM *atomsT, int *numNodes )
{
   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z, maxR = atoms[ atomsStartID ].r;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      
      
      if ( atoms[ i ].r > maxR ) maxR = atoms[ i ].r;            
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

   double cr = sqrt( r2 ) + maxR;
      
   *numNodes = 1;
      
   if ( ( atomsEndID - atomsStartID +  1 > maxLeafSize ) && ( cr > minRadius ) )
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



int clashFilter::constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms, ATOMS_OCTREE_NODE *atomsOctree )
{
   int nodeID = nextFreeNode( );
      
   atomsOctree[ nodeID ].atomsStartID = atomsStartID;
   atomsOctree[ nodeID ].atomsEndID = atomsEndID;
   
   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z, maxR = atoms[ atomsStartID ].r;
   double qSum = 0;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      

      if ( atoms[ i ].r > maxR ) maxR = atoms[ i ].r;      
      
      qSum += atoms[ i ].q;
     } 
   
   double cx = atomsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = atomsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = atomsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;
   
   double cq = atomsOctree[ nodeID ].cq = qSum;

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
   
   double cr = atomsOctree[ nodeID ].cr = sqrt( r2 ) + maxR;
         
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
           int j = constructAtomsOctree( atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atoms, atomsOctree );
           atomsOctree[ nodeID ].cPtr[ i ] = j; 
          }
        else atomsOctree[ nodeID ].cPtr[ i ] = -1;           
     }  
     
   return nodeID;  
}


void clashFilter::fillMovingAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel )
{
   if ( nodeID < 0 ) return;
   
   if ( nodeID == movingAtomsOctreeRoot ) maxSubtreeRoot = 0;
   
   if ( movingAtomsOctree[ nodeID ].leaf || ( maxNodesInLevel >= movingAtomsSubtreeRootsSize ) ) movingAtomsSubtreeRoots[ maxSubtreeRoot++ ] = nodeID;
   else 
     {
      for ( int i = 0; i < 8; i++ ) 
        if ( movingAtomsOctree[ nodeID ].cPtr[ i ] >= 0 )
           fillMovingAtomsSubtreeRootsArray( movingAtomsOctree[ nodeID ].cPtr[ i ], 8 * maxNodesInLevel );
     }
}


bool clashFilter::buildStaticAtomsOctree( void )
{  
   if ( printStatus ) printf( (char *)"\nbuilding static atoms octree... " );
   
   double startT = getTime( );
   
   ATOM *atomsT;         
   atomsT = ( ATOM * ) malloc( nStaticAtoms * sizeof( ATOM ) );
   
   if ( atomsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for static atoms!" );
      if ( !staticAtomsOctreeBuilt ) exit( 1 );
      return false;
     }

   countAtomsOctreeNodesAndSortAtoms( staticAtoms, 0, nStaticAtoms - 1, atomsT, &numStaticAtomsOctreeNodes );
   freeMem( atomsT );
   
   ATOMS_OCTREE_NODE *atomsOctreeT;
   
   atomsOctreeT = ( ATOMS_OCTREE_NODE * ) malloc( numStaticAtomsOctreeNodes * sizeof( ATOMS_OCTREE_NODE ) );

   if ( atomsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s static atoms octree - memory allocation failed!", ( staticAtomsOctreeBuilt ) ? (char *)"rebuild" : (char *)"build" );
      if ( !staticAtomsOctreeBuilt ) exit( 1 );
      return false;
     }
 
   freeMem( staticAtomsOctree );
   staticAtomsOctree = atomsOctreeT; 

   initFreeNodeServer( nStaticAtoms );
   
   staticAtomsOctreeRoot = constructAtomsOctree( 0, nStaticAtoms - 1, staticAtoms, staticAtomsOctree );

   staticAtomsOctreeBuilt = true;

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}



bool clashFilter::buildMovingAtomsOctree( void )
{  
   if ( printStatus ) printf( (char *)"\nbuilding moving atoms octree... " );
   
   double startT = getTime( );
   
   ATOM *atomsT;         
   atomsT = ( ATOM * ) malloc( nMovingAtoms * sizeof( ATOM ) );
   
   if ( atomsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for moving atoms!" );
      if ( !movingAtomsOctreeBuilt ) exit( 1 );
      return false;
     }

   countAtomsOctreeNodesAndSortAtoms( movingAtoms, 0, nMovingAtoms - 1, atomsT, &numMovingAtomsOctreeNodes );
   freeMem( atomsT );
   
   ATOMS_OCTREE_NODE *atomsOctreeT;
   
   atomsOctreeT = ( ATOMS_OCTREE_NODE * ) malloc( numMovingAtomsOctreeNodes * sizeof( ATOMS_OCTREE_NODE ) );

   if ( atomsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s moving atoms octree - memory allocation failed!", ( movingAtomsOctreeBuilt ) ? (char *)"rebuild" : (char *)"build" );
      if ( !movingAtomsOctreeBuilt ) exit( 1 );
      return false;
     }
 
   freeMem( movingAtomsOctree );
   movingAtomsOctree = atomsOctreeT;  

   initFreeNodeServer( nMovingAtoms );
   
   movingAtomsOctreeRoot = constructAtomsOctree( 0, nMovingAtoms - 1, movingAtoms, movingAtomsOctree );
   
   movingAtomsOctreeBuilt = true;
   
   fillMovingAtomsSubtreeRootsArray( movingAtomsOctreeRoot, 1 );   

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}


bool clashFilter::buildOctrees( void )
{
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !staticAtomsOctreeBuilt ) ) buildStaticAtomsOctree( );  
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !movingAtomsOctreeBuilt ) ) buildMovingAtomsOctree( );
   minRadiusUsed = minRadius;
   maxLeafSizeUsed = maxLeafSize;
}


void clashFilter::approximateInteractions( Matrix transMat, int nodeS, int nodeM, int *nClashes, int *nSevereClashes, double *interactionValue )
{
   double sumRad = staticAtomsOctree[ nodeS ].cr + movingAtomsOctree[ nodeM ].cr;
   double sumRad2 = sumRad * sumRad;   
   
   Vector oldPos( movingAtomsOctree[ nodeM ].cx, movingAtomsOctree[ nodeM ].cy, movingAtomsOctree[ nodeM ].cz, 1.0 );
   Vector newPos = transMat * oldPos;          
   double cxM = newPos[ 0 ], cyM = newPos[ 1 ], czM = newPos[ 2 ];      
    
   double dx = staticAtomsOctree[ nodeS ].cx - cxM,
          dy = staticAtomsOctree[ nodeS ].cy - cyM,
          dz = staticAtomsOctree[ nodeS ].cz - czM;
   double d2 = dx * dx + dy * dy + dz * dz;
   
   bool farEnough = false;
   
   if ( ( d2 > 4 * sumRad2 ) && ( d2 > ( sumRad2 / ( epsilon * epsilon ) ) ) ) farEnough = true;
   
   *nClashes = *nSevereClashes = 0;
   *interactionValue = 0;
   
   if ( farEnough ) 
      {
        double qiqj = staticAtomsOctree[ nodeS ].cq * movingAtomsOctree[ nodeM ].cq;
        *interactionValue = qiqj / d2; 
      }  
   else
      {
       if ( staticAtomsOctree[ nodeS ].leaf && movingAtomsOctree[ nodeM ].leaf )
         {            
           for ( int i = movingAtomsOctree[ nodeM ].atomsStartID; i <= movingAtomsOctree[ nodeM ].atomsEndID; i++ )
             {
              Vector oldPos( movingAtoms[ i ].x, movingAtoms[ i ].y, movingAtoms[ i ].z, 1.0 );
              Vector newPos = transMat * oldPos;          
              double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];      
              
              bool clash = false, severeClash = false;
             
              for ( int j = staticAtomsOctree[ nodeS ].atomsStartID; j <= staticAtomsOctree[ nodeS ].atomsEndID; j++ )
                {                            
                 dx = staticAtoms[ j ].x - xM;
                 dy = staticAtoms[ j ].y - yM;
                 dz = staticAtoms[ j ].z - zM;                     
                                                
                 d2 = dx * dx + dy * dy + dz * dz;       

                 double rSum = movingAtoms[ i ].r + staticAtoms[ j ].r;
                 double qiqj = movingAtoms[ i ].q * staticAtoms[ j ].q;
         
                 if ( d2 <= ( clashFrac * rSum ) * ( clashFrac * rSum ) ) clash = true; //( *nClashes )++;
                 if ( d2 <= ( severeClashFrac * rSum ) * ( severeClashFrac * rSum ) ) severeClash = true; //( *nSevereClashes )++;
         
                 if ( d2 < ( fuzzyFrac * rSum ) * ( fuzzyFrac * rSum ) ) d2 = ( fuzzyFrac * rSum ) * ( fuzzyFrac * rSum );     
                
                 ( *interactionValue ) += ( qiqj / d2 );                         
                }                 
                
              if ( clash ) ( *nClashes )++;
              if ( severeClash ) ( *nSevereClashes )++;
             }   
         }
       else if ( !staticAtomsOctree[ nodeS ].leaf && !movingAtomsOctree[ nodeM ].leaf )         
              {
                for ( int i = 0; i < 8; i++ )
                  if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 ) 
                         {
                           int nC, nSC;
                           double intVal;
                           
                           approximateInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], movingAtomsOctree[ nodeM ].cPtr[ j ], &nC, &nSC, &intVal );
                           
                           ( *nClashes ) += nC;
                           ( *nSevereClashes ) += nSC;
                           ( *interactionValue ) += intVal;                         
                         }  
              }
            else if ( !staticAtomsOctree[ nodeS ].leaf )         
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                         {
                           int nC, nSC;
                           double intVal;
                           
                           approximateInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], nodeM, &nC, &nSC, &intVal );
                           
                           ( *nClashes ) += nC;
                           ( *nSevereClashes ) += nSC;
                           ( *interactionValue ) += intVal;
                         }  
                   }
                 else 
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 )
                         {
                           int nC, nSC;
                           double intVal;
                           
                           approximateInteractions( transMat, nodeS, movingAtomsOctree[ nodeM ].cPtr[ j ], &nC, &nSC, &intVal );
                           
                           ( *nClashes ) += nC;
                           ( *nSevereClashes ) += nSC;
                           ( *interactionValue ) += intVal;
                         }  
                   }              
      }      
}



bool clashFilter::computeInteractions( Matrix transMat, int *nClashes, int *nSevereClashes, double *interactionValue )
{
   if ( !staticAtomsOctreeBuilt || !movingAtomsOctreeBuilt ) return false;

   if ( printStatus ) printf( (char *)"\napproximating interactions and computing clashes... " );

   double startT = getTime( );
   
   if ( numThreads == 1 ) approximateInteractions( transMat, staticAtomsOctreeRoot, movingAtomsOctreeRoot, nClashes, nSevereClashes, interactionValue );
   else
      {
       initSubtreeRootServer( );
       initFreeNodeServer( numThreads );      
             
       pthread_t p[ numThreads ];         
       THREAD_RESULT threadResults[ numThreads ]; 
       
       for ( int i = 0; i < numThreads; i++ )
         {
           threadResults[ i ].cF = this;
           threadResults[ i ].transMat = &transMat;
           pthread_create( &p[ i ], NULL, approximateInteractionsThread, ( void * ) &threadResults[ i ] );
         }  
          
       for ( int i = 0; i < numThreads; i++ )          
          pthread_join( p[ i ], NULL );          
          
       *nClashes = *nSevereClashes = 0;
       *interactionValue = 0;
       
       for ( int i = 0; i < numThreads; i++ )          
         {
           ( *nClashes ) += threadResults[ i ].nClashes;
           ( *nSevereClashes ) += threadResults[ i ].nSevereClashes;           
           ( *interactionValue ) += threadResults[ i ].interactionValue;           
         }  
      }
         
   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}



bool clashFilter::computeInteractions( int *nClashes, int *nSevereClashes, double *interactionValue )
{
   return computeInteractions( transMatrix, nClashes, nSevereClashes, interactionValue );
}


bool clashFilter::computeInteractionsNaively( Matrix transMat, int *nClashes, int *nSevereClashes, double *interactionValue )
{
   if ( printStatus ) printf( (char *)"\ncomputing steric interactions naively... " );

   double startT = getTime( );
   
   *nClashes = *nSevereClashes = 0;
   *interactionValue = 0.0; 
   
//   printf( "\nnStaticAtoms = %d\n", nStaticAtoms );
//
//   for ( int i = 0; i < nStaticAtoms; i++ )
//     printf( "\n< %lf, %lf, %lf, %lf, %lf >\n", staticAtoms[ i ].x, staticAtoms[ i ].y, staticAtoms[ i ].z, staticAtoms[ i ].q, staticAtoms[ i ].r );
//   
//   printf( "\nnMovingAtoms = %d\n", nMovingAtoms );
//
//   for ( int i = 0; i < nMovingAtoms; i++ )
//     printf( "\n< %lf, %lf, %lf, %lf, %lf >\n", movingAtoms[ i ].x, movingAtoms[ i ].y, movingAtoms[ i ].z, movingAtoms[ i ].q, movingAtoms[ i ].r );
//   
//   printf( "\n\n" );

   double frac = 1.5;
   double step = 1.5;
   
   int numClose = 0;
   for ( int i = 0; i < nMovingAtoms; i++ )
     {
      Vector oldPos( movingAtoms[ i ].x, movingAtoms[ i ].y, movingAtoms[ i ].z, 1.0 );
      Vector newPos = transMat * oldPos;
          
      double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];      
     
      int k = 0;
      
      for ( int j = 0; j < nStaticAtoms; j++ )
        {                            
         double dx = staticAtoms[ j ].x - xM;
         double dy = staticAtoms[ j ].y - yM;
         double dz = staticAtoms[ j ].z - zM;                     
                
         double d2 = dx * dx + dy * dy + dz * dz;
         double rSum = movingAtoms[ i ].r + staticAtoms[ j ].r;
         double qiqj = movingAtoms[ i ].q * staticAtoms[ j ].q;
         
         if ( d2 <= ( clashFrac * rSum ) * ( clashFrac * rSum ) ) ( *nClashes )++;
         if ( d2 <= ( severeClashFrac * rSum ) * ( severeClashFrac * rSum ) ) ( *nSevereClashes )++;
         
//         if ( d2 < ( fuzzyFrac * rSum ) * ( fuzzyFrac * rSum ) ) d2 = ( fuzzyFrac * rSum ) * ( fuzzyFrac * rSum );     

         double minD2 = ( frac * rSum ) * ( frac * rSum );

         if ( d2 < ( 1 / ( frac * frac ) ) * minD2 ) k++;
         
         if ( d2 < minD2 ) d2 = minD2; 
         
//         d2 = ceil( d2 / ( step * step ) );
//         d2 = ( step * step ) * d2;
               
         ( *interactionValue ) += ( qiqj / d2 );       
//         ( *interactionValue ) += ( minD2 / d2 );// ( ( minD2 * minD2 * minD2 ) / ( d2 * d2 * d2 ) );       
        }           
        
      if ( k > 0 ) numClose++;                
     }   
     
//  ( *interactionValue ) /= numClose;  

   double endT = getTime( );
   
   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
   
   return true;
}


bool clashFilter::computeInteractionsNaively( int *nClashes, int *nSevereClashes, double *interactionValue )
{
   return computeInteractionsNaively( transMatrix, nClashes, nSevereClashes, interactionValue );
}
