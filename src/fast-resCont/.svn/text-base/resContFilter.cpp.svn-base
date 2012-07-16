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

#include "resContFilter.h"


void resContFilter::printError( char *format, ... )
{
   char eMsg[ 500 ];
   va_list args;

   va_start( args, format );

   vsprintf( eMsg, format, args );

   va_end( args );

   printf( "\nError: %s\n\n", eMsg );
}


double resContFilter::getTime( void )
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


void resContFilter::freeMemory( void )
{
   freeMem( staticAtoms );
   freeMem( movingAtoms );
   freeMem( staticAtomsOctree );
   freeMem( movingAtomsOctree );
   freeMem( movingAtomsSubtreeRoots );
}


bool resContFilter::allocateMovingAtomsSubtreeRootsArray( int nThreads )
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


void resContFilter::setDefaults( void )
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

   distCutoff = 10.0;
   distCutoffSQ = distCutoff * distCutoff;
   distCutoffUsed = -1;

   staticAtomsOctreeBuilt = false;
   movingAtomsOctreeBuilt = false;

   priorComputationCleared = true;

   curNode = maxNode = 0;
   curSubtreeRoot = maxSubtreeRoot = 0;

   staticAtomsOctreeRoot = -1;
   movingAtomsOctreeRoot = -1;

   interactionMatrixRead = false;

   numThreads = 1;

   transMatrix.reset( );

   printStatus = true;

   pthread_mutex_init( &nodesLock, NULL );
}


void resContFilter::printCurrentSettings( void )
{
   printf( "\nCurrent Parameter Settings:\n" );
   printf( "\tminRadius = %lf, maxLeafSize = %d, distCutoff = %lf, numThreads = %d\n", minRadius, maxLeafSize, distCutoff, numThreads );
}


resContFilter::resContFilter( int numStaticAtoms, double *stAtoms, int numMovingAtoms, double *mvAtoms, char *intFile, bool printStat )
{
   setDefaults( );

   if ( intFile == NULL )
     {
       for ( int i = 0; i < NUM_RESIDUE_TYPES; i++ )
         for ( int j = 0; j < NUM_RESIDUE_TYPES; j++ )
           intVal[ i ][ j ] = intValDefault[ i ][ j ];
     }
   else if ( !readInteractionMatrix( intFile ) ) exit( 1 );

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



resContFilter::~resContFilter( )
{
   freeMemory( );
   if ( printStatus ) printf( "\n" );
}


bool resContFilter::readInteractionMatrix( char *intFile )
{
   FILE *fp;

   fp = fopen( intFile, "rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open interaction matrix file (%s)!", intFile );
        return false;
     }

   interactionMatrixRead = false;

   for ( int i = 0; i < NUM_RESIDUE_TYPES; i++ )
    {
      for ( int j = 0; j < NUM_RESIDUE_TYPES; j++ )
        {
          double v;

          if ( fscanf( fp, "%lf", &v ) != 1 )
            {
              fclose( fp );
              return false;
            }

          intVal[ i ][ j ] = v;
        }        
    }    
    
   fflush( stdout ); 

   fclose( fp );

   interactionMatrixRead = true;

   return true;
}


bool resContFilter::setMinRadius( double minRad )
{
   if ( minRad < 0 )
     {
      printError( (char *)"minRadius must be a non-negative real number!" );
      return false;
     }

   minRadius = minRad;

   buildOctrees( );

   if ( printStatus ) printf( "\nminRadius is set to %lf\n", minRad );

   return true;
}


bool resContFilter::setMaxLeafSize( int maxLfSize )
{
   if ( maxLfSize <= 0 )
     {
      printError( (char *)"maxLeafSize must be a positive integer!" );
      return false;
     }

   maxLeafSize = maxLfSize;

   buildOctrees( );

   if ( printStatus ) printf( "\nmaxLeafSize is set to %d\n", maxLfSize );

   return true;
}


bool resContFilter::setDistanceCutoff( double distanceCutoff )
{
   if ( distanceCutoff < 0 )
     {
      printError( (char *)"distCutoff must be a non-negative real number!" );
      return false;
     }

   distCutoff = distanceCutoff;
   distCutoffSQ = distCutoff * distCutoff;

   buildOctrees( );

   if ( printStatus ) printf( "\ndistCutoff is set to %lf\n", distanceCutoff );

   return true;
}



bool resContFilter::setNumThreads( int nThreads )
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
      if ( printStatus ) printf( "\nnumThreads is set to %d\n", numThreads );
      return true;
     }
   else
     {
      movingAtomsSubtreeRoots = movingAtomsSubtreeRootsT;
      if ( printStatus ) printf( "\nnumThreads remains unchanged ( %d )\n", numThreads );
      return false;
     }
}



bool resContFilter::setTransformationMatrix( Matrix transMat )
{
   transMatrix.set( transMat );

   if ( printStatus )
     {
       printf( "\nTransformation Matrix is set to:\n" );
       transMatrix.print( );
     }

   return true;
}


void resContFilter::setPrintStatus( bool printStat )
{
   printStatus = printStat;

   if ( printStatus ) printf( "\nprintStatus is set to true\n" );
}


int resContFilter::nextFreeNode( void )
{
   int nextNode = -1;

   pthread_mutex_lock( &nodesLock );
   if ( curNode < maxNode ) nextNode = curNode++;
   pthread_mutex_unlock( &nodesLock );

   return nextNode;
}


void resContFilter::initFreeNodeServer( int numNodes )
{
//   pthread_mutex_init( &nodesLock, NULL );
   curNode = 0;
   maxNode = numNodes;
}


int resContFilter::nextSubtreeRoot( void )
{
   int nextRoot = -1;

   pthread_mutex_lock( &subtreeRootsLock );
   if ( curSubtreeRoot < maxSubtreeRoot ) nextRoot = curSubtreeRoot++;
   pthread_mutex_unlock( &subtreeRootsLock );

   return ( nextRoot == -1 ) ? ( - 1 ) : movingAtomsSubtreeRoots[ nextRoot ];
}


void resContFilter::initSubtreeRootServer( void )
{
   pthread_mutex_init( &subtreeRootsLock, NULL );
   curSubtreeRoot = 0;
}


bool resContFilter::copyAtomsFromArray( int numAtomsSrc, double *atmsSrc, int *numAtomsDest, ATOM **atmsDest )
{
   if ( printStatus ) printf( "\ncopying atoms from array... " );

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

      ( *atmsDest )[ i ].rID = atmsSrc[ 5 * i + 3 ];

      ( *atmsDest )[ i ].rNum = atmsSrc[ 5 * i + 4 ];
     }

   double endT = getTime( );

   if ( printStatus ) printf( "done ( %lf sec, copied %d atoms )\n", endT - startT, numAtomsSrc );

   return true;
}


void resContFilter::countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID,
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

   double cr = sqrt( r2 );

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



int resContFilter::constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms, ATOMS_OCTREE_NODE *atomsOctree )
{
   int nodeID = nextFreeNode( );

   atomsOctree[ nodeID ].atomsStartID = atomsStartID;
   atomsOctree[ nodeID ].atomsEndID = atomsEndID;

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
           int j = constructAtomsOctree( atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atoms, atomsOctree );
           atomsOctree[ nodeID ].cPtr[ i ] = j;
          }
        else atomsOctree[ nodeID ].cPtr[ i ] = -1;
     }

   return nodeID;
}


void resContFilter::fillMovingAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel )
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


bool resContFilter::buildStaticAtomsOctree( void )
{
   if ( printStatus ) printf( "\nbuilding static atoms octree... " );

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
      printError( (char *)"Unable to %s static atoms octree - memory allocation failed!", ( staticAtomsOctreeBuilt ) ? "rebuild" : "build" );
      if ( !staticAtomsOctreeBuilt ) exit( 1 );
      return false;
     }

   freeMem( staticAtomsOctree );
   staticAtomsOctree = atomsOctreeT;

   initFreeNodeServer( nStaticAtoms );

   staticAtomsOctreeRoot = constructAtomsOctree( 0, nStaticAtoms - 1, staticAtoms, staticAtomsOctree );

   staticAtomsOctreeBuilt = true;

   double endT = getTime( );

   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );

   return true;
}



bool resContFilter::buildMovingAtomsOctree( void )
{
   if ( printStatus ) printf( "\nbuilding moving atoms octree... " );

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
      printError( (char *)"Unable to %s moving atoms octree - memory allocation failed!", ( movingAtomsOctreeBuilt ) ? "rebuild" : "build" );
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

   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );

   return true;
}


bool resContFilter::buildOctrees( void )
{
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed )
     || ( distCutoff != distCutoffUsed ) || ( !staticAtomsOctreeBuilt ) ) buildStaticAtomsOctree( );

   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed )
     || ( distCutoff != distCutoffUsed ) || ( !movingAtomsOctreeBuilt ) ) buildMovingAtomsOctree( );

   minRadiusUsed = minRadius;
   maxLeafSizeUsed = maxLeafSize;
   distCutoffUsed = distCutoff;
}


void resContFilter::computeResResInteractions( Matrix transMat, int nodeS, int nodeM, double *interactionValuePos, double *interactionValueNeg )
{
   double sumRad = staticAtomsOctree[ nodeS ].cr + movingAtomsOctree[ nodeM ].cr;

   Vector oldPos( movingAtomsOctree[ nodeM ].cx, movingAtomsOctree[ nodeM ].cy, movingAtomsOctree[ nodeM ].cz, 1.0 );
   Vector newPos = transMat * oldPos;
   double cxM = newPos[ 0 ], cyM = newPos[ 1 ], czM = newPos[ 2 ];

   double dx = staticAtomsOctree[ nodeS ].cx - cxM,
          dy = staticAtomsOctree[ nodeS ].cy - cyM,
          dz = staticAtomsOctree[ nodeS ].cz - czM;
   double d2 = dx * dx + dy * dy + dz * dz;

   bool farEnough = false;

   if ( d2 > ( sumRad + distCutoff ) * ( sumRad + distCutoff ) ) farEnough = true;

   *interactionValuePos = *interactionValueNeg = 0;

   if ( farEnough ) return;
   else
      {
       if ( staticAtomsOctree[ nodeS ].leaf && movingAtomsOctree[ nodeM ].leaf )
         {
           for ( int i = movingAtomsOctree[ nodeM ].atomsStartID; i <= movingAtomsOctree[ nodeM ].atomsEndID; i++ )
             {
              Vector oldPos( movingAtoms[ i ].x, movingAtoms[ i ].y, movingAtoms[ i ].z, 1.0 );
              Vector newPos = transMat * oldPos;
              double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];
              int rIDM = movingAtoms[ i ].rID;

              for ( int j = staticAtomsOctree[ nodeS ].atomsStartID; j <= staticAtomsOctree[ nodeS ].atomsEndID; j++ )
                {
                 dx = staticAtoms[ j ].x - xM;
                 dy = staticAtoms[ j ].y - yM;
                 dz = staticAtoms[ j ].z - zM;

                 int rIDS = staticAtoms[ j ].rID;

                 d2 = dx * dx + dy * dy + dz * dz;

                 if ( d2 <= distCutoffSQ )
                   {
                     if ( intVal[ rIDM ][ rIDS ] > 0.0 ) ( *interactionValuePos ) += intVal[ rIDM ][ rIDS ];
                     else if ( intVal[ rIDM ][ rIDS ] < -3.0 ) ( *interactionValueNeg ) += ( - intVal[ rIDM ][ rIDS ] );
                   }
                }
             }
         }
       else if ( !staticAtomsOctree[ nodeS ].leaf && !movingAtomsOctree[ nodeM ].leaf )
              {
                for ( int i = 0; i < 8; i++ )
                  if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 )
                         {
                           double intValP, intValN;

                           computeResResInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], movingAtomsOctree[ nodeM ].cPtr[ j ], &intValP, &intValN );

                           ( *interactionValuePos ) += intValP;
                           ( *interactionValueNeg ) += intValN;
                         }
              }
            else if ( !staticAtomsOctree[ nodeS ].leaf )
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                         {
                           double intValP, intValN;

                           computeResResInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], nodeM, &intValP, &intValN );

                           ( *interactionValuePos ) += intValP;
                           ( *interactionValueNeg ) += intValN;
                         }
                   }
                 else
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 )
                         {
                           double intValP, intValN;

                           computeResResInteractions( transMat, nodeS, movingAtomsOctree[ nodeM ].cPtr[ j ], &intValP, &intValN );

                           ( *interactionValuePos ) += intValP;
                           ( *interactionValueNeg ) += intValN;
                         }
                   }
      }
}



bool resContFilter::computeInteractions( Matrix transMat, double *interactionValuePos, double *interactionValueNeg )
{
   if ( !staticAtomsOctreeBuilt || !movingAtomsOctreeBuilt ) return false;

   if ( printStatus ) printf( "\ncomputing residue-residue interactions... " );

   double startT = getTime( );

   if ( numThreads == 1 ) computeResResInteractions( transMat, staticAtomsOctreeRoot, movingAtomsOctreeRoot, interactionValuePos, interactionValueNeg );
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
           pthread_create( &p[ i ], NULL, computeResResInteractionsThread, ( void * ) &threadResults[ i ] );
         }

       for ( int i = 0; i < numThreads; i++ )
          pthread_join( p[ i ], NULL );

       *interactionValuePos = *interactionValueNeg = 0;       

       for ( int i = 0; i < numThreads; i++ )
         {
           ( *interactionValuePos ) += threadResults[ i ].interactionValuePos;
           ( *interactionValueNeg ) += threadResults[ i ].interactionValueNeg;           
         }  
      }

   double endT = getTime( );

   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );

   return true;
}



bool resContFilter::computeInteractions( double *interactionValuePos, double *interactionValueNeg )
{
   return computeInteractions( transMatrix, interactionValuePos, interactionValueNeg );
}


bool resContFilter::computeInteractionsNaively( Matrix transMat, double *interactionValuePos, double *interactionValueNeg )
{
   if ( printStatus ) printf( "\ncomputing residue-residue interactions naively... " );

   double startT = getTime( );

   *interactionValuePos = *interactionValueNeg = 0.0;

   int numClose = 0;
   for ( int i = 0; i < nMovingAtoms; i++ )
     {
      Vector oldPos( movingAtoms[ i ].x, movingAtoms[ i ].y, movingAtoms[ i ].z, 1.0 );
      Vector newPos = transMat * oldPos;

      double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];

      int rIDM = movingAtoms[ i ].rID;

      for ( int j = 0; j < nStaticAtoms; j++ )
        {
         double dx = staticAtoms[ j ].x - xM;
         double dy = staticAtoms[ j ].y - yM;
         double dz = staticAtoms[ j ].z - zM;

         int rIDS = staticAtoms[ j ].rID;

         double d2 = dx * dx + dy * dy + dz * dz;

         if ( d2 <= distCutoffSQ ) 
           {
             if ( intVal[ rIDM ][ rIDS ] > 0 ) ( *interactionValuePos ) += intVal[ rIDM ][ rIDS ];
             else if ( intVal[ rIDM ][ rIDS ] < -3.0 ) ( *interactionValueNeg ) += ( - intVal[ rIDM ][ rIDS ] );
           }
        }
     }

   double endT = getTime( );

   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );

   return true;
}


bool resContFilter::computeInteractionsNaively( double *interactionValuePos, double *interactionValueNeg )
{
   return computeInteractionsNaively( transMatrix, interactionValuePos, interactionValueNeg );
}
