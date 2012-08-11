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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef FAST_BORN_RADIUS_H

#define FAST_BORN_RADIUS_H

#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <pthread.h>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include <time.h>

#ifdef _WIN32
   #include <sys/types.h>
   #include <sys/timeb.h>
#else
   #include <sys/time.h>
#endif

#include "../utils/utils.h"

#ifdef COMBINE_INTEGRALS
   #undef COMBINE_INTEGRALS
#endif

#ifdef BRR4
   #define COMBINE_INTEGRALS( a, b ) ( ( 1 - INV_SQRT_TWO ) * ( a ) / ( 4 * M_PI ) + pow( ( b ) / ( 16 * M_PI ), 0.5 ) )
#else   
//   #define COMBINE_INTEGRALS( a, b ) ( pow( ( 3 * ( a ) ) / ( 20 * M_PI ), 1.0 / 3.0 ) )
   #define COMBINE_INTEGRALS( a, b ) ( pow( ( a ) / ( 4 * M_PI ), 1.0 / 3.0 ) )
#endif   


namespace fastGB {

class fastBornRadius
{
 private:

   enum atomType { N = 0, C, O, H, S, P };
   double vdwRad[ 6 /*sizeof( atomType )*/ ]; 

   typedef struct
    {
     int id;
     double x, y, z;
     double q, r;
     double invR, invCR, R, naiveR;
    } ATOM;

   typedef struct
    {
     double x, y, z;
     double nx, ny, nz;
     double w;     
    } QPOINT;
    
   typedef struct
    {
     double cx, cy, cz;
     double wnx, wny, wnz;
     double cr;
     bool leaf;
     int cPtr[ 8 ];
     int qPtsStartID, qPtsEndID;
    } QPOINTS_OCTREE_NODE;    

   typedef struct
    {
     double cx, cy, cz;
     double cr;
     double invR, invCR;
     bool leaf;
     int cPtr[ 8 ];
     int atomsStartID, atomsEndID;
    } ATOMS_OCTREE_NODE;    

   QPOINT *qPoints;      
   int nQPoints;
   
   int numQPointsOctreeNodes;
   bool qPointsOctreeBuilt;
   QPOINTS_OCTREE_NODE *qPointsOctree;
   int qPointsOctreeRoot;
   
   ATOM *atoms;
   int nAtoms;      

   int numAtomsOctreeNodes;
   bool atomsOctreeBuilt;
   ATOMS_OCTREE_NODE *atomsOctree;
   int atomsOctreeRoot;
   
   bool priorComputationCleared;
   
   bool computedBornRadiiNaively;

   double minRadius, minRadiusUsed;
   int maxLeafSize, maxLeafSizeUsed;
   
   double epsilon, epsilonUsed;
   
   double maxBornRadius, maxBornRadiusInv;
      
   pthread_mutex_t nodesLock;
   int curNode, maxNode;   

   pthread_mutex_t subtreeRootsLock;
   int curSubtreeRoot, maxSubtreeRoot;   
   int *atomsSubtreeRoots;
   int atomsSubtreeRootsSize;
   
   double fastBRTime, naiveBRTime;
   
   int numThreads;
   
   bool printStatus;

   void freeMemory( void );
   void setDefaults( void );
   bool allocateAtomsSubtreeRootsArray( int nThreads );   
   bool readQPoints( char *qPtsFile );   
   bool copyQPointsFromArray( int numQPoints, double *qPts );
   bool readAtomsFromPQR( char *atomsFile );
   bool readAtoms( char *atomsFile ); 
   bool copyAtomsFromArray( int numAtoms, double *atms );   
   void initFreeNodeServer( int numNodes );   
   int nextFreeNode( void );
   void initSubtreeRootServer( );   
   int nextSubtreeRoot( void );   
   void countQPointsOctreeNodesAndSortQPoints( QPOINT *qPts, int qPtsStartID, int qPtsEndID, QPOINT *qPtsT, int *numNodes );   
   int constructQPointsOctree( int qPtsStartID, int qPtsEndID, QPOINT *qPts );
   bool buildQPointsOctree( void );   
   void countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID, ATOM *atomsT, int *numNodes );
   int constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms );      
   bool buildAtomsOctree( void );            
   void fillAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel );   
   void pushIntegralsAndFillAtomsSubtreeRootsArray( int nodeID, double intVal, double intCVal, int maxNodesInLevel );
   void transformAtomsOctree( int nodeID, double *transMatrix );
   void cleanupAtomsOctree( int nodeID );
   void cleanupPriorComputation( void );
   void approximateIntegrals( int nodeA, int nodeQ );   
   void pushIntegralsToAtoms( int nodeID, double intVal, double intCVal );      

   static void *approximateIntegralsThread( void *v )
    {
     fastGB::fastBornRadius *fBR = ( fastGB::fastBornRadius * ) v;
     
     while ( 1 )
       {
        int nextRoot = fBR->nextSubtreeRoot( );     
        if ( nextRoot < 0 ) break;      
        fBR->approximateIntegrals( nextRoot, fBR->qPointsOctreeRoot );
       }
    }

   static void *pushIntegralsToAtomsThread( void *v )
    {
     fastGB::fastBornRadius *fBR = ( fastGB::fastBornRadius * ) v;
     
     while ( 1 )
       {
        int nextRoot = fBR->nextSubtreeRoot( );     
        if ( nextRoot < 0 ) break;      
        fBR->pushIntegralsToAtoms( nextRoot, 0, 0 );
       }
    }
    
 public:

   fastBornRadius( char *qPtsFile, char *atmsPQRFile, bool printStat );   
   fastBornRadius( char *qPtsFile, int numAtoms, double *atms, bool printStat );      
   fastBornRadius( int numQPoints, double *qPts, int numAtoms, double *atms, bool printStat );   
   ~fastBornRadius( );
   
   bool setMinRadius( double minRad );
   bool setMaxLeafSize( int maxLfSize );
   bool setEpsilon( double eps );
   bool setMaxBornRadius( double maxBornRad );   
   bool setNumThreads( int nThreads );

   void setPrintStatus( bool printStat );
      
   void printCurrentSettings( void );   
   
   void transformAtoms( double *transMatrix );   
   
   bool buildOctrees( void );      

   bool computeBornRadii( void );
   bool computeBornRadiiNaively( void );
 
   bool getQPoints( int *numQPoints, double **qPts );
   bool getAtomsPQR( int *numAtoms, double **atomsPQR );
   bool getAtomsPQRR( int *numAtoms, double **atomsPQRR );
   bool getAtomsPQRR( int *numAtoms, double **atomsPQRR, bool getNaiveR );   
   
   double *getBornRadii( void );
   double *getIntegrals( void );
   
   double getDispersionEnergy( void );

   double *getNaiveBornRadii( void );
   
   bool writeBornRadiiToFile( char *outFile );
   bool writeIntegralsToFile( char *outFile ); 

   bool writePQRRFile( char *outFile );
   
   bool writeFastAndNaiveBornRadiiToFile( char *outFile );
};

};

#endif
