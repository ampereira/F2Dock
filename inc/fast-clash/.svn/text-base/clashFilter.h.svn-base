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

#ifndef CLASH_FILTER_H

#define CLASH_FILTER_H

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
#include <stdarg.h>
#include "../math/Matrix.h"
#include "../math/Vector.h"

#ifdef _WIN32
   #include <sys/types.h>
   #include <sys/timeb.h>
#else
   #include <sys/time.h>
#endif

#ifdef freeMem
   #undef freeMem
#endif
#define freeMem( ptr ) { if ( ptr != NULL ) free( ptr ); }

#ifdef zeroIfLess   
   #undef zeroIfLess
#endif
#define zeroIfLess( a, b ) ( ( ( a ) < ( b ) ) ? 0 : 1 )


using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;


class clashFilter
{
 private:

   typedef struct
    {
     int id;
     double x, y, z;
     double q, r;
    } ATOM;

   typedef struct
    {
     double cx, cy, cz;
     double cq, cr;
     bool leaf;
     int cPtr[ 8 ];
     int atomsStartID, atomsEndID;
    } ATOMS_OCTREE_NODE;    

   typedef struct
    {
     clashFilter *cF;
     Matrix *transMat;
     int nClashes, nSevereClashes;
     double interactionValue;
    } THREAD_RESULT;

   ATOM *staticAtoms;
   int nStaticAtoms;      

   int numStaticAtomsOctreeNodes;
   bool staticAtomsOctreeBuilt;
   ATOMS_OCTREE_NODE *staticAtomsOctree;
   int staticAtomsOctreeRoot;

   ATOM *movingAtoms;
   int nMovingAtoms;      

   int numMovingAtomsOctreeNodes;
   bool movingAtomsOctreeBuilt;
   ATOMS_OCTREE_NODE *movingAtomsOctree;
   int movingAtomsOctreeRoot;
   
   bool priorComputationCleared;
   
   bool computedInteractionNaively;

   double minRadius, minRadiusUsed;
   int maxLeafSize, maxLeafSizeUsed;
   
   double epsilon;
   
   pthread_mutex_t nodesLock;
   int curNode, maxNode;   

   pthread_mutex_t subtreeRootsLock;
   int curSubtreeRoot, maxSubtreeRoot;   
   int *movingAtomsSubtreeRoots;
   int movingAtomsSubtreeRootsSize;
   
   Matrix transMatrix;   
   double clashFrac, severeClashFrac, fuzzyFrac;     
      
   int numThreads;   
   
   bool printStatus;

   void printError( char *format, ... );
   double getTime( void );   
   void freeMemory( void );
   void setDefaults( void );
   bool allocateMovingAtomsSubtreeRootsArray( int nThreads );   
   void initFreeNodeServer( int numNodes );   
   int nextFreeNode( void );
   void initSubtreeRootServer( );   
   int nextSubtreeRoot( void );   
   bool copyAtomsFromArray( int numAtomsSrc, double *atmsSrc, int *numAtomsDest, ATOM **atmsDest );   
   void countAtomsOctreeNodesAndSortAtoms( ATOM *sAtoms, int sAtomsStartID, int sAtomsEndID, ATOM *sAtomsT, int *numNodes );      
   int constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms, ATOMS_OCTREE_NODE *atomsOctree );
   bool buildStaticAtomsOctree( void );   
   bool buildMovingAtomsOctree( void );            
   bool buildOctrees( void );         
   void fillMovingAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel );   
   void approximateInteractions( Matrix transMat, int nodeS, int nodeM, int *nClashes, int *nSevereClashes, double *interactionValue );

   static void *approximateInteractionsThread( void *v )
    {
     THREAD_RESULT *tR = ( THREAD_RESULT * ) v;     
     
     tR->nClashes = tR->nSevereClashes = 0; 
     tR->interactionValue = 0;
     
     while ( 1 )
       {
        int nextRoot = tR->cF->nextSubtreeRoot( );     
        if ( nextRoot < 0 ) break;      
        
        int nClashes, nSevereClashes;
        double interactionValue;
        
        tR->cF->approximateInteractions( *( tR->transMat ), tR->cF->staticAtomsOctreeRoot, nextRoot, &nClashes, &nSevereClashes, &interactionValue );
        
        tR->nClashes += nClashes;
        tR->nSevereClashes += nSevereClashes;
        tR->interactionValue += interactionValue;
       }
    }

    
 public:

   clashFilter( int numStaticAtoms, double *staticAtoms, int numMovingAtoms, double *movingAtoms, bool printStat ); 
   clashFilter( int numStaticAtoms, double *staticAtoms, int numMovingAtoms, double *movingAtoms, bool printStat , double cF); 	  
   clashFilter( int numStaticAtoms, double *staticAtoms, int numMovingAtoms, double *movingAtoms );      
   ~clashFilter( );
   
   bool setMinRadius( double minRad );
   bool setMaxLeafSize( int maxLfSize );
   bool setEpsilon( double eps );
   bool setNumThreads( int nThreads );

   bool setProximityFactors( double clashFactor, double severeClashFactor, double fuzzyFactor );
   bool setTransformationMatrix( Matrix transMat );

   void setPrintStatus( bool printStat );
      
   void printCurrentSettings( void );   
   
   bool computeInteractions( int *nClashes, int *nSevereClashes, double *interactionValue );
   bool computeInteractions( Matrix transMat, int *nClashes, int *nSevereClashes, double *interactionValue );
   
   bool computeInteractionsNaively( int *nClashes, int *nSevereClashes, double *interactionValue );
   bool computeInteractionsNaively( Matrix transMat, int *nClashes, int *nSevereClashes, double *interactionValue );   
};

#endif
