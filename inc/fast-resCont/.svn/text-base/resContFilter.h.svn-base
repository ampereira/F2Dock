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


#ifndef RESCONT_FILTER_H

#define RESCONT_FILTER_H

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

#define NUM_RESIDUE_TYPES 20


using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;


// The following table is from ( Table III, page 94 ): 
//     Fabian Glaser, David M. Steinberg, Ilya A. Vakser, and Nir Ben-Tal,
//     "Residue Frequencies and Pairing Preferences at Protein–Protein Interfaces",
//     PROTEINS: Structure, Function, and Genetics 43:89–102 (2001)

// residue order in intValDefault: ILE, VAL, LEU, PHE, CYS, MET, ALA, GLY, THR, SER, TRP, TYR, PRO, HIS, GLU, GLN, ASP, ASN, LYS, ARG
const double intValDefault[ NUM_RESIDUE_TYPES ][ NUM_RESIDUE_TYPES ] 
 = { { -0.67,  1.07,  0.03,  0.40, -1.21,  0.68,  0.95,  0.33, -0.27, -1.26,  0.49,  0.49,  0.09, -0.97, -0.80, -0.60, -0.92, -1.87, -1.70, -1.53 }, 
     {  1.07,  0.62,  0.36,  0.48,  0.64,  0.53,  1.40, -0.13,  0.24, -0.12, -2.11, -0.45,  0.45, -0.41, -0.06, -0.26, -0.57, -1.37,  0.23, -0.43 }, 
     {  0.03,  0.36, -0.53, -0.07, -0.05,  0.76,  0.88, -0.82, -1.25, -0.85,  0.02, -0.93, -0.67,  0.53, -0.89, -0.74, -1.82, -1.14, -1.78, -0.34 }, 
     {  0.40,  0.48, -0.07,  0.04,  0.34,  0.35,  0.74, -0.68, -0.34, -0.89, -0.29,  0.33,  0.70, -1.25, -1.51, -0.32, -2.61, -0.71, -1.73, -1.21 }, 
     { -1.21,  0.64, -0.05,  0.34,  6.26, -1.13,  1.16,  0.90, -0.70,  1.81, -2.02, -1.06,  1.16,  1.36,  0.10, -1.28, -1.39, -2.28, -1.29, -0.93 }, 
     {  0.68,  0.53,  0.76,  0.35, -1.13,  1.45,  0.41,  0.47, -1.23, -0.65, -0.86, -0.32,  0.21,  0.30, -0.13, -0.02, -2.86, -1.15, -1.01, -1.71 }, 
     {  0.95,  1.40,  0.88,  0.74,  1.16,  0.41,  0.26,  0.46,  0.57,  0.80,  0.29,  0.02,  0.72,  0.92,  0.38,  0.19,  0.59,  0.91, -0.13, -0.76 }, 
     {  0.33, -0.13, -0.82, -0.68,  0.90,  0.47,  0.46, -0.72,  1.01,  0.33, -0.21,  0.25,  0.44,  0.85, -0.77,  0.62,  0.82,  0.13,  0.52,  0.38 }, 
     { -0.27,  0.24, -1.25, -0.34, -0.70, -1.23,  0.57,  1.01, -0.79,  0.89,  0.62, -0.73,  0.72, -0.39,  0.12, -1.13,  1.91,  0.32, -0.02, -0.32 }, 
     { -1.26, -0.12, -0.85, -0.89,  1.81, -0.65,  0.80,  0.33,  0.89, -0.06, -0.58, -0.52,  0.46, -1.25,  0.90,  0.10,  2.02,  0.62,  0.10, -0.21 }, 
     {  0.49, -2.11,  0.02, -0.29, -2.02, -0.86,  0.29, -0.21,  0.62, -0.58, -1.09, -0.12,  3.51,  0.92, -3.99, -4.02, -1.79, -1.10, -0.36,  2.05 }, 
     {  0.49, -0.45, -0.93,  0.33, -1.06, -0.32,  0.02,  0.25, -0.73, -0.52, -0.12,  0.25,  0.49,  1.14, -0.02, -2.71, -2.02, -0.35, -0.23, -0.61 }, 
     {  0.09,  0.45, -0.67,  0.70,  1.16,  0.21,  0.72,  0.44,  0.72,  0.46,  3.51,  0.49, -1.18, -0.07,  0.55,  0.69, -0.37,  1.03,  0.21,  0.05 }, 
     { -0.97, -0.41,  0.53, -1.25,  1.36,  0.30,  0.92,  0.85, -0.39, -1.25,  0.92,  1.14, -0.07,  1.24, -1.49,  0.02,  2.20, -0.85, -2.00, -0.22 }, 
     { -0.80, -0.06, -0.89, -1.51,  0.10, -0.13,  0.38, -0.77,  0.12,  0.90, -3.99, -0.02,  0.55, -1.49, -1.80, -1.69, -2.59, -0.22,  0.95,  0.97 }, 
     { -0.60, -0.26, -0.74, -0.32, -1.28, -0.02,  0.19,  0.62, -1.13,  0.10, -4.02, -2.71,  0.69,  0.02, -1.69, -1.00,  0.41,  0.36, -1.07, -0.47 }, 
     { -0.92, -0.57, -1.82, -2.60, -1.39, -2.86,  0.59,  0.82,  1.91,  2.02, -1.79, -2.02, -0.37,  2.20, -2.59,  0.41, -1.74,  1.74,  0.31,  0.95 }, 
     { -1.87, -1.37, -1.14, -0.71, -2.28, -1.15,  0.91,  0.13,  0.32,  0.62, -1.10, -0.35,  1.03, -0.85, -0.22,  0.36,  1.74,  0.58, -0.66, -0.38 }, 
     { -1.70,  0.23, -1.78, -1.73, -1.29, -1.01, -0.13,  0.52, -0.02,  0.10, -0.36, -0.23,  0.21, -2.00,  0.95, -1.07,  0.31, -0.66, -2.06, -3.41 }, 
     { -1.53, -0.43, -0.34, -1.21, -0.93, -1.71, -0.76,  0.38, -0.32, -0.21,  2.05, -0.61,  0.05, -0.22,  0.97, -0.47,  0.95, -0.38, -3.41, -3.23 } };


class resContFilter
{
 private:

   double intVal[ NUM_RESIDUE_TYPES ][ NUM_RESIDUE_TYPES ];

   typedef struct
    {
     int id;
     double x, y, z;
     int rID, rNum;
    } ATOM;

   typedef struct
    {
     double cx, cy, cz;
     double cr;
     bool leaf;
     int cPtr[ 8 ];
     int atomsStartID, atomsEndID;
    } ATOMS_OCTREE_NODE;    

   typedef struct
    {
     resContFilter *cF;
     Matrix *transMat;
     double interactionValuePos, interactionValueNeg;
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
   
   bool interactionMatrixRead;
   
   bool priorComputationCleared;
   
   bool computedInteractionNaively;

   double minRadius, distCutoff, distCutoffSQ;
   int maxLeafSize, minRadiusUsed, maxLeafSizeUsed, distCutoffUsed;
   
   pthread_mutex_t nodesLock;
   int curNode, maxNode;   

   pthread_mutex_t subtreeRootsLock;
   int curSubtreeRoot, maxSubtreeRoot;   
   int *movingAtomsSubtreeRoots;
   int movingAtomsSubtreeRootsSize;
   
   Matrix transMatrix;   
      
   int numThreads;   
   
   bool printStatus;

   void printError( char *format, ... );
   double getTime( void );   
   void freeMemory( void );
   void setDefaults( void );
   bool readInteractionMatrix( char *intFile );   
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
   void computeResResInteractions( Matrix transMat, int nodeS, int nodeM, double *interactionValuePos, double *interactionValueNeg );

   static void *computeResResInteractionsThread( void *v )
    {
     THREAD_RESULT *tR = ( THREAD_RESULT * ) v;     
     
     tR->interactionValuePos = tR->interactionValueNeg = 0;
     
     while ( 1 )
       {
        int nextRoot = tR->cF->nextSubtreeRoot( );     
        if ( nextRoot < 0 ) break;      
        
        double interactionValuePos, interactionValueNeg;
        
        tR->cF->computeResResInteractions( *( tR->transMat ), tR->cF->staticAtomsOctreeRoot, nextRoot, &interactionValuePos, &interactionValueNeg );
        
        tR->interactionValuePos += interactionValuePos;
        tR->interactionValueNeg += interactionValueNeg;        
       }
    }

    
 public:

   resContFilter( int numStaticAtoms, double *staticAtoms, int numMovingAtoms, double *movingAtoms, char *intFile, bool printStat );   
   ~resContFilter( );
   
   bool setMinRadius( double minRad );
   bool setMaxLeafSize( int maxLfSize );
   bool setDistanceCutoff( double distanceCutoff );
   bool setNumThreads( int nThreads );

   bool setTransformationMatrix( Matrix transMat );

   void setPrintStatus( bool printStat );
      
   void printCurrentSettings( void );   
   
   bool computeInteractions( double *interactionValuePos, double *interactionValueNeg );
   bool computeInteractions( Matrix transMat, double *interactionValuePos, double *interactionValueNeg );
   
   bool computeInteractionsNaively( double *interactionValuePos, double *interactionValueNeg );
   bool computeInteractionsNaively( Matrix transMat, double *interactionValuePos, double *interactionValueNeg );   
};

#endif
