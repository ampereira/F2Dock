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


#ifndef FAST_G_POL_H

#define FAST_G_POL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <pthread.h>

#ifdef _WIN32
   #include <sys/types.h>
   #include <sys/timeb.h>
#else
   #include <sys/time.h>
#endif

#ifdef USE_SSE
   #include "SSEApproxMath.h"
#endif

#include "../utils/utils.h"


namespace fastGB {

class fastGpol
{
 private:

   typedef struct
     {
       double cx, cy, cz;
       double cr;
       bool leaf;
       int cPtr[ 8 ];
       int atomStartID, atomEndID;
     } OCTREE_NODE;
        
   typedef struct
     {
       double x, y, z;
       double r, q;
       double R;
       int groupR;
     } ATOM;
        
   typedef union 
     {
       unsigned long ix;
       float fx;
     } UL_F_union;
        
   typedef struct
     {
       int nodeI, nodeJ;
     } JOB;
    
   typedef struct
     {
       int threadID; // unique integer ID of the thread
        
       OCTREE_NODE *octree;
       double *qSum;
       double *approxR;
       int numGroupR;
       ATOM *atom;
       double eps;
       
       double Gpol;	     	    
       
       fastGpol *thisPtr;
     } PARAMS;

#ifdef USE_SSE
   v4sf ZERO, FOUR; 
#endif

   std::vector< JOB > jobs;
   pthread_mutex_t jobLock;
   
   int numThreads;
   double epsilon;
   double minRad;
   int maxLeafAtoms;   
   bool useApproxMath;
   
   int numAtoms;
   ATOM *atoms;
   OCTREE_NODE *octree;
   double *approxR;
   double *qSum;   

   int maxGroupR;
   
   double e_in, e_out;   
   
   bool printStatus;

   void setDefaults( void );   
   void freeMemory( void );
   bool readFromPQRRFile( char *fname );
   bool readFromPQRRArray( int nAtoms, double *pqrR );
   bool initDataStructures( void );

   int nextJob( int &nodeI, int &nodeJ );
   void initJobServer( std::vector< JOB > &leafJobs, std::vector< JOB > &nonLeafJobs );
   void generateJobs( OCTREE_NODE *octree, int nodeI, int nodeJ, double eps, std::vector< JOB > &leafJobs, std::vector< JOB > &nonLeafJobs );

   inline float invSqrt( float x );
   inline float fastExp( float x );

#ifdef USE_SSE
   inline float vectorGB( v4sf xi, v4sf yi, v4sf zi, v4sf Ri, v4sf qi,
                          v4sf xj, v4sf yj, v4sf zj, v4sf Rj, v4sf qj );
   inline float vectorGB( float rij2f, v4sf Ri, v4sf qi, v4sf Rj, v4sf qj );
#endif

   void countOctreeNodesAndSortAtoms( ATOM *atom, int atomStartID, int atomEndID, 
                                      double minRad2, int maxLeafAtoms,
                                      ATOM *atomT, int *numNodes );
   void constructOctree( OCTREE_NODE *octree, int nodeID, double *qSum, 
                         int numGroupR, double minRad, int maxLeafAtoms, 
                         ATOM *atom, int *freeID );
   void recursiveFastGpol( OCTREE_NODE *octree, int nodeI, int nodeJ, 
                           double *qSum, double *approxR, int numGroupR, ATOM *atom, 
                           double eps, double *Gpol );
   void computeThreadedFastGpol( int threadID, OCTREE_NODE *octree,
                                 double *qSum, double *approxR, int numGroupR, ATOM *atom, 
                                 double eps, double *Gpol );
                                 
   static void *startThreadedFastGpol( void *v )
      {
        fastGB::fastGpol::PARAMS *pr = ( fastGB::fastGpol::PARAMS * ) v;

        pr->thisPtr->computeThreadedFastGpol( pr->threadID, pr->octree, pr->qSum, pr->approxR, pr->numGroupR, pr->atom, pr->eps, &( pr->Gpol ) );
      }
                                 
   
 public:

   fastGpol( char *pqrRFile );
   fastGpol( char *pqrRFile, bool printStat );   
   fastGpol( int nAtoms, double *pqrR );
   fastGpol( int nAtoms, double *pqrR, bool printStat );   
   ~fastGpol( );
    
   bool setMinRadius( double minRadius ); 
   bool setMaxLeafSize( int maxLfSize ); 
   bool setEpsilon( double eps );
   bool setNumThreads( int nThreads ); 
   bool setIonDielectric( double e_ion );
   bool setSolventDielectric( double e_sol );
   void useApproxMathFunctions( bool approxMath );
   
   void setPrintStatus( bool printStat );   
   
   void printCurrentSettings( void );
      
   void computeQuadGpol( double *Gpol );   
   void computeFastGpol( double *Gpol );
};

};

#endif

