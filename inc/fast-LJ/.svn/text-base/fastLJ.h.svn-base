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


#ifndef FAST_LJ_H

#define FAST_LJ_H

#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cstring>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include <time.h>
#include <stdarg.h>

#ifdef _WIN32
   #include <sys/types.h>
   #include <sys/timeb.h>
#else
   #include <sys/time.h>
#endif

#include "../utils/utils.h"

class fastLJ
{
 private:

   enum ATOM_TYPES { C, H, N, O, P, S, numAtomType };

   typedef struct
    {
     int id;
     double x, y, z;
    } ATOM;

   typedef struct
    {
     double cx, cy, cz;
     double cr;     
     bool leaf;
     int cPtr[ 8 ];
     int numAtoms[ numAtomType ], atomsStartID[ numAtomType ], atomsEndID[ numAtomType ];
    } ATOMS_OCTREE_NODE;    

   ATOM *staticAtoms[ numAtomType ];
   int nStaticAtoms[ numAtomType ];      

   int numStaticAtomsOctreeNodes;
   bool staticAtomsOctreeBuilt;
   ATOMS_OCTREE_NODE *staticAtomsOctree;
   int staticAtomsOctreeRoot;

   ATOM *movingAtoms[ numAtomType ], *movingAtomsOrig[ numAtomType ];
   int nMovingAtoms[ numAtomType ];      

   int numMovingAtomsOctreeNodes;
   bool movingAtomsOctreeBuilt;
   ATOMS_OCTREE_NODE *movingAtomsOctree;
   int movingAtomsOctreeRoot;

   double r_eqm_XY[ numAtomType ][ numAtomType ];
   double A[ numAtomType ][ numAtomType ];
   double B[ numAtomType ][ numAtomType ];
   double vdWEqmRadScale;
   
   double minRadius, minRadiusUsed;
   int maxLeafSize, maxLeafSizeUsed;
   double minInterAtomDist, minInterAtomDistUsed;   
   double minD2;
   
   double epsilon;
   
   int curNode, maxNode;   
   
   double transMatrix[ 16 ];   
      
   int numThreads;   
   
   bool printStatus;
   
#ifdef USE_SSE
   v4sf MIND2; 
#endif

   bool useSSEFunctions;
   
   void freeMemory( void );
   void computeConstants( void );   
   void setDefaults( void );
   void transformMovingAtoms( int threadID, double *transMat );   
   void initFreeNodeServer( int numNodes );   
   int nextFreeNode( void );
   bool getParamsFromFile( char *paramFile, char **staticAtomsFile, char **movingAtomsFile );   
   bool readAtomsFromPDB( char *atomsFile, int *numAtoms, ATOM **atms );   
   bool readAtomsFromPQR( char *atomsFile, int *numAtoms, ATOM **atms );
   bool readAtomsFromFile( char *atomsFile, int *numAtoms, ATOM **atms );   
   bool copyAtomsFromArray( int numAtomsSrc, double *atmsSrcXYZ, char *atmsSrcType, int *numAtomsDest, ATOM **atmsDest );   
   void countAtomsOctreeNodesAndSortAtoms( ATOM **sAtoms, int *sAtomsStartID, int *sAtomsEndID, ATOM **sAtomsT, int *numNodes );      
   int constructAtomsOctree( int *atomsStartID, int *atomsEndID, ATOM **atoms, ATOMS_OCTREE_NODE *atomsOctree );
   bool buildStaticAtomsOctree( void );   
   bool buildMovingAtomsOctree( void );            
   bool buildOctrees( void );         
   
#ifdef USE_SSE
   inline float vectorLJ( v4sf xi, v4sf yi, v4sf zi, v4sf xj, v4sf yj, v4sf zj, v4sf Aij, v4sf Bij );
#endif
   
   void approximatePotential( int threadID, double *transMat, int nodeS, int nodeM, double *LJPot );
    
 public:

   fastLJ( int numStaticAtoms, double *stAtomsXYZ, char *stAtomsType, int numMovingAtoms, double *mvAtomsXYZ, char *mvAtomsType, int numThreads, bool printStat );
   fastLJ( char *staticAtomsFile, char *movingAtomsFile, int numThreads, bool printStat );   
   fastLJ( char *paramFile );            
   ~fastLJ( );
   
   bool setMinRadius( double minRad );
   bool setMaxLeafSize( int maxLfSize );
   bool setMinInterAtomDist( double minAtomDist );   
   bool setEpsilon( double eps );

   void useSSE( bool useSSEFuncs );

   void setTransformationMatrix( double *transMat );

   void setPrintStatus( bool printStat );
      
   void printCurrentSettings( void );   
   
   bool computePotential( int threadID, double *LJPot );
   bool computePotential( int threadID, double *transMat, double *LJPot );
   
   void computePotentialNaively( int threadID, double *LJPot );
   void computePotentialNaively( int threadID, double *transMat, double *LJPot );   
};

#endif
