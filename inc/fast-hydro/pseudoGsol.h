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


#ifndef PSEUDO_GSOL_H

#define PSEUDO_GSOL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#include "utils/utils.h"
#include "f2dock/ElementInformation.h"
#include "PG-range/PG.h"

#ifdef LID
   #undef LID
#endif

#define LID( r, c ) ( ( ( r ) * 4 ) + ( c ) )

class pseudoGsol
{
 private:

    typedef struct
      {
        char *staticMoleculePQR;
        char *movingMoleculePQR;        
        char *staticMoleculeQUAD;
        char *movingMoleculeQUAD;    
        double distanceCutoff;           
        bool useInterfacePropensity;
        bool perResidueHydrophobicity;
        int numThreads;
      } PARAMS_IN;
    
    typedef struct
      {
       int id;
       double x, y, z;
       double r, h;
      } ATOM;
    
    typedef struct
      {
       double x, y, z;
       double w, h;     
      } QPOINT;
      
    typedef struct
      {
       double cx, cy, cz;
       double cr;
       bool leaf;
       int cPtr[ 8 ];
       int qPtsStartID, qPtsEndID;
      } QPOINTS_OCTREE_NODE;    
    
    typedef struct
      {
       double cx, cy, cz;
       double cr;
       bool leaf;
       int cPtr[ 8 ];
       int atomsStartID, atomsEndID;
      } ATOMS_OCTREE_NODE;    

    
    PARAMS_IN params;
    
    int numStaticAtoms, numMovingAtoms;
    
    ATOM *staticAtoms;    
    ATOM *movingAtoms;    

    int numStaticQPoints, numMovingQPoints;
                    
    QPOINT *staticQPoints;    
    QPOINT *movingQPoints;    
                
    Point *staticQPointsPG;
    Point *movingQPointsPG;
               
    PG *staticPG;
    PG *movingPG;
    
    int numStaticAtomsOctreeNodes;
    ATOMS_OCTREE_NODE *staticAtomsOctree;
    int staticAtomsOctreeRoot;

    int numMovingAtomsOctreeNodes;
    ATOMS_OCTREE_NODE *movingAtomsOctree;
    int movingAtomsOctreeRoot;

    int numStaticQPointsOctreeNodes;
    QPOINTS_OCTREE_NODE *staticQPointsOctree;
    int staticQPointsOctreeRoot;
    bool *staticQPointsOctreeFlags;

    int numMovingQPointsOctreeNodes;
    QPOINTS_OCTREE_NODE *movingQPointsOctree;
    int movingQPointsOctreeRoot;
    bool *movingQPointsOctreeFlags;    

    int *hydroIndex;
    int nHydroIndex;    
    
    int numThreads;

    pthread_mutex_t nodesLock;
    int curNode, maxNode;   
    
    double minRadius;
    int maxLeafSize;    
      
    void freeMemory( void );
    void setDefaults( void );
    void processQPoints( void );
    bool getParamsFromFile( PARAMS_IN *p, char *paramFile, bool atomsFromFile );
    double determinant( double *trans );
    void invert( double *trans, double *transI );    
    inline void transformPoint( double x, double y, double z, double *transMat, double *nx, double *ny, double *nz );
    inline void transformPoint( Point p, double *transMat, Point *np );
    double computeXlateForPG( int numStQPoints, QPOINT *stQPoints, int numMvQPoints, QPOINT *mvQPoints );
    void preprocessElementInformationTable( int **index, int *nIndex );
    int strcmp_nospace( char *s1, char *s2 );
//    float getHydrophobicity( char *atomName, char *residueName );
    float getHydrophobicity( char *atomName, char *residueName, int *index, int nIndex, bool useInterfacePropensity, bool perResidueHydrophobicity );
    bool readQPoints( char *qPtsFile, int *nQPoints, QPOINT **qPoints );
    bool copyAtomsFromArray( int numAtoms, double *atms, int *nAtoms, ATOM **atoms );
    bool readAtoms( char *atomsFile, int *nAtoms, ATOM **atoms );    
    void initFreeNodeServer( int numNodes );
    int nextFreeNode( void );
    void countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID, ATOM *atomsT, int *numNodes );
    int constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms, ATOMS_OCTREE_NODE *atomsOctree );
    bool buildAtomsOctree( int nAtoms, ATOM *atoms, int *numAtomsOctreeNodes, ATOMS_OCTREE_NODE **atomsOctree, int *atomsOctreeRoot );
    void countQPointsOctreeNodesAndSortQPoints( QPOINT *qPts, int qPtsStartID, int qPtsEndID, 
                                                QPOINT *qPtsT, int *numNodes );
    int constructQPointsOctree( int qPtsStartID, int qPtsEndID, QPOINT *qPts, QPOINTS_OCTREE_NODE *qPointsOctree );
    bool buildQPointsOctree( int nQPoints, QPOINT *qPoints, int *numQPointsOctreeNodes, QPOINTS_OCTREE_NODE **qPointsOctree, int *qPointsOctreeRoot );
    void assignHydrophobicityToQPoints( ATOMS_OCTREE_NODE *atomsOctree, int nodeA, ATOM *atoms,
                                        QPOINTS_OCTREE_NODE *qPointsOctree, int nodeQ, QPOINT *qPoints,
                                        double farDist, double rangeExt );
    void initOctreeFlags( int threadID );
    void collectPseudoGsol( int threadID, double *trans, double *transI, double *pGsol, 
                            double *pGsolHStaticPos, double *pGsolHStaticNeg, double *pGsolHMovingPos, double *pGsolHMovingNeg );
    void markPotentialQPoints( int threadID, int nodeS, int nodeM, double *trans );
                                                
 public: 
 
    pseudoGsol( char *paramFile, int nStAtoms, double *stAtoms, int nMvAtoms, double *mvAtoms, int nThreads );
    pseudoGsol( char *paramFile, int nThreads );    
    ~pseudoGsol( );

    void getPseudoGsol( int threadID, double *trans, double *pGsol,
                                    double *pGsolHStaticPos, double *pGsolHStaticNeg, double *pGsolHMovingPos, double *pGsolHMovingNeg );
    void getPseudoGsol( int threadID, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, 
                                    double *pGsolHMovingPos, double *pGsolHMovingNeg );
    void getPseudoGsol( int threadID, double *trans, double *pGsol, double *pGsolH );
    void getPseudoGsol( int threadID, double *pGsol, double *pGsolH );
    
    void printGsolParamters( FILE* fp );        
};


#endif
