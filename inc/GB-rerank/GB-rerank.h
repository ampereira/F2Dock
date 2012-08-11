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


#ifndef GB_RERANK_H

#define GB_RERANK_H

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

#include "fast-GB/fastBornRadius.h"
#include "fast-GB/fastGpol.h"
#include "utils/utils.h"
#include "GB-rerank/pairingHeap.h"
#include "PG-range/PG.h"

typedef struct
  {
    char *staticMoleculePQR;
    char *movingMoleculePQR;
    
    char *staticMoleculeQUAD;
    char *movingMoleculeQUAD;

    char *F2DockOutputFile;
    char *rerankedOutputFile;
    
    double F2DockScoreWeight;
    double GpolWeight; 
    double GnonpolWeight;
    
    double distanceCutoff;   
    
    int numThreadsBR, numThreadsGpol;
    double epsilonBR, epsilonGpol;
    bool useApproxMath;
    
    int numBands;
    int *bands;
    
    int numSol;
    
  } PARAMS_IN;


typedef struct
  {
    int numAtoms, numQPoints;    
    
    double *atoms;
    double *qPoints;    
    
    double Gpol;
            
  } MOLECULE_INFO;


typedef struct
  {
    char *sol;        
    int F2DockRank;
    double rmsd;
    double DelGpol;
    double areaProp;
    double newScore;
    
    int nextIndex;
                
  } SOLUTION_INFO;

#endif
