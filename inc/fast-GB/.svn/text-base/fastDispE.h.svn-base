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


#ifndef FAST_DISP_E_H

#define FAST_DISP_E_H

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

#include "fastBornRadius.h"
#include "../utils/utils.h"


namespace fastGB {

class fastDispE
{
 private:

   typedef struct
     {
       char *staticMoleculeQUAD;
       char *movingMoleculeQUAD;
        
       int numThreads;
       int numThreadsPerThread;
       
       double epsilonBR;
        
     } PARAMS_IN;

   PARAMS_IN params;

   typedef struct
     {
       int numAtoms, numQPoints;    
        
       double *atoms;
       double *qPoints;    
        
       double dispE;
                
     } MOLECULE_INFO;

   MOLECULE_INFO staticMol, movingMol;  
   
   int numAtoms, numQPoints;
   
   double *atoms;
   double *qPoints;

   void freeMemory( void );
   void getAtomsQuadsAndDispE( char *quadFile, int numAtoms, double *atomsPQR, MOLECULE_INFO *mol );
   bool getParamsFromFile( PARAMS_IN *p, char *paramFile );
   void copyStaticMoleculeToAllThreads( void );
   void transformAndCopyMovingMolecule( int threadID, double *transMat );
   
 public:

   fastDispE( char *paramFile, int numStaticAtoms, double *staticAtomsPQR, int numMovingAtoms, double *movingAtomsPQR,
              int numThreads, int numThreadsPerThread );
   ~fastDispE( );

   double getDispE( int threadID, double *transMat );
   double getDispE( int threadID );
   
   double getDelDispE( int threadID, double *transMat );
   double getDelDispE( int threadID );   
};

};

#endif

