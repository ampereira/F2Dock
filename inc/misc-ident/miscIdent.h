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


#ifndef MISC_IDENT_H

#define MISC_IDENT_H

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <pthread.h>
#include <string.h>

#ifdef _WIN32
   #include <sys/types.h>
   #include <sys/timeb.h>
#else
   #include <sys/time.h>
#endif

#include "../utils/utils.h"
#include "../fast-clash/clashFilter.h"
#include "../fast-resCont/resContFilter.h"

//enum resType { NONE = 0, ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL };
enum resType { NONE = 0, ILE, VAL, LEU, PHE, CYS, MET, ALA, GLY, THR, SER, TRP, TYR, PRO, HIS, GLU, GLN, ASP, ASN, LYS, ARG };

typedef struct
  {
    int resNum, resID;
    int chainID;                
  } RESIDUE;


int getResidueID( char *resName );
bool countResidues( char *pqrFile, int res, int *count, int *total );

bool isGXY( int nRes, RESIDUE *res, int i );
bool isYXG( int nRes, RESIDUE *res, int i );

bool readResidues( char *pqrFile, int *nRes, RESIDUE **res, int *nChn, int **chn );
bool readAtomsAndResidues( char *pqrFile, int *nAtm, double **atm, int *nRes, RESIDUE **res, int *nChn, int **chn );
bool readAtomsOnly( char *pqrFile, int *nAtm, double **atm );
bool readAtomsWithResidueInfo( char *pqrFile, int *nAtm, double **atm );
bool readGlycines( char *pqrFile, int *nAtm, double **atm );

bool getTotalCharge( char *pqrFile, double *tCharge );

bool CDR_L1_start( RESIDUE *res, int i );
bool CDR_L1_end( RESIDUE *res, int i, int chainEnd );
bool CDR_L2_start( RESIDUE *res, int i );
bool CDR_L3_start( RESIDUE *res, int i );
bool CDR_L3_end( RESIDUE *res, int i, int chainEnd );
bool CDR_H1_start( RESIDUE *res, int i );
bool CDR_H1_end( RESIDUE *res, int i, int chainEnd );
bool CDR_H2_end( RESIDUE *res, int i, int chainEnd );
bool CDR_H3_start( RESIDUE *res, int i );
bool CDR_H3_end( RESIDUE *res, int i, int chainEnd );
bool CDR_L3_identified( RESIDUE *res, int c, int l, int s, int *id );
bool CDR_L2_L3_identified( RESIDUE *res, int c, int l, int s, int *id );
bool CDR_L_identified( RESIDUE *res, int c, int l, int *id );
bool CDR_H3_identified( RESIDUE *res, int c, int l, int s, int *id );
bool CDR_H2_H3_identified( RESIDUE *res, int c, int l, int s, int *id );
bool CDR_H_identified( RESIDUE *res, int c, int l, int *id );

bool isAntibody( char *pqrFile );
bool getAntibodyBindingSite( char *pqrFile, int *numAtoms, double **atm, bool low, bool high, int l1, int l2 );

bool countGXYandYXG( char *pqrFile, int *count, int *total, int *nGly );

bool initAntibodyClashFilter( char *staticPQR, char *movingPQR, bool low, bool high, int l1, int l2, clashFilter **cFilter );
bool initEnzymeClashFilter( char *staticPQR, char *movingPQR, clashFilter **cFilter );
bool initResContFilter( char *staticPQR, char *movingPQR, char *resContFile, resContFilter **cFilter );

#endif
