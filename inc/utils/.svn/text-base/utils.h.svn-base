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


#ifndef UTILS_H

#define UTILS_H

#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <pthread.h>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include <stdarg.h>
#include <time.h>

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

#define transform( ox, oy, oz, M, nx, ny, nz ) {                                                                      \
                                                 nx = ( ox ) * M[ 0 ] + ( oy ) * M[ 1 ] + ( oz ) * M[  2 ] + M[  3 ]; \
                                                 ny = ( ox ) * M[ 4 ] + ( oy ) * M[ 5 ] + ( oz ) * M[  6 ] + M[  7 ]; \
                                                 nz = ( ox ) * M[ 8 ] + ( oy ) * M[ 9 ] + ( oz ) * M[ 10 ] + M[ 11 ]; \
                                               }

#ifndef M_PI
   #define M_PI 3.1415926535897932384626433832795
#endif

#ifndef INV_SQRT_TWO
   #define INV_SQRT_TWO 0.70710678118654752440084436210485
#endif

void printError( char *format, ... );
void flushPrint( char *format, ... );
void f_printf( FILE *fp, char *format, ... );
double getTime( void );
int skipWhiteSpaces( char *buf, int i );
int skipInitial( char *s1, char *s2, char *p );
int getInt( char *buf, int i, int *v );
int getDouble( char *buf, int i, double *v );
bool getDoublesInRange( char *buf, int i1, int i2, double *v );
int getAlphaString( char *buf, int i, char *s );
int getString( char *buf, int i, char *s );

#endif
