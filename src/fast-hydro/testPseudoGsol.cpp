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


#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <stdio.h>

#include "../utils/utils.h"
#include "pseudoGsol.h"

bool readAtoms( char *atomsFile, int *nAtoms, double **atoms )
{
   FILE *fp;
   
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open atoms file (%s)!", atomsFile );
      return false;
     }
   
   *nAtoms = 0;
   double x, y, z, r, h;
          
   while ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf", &x, &y, &z, &r, &h ) == 5 ) ( *nAtoms )++; 

   fclose( fp );  
   
   ( *atoms ) = ( double * ) malloc( 5 * ( *nAtoms ) * sizeof( double ) );
   
   if ( *atoms == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }
     
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open atoms file (%s)!", atomsFile );
      return false;
     }

   for ( int i = 0; i < *nAtoms; i++ )   
     {          
      if ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf", &x, &y, &z, &r, &h ) != 5 )
        {
         printError( (char *)"Failed to read the atoms file (%s)!", atomsFile );
         return false;
        }
      
      ( *atoms )[ 5 * i + 0 ] = x;         
      ( *atoms )[ 5 * i + 0 ] = y;         
      ( *atoms )[ 5 * i + 0 ] = z;         
      ( *atoms )[ 5 * i + 0 ] = r;         
      ( *atoms )[ 5 * i + 0 ] = h;                                 
     }    
     
   fclose( fp );
         
   return true;
}


int main( int argc, char *argv[ ] )
{
  if ( argc < 4 )
    {
      printError( (char *)"Input files not specified ( #1: parameter file, #2: static atoms file, #3: moving atoms file )!" );
      return 1;
    }
    
  int nStaticAtoms, nMovingAtoms;  
  double *stAtoms, *mvAtoms;  
    
  if ( !readAtoms( argv[ 2 ], &nStaticAtoms, &stAtoms )
    || !readAtoms( argv[ 3 ], &nMovingAtoms, &mvAtoms ) )  
     return 1;
    
  pseudoGsol pGsol( argv[ 1 ], nStaticAtoms, stAtoms, nMovingAtoms, mvAtoms, 1 );  
     
  double trans[ ] = { 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 0, 0 };

  double u, v;
  
  pGsol.getPseudoGsol( 0, trans, &u, &v );

  printf( (char *)"\npseudo G_sol = %lf, pseudo G_sol with hydrophobicity = %lf\n", u, v );   
    
  return 0;
}
