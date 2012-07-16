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


#include <stdio.h>
#include "fastLJ.h"

int main( int argc, char *argv[ ] )
{
  if ( argc < 2 )
    {
      printError( "Input text file not specified!" );
      return 1;
    }
     
  fastLJ testLJ( argv[ 1 ] );   

  testLJ.setPrintStatus( true );     
  
  double fastP = 0, naiveP = 0;
  
  testLJ.computePotential( 0, &fastP );
  testLJ.computePotentialNaively( 0, &naiveP );
  
  printf( "\nfastP = %lf, naiveP = %lf\n", fastP, naiveP );
    
  return 0;
}
