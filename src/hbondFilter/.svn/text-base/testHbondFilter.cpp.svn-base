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
#include <hbondFilter/hbondFilter.h>

int main( int argc, char* argv[ ] )
{
  if ( argc < 10 )
    {
      printf( "Input files not specified (  staticPQR,  movingPQR,  staticPSF,  movingPSF,  staticMol2,  movingMol2,  rtfFile,  prmFile,  aprmFile )! \n" );
      return 1;
    }
    
  hbondFilter hf(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
  hf.initializeFilter();
       
  double trans[ ] = { 1, 0, 0, 1,
                      0, 1, 0, 0,
                      0, 0, 1, 0 };	// small translation
  double e;
  hf.getEnergy( trans, &e );

  printf( "\nhbondEnergy = %lf\n", e );   
    
  return 0;
}
