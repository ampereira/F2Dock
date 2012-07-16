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


#ifndef HBOND_FILTER_H
#define HBOND_FILTER_H

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <unistd.h> 
#include <math.h> 
#include <time.h>
#include <libmol/libmol.h>

//using namespace LibMol; 

#define MAX_ITER 1000

class hbondFilter
{
	private:
		int* fix1;
		int* fix2; 
		int nfix1, nfix2; 

        	double eface; 
		double efacv; 

		struct LibMol::prm* aprm;
		struct LibMol::atomgrp* staticAG;
		struct LibMol::atomgrp* movingAG;
		struct LibMol::agsetup* ags1;
		struct LibMol::agsetup* ags2;  
		int ndim1, ndim2;
	 	double* startag1;
	 	double* startag2;

	        LibMol::OCTREE staticOctree;
	        LibMol::OCTREE movingOctree;
		LibMol::OCTREE_PARAMS static_eng_params;
		LibMol::OCTREE_PARAMS moving_eng_params;
		LibMol::OCTREE_PARAMS complex_eng_params;

		double staticHbondEnergy;
		double movingHbondEnergy;
		double complexHbondEnergy;

		void read_fix(char *ffile, int *nfix, int **fix);

		
	public:
		hbondFilter(char* staticPQR, char* movingPQR, char* staticPSF, char* movingPSF, char* staticMol2, char* movingMol2, char* rtfFile, char* prmFile, char* aprmFile);
		~hbondFilter();
		bool initializeFilter();
		bool getEnergy(double* trans, double* en);
};

#endif
