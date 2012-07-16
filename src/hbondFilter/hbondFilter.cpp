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


#include <hbondFilter/hbondFilter.h>

using namespace LibMol; 


hbondFilter::hbondFilter(char* staticPQR, char* movingPQR, char* staticPSF, char* movingPSF, char* staticMol2, char* movingMol2, char* rtfFile, char* prmFile, char* aprmFile)
{
        eface=1.0; 
	efacv=1.0;

	aprm = read_prm (aprmFile, "0.0.6");        
	
	staticAG = read_pdb (staticPQR, aprm); 
	
	fix1 = (int*)mymalloc(staticAG->natoms*sizeof(int)); 
        read_ff_charmm(staticPSF, prmFile, rtfFile, staticAG); 

        if ( !read_hybridization_states_from_mol2( staticMol2, staticPQR, staticAG ) )
	{
		printf("Error reading hybridization states\n");
		return; 
	}       
          
        fix_acceptor_bases( staticAG, aprm );

	//Read in fixed part of the molecule 
	read_fix(staticPQR, &nfix1, &fix1); 
	fixed_init(staticAG); 
        fixed_update(staticAG, nfix1, fix1); 
        zero_grads(staticAG); 
	fill_ingrp(staticAG);

        //Setup nblist information
        ags1 = (agsetup*)mymalloc(sizeof(struct agsetup)); 
        init_nblst(staticAG,ags1); 
	update_nblst(staticAG,ags1); 
	
        assign_combined_residue_sequence( staticAG );	
        ndim1 = staticAG->nactives*3; 
 	startag1 = (double*)mymalloc(ndim1*sizeof(double) );
	ag2array(startag1,staticAG);


	movingAG = read_pdb (staticPQR, aprm); 
	
	fix2 = (int*)mymalloc(movingAG->natoms*sizeof(int)); 
        read_ff_charmm(movingPSF, prmFile, rtfFile, movingAG); 

        if ( !read_hybridization_states_from_mol2( movingMol2, movingPQR, movingAG ) )
	{
		printf("Error reading hybridization states 2\n");
		return; 
	}       
          
        fix_acceptor_bases( movingAG, aprm );

	//Read in fixed part of the molecule 
	read_fix(movingPQR, &nfix2, &fix2); 
	fixed_init(movingAG); 
        fixed_update(movingAG, nfix2, fix2); 
        zero_grads(movingAG); 
	fill_ingrp(movingAG);

        //Setup nblist information
        ags2 = (agsetup*)mymalloc(sizeof(struct agsetup)); 
        init_nblst(movingAG,ags2); 
	update_nblst(movingAG,ags2); 
	
        assign_combined_residue_sequence( movingAG );	
        ndim2 = movingAG->nactives*3; 
 	startag2 = (double*)mymalloc(ndim2*sizeof(double) );
	ag2array(startag2,movingAG);
}

hbondFilter::~hbondFilter()
{
	free(fix1);
	free(fix2);
	free(startag1);
	free(startag2);
}


bool hbondFilter::initializeFilter()
{
	//building octrees        
        if ( !build_octree( &staticOctree, max( 60, ags1->nblst->npairs / staticAG->natoms ), 6, 1.0, staticAG ) )        
	{
		printf("Failed to build octree!\n");
		return false;
	}

        if ( !build_octree( &movingOctree, max( 60, ags2->nblst->npairs / movingAG->natoms ), 6, 1.0, movingAG ) )        
	{
		printf("Failed to build octree!\n");
		return false;
	}
             
        complex_eng_params.octree_static = &staticOctree;
        complex_eng_params.octree_moving = &movingOctree;
        moving_eng_params.octree_static = &movingOctree;
        moving_eng_params.octree_moving = &movingOctree;
        static_eng_params.octree_static = &staticOctree;
        static_eng_params.octree_moving = &staticOctree;

        complex_eng_params.ags = static_eng_params.ags = moving_eng_params.ags = ags1; 
        complex_eng_params.eps = static_eng_params.eps = moving_eng_params.eps = eface;   
        complex_eng_params.hdist_cutoff = static_eng_params.hdist_cutoff = moving_eng_params.hdist_cutoff = MAX_R;
        complex_eng_params.fixed_cull = static_eng_params.fixed_cull = moving_eng_params.fixed_cull = 1;
        complex_eng_params.trans = static_eng_params.trans = moving_eng_params.trans = NULL;
        complex_eng_params.engcat = static_eng_params.engcat = moving_eng_params.engcat = NULL;        

        complex_eng_params.proc_func_params = &complex_eng_params;
        static_eng_params.proc_func_params = &static_eng_params;
        moving_eng_params.proc_func_params = &moving_eng_params;

	//computing individual hbond energies
	staticHbondEnergy = octree_accumulation_excluding_far( &staticOctree, &staticOctree, MAX_R, MAX_R, 1, NULL, ( void * ) ( &static_eng_params ), hbondeng_octree_single_mol );
	movingHbondEnergy = octree_accumulation_excluding_far( &movingOctree, &movingOctree, MAX_R, MAX_R, 1, NULL, ( void * ) ( &moving_eng_params ), hbondeng_octree_single_mol );

	return true;
}


bool hbondFilter::getEnergy(double* trans, double* en)
{
	complexHbondEnergy = octree_accumulation_excluding_far( &staticOctree, &movingOctree, MAX_R, MAX_R, 1, trans, ( void * ) ( &complex_eng_params ), hbondeng_octree_single_mol );
	
	*en = complexHbondEnergy - movingHbondEnergy - staticHbondEnergy;

	return true;
}


void hbondFilter::read_fix(char *ffile, int *nfix, int **fix) 
{ 
	int linesz=91; 
	char *buffer= (char*)mymalloc(sizeof(char)*linesz); 
	*nfix=0; 
	FILE* fp = fopen (ffile, "r"); 
	while (fgets(buffer, linesz-1, fp)!=NULL) 
	{ 
		if(!strncmp(buffer,"ATOM",4))(*nfix)++; 
	} 
	fclose(fp); 
	*fix=(int*)mymalloc(*nfix*sizeof(int)); 
	fp = fopen (ffile, "r"); 
	int na=0; 
	while(fgets(buffer, linesz-1, fp)!=NULL) 
	{ 
		if(!strncmp(buffer,"ATOM",4)) 
		{ 
			(*fix)[na]=atoi(buffer+4)-1; 
			na++; 
		} 
	} 
	free(buffer); 
	fclose(fp); 
} 
