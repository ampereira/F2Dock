#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <unistd.h> 
#include <math.h> 
#include <time.h>
#include <libmol/libmol.h>

using namespace LibMol; 

#define MAX_ITER 1000
 
void read_fix(char *ffile, int *nfix, int **fix) 
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
 //  free(fix); 
   fclose(fp); 
} 


void print_short_args (char* app_name) 
{ 
        fprintf (stderr, "usage: %s PDBFILE OUTPUTPDB PSF PRM RTF FIXEEDATOMS ATOMSPRMFILE MOL2FILE\n", app_name); 
} 
 
 
int  main(int argc, char* argv[])  
{ 
	char* app_name = argv[0]; 
	if (argc != 9) 
	{ 
		print_short_args (app_name); 
		exit (EXIT_FAILURE); 
	} 
 
	char* ifile; 
        char* ofile;
	char* ofilefull;
	char* psffile; 
        char* prmfile; 
        char* rtffile; 
	char* ffile; 
	char* aprmfile;
	char* mol2file;
	 
	int *fix; 
        int nfix; 
	 
	ifile    = argv[1]; 
        ofile    = argv[2]; 
	psffile  = argv[3]; 
        prmfile  = argv[4]; 
        rtffile  = argv[5]; 
 	ffile    = argv[6]; 
 	aprmfile = argv[7]; 
 	mol2file = argv[8];  	

        double eface=1.0; 
	double efacv=1.0; 
        ofilefull=(char*)mymalloc(100*sizeof(char));
        struct prm* aprm = read_prm (aprmfile, "0.0.6");        

	//Read in file and it's forcefield parameters
	struct atomgrp* ag = read_pdb (ifile, aprm); 
	
	fix=(int*)mymalloc(ag->natoms*sizeof(int)); 
        read_ff_charmm(psffile, prmfile, rtffile, ag); 

        if ( !read_hybridization_states_from_mol2( mol2file, ifile, ag ) )
          {
	    exit (EXIT_FAILURE);             
          }       
          
        fix_acceptor_bases( ag, aprm );

	//Read in fixed part of the molecule 
	read_fix(ffile, &nfix, &fix); 
	
        fixed_init(ag); 
        fixed_update(ag,nfix,fix); 
	//Finish structure initialisation 
        zero_grads(ag); 
	fill_ingrp(ag);

        //Setup nblist information
	struct agsetup* ags; 
        ags=(agsetup*)mymalloc(sizeof(struct agsetup)); 
        init_nblst(ag,ags); 
        //Update nblst 
	update_nblst(ag,ags); 
	
        assign_combined_residue_sequence( ag );	
        int ndim = ag->nactives*3; 
 	double* startag = (double*)mymalloc(ndim*sizeof(double) );
	//Storing current active list of atoms in array
	ag2array(startag,ag);
                
        OCTREE octree;
        
        if ( !build_octree( &octree, max( 60, ags->nblst->npairs / ag->natoms ), 6, 1.0, ag ) )        
          {
            print_error( "Failed to build octree!" );
            return 1;
          }        

        OCTREE_PARAMS eng_params;        
        eng_params.ags = ags; 
        eng_params.eps = eface;
        
        eng_params.octree_static = &octree;
        eng_params.octree_moving = &octree;
        eng_params.hdist_cutoff = MAX_R;
        eng_params.fixed_cull = 1;
        eng_params.trans = NULL;
        eng_params.engcat = NULL;        
        eng_params.proc_func_params = &eng_params;

        printf( "ags->nblst->nbcof = %lf\n", ags->nblst->nbcof );
        
        int k;
 
        double start, end, elapsed; 
        double en; 
       
        start = clock();         
        for ( k = 0; k < MAX_ITER; k++ )
          {
            zero_grads( ag );           
            en = 0;
            hbondeng( ag, &en, ags->nblst ); 
          }  
        end = clock();
        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("nblist_hbond_energy = %lf, time = %lf sec\n", en, elapsed/MAX_ITER);

        double *engcat = alloc_categorized_hbondeng( );

        init_categorized_hbondeng( engcat );
        hbondengcat(ag, engcat, ags->nblst);
        
        double bb_bb_sr, bb_bb_lr, bb_sc, sc_sc;        
        get_categorized_hbondeng( engcat, &bb_bb_sr, &bb_bb_lr, &bb_sc, &sc_sc );
        printf( "nblist_hond_energy_cat = < BB-BB-SR: %lf, BB-BB-LR: %lf, BB-SC: %lf, SC-SC: %lf >\n", bb_bb_sr, bb_bb_lr, bb_sc, sc_sc );

        start = clock();     
        for ( k = 0; k < MAX_ITER; k++ )
          {
            zero_grads( ag );                     
            en = octree_accumulation_excluding_far( &octree, &octree, MAX_R, MAX_R, 1, NULL, ( void * ) ( &eng_params ), hbondeng_octree_single_mol );
          }  
        end = clock();
        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("octree_hbond_energy = %lf, time = %lf sec\n", en, elapsed/MAX_ITER);

        init_categorized_hbondeng( engcat );	
        eng_params.engcat = engcat;	
        en = octree_accumulation_excluding_far( &octree, &octree, MAX_R, MAX_R, 1, NULL, ( void * ) ( &eng_params ), hbondeng_octree_single_mol );
        
        get_categorized_hbondeng( engcat, &bb_bb_sr, &bb_bb_lr, &bb_sc, &sc_sc );
        printf( "octree_hond_energy_cat = < BB-BB-SR: %lf, BB-BB-LR: %lf, BB-SC: %lf, SC-SC: %lf >\n", bb_bb_sr, bb_bb_lr, bb_sc, sc_sc );


        FLOW_STRUCT fs;
        init_flow_struct( &fs );
        ag->flow_struct = ( void * ) ( &fs );	

//        start = clock();         
//        for ( k = 0; k < MAX_ITER; k++ )
//          {
            zero_grads( ag );           
            en = 0;
            flow_hbondeng( ag, &en, ags->nblst ); 
//          }  
//        end = clock();
//        elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
//	  printf("flow_energy = %lf, time = %lf sec\n", en, elapsed/MAX_ITER);
	printf("flow_hbond_energy = %lf\n", en);
        
        free_flow_struct( &fs );

        
        return 0; 
} 
 

