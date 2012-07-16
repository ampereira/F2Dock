/**-----------------------------------------------------------------------------                                                                                                  
**                                                                                                                                                                               
**  Copyright (C) : Structural Bioinformatics Laboratory, Boston University.                                                                                                                        
**                                                                                                                                                                               
**  This software was developed at the Boston University 2006-2011, by                                                                                                      
**  Structural Bioinformatics Laboratory, as part of NIH funded research.                                                                                                                      
**                                                                                                                                                                               
**  Explicit permission is hereby granted to US Universities and US                                                                                                     
**  Government supported Research Institutions to copy and modify this                                                                                                           
**  software for educational and research purposes, provided copies include                                                                                                      
**  this notice. This software (or modified copies thereof) may not be                                                                                                           
**  distributed to any other institution without express permission from the                                                                                                     
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 **  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <libmol/libmol.h>

namespace LibMol{


File_Type file_ext (const char* path)
{
	char* ri = (char *)rindex (path, '.');

	if (ri == NULL)
		return FILE_UNKNOWN;
	if (strncmp (ri, ".pdb", 4) == 0)
		return FILE_PDB;
	if (strncmp (ri, ".xyz", 4) == 0)
		return FILE_XYZ;
	if (strncmp (ri, ".ms", 3) == 0)
		return FILE_MS;
	if (strncmp (ri, ".dx", 3) == 0)
		return FILE_DX;
	if (strncmp (ri, ".pot", 4) == 0)
		return FILE_POT;

    return FILE_UNKNOWN;
}

struct atomgrp* read_file_atomgrp (const char* path, struct prm* prm, float msur_k)
{
	struct atomgrp* ag = NULL;
	if (file_ext (path) == FILE_PDB)
	{
		ag=read_pdb (path, prm);
		msur (ag, prm, msur_k);
		return ag;
	}

	if (file_ext (path) == FILE_XYZ)
		return read_oldxyz (path, prm);

	if (file_ext (path) == FILE_MS)
		return read_ms (path, prm);

	// file type unknown
	fprintf (stderr, "file ext of %s is not a recognized file ext\n", path);
	exit (EXIT_FAILURE);

	return NULL;
}
/*
struct grid* read_file_grid (const char* path)
{
	if (file_ext (path) == FILE_DX)
		return opendx_read_grid (path);

	if (file_ext (path) == FILE_POT)
		return read_pot_grid (path);

	return NULL;
}*/

void fprint_file_atomgrp (const char* path, struct atomgrp* ag, struct prm* prm)
{
	if (file_ext (path) == FILE_PDB)
		return fprint_pdb (ag, prm, path);

	if (file_ext (path) == FILE_XYZ)
		return fprint_oldxyz (ag, prm, path);

	if (file_ext (path) == FILE_MS)
		return fprint_ms (ag, prm, path);

	// file type unknown
	fprintf (stderr, "file ext of %s is not a recognized file ext\n", path);
	exit (EXIT_FAILURE);
}

void read_mod_vdw(char *mfile, int *nmod, int **mod, double **modeps, double **modrminh)
{
   int linesz=91;
   char *buffer=(char*)mymalloc(sizeof(char)*linesz);
   *nmod=0;
   FILE* fp = myfopen (mfile, "r");
   while (fgets(buffer, linesz-1, fp)!=NULL)
   {
      if(!strncmp(buffer,"ATOM",4))(*nmod)++;
   }
   fclose(fp);
   *mod=(int*)mymalloc(*nmod*sizeof(int));
   *modeps=(double*)mymalloc(*nmod*sizeof(double));
   *modrminh=(double*)mymalloc(*nmod*sizeof(double));
   fp = fopen (mfile, "r");
   int na=0;
   while(fgets(buffer, linesz-1, fp)!=NULL)
   {
      if(!strncmp(buffer,"ATOM",4))
      {
         (*mod)[na]=atoi(buffer+4)-1;
	 (*modeps)[na]=atof(buffer+54);
	 (*modrminh)[na]=atof(buffer+60);
         na++;
      }
   }
   free(buffer);
   fclose(fp);
}

};
