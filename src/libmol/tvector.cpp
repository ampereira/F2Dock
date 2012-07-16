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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

mol_tvector *
mol_tvector_create ()
{
	mol_tvector *v = (mol_tvector*)_mol_malloc (sizeof (mol_tvector));
	return v;
}

void
mol_tvector_destroy (mol_tvector *v)
{
	free (v);
}

void
mol_tvector_init (mol_tvector *v, float x, float y, float z)
{
	v->X = x;
	v->Y = y;
	v->Z = z;
}

void
mol_tvector_copy (mol_tvector* src, mol_tvector* dst)
{
	assert (src != NULL);
	assert (dst != NULL);

	memcpy ( dst, src, sizeof(mol_tvector) );
}

mol_tvector *
mol_tvector_clone (mol_tvector *a)
{
	mol_tvector *b;

	assert (a != NULL);
	b = mol_tvector_create ();
	mol_tvector_copy (a, b);

	return b;
}

void
mol_tvector_add_inplace (mol_tvector *a, mol_tvector *b)
{
	a->X += b->X;
	a->Y += b->Y;
	a->Z += b->Z;
}

mol_tvector *
mol_tvector_add (mol_tvector *a, mol_tvector *b)
{
	mol_tvector *c;

	c = mol_tvector_clone (a);
	mol_tvector_add_inplace (c, b);

	return c;
}


int
mol_tvector_normalize (mol_tvector *a)
{
	float d = a->X * a->X + a->Y * a->Y + a->Z * a->Z;

	if ( d == 0 ) return 0;

	d = sqrtf( d );

	a->X /= d;
	a->Y /= d;
	a->Z /= d;

	return 1;
}


void
mol_tvector_print (mol_tvector* t)
{
	printf ("%.4f, %.4f, %.4f\n", t->X, t->Y, t->Z);
}

};
