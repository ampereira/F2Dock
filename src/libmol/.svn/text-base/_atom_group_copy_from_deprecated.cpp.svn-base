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
#include <assert.h>
#include <stdlib.h>

#include <libmol/libmol.h>

namespace LibMol{

void
_mol_atom_create_bond_indices (mol_atom* a, int nbondis)
{
	assert (a != NULL);

	a->nbondis = nbondis;
	if (a->nbondis > 0)
		a->bondis = (int*) _mol_malloc (a->nbondis * sizeof (int));
	else
		a->bondis = NULL;

	return;
}

void
_mol_atom_copy_bond_ptrs_to_bond_indices (mol_atom* atom, int nbonds, mol_bond* bonds)
{
	int bondssi, bondsi;
	mol_bond** bondss;

	atom->nbondis = atom->nbonds;
	bondss = atom->bonds;

	for (bondssi = 0; bondssi < atom->nbonds; bondssi++)
	{
		int target_count = 0;
		mol_bond* bond_target = bondss[bondssi];

		for (bondsi = 0; bondsi < nbonds; bondsi++)
		{
			if (&bonds[bondsi] == bond_target)
			{
				atom->bondis[bondssi] = bondsi;
				target_count++;
			}
		}

		assert (target_count == 1);
	}
}

void
_mol_bond_copy_atom_ptrs_to_atom_indices (mol_bond* b)
{
	b->ai = b->a0->ingrp;
	b->aj = b->a1->ingrp;
}

void
_mol_atom_group_copy_from_deprecated (mol_atom_group* ag)
{
	int atomsi, bondsi;

	for (atomsi = 0; atomsi < ag->natoms; atomsi++)
	{
		_mol_atom_copy_bond_ptrs_to_bond_indices (&ag->atoms[atomsi],
				ag->nbonds, ag->bonds);
	}

	for (bondsi = 0; bondsi < ag->nbonds; bondsi++)
	{
		_mol_bond_copy_atom_ptrs_to_atom_indices (&ag->bonds[bondsi]);
	}
}
};
