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
#ifndef _MOL_ATOM_GROUP_COPY_FROM_DEPRECATED_H_
#define _MOL_ATOM_GROUP_COPY_FROM_DEPRECATED_H_

namespace LibMol{

void
_mol_atom_create_bond_indices (mol_atom* a, int nbondis);

void
_mol_atom_copy_bond_ptrs_to_bond_indices (mol_atom* atom, int nbonds, mol_bond*
		bonds);

void
_mol_bond_copy_atom_ptrs_to_atom_indices (mol_bond* b);

void
_mol_atom_group_copy_from_deprecated (mol_atom_group* ag);

};
#endif
