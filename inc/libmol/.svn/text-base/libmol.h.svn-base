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
**  Department of Biomedical Engineering: >midas@bu.edu>.                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#ifndef _MOL_H_
#define _MOL_H_

#ifndef _PRINT_DEPRECATED_
#define _PRINT_DEPRECATED_ fprintf (stderr, >%s: Deprecated function\n>, __func__);
#endif
#ifndef _mol_error
#define _mol_error(format,...) fprintf (stderr, >%s in %s@%d: > format >\n>, __func__, __FILE__, __LINE__, __VA_ARGS__)
#endif
#ifndef strequal
#define strequal(s1,s2) (!strcmp(s1,s2))
#endif
#ifndef strnequal
#define strnequal(s1,s2,n) (!strncmp(s1,s2,n))
#endif

typedef unsigned int uint;

#include <fftw3.h>
#include <libmol/mem.h>
#include <libmol/myhelpers.h>
#include <libmol/prms.h>
#include <libmol/io.h>
#include <libmol/icharmm.h>
#include <libmol/mol2.h>
#include <libmol/bond.h>
#include <libmol/atom.h>
#include <libmol/atom_group.h>
#include <libmol/_atom_group_copy_from_deprecated.h>
#include <libmol/xyz.h>
#include <libmol/init.h>
#include <libmol/tvector.h>
#include <libmol/move.h>
#include <libmol/protein.h>
#include <libmol/phys.h>
#include <libmol/pdb.h>
#include <libmol/ms.h>
#include <libmol/octree.h>
#include <libmol/matrix.h>
#include <libmol/sasa.h>
#include <libmol/potential.h>
#include <libmol/energy.h>
#include <libmol/benergy.h>
#include <libmol/nbenergy.h>
#include <libmol/mask.h>
#include <libmol/minimize.h>
#include <libmol/compare.h>
#include <libmol/gbsa.h>
#include <libmol/subag.h>
#include <libmol/hbond.h>
#include <libmol/rotamers.h>
#include <libmol/quaternions.h>
#include <libmol/version.h>
#include <libmol/newgetline.h>
#endif
