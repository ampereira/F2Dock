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
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

#define linesize 128

void read_ff_charmm(const char *psffile, char *prmfile, char *rtffile,
		    struct atomgrp *ag)
{
	int i, i1, i2, i3, i4;
	int nat, nb, nang, ndih0, ndih, nimp, ndon, nacc;
	int *atoms, *bonds, *angs, *dihs, *imps, *dons, *acpts;
	char *atnam;
	char **atom_names;	//CA, CB, etc.
	float *lbond, *kbond, *fang, *kang, *fdih, *kdih;
	float *fimp, *kimp, *cha, *siga, *epa, *acevolumes;
	int *pdih, *tdih, ti;
	int nres, *ires;
	char *res_name;
	char **dres = (char**)_mol_malloc(ag->natoms * sizeof(char *));

	read_psf(psffile, &nat, &atoms, &atom_names, &nb, &bonds, &nang, &angs,
		 &ndih0, &dihs, &nimp, &imps, &cha, &nres, &ires, dres,
		 &ndon, &dons, &nacc, &acpts);

// residues
	ag->nres = nres;
	if (nres > 0) {
		ag->iares = (int*)_mol_malloc((nres + 1) * sizeof(int));
		ag->idres = (char**)_mol_malloc(nres * sizeof(char *));
		ag->res_type = (enum mol_res_type*)_mol_malloc(nres * sizeof(enum mol_res_type));
	}
	for (i = 0; i < nres; i++) {
		ag->iares[i] = ires[i];
		ag->idres[i] = dres[i];
	}
	ag->iares[nres] = nat;

	for (i = 0; i < nres; i++) {
		res_name = (ag->idres[i]) + 10; //location of residue name on psf line

		if (strncmp(res_name, "ALA ", 4) == 0)
			ag->res_type[i] = ALA;
		else if (strncmp(res_name, "ARG ", 4) == 0)
			ag->res_type[i] = ARG;
		else if (strncmp(res_name, "ARGN", 4) == 0)
			ag->res_type[i] = ARGN;
		else if (strncmp(res_name, "ASN ", 4) == 0)
			ag->res_type[i] = ASN;
		else if (strncmp(res_name, "ASP ", 4) == 0)
			ag->res_type[i] = ASP;
		else if (strncmp(res_name, "ASPH", 4) == 0)
			ag->res_type[i] = ASPH;
		else if (strncmp(res_name, "CYS ", 4) == 0)
			ag->res_type[i] = CYS;
		else if (strncmp(res_name, "GLN ", 4) == 0)
			ag->res_type[i] = GLN;
		else if (strncmp(res_name, "GLU ", 4) == 0)
			ag->res_type[i] = GLU;
		else if (strncmp(res_name, "GLUH", 4) == 0)
			ag->res_type[i] = GLUH;
		else if (strncmp(res_name, "GLY ", 4) == 0)
			ag->res_type[i] = GLY;
		else if (strncmp(res_name, "HIS ", 4) == 0)
			ag->res_type[i] = HIS;
		else if (strncmp(res_name, "ILE ", 4) == 0)
			ag->res_type[i] = ILE;
		else if (strncmp(res_name, "LEU ", 4) == 0)
			ag->res_type[i] = LEU;
		else if (strncmp(res_name, "LYS ", 4) == 0)
			ag->res_type[i] = LYS;
		else if (strncmp(res_name, "LYSN", 4) == 0)
			ag->res_type[i] = LYSN;
		else if (strncmp(res_name, "MET ", 4) == 0)
			ag->res_type[i] = MET;
		else if (strncmp(res_name, "PHE ", 4) == 0)
			ag->res_type[i] = PHE;
		else if (strncmp(res_name, "PRO ", 4) == 0)
			ag->res_type[i] = PRO;
		else if (strncmp(res_name, "SER ", 4) == 0)
			ag->res_type[i] = SER;
		else if (strncmp(res_name, "SEP ", 4) == 0)
			ag->res_type[i] = SEP;
		else if (strncmp(res_name, "THR ", 4) == 0)
			ag->res_type[i] = THR;
		else if (strncmp(res_name, "TRP ", 4) == 0)
			ag->res_type[i] = TRP;
		else if (strncmp(res_name, "TYR ", 4) == 0)
			ag->res_type[i] = TYR;
		else if (strncmp(res_name, "VAL ", 4) == 0)
			ag->res_type[i] = VAL;
		else if (strncmp(res_name, "HSC ", 4) == 0)
			ag->res_type[i] = HSC;
		else if (strncmp(res_name, "HSD ", 4) == 0)
			ag->res_type[i] = HSD;
		else if (strncmp(res_name, "ACE ", 4) == 0)
			ag->res_type[i] = ACE;
		else if (strncmp(res_name, "PCA ", 4) == 0)
			ag->res_type[i] = PCA;
		else if (strncmp(res_name, "HYL ", 4) == 0)
			ag->res_type[i] = HYL;
		else if (strncmp(res_name, "HYP ", 4) == 0)
			ag->res_type[i] = HYP;
		else if (strncmp(res_name, "HSE ", 4) == 0)
			ag->res_type[i] = HSE;
		else if (strncmp(res_name, "ORN ", 4) == 0)
			ag->res_type[i] = ORN;
		else if (strncmp(res_name, "PEN ", 4) == 0)
			ag->res_type[i] = PEN;
		else if (strncmp(res_name, "ALB ", 4) == 0)
			ag->res_type[i] = ALB;
		else if (strncmp(res_name, "ABU ", 4) == 0)
			ag->res_type[i] = ABU;
		else if (strncmp(res_name, "ST2 ", 4) == 0)
			ag->res_type[i] = ST2;
		else if (strncmp(res_name, "TIP3", 4) == 0)
			ag->res_type[i] = TIP3;
		else if (strncmp(res_name, "OH2 ", 4) == 0)
			ag->res_type[i] = OH2;
		else if (strncmp(res_name, "HOH ", 4) == 0)
			ag->res_type[i] = HOH;
		else if (strncmp(res_name, "WAT ", 4) == 0)
			ag->res_type[i] = WAT;
		else
			ag->res_type[i] = UNK;
	}

	atnam = get_atnam_rtf(rtffile, atoms, nat);
	ndih = ndih0;
	read_par(prmfile, nat, atnam, nb, bonds, &kbond, &lbond,
		 nang, angs, &kang, &fang,
		 &ndih, dihs, &kdih, &fdih, &pdih, &tdih,
		 nimp, imps, &kimp, &fimp, &epa, &siga, &acevolumes);

//        printf("DIHEDRAL ndih0, ndih %d %d\n", ndih0, ndih);
	if ((ag->natoms) != nat) {
		fprintf(stderr,
			"different number of atoms after reading protein structure file (PSF)\n");
		exit(EXIT_FAILURE);
	}
// vdw and charges, types and acevolumes
	for (i = 0; i < nat; i++) {
		ag->atoms[i].atom_ftypen = atoms[i];
		ag->atoms[i].name = atom_names[i];
		if (acevolumes[i] < 0.0) {
			printf("Warning: Ace Volume not set for atom %d", i);
			printf("Setting volume to %f", 18.0);
			ag->atoms[i].acevolume = 18.0;
		} else {
			ag->atoms[i].acevolume = acevolumes[i];
		}

		ag->atoms[i].chrg = cha[i];
		if (epa[2 * i] <= 0)
			ag->atoms[i].eps = sqrt(-epa[2 * i]);
		else {
			printf("error on pdb atom %d: positive vdw eps for atom type %.4s",
			     i, &atnam[4 * i]);
			exit(EXIT_FAILURE);
		}
		ag->atoms[i].rminh = siga[2 * i];

		if (epa[2 * i + 1] <= 0) {
			ag->atoms[i].eps03 = sqrt(-epa[2 * i + 1]);
			ag->atoms[i].rminh03 = siga[2 * i + 1];
		} else {
			ag->atoms[i].eps03 = ag->atoms[i].eps;
			ag->atoms[i].rminh03 = ag->atoms[i].rminh;
		}
		ag->atoms[i].nbonds = 0;
		ag->atoms[i].nangs = 0;
		ag->atoms[i].ntors = 0;
		ag->atoms[i].nimps = 0;
	}
// bonds
	ag->nbonds = nb;
	ag->bonds = (struct atombond*)_mol_malloc(nb * sizeof(struct atombond));
	for (i = 0; i < nb; i++) {
		i1 = bonds[2 * i] - 1;
		i2 = bonds[2 * i + 1] - 1;
		ag->bonds[i].a0 = &(ag->atoms[i1]);
		ag->bonds[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i1].nbonds)++;
		(ag->atoms[i2].nbonds)++;
		ag->bonds[i].l0 = lbond[i];
		ag->bonds[i].k = kbond[i];
	}
//angles
	ag->nangs = nang;
	ag->angs = (struct atomangle*)_mol_malloc(nang * sizeof(struct atomangle));
	for (i = 0; i < nang; i++) {
		i1 = angs[3 * i] - 1;
		i2 = angs[3 * i + 1] - 1;
		i3 = angs[3 * i + 2] - 1;
		ag->angs[i].a0 = &(ag->atoms[i1]);
		ag->angs[i].a1 = &(ag->atoms[i2]);
		ag->angs[i].a2 = &(ag->atoms[i3]);
		ag->atoms[i1].nangs++;
		ag->atoms[i2].nangs++;
		ag->atoms[i3].nangs++;
		ag->angs[i].th0 = fang[i];
		ag->angs[i].k = kang[i];
	}
//impropers
	ag->nimps = nimp;
	ag->imps = (struct atomimproper*)_mol_malloc(nimp * sizeof(struct atomimproper));
	for (i = 0; i < nimp; i++) {
		i1 = imps[4 * i] - 1;

		i2 = imps[4 * i + 1] - 1;
		i3 = imps[4 * i + 2] - 1;
		i4 = imps[4 * i + 3] - 1;
		ag->imps[i].a0 = &(ag->atoms[i1]);
		ag->imps[i].a1 = &(ag->atoms[i2]);
		ag->imps[i].a2 = &(ag->atoms[i3]);
		ag->imps[i].a3 = &(ag->atoms[i4]);
		ag->atoms[i1].nimps++;
		ag->atoms[i2].nimps++;
		ag->atoms[i3].nimps++;
		ag->atoms[i4].nimps++;
		ag->imps[i].k = kimp[i];
		ag->imps[i].psi0 = fimp[i];
	}
// hbond donors
	for (i = 0; i < ndon; i++) {
		int donor_i = dons[2 * i] - 1;
		int hydro_i = dons[2 * i + 1] - 1;

		mol_atom *hydro = &(ag->atoms[hydro_i]);
		mol_atom *donor = &(ag->atoms[donor_i]);

		hydro->hprop |= DONATABLE_HYDROGEN;
		donor->hprop |= HBOND_DONOR;

		hydro->base = donor_i;
	}
// hbond acceptors
	for (i = 0; i < nacc; i++) {
		int acc_i = acpts[2 * i] - 1;
		int base_i = acpts[2 * i + 1] - 1;

		mol_atom *acc = &(ag->atoms[acc_i]);

		acc->hprop |= HBOND_ACCEPTOR;
		acc->base = base_i;
	}
//torsions
	ag->tors = (struct atomtorsion*)_mol_malloc(ndih * sizeof(struct atomtorsion));
	int ltors = 0;
	for (i = 0; i < ndih; i++) {
		ti = tdih[i];
		if (ti < 0) {
			printf("read_ff_charmm WARNING:\n");
			printf("no parameters for the torsion\n");
			printf("%d %d %d %d\n", dihs[4 * i], dihs[4 * i + 1],
			       dihs[4 * i + 2], dihs[4 * i + 3]);
			continue;
		}
		i1 = dihs[4 * ti] - 1;
		i2 = dihs[4 * ti + 1] - 1;
		i3 = dihs[4 * ti + 2] - 1;
		i4 = dihs[4 * ti + 3] - 1;
		ag->tors[ltors].a0 = &(ag->atoms[i1]);
		ag->tors[ltors].a1 = &(ag->atoms[i2]);
		ag->tors[ltors].a2 = &(ag->atoms[i3]);
		ag->tors[ltors].a3 = &(ag->atoms[i4]);
		ag->atoms[i1].ntors++;
		ag->atoms[i2].ntors++;
		ag->atoms[i3].ntors++;
		ag->atoms[i4].ntors++;
		ag->tors[ltors].k = kdih[i];
		ag->tors[ltors].d = fdih[i];
		ag->tors[ltors].n = pdih[i];
		ltors++;
	}
	ag->ntors = ltors;
//allocate atom arrays of pointers to parameters
	for (i = 0; i < nat; i++) {
		i1 = ag->atoms[i].nbonds;
		ag->atoms[i].bonds = (struct atombond **)
		    _mol_malloc(i1 * sizeof(struct atombond *));
		ag->atoms[i].nbonds = 0;
		i1 = ag->atoms[i].nangs;
		ag->atoms[i].angs = (struct atomangle **)
		    _mol_malloc(i1 * sizeof(struct atomangle *));
		ag->atoms[i].nangs = 0;
		i1 = ag->atoms[i].ntors;
		ag->atoms[i].tors = (struct atomtorsion **)
		    _mol_malloc(i1 * sizeof(struct atomtorsion *));
		ag->atoms[i].ntors = 0;
		i1 = ag->atoms[i].nimps;
		ag->atoms[i].imps = (struct atomimproper **)
		    _mol_malloc(i1 * sizeof(struct atomimproper *));
		ag->atoms[i].nimps = 0;
	}
	struct atom *atm;
//fill bonds
	for (i = 0; i < ag->nbonds; i++) {
		atm = ag->bonds[i].a0;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
		atm = ag->bonds[i].a1;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
	}
//fill angles
	for (i = 0; i < ag->nangs; i++) {
		atm = ag->angs[i].a0;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a1;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a2;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
	}
//fill torsions
	for (i = 0; i < ag->ntors; i++) {
		atm = ag->tors[i].a0;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a1;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a2;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a3;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
	}
//fill impropers
	for (i = 0; i < ag->nimps; i++) {
		atm = ag->imps[i].a0;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a1;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a2;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a3;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
	}
//atom indices in the group
	fill_ingrp(ag);

	ag->is_psf_read = true;

	free(ires);
	free(dres);
	free(atnam);
	free(atoms);
	free(bonds);
	free(angs);
	free(dihs);
	free(imps);
	free(dons);
	free(acpts);
	free(kbond);
	free(lbond);
	free(kang);
	free(fang);
	free(kdih);
	free(fdih);
	free(pdih);
	free(tdih);
	free(kimp);
	free(fimp);
	free(cha);
	free(epa);
	free(siga);
	free(acevolumes);

	// copy vals from deprecated to new data structures
	int atomsi;
	for (atomsi = 0; atomsi < ag->natoms; atomsi++) {
		_mol_atom_create_bond_indices(&ag->atoms[atomsi],
					      ag->atoms[atomsi].nbonds);
	}
	_mol_atom_group_copy_from_deprecated(ag);
}

char *get_atnam_rtf(char *rtffile, int *atoms, int nat)
{
	char *buffer = (char*)_mol_malloc(sizeof(char) * linesize);
	char *atnam = (char*)_mol_malloc(sizeof(char) * nat * 4);
	char *p;
	int it, j, i;

	FILE *fp = myfopen(rtffile, "r");

	while (fgets(buffer, linesize, fp) != NULL) {
		if ((p = strchr(buffer, '!')))
			*p = '\0';
		if (strstr(buffer, "MASS") != NULL) {
			nthv(&p, buffer, ' ', 1);
			it = atoi(p);
			nthv(&p, buffer, ' ', 2);
			for (j = 0; j < nat; j++) {
				if (atoms[j] == it) {
					for (i = 0; i < 4; i++)
						atnam[4 * j + i] = *(p + i);
					continue;
				}
			}
		}
	}
	free(buffer);
	myfclose(fp);
	return atnam;
}

int nthv(char **pptext, char *ptext, char ch, int n)
{
	int r = 1, i;
	char *p = ptext;
	if (ptext[0] == ch)
		n++;
	for (i = 0; i < n; i++) {
		p = strchr(p, ch);
		if (p == NULL)
			return 0;
		while (*p == ch)
			p++;
	}
	*pptext = p;
	return r;
}

void read_psf(const char *psffile, int *o_nat, int **o_atoms,
	      char ***o_atom_names, int *o_nb, int **o_bonds, int *o_nang,
	      int **o_angs, int *o_ndih, int **o_dihs, int *o_nimp,
	      int **o_imps, float **o_cha, int *o_nres, int **o_ires,
	      char **o_dres, int *o_ndon, int **o_dons, int *o_nacc,
	      int **o_acpts)
{
	char *buffer = (char*)_mol_malloc(sizeof(char) * linesize);

	FILE *fp = myfopen(psffile, "r");

	int nat = 0, nb = 0, nang = 0, ndih = 0, nimp = 0, i, j, nres =
	    0, ndon = 0, nacc = 0;
	int *bonds = NULL, *atoms = NULL, *angs = NULL, *dihs = NULL, *imps =
	    NULL, *dons = NULL, *acpts = NULL;
	char **atom_names = NULL;
	int *ires = NULL;
	char *p, *p1, *p2, seg0[16] = "---------------", seg[16] =
	    "---------------";
	float *cha = NULL;
	while (fgets(buffer, linesize, fp) != NULL) {
		if (strstr(buffer, "NATOM") != NULL) {
			nat = atoi(buffer);
			atoms = (int*)_mol_malloc(sizeof(int) * nat);
			atom_names = (char**)_mol_malloc(sizeof(char *) * nat);
			cha = (float*)_mol_malloc(sizeof(float) * nat);
			ires = (int*)_mol_malloc(sizeof(int) * nat);

			for (i = 0; i < nat; i++) {
				if (fgets(buffer, linesize, fp) == NULL ||
				    nthv(&p1, buffer, ' ', 1) == 0 ||
				    nthv(&p, buffer, ' ', 2) == 0 ||
				    nthv(&p2, buffer, ' ', 4) == 0 ||
				    nthv(&p, buffer, ' ', 5) == 0 ||
				    (j = atoi(p)) == 0
				    || nthv(&p, buffer, ' ', 6) == 0) {
					printf
					    ("error: reading atoms from psf\n");
					exit(EXIT_FAILURE);
				}

				strncpy(seg, p1, 15);
				if (strcmp(seg, seg0)) {
					ires[nres] = i;
					strcpy(seg0, seg);
					o_dres[nres] = (char*)
					    _mol_malloc(sizeof(char) * 16);
					strcpy(o_dres[nres], seg);
					nres++;
				}

				atom_names[i] = (char*)_mol_calloc(sizeof(char), 5);
				strncpy(atom_names[i], p2, 4);
				rstrip(atom_names[i]);

				atoms[i] = j;
				cha[i] = atof(p);
			}
			continue;
		}
		if (strstr(buffer, "NBOND") != NULL) {
			nb = atoi(buffer);
			bonds = (int*)_mol_malloc(sizeof(int) * 2 * nb);
			for (i = 0; i < nb; i++) {
				fscanf(fp, "%d %d", &bonds[2 * i],
				       &bonds[2 * i + 1]);
			}
			continue;
		}
		if (strstr(buffer, "NTHETA") != NULL) {
			nang = atoi(buffer);
			angs = (int*)_mol_malloc(sizeof(int) * 3 * nang);
			for (i = 0; i < nang; i++) {
				fscanf(fp, "%d %d %d", &angs[3 * i],
				       &angs[3 * i + 1], &angs[3 * i + 2]);
			}
			continue;
		}
		if (strstr(buffer, "NPHI") != NULL) {
			ndih = atoi(buffer);
			dihs = (int*)_mol_malloc(sizeof(int) * 4 * ndih);
			for (i = 0; i < ndih; i++) {
				fscanf(fp, "%d %d %d %d", &dihs[4 * i],
				       &dihs[4 * i + 1], &dihs[4 * i + 2],
				       &dihs[4 * i + 3]);
			}
			continue;
		}
		if (strstr(buffer, "NIMPHI") != NULL) {
			nimp = atoi(buffer);
			imps = (int*)_mol_malloc(sizeof(int) * 4 * nimp);
			for (i = 0; i < nimp; i++) {
				fscanf(fp, "%d %d %d %d", &imps[4 * i],
				       &imps[4 * i + 1], &imps[4 * i + 2],
				       &imps[4 * i + 3]);
			}
			continue;
		}
		if (strstr(buffer, "NDON") != NULL) {
			ndon = atoi(buffer);
			dons = (int*)_mol_malloc(sizeof(int) * 2 * ndon);
			for (i = 0; i < ndon; i++) {
				fscanf(fp, "%d %d", &dons[2 * i],
				       &dons[2 * i + 1]);
			}
			continue;
		}
		if (strstr(buffer, "NACC") != NULL) {
			nacc = atoi(buffer);
			acpts = (int*)_mol_malloc(sizeof(int) * 2 * nacc);
			for (i = 0; i < nacc; i++) {
				fscanf(fp, "%d %d", &acpts[2 * i],
				       &acpts[2 * i + 1]);
			}
			continue;
		}
	}
	myfclose(fp);
	free(buffer);
	*o_nat = nat;
	*o_atoms = atoms;
	*o_atom_names = atom_names;
	*o_nb = nb;
	*o_bonds = bonds;
	*o_nang = nang;
	*o_angs = angs;
	*o_ndih = ndih;
	*o_dihs = dihs;
	*o_nimp = nimp;
	*o_imps = imps;
	*o_cha = cha;
	*o_nres = nres;
	*o_ires = ires;
	*o_ndon = ndon;
	*o_dons = dons;
	*o_nacc = nacc;
	*o_acpts = acpts;
}

void read_par(char *prmfile, int nat, char *atnam,
	      int nb, int *bonds, float **o_kbond, float **o_lbond,
	      int nang, int *angles, float **o_kangle, float **o_fangle,
	      int *ndih, int *dihs, float **o_kdih, float **o_fdih,
	      int **o_pdih, int **o_tdih, int nimp, int *imps, float **o_kimp,
	      float **o_fimp, float **o_epa, float **o_siga,
	      float **o_acevolumes)
{
	int i, ndih0 = *ndih;
	char *buffer = (char*)_mol_malloc(sizeof(char) * linesize);
	char *p;

	float *kbond = (float*)_mol_malloc(sizeof(float) * nb);
	for (i = 0; i < nb; i++)
		kbond[i] = -1.0;

	float *lbond = (float*)_mol_malloc(sizeof(float) * nb);
	for (i = 0; i < nb; i++)
		lbond[i] = -1.0;

	float *kangle = (float*)_mol_malloc(sizeof(float) * nang);
	for (i = 0; i < nang; i++)
		kangle[i] = -1.0;

	float *fangle = (float*)_mol_malloc(sizeof(float) * nang);
	for (i = 0; i < nang; i++)
		fangle[i] = -1.0;

	float *kdih = (float*)_mol_malloc(sizeof(float) * 3 * ndih0);
	for (i = 0; i < 3 * ndih0; i++)
		kdih[i] = -1.0;

	float *fdih = (float*)_mol_malloc(sizeof(float) * 3 * ndih0);
	for (i = 0; i < 3 * ndih0; i++)
		fdih[i] = -1.0;

	int *pdih = (int*)_mol_malloc(sizeof(int) * 3 * ndih0);
	for (i = 0; i < 3 * ndih0; i++)
		pdih[i] = -1;

	int *wdih = (int*)_mol_malloc(sizeof(int) * ndih0);
	for (i = 0; i < ndih0; i++)
		wdih[i] = -1;

	int *tdih = (int*)_mol_malloc(sizeof(int) * 3 * ndih0);
	for (i = 0; i < 3 * ndih0; i++)
		tdih[i] = -1;

	float *kimp = (float*)_mol_malloc(sizeof(float) * nimp);
	for (i = 0; i < nimp; i++)
		kimp[i] = -1.0;

	float *fimp = (float*)_mol_malloc(sizeof(float) * nimp);
	for (i = 0; i < nimp; i++)
		fimp[i] = -1.0;

	float *epa = (float*)_mol_malloc(2 * sizeof(float) * nat);
	for (i = 0; i < 2 * nat; i++)
		epa[i] = 1.0;

	float *siga = (float*)_mol_malloc(2 * sizeof(float) * nat);
	for (i = 0; i < 2 * nat; i++)
		siga[i] = -1.0;

	float *acevolumes = (float*)_mol_malloc(sizeof(float) * nat);
	for (i = 0; i < nat; i++)
		acevolumes[i] = -1.0;

	FILE *fp = myfopen(prmfile, "r");
	double efac = 0.0;

	int read_mode = 0;
	while (fgets(buffer, linesize, fp) != NULL) {
		if ((p = strstr(buffer, "e14fac"))
		    || (p = strstr(buffer, "E14FAC")))
			efac = atof(p + 6);
		if ((p = strchr(buffer, '!')))
			*p = '\0';
		if ((strstr(buffer, "NBON") != NULL)
		    || (strstr(buffer, "NONB") != NULL)) {
			read_mode = 5;
			continue;
		} else if (strstr(buffer, "HBOND") != NULL) {
			read_mode = 7;
			continue;
		} else if (strstr(buffer, "BOND") != NULL) {
			read_mode = 1;
			continue;
		} else if (strstr(buffer, "ANGL") != NULL
			   || strstr(buffer, "THETA") != NULL) {
			read_mode = 2;
			continue;
		} else if (strstr(buffer, "IMPR") != NULL
			   || strstr(buffer, "IMPHI") != NULL) {
			read_mode = 4;
			continue;
		} else if (strstr(buffer, "DIHE") != NULL
			   || strstr(buffer, "PHI") != NULL) {
			read_mode = 3;
			continue;
		} else if (strstr(buffer, "NBFIX") != NULL) {
			read_mode = 6;
			continue;
		} else if (strstr(buffer, "VOLUME") != NULL) {
			read_mode = 8;
			continue;
		}
		if (read_mode == 1)
			chkbond(buffer, nb, bonds, atnam, kbond, lbond);
		if (read_mode == 2)
			chkangle(buffer, nang, angles, atnam, kangle, fangle);
		if (read_mode == 3)
			chkdih(buffer, ndih, ndih0, dihs, atnam, kdih, fdih,
			       pdih, tdih, wdih);
		if (read_mode == 4)
			chkimp(buffer, nimp, imps, atnam, kimp, fimp);
		if (read_mode == 5)
			chkvdw(buffer, nat, atnam, epa, siga);
		if (read_mode == 8)
			chkacevolumes(buffer, nat, atnam, acevolumes);
	}

	myfclose(fp);
	free(buffer);
	free(wdih);
	*o_kbond = kbond;
	*o_lbond = lbond;
	*o_kangle = kangle;
	*o_fangle = fangle;
	*o_kdih = kdih;
	*o_fdih = fdih;
	*o_pdih = pdih;
	*o_kimp = kimp;
	*o_fimp = fimp;
	*o_epa = epa;
	*o_siga = siga;
	*o_tdih = tdih;
	*o_acevolumes = acevolumes;
}

void chkbond(char *buffer, int nb, int *bonds, char *atnam, float *kbond,
	     float *lbond)
{
	int i;
	char *p, match[10] = "         ", tag[20];
	float l, k;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 1) == 0)
		return;
	strncpy(&match[5], p, 4);
	if (nthv(&p, buffer, ' ', 2) == 0)
		return;
	if ((k = atof(p)) == 0)
		return;
	if (nthv(&p, buffer, ' ', 3) == 0)
		return;
	if ((l = atof(p)) == 0)
		return;
	for (i = 0; i < nb; i++) {
		bondtag(i, bonds, atnam, tag);
		if (strstr(tag, match) != NULL) {
			kbond[i] = k;
			lbond[i] = l;
		}
	}
}

void bondtag(int i, int *bonds, char *atnam, char *p)
{
	int j;
	for (j = 0; j < 4; j++) {
		p[j] = atnam[4 * (bonds[2 * i] - 1) + j];
		p[j + 5] = atnam[4 * (bonds[2 * i + 1] - 1) + j];
		p[j + 10] = p[j + 5];
		p[j + 15] = p[j];
	}
	p[19] = '\0';
	p[4] = ' ';
	p[9] = '+';
	p[14] = ' ';
}

void chkangle(char *buffer, int nang, int *angles,
	      char *atnam, float *kangle, float *fangle)
{
	int i;
	char *p, match[15] = "              ", tag[30];
	float l, k;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 1) == 0)
		return;
	strncpy(&match[5], p, 4);
	if (nthv(&p, buffer, ' ', 2) == 0)
		return;
	strncpy(&match[10], p, 4);
	if (nthv(&p, buffer, ' ', 3) == 0)
		return;
	if ((k = atof(p)) == 0)
		return;
	if (nthv(&p, buffer, ' ', 4) == 0)
		return;
	if ((l = atof(p)) == 0)
		return;
	for (i = 0; i < nang; i++) {
		angletag(i, angles, atnam, tag);
		if (strstr(tag, match) != NULL) {
			kangle[i] = k;
			fangle[i] = l;
		}
	}
}

void angletag(int i, int *angles, char *atnam, char *p)
{
	int j;
	for (j = 0; j < 4; j++) {
		p[j] = atnam[4 * (angles[3 * i] - 1) + j];
		p[j + 5] = atnam[4 * (angles[3 * i + 1] - 1) + j];
		p[j + 10] = atnam[4 * (angles[3 * i + 2] - 1) + j];
		p[j + 15] = p[j + 10];
		p[j + 20] = p[j + 5];
		p[j + 25] = p[j];
	}
	p[29] = '\0';
	p[4] = ' ';
	p[9] = ' ';
	p[14] = '+';
	p[19] = ' ';
	p[24] = ' ';
}

/* modes are following
0) A - B - C - D
1) A - X - X - D
2) X - B - C - D
3) X - B - C - X
4) X - X - C - D
*/
void dihtag(int i, int *dihs, char *atnam, int mode, char *p)
{
	int j;
	for (j = 0; j < 4; j++) {
		p[j] = atnam[4 * (dihs[4 * i] - 1) + j];
		p[j + 5] = atnam[4 * (dihs[4 * i + 1] - 1) + j];
		p[j + 10] = atnam[4 * (dihs[4 * i + 2] - 1) + j];
		p[j + 15] = atnam[4 * (dihs[4 * i + 3] - 1) + j];
	}
	if (mode == 1) {
		putXtoP(p, 5);
		putXtoP(p, 10);
	} else if (mode == 2)
		putXtoP(p, 0);
	else if (mode == 3) {
		putXtoP(p, 0);
		putXtoP(p, 15);
	} else if (mode == 4) {
		putXtoP(p, 0);
		putXtoP(p, 5);
	}
	for (j = 0; j < 4; j++) {
		p[j + 20] = p[j + 15];
		p[j + 25] = p[j + 10];
		p[j + 30] = p[j + 5];
		p[j + 35] = p[j];
	}
	p[39] = '\0';
	p[4] = ' ';
	p[9] = ' ';
	p[14] = ' ';
	p[19] = '+';
	p[24] = ' ';
	p[29] = ' ';
	p[34] = ' ';
}

void chkdih(char *buffer, int *ndih, int ndih0, int *dihs,
	    char *atnam, float *kdih, float *fdih, int *pdih, int *tdih,
	    int *wdih)
{
	int i, ndih03 = 3 * ndih0;
	char *p, match[20] = "                   ";
	char tag0[40], tag3[40];
	float l, k;
	int n;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 1) == 0)
		return;
	strncpy(&match[5], p, 4);
	if (nthv(&p, buffer, ' ', 2) == 0)
		return;
	strncpy(&match[10], p, 4);
	if (nthv(&p, buffer, ' ', 3) == 0)
		return;
	strncpy(&match[15], p, 4);
	if (nthv(&p, buffer, ' ', 4) == 0)
		return;
	k = atof(p);
	if (nthv(&p, buffer, ' ', 5) == 0)
		return;
	n = atoi(p);
	if (nthv(&p, buffer, ' ', 6) == 0)
		return;
	l = atof(p);
	for (i = 0; i < ndih0; i++) {
		dihtag(i, dihs, atnam, 0, tag0);
		dihtag(i, dihs, atnam, 3, tag3);
		if (strstr(tag0, match) != NULL) {
			if (wdih[i] < 1) {
				wdih[i] = 1;
				kdih[i] = k;
				fdih[i] = l;
				pdih[i] = n;
				tdih[i] = i;
				continue;
			} else {
				if ((*ndih) < ndih03) {
					kdih[(*ndih)] = k;
					fdih[(*ndih)] = l;
					pdih[(*ndih)] = n;
					tdih[(*ndih)] = i;
					(*ndih)++;
					continue;
				} else {
					printf
					    ("error: need more memory for dihedrals\n");
					exit(1);
				}
			}
		} else if (strstr(tag3, match) != NULL) {
			if (wdih[i] < 0) {
				wdih[i] = 0;
				kdih[i] = k;
				fdih[i] = l;
				pdih[i] = n;
				tdih[i] = i;
				continue;
			} else if (wdih[i] == 0) {
				if ((*ndih) < ndih03) {
					kdih[(*ndih)] = k;
					fdih[(*ndih)] = l;
					pdih[(*ndih)] = n;
					tdih[(*ndih)] = i;
					(*ndih)++;
					continue;
				}
			}
		}
	}
}

void putXtoP(char *p, int n)
{
	int i;
	p[n] = 'X';
	for (i = 1; i < 4; i++)
		p[n + i] = ' ';
}

void chkimp(char *buffer, int nimp, int *imps,
	    char *atnam, float *kimp, float *fimp)
{
	int i;
	char *p, match[20] = "                   ";
	char tag0[40], tag1[40], tag2[40], tag3[40], tag4[40];
	float l, k;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 1) == 0)
		return;
	strncpy(&match[5], p, 4);
	if (nthv(&p, buffer, ' ', 2) == 0)
		return;
	strncpy(&match[10], p, 4);
	if (nthv(&p, buffer, ' ', 3) == 0)
		return;
	strncpy(&match[15], p, 4);
	if (nthv(&p, buffer, ' ', 4) == 0)
		return;
	if ((k = atof(p)) == 0)
		return;
	if (nthv(&p, buffer, ' ', 6) == 0)
		return;
	l = atof(p);
	for (i = 0; i < nimp; i++) {
		dihtag(i, imps, atnam, 0, tag0);
		dihtag(i, imps, atnam, 1, tag1);
		dihtag(i, imps, atnam, 2, tag2);
		dihtag(i, imps, atnam, 3, tag3);
		dihtag(i, imps, atnam, 4, tag4);
		if (strstr(tag0, match) != NULL || strstr(tag1, match) != NULL
		    || strstr(tag2, match) != NULL
		    || strstr(tag3, match) != NULL
		    || strstr(tag4, match) != NULL) {
			kimp[i] = k;
			fimp[i] = l;
			continue;
		}
	}
}

void chkvdw(char *buffer, int nat, char *atnam, float *epa, float *siga)
{
	int i, j, nm;
	char *p, *p1, match[5] = "    ";
	char tag[5];
	float e, s, e03, s03;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 2) == 0)
		return;
	e = atof(p);
	if (nthv(&p, buffer, ' ', 3) == 0)
		return;
	s = atof(p);
	e03 = 1.0;
	s03 = -1.0;
	if (nthv(&p, buffer, ' ', 5) != 0) {
		if (p[0] == '-') {
			if (nthv(&p1, buffer, ' ', 6) != 0) {
				e03 = atof(p);
				s03 = atof(p1);
			}
		}
	}
	nm = 4;
	for (i = 0; i < 4; i++) {
		if (match[i] == '*') {
			nm = i;
			break;
		}
	}
	if (nm < 4) {
		for (i = 0; i < nat; i++) {
			for (j = 0; j < nm; j++)
				tag[j] = atnam[4 * i + j];
			if (strncmp(tag, match, nm) == 0) {

				if (epa[2 * i] > 0) {
					epa[2 * i] = e;
					siga[2 * i] = s;
					epa[2 * i + 1] = e03;
					siga[2 * i + 1] = s03;
				}
			}
		}
	} else {
		for (i = 0; i < nat; i++) {
			for (j = 0; j < 4; j++)
				tag[j] = atnam[4 * i + j];
			if (strncmp(tag, match, 4) == 0) {
				epa[2 * i] = e;
				siga[2 * i] = s;
				epa[2 * i + 1] = e03;
				siga[2 * i + 1] = s03;
			}
		}
	}
}

void chkacevolumes(char *buffer, int nat, char *atnam, float *acevolumes)
{
	int i, j, nm;
	char *p, match[5] = "    ";
	char tag[5];
	float volume;
	if (nthv(&p, buffer, ' ', 0) == 0)
		return;
	strncpy(match, p, 4);
	if (nthv(&p, buffer, ' ', 1) == 0)
		return;
	volume = atof(p);
	nm = 4;
	for (i = 0; i < 4; i++) {
		if (match[i] == '*') {
			nm = i;
			break;
		}
	}
	if (nm < 4) {
		for (i = 0; i < nat; i++) {
			for (j = 0; j < nm; j++)
				tag[j] = atnam[4 * i + j];
			if (strncmp(tag, match, nm) == 0) {
				acevolumes[i] = volume;
			}
		}
	} else {
		for (i = 0; i < nat; i++) {
			for (j = 0; j < 4; j++)
				tag[j] = atnam[4 * i + j];
			if (strncmp(tag, match, 4) == 0) {
				acevolumes[i] = volume;
			}
		}
	}
}

void fixed_init(struct atomgrp *ag)
{
	ag->nimpact = 0;
	ag->ntoract = 0;
	ag->nangact = 0;
	ag->nbact = 0;
	ag->nactives = 0;
}

void fixed_update(struct atomgrp *ag, int nlist, int *list)
{

	int i, m, n;
	struct atom *a0, *a1, *a2, *a3;
	if (ag->nbact > 0) {
		free(ag->bact);
		ag->nbact = 0;
	}
	if (ag->nangact > 0) {
		free(ag->angact);
		ag->nangact = 0;
	}

	if (ag->ntoract > 0) {
		free(ag->toract);
		ag->ntoract = 0;
	}
	if (ag->nimpact > 0) {
		free(ag->impact);
		ag->nimpact = 0;
	}
	if (ag->nactives > 0) {
		free(ag->activelist);
		ag->nactives = 0;
	}
// atoms
	for (i = 0; i < ag->natoms; i++)
		ag->atoms[i].fixed = 0;
	for (i = 0; i < nlist; i++)
		ag->atoms[list[i]].fixed = 1;
	int ci;
	ag->activelist = (int *)_mol_malloc((ag->natoms - nlist) * sizeof(int));
	ci = 0;
	for (i = 0; i < ag->natoms; i++) {
		if (ag->atoms[i].fixed == 0)
			ag->activelist[ci++] = i;
	}
	ag->nactives = ci;
// bonds
	n = ag->nbonds;
	m = 0;
	ag->bact =(struct atombond**)_mol_malloc(n * sizeof(struct atombond *));
	for (i = 0; i < n; i++) {
		a0 = ag->bonds[i].a0;
		a1 = ag->bonds[i].a1;
		if (a0->fixed == 1 && a1->fixed == 1)
			continue;
		ag->bact[m++] = &(ag->bonds[i]);
	}
	ag->nbact = m;
	if (m > 0)
		ag->bact = (struct atombond **)_mol_realloc(ag->bact,
						     m *
						     sizeof(struct atombond *));
	else
		free(ag->bact);
// angles
	n = ag->nangs;
	m = 0;
	ag->angact = (struct atomangle **)_mol_malloc(n * sizeof(struct atomangle *));
	for (i = 0; i < n; i++) {
		a0 = ag->angs[i].a0;
		a1 = ag->angs[i].a1;
		a2 = ag->angs[i].a2;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1)
			continue;
		ag->angact[m++] = &(ag->angs[i]);
	}
	ag->nangact = m;
	if (m > 0)
		ag->angact =
		    (struct atomangle **)_mol_realloc(ag->angact,
						      m *
						      sizeof(struct atomangle
							     *));
	else
		free(ag->angact);
// dihedrals
	n = ag->ntors;
	m = 0;
	ag->toract = (struct atomtorsion **)_mol_malloc(n * sizeof(struct atomtorsion *));
	for (i = 0; i < n; i++) {
		a0 = ag->tors[i].a0;
		a1 = ag->tors[i].a1;
		a2 = ag->tors[i].a2;
		a3 = ag->tors[i].a3;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1
		    && a3->fixed == 1)
			continue;
		ag->toract[m++] = &(ag->tors[i]);
	}
	ag->ntoract = m;
	if (m > 0)
		ag->toract =
		    (struct atomtorsion **)_mol_realloc(ag->toract,
							m *
							sizeof(struct
							       atomtorsion *));
	else
		free(ag->toract);
// impropers
	n = ag->nimps;
	m = 0;
	ag->impact = (struct atomimproper **)_mol_malloc(n * sizeof(struct atomimproper *));
	for (i = 0; i < n; i++) {
		a0 = ag->imps[i].a0;
		a1 = ag->imps[i].a1;
		a2 = ag->imps[i].a2;
		a3 = ag->imps[i].a3;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1
		    && a3->fixed == 1)
			continue;
		ag->impact[m++] = &(ag->imps[i]);
	}
	ag->nimpact = m;
	if (m > 0)
		ag->impact =
		    (struct atomimproper **)_mol_realloc(ag->impact,
							 m *
							 sizeof(struct
								atomimproper
								*));
	else
		free(ag->impact);
}

/*Light version to read only bond info*/
void read_ff_charmm_light(const char *psffile, struct atomgrp *ag)
{
	int i;
	int nat, nb, nang, ndih0, nimp, ndon, nacc;
	int *atoms, *bonds, *angs, *dihs, *imps, *dons, *acpts;
	char **atom_names;	//CA, CB, etc.
	float *cha;
	int nres, *ires;
	char **dres = (char **)_mol_malloc(ag->natoms * sizeof(char *));
//        printf("PSF INn\n");

	read_psf(psffile, &nat, &atoms, &atom_names, &nb, &bonds, &nang, &angs,
		 &ndih0, &dihs, &nimp, &imps, &cha, &nres, &ires, dres,
		 &ndon, &dons, &nacc, &acpts);
//      printf("PSF INn\n");
// residues
	ag->nres = nres;
	if (nres > 0) {
		ag->iares = (int*)_mol_malloc((nres + 1) * sizeof(int));
		ag->idres = (char **)_mol_malloc(nres * sizeof(char *));
	}
	for (i = 0; i < nres; i++) {
		ag->iares[i] = ires[i];
		ag->idres[i] = dres[i];
	}
	ag->iares[nres] = nat;

	//atnam=get_atnam_rtf(rtffile, atoms, nat);
/*	read_par(prmfile, nat, atnam, nb, bonds, &kbond, &lbond,
					nang, angs, &kang, &fang,
					&ndih, dihs, &kdih, &fdih, &pdih, &tdih,
					nimp, imps, &kimp, &fimp,
		 &epa, &siga,&acevolumes);
*/
//        printf("DIHEDRAL ndih0, ndih %d %d\n", ndih0, ndih);
	if ((ag->natoms) != nat) {
		fprintf(stderr,
			"different number of atoms after reading protein structure file (PSF)");
		exit(EXIT_FAILURE);
	}
// vdw and charges, types and acevolumes
//printf("Atoms INn\n");
	for (i = 0; i < nat; i++) {
		ag->atoms[i].atom_ftypen = atoms[i];
		ag->atoms[i].name = atom_names[i];
		ag->atoms[i].nbonds = 0;
		ag->atoms[i].nangs = 0;
		ag->atoms[i].ntors = 0;
		ag->atoms[i].nimps = 0;
	}
// bonds
//printf("Bonds  INn\n");
	ag->nbonds = nb;
	ag->bonds = (struct atombond*)_mol_malloc(nb * sizeof(struct atombond));
	for (i = 0; i < nb; i++) {
		int i1 = bonds[2 * i] - 1;
		int i2 = bonds[2 * i + 1] - 1;
		ag->bonds[i].a0 = &(ag->atoms[i1]);
		ag->bonds[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i1].nbonds)++;
		(ag->atoms[i2].nbonds)++;

	}
//printf("Angles  INn\n");
//angles
	ag->nangs = nang;
	ag->angs = (struct atomangle*)_mol_malloc(nang * sizeof(struct atomangle));
	for (i = 0; i < nang; i++) {
		int i1 = angs[3 * i] - 1;
		int i2 = angs[3 * i + 1] - 1;
		int i3 = angs[3 * i + 2] - 1;
		ag->angs[i].a0 = &(ag->atoms[i1]);
		ag->angs[i].a1 = &(ag->atoms[i2]);
		ag->angs[i].a2 = &(ag->atoms[i3]);
		ag->atoms[i1].nangs++;
		ag->atoms[i2].nangs++;
		ag->atoms[i3].nangs++;
	}
//impropers
	ag->nimps = nimp;
	ag->imps = (struct atomimproper*)_mol_malloc(nimp * sizeof(struct atomimproper));
	for (i = 0; i < nimp; i++) {
		int i1 = imps[4 * i] - 1;

		int i2 = imps[4 * i + 1] - 1;
		int i3 = imps[4 * i + 2] - 1;
		int i4 = imps[4 * i + 3] - 1;
		ag->imps[i].a0 = &(ag->atoms[i1]);
		ag->imps[i].a1 = &(ag->atoms[i2]);
		ag->imps[i].a2 = &(ag->atoms[i3]);
		ag->imps[i].a3 = &(ag->atoms[i4]);
		ag->atoms[i1].nimps++;
		ag->atoms[i2].nimps++;
		ag->atoms[i3].nimps++;
		ag->atoms[i4].nimps++;
	}
// hbond donors
	for (i = 0; i < ndon; i++) {
		int donor_i = dons[2 * i] - 1;
		int hydro_i = dons[2 * i + 1] - 1;

		mol_atom *hydro = &(ag->atoms[hydro_i]);
		mol_atom *donor = &(ag->atoms[donor_i]);

		hydro->hprop |= DONATABLE_HYDROGEN;
		donor->hprop |= HBOND_DONOR;

		hydro->base = donor_i;
	}
// hbond acceptors
	for (i = 0; i < nacc; i++) {
		int acc_i = acpts[2 * i] - 1;
		int base_i = acpts[2 * i + 1] - 1;

		mol_atom *acc = &(ag->atoms[acc_i]);

		acc->hprop |= HBOND_ACCEPTOR;
		acc->base = base_i;
	}
//torsions
	ag->tors = (struct atomtorsion*)_mol_malloc(ndih0 * sizeof(struct atomtorsion));
	int ltors = 0;
	for (i = 0; i < ndih0; i++) {
		int i1 = dihs[4 * i] - 1;
		int i2 = dihs[4 * i + 1] - 1;
		int i3 = dihs[4 * i + 2] - 1;
		int i4 = dihs[4 * i + 3] - 1;
		ag->tors[ltors].a0 = &(ag->atoms[i1]);
		ag->tors[ltors].a1 = &(ag->atoms[i2]);
		ag->tors[ltors].a2 = &(ag->atoms[i3]);
		ag->tors[ltors].a3 = &(ag->atoms[i4]);
		ag->atoms[i1].ntors++;
		ag->atoms[i2].ntors++;
		ag->atoms[i3].ntors++;
		ag->atoms[i4].ntors++;
		ltors++;
	}
	ag->ntors = ltors;

//printf("All  INn\n");
//allocate atom arrays of pointers to parameters
	for (i = 0; i < nat; i++) {
		int i1 = ag->atoms[i].nbonds;
		ag->atoms[i].bonds = (struct atombond **)
		    _mol_malloc(i1 * sizeof(struct atombond *));
		ag->atoms[i].nbonds = 0;
		i1 = ag->atoms[i].nangs;
		ag->atoms[i].angs = (struct atomangle **)
		    _mol_malloc(i1 * sizeof(struct atomangle *));
		ag->atoms[i].nangs = 0;
		i1 = ag->atoms[i].ntors;
		ag->atoms[i].tors = (struct atomtorsion **)
		    _mol_malloc(i1 * sizeof(struct atomtorsion *));
		ag->atoms[i].ntors = 0;
		i1 = ag->atoms[i].nimps;
		ag->atoms[i].imps = (struct atomimproper **)
		    _mol_malloc(i1 * sizeof(struct atomimproper *));
		ag->atoms[i].nimps = 0;
	}
	struct atom *atm;
//fill bonds
	for (i = 0; i < ag->nbonds; i++) {
		atm = ag->bonds[i].a0;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
		atm = ag->bonds[i].a1;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
	}
//fill angles
	for (i = 0; i < ag->nangs; i++) {
		atm = ag->angs[i].a0;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a1;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a2;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
	}
//fill torsions
	for (i = 0; i < ag->ntors; i++) {
		atm = ag->tors[i].a0;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a1;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a2;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a3;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
	}
//fill impropers
	for (i = 0; i < ag->nimps; i++) {
		atm = ag->imps[i].a0;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a1;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a2;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a3;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
	}

//atom indices in the group
	fill_ingrp(ag);

	ag->is_psf_read = true;

	// copy vals from deprecated to new data structures
	int atomsi;
	for (atomsi = 0; atomsi < ag->natoms; atomsi++) {
		_mol_atom_create_bond_indices(&ag->atoms[atomsi],
					      ag->atoms[atomsi].nbonds);
	}
	_mol_atom_group_copy_from_deprecated(ag);

	free(ires);
	free(dres);
	free(atoms);
	free(bonds);
	free(angs);
	free(dihs);
	free(imps);
	free(cha);
}

};
