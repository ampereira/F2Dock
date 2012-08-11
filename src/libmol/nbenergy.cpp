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
#include <math.h>
#include <errno.h>

#include <libmol/libmol.h>

namespace LibMol{

#define MAXPAIR 36
#define CCELEC 332.0716


void test_nbgrads(struct atomgrp *ag, double d, struct nblist *nblst,
                  int n03, int* list03)
{
int n=ag->natoms, i;
        double en, en1, t, f=0.5;
	double *fs=(double*)_mol_malloc(3*n*sizeof(double));
//en0
        en=0;
        vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
//        vdweng03(f, ag, &en, n03, list03);
//        eleng03(f, ag, 1.0, &en, n03, list03);
//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
//        vdweng(ag, &en, nblst);
//        eleng(ag, 1.0, &en, nblst);

        for(i=0; i<n; i++)
        {
//x
                en1=0;
                t=ag->atoms[i].X;
                ag->atoms[i].X=d+t;
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
//                vdweng03(f, ag, &en1, n03, list03);
//                eleng03(f, ag, 1.0, &en1, n03, list03);
//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
//                vdweng(ag, &en1, nblst);
//                eleng(ag, 1.0, &en1, nblst);
                ag->atoms[i].X=t;
                fs[3*i]=(en-en1)/d;
//y
                en1=0;
                t=ag->atoms[i].Y;
                ag->atoms[i].Y=d+t;
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
//                vdweng03(f, ag, &en1, n03, list03);
//                eleng03(f, ag, 1.0, &en1, n03, list03);
//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
//                vdweng(ag, &en1, nblst);
//                eleng(ag, 1.0, &en1, nblst);
                ag->atoms[i].Y=t;
                fs[3*i+1]=(en-en1)/d;
//z
                en1=0;
                t=ag->atoms[i].Z;
                ag->atoms[i].Z=d+t;
                vdwengs03(f, nblst->nbcof, ag, &en1, n03, list03);
//                vdweng03(f, ag, &en1, n03, list03);
//                eleng(f, ag, 1.0, &en1, n03, list03);
//                elengs03(f, nblst->nbcof, ag, 1.0, &en1, n03, list03);
//                vdweng(ag, &en1, nblst);
//                eleng(ag, 1.0, &en1, nblst);
                ag->atoms[i].Z=t;
                fs[3*i+2]=(en-en1)/d;
        }
        en=0;
        zero_grads(ag);
        vdwengs03(f, nblst->nbcof, ag, &en, n03, list03);
//        vdweng03(f, ag, &en, n03, list03);
//        eleng(f, ag, 1.0, &en, n03, list03);
//        elengs03(f, nblst->nbcof, ag, 1.0, &en, n03, list03);
//        vdweng(ag, &en, nblst);
//        eleng(ag, 1.0, &en, nblst);
        for(i=0; i<n; i++)
        {
                printf("%d calculated: %lf %lf %lf\n",
                i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
                printf("%d numerical : %lf %lf %lf\n",
                i,fs[3*i],fs[3*i+1],fs[3*i+2]);
        }
        free(fs);
}

void vdweng03(double f, struct atomgrp *ag, double* ven, int n03, int* list03)
{
        int i, i1, i2;
        struct atom *a1, *a2;
        double ei, ej, eij, ri, rj, rij;
        double dx, dy, dz, d2;
        double Rd6, Rd12, dven, g;
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           ei=f*(a1->eps03);
           ri=a1->rminh03;

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           ej=a2->eps03;
           rj=a2->rminh03;

           eij=ei*ej;
           rij=ri+rj;
           rij*=rij;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           Rd6=rij/d2;
           Rd6=Rd6*Rd6*Rd6;
           Rd12=Rd6*Rd6;
           (*ven)+=eij*(Rd12-2*Rd6);
           dven=-eij*12*(-Rd12+Rd6)/d2;
           g=dven*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=dven*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=dven*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
        }
}

void vdwengs03(double f, double rc, struct atomgrp *ag, double* ven,
               int n03, int* list03)
{
        int i, i1, i2;
        double ei, ri, ej, rj, dx, dy, dz, eij;
        double d2, Rd12, Rd6, Rr12, Rr6, dr6;
        double rij, dven, g;
        struct atom *a1, *a2;
        double rc2=rc*rc;

        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           i2=list03[2*i+1];
           if(i1 == i2)continue;
           a1=&(ag->atoms[i1]);
           ei=f*(a1->eps03);
           ri=a1->rminh03;

           a2=&(ag->atoms[i2]);
           ej=a2->eps03;
           rj=a2->rminh03;

           eij=ei*ej;
           rij=ri+rj;
           rij*=rij;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;

           Rd6=rij/d2;
           Rd6=Rd6*Rd6*Rd6;
           Rd12=Rd6*Rd6;
           Rr6=rij/rc2;
           Rr6=Rr6*Rr6*Rr6;
           Rr12=Rr6*Rr6;
           dr6=d2/rc2;
           dr6=dr6*dr6*dr6;

           (*ven)+=eij*(Rd12-2*Rd6+Rr6*(4.0-2*dr6)+Rr12*(2*dr6-3.0));
           dven=-eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;
           g=dven*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=dven*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=dven*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
	}
}

void vdweng(struct atomgrp *ag, double* ven, struct nblist *nblst)
{
        int i, i1, n2, *p, j, i2;
        double ei, ri, ej, rj, x1, y1, z1, dx, dy, dz, eij;
        double d2, Rd12, Rd6, Rr12, Rr6, dr6;
        double rij, dven, g;
        struct atom *a1, *a2;
        double rc=nblst->nbcof;
        double rc2=rc*rc;

        for(i=0; i<nblst->nfat; i++)
        {
           i1=nblst->ifat[i];
           a1=&(ag->atoms[i1]);
           ei=a1->eps;
           ri=a1->rminh;
           x1=ag->atoms[i1].X;
           y1=ag->atoms[i1].Y;
           z1=ag->atoms[i1].Z;
           n2=nblst->nsat[i];
           p=nblst->isat[i];
           for(j=0; j<n2; j++)
           {
              i2=p[j];
              a2=&(ag->atoms[i2]);
              rj=a2->rminh;
              ej=a2->eps;
              dx=x1-ag->atoms[i2].X;
              dy=y1-ag->atoms[i2].Y;
              dz=z1-ag->atoms[i2].Z;
              d2=dx*dx+dy*dy+dz*dz;
              if(d2<rc2)
              {
                 eij=ei*ej;
                 rij=ri+rj;
                 rij*=rij;

                 Rd6=rij/d2;
                 Rd6=Rd6*Rd6*Rd6;
                 Rd12=Rd6*Rd6;
                 Rr6=rij/rc2;
                 Rr6=Rr6*Rr6*Rr6;
                 Rr12=Rr6*Rr6;
                 dr6=d2/rc2;
                 dr6=dr6*dr6*dr6;


                 (*ven)+=eij*(Rd12-2*Rd6+Rr6*(4.0-2*dr6)+Rr12*(2*dr6-3.0));
                 dven=-eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;
                 g=dven*dx;
                 (a1->GX)+=g;
                 (a2->GX)-=g;
                 g=dven*dy;
                 (a1->GY)+=g;
                 (a2->GY)-=g;
                 g=dven*dz;
                 (a1->GZ)+=g;
                 (a2->GZ)-=g;
              }
           }
        }
}



void eleng03(double f, struct atomgrp *ag, double eps, double* een,
             int n03, int* list03)
{
        int i, i1, i2;
        struct atom *a1, *a2;
        double chi, chj, chij, ee, deen, g;
        double dx, dy, dz, d, d2;
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           chi=a1->chrg;

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           chj=a2->chrg;

           chij=f*CCELEC*chi*chj/eps;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           d=sqrt(d2);
           ee=chij/d;
           (*een)+=ee;
           deen=ee/d2;
           g=deen*dx;
           (a1->GX)+=g;
           (a2->GX)-=g;
           g=deen*dy;
           (a1->GY)+=g;
           (a2->GY)-=g;
           g=deen*dz;
           (a1->GZ)+=g;
           (a2->GZ)-=g;
        }
}

void elengs03(double f, double rc, struct atomgrp *ag, double eps, double* een,
              int n03, int* list03)
{
        int i, i1, i2;
        double dx, dy, dz;
        double d2, d1, g, ch1, ch2;
        double esh, desh;
        struct atom *a1, *a2;

        double pf=f*CCELEC/eps;
        double rc2=1.0/(rc*rc);
        for(i=0; i<n03; i++)
        {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           ch1=pf*(a1->chrg);

           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           ch2=a2->chrg;

           dx=ag->atoms[i1].X-ag->atoms[i2].X;
           dy=ag->atoms[i1].Y-ag->atoms[i2].Y;
           dz=ag->atoms[i1].Z-ag->atoms[i2].Z;
           d2=dx*dx+dy*dy+dz*dz;
           d1=sqrt(d2);
           if(d1<rc)
           {
              esh=1.0-d1/rc;
              esh*=esh/d1;
              desh=(1.0/d2-rc2)/d1;
              ch2=ch2*ch1;
              (*een)+=ch2*esh;
              g=ch2*desh*dx;
              (a1->GX)+=g;
              (a2->GX)-=g;
              g=ch2*desh*dy;
              (a1->GY)+=g;
              (a2->GY)-=g;
              g=ch2*desh*dz;
              (a1->GZ)+=g;
              (a2->GZ)-=g;
           }
        }
}


void eleng(struct atomgrp *ag, double eps, double* een, struct nblist *nblst)
{
        int i, j, i1, n2, i2, *p;
        double dx, dy, dz;
        double d2, d1, g, ch1, ch2, x1, y1, z1;
        double esh, desh;
        struct atom *a1, *a2;

        double pf=CCELEC/eps;
        double rc=nblst->nbcof;
        double rc2=1.0/(rc*rc);
	for(i=0; i<nblst->nfat; i++)
	{
	   i1=nblst->ifat[i];
           a1=&(ag->atoms[i1]);
           ch1=pf*(a1->chrg);
           x1=ag->atoms[i1].X;
           y1=ag->atoms[i1].Y;
           z1=ag->atoms[i1].Z;
           n2=nblst->nsat[i];
           p=nblst->isat[i];
           for(j=0; j<n2; j++)
	   {
	      i2=p[j];
              a2=&(ag->atoms[i2]);
              ch2=a2->chrg;
              dx=x1-ag->atoms[i2].X;
              dy=y1-ag->atoms[i2].Y;
              dz=z1-ag->atoms[i2].Z;
              d2=dx*dx+dy*dy+dz*dz;
              d1=sqrt(d2);
              if(d1<rc)
              {
                 esh=1.0-d1/rc;
                 esh*=esh/d1;
                 desh=(1.0/d2-rc2)/d1;
                 ch2=ch2*ch1;
                 (*een)+=ch2*esh;
                 g=ch2*desh*dx;
                 (a1->GX)+=g;
                 (a2->GX)-=g;
                 g=ch2*desh*dy;
                 (a1->GY)+=g;
                 (a2->GY)-=g;
                 g=ch2*desh*dz;
                 (a1->GZ)+=g;
                 (a2->GZ)-=g;
              }
	   }
	}
}

void destroy_nblist(struct nblist *nblst)
{
        for(int i=0; i<nblst->nfat; i++)
        {
           if(nblst->nsat[i]>0)
              free(nblst->isat[i]);
        }
        free(nblst->isat);
        free(nblst->nsat);
        free(nblst->ifat);
        free(nblst->crds);
}
void free_nblist(struct nblist *nblst)
{
	destroy_nblist(nblst);
	free(nblst);
}

void gen_nblist_orig(struct atomgrp *ag, struct cubeset *cust, struct clusterset *clst,
                     int *excl_list, int **pd1, int **pd2, int ndm, struct nblist *nblst)
{
        int i, j, ic1, ic2, icl1, icl2;
        int k1, k2, ak1, ak2, ka1, ka2, ka3, ip;
        double mdc1, mdc2, dna, dda, dx, dy, dz, d;
        double x1, x2, y1, y2, z1, z2, nbcut;
        double dd, dista;

        int natoms=ag->natoms;

//prepare nonbond list arrays
        int *p;
        for(i=0; i<nblst->nfat; i++)
        {
           if((nblst->nsat[i])>0)
               free(nblst->isat[i]);
        }
        nblst->nfat=0;
        nblst->ifat=(int*)_mol_realloc(nblst->ifat,natoms*sizeof(int));
        nblst->nsat=(int*)_mol_realloc(nblst->nsat,natoms*sizeof(int));
        nblst->isat=(int**)_mol_realloc(nblst->isat,natoms*sizeof(int*));
//ifat is temporal
        int *ifat=(int*)_mol_malloc(natoms*sizeof(int));
        for(i=0; i<natoms; i++)
        {
           nblst->nsat[i]=0;
           nblst->ifat[i]=-1;
           ifat[i]=-1;
        }

        nbcut=nblst->nbcut;
        double nbcuts=nbcut*nbcut;

//loop over the filled cubes
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
//loop over all clusters in cube1
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
//loop over neighboring cubes including cube1
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
//second loop over clusters
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
//for clusters icl1 and icl2 go through all pairs of their atoms
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
//check the exclusion table
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
//check distances of distant atom pairs (d>dna)
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
//write a successful atom pair to a nonbond list structure
                             if(ifat[ka1]==-1)
                             {
                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
                                nblst->isat[nblst->nfat]=(int*)_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
                             p=nblst->isat[ip];
                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}
//cleanup, total pair calculation
        int ntotp=0;
        free(ifat);
        if(nblst->nfat > 0)  {
            nblst->ifat=(int*)_mol_realloc(nblst->ifat,(nblst->nfat)*sizeof(int));
            nblst->nsat=(int*)_mol_realloc(nblst->nsat,(nblst->nfat)*sizeof(int));
            nblst->isat=(int**)_mol_realloc(nblst->isat,(nblst->nfat)*sizeof(int*));
        }

        for(i=0; i<nblst->nfat; i++)
        {
           ntotp+=nblst->nsat[i];
           nblst->isat[i]=(int*)_mol_realloc(nblst->isat[i],(nblst->nsat[i])*sizeof(int));
        }
        nblst->npairs=ntotp;
}



void gen_nblist(struct atomgrp *ag, struct cubeset *cust, struct clusterset *clst,
                int *excl_list, int **pd1, int **pd2, int ndm, struct nblist *nblst)
{
        int i, j, ic1, ic2, icl1, icl2;
        int k1, k2, ak1, ak2, ka1, ka2, ka3, ip;
        double mdc1, mdc2, dna, dda, dx, dy, dz, d;
        double x1, x2, y1, y2, z1, z2, nbcut;
        double dd, dista;

        int natoms=ag->natoms;

//prepare nonbond list arrays
        int *p;
        for(i=0; i<nblst->nfat; i++)
        {
           if((nblst->nsat[i])>0)
               free(nblst->isat[i]);
        }
        nblst->nfat=0;
        nblst->ifat=(int*)_mol_realloc(nblst->ifat,natoms*sizeof(int));
        nblst->nsat=(int*)_mol_realloc(nblst->nsat,natoms*sizeof(int));
        nblst->isat=(int**)_mol_realloc(nblst->isat,natoms*sizeof(int*));
//ifat is temporal
        int *ifat=(int*)_mol_malloc(natoms*sizeof(int));
        for(i=0; i<natoms; i++)
        {
           nblst->nsat[i]=0;
           nblst->ifat[i]=-1;
           ifat[i]=-1;
        }

        nbcut=nblst->nbcut;
        double nbcuts=nbcut*nbcut;

//loop over the filled cubes
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
//loop over all clusters in cube1
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
//loop over neighboring cubes including cube1
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
//second loop over clusters
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
//for clusters icl1 and icl2 go through all pairs of their atoms
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
//check the exclusion table
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
//check distances of distant atom pairs (d>dna)
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
//write a successful atom pair to a nonbond list structure
                             if(ifat[ka1]==-1)
                             {
//                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
//                                nblst->isat[nblst->nfat]=(int*)_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
//                             p=nblst->isat[ip];
//                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}

        for(i=0; i<natoms; i++)
        {
           ip=ifat[i];
	   if ( ip > -1 )
              {
                if ( nblst->nsat[ip] > 0 ) nblst->isat[ip]=(int*)_mol_malloc(nblst->nsat[ip]*sizeof(int));
                nblst->nsat[ip]=0;
              } 
           nblst->ifat[i]=-1;           
           ifat[i]=-1;
        }
        nblst->nfat=0;	
	
//loop over the filled cubes
	for (i=0; i<cust->nfcubes; i++)
	{
           ic1=cust->ifcubes[i];
           icl1=cust->cubes[ic1].hstincube;
//loop over all clusters in cube1
	   while(icl1>=0)
           {
              mdc1=clst->clusters[icl1].mdc;
              x1=clst->clusters[icl1].gcent[0];
              y1=clst->clusters[icl1].gcent[1];
              z1=clst->clusters[icl1].gcent[2];
//loop over neighboring cubes including cube1
              for(j=-1; j<cust->cubes[ic1].nncubes; j++)
              {
                 if(j==-1)                                  // cube-self part
                    icl2=clst->clusters[icl1].nextincube;
                 else                                       // cube-neighbor cube part
                 {
                    ic2=cust->cubes[ic1].icubes[j];
                    icl2=cust->cubes[ic2].hstincube;
                 }
//second loop over clusters
                 while(icl2>=0)
                 {
                    mdc2=clst->clusters[icl2].mdc;
                    x2=clst->clusters[icl2].gcent[0];
                    y2=clst->clusters[icl2].gcent[1];
                    z2=clst->clusters[icl2].gcent[2];
                    dx=x1-x2;
                    dx*=dx;
                    dy=y1-y2;
                    dy*=dy;
                    dz=z1-z2;
                    dz*=dz;
                    d=dx+dy+dz;
                    dna=nbcut-mdc1-mdc2;
                    dna*=dna;
                    dda=nbcut+mdc1+mdc2;
                    dda*=dda;
                    if(d<=dda)
                    {
//for clusters icl1 and icl2 go through all pairs of their atoms
                       for(k1=0; k1<clst->clusters[icl1].natoms; k1++)
                       {
                          ak1=clst->clusters[icl1].iatom[k1];
                          for(k2=0; k2<clst->clusters[icl2].natoms; k2++)
                          {
                             ak2=clst->clusters[icl2].iatom[k2];
                             if(ag->atoms[ak1].fixed+ag->atoms[ak2].fixed==2)continue;
                             if(ak1<ak2)
                             {
                                ka1=ak1;
                                ka2=ak2;
                             }
                             else
                             {
                                ka1=ak2;
                                ka2=ak1;
                             }
//check the exclusion table
                             ka3=exta(ka1, ka2, excl_list, pd1, pd2, ndm);
                             if(ka3>0)continue;
//check distances of distant atom pairs (d>dna)
                             if(d>dna)
                             {
                                dista=0.0;
                                dd=ag->atoms[ka1].X-ag->atoms[ka2].X;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Y-ag->atoms[ka2].Y;
                                dista+=dd*dd;
                                dd=ag->atoms[ka1].Z-ag->atoms[ka2].Z;
                                dista+=dd*dd;
                                if(dista>nbcuts)continue;
                             }
//write a successful atom pair to a nonbond list structure
                             if(ifat[ka1]==-1)
                             {
                                nblst->ifat[nblst->nfat]=ka1;
                                ifat[ka1]=nblst->nfat;
//                                nblst->isat[nblst->nfat]=(int*)_mol_malloc((natoms-ka1-1)*sizeof(int));
                                (nblst->nfat)++;
                             }
                             ip=ifat[ka1];
                             p=nblst->isat[ip];
                             p[nblst->nsat[ip]]=ka2;
                             (nblst->nsat[ip])++;
                          }
                       }
                    }
                    icl2=clst->clusters[icl2].nextincube;
                 }
              }
              icl1=clst->clusters[icl1].nextincube;
           }
	}	
//cleanup, total pair calculation
        int ntotp=0;
        free(ifat);
        if(nblst->nfat > 0)  {
            nblst->ifat=(int*)_mol_realloc(nblst->ifat,(nblst->nfat)*sizeof(int));
            nblst->nsat=(int*)_mol_realloc(nblst->nsat,(nblst->nfat)*sizeof(int));
            nblst->isat=(int**)_mol_realloc(nblst->isat,(nblst->nfat)*sizeof(int*));
        }

        for(i=0; i<nblst->nfat; i++)
        {
           ntotp+=nblst->nsat[i];
//           nblst->isat[i]=_mol_realloc(nblst->isat[i],(nblst->nsat[i])*sizeof(int));
        }
        nblst->npairs=ntotp;
}


void free_cubeset(struct cubeset *cust)
{
        int i,ic;

        for (i=0; i<(cust->nfcubes); i++)
        {
             ic=cust->ifcubes[i];
             if (cust->cubes[ic].nncubes>0)
                   free(cust->cubes[ic].icubes);
        }
        free(cust->cubes);
        free(cust->ifcubes);
}

void gen_cubeset(double nbcut, struct clusterset *clst, struct cubeset *cust)
{
	double xmin, xmax, ymin, ymax, zmin, zmax, cubel;
        int i, nx, ny, nz, nxy, n, j, k, l;
        int    ix, iy, iz, ic, ic1;
        int    ix1, ix2, iy1, iy2, iz1, iz2;
        int    nnbr, temp[27];

        cubel=nbcut+clst->marg;
        cust->cubel=cubel;
        xmin=clst->clusters[0].gcent[0];
        xmax=xmin;
        ymin=clst->clusters[0].gcent[1];
        ymax=ymin;
        zmin=clst->clusters[0].gcent[2];
        zmax=zmin;
        for(i=1; i<clst->nclusters; i++)
	{
 	   if(xmin>clst->clusters[i].gcent[0])xmin=clst->clusters[i].gcent[0];
           if(xmax<clst->clusters[i].gcent[0])xmax=clst->clusters[i].gcent[0];
           if(ymin>clst->clusters[i].gcent[1])ymin=clst->clusters[i].gcent[1];
           if(ymax<clst->clusters[i].gcent[1])ymax=clst->clusters[i].gcent[1];
           if(zmin>clst->clusters[i].gcent[2])zmin=clst->clusters[i].gcent[2];
           if(zmax<clst->clusters[i].gcent[2])zmax=clst->clusters[i].gcent[2];
	}
        nx=1+(int)((xmax-xmin)/cubel);
        ny=1+(int)((ymax-ymin)/cubel);
        nz=1+(int)((zmax-zmin)/cubel);
        nxy=nx*ny;
        n=nxy*nz;
        cust->ncubes=n;
        cust->nfcubes=0;
        cust->cubes=(struct cube*)_mol_malloc(n*sizeof(struct cube));
        cust->ifcubes=(int*)_mol_malloc(n*sizeof(int));
        for(i=0; i<n; i++)
          {
	   cust->cubes[i].hstincube=-1;
	   cust->cubes[i].nncubes=0;
	  }
        for(i=0; i<clst->nclusters; i++)
        {
	   ix=(int)((clst->clusters[i].gcent[0]-xmin)/cubel);
           iy=(int)((clst->clusters[i].gcent[1]-ymin)/cubel);
           iz=(int)((clst->clusters[i].gcent[2]-zmin)/cubel);
           ic=ix+iy*nx+iz*nxy;
           if(cust->cubes[ic].hstincube==-1)
           {
              cust->cubes[ic].ix=ix;
              cust->cubes[ic].iy=iy;
              cust->cubes[ic].iz=iz;
              cust->ifcubes[(cust->nfcubes)++]=ic;
           }
	   clst->clusters[i].nextincube=cust->cubes[ic].hstincube;
	   cust->cubes[ic].hstincube=i;
        }
        for(i=0; i<cust->nfcubes; i++)
        {
	   ic=cust->ifcubes[i];
           if(cust->cubes[ic].hstincube==-1)continue;
           ix=cust->cubes[ic].ix;
              ix1=ix-1;
              if(ix1<0)ix1=0;
              ix2=ix+1;
              if(ix2>nx-1)ix2=nx-1;
           iy=cust->cubes[ic].iy;
              iy1=iy-1;
              if(iy1<0)iy1=0;
              iy2=iy+1;
              if(iy2>ny-1)iy2=ny-1;
           iz=cust->cubes[ic].iz;
              iz1=iz-1;
              if(iz1<0)iz1=0;
              iz2=iz+1;
              if(iz2>nz-1)iz2=nz-1;
           nnbr=0;
           for(j=iz1; j<=iz2; j++)
           {
              for(k=iy1; k<=iy2; k++)
              {
                 for(l=ix1; l<=ix2; l++)
                 {
                    ic1=l+k*nx+j*nxy;
                    if(ic1<=ic)continue;
                    if(cust->cubes[ic1].hstincube==-1)continue;
                    temp[nnbr++]=ic1;
                 }
              }
           }
           cust->cubes[ic].nncubes=nnbr;
           if(nnbr>0)
           {
              cust->cubes[ic].icubes=(int*)_mol_malloc(nnbr*sizeof(int));
              for(j=0; j<nnbr; j++)cust->cubes[ic].icubes[j]=temp[j];
           }
	}
}

void findmarg(struct atomgrp *ag, struct clusterset *clst)
{
	int i, j, k, a;
        double dist1=0, dist2=0;
        double dx, dy, dz, d, md, gc[3];
        for(i=0; i<clst->nclusters; i++)
	{
           md=0.0;
           for(j=0; j<3; j++)gc[j]=0.0;
           a=clst->clusters[i].natoms;
           for(j=0; j<a; j++)
           {
              k=clst->clusters[i].iatom[j];
              gc[0]+=(ag->atoms[k].X)/a;
              gc[1]+=(ag->atoms[k].Y)/a;
              gc[2]+=(ag->atoms[k].Z)/a;
           }
           for(j=0; j<a; j++)
           {
              k=clst->clusters[i].iatom[j];
              dx=ag->atoms[k].X-gc[0];
              dx*=dx;
              dy=ag->atoms[k].Y-gc[1];
              dy*=dy;
              dz=ag->atoms[k].Z-gc[2];
              dz*=dz;
              d=sqrt(dx+dy+dz);
              if(d>md)md=d;
              if(d>dist1)
              {
                 dist2=dist1;
                 dist1=d;
              }
              else if(d>dist2)dist2=d;
           }
           clst->clusters[i].mdc=md;
           for(j=0; j<3; j++)clst->clusters[i].gcent[j]=gc[j];
	}
        clst->marg=dist1+dist2;
}

void free_clset(struct clusterset *clst)
{
        int i;
        for(i=0; i<clst->nclusters; i++)
            free(clst->clusters[i].iatom);
        free(clst->clusters);
}

void gen_clset(int natoms, struct clusterset *clst,
               int nclust, int *clust)
{
	clst->nclusters=nclust;
        clst->clusters=(struct cluster*)_mol_malloc(nclust*sizeof(struct cluster));
        int i, j;
        for(i=0; i<nclust; i++)
            clst->clusters[i].natoms=0;
        for(i=0; i<natoms; i++)
	    (clst->clusters[clust[i]].natoms)++;
        for(i=0; i<nclust; i++)
	{
	    clst->clusters[i].iatom=(int*)_mol_malloc((clst->clusters[i].natoms)*sizeof(int));
            clst->clusters[i].natoms=0;
	}
        for(i=0; i<natoms; i++)
        {
            j=clst->clusters[clust[i]].natoms;
            clst->clusters[clust[i]].iatom[j++]=i;
            clst->clusters[clust[i]].natoms=j;
        }
}

void clust14(struct atomgrp* ag, int *nclust, int *clust)
{

	*nclust=0;
	int npair;
	int ma, mapair, ipair;
	int *pairs=(int*)_mol_malloc(MAXPAIR*2*sizeof(int));
	int i, j, k, l, m, n=ag->natoms;
        struct atombond *b;
	struct atom *a0, *a1, *a2;
	for(i=0; i<n; i++)
	{
		clust[i]=-1;
	}
        for(i=0; i<n; i++)
        {
	   if(clust[i]==-1)
	   {
		npair=0;
		a0=&(ag->atoms[i]);
		for(j=0; j<(a0->nbonds); j++)
		{
	  	   b=a0->bonds[j];
		   a1=b->a1;
		   if(a1==a0)a1=b->a0;
		   k=(a1->ingrp);
		   if(clust[k]==-1)
		   {
			addpair(i, k, &npair, pairs);
			for(l=0; l<(a1->nbonds); l++)
			{
			   b=a1->bonds[l];
			   a2=b->a1;
			   if(a2==a1)a2=b->a0;
			   if(a2==a0)continue;
                           m=(a2->ingrp);
                           if(clust[m]==-1)
				addpair(k, m, &npair, pairs);
			}
		   }
		}
		if(npair==0)clust[i]=*nclust;
		else
		{
		   mapair=2;
		   ipair=0;
		   for(l=0; l<npair; l++)
		   {
			ma=natpair(ag, &pairs[2*l], clust);
			if(ma>mapair)
			{
			   ipair=l;
			   mapair=ma;
			}
		   }
		   addclust(ag, &pairs[2*ipair], clust, *nclust);
		}
		(*nclust)++;
	   }
        }
	free(pairs);
}

void addpair(int i, int k, int *npair, int* pairs)
{
   int l=*npair;
   if(l+1>MAXPAIR)
   {
       printf("Increase MAXPAIR in clust14\n");
       exit (EXIT_FAILURE);
   }
   pairs[2*l]=i;
   pairs[2*l+1]=k;
   (*npair)+=1;
}

int natpair(struct atomgrp* ag, int *pair, int *clust)
{
	int n=0, j, k;
	struct atom *a1=&(ag->atoms[pair[0]]);
	struct atom *a2=&(ag->atoms[pair[1]]);
	struct atom *a;
	struct atombond *b;

        for(j=0; j<(a1->nbonds); j++)
        {
           b=a1->bonds[j];
           a=b->a1;
           if(a==a1)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)n++;
	}
        for(j=0; j<(a2->nbonds); j++)
        {
           b=a2->bonds[j];
           a=b->a1;
           if(a==a2)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)n++;
        }
	return n--;
}

void addclust(struct atomgrp* ag, int *pair, int *clust, int nclust)
{
        int j, k;
        struct atom *a1=&(ag->atoms[pair[0]]);
        struct atom *a2=&(ag->atoms[pair[1]]);
        struct atom *a;
        struct atombond *b;

        for(j=0; j<(a1->nbonds); j++)
        {
           b=a1->bonds[j];
           a=b->a1;
           if(a==a1)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)clust[k]=nclust;
        }
        for(j=0; j<(a2->nbonds); j++)
        {
           b=a2->bonds[j];
           a=b->a1;
           if(a==a2)a=b->a0;
           k=(a->ingrp);
           if(clust[k]==-1)clust[k]=nclust;
        }
}

void comp_list01(struct atomgrp* ag, int *list01, int *na01, int **pna01)
{
	int i, j, k, nb;
	struct atom *a, *a0;
	struct atombond *b;
	for(i=0; i<ag->natoms; i++)
	{
		na01[i]=0;
		a0=&(ag->atoms[i]);
		nb=a0->nbonds;
		pna01[i]=list01;
		for(j=0; j<nb; j++)
		{
			b=a0->bonds[j];
			a=b->a1;
			if(a==a0)a=b->a0;
			k=(a->ingrp);
			*list01=k;
			na01[i]++;
			list01++;
		}
		if(na01[i]==0)pna01[i]=NULL;
	}
}

void comp_n23(int natoms, int *na01, int **pna01, int *n02, int *n03)
{
	int i, j, k, l, m, *p;
	*n02=0;
	*n03=0;
	for(i=0; i<natoms; i++)
        {
	   j=na01[i];
	   if(j>1)
	   {
		(*n02)+=j*(j-1)/2;
		p=pna01[i];
		for(k=0; k<j; k++)
		{
		   l=p[k];
		   if(l>i)
		   {
			m=na01[l];
			if(m>1)(*n03)+=(j-1)*(m-1);
		   }
		}
	   }
	}
}

int trim_comp(const void *s1, const void *s2)
{
     int *v1, *v2;
     v1=(int *)s1;
     v2=(int *)s2;

     int i11=*v1;
     int i21=*v2;
     if(i11<i21)
        return -1;
     if(i11>i21)
        return 1;

     int i12=*(v1+1);
     int i22=*(v2+1);
     if(i12<i22)
        return -1;
     if(i12>i22)
        return 1;
     return 0;
}

void comp_list02(int natoms, int *na01, int **pna01,
                 int *na02, int **pna02, int *n02, int *list02)
{
	int i,j,k,l,m,n, *p;
        int *list=list02;

	for(i=0; i<natoms; i++)
        {
           j=na01[i];
           if(j>1)
           {
		p=pna01[i];
		for(k=0; k<j-1; k++)
		{
			n=p[k];
			for(l=k+1; l<j; l++)
			{
			   m=p[l];
			   if(m>n)
			   {
				*list02=n;
				list02++;
				*list02=m;
				list02++;
			   }
			   else
                           {
                                *list02=m;
                                list02++;
                                *list02=n;
                                list02++;
                           }
			}
		}
	   }
	}

        for(i=0; i<*n02; i++)
        {
           j=list[2*i];
           k=list[2*i+1];
           l=na01[j];
           if(l>0)
           {
              p=pna01[j];
              for(m=0; m<l; m++)
              {
                 n=p[m];
                 if(n==k)
                 {
                    list[2*i]=natoms;
                    list[2*i+1]=natoms;
                    break;
                 }
              }
           }
        }
        qsort(list, *n02, 2*sizeof(list[0]), trim_comp);
        j=(*n02)-1;
        for(i=j; i>0; i--)
        {
           if(list[2*i]==natoms)
               (*n02)--;
           else
               break;
        }
        for(i=0; i<natoms; i++)
        {
                na02[i]=0;
                pna02[i]=NULL;
        }
        for(i=0; i<*n02; i++)
        {
                j=list[2*i];
                (na02[j])++;
                if(i==0||j>list[2*i-2])
                   pna02[j]=&(list[2*i]);
        }
}

void comp_list03(int natoms, int *na01, int **pna01, int *list03)
{
	int i,j,k,l,m,n,q,r,s, *p, *o;
	for(i=0; i<natoms; i++)
	{
	   j=na01[i];
	   if(j>1)
	   {
		p=pna01[i];
		for(k=0; k<j; k++)
                {
		   l=p[k];
		   if(l>i)
		   {
			m=na01[l];
			if(m>1)
			{
			   o=pna01[l];
			   for(n=0; n<m; n++)
			   {
			      q=o[n];
			      if(q==i)continue;
			      for(r=0; r<j; r++)
			      {
			         if(r==k)continue;
			         s=p[r];
	                         if(q>s)
                           	 {
                                    *list03=s;
                                    list03++;
                                    *list03=q;
                                    list03++;
                           	 }
                                 else
                                 {
                                    *list03=q;
                                    list03++;
                                    *list03=s;
                                    list03++;
                                 }
			      }
			   }
			}
		   }
		}
	   }
	}
}

void fix_list03(struct atomgrp *ag, int n03, int* list03, int *nf03, int *listf03)
{
    int i, i1, i2;
    int n=0;
    struct atom *a1, *a2;
    for(i=0; i<n03; i++)
    {
           i1=list03[2*i];
           a1=&(ag->atoms[i1]);
           i2=list03[2*i+1];
           a2=&(ag->atoms[i2]);
           if(a2->fixed == 1 && a1->fixed == 1)continue;
           listf03[2*n]=i1;
           listf03[2*n+1]=i2;
           n++;
    }
    *nf03=n;
}

void trim_list03(int natoms,
                 int *na01, int **pna01,
                 int *na02, int **pna02,
                 int *n03, int *list03)
{
	qsort(list03, *n03, 2*sizeof(list03[0]), trim_comp);
        int i, j, k, l, m, n, o, *p;
        int *list=(int*)_mol_malloc(*n03*sizeof(int));
        for(i=1; i<*n03; i++)
        {
           j=2*i;
           list[i]=0;
           k=list03[j];
           l=list03[j+1];
// redundant pairs
           if(k==list03[j-2] && l==list03[j-1])
           {
              list[i]=1;
              continue;
           }
// directly bonded pairs
           m=na01[k];
           if(m>0)
           {
              p=pna01[k];
              for(n=0; n<m; n++)
              {
                 o=p[n];
                 if(l==o)
                 {
                    list[i]=1;
                    continue;
                 }
              }
           }
// bonded through an angle pairs
           m=na02[k];
           if(m>0)
           {
              p=pna02[k];
              for(n=0; n<m; n++)
              {
                 o=p[2*n+1];
                 if(l==o)
                 {
                    list[i]=1;
                    continue;
                 }
              }
           }
        }
        for(i=1; i<*n03; i++)
        {
           if(list[i]==1)
           {
               j=2*i;
               list03[j]=natoms;
               list03[j+1]=natoms;
           }
        }
        qsort(list03, *n03, 2*sizeof(list03[0]), trim_comp);
        j=*n03-1;
        for(i=j; i>0; i--)
        {
           if(list03[2*i]==natoms)(*n03)--;
           else break;
        }
        free(list);
//        for(i=0; i<*n03; i++)
//           printf("list03 %d, %d, %d, %d\n", i, *n03, list03[2*i], list03[2*i+1]);
//        exit(0);
}

void excl_dims(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int *nd1, int *nd2, int *ndm, int *atmind)
{
	int i,  j,  k, l, *p;
	for(i=0; i<2*natoms; i++)atmind[i]=0;
	*ndm=20;
	*nd2=0;
	*nd1=0;
	for(i=0; i<natoms; i++)
	{
		p=pna01[i];
		for(k=0; k<na01[i]; k++)
		{
		   j=p[k];
		   l=j-i;
		   if(l>*ndm)*ndm=l;
		   l--;
		   if(l>19 && atmind[2*i+1]==0)
		   {
			atmind[2*i+1]=1;
			(*nd2)++;
		   }
		   else if(l>9 && atmind[2*i]==0)
		   {
			atmind[2*i]=1;
			(*nd1)++;
		   }
		}
	}
	for(i=0; i<n02; i++)
	{
		j=list02[2*i];
		k=list02[2*i+1];
		l=k-j;
		if(l>*ndm)*ndm=l;
		l--;
                if(l>19 && atmind[2*j+1]==0)
                {
                     atmind[2*j+1]=1;
                     (*nd2)++;
                }
                else if(l>9 && atmind[2*j]==0)
                {
                     atmind[2*j]=1;
                     (*nd1)++;
                }
	}
        for(i=0; i<n03; i++)
        {
                j=list03[2*i];
                k=list03[2*i+1];
                l=k-j;
                if(l>*ndm)
                {
                   *ndm=l;
//                   printf("ndm %d %d %d\n",j,k,*ndm);
                }
                l--;
                if(l>19 && atmind[2*j+1]==0)
                {
                     atmind[2*j+1]=1;
                     (*nd2)++;
                }
                else if(l>9 && atmind[2*j]==0)
                {
                     atmind[2*j]=1;
                     (*nd1)++;
                }
        }
}

int exta(int a1, int a2, int *excl_list, int **pd1, int **pd2, int ndm)
{
	int i=a2-a1-1;
        if(i<0)
        {
           printf("exta: ERROR a1<=a2\n");
           exit (EXIT_FAILURE);
        }
        if(i>=ndm)
           return 0;
        if(i<10)
           return excl_list[10*a1+i];
        int *p;
        if(i<20)
        {
           p=pd1[a1];
           return p[i-10];
        }
        p=pd2[a1];
        return p[i-20];
}


void excl_tab(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int nd1, int nd2, int ndm, int *atmind,
               int *excl_list, int **pd1, int **pd2)
{
        int i,  j,  k, l, *p, *p1;
        int md1=-1, md2=-1;
// zero the exclusion list
        i=(natoms+nd1+1)*10+(nd2+1)*(ndm-20);
        for(j=0; j<i; j++)excl_list[j]=0;
// populate atom pointers to the exclusion list
        for(i=0; i<natoms; i++)
        {
           if(atmind[2*i]==1)
              pd1[i]=&(excl_list[10*(++md1+natoms)]);
           else
              pd1[i]=&(excl_list[10*(nd1+natoms)]);
           if(atmind[2*i+1]==1)
              pd2[i]=&(excl_list[10*(natoms+nd1+1)+(ndm-20)*(++md2)]);
           else
              pd2[i]=&(excl_list[10*(natoms+nd1+1)+(ndm-20)*nd2]);
        }
        for(i=0; i<natoms; i++)
        {
                p=pna01[i];
                for(k=0; k<na01[i]; k++)
                {
                   j=p[k];
                   l=j-i-1;
                   if(l>19)
                   {
                      p1=pd2[i];
                      p1[l-20]=1;
                   }
                   else if(l>9)
                   {
                      p1=pd1[i];
                      p1[l-10]=1;
                   }
                   else if(l>=0)
                   {
                      excl_list[10*i+l]=1;
                   }
                }
        }
        for(i=0; i<n02; i++)
        {
                j=list02[2*i];
                k=list02[2*i+1];
                l=k-j-1;
                if(l>19)
                {
                   p1=pd2[j];
                   p1[l-20]=2;
                }
                else if(l>9)
                {
                   p1=pd1[j];
                   p1[l-10]=2;
                }
                else if(l>=0)
                {
                   excl_list[10*j+l]=2;
                }
        }
        for(i=0; i<n03; i++)
        {
                j=list03[2*i];
                k=list03[2*i+1];
                l=k-j-1;
                if(l>19)
                {
                   p1=pd2[j];
                   p1[l-20]=3;
                }
                else if(l>9)
                {
                   p1=pd1[j];
                   p1[l-10]=3;
                }
                else if(l>=0)
                {
                   excl_list[10*j+l]=3;
                }
        }
}

//free nblist arrays in agsetup
void destroy_agsetup(struct agsetup* ags)
{
    free(ags->excl_list);
    free(ags->pd1);
    free(ags->pd2);
    free(ags->listf03);
    free(ags->list03);
    free(ags->list02);
    free_nblist(ags->nblst);
    free_clset(ags->clst);
    free(ags->clst);
}
void free_agsetup(struct agsetup* ags)
{
    destroy_agsetup(ags);
    free(ags);
}


//Setup nblist calculation
void init_nblst(struct atomgrp* ag, struct agsetup* ags)
{
        int *list01=(int*)_mol_malloc(2*(ag->nbonds)*sizeof(int));
        int *na01=(int*)_mol_malloc((ag->natoms)*sizeof(int));
	int **pna01=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
        int *na02=(int*)_mol_malloc((ag->natoms)*sizeof(int));
        int **pna02=(int**)_mol_malloc((ag->natoms)*sizeof(int*));

        comp_list01(ag, list01, na01, pna01);
        int i,j;

        int n02, n03,nf03;
	comp_n23(ag->natoms, na01, pna01, &n02, &n03);
        int *list02=(int*)_mol_malloc(2*n02*sizeof(int));
        int *list03=(int*)_mol_malloc(2*n03*sizeof(int));
	int *listf03=(int*)_mol_malloc(2*n03*sizeof(int));
	comp_list02(ag->natoms, na01, pna01, na02, pna02, &n02, list02);

        comp_list03(ag->natoms, na01, pna01, list03);
        trim_list03(ag->natoms, na01, pna01, na02, pna02, &n03, list03);
	fix_list03(ag,n03,list03,&nf03,listf03);
	if(nf03>0)
           listf03=(int*)_mol_realloc(listf03,2*nf03*sizeof(int));
        else
           listf03=(int*)_mol_realloc(listf03,2*sizeof(int));
        int nd1,nd2,ndm;
	int *atmind=(int*)_mol_malloc(2*(ag->natoms)*sizeof(int));
	excl_dims(ag->natoms, na01, pna01,
                       n02, list02, n03, list03,
                       &nd1, &nd2, &ndm, atmind);

	i=10*((ag->natoms)+nd1+1)+(nd2+1)*(ndm-20);
	int *excl_list=(int*)_mol_malloc(i*sizeof(int));
        int **pd1=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
        int **pd2=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
        excl_tab(ag->natoms, na01, pna01,
               n02, list02, n03, list03,
               nd1, nd2, ndm, atmind,
               excl_list, pd1, pd2);
	ags->list02=list02;
	ags->list03=list03;
	ags->listf03=listf03;
	ags->n02=n02;
	ags->n03=n03;
	ags->nf03=nf03;
	ags->ndm=ndm;
	ags->excl_list=excl_list;
	ags->pd1=pd1;
	ags->pd2=pd2;
	free(pna01);
        free(na01);
        free(list01);

        free(pna02);
        free(na02);
        free(atmind);
//***************Generating cluster------------------------
	int *clust=(int*)_mol_malloc((ag->natoms)*sizeof(int));
        int nclust;
        clust14(ag, &nclust, clust);

        struct clusterset *clst=(struct clusterset*)_mol_malloc(sizeof(struct clusterset));

        gen_clset(ag->natoms, clst, nclust, clust);
	ags->clst=clst;
//**********************free****************************************
        free(clust);

        struct nblist *nblst=(struct nblist*)_mol_malloc(sizeof(struct nblist));
        double nbcut=13.0;
        double nbcof=12.0;
        nblst->nbcut=nbcut;
        nblst->nbcof=nbcof;
        nblst->crds=(float*)_mol_malloc(3*(ag->natoms)*sizeof(float));
        nblst->nfat=0;
        nblst->ifat=(int*)_mol_malloc((ag->natoms)*sizeof(int));
        nblst->nsat=(int*)_mol_malloc((ag->natoms)*sizeof(int));
        for(j=0; j<(ag->natoms); j++)
                            nblst->nsat[j]=0;
        nblst->isat=(int**)_mol_malloc((ag->natoms)*sizeof(int*));
	ags->nblst=nblst;
}


void update_nblst(struct atomgrp* ag, struct agsetup* ags){
	int i;
        for(i=0; i<ag->natoms; i++)
        {
          ags->nblst->crds[3*i]=ag->atoms[i].X;
	  ags->nblst->crds[3*i+1]=ag->atoms[i].Y;
	  ags->nblst->crds[3*i+2]=ag->atoms[i].Z;
        }
        findmarg(ag, ags->clst);

        struct cubeset *cust = (struct cubeset*) _mol_malloc(sizeof(struct cubeset));
        gen_cubeset(ags->nblst->nbcut, ags->clst, cust);
        gen_nblist(ag, cust, ags->clst, ags->excl_list, ags->pd1, ags->pd2, ags->ndm, ags->nblst);
	free_cubeset(cust);
        free(cust);

}
//Checks wether cluster needs update
int check_clusterupdate(struct atomgrp* ag,struct agsetup* ags){
 int i=0,j;
 short flag=0;
 double delta0,delta;
 delta0=fabs(ags->nblst->nbcut-ags->nblst->nbcof)/2.0;
 delta0=delta0*delta0;
    for (i=0;i<ag->nactives;i++){
        j=ag->activelist[i];
        delta=(ags->nblst->crds[3*j]-ag->atoms[j].X)*(ags->nblst->crds[3*j]-ag->atoms[j].X)+
	   (ags->nblst->crds[3*j+1]-ag->atoms[j].Y)*(ags->nblst->crds[3*j+1]-ag->atoms[j].Y)+
	   (ags->nblst->crds[3*j+2]-ag->atoms[j].Z)*(ags->nblst->crds[3*j+2]-ag->atoms[j].Z);
        delta=delta*delta;
	if (delta>delta0){
	    flag=1;
	    break;
	}
	}
    if (flag){
	update_nblst(ag,ags);
	return 1;
	    } else {return 0;}
}

void give_012(struct atomgrp* ag, struct agsetup* ags,
              int* na012, int** la012, int* nf012, int** lf012)
{
     struct atom *a0, *a1;
     int ia0, ia1;
     int n01=ag->nbonds;
     int n02=ags->n02;
     int n=n01+n02;
     int na=0;
     int nf=0;

     if(n>0)
     {
        *la012=(int*)_mol_malloc(2*n*sizeof(int));
        *lf012=(int*)_mol_malloc(2*n*sizeof(int));
        int i;
        for(i=0; i<n01; i++)
        {
           a0=ag->bonds[i].a0;
           a1=ag->bonds[i].a1;
           if(a0->fixed == 0 || a1->fixed == 0)
           {
              *la012[2*na]=a0->ingrp;
              *la012[2*na+1]=a1->ingrp;
              na++;
           }
           else
           {
              *lf012[2*nf]=a0->ingrp;
              *lf012[2*nf+1]=a1->ingrp;
              nf++;
           }
        }
        for(i=0; i<n02; i++)
        {
           ia0=ags->list02[2*i];
           ia1=ags->list02[2*i+1];
           if(ag->atoms[ia0].fixed == 0 || ag->atoms[ia1].fixed == 0)
           {
              *la012[2*na]=ia0;
              *la012[2*na+1]=ia1;
              na++;
           }
           else
           {
              *lf012[2*nf]=ia0;
              *lf012[2*nf+1]=ia1;
              nf++;
           }
        }
        if(na>0)*la012=(int*)_mol_realloc(*la012, 2*na*sizeof(int));
        else free(*la012);
        if(nf>0)*lf012=(int*)_mol_realloc(*lf012, 2*nf*sizeof(int));
        else free(*lf012);
     }
     *na012=na;
     *nf012=nf;
}

void springeng(struct springset *sprst, double* spren)
{
     int i,j,nat,lj;
     double xtot,ytot,ztot,fk;
     struct atomgrp *ag;
     for (i=0; i<sprst->nsprings; i++)
     {
         nat=sprst->springs[i].naspr;
         if(nat>0)
         {
            xtot=0.0;
            ytot=0.0;
            ztot=0.0;
            ag =sprst->springs[i].agp;
            for (j=0; j<nat; j++)
            {
               lj=sprst->springs[i].laspr[j];
               xtot+=ag->atoms[lj].X;
               ytot+=ag->atoms[lj].Y;
               ztot+=ag->atoms[lj].Z;
            }
            xtot=xtot/nat-sprst->springs[i].X0;
            ytot=ytot/nat-sprst->springs[i].Y0;
            ztot=ztot/nat-sprst->springs[i].Z0;
            fk=sprst->springs[i].fkspr;
            *spren+=fk*(xtot*xtot+ytot*ytot+ztot*ztot);
            fk=2*fk/nat;
            for (j=0; j<nat; j++)
            {
                lj=sprst->springs[i].laspr[j];
                ag->atoms[lj].GX -= xtot*fk;
                ag->atoms[lj].GY -= ytot*fk;
                ag->atoms[lj].GZ -= ztot*fk;
            }
         }
     }
}

void modify_vdw_save(struct atomgrp* ag, int nmod, int* mod, double* enew, double* rnew,
		                                             double* eold, double* rold)
{
   int i;
   for(i=0; i<nmod; i++)
   {
      eold[i]=ag->atoms[mod[i]].eps;
      if(enew[i]<=0)ag->atoms[mod[i]].eps=sqrt(-enew[i]);
      rold[i]=ag->atoms[mod[i]].rminh;
      ag->atoms[mod[i]].rminh=rnew[i];
   }

}

void modify_vdw_back(struct atomgrp* ag, int nmod, int* mod, double* eold, double* rold)
{
   int i;
   for(i=0; i<nmod; i++)
   {
      ag->atoms[mod[i]].eps=eold[i];
      ag->atoms[mod[i]].rminh=rold[i];
   }

}

void get_mod_vdw_all(double lambe, double lambr, struct atomgrp* ag, int *nmod, int **mod, double **modeps, double **modrminh)
{
   int i;
   *nmod=ag->natoms;
   *mod=(int*)mymalloc(*nmod*sizeof(int));
   *modeps=(double*)mymalloc(*nmod*sizeof(double));
   *modrminh=(double*)mymalloc(*nmod*sizeof(double));
   for(i=0; i<*nmod; i++)
   {
      (*mod)[i]=i;
      (*modeps)[i]=lambe*(ag->atoms[i].eps);
      (*modrminh)[i]=lambr*(ag->atoms[i].rminh);
   }
}



/* ####################################################################### */
/* #                       START: ADDED BY REZAUL                        # */
/* ####################################################################### */

inline int bonded( mol_atom *ai, mol_atom *aj )
{
   if ( ai->nbonds > aj->nbonds )
     {
       mol_atom *t = ai;
       ai = aj;
       aj = t;
     }

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *b = ai->bonds[ i ];

       if ( ( b->a0 == aj ) || ( b->a1 == aj ) ) return 1;
     }

   return 0;
}


inline int bonded2( mol_atom *ai, mol_atom *aj )
{
   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       for ( int j = 0; j < aj->nbonds; j++ )
         {
           mol_bond *bj = aj->bonds[ j ];

           if ( ( bi->a0 == bj->a0 ) || ( bi->a0 == bj->a1 )
             || ( bi->a1 == bj->a0 ) || ( bi->a1 == bj->a1 ) ) return 1;
         }
     }

   return 0;
}


inline int within_3_bonds( mol_atom *ai, mol_atom *aj )
{
   if ( ai->nbonds > aj->nbonds )
     {
       mol_atom *t = ai;
       ai = aj;
       aj = t;
     }

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       if ( ap == aj ) return 1;

       for ( int p = 0; p < ap->nbonds; p++ )
         {
           mol_bond *bp = ap->bonds[ p ];

           mol_atom *aq = bp->a0;
           if ( aq == ap ) aq = bp->a1;
           if ( aq == ai ) continue;

           if ( aq == aj ) return 1;

           for ( int q = 0; q < aq->nbonds; q++ )
             {
               mol_bond *bq = aq->bonds[ q ];

               mol_atom *ar = bq->a0;
               if ( ar == aq ) ar = bq->a1;
               if ( ( ar == ai ) || ( ar == ap ) ) continue;

               if ( ar == aj ) return 1;
             }
         }
     }

   return 0;
}



inline int within_3_bonds_bak( mol_atom *ai, mol_atom *aj )
{
   if ( ai->nbonds > aj->nbonds )
     {
       mol_atom *t = ai;
       ai = aj;
       aj = t;
     }

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       if ( ap == aj ) return 1;
     }

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       for ( int p = 0; p < ap->nbonds; p++ )
         {
           mol_bond *bp = ap->bonds[ p ];

           mol_atom *aq = bp->a0;
           if ( aq == ap ) aq = bp->a1;
           if ( aq == ai ) continue;

           if ( aq == aj ) return 1;
         }
     }

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       for ( int p = 0; p < ap->nbonds; p++ )
         {
           mol_bond *bp = ap->bonds[ p ];

           mol_atom *aq = bp->a0;
           if ( aq == ap ) aq = bp->a1;
           if ( aq == ai ) continue;

           for ( int q = 0; q < aq->nbonds; q++ )
             {
               mol_bond *bq = aq->bonds[ q ];

               mol_atom *ar = bq->a0;
               if ( ar == aq ) ar = bq->a1;
               if ( ( ar == ai ) || ( ar == ap ) ) continue;

               if ( ar == aj ) return 1;
             }
         }
     }


   return 0;
}



inline void mark_atoms_within_3_bonds( mol_atom *ai )
{
   ai->octree_ptr = - ai->octree_ptr;

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       if ( ap->octree_ptr < 0 ) continue;

       ap->octree_ptr = - ap->octree_ptr;

       for ( int p = 0; p < ap->nbonds; p++ )
         {
           mol_bond *bp = ap->bonds[ p ];

           mol_atom *aq = bp->a0;
           if ( aq == ap ) aq = bp->a1;

           if ( aq->octree_ptr < 0 ) continue;

           aq->octree_ptr = - aq->octree_ptr;

           for ( int q = 0; q < aq->nbonds; q++ )
             {
               mol_bond *bq = aq->bonds[ q ];

               mol_atom *ar = bq->a0;
               if ( ar == aq ) ar = bq->a1;

               if ( ar->octree_ptr < 0 ) continue;

               ar->octree_ptr = - ar->octree_ptr;
             }
         }
     }
}


inline void unmark_atoms_within_3_bonds( mol_atom *ai )
{
   ai->octree_ptr = - ai->octree_ptr;

   for ( int i = 0; i < ai->nbonds; i++ )
     {
       mol_bond *bi = ai->bonds[ i ];

       mol_atom *ap = bi->a0;
       if ( ap == ai ) ap = bi->a1;

       if ( ap->octree_ptr >= 0 ) continue;

       ap->octree_ptr = - ap->octree_ptr;

       for ( int p = 0; p < ap->nbonds; p++ )
         {
           mol_bond *bp = ap->bonds[ p ];

           mol_atom *aq = bp->a0;
           if ( aq == ap ) aq = bp->a1;

           if ( aq->octree_ptr >= 0 ) continue;

           aq->octree_ptr = - aq->octree_ptr;

           for ( int q = 0; q < aq->nbonds; q++ )
             {
               mol_bond *bq = aq->bonds[ q ];

               mol_atom *ar = bq->a0;
               if ( ar == aq ) ar = bq->a1;

               if ( ar->octree_ptr >= 0 ) continue;

               ar->octree_ptr = - ar->octree_ptr;
             }
         }
     }
}



inline int exclude( int a1, int a2, int *excl_list, int **pd1, int **pd2, int ndm )
{
   if ( a1 <= a2 )
     {
       int i = a2 - a1 - 1;

       if ( i >= ndm ) return 0;
       if ( i < 10 ) return ( excl_list[ 10 * a1 + i ] );

       int *p;

       if ( i < 20 )
         {
           p = pd1[ a1 ];
           return ( p[ i - 10 ] );
         }

       p = pd2[ a1 ];
       return ( p[ i - 20 ] );
     }
   else
     {
       int i = a1 - a2 - 1;

       if ( i >= ndm ) return 0;
       if ( i < 10 ) return ( excl_list[ 10 * a2 + i ] );

       int *p;

       if ( i < 20 )
         {
           p = pd1[ a2 ];
           return ( p[ i - 10 ] );
         }

       p = pd2[ a2 ];
       return ( p[ i - 20 ] );
     }
}


inline double get_pairwise_vdweng_near( int ai, mol_atom *atom_i, double xi, double yi, double zi, double e_i, double r_i,
				        OCTREE *octree_static, int aj, double rc2, double irc2, struct agsetup *ags )
{
   mol_atom *atom_j = &( octree_static->atoms[ aj ] );

   double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

   double d2 = dx * dx + dy * dy + dz * dz;

   if ( ( d2 < 0.0001 ) || ( d2 >= rc2 ) ) return 0;

//   if ( ( d2 < 20 ) && exclude( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm ) ) return 0;

   if ( d2 < 20 )
     {
       int ka1, ka2, ka3;

       if ( ai < aj )
         {
           ka1 = ai;
           ka2 = aj;
         }
       else
         {
           ka1 = aj;
           ka2 = ai;
         }

       ka3 = exta( ka1, ka2, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );

       if ( ka3 > 0 ) return 0;
     }

//   if ( ( d2 < 20 ) && within_3_bonds( atom_i, atom_j ) ) return 0;

   double id2 = 1 / d2;

   double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
   r_ij *= r_ij;

   double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;

   Rd6 = Rd6 * Rd6 * Rd6;
   Rr6 = Rr6 * Rr6 * Rr6;
   dr6 = dr6 * dr6 * dr6;

   double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;

   double en = e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );

   double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;

   double g = dven * dx;

   ( atom_i->GX ) += g;
   ( atom_j->GX ) -= g;

   g = dven * dy;

   ( atom_i->GY ) += g;
   ( atom_j->GY ) -= g;

   g = dven * dz;

   ( atom_i->GZ ) += g;
   ( atom_j->GZ ) -= g;

   return en;
}


inline double get_pairwise_vdweng_far( int ai, mol_atom *atom_i, double xi, double yi, double zi, double e_i, double r_i,
				       OCTREE *octree_static, int aj, double rc2, double irc2 )
{
   mol_atom *atom_j = &( octree_static->atoms[ aj ] );

   double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

   double d2 = dx * dx + dy * dy + dz * dz;

   if ( d2 >= rc2 ) return 0;

   double id2 = 1 / d2;

   double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
   r_ij *= r_ij;

   double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;

   Rd6 = Rd6 * Rd6 * Rd6;
   Rr6 = Rr6 * Rr6 * Rr6;
   dr6 = dr6 * dr6 * dr6;

   double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;

   double en = e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );

   double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;

   double g = dven * dx;

   ( atom_i->GX ) += g;
   ( atom_j->GX ) -= g;

   g = dven * dy;

   ( atom_i->GY ) += g;
   ( atom_j->GY ) -= g;

   g = dven * dz;

   ( atom_i->GZ ) += g;
   ( atom_j->GZ ) -= g;

   return en;
}



double get_pairwise_vdweng_near_node( int ai, mol_atom *atom_i, double xi, double yi, double zi, double e_i, double r_i,
				      OCTREE *octree_static, OCTREE_NODE *snode, double rc2, double irc2, struct agsetup *ags )
{
   double en = 0, gx = 0, gy = 0, gz = 0;
   int nf = snode->nfixed;

   for ( int j = 0; j < snode->n; j++ )
      {
        int aj = snode->indices[ j ];

        if ( ( j >= nf ) && ( aj <= ai ) ) continue;

        mol_atom *atom_j = &( octree_static->atoms[ aj ] );

        double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

        double d2 = dx * dx + dy * dy + dz * dz;

        if ( d2 >= rc2 ) continue;

//        if ( d2 < 20 )
          {
            int k;

            if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
            else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );

            if ( k > 0 ) continue;
          }

        double id2 = 1 / d2;

        double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
        r_ij *= r_ij;

        double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;

        Rd6 = Rd6 * Rd6 * Rd6;
        Rr6 = Rr6 * Rr6 * Rr6;
        dr6 = dr6 * dr6 * dr6;

        double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;

        en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );

        double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;

        double g = dven * dx;

        gx += g;
        ( atom_j->GX ) -= g;

        g = dven * dy;

        gy += g;
        ( atom_j->GY ) -= g;

        g = dven * dz;

        gz += g;
        ( atom_j->GZ ) -= g;
      }

   ( atom_i->GX ) += gx;
   ( atom_i->GY ) += gy;
   ( atom_i->GZ ) += gz;

   return en;
}



double get_pairwise_vdweng_far_node( int ai, mol_atom *atom_i, double xi, double yi, double zi, double e_i, double r_i,
				     OCTREE *octree_static, OCTREE_NODE *snode, double rc2, double irc2 )
{
   double en = 0, gx = 0, gy = 0, gz = 0;
   int nf = snode->nfixed;

   for ( int j = 0; j < snode->n; j++ )
      {
        int aj = snode->indices[ j ];

        if ( ( j >= nf ) && ( aj <= ai ) ) continue;

        mol_atom *atom_j = &( octree_static->atoms[ aj ] );

        double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

        double d2 = dx * dx + dy * dy + dz * dz;

        if ( d2 >= rc2 ) continue;

        double id2 = 1 / d2;

        double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
        r_ij *= r_ij;

        double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;

        Rd6 = Rd6 * Rd6 * Rd6;
        Rr6 = Rr6 * Rr6 * Rr6;
        dr6 = dr6 * dr6 * dr6;

        double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;

        en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );

        double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;

        double g = dven * dx;

        gx += g;
        ( atom_j->GX ) -= g;

        g = dven * dy;

        gy += g;
        ( atom_j->GY ) -= g;

        g = dven * dz;

        gz += g;
        ( atom_j->GZ ) -= g;
      }

   ( atom_i->GX ) += gx;
   ( atom_i->GY ) += gy;
   ( atom_i->GZ ) += gz;

   return en;
}



//#define DEBUG_VDW

void vdweng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy )
{
    OCTREE *octree_static = octpar->octree_static;
    OCTREE *octree_moving = octpar->octree_moving;
    double *trans_mat = octpar->trans;

    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) octpar->proc_func_params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ octpar->node_static ] );
    OCTREE_NODE *mnode = &( octree_moving->nodes[ octpar->node_moving ] );

    double rc2 = octpar->dist_cutoff * octpar->dist_cutoff;
    double irc2 = 1 / rc2;

    int nf = snode->nfixed;

    double en = 0;

#ifdef DEBUG_VDW
    int m = 0, n = 0;
    double minD2 = 100000000000.0;
    double maxD2 = 0;
#endif

    if ( ( trans_mat != NULL ) || ( mnode->n - mnode->nfixed <= snode->n ) )
      {
        for ( int i = mnode->nfixed; i < mnode->n; i++ )
          {
            int ai = mnode->indices[ i ];
            mol_atom *atom_i = &( octree_moving->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h
      
            double d2 = min_pt2bx_dist2( snode->lx, snode->ly, snode->lz, snode->dim, x, y, z );

#ifdef DEBUG_VDW                        
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_VDW              
                m++;
#endif      
                double e_i = atom_i->eps, r_i = atom_i->rminh;
      
                double gx = 0, gy = 0, gz = 0;
                      
                for ( int j = 0; j < snode->n; j++ )
                  {
                    int aj = snode->indices[ j ];
      
                    if ( ( j >= nf ) && ( aj <= ai ) ) continue;
      
                    mol_atom *atom_j = &( octree_static->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                    d2 = dx * dx + dy * dy + dz * dz;
                    
                    if ( d2 >= rc2 ) continue;
      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;
                      }

#ifdef DEBUG_VDW      
                    n++;
#endif      
                    double id2 = 1 / d2;
      
                    double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
                    r_ij *= r_ij;
      
                    double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;
      
                    Rd6 = Rd6 * Rd6 * Rd6;
                    Rr6 = Rr6 * Rr6 * Rr6;
                    dr6 = dr6 * dr6 * dr6;
      
                    double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;
      
                    en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );
      
                    double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;
      
                    double g = dven * dx;
      
                    gx += g;
                    ( atom_j->GX ) -= g;
      
                    g = dven * dy;
      
                    gy += g;
                    ( atom_j->GY ) -= g;
      
                    g = dven * dz;
      
                    gz += g;
                    ( atom_j->GZ ) -= g;
                  }
      
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }
    else
      {
        for ( int i = 0; i < snode->n; i++ )
          {
            int ai = snode->indices[ i ];
            mol_atom *atom_i = &( octree_static->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            double d2 = min_pt2bx_dist2( mnode->lx, mnode->ly, mnode->lz, mnode->dim, x, y, z );
      
#ifdef DEBUG_VDW            
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_VDW              
                m++;
#endif                
      
                double e_i = atom_i->eps, r_i = atom_i->rminh;
      
                double gx = 0, gy = 0, gz = 0;
      
                for ( int j = mnode->nfixed; j < mnode->n; j++ )
                  {
                    int aj = mnode->indices[ j ];
      
                    if ( ( i >= nf ) && ( aj >= ai ) ) continue;
      
                    mol_atom *atom_j = &( octree_moving->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                    d2 = dx * dx + dy * dy + dz * dz;
                    
                    if ( d2 >= rc2 ) continue;
      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;
                      }

#ifdef DEBUG_VDW      
                    n++;
#endif                    
      
                    double id2 = 1 / d2;
      
                    double e_ij = e_i * atom_j->eps, r_ij = r_i + atom_j->rminh;
                    r_ij *= r_ij;
      
                    double Rd6 = r_ij * id2, Rr6 = r_ij * irc2, dr6 = d2 * irc2;
      
                    Rd6 = Rd6 * Rd6 * Rd6;
                    Rr6 = Rr6 * Rr6 * Rr6;
                    dr6 = dr6 * dr6 * dr6;
      
                    double Rd12 = Rd6 * Rd6, Rr12 = Rr6 * Rr6;
      
                    en += e_ij * ( Rd12 - 2 * Rd6 + Rr6 * ( 4.0 - 2 * dr6 ) + Rr12 * ( 2 * dr6 - 3.0 ) );
      
                    double dven = e_ij * 12 * ( -Rd12 + Rd6 + dr6 * ( Rr12 - Rr6 ) ) * id2;
      
                    double g = dven * dx;
      
                    gx += g;
                    ( atom_j->GX ) -= g;
      
                    g = dven * dy;
      
                    gy += g;
                    ( atom_j->GY ) -= g;
      
                    g = dven * dz;
      
                    gz += g;
                    ( atom_j->GZ ) -= g;
                  }
      
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }  

#ifdef DEBUG_VDW
    int d = within_distance_cutoff( snode->lx, snode->ly, snode->lz, snode->dim,
                                    mnode->lx, mnode->ly, mnode->lz, mnode->dim,
                                    octpar->dist_cutoff );                                                           
                                         
    printf( "< %2d %2d %9lf > < %2d %2d %9lf > %2d %3d %15.12lf %9lf %9lf %d ( %9lf, %9lf, %9lf | %9lf %9lf %9lf )\n", 
            mnode->nfixed, mnode->n, mnode->dim, snode->nfixed, snode->n, snode->dim, m, n, en, sqrt( minD2 ), sqrt( maxD2 ), 
            d, snode->lx, snode->ly, snode->lz, mnode->lx, mnode->ly, mnode->lz );
#endif

    *energy = en;
}



inline double get_pairwise_eleng_near( int ai, mol_atom *atom_i, double xi, double yi, double zi, double c_i,
				       OCTREE *octree_static, int aj, double rc, double irc2, struct agsetup *ags )
{
   mol_atom *atom_j = &( octree_static->atoms[ aj ] );

   double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

   double d2 = dx * dx + dy * dy + dz * dz;
   double rc2 = rc * rc;

   if ( ( d2 < 0.0001 ) || ( d2 >= rc2 ) ) return 0;

//   if ( ( d2 < 20 ) && exclude( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm ) ) return 0;

   if ( d2 < 20 )
     {
       int ka1, ka2, ka3;

       if ( ai < aj )
         {
           ka1 = ai;
           ka2 = aj;
         }
       else
         {
           ka1 = aj;
           ka2 = ai;
         }

       ka3 = exta( ka1, ka2, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );

       if ( ka3 > 0 ) return 0;
     }

//   if ( ( d2 < 20 ) && within_3_bonds( atom_i, atom_j ) ) return 0;

   double d1 = sqrt( d2 );
   double id2 = 1 / d2;

   double c_ij = c_i * atom_j->chrg;

   double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;

   esh *= ( esh / d1 );

   double en = c_ij * esh;

   double g = c_ij * desh * dx;

   ( atom_i->GX ) += g;
   ( atom_j->GX ) -= g;

   g = c_ij * desh * dy;

   ( atom_i->GY ) += g;
   ( atom_j->GY ) -= g;

   g = c_ij * desh * dz;

   ( atom_i->GZ ) += g;
   ( atom_j->GZ ) -= g;

   return en;
}



inline double get_pairwise_eleng_far( int ai, mol_atom *atom_i, double xi, double yi, double zi, double c_i,
				       OCTREE *octree_static, int aj, double rc, double irc2 )
{
   mol_atom *atom_j = &( octree_static->atoms[ aj ] );

   double dx = atom_j->X - xi, dy = atom_j->Y - yi, dz = atom_j->Z - zi;

   double d2 = dx * dx + dy * dy + dz * dz;
   double rc2 = rc * rc;

   if ( d2 >= rc2 ) return 0;

   double d1 = sqrt( d2 );
   double id2 = 1 / d2;

   double c_ij = c_i * atom_j->chrg;

   double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;

   esh *= ( esh / d1 );

   double en = c_ij * esh;

   double g = c_ij * desh * dx;

   ( atom_i->GX ) += g;
   ( atom_j->GX ) -= g;

   g = c_ij * desh * dy;

   ( atom_i->GY ) += g;
   ( atom_j->GY ) -= g;

   g = c_ij * desh * dz;

   ( atom_i->GZ ) += g;
   ( atom_j->GZ ) -= g;

   return en;
}




void eleng_octree_single_mol_bak( OCTREE *octree_static, int node_static, OCTREE *octree_moving, int node_moving,
                                  void *params, double dist_cutoff, double *trans_mat, double *energy )
{
    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ node_static ] );
    OCTREE_NODE *mnode = &( octree_moving->nodes[ node_moving ] );

    if ( ( snode->nfixed == snode->n ) && ( mnode->nfixed == mnode->n ) ) return;

    double pf = CCELEC / prms->eps;
    double rc = dist_cutoff;
    double rc2 = rc * rc;
    double irc2 = 1 / rc2;

    double half_dim = 0.5 * snode->dim;
    double sx = snode->lx + half_dim,
           sy = snode->ly + half_dim,
           sz = snode->lz + half_dim;
    double hd2 = ( dist_cutoff + HALF_SQRT_THREE * snode->dim ) * ( dist_cutoff + HALF_SQRT_THREE * snode->dim );
    double ld2 = ( 4.5 + HALF_SQRT_THREE * snode->dim ) * ( 4.5 + HALF_SQRT_THREE * snode->dim );

    *energy = 0;

    for ( int i = mnode->nfixed; i < mnode->n; i++ )
      {
        int ai = mnode->indices[ i ];
        mol_atom *atom_i = &( octree_moving->atoms[ ai ] );

        double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;

        if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h

        double dx = sx - x, dy = sy - y, dz = sz - z;
        double d2 = dx * dx + dy * dy + dz * dz;

        if ( d2 < hd2 )
          {
            double c_i = pf * atom_i->chrg;

            if ( d2 < ld2 )
              {
                for ( int j = 0; j < snode->nfixed; j++ )
                   {
                     int aj = snode->indices[ j ];
                     ( *energy ) += get_pairwise_eleng_near( ai, atom_i, x, y, z, c_i, octree_static, aj, rc, irc2, ags );
                   }

                for ( int j = snode->nfixed; j < snode->n; j++ )
                   {
                     int aj = snode->indices[ j ];

                     if ( aj > ai ) ( *energy ) += get_pairwise_eleng_near( ai, atom_i, x, y, z, c_i, octree_static, aj, rc, irc2, ags );
                   }
              }
            else
              {
                for ( int j = 0; j < snode->nfixed; j++ )
                   {
                     int aj = snode->indices[ j ];
                     ( *energy ) += get_pairwise_eleng_far( ai, atom_i, x, y, z, c_i, octree_static, aj, rc, irc2 );
                   }

                for ( int j = snode->nfixed; j < snode->n; j++ )
                   {
                     int aj = snode->indices[ j ];

                     if ( aj > ai ) ( *energy ) += get_pairwise_eleng_far( ai, atom_i, x, y, z, c_i, octree_static, aj, rc, irc2 );
                   }
              }
          }
      }
}


//#define DEBUG_ELENG

void eleng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy )
{
    OCTREE *octree_static = octpar->octree_static;
    OCTREE *octree_moving = octpar->octree_moving;
    double *trans_mat = octpar->trans;

    OCTREE_PARAMS *prms = ( OCTREE_PARAMS * ) octpar->proc_func_params;
    struct agsetup *ags = prms->ags;

    OCTREE_NODE *snode = &( octree_static->nodes[ octpar->node_static ] );
    OCTREE_NODE *mnode = &( octree_moving->nodes[ octpar->node_moving ] );

    double pf = CCELEC / prms->eps;
    double rc = octpar->dist_cutoff;
    double rc2 = rc * rc;
    double irc2 = 1 / rc2;

    int nf = snode->nfixed;

    double en = 0;

#ifdef DEBUG_ELENG
    int m = 0, n = 0;
    double minD2 = 100000000000.0;
    double maxD2 = 0;
#endif

    if ( ( trans_mat != NULL ) || ( mnode->n - mnode->nfixed <= snode->n ) )
      {
        for ( int i = mnode->nfixed; i < mnode->n; i++ )
          {
            int ai = mnode->indices[ i ];
            mol_atom *atom_i = &( octree_moving->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            if ( trans_mat != NULL ) transform_point( x, y, z, trans_mat, &x, &y, &z );    // defined in octree.h
      
            double d2 = min_pt2bx_dist2( snode->lx, snode->ly, snode->lz, snode->dim, x, y, z );

#ifdef DEBUG_ELENG                        
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_ELENG              
                m++;
#endif      
                double c_i = pf * atom_i->chrg;
      
                double gx = 0, gy = 0, gz = 0;
                      
                for ( int j = 0; j < snode->n; j++ )
                  {
                    int aj = snode->indices[ j ];
      
                    if ( ( j >= nf ) && ( aj <= ai ) ) continue;
      
                    mol_atom *atom_j = &( octree_static->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                    d2 = dx * dx + dy * dy + dz * dz;
                    
                    if ( d2 >= rc2 ) continue;
      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;
                      }

#ifdef DEBUG_ELENG      
                    n++;
#endif      
                    double d1 = sqrt( d2 );
                    double id2 = 1 / d2;
                    
                    double c_ij = c_i * atom_j->chrg;
                    
                    double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;
                    
                    esh *= ( esh / d1 );
                    
                    en += c_ij * esh;
                    
                    double g = c_ij * desh * dx;
                    
                    gx += g;
                    ( atom_j->GX ) -= g;
                    
                    g = c_ij * desh * dy;
                    
                    gy += g;
                    ( atom_j->GY ) -= g;
                    
                    g = c_ij * desh * dz;
                    
                    gz += g;
                    ( atom_j->GZ ) -= g;
                  }
      
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }
    else
      {
        for ( int i = 0; i < snode->n; i++ )
          {
            int ai = snode->indices[ i ];
            mol_atom *atom_i = &( octree_static->atoms[ ai ] );
      
            double x = atom_i->X, y = atom_i->Y, z = atom_i->Z;
      
            double d2 = min_pt2bx_dist2( mnode->lx, mnode->ly, mnode->lz, mnode->dim, x, y, z );
      
#ifdef DEBUG_ELENG            
            if ( d2 < minD2 ) minD2 = d2;
            if ( d2 > maxD2 ) maxD2 = d2;                
#endif
      
            if ( d2 < rc2 )
              {
#ifdef DEBUG_ELENG              
                m++;
#endif                      
                double c_i = pf * atom_i->chrg;
      
                double gx = 0, gy = 0, gz = 0;
      
                for ( int j = mnode->nfixed; j < mnode->n; j++ )
                  {
                    int aj = mnode->indices[ j ];
      
                    if ( ( i >= nf ) && ( aj >= ai ) ) continue;
      
                    mol_atom *atom_j = &( octree_moving->atoms[ aj ] );
      
                    double dx = atom_j->X - x, 
                           dy = atom_j->Y - y,
                           dz = atom_j->Z - z;
      
                    d2 = dx * dx + dy * dy + dz * dz;
                    
                    if ( d2 >= rc2 ) continue;
      
                    if ( d2 < 20 )
                      {
                        int k;
      
                        if ( ai < aj ) k = exta( ai, aj, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
                        else k = exta( aj, ai, ags->excl_list, ags->pd1, ags->pd2, ags->ndm );
      
                        if ( k > 0 ) continue;
                      }

#ifdef DEBUG_ELENG      
                    n++;
#endif                    

                    double d1 = sqrt( d2 );
                    double id2 = 1 / d2;
                    
                    double c_ij = c_i * atom_j->chrg;
                    
                    double esh = 1.0 - d1 / rc, desh = ( id2 - irc2 ) / d1;
                    
                    esh *= ( esh / d1 );
                    
                    en += c_ij * esh;
                    
                    double g = c_ij * desh * dx;
                    
                    gx += g;
                    ( atom_j->GX ) -= g;
                    
                    g = c_ij * desh * dy;
                    
                    gy += g;
                    ( atom_j->GY ) -= g;
                    
                    g = c_ij * desh * dz;
                    
                    gz += g;
                    ( atom_j->GZ ) -= g;                          
                  }
      
                ( atom_i->GX ) += gx;
                ( atom_i->GY ) += gy;
                ( atom_i->GZ ) += gz;
              }
          }
      }  

#ifdef DEBUG_ELENG
    int d = within_distance_cutoff( snode->lx, snode->ly, snode->lz, snode->dim,
                                    mnode->lx, mnode->ly, mnode->lz, mnode->dim,
                                    octpar->dist_cutoff );                                                           
                                         
    printf( "< %2d %2d %9lf > < %2d %2d %9lf > %2d %3d %15.12lf %9lf %9lf %d ( %9lf, %9lf, %9lf | %9lf %9lf %9lf )\n", 
            mnode->nfixed, mnode->n, mnode->dim, snode->nfixed, snode->n, snode->dim, m, n, en, sqrt( minD2 ), sqrt( maxD2 ), 
            d, snode->lx, snode->ly, snode->lz, mnode->lx, mnode->ly, mnode->lz );
#endif

    *energy = en;
}


/* ####################################################################### */
/* #                        END: ADDED BY REZAUL                         # */
/* ####################################################################### */

};
