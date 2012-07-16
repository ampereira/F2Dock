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

#define small 0.0000001

void beng(struct atomgrp *ag, double* en)
{
	int i;
	struct atom *a0, *a1;
	double dx, dy, dz, l, dl, e1;
        struct atombond* bp;
	for(i=0; i<ag->nbact; i++)
	{
                bp=ag->bact[i];
                a0=bp->a0;
                a1=bp->a1;
//                if(a0->fixed == 1 && a1->fixed == 1)continue;
                dx=a0->X-a1->X;
		dy=a0->Y-a1->Y;
		dz=a0->Z-a1->Z;
		l=sqrt(dx*dx+dy*dy+dz*dz);
		dl=l-(bp->l0);
		e1=(bp->k)*dl;
		(*en)+=e1*dl;
		e1=-2*e1/l;
		dl=e1*dx;
		(a0->GX)+=dl;
		(a1->GX)-=dl;
		dl=e1*dy;
		(a0->GY)+=dl;
		(a1->GY)-=dl;
		dl=e1*dz;
                (a0->GZ)+=dl;
                (a1->GZ)-=dl;
	}
}

void aeng(struct atomgrp *ag, double* en)
{
        double DEGRA = acos(-1.0)/180.0;
        int i;
        struct atom *a0, *a1, *a2;
	double dx10, dx12, dy10, dy12, dz10, dz12;
	double dsx10, dsy10, dsz10, ls10, l10;
	double dsx12, dsy12, dsz12, ls12, l12, l02;
	double dth, sth, e1, e2, e3;
	double yzp, yzm, xzp, xzm, xyp, xym;
        struct atomangle* ap;
        for(i=0; i<ag->nangact; i++)
        {
                ap=ag->angact[i];
                a0=ap->a0;
                a1=ap->a1;
                a2=ap->a2;
//                if(a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1)continue;
		dx10=a0->X-a1->X;
		dx12=a2->X-a1->X;
                dy10=a0->Y-a1->Y;
                dy12=a2->Y-a1->Y;
                dz10=a0->Z-a1->Z;
                dz12=a2->Z-a1->Z;

		dsx10=dx10*dx10;
		dsy10=dy10*dy10;
		dsz10=dz10*dz10;

		ls10=dsx10+dsy10+dsz10;
		l10=sqrt(ls10);

		dsx12=dx12*dx12;
		dsy12=dy12*dy12;
		dsz12=dz12*dz12;

		ls12=dsx12+dsy12+dsz12;
		l12=sqrt(ls12);

		l02=l10*l12;
                dth=acos((dx10*dx12+dy10*dy12+dz10*dz12)/l02);
		sth=sin(dth);

		if(sth < small)
		{
			fprintf (stderr, 
			"angle %d (%d, %d, %d) is close to linear\n", i, a0->ingrp, a1->ingrp, a2->ingrp);
                        fprintf (stderr,
                        "accuracy of forces will be compromised\n");
			sth=small;
		}

		dth-=DEGRA*(ap->th0);
		e1=(ap->k)*dth;
                (*en)+=e1*dth;
		e1*=2.0/sth/l02;
// 10 bond
		e2=e1/ls10;
		
		yzp=dsy10+dsz10;
		yzm=-dy10*dz10;
		xzp=dsx10+dsz10;
		xzm=-dx10*dz10;
		xyp=dsx10+dsy10;
		xym=-dx10*dy10;

		e3=e2*(yzp*dx12+xym*dy12+xzm*dz12);
		(a0->GX)+=e3;
		(a1->GX)-=e3;

		e3=e2*(xym*dx12+xzp*dy12+yzm*dz12);
		(a0->GY)+=e3;
		(a1->GY)-=e3;

		e3=e2*(xzm*dx12+yzm*dy12+xyp*dz12);
		(a0->GZ)+=e3;
                (a1->GZ)-=e3;
// 12 bond
                e2=e1/ls12;

                yzp=dsy12+dsz12;
                yzm=-dy12*dz12;
                xzp=dsx12+dsz12;
                xzm=-dx12*dz12;
                xyp=dsx12+dsy12;
                xym=-dx12*dy12;

                e3=e2*(yzp*dx10+xym*dy10+xzm*dz10);
                (a2->GX)+=e3;
                (a1->GX)-=e3;

                e3=e2*(xym*dx10+xzp*dy10+yzm*dz10);
                (a2->GY)+=e3;
                (a1->GY)-=e3;

                e3=e2*(xzm*dx10+yzm*dy10+xyp*dz10);
                (a2->GZ)+=e3;
                (a1->GZ)-=e3;
        }
}

void ieng(struct atomgrp *ag, double* en)
{
        double PI=acos(-1.0);
        double PI2=2*PI;
        double DEGRA = PI/180.0;
   
        int i;
        struct atom *a0, *a1, *a2, *a3;
	double dx01, dy01, dz01;
        double dx12, dy12, dz12;
        double dx23, dy23, dz23;
        double dx02, dy02, dz02;
        double dx13, dy13, dz13;
	double vx02, vy02, vz02;
        double vx13, vy13, vz13;
        double vx03, vy03, vz03;
	double ds02, ds13, d12;
	double xco, ysi, impan, imp0, dimp, e1, e2;
	double x1, y1, z1, x2, y2, z2, dx, dy, dz;
        struct atomimproper* ip;
	
        for(i=0; i<ag->nimpact; i++)
        {
                ip=ag->impact[i];
                a0=ip->a0;
                a1=ip->a1;
                a2=ip->a2;
		a3=ip->a3;
//                if(a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1 && a3->fixed == 1)continue;
		dx01=a1->X-a0->X;
		dy01=a1->Y-a0->Y;
                dz01=a1->Z-a0->Z;
                dx12=a2->X-a1->X;
                dy12=a2->Y-a1->Y;
                dz12=a2->Z-a1->Z;
                dx23=a3->X-a2->X;
                dy23=a3->Y-a2->Y;
                dz23=a3->Z-a2->Z;
                dx02=a2->X-a0->X;
                dy02=a2->Y-a0->Y;
                dz02=a2->Z-a0->Z;
                dx13=a3->X-a1->X;
                dy13=a3->Y-a1->Y;
                dz13=a3->Z-a1->Z;
// 01x12
		vx02=dy01*dz12-dy12*dz01;
		vy02=dz01*dx12-dz12*dx01;
		vz02=dx01*dy12-dx12*dy01;
// 12x23
		vx13=dy12*dz23-dy23*dz12;
                vy13=dz12*dx23-dz23*dx12;
                vz13=dx12*dy23-dx23*dy12;
// (01x12)x(12x23)
		vx03=vy02*vz13-vy13*vz02;
		vy03=vz02*vx13-vz13*vx02;
		vz03=vx02*vy13-vx13*vy02;
// lengths
		ds02=vx02*vx02+vy02*vy02+vz02*vz02;
		ds13=vx13*vx13+vy13*vy13+vz13*vz13;
		d12=sqrt(dx12*dx12+dy12*dy12+dz12*dz12);

		xco=vx02*vx13+vy02*vy13+vz02*vz13;
		ysi=(dx12*vx03+dy12*vy03+dz12*vz03)/d12;
		impan=atan2(ysi,xco);
		imp0=DEGRA*(ip->psi0);
		dimp=impan-imp0;
		while(dimp>PI)dimp-=PI2;
		while(dimp<-PI)dimp+=PI2;
		e1=(ip->k)*dimp;
		(*en)+=e1*dimp;

		e1*=2.0/d12;
		e2=e1/ds13;
		e1/=(-ds02);

		x1=e1*(vy02*dz12-dy12*vz02);
		y1=e1*(vz02*dx12-dz12*vx02);
		z1=e1*(vx02*dy12-dx12*vy02);

		x2=e2*(vy13*dz12-dy12*vz13);
                y2=e2*(vz13*dx12-dz12*vx13);
                z2=e2*(vx13*dy12-dx12*vy13);

		dx=dz12*y1-dy12*z1;
		dy=dx12*z1-dz12*x1;
		dz=dy12*x1-dx12*y1;

		(a0->GX)+=dx;
                (a0->GY)+=dy;
                (a0->GZ)+=dz;

		dx=dy02*z1-dz02*y1+dz23*y2-dy23*z2;
		dy=dz02*x1-dx02*z1+dx23*z2-dz23*x2;
		dz=dx02*y1-dy02*x1+dy23*x2-dx23*y2;

                (a1->GX)+=dx;
                (a1->GY)+=dy;
                (a1->GZ)+=dz;

		dx=dz01*y1-dy01*z1+dy13*z2-dz13*y2;
		dy=dx01*z1-dz01*x1+dz13*x2-dx13*z2;
		dz=dy01*x1-dx01*y1+dx13*y2-dy13*x2;

                (a2->GX)+=dx;
                (a2->GY)+=dy;
                (a2->GZ)+=dz;

		dx=dz12*y2-dy12*z2;
                dy=dx12*z2-dz12*x2;
                dz=dy12*x2-dx12*y2;

                (a3->GX)+=dx;
                (a3->GY)+=dy;
                (a3->GZ)+=dz;
	}
}

void teng(struct atomgrp *ag, double* en)
{
        long double PI=acos(-1.0);
        long double DEGRA=PI/180.0;
        int i, n;
        struct atom *a0, *a1, *a2, *a3;
        double dx01, dy01, dz01;
        double dx12, dy12, dz12;
        double dx23, dy23, dz23;
        double dx02, dy02, dz02;
        double dx13, dy13, dz13;
        double vx02, vy02, vz02;
        double vx13, vy13, vz13;
        double vx03, vy03, vz03;
        double ds02, ds13, d12;
        double xco, ysi, toran, tor0, dtor, e1, e2;
        double x1, y1, z1, x2, y2, z2, dx, dy, dz;
        struct atomtorsion* tp;

        for(i=0; i<ag->ntoract; i++)
        {
                tp=ag->toract[i];
                a0=tp->a0;
                a1=tp->a1;
                a2=tp->a2;
                a3=tp->a3;
//                if(a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1 && a3->fixed == 1)continue;
                dx01=a1->X-a0->X;
                dy01=a1->Y-a0->Y;
                dz01=a1->Z-a0->Z;
                dx12=a2->X-a1->X;
                dy12=a2->Y-a1->Y;
                dz12=a2->Z-a1->Z;
                dx23=a3->X-a2->X;
                dy23=a3->Y-a2->Y;
                dz23=a3->Z-a2->Z;
                dx02=a2->X-a0->X;
                dy02=a2->Y-a0->Y;
                dz02=a2->Z-a0->Z;
                dx13=a3->X-a1->X;
                dy13=a3->Y-a1->Y;
                dz13=a3->Z-a1->Z;
// 01x12
                vx02=dy01*dz12-dy12*dz01;
                vy02=dz01*dx12-dz12*dx01;
                vz02=dx01*dy12-dx12*dy01;
// 12x23
                vx13=dy12*dz23-dy23*dz12;
                vy13=dz12*dx23-dz23*dx12;
                vz13=dx12*dy23-dx23*dy12;
// (01x12)x(12x23)
                vx03=vy02*vz13-vy13*vz02;
                vy03=vz02*vx13-vz13*vx02;
                vz03=vx02*vy13-vx13*vy02;
// lengths
                ds02=vx02*vx02+vy02*vy02+vz02*vz02;
                ds13=vx13*vx13+vy13*vy13+vz13*vz13;
                d12=sqrt(dx12*dx12+dy12*dy12+dz12*dz12);

                xco=vx02*vx13+vy02*vy13+vz02*vz13;
                ysi=(dx12*vx03+dy12*vy03+dz12*vz03)/d12;
                toran=atan2(ysi,xco);
                tor0=DEGRA*(tp->d);
		n=tp->n;
                dtor=n*toran-tor0;
                e1=tp->k;
                (*en)+=e1*(1.0+cos(dtor));

                e1*=(-n*sin(dtor)/d12);
                e2=e1/ds13;
                e1/=(-ds02);

                x1=e1*(vy02*dz12-dy12*vz02);
                y1=e1*(vz02*dx12-dz12*vx02);
                z1=e1*(vx02*dy12-dx12*vy02);

                x2=e2*(vy13*dz12-dy12*vz13);
                y2=e2*(vz13*dx12-dz12*vx13);
                z2=e2*(vx13*dy12-dx12*vy13);

                dx=dz12*y1-dy12*z1;
                dy=dx12*z1-dz12*x1;
                dz=dy12*x1-dx12*y1;

                (a0->GX)+=dx;
                (a0->GY)+=dy;
                (a0->GZ)+=dz;

                dx=dy02*z1-dz02*y1+dz23*y2-dy23*z2;
                dy=dz02*x1-dx02*z1+dx23*z2-dz23*x2;
                dz=dx02*y1-dy02*x1+dy23*x2-dx23*y2;

                (a1->GX)+=dx;
                (a1->GY)+=dy;
                (a1->GZ)+=dz;

                dx=dz01*y1-dy01*z1+dy13*z2-dz13*y2;
                dy=dx01*z1-dz01*x1+dz13*x2-dx13*z2;
                dz=dy01*x1-dx01*y1+dx13*y2-dy13*x2;

                (a2->GX)+=dx;
                (a2->GY)+=dy;
                (a2->GZ)+=dz;

                dx=dz12*y2-dy12*z2;
                dy=dx12*z2-dz12*x2;
                dz=dy12*x2-dx12*y2;

                (a3->GX)+=dx;
                (a3->GY)+=dy;
                (a3->GZ)+=dz;
        }
}

void zero_grads(struct atomgrp *ag)
{
	int i;
	for(i=0; i<ag->natoms; i++)
	{
		ag->atoms[i].GX=0.0;
		ag->atoms[i].GY=0.0;
		ag->atoms[i].GZ=0.0;
	}	
}

void check_b_grads(struct atomgrp *ag, double d, void (*efun)(struct atomgrp *, double*))
{
        int n=ag->natoms, i;
        double en, en1, t;
        double *fs=(double*)_mol_malloc(3*n*sizeof(double));
//en0
        en=0;
        (*efun)(ag, &en);

        for(i=0; i<n; i++)
        {
//x
                en1=0;
                t=ag->atoms[i].X;
                ag->atoms[i].X=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].X=t;
                fs[3*i]=(en-en1)/d;
//y
                en1=0;
                t=ag->atoms[i].Y;
                ag->atoms[i].Y=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].Y=t;
                fs[3*i+1]=(en-en1)/d;
//z
                en1=0;
                t=ag->atoms[i].Z;
                ag->atoms[i].Z=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].Z=t;
                fs[3*i+2]=(en-en1)/d;
        }
        en=0;
        zero_grads(ag);
        (*efun)(ag, &en);
        for(i=0; i<n; i++)
        {
                printf("%d calculated: %lf %lf %lf\n",
                i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
                printf("%d numerical : %lf %lf %lf\n",
                i,fs[3*i],fs[3*i+1],fs[3*i+2]);
        }
        free(fs);
}

void check_epeng_grads(struct atomgrp *ag, struct grid* pot, double d, 
		void (*efun)(double, struct atomgrp *, double*, struct grid*))
{
	int n=ag->natoms, i;
	double en, en1, t;
	double *fs=(double*)_mol_malloc(3*n*sizeof(double));
//en0
	en=0;
//	(*efun)(ag, &en);
        (*efun)(1.0,ag, &en, pot);

        for(i=0; i<n; i++)
	{
//x
		en1=0;
		t=ag->atoms[i].X;
		ag->atoms[i].X=d+t;
//		(*efun)(ag, &en1);
                (*efun)(1.0,ag, &en1, pot);
		ag->atoms[i].X=t;
		fs[3*i]=(en-en1)/d;
//y
		en1=0;
                t=ag->atoms[i].Y;
                ag->atoms[i].Y=d+t;
//                (*efun)(ag, &en1);
                  (*efun)(1.0,ag, &en1, pot);
                ag->atoms[i].Y=t;
                fs[3*i+1]=(en-en1)/d;
//z
                en1=0;
                t=ag->atoms[i].Z;
                ag->atoms[i].Z=d+t;
//                (*efun)(ag, &en1);
                  (*efun)(1.0,ag, &en1, pot);
                ag->atoms[i].Z=t;
                fs[3*i+2]=(en-en1)/d;
	}
	en=0;
	zero_grads(ag);
//	(*efun)(ag, &en);
        (*efun)(1.0,ag, &en, pot);
	for(i=0; i<n; i++)
        {
		printf("%d calculated: %lf %lf %lf\n",
		i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
		printf("%d numerical : %lf %lf %lf\n",
                i,fs[3*i],fs[3*i+1],fs[3*i+2]);
	}
	free(fs);
}

void check_speng_grads(int nstart, int nend,
                struct atomgrp *ag, double d, double stens,
                double* hx0, double* hy0, double* hz0, 
                int nx, int ny, int nz,
                double dcx, double dcy, double dcz,
                double cx, double cy, double cz, double w,
                void (*efun)(double, double, struct atomgrp *, double*, 
                double*, double*, double*, 
                int,int,int,
                double,double,double,
                double,double,double,double))
{
        int n=ag->natoms, i;
        double en, en1, t;
        double *fs=(double*)_mol_malloc(3*n*sizeof(double));
//en0
        en=0;
//      (*efun)(ag, &en);
//        (*efun)(1.0,ag, &en, pot);
        (*efun)( 1.0, stens, ag, &en,
                 hx0,hy0,hz0,
                 nx,ny,nz,
                 dcx, dcy, dcz,
                 cx,cy,cz,w);

        for(i=nstart; i<=nend; i++)
        {
                printf("i= %d\n",i);
//x
                en1=0;
                t=ag->atoms[i].X;
                ag->atoms[i].X=d+t;
//              (*efun)(ag, &en1);
//              (*efun)(1.0,ag, &en1, pot);
                (*efun)( 1.0, stens, ag, &en1,
                 hx0,hy0,hz0,
                 nx,ny,nz,
                 dcx, dcy, dcz,
                 cx,cy,cz,w);
                ag->atoms[i].X=t;
                fs[3*i]=(en-en1)/d;
//y
                en1=0;
                t=ag->atoms[i].Y;
                ag->atoms[i].Y=d+t;
//                (*efun)(ag, &en1);
//                (*efun)(1.0,ag, &en1, pot);
                (*efun)( 1.0, stens, ag, &en1,
                 hx0,hy0,hz0,
                 nx,ny,nz,
                 dcx, dcy, dcz,
                 cx,cy,cz,w);
                ag->atoms[i].Y=t;
                fs[3*i+1]=(en-en1)/d;
//z
                en1=0;
                t=ag->atoms[i].Z;
                ag->atoms[i].Z=d+t;
//                (*efun)(ag, &en1);
//                (*efun)(1.0,ag, &en1, pot);
                (*efun)( 1.0, stens, ag, &en1,
                 hx0,hy0,hz0,
                 nx,ny,nz,
                 dcx, dcy, dcz,
                 cx,cy,cz,w);
                ag->atoms[i].Z=t;
                fs[3*i+2]=(en-en1)/d;
        }
        en=0;
        zero_grads(ag);
//      (*efun)(ag, &en);
//      (*efun)(1.0,ag, &en, pot);
        (*efun)( 1.0, stens, ag, &en,
                 hx0,hy0,hz0,
                 nx,ny,nz,
                 dcx, dcy, dcz,
                 cx,cy,cz,w);
        for(i=nstart; i<=nend; i++)
        {
                printf("%d calculated: %lf %lf %lf\n",
                i,ag->atoms[i].GX,ag->atoms[i].GY,ag->atoms[i].GZ);
                printf("%d numerical : %lf %lf %lf\n",
                i,fs[3*i],fs[3*i+1],fs[3*i+2]);
        }
        free(fs);
}

};
