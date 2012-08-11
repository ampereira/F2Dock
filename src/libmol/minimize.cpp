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
#include <string.h>
#include <math.h>
#include <float.h>
#include <libmol/libmol.h>
#include <libmol/minimize.h>

namespace LibMol{

//Powell/dirPowell is below


void bracket(double* orig, double* dir, double step,
             int ndim, void* prms,
             void (*egfun)(int , double* , void* , double*, double* ),
             double *fb,  double *la, double *lb, double *lc)
{
         int mindim=ndim;
         double* newD=(double*)_mol_malloc(mindim*sizeof(double) );

         double fa,fc,fnew,lnew;
         double G = 1.618034;
         double TooSmall=pow(10, -10);

         int i, j;

         *la=0;
         *lb=step;
         egfun(ndim, orig, prms, &fa, NULL);
         for (i =0; i<mindim ; i++)
                     newD[i] = orig[i]+step*dir[i];
         egfun(ndim, newD, prms, fb, NULL);
         if(fa<*fb)
         {
            double temp;
            temp = fa;
            fa = *fb;
            *fb = temp;
            temp =*la;
            *la=*lb;
            *lb=temp;
         }

         *lc=*lb+G*(*lb-*la);
//       printf("la %f lb %f lc %f\n", *la, *lb, *lc);
         for  (i =0; i<mindim ; i++)
             newD[i] = orig[i]+*lc*dir[i];
            
         egfun(ndim, newD, prms, &fc, NULL);
        
         double parabp1, parabp2;

         while( (fc<=*fb))
         {
//       printf("lb=%f", *lb);

             parabp1 = (*lb -*la)*(*fb-fc);
             parabp2 = (*lb -*lc)*(*fb-fa);
             if( fabs(parabp1-parabp2)>TooSmall)
             {
                 lnew=*lb+((*lb-*lc)*parabp2-(*lb-*la)*parabp1)/(2*(parabp1-parabp2));
                 if(((*lb-lnew)*(lnew-*lc))>0)
                 {
                     for(j =0; j<mindim ; j++)
                        newD[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, newD, prms, &fnew, NULL);
                     if(fnew<fc)
                     {
                        *la=*lb;
                        *fb=fnew;
                        *lb=lnew;
                        free(newD);
                        return;
                     }
                     else if(fnew>*fb)
                     {
                        *lc=lnew;
                        free(newD);
                        return;
                     }
                     lnew=*lc+G*(*lc-*lb);
                     for(j =0; j<mindim ; j++)
                         newD[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, newD, prms, &fnew, NULL);
                 }
                 else
                 {
                     for(j =0; j<mindim ; j++)
                            newD[j] = orig[j] + lnew*(dir[j]);
                     egfun(ndim, newD, prms, &fnew, NULL);
                     if(fnew<fc)
                     {
                         *lb=*lc;
                         *lc=lnew;
                         lnew=*lc+G*(*lc-*lb);
                         *fb=fc;
                         fc=fnew;
                         for(j =0; j<mindim ; j++)
                                newD[j] = orig[j] + lnew*(dir[j]);
                         egfun(ndim, newD, prms, &fnew, NULL);
                     } 
                 }
             }
             else
             {
                 lnew=*lc+G*(*lc-*lb);
                 for(j =0; j<mindim ; j++)
                      newD[j] = orig[j] + lnew*(dir[j]);
                 egfun(ndim, newD, prms, &fnew, NULL);
             }
             *la=*lb;
             *lb=*lc;
             *lc=lnew;
             fa=*fb;
             *fb=fc;
             fc=fnew;
         }

         free(newD);
}
void old_brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, int maxtimes,
           double*  min, double* fmim)
{
           int j;
           double* newD=(double*)_mol_malloc(ndim*sizeof(double) );
           double s=pow(10,-10);
           double CGOLD=0.3819660;
           int iter;
           double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
           double e=0.0; 
           a=(la < lc ? la : lc); 
           b=(la > lc ? la : lc); 
           x=w=v=lb; 
//           for(j=0; j<ndim; j++)
//                 newD[j]=orig[j]+x*dir[j]; 
//           egfun(ndim, newD, prms, &fx, NULL); 
           fw=fv=fx=fb;
           for (iter=1;iter<=maxtimes;iter++) 
           { 
                xm=0.5*(a+b);
                tol2=2.0*(tol1=tol*fabs(x)+s);
                if (fabs(x-xm) <= (tol2-0.5*(b-a))) 
                {
                    for(j=0; j<ndim; j++)
                        min[j]=orig[j]+x*dir[j];
                    *fmim=fx; 
                    free(newD);
                    return ;
                }
                if (fabs(e) > tol1) 
                {
                    r=(x-w)*(fx-fv);
                    q=(x-v)*(fx-fw);
                    p=(x-v)*q-(x-w)*r;
                    q=2.0*(q-r);
                    if (q > 0.0) p = -p;
                    q=fabs(q);
                    etemp=e;
                    e=d;
                    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                                                    d=CGOLD*(e=(x >= xm ? a-x : b-x));
                    else 
                    {
                          d=p/q; 
                          u=x+d;
                          if (u-a < tol2 || b-u < tol2)
                                  d=((xm-x) >= 0.0 ? fabs(tol1) : -fabs(tol1));
                    }
                } 
                else 
                {
                    d=CGOLD*(e=(x >= xm ? a-x : b-x));
                }
                u=(fabs(d) >= tol1 ? x+d : x+((d) >= 0.0 ? fabs(tol1) : -fabs(tol1)));
                for(j=0; j<ndim; j++)
                    newD[j]=orig[j]+u*dir[j];
                egfun(ndim, newD, prms, &fu, NULL);
                if (fu <= fx) 
                { 
                      if(u >= x) a=x; 
                      else b=x; 
                      v=w;w=x;x=u;
                      fv=fw;fw=fx;fx=fu;
                } 
                else 
                {
                    if (u < x) a=u; 
                    else b=u;
                    if (fu <= fw || w == x) 
                    {
                        v=w;
                        w=u;
                        fv=fw;
                        fw=fu;
                    } 
                    else if (fu <= fv || v == x || v == w) 
                    {
                       v=u;
                       fv=fu;
                    }
                } 
           } 
           printf("Too many iterations in brent\n");
           for(j=0; j<ndim; j++)
                        min[j]=orig[j]+x*dir[j];
           *fmim=fx;
           free(newD); 
}



void brent(double* orig, double* dir, 
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, unsigned int maxtimes,
           double*  min, double* fmim)
     {
        double lbound1,lbound2;
        unsigned int numIt = 0;
        int doGolden,i;
        int numAct=ndim;
        if( la<lc)
        {
            lbound1=la;
            lbound2=lc;
        }
        else
        {
            lbound1=lc;
            lbound2=la;
        }
        double lcurmin=lb;
        double lprevmin=lb;
        double lprevmin2=lb;
        double fcurmin=fb;
        double fprevmin=fb;
        double fprevmin2=fb;
        double dMoved=0;
        double dMoved2=0;
        double s= pow(10,-10);
        double lim;
        double par1,par2,lpartial,lnew,fnew;
        double GOLD=0.381966;
        double*  newD=(double*)_mol_malloc(numAct*sizeof(double) );
        double denom;
        double lmidpoint = 0.5*(lbound1+lbound2);
        double tempdMoved;   

        lim=tol*fabs(lcurmin)+s;
        while(fabs(lcurmin -lmidpoint)>(2*lim - 0.5*(lbound2-lbound1)) && numIt<maxtimes)
        { 
           lmidpoint = 0.5*(lbound1+lbound2);
           lim=tol*fabs(lcurmin)+s;
           doGolden=1;
           if(fabs(dMoved2)>lim) 
           {
              par1=(lcurmin-lprevmin)*(fcurmin-fprevmin2);
              par2=(lcurmin-lprevmin2)*(fcurmin-fprevmin);
              lpartial=((lcurmin-lprevmin)*par1-(lcurmin-lprevmin2)*par2);
              denom=2*(par1-par2);
    
              if(denom>0)
                 lpartial*=-1;
              denom=fabs(denom);
              tempdMoved=dMoved2;
              if(fabs(lpartial)<abs(0.5*denom*tempdMoved)
                 && lpartial>denom* (lbound1-lcurmin)
                 && lpartial<denom*(lbound2-lcurmin))
              {
                 dMoved=lpartial/denom;
                 lnew=lcurmin+dMoved;
                 if(lnew-lbound1<2*lim||lbound2-lnew<2*lim)
                 {
                    doGolden=0;
                    if(lmidpoint-lcurmin !=0)
                        dMoved=lim*(lmidpoint-lcurmin)/fabs(lmidpoint-lcurmin);
                    else
                        dMoved=lim;
                 }
              }
           }

           if( doGolden != 0)
           {
              if(lcurmin>=lmidpoint)
              dMoved2=lbound1-lcurmin;
              else
              dMoved2=lbound2-lcurmin;
              dMoved=GOLD*dMoved2;
           }
      
           if (fabs(dMoved)>lim)
               lnew=lcurmin+dMoved;
           else
           {
              if(dMoved != 0)
              lnew=lcurmin+lim*(dMoved/fabs(dMoved));
              else
              lnew = lcurmin;
           }
           for(i=0; i<numAct; i++)
               newD[i] = orig[i] + lnew*dir[i];
           egfun(ndim, newD, prms, &fnew, NULL);
           if(fnew<fcurmin)
           {
              if (lnew>=lcurmin)
                 lbound1=lcurmin;
              else
                 lbound2=lcurmin;
                 lprevmin2=lprevmin;
                 fprevmin2=fprevmin;
            
                 lprevmin=lcurmin;
                 fprevmin=fcurmin;
              
                 lcurmin=lnew;
                 fcurmin=fnew;
           }
           else
           {
              if(lnew<lcurmin)
                 lbound1=lnew;
              else
                 lbound2=lnew;
              if(fnew<fprevmin || lprevmin == fcurmin)
              {
                 lprevmin2=lprevmin;
                 fprevmin2=fprevmin;
                 lprevmin=lnew;
                 fprevmin=fnew;
              }
              else if(fnew<fprevmin2 || lprevmin2==lcurmin||lprevmin2==lprevmin)
              {
                 lprevmin2=lnew;
                 fprevmin2=fnew;
              }
           }
           dMoved2=dMoved;
           dMoved=lcurmin-lprevmin;
           numIt++;
        }       
        if(numIt==maxtimes)
             printf("MAXIMUM ITERATIONS REACHED. RETURNING CURRENT BEST\n");
        for(i=0; i<numAct; i++)
           min[i]=orig[i]+lcurmin*dir[i];
  
        *fmim=fcurmin;
        free(newD);
}

void dirbrent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, int maxtimes,
           double*  min, double* fmim, double*  grad)
{
       double lbound1,lbound2;
       int numIt=0;
       int doBisect, i, canUse1, canUse2;
       int numAct=ndim;

       if( la<lc)
       {
           lbound1=la;
           lbound2=lc;
       }
       else
       {
           lbound1=lc;
           lbound2=la;
       }
    
       double db=0;
        
       for(i=0; i<numAct; i++)
            db+=grad[i]*dir[i];

       double lcurmin=lb;
       double lprevmin=lb;
       double lprevmin2=lb;
       double dcurmin= db;
       double dprevmin= db;
       double dprevmin2= db;
       double fcurmin=fb;
       double fprevmin=fb;
       double fprevmin2=fb;
       double distMoved=0;
       double distMoved2=0;
       double s= pow(10,-10);
       double lim;
       double sec1=0,sec2=0,lnew,fnew,dnew;
       double*  newD=(double*)_mol_malloc(numAct*sizeof(double) );
       double lmidpoint;
       double tempdistMoved;

       for(numIt=0; numIt<maxtimes; numIt++ )
       {
           lmidpoint = 0.5*(lbound1+lbound2);
           lim=tol*fabs(lcurmin)+s;

           if    (fabs(lcurmin -lmidpoint)<=(2*lim - 0.5*(lbound2-lbound1))) 
           {
                for(i=0; i<numAct; i++)
                       min[i]=orig[i]+lcurmin*dir[i];
                *fmim=fcurmin;
                free(newD);
                return;
           }     

           doBisect=1;
           if(fabs(distMoved2)>lim) 
           {
               canUse1=0;
               canUse2=0;
               if (dprevmin != dcurmin)
               {
                   sec1=(lprevmin-lcurmin)*dcurmin/(dcurmin-dprevmin);
                   canUse1=1;
               }
   
               if (dprevmin2 != dcurmin)
               {
                   sec2=(lprevmin2-lcurmin)*dcurmin/(dcurmin-dprevmin2);
                   canUse2=1;
               }

               tempdistMoved=distMoved2;
               distMoved2=distMoved;

               if (canUse1 !=0)
               {
                   if(((lbound1-lcurmin-sec1) * (lcurmin+sec1-lbound2)<=0)||(dcurmin*sec1>0))
                             canUse1=0;
               }
        
               if (canUse2 !=0)
               {    
                   if(((lbound1-lcurmin-sec2) * (lcurmin+sec2-lbound2)<=0)||(dcurmin*sec2>0))
                                                 canUse2=0;
               }   

               if(canUse2 + canUse1 != 0)
               {
                   if(canUse1+canUse2 == 2)
                   {
                       if(fabs(sec1)<fabs(sec2))
                           distMoved=sec1;
                       else
                           distMoved=sec2;
                   }
                   else if(canUse1 != 0)
                       distMoved=sec1;
                   else if(canUse2 != 0)
                       distMoved=sec2;
                   doBisect=0;
                   if(fabs(distMoved)<= fabs(0.5*tempdistMoved))
                   {
                       lnew=lcurmin+distMoved;
                       if(((lnew-lbound1)<2*lim)||((lbound2-lnew)<2*lim))
                       {
                           if(lcurmin<lmidpoint)
                              distMoved=lim;
                           else 
                              distMoved=-lim;
                              doBisect=0;
                       }
                   }
                   else 
                      doBisect=1;                    
               }    
                     else doBisect=1;
           }
           if( doBisect != 0)
           {
              if(dcurmin>0)
                  distMoved2=lbound1-lcurmin;
              else
                  distMoved2=lbound2-lcurmin;
              distMoved=0.5*distMoved2;
           } 
              
           if (fabs(distMoved)>=lim)
           {
              lnew=lcurmin+distMoved;
              for(i=0; i<numAct; i++)
                  newD[i] = orig[i] + lnew*dir[i];
              egfun(ndim, newD, prms, &fnew, grad);
           }
           else
           {
              if(distMoved >=0)
                   lnew=lcurmin+lim;
              else
                   lnew = lcurmin-lim;
              for(i=0; i<numAct; i++)
                  newD[i] = orig[i] + lnew*dir[i];
              egfun(ndim, newD, prms, &fnew, grad);
             
              if (fnew>fcurmin)
              {
                 for(i=0; i<numAct; i++)
                     min[i]=orig[i]+lcurmin*dir[i];
                 *fmim=fcurmin;
                 free(newD);
                 return;
              }
           }  
           dnew=0;
           for(i=0; i<numAct; i++)
               dnew+=dir[i]*grad[i];      
           if((fnew<=fcurmin))
           {
               if (lnew>=lcurmin)
                   lbound1=lcurmin;
               else
                   lbound2=lcurmin;

               lprevmin2=lprevmin;
               fprevmin2=fprevmin;
               dprevmin2=dprevmin;

               lprevmin=lcurmin;
               fprevmin=fcurmin;
               dprevmin=dcurmin;

               lcurmin=lnew;
               fcurmin=fnew;
               dcurmin=dnew;
           }
           else
           {
              if(lnew<lcurmin)
                    lbound1=lnew;
              else
                    lbound2=lnew;
              if(fnew<=fprevmin || lprevmin ==lcurmin)
              {
                    lprevmin2=lprevmin;
                    fprevmin2=fprevmin;
                    dprevmin2=dprevmin;
                    lprevmin=lnew;
                    fprevmin=fnew;
                    dprevmin=dnew;
              }
              else if(fnew<=fprevmin2 || lprevmin2==lcurmin||lprevmin2==lprevmin)
              {
                    lprevmin2=lnew;
                    fprevmin2=fnew;
                    dprevmin2=dnew;
              }
           }
       }
       if(numIt==maxtimes)
             printf("MAXIMUM ITERATIONS REACHED. RETURNING CURRENT BEST\n");
       for(i=0; i<numAct; i++)
             min[i]=orig[i]+lcurmin*dir[i];
       *fmim=fcurmin;
       free(newD);
}


void powell(double* orig, double* directions, unsigned int maxIt, double tol,
            int ndim, void* prms,
            void (*egfun)(int , double* , void* , double* , double*),
            double* min, double* fmim) 
{
    int 
        j, k, //variables that will be used in for loops
        jmostdec=0,//direction of  the biggest  decrease withing a single directions set
        numAct=ndim;
    unsigned int i;
    double mostdec=0;
    double fbrac,fprev, fprev2=0, val;
    double    la=0,lb=0,lc=0;
    double    test;
    double* pcur=(double*)_mol_malloc(numAct*sizeof(double));//point of current minimum
    double* prev=(double*)_mol_malloc(numAct*sizeof(double));//previous minimum
    double* newdir=(double*)_mol_malloc(numAct*sizeof(double));//average direction moved in one iteration

    
    for(j=0; j<numAct;j++)
    {
        prev[j]=orig[j];
        pcur[j]=orig[j];
    }
    egfun(ndim, pcur, prms, &val, NULL);

    for(i=0; i<maxIt; i++)//loop through once for each direction set(break out if close enough)
    {
        fprev=val;
        mostdec=0;
        jmostdec=0;
        for(j=0; j<numAct ;j++)//loop through once for each direction
        {
            for(k=0;k<numAct;k++)
               newdir[k]=directions[numAct*j+k];//get direction
            fprev2=val;
            bracket(pcur, newdir, 1, ndim, prms, egfun, &fbrac,  &la, &lb, &lc );
            brent(pcur,  newdir,  
                  fbrac, la, lb, lc,
                  ndim, prms, egfun,
                  pow(10,-5), 100,    
                  pcur, &val);
            if(fabs(fprev2-val)>mostdec)//get direction of greatest decrease and greatest decrease
            {
                mostdec=fabs(fprev2-val);
                jmostdec=j;
            }
        }

        if( 2*fabs(fprev-val)<=tol*(fabs(fprev)+fabs(val)))//if you reach the desired tolerance end powell
        { 
            for(j=0; j<numAct; j++)
                min[j]=pcur[j];
            
            *fmim=val;
            free(pcur);
            free(prev);
            free(newdir);
            return;
        }


        for(j=0; j<numAct; j++)
            newdir[j]=2*pcur[j]-prev[j];//get new direction

        egfun(ndim, newdir, prms, &fprev2, NULL);

        for(j=0; j<numAct; j++)
        {
            newdir[j]=pcur[j]-prev[j];
            prev[j]=pcur[j];
        }
        
        if(fprev2<fprev)
        {
            test=2.*(fprev-2.*val+fprev2)*pow((fprev-val-mostdec),2)-mostdec*pow((fprev-fprev2),2);
            if(test<0)
            {//if it makes sence to the direction set
                bracket(pcur, newdir, 1, ndim, prms, egfun,
                        &fbrac,  &la, &lb, &lc );//minimize in new direction
                brent(pcur, newdir,
                      fbrac, la, lb, lc,
                      ndim, prms, egfun,
                      pow(10,-5), 100,
                      pcur, &val);
                for(k=0; k<numAct; k++)
                {
                      directions[jmostdec*numAct+k]=directions[(numAct-1)*(numAct)+k];
                      directions[(numAct-1)*(numAct)+k]= newdir[k];
                }
            }
        }
    }         
    printf("warning maximum number of iterations reached. using current best\n");
        
    for(j=0; j<numAct; j++)
                min[j]=pcur[j];

    *fmim=val;
    free(pcur);
    free(prev);
    free(newdir);
}

void dirMin(double* orig,  unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim)
{
     int numAct = ndim;
     int i,j;
     unsigned int it;

     double* dir=(double*)malloc(numAct*sizeof(double)); 
     double* interm=(double*)malloc(numAct*sizeof(double)); 
     double* temp=(double*)malloc(numAct*sizeof(double));
     double* curmin=(double*)malloc(numAct*sizeof(double));
     double* grad = (double*)malloc(numAct*sizeof(double));
     double* mvec = (double*)malloc(numAct*sizeof(double));

     for(i=0; i<numAct; i++)
            curmin[i]=orig[i];
        
     double val,t1,t2,t3,la,lb,lc,fbrac ;
     double s = pow(10,-5);
     double fprev;
     egfun(ndim, orig, prms, &fprev, grad);
     double dnorm=0, maxnorm=1.5;

     for(i=0; i<numAct; i++)
     {
        interm[i] =-grad[i] ;
        temp[i]=interm[i];
        dir[i]=interm[i];
     }

     for(it=0; it<maxIt; it++)
     {
        dnorm=0;
        for(i=0; i<numAct; i++)
            dnorm+= dir[i]*dir[i];
        dnorm=sqrt(dnorm);
        if(dnorm>maxnorm)
            for(i=0; i<numAct; i++)
                dir[i]*=maxnorm/dnorm;
                 
        bracket(curmin,  dir, 1,
                ndim, prms, egfun,
                &fbrac,  &la, &lb, &lc);
        
        dirbrent(curmin,dir,  fbrac, la, lb, lc,
                 ndim, prms, egfun,
                 pow(10,-10), 100,
                 mvec, &val, grad);

        if(2*fabs(val-fprev)<= tol*(fabs(fprev)+fabs(val)+s))
        {
           for(j=0; j<numAct; j++)
               min[j]=mvec[j];
           *fmim=val;
           free(dir);
           free(mvec);
           free(interm);
           free(temp);
           return;
        }
    
        for(j=0; j<numAct; j++)
        {
           curmin[j]=mvec[j];
           dir[j]=grad[j];
        }    
        fprev =val;
        t1=0;
        t2=0;
        for(j=0; j<numAct; j++)
        {
            t1+= pow(interm[j],2);
            t2+= (dir[j]+interm[j])*dir[j] ;
        }
      
        if(t1 ==0)
        {
            for(j=0; j<numAct; j++)
                min[j]=curmin[j];
            *fmim=val;
            free(dir);
            free(interm);
            free(temp);
            free(mvec);
            return;
        }
        
        t3=t2/t1;

        for(j=0; j<numAct; j++)
        {
            interm[j] =-dir[j];
            temp[j]=interm[j]+t3*temp[j];
            dir[j]=temp[j];
        }
     }
     printf("maximum num of iterations reached displaying current min\n");
     for(j=0; j<numAct; j++)
                min[j]=curmin[j];
     *fmim=val;
     free(dir);
     free(interm);
     free(temp);
     free(mvec);
     return;
}

void limin(double* orig, double* dir, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim)
{
        double fbrac, la, lb, lc;
        bracket(orig,  dir, 1,
                ndim, prms, egfun,
                &fbrac,  &la, &lb, &lc);
        brent(orig,  dir,
                  fbrac, la, lb, lc,
                  ndim, prms, egfun,
                  tol, maxIt,
                  min, fmim);
}


//Powell/dirPowell is above
 

double
max3 (double x, double y, double z)
{
  return x < y ? (y < z ? z : y) : (x < z ? z : x);
}

double
max (double x, double y)
{
  return x < y ? y : x;
}

double
min (double x, double y)
{
  return x < y ? x : y;
}

double
sqr (double x)
{
  return x * x;
}
void

linestep (double *stx, double *fx, double *dx, double *sty, double *fy,
	  double *dy, double *stp, double fp, double dp, short *brackt,
	  double stpmin, double stpmax, int *info)
{

  double xsafe = 0.001;
  short bound;
  double lsgamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;

  info[0] = 0;

  if ((brackt[0]
       && (stp[0] <= min (stx[0], sty[0]) || stp[0] >= max (stx[0], sty[0])))
      || dx[0] * (stp[0] - stx[0]) >= 0.0 || stpmax < stpmin)
    return;

  // Determine if the derivatives have opposite sign.

  sgnd = dp * (dx[0] / fabs (dx[0]));

  if (fp > fx[0])
    {
      // First case. A higher function value.
      // The minimum is bracketed. If the cubic step is closer
      // to stx than the quadratic step, the cubic step is taken,
      // else the average of the cubic and quadratic steps is taken.

      info[0] = 1;
      bound = 1;
      theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
      s = max3 (fabs (theta), fabs (dx[0]), fabs (dp));
      lsgamma = s * sqrt (sqr (theta / s) - (dx[0] / s) * (dp / s));
      if (stp[0] < stx[0])
	lsgamma = -lsgamma;
      p = (lsgamma - dx[0]) + theta;
      q = ((lsgamma - dx[0]) + lsgamma) + dp;
      r = p / q;
      stpc = stx[0] + r * (stp[0] - stx[0]);
      stpq =
	stx[0] +
	((dx[0] / ((fx[0] - fp) / (stp[0] - stx[0]) + dx[0])) / 2) * (stp[0] -
								      stx[0]);
      if (fabs ((stpc - stx[0])) < fabs ((stpq - stx[0])))
	{
	  stpf = stpc;
	}
      else
	{
	  stpf = stpc + (stpq - stpc) / 2;
	}
	   
        if ( stp[0] > stx[0])
                stpf = max (stx[0]+xsafe*(stp[0]-stx[0]),stpf);
        else
                stpf = min (stx[0]+xsafe*(stp[0]-stx[0]),stpf);

      brackt[0] = 1;
    }
  else if (sgnd < 0.0)
    {
      // Second case. A lower function value and derivatives of
      // opposite sign. The minimum is bracketed. If the cubic
      // step is closer to stx than the quadratic (secant) step,
      // the cubic step is taken, else the quadratic step is taken.

      info[0] = 2;
      bound = 0;
      theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
      s = max3 (fabs (theta), fabs (dx[0]), fabs (dp));
      lsgamma = s * sqrt (sqr (theta / s) - (dx[0] / s) * (dp / s));
      if (stp[0] > stx[0])
	lsgamma = -lsgamma;
      p = (lsgamma - dp) + theta;
      q = ((lsgamma - dp) + lsgamma) + dx[0];
      r = p / q;
      stpc = stp[0] + r * (stx[0] - stp[0]);
      stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
      if (fabs (stpc - stp[0]) > fabs (stpq - stp[0]))
	{
	  stpf = stpc;
	}
      else
	{
	  stpf = stpq;
	}
      brackt[0] = 1;
    }
  else if (fabs (dp) < fabs (dx[0]))
    {
      // Third case. A lower function value, derivatives of the
      // same sign, and the magnitude of the derivative decreases.
      // The cubic step is only used if the cubic tends to infinity
      // in the direction of the step or if the minimum of the cubic
      // is beyond stp. Otherwise the cubic step is defined to be
      // either stpmin or stpmax. The quadratic (secant) step is also
      // computed and if the minimum is bracketed then the the step
      // closest to stx is taken, else the step farthest away is taken.

      info[0] = 3;
      bound = 1;
      theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
      s = max3 (fabs (theta), fabs (dx[0]), fabs (dp));
      lsgamma = s * sqrt (max (0, sqr (theta / s) - (dx[0] / s) * (dp / s)));
      if (stp[0] > stx[0])
	lsgamma = -lsgamma;
      p = (lsgamma - dp) + theta;
      q = (lsgamma + (dx[0] - dp)) + lsgamma;
      r = p / q;
      if (r < 0.0 && lsgamma != 0.0)
	{
	  stpc = stp[0] + r * (stx[0] - stp[0]);
	}
      else if (stp[0] > stx[0])
	{
	  stpc = stpmax;
	}
      else
	{
	  stpc = stpmin;
	}
      stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
      if (brackt[0])
	{
	  if (fabs (stp[0] - stpc) < fabs (stp[0] - stpq))
	    {
	      stpf = stpc;
	    }
	  else
	    {
	      stpf = stpq;
	    }
	}
      else
	{
	  if (fabs (stp[0] - stpc) > fabs (stp[0] - stpq))
	    {
	      stpf = stpc;
	    }
	  else
	    {
	      stpf = stpq;
	    }
	}
    }
  else
    {
      // Fourth case. A lower function value, derivatives of the
      // same sign, and the magnitude of the derivative does
      // not decrease. If the minimum is not bracketed, the step
      // is either stpmin or stpmax, else the cubic step is taken.

      info[0] = 4;
      bound = 0;
      if (brackt[0])
	{
	  theta = 3 * (fp - fy[0]) / (sty[0] - stp[0]) + dy[0] + dp;
	  s = max3 (fabs (theta), fabs (dy[0]), fabs (dp));
	  lsgamma = s * sqrt (sqr (theta / s) - (dy[0] / s) * (dp / s));
	  if (stp[0] > sty[0])
	    lsgamma = -lsgamma;
	  p = (lsgamma - dp) + theta;
	  q = ((lsgamma - dp) + lsgamma) + dy[0];
	  r = p / q;
	  stpc = stp[0] + r * (sty[0] - stp[0]);
	  stpf = stpc;
	}
      else if (stp[0] > stx[0])
	{
	  stpf = stpmax;
	}
      else
	{
	  stpf = stpmin;
	}
    }

  // Update the interval of uncertainty. This update does not
  // depend on the new step or the case analysis above.

  if (fp > fx[0])
    {
      sty[0] = stp[0];
      fy[0] = fp;
      dy[0] = dp;
    }
  else
    {
      if (sgnd < 0.0)
	{
	  sty[0] = stx[0];
	  fy[0] = fx[0];
	  dy[0] = dx[0];
	}
      stx[0] = stp[0];
      fx[0] = fp;
      dx[0] = dp;
    }

  // Compute the new step and safeguard it.

  stpf = min (stpmax, stpf);
  stpf = max (stpmin, stpf);
  stp[0] = stpf;

  if (brackt[0] && bound)
    {
      if (sty[0] > stx[0])
	{
	  stp[0] = min (stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
	}
      else
	{
	  stp[0] = max (stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
	}
    }

  return;
}



double gtol = 0.9;
int infoc[1], j = 0;
double stpmin = 1e-20;
double stpmax = 1e20;
short brackt[1], stage1 = 0;
  double dg = 0, dgm = 0, dginit = 0, dgtest =
    0, dgx[1], dgxm[1], dgy[1], dgym[1], finit = 0, ftest1 = 0, fm =
    0, fx[1], fxm[1], fy[1], fym[1], p5 = 0, p66 = 0, stx[1], sty[1], stmin =
    0, stmax = 0, width = 0, width1 = 0, xtrapf = 0;
void
linesearch (int n, double *x, double f, double *g, double *s, int is0,
	    double *stp, double ftol, double xtol, int maxfev, int *info,
	    int *nfev, double *wa)
{
    /* double stpmin = 1e-20;
  double stpmax = 1e20;
  short brackt[1], stage1 = 0;
  double dg = 0, dgm = 0, dginit = 0, dgtest =
    0, dgx[1], dgxm[1], dgy[1], dgym[1], finit = 0, ftest1 = 0, fm =
    0, fx[1], fxm[1], fy[1], fym[1], p5 = 0, p66 = 0, stx[1], sty[1], stmin =
    0, stmax = 0, width = 0, width1 = 0, xtrapf = 0;*/
  p5 = 0.5;
  p66 = 0.66;
  xtrapf = 4;

  if (info[0] != -1)
    {
      infoc[0] = 1;
      if (n <= 0 || stp[0] <= 0 || ftol < 0 || gtol < 0 || xtol < 0
	  || stpmin < 0 || stpmax < stpmin || maxfev <= 0)
	return;
      dginit = 0;
      for (j = 1; j <= n; j += 1)
	{
	  dginit = dginit + g[j - 1] * s[is0 + j - 1];
	}
      if (dginit >= 0)
	{
	  printf ("Bad search direction.\n");
	  return;
	}
      brackt[0] = 0;
      stage1 = 1;
      nfev[0] = 0;
      finit = f;
      dgtest = ftol * dginit;
      width = stpmax - stpmin;
      width1 = width / p5;
      for (j = 1; j <= n; j += 1)
	{
	  wa[j - 1] = x[j - 1];
	}

      // function, and derivative at the other endpoint of
      // the interval of uncertainty.
      // The variables stp, f, dg contain the values of the step,
      // function, and derivative at the current step.

      stx[0] = 0;
      fx[0] = finit;
      dgx[0] = dginit;
      sty[0] = 0;
      fy[0] = finit;
      dgy[0] = dginit;
    }

  while (1)
    {
      if (info[0] != -1)
	{
	  // Set the minimum and maximum steps to correspond
	  // to the present interval of uncertainty.

	  if (brackt[0])
	    {
	      stmin = min (stx[0], sty[0]);
	      stmax = max (stx[0], sty[0]);
	    }
	  else
	    {
	      stmin = stx[0];
	      stmax = stp[0] + xtrapf * (stp[0] - stx[0]);
	    }

	  // Force the step to be within the bounds stpmax and stpmin.

	  stp[0] = max (stp[0], stpmin);
	  stp[0] = min (stp[0], stpmax);

	  // If an unusual termination is to occur then let
	  // stp be the lowest point obtained so far.

	  if ((brackt[0] && (stp[0] <= stmin || stp[0] >= stmax))
	      || nfev[0] >= maxfev - 1 || infoc[0] == 0 || (brackt[0]
							    && stmax -
							    stmin <=
							    xtol * stmax))
	    stp[0] = stx[0];
	  // Evaluate the function and gradient at stp
	  // and compute the directional derivative.
	  // We return to main program to obtain F and G.
	  for (j = 1; j <= n; j += 1)
	    {
	      x[j - 1] = wa[j - 1] + stp[0] * s[is0 + j - 1];
	    }
	  info[0] = -1;
	  return;
	}
      info[0] = 0;
      nfev[0] = nfev[0] + 1;
      dg = 0;
      for (j = 1; j <= n; j += 1)
	{
	  dg = dg + g[j - 1] * s[is0 + j - 1];
	}
      ftest1 = finit + stp[0] * dgtest;
      // Test for convergence.
      if ((brackt[0] && (stp[0] <= stmin || stp[0] >= stmax))
	  || infoc[0] == 0)
	info[0] = 6;
      if (stp[0] == stpmax && f <= ftest1 && dg <= dgtest)
	info[0] = 5;
      if (stp[0] == stpmin && (f > ftest1 || dg >= dgtest))
	info[0] = 4;
      if (nfev[0] >= maxfev)
	info[0] = 3;
      if (brackt[0] && stmax - stmin <= xtol * stmax)
	info[0] = 2;
      if (f <= ftest1 && fabs (dg) <= gtol * (-dginit))
	info[0] = 1;
      // Check for termination.
      if (info[0] != 0)
	return;
      // In the first stage we seek a step for which the modified
      // function has a nonpositive value and nonnegative derivative.
      if (stage1 && f <= ftest1 && dg >= min (ftol, gtol) * dginit)
	stage1 = 0;
      // A modified function is used to predict the step only if
      // we have not obtained a step for which the modified
      // function has a nonpositive function value and nonnegative
      // derivative, and if a lower function value has been
      // obtained but the decrease is not sufficient.
      if (stage1 && f <= fx[0] && f > ftest1)
	{
	  // Define the modified function and derivative values.
	  fm = f - stp[0] * dgtest;
	  fxm[0] = fx[0] - stx[0] * dgtest;
	  fym[0] = fy[0] - sty[0] * dgtest;
	  dgm = dg - dgtest;
	  dgxm[0] = dgx[0] - dgtest;
	  dgym[0] = dgy[0] - dgtest;
	  // Call cstep to update the interval of uncertainty
	  // and to compute the new step.
	  linestep (stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt,
		    stmin, stmax, infoc);
	  // Reset the function and gradient values for f.
	  fx[0] = fxm[0] + stx[0] * dgtest;
	  fy[0] = fym[0] + sty[0] * dgtest;
	  dgx[0] = dgxm[0] + dgtest;
	  dgy[0] = dgym[0] + dgtest;
	}
      else
	{
	  // Call mcstep to update the interval of uncertainty
	  // and to compute the new step.
	  linestep (stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin,
		    stmax, infoc);
	}

      // Force a sufficient decrease in the size of the
      // interval of uncertainty.

      if (brackt[0])
	{
	  if (fabs (sty[0] - stx[0]) >= p66 * width1)
	    stp[0] = stx[0] + p5 * (sty[0] - stx[0]);
	  width1 = width;
	  width = fabs (sty[0] - stx[0]);
	}
    }
}

	/** Print debugging and status messages for <code>lbfgs</code>.
	  * Depending on the parameter <code>iprint</code>, this can include 
	  * number of function evaluations, current function value, etc.
	  * The messages are output to <code>System.err</code>.
	  *
	  * @param iprint Specifies output generated by <code>lbfgs</code>.<p>
	  *		<code>iprint[0]</code> specifies the frequency of the output:
	  *		<ul>
	  *		<li> <code>iprint[0] &lt; 0</code>: no output is generated,
	  *		<li> <code>iprint[0] = 0</code>: output only at first and last iteration,
	  *		<li> <code>iprint[0] &gt; 0</code>: output every <code>iprint[0]</code> iterations.
	  *		</ul><p>
	  *
	  *		<code>iprint[1]</code> specifies the type of output generated:
	  *		<ul>
	  *		<li> <code>iprint[1] = 0</code>: iteration count, number of function 
	  *			evaluations, function value, norm of the gradient, and steplength,
	  *		<li> <code>iprint[1] = 1</code>: same as <code>iprint[1]=0</code>, plus vector of
	  *			variables and  gradient vector at the initial point,
	  *		<li> <code>iprint[1] = 2</code>: same as <code>iprint[1]=1</code>, plus vector of
	  *			variables,
	  *		<li> <code>iprint[1] = 3</code>: same as <code>iprint[1]=2</code>, plus gradient vector.
	  *		</ul>
	  * @param iter Number of iterations so far.
	  * @param nfun Number of function evaluations so far.
	  * @param gnorm Norm of gradient at current solution <code>x</code>.
	  * @param n Number of free parameters.
	  * @param m Number of corrections kept.
	  * @param x Current solution.
	  * @param f Function value at current solution.
	  * @param g Gradient at current solution <code>x</code>.
	  * @param stp Current stepsize.
	  * @param finish Whether this method should print the ``we're done'' message.
	  */
void

lb1 (int *iprint, int *iter, int *nfun, double gnorm, int n, int m, double *x,
     double f, double *g, double *stp, short finish)
{
  int i;

  if (iter[0] == 0)
    {
      printf ("*************************************************\n");
      printf
	("  n = %d   number of corrections = %d\n       initial values\n", n,
	 m);
      printf (" f =  %.3f   gnorm =  %.3f\n", f, gnorm);
      if (iprint[1] >= 1)
	{
	  printf (" vector x =\n");
	  for (i = 0; i <= n - 1; i++)
	    printf ("  %.3f", x[i]);
	  printf ("\n");

	  printf (" gradient vector g =");
	  for (i = 0; i <= n - 1; i++)
	    printf ("  %.3f", x[i]);
	  printf ("\n");
	}
      printf ("*************************************************\n");
      printf ("\ti\tnfn\tfunc\tgnorm\tsteplength\n");
    }
  else
    {
      if ((iprint[0] == 0) && (iter[0] != 1 && !finish))
	return;
      if (iprint[0] != 0)
	{
	  if ((iter[0] - 1) % iprint[0] == 0 || finish)
	    {
	      if (iprint[1] > 1 && iter[0] > 1)
		printf ("\ti\tnfn\tfunc\tgnorm\tsteplength\n");
	      printf ("\t %d\t %d\t %.12f\t %.12f\t%.12f\n", iter[0], nfun[0], f,
		      gnorm, stp[0]);
	    }
	  else
	    {
	      return;
	    }
	}
      else
	{
	  if (iprint[1] > 1 && finish)
	    printf ("\ti\tnfn\tfunc\tgnorm\tsteplength\n");
	  printf ("\t %d\t %d\t %.12f\t %.12f\t%.12f\n", iter[0], nfun[0], f,
		  gnorm, stp[0]);

	}
      if (iprint[1] == 2 || iprint[1] == 3)
	{
	  if (finish)
	    {
	      printf (" final point x =\n");
	    }
	  else
	    {
	      printf (" vector x =  \n");
	    }
	  for (i = 0; i <= n - 1; i++)
	    printf ("%.3f  ", x[i - 1]);
	  printf ("\n");
	  if (iprint[1] == 3)
	    {
	      printf (" gradient vector g =");
	      for (i = 0; i <= n - 1; i++)
		printf ("  %.3f", g[i]);
	      printf ("\n");
	    }
	}
      if (finish)
	printf
	  (" The minimization terminated without detecting errors. iflag = 0\n");
    }
  return;
}


double  stp1 = 0, ftol = 0, stp[1], ys = 0, yy = 0, sq = 0, yr =
    0, beta = 0, xnorm = 0;
  int point = 0, ispt = 0, iypt = 0, maxfev = 0,info[1], bound = 0, npt =
    0, cp = 0, nfev[1], inmc = 0, iycn = 0, iscn = 0;
 
/* worker array w needs to be initialized*/
void
lbfgs (int n, int m, double *x, double f, double *g, short diagco,
       double *diag, int *iprint, double eps, double xtol, double *w,
       int *iflag, int *iter, int *nfun)
{
/*  double gnorm = 0, stp1 = 0, ftol = 0, stp[1], ys = 0, yy = 0, sq = 0, yr =
    0, beta = 0, xnorm = 0;
  int point = 0, ispt = 0, iypt = 0, maxfev = 0,info[1], bound = 0, npt =
  0, cp = 0, i = 0, nfev[1], inmc = 0, iycn = 0, iscn = 0;*/
    short finish = 0;
  short execute_entire_while_loop = 0;
  int i=0;
  double gnorm=0;
  
  if (iflag[0] == 0)
    {
      // Initialize.

      iter[0] = 0;

      if (n <= 0 || m <= 0)
	{
	  iflag[0] = -3;
	  printf
	    ("Improper input parameters  (n or m are not positive.). Error %d\n",
	     iflag[0]);
	  return;
	}

      if (gtol <= 0.0001)
	{
	  printf
	    ("LBFGS.lbfgs: gtol is less than or equal to 0.0001. It has been reset to 0.9.\n");
	  gtol = 0.9;
	}

      nfun[0] = 1;
      point = 0;
      finish = 0;

      if (diagco)
	{
	  for (i = 1; i <= n; i += 1)
	    {
	      if (diag[i - 1] <= 0)
		{
		  iflag[0] = -2;
		  printf
		    ("The %d-th diagonal element of the inverse hessian approximation is not positive. Error %d\n",
		     i, iflag[0]);
		  return;
		}
	    }
	}
      else
	{
	  for (i = 1; i <= n; i += 1)
	    {
	      diag[i - 1] = 1;
	    }
	}
      ispt = n + 2 * m;
      iypt = ispt + n * m;

      for (i = 1; i <= n; i += 1)
	{
	  w[ispt + i - 1] = -g[i - 1] * diag[i - 1];
	}

      gnorm = sqrt (ddot (n, g, 0, 1, g, 0, 1));
      stp1 = 1 / gnorm;
      ftol = 0.0001;
      maxfev = 20;

      if (iprint[1 - 1] >= 0)
	lb1 (iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);

      execute_entire_while_loop = 1;
    }

  while (1)
    {
      if (execute_entire_while_loop)
	{
	  iter[0] = iter[0] + 1;
	  info[0] = 0;
	  bound = iter[0] - 1;
	  if (iter[0] != 1)
	    {
	      if (iter[0] > m)
		bound = m;
	      ys = ddot (n, w, iypt + npt, 1, w, ispt + npt, 1);
	      if (!diagco)
		{
		  yy = ddot (n, w, iypt + npt, 1, w, iypt + npt, 1);

		  for (i = 1; i <= n; i += 1)
		    {
		      diag[i - 1] = ys / yy;
		    }
		}
	      else
		{
		  iflag[0] = 2;
		  return;
		}
	    }
	}

      if (execute_entire_while_loop || iflag[0] == 2)
	{
	  if (iter[0] != 1)
	    {
	      if (diagco)
		{
		  for (i = 1; i <= n; i += 1)
		    {
		      if (diag[i - 1] <= 0)
			{
			  iflag[0] = -2;
			  printf
			    ("The %d-th diagonal element of the inverse hessian approximation is not positive. Error %d\n",
			     i, iflag[0]);
			  return;
			}
		    }
		}
	      cp = point;
	      if (point == 0)
		cp = m;
	      w[n + cp - 1] = 1 / ys;

	      for (i = 1; i <= n; i += 1)
		{
		  w[i - 1] = -g[i - 1];
		}

	      cp = point;

	      for (i = 1; i <= bound; i += 1)
		{
		  cp = cp - 1;
		  if (cp == -1)
		    cp = m - 1;
		  sq = ddot (n, w, ispt + cp * n, 1, w, 0, 1);
		  inmc = n + m + cp + 1;
		  iycn = iypt + cp * n;
		  w[inmc - 1] = w[n + cp + 1 - 1] * sq;
		  daxpy (n, -w[inmc - 1], w, iycn, 1, w, 0, 1);
		}

	      for (i = 1; i <= n; i += 1)
		{
		  w[i - 1] = diag[i - 1] * w[i - 1];
		}

	      for (i = 1; i <= bound; i += 1)
		{
		  yr = ddot (n, w, iypt + cp * n, 1, w, 0, 1);
		  beta = w[n + cp + 1 - 1] * yr;
		  inmc = n + m + cp + 1;
		  beta = w[inmc - 1] - beta;
		  iscn = ispt + cp * n;
		  daxpy (n, beta, w, iscn, 1, w, 0, 1);
		  cp = cp + 1;
		  if (cp == m)
		    cp = 0;
		}

	      for (i = 1; i <= n; i += 1)
		{
		  w[ispt + point * n + i - 1] = w[i - 1];
		}
	    }

	  nfev[0] = 0;
	  stp[0] = 1;
	  if (iter[0] == 1)
	    stp[0] = stp1;

	  for (i = 1; i <= n; i += 1)
	    {
	      w[i - 1] = g[i - 1];
	    }
	}

      linesearch (n, x, f, g, w, ispt + point * n, stp, ftol, xtol, maxfev,
		  info, nfev, diag);

      if (info[0] == -1)
	{
	  iflag[0] = 1;
	  return;
	}

      if (info[0] != 1)
	{
	  iflag[0] = -1;
	  if (iprint[1 - 1] >= 0){
	  printf
	    ("Line search failed. See documentation of routine mcsrch. Error return of line search: info = %d Possible causes: function or gradient are incorrect, or incorrect tolerances.\n",
	     info[0]);
	  }
	  return;
	}

      nfun[0] = nfun[0] + nfev[0];
      npt = point * n;

      for (i = 1; i <= n; i += 1)
	{
	  w[ispt + npt + i - 1] = stp[0] * w[ispt + npt + i - 1];
	  w[iypt + npt + i - 1] = g[i - 1] - w[i - 1];
	}

      point = point + 1;
      if (point == m)
	point = 0;

      gnorm = sqrt (ddot (n, g, 0, 1, g, 0, 1));
      xnorm = sqrt (ddot (n, x, 0, 1, x, 0, 1));
      xnorm = max (1.0, xnorm);

      if (gnorm / xnorm <= eps)
	finish = 1;

      if (iprint[1 - 1] >= 0)
	lb1 (iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
      if (finish)
	{
	  iflag[0] = 0;
	  return;
	}

      execute_entire_while_loop = 1;	// from now on, execute whole loop
    }
}

	/** Compute the sum of a vector times a scalara plus another vector.
	  * Adapted from the subroutine <code>daxpy</code> in <code>lbfgs.f</code>.
	  * There could well be faster ways to carry out this operation; this
	  * code is a straight translation from the Fortran.
	  */
void
daxpy (int n, double da, double *dx, int ix0, int incx, double *dy, int iy0,
       int incy)
{
  int i, ix, iy, m, mp1;

  if (n <= 0)
    return;

  if (da == 0)
    return;

  if (!(incx == 1 && incy == 1))
    {
      ix = 1;
      iy = 1;

      if (incx < 0)
	ix = (-n + 1) * incx + 1;
      if (incy < 0)
	iy = (-n + 1) * incy + 1;

      for (i = 1; i <= n; i += 1)
	{
	  dy[iy0 + iy - 1] = dy[iy0 + iy - 1] + da * dx[ix0 + ix - 1];
	  ix = ix + incx;
	  iy = iy + incy;
	}

      return;
    }

  m = n % 4;
  if (m != 0)
    {
      for (i = 1; i <= m; i += 1)
	{
	  dy[iy0 + i - 1] = dy[iy0 + i - 1] + da * dx[ix0 + i - 1];
	}

      if (n < 4)
	return;
    }

  mp1 = m + 1;
  for (i = mp1; i <= n; i += 4)
    {
      dy[iy0 + i - 1] = dy[iy0 + i - 1] + da * dx[ix0 + i - 1];
      dy[iy0 + i + 1 - 1] = dy[iy0 + i + 1 - 1] + da * dx[ix0 + i + 1 - 1];
      dy[iy0 + i + 2 - 1] = dy[iy0 + i + 2 - 1] + da * dx[ix0 + i + 2 - 1];
      dy[iy0 + i + 3 - 1] = dy[iy0 + i + 3 - 1] + da * dx[ix0 + i + 3 - 1];
    }
  return;
}

	/** Compute the dot product of two vectors.
	  * Adapted from the subroutine <code>ddot</code> in <code>lbfgs.f</code>.
	  * There could well be faster ways to carry out this operation; this
	  * code is a straight translation from the Fortran.
	  */
double
ddot (int n, double *dx, int ix0, int incx, double *dy, int iy0, int incy)
{
  double dtemp;
  int i, ix, iy, m, mp1;

  dtemp = 0;

  if (n <= 0)
    return 0;

  if (!(incx == 1 && incy == 1))
    {
      ix = 1;
      iy = 1;
      if (incx < 0)
	ix = (-n + 1) * incx + 1;
      if (incy < 0)
	iy = (-n + 1) * incy + 1;
      for (i = 1; i <= n; i += 1)
	{
	  dtemp = dtemp + dx[ix0 + ix - 1] * dy[iy0 + iy - 1];
	  ix = ix + incx;
	  iy = iy + incy;
	}
      return dtemp;
    }

  m = n % 5;
  if (m != 0)
    {
      for (i = 1; i <= m; i += 1)
	{
	  dtemp = dtemp + dx[ix0 + i - 1] * dy[iy0 + i - 1];
	}
      if (n < 5)
	return dtemp;
    }

  mp1 = m + 1;
  for (i = mp1; i <= n; i += 5)
    {
      dtemp =
	dtemp + dx[ix0 + i - 1] * dy[iy0 + i - 1] + dx[ix0 + i + 1 -1] * dy[iy0 + i + 1 -1] + dx[ix0 + i + 2 - 1] *
	dy[iy0 + i + 2 - 1] + dx[ix0 + i + 3 - 1] * dy[iy0 + i + 3 - 1] +
	dx[ix0 + i + 4 - 1] * dy[iy0 + i + 4 - 1];
    }

  return dtemp;
}

lbfgs_t *
lbfgs_create (int n, int m, double eps)
{
  lbfgs_t *opt = (lbfgs_t *) malloc (sizeof (lbfgs_t));

  if (!opt)
    return 0;
  opt->w = (double *) malloc (sizeof (double) * (n * (2 * m + 1) + 2 * m));
  if (!opt->w)
    {
      free (opt);
      return 0;
    }
  opt->diag = (double *) malloc (sizeof (double) * n);
  if (!opt->diag)
    {
      free (opt->w);
      free (opt);
      return 0;
    }

  opt->n = n;
  opt->m = m;
  opt->eps = eps;
  opt->xtol = DBL_EPSILON;
  opt->diagco = 0;		/* by default we do not provide diagonal matrices */
  opt->iflag = 0;
  opt->iprint[0] = -1;		/* by default print nothing                       */
  opt->iprint[1] = 0;
  opt->niter = 0;
  opt->nfuns = 0;

  return opt;
}

/* free all the memory used by the optimizer */
void
lbfgs_destroy (lbfgs_t * opt)
{
  free (opt->diag);
  free (opt->w);
  free (opt);
}


/* x is the n-dimension variable
 * f is the current value of objective function
 * g is the n-dimension gradient of f
 *
 * return value:
 *       = 0: success
 *       < 0: some error occur, see comment in lbfgs.f for detail
 *       = 1: user must evaluate F and G
 *       = 2: user must provide the diagonal matrix Hk0.
 */
int
lbfgs_run (lbfgs_t * opt, double *x, double f, double *g)
{
  lbfgs (opt->n, opt->m, x, f, g, opt->diagco, opt->diag,
	 opt->iprint, opt->eps, opt->xtol, opt->w, &opt->iflag,
	 &opt->niter, &opt->nfuns);
//    printf("%d\n",opt->iflag);
  return opt->iflag;
}

//LBFGS minimizer main loop
void lbfgsMin(double* orig, unsigned int maxIt, double tol, int ndim,  void* prms,void (*egfun)(int , double* , void* , double* , double*),double* minv, double* fmim ){
    lbfgs_t* opt;
    int i;
    unsigned int nsteps;
    double* grad = (double*) _mol_malloc(ndim*sizeof(double));
    opt = lbfgs_create(ndim, 5, tol);
    opt->iprint[0] = 0;
    opt->iprint[1] = 0;
    opt->diagco = 0;
    nsteps = 0;
    for(i=0; i<ndim; i++)
            minv[i]=orig[i];
    while (1) {
	 egfun(ndim, minv, prms, fmim, grad); 
         if ((lbfgs_run(opt, minv, (*fmim), grad) <= 0)|| (++nsteps >maxIt))
	               break;
	    }
    lbfgs_destroy(opt);
    free(grad);
}



/* Minimizes atomgroup with specified parameter set
 minimization types; 0 = LBFGS, 1- Conjugate gradients, 2 -Powell 
*/
void minimize_ag(mol_min_method min_type, unsigned int maxIt, double tol,struct atomgrp* ag, void* minprms, void (*egfun)(int , double* , void* , double* , double*)){
     int ndim = ag->nactives*3;
     double* xyz = (double*)_mol_malloc(ndim*sizeof(double) ); 
     double* minv = (double*)_mol_malloc(ndim*sizeof(double) );
       	ag2array(xyz, ag);
	double fmim; 
        /*my_en_grad(ndim, xyz, minprms, &fmim, NULL); 
	 printf("Start E-value=%f\n" , fmim );
	*/
	if (min_type == MOL_LBFGS) lbfgsMin(xyz,maxIt,tol,ndim,minprms, egfun,minv,&fmim);
	if (min_type == MOL_CONJUGATE_GRADIENTS) dirMin(xyz,maxIt,tol,ndim,minprms, egfun,minv,&fmim);
	if (min_type == MOL_POWELL) {
        double* directions = (double*)_mol_malloc(ndim*ndim*sizeof(double) );
        int i;
        for(i=0; i<ndim*ndim;i++) directions[i]=0;
        for(i=0; i<ndim; i++) directions[i*ndim+i]=1;
        powell(xyz,directions,maxIt,tol,ndim,minprms, egfun,minv,&fmim);
	free(directions);
	}
	/*
        my_en_grad(ndim, minv, minprms, &fmim, NULL); 
        printf("fmin=%f\n" , fmim );
        */
	array2ag(minv, ag);
	free(minv);
	free(xyz); 
}

};
