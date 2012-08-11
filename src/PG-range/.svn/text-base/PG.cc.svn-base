/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/


#include "PG.h"
#include <cmath>
#include <time.h>
#include <sys/time.h>

using namespace std;

//#define TRANSLATE 2000.0

double gtod_sec = 0.0E0;

double gtod_timer()
{
   struct timeval tv;
   struct timezone Tzp;
   double sec;

   gettimeofday(&tv, &Tzp);

   if(gtod_sec == 0.0E0) 
      gtod_sec = (double)tv.tv_sec;
   sec = (double)tv.tv_sec - gtod_sec;

   return sec + 1.0E-06*(double)tv.tv_usec;
}




vector<Point*> PG::range(Point *q, double delta) 
{
  Point p = *q;

  cout<<"Point = "<<q->x<<" "<<q->y<<" "<<q->z<<" "<<endl;

  p.x += TRANSLATE;
  p.y += TRANSLATE;
  p.z += TRANSLATE;

  int l = (int)(( p.z - delta) / DIM);
  int h = (int)((p.z + delta) / DIM);

  int i,j,k,size,size1,size2,m;
  
  vector<Point*> result;
  
  vector<tuple<line*> >  temp1;
  vector<tuple<gridcell*> >  temp2, S0, S1;
  vector<tuple<gridcell*> >::iterator start2, end2; 
  
  

  int tts;


  vector<tuple<plane*> > S2 = g.RR.report(l,h);
  rangeCount++;

  if(S2.empty())
  {
    return result;
  }
  size2 = (int)S2.size();


  for(i = 0; i < size2; i++) 
  {
    rangeCount++;
   
    l = (int)((p.y - delta)/DIM);
    h = (int)((p.y + delta)/DIM);
    temp1 = (S2[i].ptr)->RR.report(l,h);
    
    size1 = (int) temp1.size();
    
    for(j = 0; j < size1; j++) 
    {
      rangeCount++;

      l = (int)((p.x - delta)/DIM);
      h = (int)((p.x + delta)/DIM);


      if((temp1[j].ptr)->RR.getn())
      {
	temp2 = (temp1[j].ptr)->RR.report(l,h);
        
	size = (int)temp2.size();
	int pbcount;

	Point *oa;

	double delsq = delta*delta;
	for(int mn=0;mn<size;mn++)
	{
		int atomsInCell = temp2[mn].ptr->balls.size();
		for(int pq=0;pq<atomsInCell;pq++)
		{
			oa = temp2[mn].ptr->balls[pq];
			if(oa->distsq(*oa, *q) <= delsq)
				result.push_back(oa);
		}
	}
	
      	temp2.clear();
	
      }
    }
    temp1.clear();
  }

  return result;
}



bool PG::pointsWithinRange(Point *q, double delta) 
{
  Point p = *q;

  p.x += TRANSLATE;
  p.y += TRANSLATE;
  p.z += TRANSLATE;

  int l = (int)(( p.z - delta) / DIM);
  int h = (int)((p.z + delta) / DIM);

  int i,j,k,size,size1,size2,m;
  
  vector<tuple<line*> >  temp1;
  vector<tuple<gridcell*> >  temp2, S0, S1;
  vector<tuple<gridcell*> >::iterator start2, end2; 
  
  int tts;


  vector<tuple<plane*> > S2 = g.RR.report(l,h);
  rangeCount++;

  if(S2.empty())
  {
    return false;
  }
  size2 = (int)S2.size();


  for(i = 0; i < size2; i++) 
  {
    rangeCount++;
   
    l = (int)((p.y - delta)/DIM);
    h = (int)((p.y + delta)/DIM);
    temp1 = (S2[i].ptr)->RR.report(l,h);
    
    size1 = (int) temp1.size();
    
    for(j = 0; j < size1; j++) 
    {
      rangeCount++;

      l = (int)((p.x - delta)/DIM);
      h = (int)((p.x + delta)/DIM);


      if((temp1[j].ptr)->RR.getn())
      {
	temp2 = (temp1[j].ptr)->RR.report(l,h);
        
	size = (int)temp2.size();
	int pbcount;

	Point *oa;

	double delsq = delta*delta;
	for(int mn=0;mn<size;mn++)
	{
		int atomsInCell = temp2[mn].ptr->balls.size();
		for(int pq=0;pq<atomsInCell;pq++)
		{
			oa = temp2[mn].ptr->balls[pq];
			if(oa->distsq(*oa, *q) <= delsq) return true;
		}
	}
	
      	temp2.clear();
	
      }
    }
    temp1.clear();
  }

  return false;
}



int PG::countPointsWithinRange(Point *q, double delta) 
{
  Point p = *q;

  p.x += TRANSLATE;
  p.y += TRANSLATE;
  p.z += TRANSLATE;

  int l = (int)(( p.z - delta) / DIM);
  int h = (int)((p.z + delta) / DIM);

  int i,j,k,size,size1,size2,m;
  
  vector<tuple<line*> >  temp1;
  vector<tuple<gridcell*> >  temp2, S0, S1;
  vector<tuple<gridcell*> >::iterator start2, end2; 
  
  int tts;

  vector<tuple<plane*> > S2 = g.RR.report(l,h);
  rangeCount++;

  if(S2.empty())
  {
    return 0;
  }
  size2 = (int)S2.size();

  int count = 0;

  for(i = 0; i < size2; i++) 
  {
    rangeCount++;
   
    l = (int)((p.y - delta)/DIM);
    h = (int)((p.y + delta)/DIM);
    temp1 = (S2[i].ptr)->RR.report(l,h);
    
    size1 = (int) temp1.size();
    
    for(j = 0; j < size1; j++) 
    {
      rangeCount++;

      l = (int)((p.x - delta)/DIM);
      h = (int)((p.x + delta)/DIM);


      if((temp1[j].ptr)->RR.getn())
      {
	temp2 = (temp1[j].ptr)->RR.report(l,h);
        
	size = (int)temp2.size();
	int pbcount;

	Point *oa;

	double delsq = delta*delta;
	for(int mn=0;mn<size;mn++)
	{
		int atomsInCell = temp2[mn].ptr->balls.size();
		for(int pq=0;pq<atomsInCell;pq++)
		{
			oa = temp2[mn].ptr->balls[pq];
			if(oa->distsq(*oa, *q) <= delsq) count++;
		}
	}
	
      	temp2.clear();
	
      }
    }
    temp1.clear();
  }

  return count;
}


/*adds a ball to the collection of balls*/
void PG::addPoint(Point *a) 
{
//  cout<<"INSERTING AN ATOM "<<endl;

  Point center = *a;

  center.x += TRANSLATE;
  center.y += TRANSLATE;
  center.z += TRANSLATE;
//  double radius = a->getr();

  int cz = (int)center.z/DIM;
  int cy = (int)center.y/DIM;
  int cx = (int)center.x/DIM;
  tuple<plane*> temp;
  plane* p;
  line* l;
  gridcell* c;

//  if(rmax<radius)
//	rmax = radius;

  vector<tuple<plane*> > P = g.RR.report(cz,cz);
  
  if(P.empty()) 
  {
    planes++;
    p = new plane(cz);
    temp = tuple<plane*>(cz, p);
    P.push_back(temp);
    g.RR.insert(cz, p);
  }
  
  vector<tuple<line*> > L = P[0].ptr->RR.report(cy,cy);
  
  if(L.empty()) 
  {
    lines++;
    l = new line(cy, cz); 
    tuple<line*> temp1(cy, l);
    L.push_back(temp1);
    P[0].ptr->RR.insert(cy, l); //same here
  }
  
  vector<tuple<gridcell*> > C = L[0].ptr->RR.report(cx,cx); //same here
  
  if(C.empty()) 
  {
    cells++;
    c = new gridcell(cx, cy, cz); //and here
    tuple<gridcell*> temp2(cx, c);
    C.push_back(temp2);
    L[0].ptr->RR.insert(cx, c); //and here
  }

  (C[0].ptr->balls).push_back(a); //and here

//cout<<"Inserted"<<endl;

return;
}

/*removes the given ball from the collection of balls*/
void PG::removePoint(Point *a) 
{
cout<<"Inside removeatom"<<endl;

  Point center = *a;

  center.x += TRANSLATE;
  center.y += TRANSLATE;
  center.z += TRANSLATE;
//  double radius = a->getr();

  cout<<"REMOVING AN ATOM "<<endl;
  int cz = (int)center.z/DIM;
  int cy = (int)center.y/DIM;
  int cx = (int)center.x/DIM;
  
  vector<tuple<plane*> > P = g.RR.report(cz,cz);
  if(P.empty()) {
    cout<<"Atom does not exist"<<endl;
    return;
  }

  vector<tuple<line*> > L = P[0].ptr->RR.report(cy,cy);
  if(L.empty()) {
    cout<<"Atom does not exist"<<endl;
    return;
  }

  vector<tuple<gridcell*> > C = L[0].ptr->RR.report(cx,cx); //same here
  if(C.empty()) {
    cout<<"Atom does not exist"<<endl;
    return;
  }

  int index;
  for(int i=0; i<(int)(C[0].ptr->balls).size(); i++) 
  {
    if(center.x == C[0].ptr->balls[i]->x && center.y == C[0].ptr->balls[i]->y && center.z == C[0].ptr->balls[i]->z) 
    {
      index = i;
      break;
    }
  }
  C[0].ptr->balls[index] = C[0].ptr->balls[C[0].ptr->balls.size()-1];
  C[0].ptr->balls.pop_back();
  
cout<<"return from removeatom"<<endl;

  return;
}


