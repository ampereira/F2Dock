/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@cs.utexas.edu>
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


#ifndef _PG_H_
#define _PG_H_

#include <vector>
#include "1D-IntegerRange.h"
#include <cmath>

//#include "grid.h"
//#include "Point.h"
//#include <arrangement_2d.h>
using namespace std;

#define INITHASHSIZE 51

//struct grid;
//template <class T> class IntegerRange;

struct Point {
  float x;
  float y;
  float z;
  Point() {}
  Point(float a, float b, float c) {
    x = a;
    y = b;
    z = c;
  }
  bool operator==(const Point& p) {
    return (x==p.x && y==p.y && z==p.z);
  }
 
  float distsq(Point a, Point b) {
  float res;
  double dx, dy, dz;
  dx = a.x - b.x;
  dy = a.y - b.y;
  dz = a.z - b.z;  
  res = dx*dx + dy*dy + dz*dz;
  return res;
}

};

struct cellID {
  int x;
  int y;
  int z;
};

struct lineID {
  int y;
  int z;
};

struct planeID {
  int z;
};

struct gridcell {
  cellID ID;
  vector<Point*> balls;
  gridcell(int a, int b, int c) {
    ID.x = a;
    ID.y = b;
    ID.z = c;
  }
};

struct line {
  lineID ID;
  IntegerRange<gridcell*> RR;
  line(int b, int c) {
    ID.y = b;
    ID.z = c;
	RR = IntegerRange<gridcell*>();
  }
};

struct plane {
  planeID ID;
  IntegerRange<line*> RR;
  plane(int c) {
    ID.z = c;
    RR = IntegerRange<line*>();
  }
};

struct grid {
  IntegerRange<plane*> RR;
  grid()
  {
	RR = IntegerRange<plane*>();
  }
//  ball* surface_root_ptr;
};

class PG {
  grid g;
  int cells;
  int lines;
  int planes;
  double rmax;
  double rangeTime, compTime, initTime;
  double minx, miny, minz, maxx, maxy, maxz;
  double nx, ny, nz;
  int rangeCount;
  double DIM;
//  dynamic_graph* dg;
 public:
  int s1size;
  int s2size;
  int checked;
  double TRANSLATE;

  PG(double D, double xlate, double r)
  {
    DIM = D;
    rmax = r;
    TRANSLATE = xlate;
    cells = 0;
    lines = 0;
    planes = 0;
    rangeTime = 0.0;
    compTime = 0.0;
    rangeCount = 0;
    initTime = 0;
}

  PG(double mnx, double mny, double mnz, double mxx, double mxy, double mxz, double divisionsize, int inithashsize = INITHASHSIZE) 
  {
    maxx = mxx;
    maxy = mxy;
    maxz = mxz;
    minx = mnx;
    miny = mny;
    minz = mnz;
    if(divisionsize!=0)
	DIM = divisionsize;
    else
	DIM = 3.5;
    rmax = 1.5;
    cells = 0;
    lines = 0;
    planes = 0;
    rangeTime = 0.0;
    compTime = 0.0;
    rangeCount = 0;
    initTime = 0;
  }

  PG(vector<Point *> *alist, double divisionsize, int inithashsize = INITHASHSIZE);

  void boundingbox(double &mnx, double &mny, double &mnz, double &mxx, double &mxy, double &mxz)
  {
	mnx = this->minx;
	mny = this->miny;
	mnz = this->minz;
	mxx = this->maxx;
	mxy = this->maxy;
	mxz = this->maxz;
  }

  /* Supported Queries*/
  vector<Point*> range(Point *, double);
  bool pointsWithinRange(Point *q, double delta); 
  int countPointsWithinRange(Point *q, double delta);
//  vector<Point*> *findintersections(Point *, double);
// vector<Point*> *findintersection(Point * a, double d) {return findintersections(a,d);}
//  vector<Point*> *findintersection(Point * a) {return findintersections(a,0.0);}
//  vector<Point*> intersect(Point, float);
//  bool exposed(Point*);
//  void surface(void);
//  void checkInsertion(FILE *fp);

  double getRangeTime()
  {
  	return rangeTime;
  }
  int getRangeCount()
  {
  	return rangeCount;
  }
  double getCompTime()
  {
  	return compTime;
  }
  double getInitTime()
  {
  	return initTime;
  }
  void resetTimes()
  {
  	rangeTime = 0;
	compTime = 0;
	rangeCount = 0;
	initTime =0;
  }


 double getdivsize() { return DIM; }
  
 int cellsstored(void) { return cells; }

 int tablesize(void) { return cells; }

 int gridsize(void) { return cells; }

  
  /* Supported updates*/
  void addPoint(Point *a);
  void addPoints(vector<Point *> *alist)
  {
	int i, size;
	size = alist->size();
	for(i=0;i<size;i++)
		addPoint(alist->at(i));
  }
  void removePoint(Point *a);
//  void removePoints(vector<Point *> *alist);
  //void move(Point, Point, float);
};
#endif
