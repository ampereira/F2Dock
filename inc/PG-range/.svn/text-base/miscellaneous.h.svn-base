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


#ifndef MISC_H
#define MISC_H
#include <vector>
using namespace std;

int binSearch(int key, const vector<int>& x, int low, int high) {
  if(high < low) return -1;
  int mid;
  while(true) {
  mid = (low+high)/2;
  if(x[mid] == key) return mid;

  else if(low == high) return -1;

  else if(low == mid) low = high;

  else if(x[mid] > key) high = mid;
  else low = mid;
  }
 /* for(int i=low;i<=high;i++)
  	if(x[i]==key)
		return i;
  return -1;*/
}
#endif
