/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
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


#ifndef _PAIR_H_

#define _PAIR_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <limits.h>

#ifndef VAL_TYPE
   #define VAL_TYPE double
#endif

#ifndef INF
  #define INF ( ( VAL_TYPE ) INT_MAX )
#endif

typedef std::vector< VAL_TYPE > data_type;

class PairingHeap
{  
 private:

  enum misc{ NIL = -1, REC_SIZE = 5 };
  enum rec_ids{ ID, KEY, CHILD, LEFT, RIGHT };

  data_type *heap;

  int n, root, auxptr, min_auxptr, freeptr, csize;
  bool use_multi_pass, use_aux_trees;
  
  void basic_init( bool use_mpass, bool use_aux )
    {
     n = 0;
     root = NIL;
     freeptr = NIL;
     auxptr = NIL;
     min_auxptr = NIL;

     heap = new data_type;
     csize = 0;

     use_multi_pass = use_mpass;
     use_aux_trees = use_aux;
    }
  
  int new_node( int x, VAL_TYPE k );
  void free_node( int xp );
  void add_to_aux_area( int xp );
  int pairing_heap_link( int xp, int yp );
  int two_pass_merge( int xp );
  int multi_pass_merge( int xp );

 public:

  PairingHeap( void )
    {
      basic_init( false, false );
    }

  PairingHeap( bool use_mpass )
    {
      basic_init( use_mpass, false );
    }

  PairingHeap( bool use_mpass, bool use_aux )
    {
      basic_init( use_mpass, use_aux );
    }

  ~PairingHeap( )
    {
      delete heap; 
    }


  int isEmpty( void )
    {
      return !n;
    }

    
  int Insert( int x, VAL_TYPE k );
  void Find_Min( int &mx, VAL_TYPE &mk );
  void Delete_Min( int &mx, VAL_TYPE &mk );
  void Decrease_Key( int xp, VAL_TYPE k );
};

#endif
