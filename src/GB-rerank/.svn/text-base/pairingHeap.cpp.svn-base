
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

#include "GB-rerank/pairingHeap.h"

int PairingHeap::new_node( int x, VAL_TYPE k )
  {
   int fp;

   if ( freeptr != NIL )
      {
       fp = freeptr;
       freeptr = ( *heap )[ freeptr + CHILD ];

       ( *heap )[ fp + ID ] = x;
       ( *heap )[ fp + KEY ] = k;
       ( *heap )[ fp + CHILD ] = NIL;
       ( *heap )[ fp + LEFT ] = NIL;
       ( *heap )[ fp + RIGHT ] = NIL;
      }
   else
      {
       fp = csize;

       heap->push_back( x );        // BACK
       heap->push_back( k );        // KEY
       heap->push_back( NIL );      // CHILD
       heap->push_back( NIL );      // LEFT
       heap->push_back( NIL );      // RIGHT

       csize += REC_SIZE;
      }

   return fp;
  }


void PairingHeap::free_node( int xp )
  {
   ( *heap )[ xp + CHILD ] = freeptr;
   freeptr = xp;
  }


void PairingHeap::add_to_aux_area( int xp )
  {
   if ( auxptr == NIL )
      min_auxptr = auxptr = ( *heap )[ xp + LEFT ] = ( *heap )[ xp + RIGHT ] = xp;
   else
     {
      int yp = ( *heap )[ auxptr + LEFT ];

      ( *heap )[ xp + LEFT ] = yp;
      ( *heap )[ xp + RIGHT ] = auxptr;

      ( *heap )[ yp + RIGHT ] = xp;
      ( *heap )[ auxptr + LEFT ] = xp;

      if ( ( *heap )[ min_auxptr + KEY ] > ( *heap )[ xp + KEY ] ) min_auxptr = xp;
     }
  }


inline int PairingHeap::pairing_heap_link( int xp, int yp )
  {
   if ( xp == NIL ) return yp;
   else if ( yp == NIL ) return xp;

   if ( ( *heap )[ xp + KEY ] > ( *heap )[ yp + KEY ] )
     {
      int t = xp;

      xp = yp; yp = t;
     }

   ( *heap )[ yp + LEFT ] = xp;
   ( *heap )[ yp + RIGHT ] = ( *heap )[ xp + CHILD ];
   if ( ( *heap )[ yp + RIGHT ] != NIL )
      ( *heap )[ ( *heap )[ yp + RIGHT ] + LEFT ] = yp;
   ( *heap )[ xp + CHILD ] = yp;

   ( *heap )[ xp + LEFT ] = ( *heap )[ xp + RIGHT ] = NIL;

   return xp;
  }


int PairingHeap::two_pass_merge( int xp )
  {
   if ( xp == NIL ) return NIL;

   int xp1, xp2, yp, zp = NIL;

   while ( xp != NIL )
     {
      xp1 = xp;

      if ( ( *heap )[ xp1 + RIGHT ] != NIL ) xp2 = ( *heap )[ xp1 + RIGHT ];
      else xp2 = NIL;

      if ( xp2 == NIL ) xp = NIL;
      else 
         {
          if ( ( *heap )[ xp2 + RIGHT ] != NIL ) xp = ( *heap )[ xp2 + RIGHT ];
          else xp = NIL;
         }

      if ( zp == NIL ) 
         {
          zp = pairing_heap_link( xp1, xp2 );
          ( *heap )[ zp + RIGHT ] = NIL;
         }
      else
         {
          yp = pairing_heap_link( xp1, xp2 );
          ( *heap )[ yp + RIGHT ] = zp;
          zp = yp;
         }
     }

   xp = NIL;
   while ( zp != NIL )
     {
      yp = zp;
      zp = ( *heap )[ zp + RIGHT ];

      xp = pairing_heap_link( xp, yp );
     }

   return xp;
  }


int PairingHeap::multi_pass_merge( int xp )
  {
   if ( xp == NIL ) return NIL;

   int xp1, xp2, yp, zp;

   while ( ( *heap )[ xp + RIGHT ] != NIL )
     {
      zp = NIL;

      while ( xp != NIL )
        {
         xp1 = xp;

         if ( ( *heap )[ xp1 + RIGHT ] != NIL ) xp2 = ( *heap )[ xp1 + RIGHT ];
         else xp2 = NIL;

         if ( xp2 == NIL ) xp = NIL;
         else 
            {
             if ( ( *heap )[ xp2 + RIGHT ] != NIL ) xp = ( *heap )[ xp2 + RIGHT ];
             else xp = NIL;
            }

         if ( zp == NIL ) 
            {
             zp = pairing_heap_link( xp1, xp2 );
             ( *heap )[ zp + RIGHT ] = NIL;
            }
         else
            {
             yp = pairing_heap_link( xp1, xp2 );
             ( *heap )[ yp + RIGHT ] = zp;
             zp = yp;
            }
        }

      xp = zp;
     }

   return xp;
  }


int PairingHeap::Insert( int x, VAL_TYPE k )
  {
   int xp = new_node( x, k );

   if ( root == NIL ) root = xp;
   else 
      {
       if ( use_aux_trees ) add_to_aux_area( xp );
       else root = pairing_heap_link( xp, root );
      }

   n++;

   return xp;
  }


void PairingHeap::Find_Min( int &mx, VAL_TYPE &mk )
  {
   if ( root == NIL ) 
       {
        mx = NIL;
        mk = INF;
       }
   else
       {
        if ( ( use_aux_trees == false ) || ( min_auxptr == NIL ) || ( ( *heap )[ min_auxptr + KEY ] > ( *heap )[ root + KEY ] ) )
          {
           mx = ( *heap )[ root + ID ];
           mk = ( *heap )[ root + KEY ];
          }
        else
          {
           mx = ( *heap )[ min_auxptr + ID ];
           mk = ( *heap )[ min_auxptr + KEY ];
          }
       }
  }


void PairingHeap::Delete_Min( int &mx, VAL_TYPE &mk )
  {
   if ( root == NIL )
       {
        mx = NIL;
        mk = INF;
       }
   else
       {
        if ( use_aux_trees && ( auxptr != NIL ) )
          {
           ( *heap )[ ( *heap )[ auxptr + LEFT ] + RIGHT ] = NIL;
           ( *heap )[ auxptr + LEFT ] = NIL;

           auxptr = multi_pass_merge( auxptr );
           root = pairing_heap_link( auxptr, root );

           auxptr = min_auxptr = NIL;
          }

        mx = ( *heap )[ root + ID ];
        mk = ( *heap )[ root + KEY ];

        int xp = ( *heap )[ root + CHILD ];

        free_node( root );

        if ( use_multi_pass ) root = multi_pass_merge( xp );
        else root = two_pass_merge( xp );

        n--;
       }
  }


void PairingHeap::Decrease_Key( int xp, VAL_TYPE k )
  {
   if ( k < ( *heap )[ xp + KEY ] )
     {
      ( *heap )[ xp + KEY ] = k;

      if ( xp != root )
        {
         if ( ( *heap )[ xp + RIGHT ] != NIL )
            ( *heap )[ ( *heap )[ xp + RIGHT ] + LEFT ] = ( *heap )[ xp + LEFT ];

         if ( ( *heap )[ ( *heap )[ xp + LEFT ] + CHILD ] == xp )
            ( *heap )[ ( *heap )[ xp + LEFT ] + CHILD ] = ( *heap )[ xp + RIGHT ];
         else
             ( *heap )[ ( *heap )[ xp + LEFT ] + RIGHT ] = ( *heap )[ xp + RIGHT ];

         ( *heap )[ xp + LEFT ] = ( *heap )[ xp + RIGHT ] = NIL;

         if ( use_aux_trees ) add_to_aux_area( xp );
         else root = pairing_heap_link( xp, root );
        }
     }
  }
