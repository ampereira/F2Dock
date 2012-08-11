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


#ifndef BINARY_SEARCH_TREE_H_
#define BINARY_SEARCH_TREE_H_
#include <vector>
#include <iostream>
#include "dsexceptions.h"

using namespace std;

// Binary node and forward declaration because g++ does
// not understand nested classes.
template <class Comparable>
class BinarySearchTree;

template <class Comparable>
class BinaryNode
{
  Comparable element;
  BinaryNode *left;
  BinaryNode *right;
  
 BinaryNode(const Comparable& theElement, BinaryNode *lt, BinaryNode *rt)
   : element(theElement), left(lt), right(rt) { }
  friend class BinarySearchTree<Comparable>;
};


// BinarySearchTree class
//
// CONSTRUCTION: with ITEM_NOT_FOUND object used to signal failed finds
//
// ******************PUBLIC OPERATIONS*********************
// void insert( x )       --> Insert x
// void remove( x )       --> Remove x
// Comparable find( x )   --> Return item that matches x
// Comparable findMin( )  --> Return smallest item
// Comparable findMax( )  --> Return largest item
// boolean isEmpty( )     --> Return true if empty; else false
// void makeEmpty( )      --> Remove all items
// void printTree( )      --> Print tree in sorted order

template <class Comparable>
class BinarySearchTree
{
  
 private:
  unsigned int num_nodes;
  BinaryNode<Comparable> *root;
  const Comparable ITEM_NOT_FOUND;
  
  /**
   * Internal method to get element field in node t.
   * Return the element field or ITEM_NOT_FOUND if t is NULL.
   */
  
  const Comparable & elementAt( BinaryNode<Comparable> *t ) const
  {
    if( t == NULL )
      return ITEM_NOT_FOUND;
    else
      return t->element;
  }
  
  /**
   * Internal method to insert into a subtree.
   * x is the item to insert.
   * t is the node that roots the tree.
   * Set the new root.
   */
  
  void insert( const Comparable & x, BinaryNode<Comparable>*& t )
  {
    if( t == NULL ) {
      t = new BinaryNode<Comparable>(x, NULL, NULL);
      num_nodes++;
    }
    else if( x < t->element )
      insert( x, t->left );
    else if( t->element < x )
      insert( x, t->right );
    return;
  }
  
  /**
   * Internal method to remove from a subtree.
   * x is the item to remove.
   * t is the node that roots the tree.
   * Set the new root.
   */
  
  
  void remove( const Comparable & x, BinaryNode<Comparable> * & t )
  {
    if( t == NULL )
      return;   // Item not found; do nothing
    if( x < t->element )
      remove( x, t->left );
    else if( t->element < x )
      remove( x, t->right );
    else if( t->left != NULL && t->right != NULL ) // Two children
      {
	t->element = findMin( t->right )->element;
	remove( t->element, t->right );
      }
    else
      {
	BinaryNode<Comparable> *oldNode = t;
	t = ( t->left != NULL ) ? t->left : t->right;
	delete oldNode;
	num_nodes--;
      }
  }
  
  
  /**
   * Internal method to find the smallest item in a subtree t.
   * Return node containing the smallest item.
   */
  
  BinaryNode<Comparable> * findMin( BinaryNode<Comparable> *t ) const
    {
      if( t == NULL )
	return NULL;
      if( t->left == NULL )
	return t;
      return findMin( t->left );
    }
  
  
  /**
   * Internal method to find the largest item in a subtree t.
   * Return node containing the largest item.
   */
  
  BinaryNode<Comparable> * findMax( BinaryNode<Comparable> *t ) const
    {
      if( t != NULL )
	while( t->right != NULL )
	  t = t->right;
      return t;
    }
  
  /**
   * Internal method to find an item in a subtree.
   * x is item to search for.
   * t is the node that roots the tree.
   * Return node containing the matched item.
   */
  
  BinaryNode<Comparable> * find( const Comparable & x, BinaryNode<Comparable> *t ) const
    {
      if( t == NULL )
	return NULL;
      else if( x < t->element )
	return find( x, t->left );
      else if( t->element < x )
	return find( x, t->right );
      else
	return t;    // Match
    }


			/**
			 * Internal method to make subtree empty.
			 */
			
			void makeEmpty( BinaryNode<Comparable> * & t )
			{
			    if( t != NULL )
			    {
				makeEmpty( t->left );
				makeEmpty( t->right );
				delete t;
			    }
			    t = NULL;
			    num_nodes = 0;
			}
	

			/**
			 * Internal method to print a subtree rooted at t in sorted order.
			 */

			void printTree( BinaryNode<Comparable> *t ) const
			{
			    if( t != NULL )
			    {
				printTree( t->left );
				cout << t->element << endl;
				printTree( t->right );
			    }
			}

			/**
			 * Internal method to clone subtree.
			 */

			BinaryNode<Comparable> * clone( BinaryNode<Comparable> * t ) const
			{
			    if( t == NULL )
				return NULL;
			    else
				return new BinaryNode<Comparable>( t->element, clone( t->left ), clone( t->right ) );
			}
			

			/*Takes a sorted list, and builds a balanced BST out of it*/

			void MakeBalanced(vector<Comparable> &sorted, int left, int right) {
			  if(left <= right) {
			    if(left == right) {
			      insert(sorted[left],root);
			      return;
			    }
			    
			    int mid = (left + right)/2;
			    insert(sorted[mid], root);
			    MakeBalanced(sorted, left, mid-1);
			    MakeBalanced(sorted, mid+1, right);
			  }	
			}
			   

			/*merges 2 sorted lists and puts result in 3rd list*/
			
			void mergeSortedLists(vector<Comparable>& t1, vector<Comparable>& t2, vector <Comparable>& merged)
			{
			  int len = (int)t1.size() + (int)t2.size();
			  merged.resize(len);
			  
			  int left_index = 0;
			  int right_index = 0;
			  for(int i = 0; i < len; i++) {
			    if(left_index < t1.size() && right_index < t2.size()) {
			      if(t1[left_index] < t2[right_index] ) {
				merged[i] = t1[left_index++];
			      }
			      else {
				merged[i] = t2[right_index++];
			      }
			    }
			    else {
			      if(left_index < t1.size()) {
				merged[i] = t1[left_index++];
			      }
			      else {
				merged[i] = t2[right_index++];
			      }
			    }
			  }
			}

			/*traverses the tree rooted at 't' in inorder and stores the sorted list in the given vector*/
			void Inorder(vector <Comparable> &inorder_vector, const BinaryNode<Comparable> *t)
			{
				if(t== NULL) return;
				Inorder(inorder_vector, t->left);
				inorder_vector.push_back(t->element);
				Inorder(inorder_vector, t->right);
			}
			

        	public:
            
			explicit BinarySearchTree( const Comparable & notFound):num_nodes(0), root(NULL), ITEM_NOT_FOUND( notFound )
        		{
        		}


			/**
			 * Copy constructor.
			 */
		
			BinarySearchTree( const BinarySearchTree<Comparable> & rhs ) :  num_nodes(rhs.num_nodes), root( NULL ), ITEM_NOT_FOUND( rhs.ITEM_NOT_FOUND )
			{ 
			    *this = rhs;
			    //num_nodes = rhs.size();
			}

			/**
			 * Destructor for the tree.
			 */
			~BinarySearchTree( )
			{
			    makeEmpty();
			}

			/**
			 * Find the smallest item in the tree.
			 * Return smallest item or ITEM_NOT_FOUND if empty.
			 */

			const Comparable & findMin( ) const
			{
			    return elementAt( findMin( root ) );
			}

			/**
			 * Find the largest item in the tree.
			 * Return the largest item of ITEM_NOT_FOUND if empty.
			 */
			
			const Comparable & findMax( ) const
			{
			    return elementAt( findMax( root ) );
			}

			/**
			 * Find item x in the tree.
			 * Return the matching item or ITEM_NOT_FOUND if not found.
			 */
			
			const Comparable & find( const Comparable & x ) const
			{
			    return elementAt( find( x, root ) );
			}

			/**
			 * Make the tree logically empty.
			 */

			void makeEmpty( )
			{
			    makeEmpty( root );
			}

			/**
			 * Test if the tree is logically empty.
			 * Return true if empty, false otherwise.
			 */

			bool isEmpty( ) const
			{
			    return root == NULL;
			}

			/**
			 * Print the tree contents in sorted order.
			 */

			void printTree( ) const
			{
			  if( isEmpty( ) ) 
			      cout << "Empty tree" << endl;
			    else
				printTree( root );
			}


		        /**
			 * Insert x into the tree; duplicates are ignored.         */
			
			void insert( const Comparable & x )
			{
			    insert( x, root );
			}

			/**
			 * Remove x from the tree. Nothing is done if x is not found.
			 */

			void remove( const Comparable & x )
			{
			    remove( x, root );
			}


			/**
			 * Deep copy.
			 */

			const BinarySearchTree<Comparable> & operator=( const BinarySearchTree<Comparable> & rhs )
			{
			    if( this != &rhs )
			    {
				makeEmpty( );
				root = clone( rhs.root );
			    }
			    return *this;
			}

		    	/* Find the successor*/
			Comparable successor(const Comparable& x) {

			  BinaryNode<Comparable> *t = find(x, root);
			  if(t -> right != NULL)
			    return findMin(t -> right)->element;

			  if(x == findMax(root)->element)
			    return ITEM_NOT_FOUND;

			  BinaryNode<Comparable> *s = NULL;
			  t = root;
			  while(t != NULL) {
			    if(x < t->element) {
			      s = t;
			      t = t->left;
			    }
			    else if(x >  t->element)
			      t = t->right;
			    else
			      break;
			  }
			  return s->element;
			}

			/*Find the predecessor*/
			Comparable predecessor(const Comparable& x) {

			  BinaryNode<Comparable> *t = find(x, root);
			  /*
			  if(t == NULL) {}

			  else if(t -> left != NULL)
			    return findMax(t -> left)->element;

			  if(x == findMin(root)->element)
			    return -1;
			  */
			  if(t) return t->element;
			  BinaryNode<Comparable> *s = NULL;
			  t = root;
			  while(t != NULL) {
			    if(x < t->element) {
			      t = t->left;
			    }
			    else if(x > t->element) {
			      s = t;
			      t = t->right;
			    }
			    else
			      break;
			  }
			  if (s) return s->element;
			  else return ITEM_NOT_FOUND;
			}


			/*Return the size of BST*/
			unsigned int size() {
			  return num_nodes;
			  //return findNumNodes(root);
			}

		    	/*Split a BST*/
			BinarySearchTree<Comparable>* split(BinarySearchTree<Comparable>* that) {
			  vector<Comparable> inorder_bst;
			  Inorder(inorder_bst, root);
			  
			  int low = 0;
			  int high = (int) inorder_bst.size() - 1;
			  int mid = (low+high)/2;

			  makeEmpty();
			  MakeBalanced(inorder_bst, low, mid);
			  that->MakeBalanced(inorder_bst, mid+1, high);

			  return that;
			}

			/*Merge two binary search tree and merged tree is stored in this; the other tree remains unaltered*/
			void merge(BinarySearchTree<Comparable>* T1) {
			  
			  vector <Comparable> inorder_bst1;
			  vector <Comparable> inorder_bst2;
			  
			  Inorder(inorder_bst1, root);
			  Inorder(inorder_bst2, T1->root);

			  vector <Comparable> inorder_merged;

			  mergeSortedLists(inorder_bst1,inorder_bst2,inorder_merged);
			  makeEmpty();
			  MakeBalanced(inorder_merged, 0, (int)inorder_merged.size()-1);
			}
        };
        

#endif
