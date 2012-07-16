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


#ifndef D_RANGE
#define D_RANGE


#include "cuckoo.h"
#include "BinarySearchTree.h"
//#include "hash3d.h"
//#include "miscellaneous.h"
#include <string>
#include <iostream>
#include <cmath>
using namespace std;

//#define W 4
//#define U 16
#define SPLIT_THRESH 2*n/W
#define MERGE_THRESH n/(2*W)


template<typename T>
class tuple {
 public:
  int id;
  T ptr;
 
  tuple(int a, T b) {
    id = a;
    ptr = b;
  }
  tuple() {
    id = -1;
    ptr = NULL;
  }

  tuple(const tuple<T>& that) {
      id = that.id;
      ptr = that.ptr;

  }

  bool operator==(const tuple<T>& that) const{
    return (id == that.id);
  }
  bool operator<(const tuple<T>& that) const{
    return id < that.id;
  }
  bool operator>(const tuple<T>& that) const{
    return id > that.id;
  }
  bool operator<=(const tuple<T>& that) const{
    return (operator<(that) || operator==(that));
  }
  bool operator>=(const tuple<T>& that) const{
    return (operator>(that) || operator==(that));
  }
  bool operator!=(const tuple<T>& that) const{
    return !(operator==(that));
  }

  /*
  ostream& operator<<(ostream& ost) {
    string end = "";
    ost << this->id;
    return ost << end;
  }
  */

};

template<typename T>
ostream& operator<<(ostream& ost, const tuple<T>& p) {
  string end = "";
  ost<<p.id<<" "<<p.ptr;
  return ost << end;
}

template <class T>
class IntegerRange 
{
	int W;
	int U;
	int n;
	vector<BinarySearchTree<tuple<T> >*> BST;
	dict_ptr D;
	vector<int> repMap; // maps cluster id with representative
	//tuple<T> notfound(-1, NULL);

	int get_pre_rep(int key) 
	{
		int m = U+key;
  
		/*if n is in hash, then key is returned*/
		if(lookup(D, m)) return key;

		int c = m;
		m = m/2;
		
		while(m) 
		{
	    		if(lookup(D, m) && lookup(D, 2*m) && c == 2*m+1)
	      			break;
	    		c = m;
			m = m/2;
		}

		if(!m)
	    		return -1;

		// cout<<"******* before "<<m<<endl;
		int level = (int) log2((double)m);
		int i = level+2;
		c=m;
		m = 2*m;
		/*
		for(; i<=W; i++) 
		{
			if(lookup(D, 2*m+1)) 	
				m = 2*m+1;
			else 
				m = 2*m;
		}*/

		while(1)
		{
			if(c==m) 
				break;

			c=m;
		
			if(lookup(D, 2*m+1)) 
				m = 2*m+1;
			else if(lookup(D, 2*m)) 
				m = 2*m;
		}

		// cout<<"******* after "<<m<<endl;
		return m-U;
	}



	/*Prints the Binary search trees*/

	void printBSTs() 
	{
		int num_clusters = (int)BST.size();
		for(int i = 0; i < num_clusters; i++) 
		{
			cout<<"Printing tree # "<<i+1<<endl;
			BST[i]->printTree(); /*cout for tuple*/
		}
		return;
	}


	/*Prints the representative element of each cluster*/

	void printReps() 
	{
		int num_clusters = (int)repMap.size();
	
		for(int i = 0; i < num_clusters; i++) 
			cout<<"Representative for tree # "<<i+1<<" is "<<repMap[i]<<endl;

		return;
	}


	/*Insert x and all its binary prefixes into the hash*/

	void insertInHash(int x) 
	{
		int k = U + x;
		while(k) 
		{
    			insertD(D, k);
    			k = k/2;
  		}
		return;
	}


	/*Remove x and all its prefixes from the hash*/

	void removeFromHash(int x) 
	{
		delete_key(D, x);
		int level = (int) log2((double)x);
		int c = x;

		for(int i = level-1; i>=0; i--) 
		{
			x = c/2;
			if(lookup(D, 2*x) || lookup(D, 2*x+1)) 
				break;
			delete_key(D, x);
			c = x;
		}
		return;
	}


	/*Given the index of a BST, 
	  splits the tree into two, if applicable
	*/

	void splitBST(int index) 
	{
		int s = (int)BST[index]->size();
  
		if(s > SPLIT_THRESH) 
		{
			tuple<T> notfound(-1, NULL);
			BinarySearchTree<tuple<T> >* t = new BinarySearchTree<tuple<T> >(notfound);
			BST[index]->split(t);
	
			//   cout<<"In split"<<endl;
			//   printReps();
			removeFromHash(repMap[index]);
			repMap[index] = (BST[index]->findMax()).id;
			vector<int>::iterator p = repMap.begin();
			p = p + (index + 1);
			repMap.insert(p, (t->findMax()).id);

			typename vector<BinarySearchTree<tuple<T> >*>::iterator q = BST.begin();
			q = q + index + 1;
			BST.insert(q, t);

			insertInHash(repMap[index]);
			insertInHash(repMap[index+1]);
			//    printReps();
		}
	
		return;
	}



	int binSearch(int key, const vector<int>& x, int low, int high) 
	{
		if(high < low) return -1;
		int mid;
		
		while(true) 
		{
			mid = (low+high)/2;
			if(x[mid] == key) 
				return mid;

			else if(low == high) 	
				return -1;

			else if(low == mid) 
				low = high;

			else if(x[mid] > key) 
				high = mid;

			else low = mid;
		}
	
		/* for(int i=low;i<=high;i++)
		  	if(x[i]==key)
				return i;
		  return -1;*/
	}



	/*Given the index of a BST,
	  merges the input tree with an adjacent tree, if applicable, and
	  returns the index of the merged tree.(-1 if no merge took place)
	*/

	int mergeBST(int index) 
	{
		int s = (int)BST[index]->size();
	
		if(s <= MERGE_THRESH) 
		{
			int merge_with = index + 1;
			if(index == ((int)BST.size()-1))
				merge_with = index-1;

			removeFromHash(repMap[merge_with]);

			BST[index]->merge(BST[merge_with]);
			removeFromHash(repMap[index]);
			repMap[index] = (BST[index]->findMax()).id;
			insertInHash(repMap[index]);

			typename vector<BinarySearchTree<tuple<T> >*>::iterator p = BST.begin();
			p = p + merge_with;
			BST.erase(p);

			vector<int>::iterator q = repMap.begin();
			q = q + merge_with;
			repMap.erase(q);

			if(merge_with < index)
				index = merge_with;
			return index;
		}
	
		else 
			return -1;
	}


	/*Given the predecessor representative of x and x,
	  returns the cluster id of the relevant cluster
	*/

	int get_cluster(int pre_rep, int x) 
	{
		if(pre_rep < 0)
			return 0;

		if(pre_rep == x)
			return binSearch(pre_rep, repMap, 0, (int)repMap.size()-1);

		return binSearch(pre_rep, repMap, 0, (int)repMap.size()-1) + 1;
	}


public:
	
	/*Constructors*/

	IntegerRange()
	{
		n = 0;
		W = 4;
		U = 1024;
		D = construct_dict(1024);
	}
  
	int getn()
	{
		return n;
	}
 

	/*Insert into S*/
	void insert(int x, T y)
	{
		//cout<<"Inserting "<<"x = "<<x<<endl;
		//printReps();
		int pre_rep = get_pre_rep(x);
		int next_clus = get_cluster(pre_rep, x);
		//cout<<" pre = "<<pre_rep<<" next = "<<next_clus<<endl;
 
		/*
		If new cluster needs to be created, then:
		1. insert representative element of new cluster into hash
		2. record the same in repMap
		3. construct a new BST for the new cluster
		*/

		int s = (int)BST.size();
		if(next_clus == s) 
		{
			//cout<<"Creating new BST"<<endl;
    
			insertInHash(x);
			repMap.push_back(x);
			tuple<T> notfound(-1, NULL);
			BST.push_back(new BinarySearchTree<tuple<T> >(notfound));
    
			//printReps();
  		}
  
		/*
		Insert x into appropriate BST, and 
		increment n
		*/

		tuple<T> item(x,y);
  
		BST[next_clus]->insert(item);

		//BST[next_clus]->printTree();
		
		n++;
		//cout<<"n = "<<n<<endl;
  
		if(BST[next_clus]->size()*2/W < 2 )
		{
			return;
		} 
  
  		splitBST(next_clus);
		//  cout<<"Splitted"<<endl;
		//  printReps();
		return;
  	}


	/*Delete from S*/

	void remove(int k)
	{
		int pre_rep = get_pre_rep(k);
		int next_clus = get_cluster(pre_rep, k);
 
		tuple<T> item(k,NULL);
		BST[next_clus]->remove(item);
		n--;
  
		/*
		If removing k from BST[next_clus] makes it empty, do the following and return:
		1. Delete BST[next_clus]
		2. delete corresponding entry from repMap
		3. remove k from hash
		*/
 
		if(BST[next_clus]->isEmpty()) 
		{
			delete BST[next_clus];
	    
			typename vector<BinarySearchTree<tuple<T> >*>::iterator p = BST.begin();
			p = p + next_clus;
			BST.erase(p);
	    
			vector<int>::iterator q  = repMap.begin();
			q  = q + next_clus;
			repMap.erase(q);

			removeFromHash(k);
			return;
	  	}
  
		/*
		If the representative element of the BST is k, then:
		1. remove k from hash
		2. set corresponding element in repMap to the new representative
		3. insert the new representative in hash
		*/
  
		if(repMap[next_clus] == k) 
		{
			removeFromHash(repMap[next_clus]);
			repMap[next_clus] = (BST[next_clus]->findMax()).id;
			insertInHash(repMap[next_clus]);
		}
  
		next_clus = mergeBST(next_clus);
		if(next_clus >= 0)
			splitBST(next_clus);
	  
		return;
	}


	/*Return the predecessor of an element*/

	tuple<T> predecessor(int x)
  	{
		int pre_rep = get_pre_rep(x);  
  		int next_clus = get_cluster(pre_rep, x);
  
  		if(next_clus == (int)BST.size())
		{   
			tuple<T> temp(pre_rep, NULL);
			return BST[next_clus-1]->find(temp);
		}

  
		tuple<T> item(x, NULL);
		tuple<T> pre = BST[next_clus]->predecessor(item);
  
		if(pre.id == -1 && next_clus > 0) 
		{
			tuple<T> pretup(pre_rep, NULL);
			return BST[next_clus-1]->find(pretup);
		}
  
		return pre;
	}


	/*Reports all integers between 'a' and 'b' (inclusive)*/
	
	vector<tuple<T> > report(int a, int b)
	{
//		cout<<"Inside RR.report"<<endl;
		  
		vector<tuple<T> > x;

//		cout<<n;
//		int abc;
//		scanf("%d", &abc);

		if(n == 0)
		{ 
//			cout<<"returning"<<endl;
			return x;
		}
		tuple<T> pre = predecessor(b);
	  
		//cout << "a = " << a << ", b = " << b << endl;
	  
		while(pre.id >= a) 
		{
			//cout << "pre.id = " << pre.id << ", pre = " << pre << endl;

			if ( pre.ptr != NULL ) 
				x.push_back(pre);
		
			b = pre.id - 1;

			if(b < 0)
				break;

			pre = predecessor(b);
		}

//		cout<<"returning with "<<(int)x.size()<<endl;
		return x;
	}


	/*Return Size of S*/
	int size() 
	{
		return n;
	}


	inline void printAll()
	{
		printBSTs();
		printReps();
		return;
	}
};






#endif

