#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>

#include "dynArray.h"
#include "../../inc/pseudoGsol.h"

/* Point_cuda is the same as the struct Point in the CPU but allows one to create new points inside the GPU */

struct Point_cuda {
  float x;
  float y;
  float z;
  
  __device__ Point_cuda() {}

  __device__ Point_cuda(float a, float b, float c) {
    x = a;
    y = b;
    z = c;
  }

};


/* ball_cuda consists on a Point with its parent plane, line, and gridcell IDs 'attached' so that the grid can be reconstructed on the GPU */

struct ball_cuda
{
	int cell_id;
	int parentLine;
	int parentPlane;

	Point ball;
};



/* gridcell_cuda for GPU contains only its id, the number of Points that it contains, and pointer to those Points */

class gridcell_cuda
{
 public:
	int id;
	int nElems;
	Point *ptr;

	__device__ gridcell_cuda(){}

	__device__ gridcell_cuda(int id_, int nElems_, Point *ptr_)
	{
		id = id_;
		nElems = nElems_;
		ptr = (Point*) malloc(nElems * sizeof(Point));
		for (int i = 0; i < nElems; i++)
		{
			ptr[i].x = ptr_[i].x;		
			ptr[i].y = ptr_[i].y;		
			ptr[i].z = ptr_[i].z;		
		}
	}	

	__device__ int getID()
	{
		return id;
	}	

	__device__ int getSize()
	{
		return nElems;
	}	

	__device__ int empty()
	{
		if (nElems == 0)
			return true;
		else
			return false;
	}	

	__device__ Point get(int i)
	{
		return ptr[i];
	}	

	__device__ void freePtr()
	{		
		free(ptr);
	}

};


/* line_cuda for GPU contains only its id, the number of gridcells that it contains, and pointer to those gridcells */

class line_cuda 
{ 
 public:
	int id;
	int nElems;
	gridcell_cuda *ptr; // contains gridcells IDs

	__device__ line_cuda(){}

	__device__ line_cuda(int id_, int nElems_, gridcell_cuda *ptr_)
	{
		id = id_;
		//printf("line id %d\n", id);
		nElems = nElems_;
		ptr = (gridcell_cuda*) malloc(nElems * sizeof(gridcell_cuda));
		for (int i = 0; i < nElems; i++)
		{
			ptr[i].id = ptr_[i].id;		
			ptr[i].nElems = ptr_[i].nElems;		
			ptr[i].ptr = ptr_[i].ptr;		
		//	printf("line_cuda i ptr_[i].id ptr[i].id %d %d %d\n", i, ptr_[i].id, ptr[i].id);
		}
	}

	__device__ int getID()
	{
		return id;
	}	

	__device__ int getSize()
	{
		return nElems;
	}	

	__device__ int empty()
	{
		if (nElems == 0)
			return true;
		else
			return false;
	}	

	__device__ gridcell_cuda get(int i)
	{
		return ptr[i];
	}	

	__device__ dynArray_cuda<gridcell_cuda> report_d(int a, int b)
	{
		dynArray_cuda<gridcell_cuda> array = dynArray_cuda<gridcell_cuda>(5);

		for ( int i = 0; i < nElems; i++) // elements are in decreasing order
		{
	//	printf("lines ptr[i].id %d\n", ptr[i].id);
			if ( ptr[i].id < a)
			{
			//		printf("lines break i a b ptr[i].d %d %d %d %d\n", i, a, b, ptr[i].id);
				break;
			}
			if ( ptr[i].id <= b){
				array.push_end(ptr[i]);
		//		printf("lines i ptr[i].i %d\n", i, ptr[i].id);
			}

		}

		return array;
	}

	__device__ void freePtr()
	{		
		free(ptr);
	}

};


/* plane_cuda for GPU contains only its id, the number of lines that it contains, and pointer to those lines */

class plane_cuda 
{
 public:

	int id;
	int nElems;
	line_cuda *ptr; // contains lines IDs

	__device__ plane_cuda(){}

	__device__ plane_cuda(int id_, int nElems_, line_cuda *ptr_)
	{	
		id = id_;
		nElems = nElems_;
		ptr = (line_cuda*) malloc(nElems * sizeof(line_cuda));
		for (int i = 0; i < nElems; i++)
		{
			ptr[i].id = ptr_[i].id;		
			ptr[i].nElems = ptr_[i].nElems;		
			ptr[i].ptr = ptr_[i].ptr;		
		//	printf("plane_cuda i ptr_[i].id ptr[i].id %d %d %d\n", i, ptr_[i].id, ptr[i].id);
		}
	}

	__device__ int getID()
	{
		return id;
	}	

	__device__ int getSize()
	{
		return nElems;
	}	

	__device__ int empty()
	{
		if (nElems == 0)
			return true;
		else
			return false;
	}	

	__device__ line_cuda get(int i)
	{
		return ptr[i];
	}	

	__device__ dynArray_cuda<line_cuda> report_d(int a, int b)
	{
		dynArray_cuda<line_cuda> array = dynArray_cuda<line_cuda>(5);

		for ( int i = 0; i < nElems; i++) // elements are in decreasing order
		{
			//printf("planes ptr[i].d %d\n", ptr[i].id);
			if ( ptr[i].id < a)
			{
			//printf("planes break i a b ptr[i].d %d %d %d %d\n", i, a, b, ptr[i].id);
				break;
			}
			if ( ptr[i].id <= b){
				array.push_end(ptr[i]);
			//printf("planes i ptr[i].i %d\n", i, ptr[i].id);
			}

		}

		return array;
	}

	__device__ void freePtr()
	{		
		free(ptr);
	}


};


/* grid_cuda for GPU contains only  the number of planes that it contains, and pointer to those planes */
  
class grid_cuda 
{
 public:

	int nElems;
	plane_cuda *ptr; // contains planes IDs

	__device__ grid_cuda(){}

	__device__ grid_cuda(int nElems_, plane_cuda *ptr_)
	{
		nElems = nElems_;
		ptr = (plane_cuda*) malloc(nElems * sizeof(plane_cuda));
		for (int i = 0; i < nElems; i++)
		{
			ptr[i].id = ptr_[i].id;		
			//printf("grid_cuda i ptr_[i].id ptr[i].id %d %d %d\n", i, ptr_[i].id, ptr[i].id);
			ptr[i].nElems = ptr_[i].nElems;		
			ptr[i].ptr = ptr_[i].ptr;		
		}
	}

	__device__ void assignValues(int nElems_, plane_cuda *ptr_)
	{
		nElems = nElems_;
		ptr = (plane_cuda*) malloc(nElems * sizeof(plane_cuda));
		for (int i = 0; i < nElems; i++)
		{
			ptr[i].id = ptr_[i].id;		
			//printf("grid_cuda i ptr_[i].id ptr[i].id %d %d %d\n", i, ptr_[i].id, ptr[i].id);
			ptr[i].nElems = ptr_[i].nElems;		
			ptr[i].ptr = ptr_[i].ptr;		
		}
	}

	__device__ int getSize()
	{
		return nElems;
	}	

	__device__ int empty()
	{
		if (nElems == 0)
			return true;
		else
			return false;
	}	

	__device__ plane_cuda get(int i)
	{
		return ptr[i];
	}	

	__device__ dynArray_cuda<plane_cuda> report_d(int a, int b)
	{
		dynArray_cuda<plane_cuda> array = dynArray_cuda<plane_cuda>(5);
		
	//	for (int i = 0; i < nElems; i++)
		//	printf("ptr[i] ptr[i].id %d %d\n", ptr[i], ptr[i].id);

		//printf("grid_cuda a b %d %d\n", a,b);
		for ( int i = 0; i < nElems; i++) // elements are in decreasing order
		{
	//		printf("grid_cuda i ptr[i].id %d %d\n", i, ptr[i].id);
			if ( ptr[i].id < a)
			{
	//			printf("grid_cuda break i a b ptr[i] %d %d %d %d\n", i, a, b, ptr[i].id);
				break;
			}
			if ( ptr[i].id <= b){
				array.push_end(ptr[i]);
		//		printf("grid_cuda push i ptr[i] %d %d\n", i, ptr[i].id);
			}


		}

		return array;
	}

	__device__ void freePtr()
	{		
		free(ptr);
	}

};




//from pseudoGsol.h

/* these structs are only used so that only the members that are used in the GPU have to be copied */

struct QPOINT_CUDA
{
	double x, y, z;
	double w, h;     
};
      
struct QPOINTS_OCTREE_NODE_CUDA
{
	int qPtsStartID, qPtsEndID;
};    


/* these structs were created so that the whole structure can be copied to the GPU at once and not each member separately */

struct static_params
{
	int offset;
	int numStaticQPointsOctreeNodes;
	double distanceCutoff;
	double TRANSLATE;
	double DIM;
	int rangeCount;
};


struct moving_params
{
	int offset;
	int numMovingQPointsOctreeNodes;
	double distanceCutoff;
	double TRANSLATE;
	double DIM;
	int rangeCount;
};


// GPU grids built/used in kern.cu

__device__ grid_cuda movingGrid;
__device__ grid_cuda staticGrid;

// function implemented in kern.cu

vector<ball_cuda> convertGrid(PG *grid_);

__device__ void printGrid(grid_cuda g); 
__global__ void copyStaticGridtoGPU(ball_cuda *balls, int *totBalls);
__global__ void copyMovingGridtoGPU(ball_cuda *balls, int *totBalls);
__device__ grid_cuda copyGrid(ball_cuda *balls, int *totBalls);
__device__ float distsq_cuda(Point_cuda *a, Point_cuda *b);
__device__ void transformPoint_cuda( Point_cuda *p, double *transMat, Point_cuda *np );

__device__ bool pointsWithinRange_cuda(Point_cuda *q, double *delta, double *TRANSLATE, double *DIM, int *rangeCount, grid_cuda *g);

__device__ void pseudoGsolStatic(int *offset, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATE, double *DIM, int *rangeCount, double *distanceCutoff);

__device__ void pseudoGsolMoving(int *offset, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes, double *distanceCutoff, double *TRANSLATE, double *DIM, int *rangeCount);

__global__ void kern1(/*static stuff*/int *offsetStatic, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATEstatic, double *DIMstatic, int *rangeCountStatic, double *distanceCutoff);// /*others*/, double *pGsolStatic_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d)

__global__ void kern2(double *distanceCutoff /*moving stuff*/,int *offsetMoving, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes,  double *TRANSLATEmoving, double *DIMmoving, int *rangeCountMoving);// /*others*/, double *pGsolMoving_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)

__device__ void warpReduce(volatile double *sdata, unsigned int tid);

__device__ void reduce4(double *g_idata, double *g_odata, int n);
