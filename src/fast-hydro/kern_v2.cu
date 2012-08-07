#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
//#include "dynArray_cpu.h"

#include "../../inc/pseudoGsol.h"
#include "pseudoGsol_cuda.h"

vector<ball_cuda> convertGrid(PG *grid_);
__device__ void printGrid(grid_cuda g);
__device__ grid_cuda copyGrid(ball_cuda *balls, int *totBalls);


__device__ void coisas(grid_cuda g) //OK
{

	for (int i = 0; i < g.nElems; i++)
			printf("	coisas i g.ptr[i].id g.ptr[i] g.ptr[i].nElems, %d %d %d\n", i, g.ptr[i].id, g.ptr[i].nElems);


	dynArray_cuda<plane_cuda> S2;
	S2 = g.report_d(20,20);

  int size2 = S2.getSize();

printf("CUDA coisas size2 %d\n", size2);


}

__device__ void printGrid(grid_cuda g) //OK
{

		plane_cuda plane;
		line_cuda line;
		gridcell_cuda cell;
		Point_cuda point;

	//	printf("CUDA PGrid nPlanes %d\n", g.nElems);


		for (int i = 0; i < g.nElems; i++)
		{	
		//	printf("i g.ptr[i].id, %d %d\n", i, g.ptr[i].id);
			plane = plane_cuda(g.ptr[i].id, g.ptr[i].nElems, g.ptr[i].ptr);
			printf("CUDA plane ID %d\n", plane.id);
		//	printf("CUDA nLines %d\n", plane.nElems);
				
			for (int j = 0; j < plane.nElems; j++)
			{
		//		printf("j plane.ptr[j].id %d %d\n", j, plane.ptr[j].id);
				line =  line_cuda(plane.ptr[j].id, plane.ptr[j].nElems, plane.ptr[j].ptr);
				printf("	CUDA line ID %d\n", line.id);
		//		printf("CUDA nCells %d\n", line.nElems);
				
				for (int k = 0; k < line.nElems; k++)
				{
				//	printf("k line.ptr[k].id %d %d\n", k, line.ptr[k].id);
					cell = gridcell_cuda(line.ptr[k].id, line.ptr[k].nElems, line.ptr[k].ptr);
					printf("		CUDA cell ID %d\n", cell.id);
					printf("			CUDA nBalls %d\n", cell.nElems);
					
					for (int l = 0; l < cell.nElems; l++)
					{
						point = Point_cuda(cell.ptr[l].x, cell.ptr[l].y, cell.ptr[l].z);
//						printf("CUDA Point xx %f\n", point.x);
//						printf("CUDA Point yy %f\n", point.y);
//						printf("CUDA Point zz %f\n", point.z);
					} // end points
					cell.freePtr();

				} // end cells
				line.freePtr();

			} // end lines
			plane.freePtr();

		}	// end planes

	
}


__device__ grid_cuda copyGrid(ball_cuda *balls, int *totBalls) //OK
{
		int curCellID, curLineID, curPlaneID;

		// to store temporary lists of planes, lines, an cells
		dynArray_cuda<plane_cuda> planes = dynArray_cuda<plane_cuda>(10); 	
		dynArray_cuda<line_cuda> lines;
		dynArray_cuda<gridcell_cuda> cells;
		dynArray_cuda<Point> particles;

		gridcell_cuda cell;
		line_cuda line;
		plane_cuda plane;


		int i = 0;
		curCellID =  balls[i].cell_id;
		curLineID =  balls[i].parentLine;
		curPlaneID =  balls[i].parentPlane;

 	
		while ( i < *totBalls )
		{
//			printf("CUDA curPlaneID %d\n", curPlaneID);
//			printf("CUDA i %d\n", i);

			lines = dynArray_cuda<line_cuda>(5);
			while( balls[i].parentPlane == curPlaneID )		
			{

				cells = dynArray_cuda<gridcell_cuda>(5);
				while ( balls[i].parentLine == curLineID && balls[i].parentPlane == curPlaneID )
				{

					particles = dynArray_cuda<Point>(10);		
					while ( balls[i].cell_id == curCellID && balls[i].parentLine == curLineID && balls[i].parentPlane == curPlaneID )
					{
						particles.push_end(balls[i].ball);
						i++;	
					}	// end cells		

	//				printf("CUDA cell id %d\n", curCellID);
//					printf("CUDA nBalls %d\n", particles.getSize());
					
					cell = gridcell_cuda(curCellID, particles.getSize(), particles.getAll());
					cells.push_end(cell);
					particles.freeElements();
					curCellID =  balls[i].cell_id;
				} // end lines

//				printf("CUDA line id %d\n", curLineID);
//				printf("CUDA nCells %d\n", cells.getSize());

				line = line_cuda(curLineID, cells.getSize(), cells.getAll());
				lines.push_end(line);
				cells.freeElements();
				curLineID =  balls[i].parentLine;
			} // end planes

//			printf("CUDA plane id %d\n", curPlaneID);
//			printf("CUDA nLines %d\n", lines.getSize());

			plane = plane_cuda(curPlaneID, lines.getSize(), lines.getAll());
		//	printf("curPlaneID plane id %d %d\n", curPlaneID, plane.id);
			planes.push_end(plane);		
			lines.freeElements();
			curPlaneID =  balls[i].parentPlane;
		} // end all
	
		grid_cuda g = grid_cuda(planes.getSize(), planes.getAll());
/*
		for (int i = 0; i < planes.getSize(); i++)
			printf(" i, planes.get(i), planes.get(i).id %d %d %d\n", i, planes.get(i), planes.get(i).id);				
		

		//printf("CUDA nPlanes %d\n", planes.getSize());


		for (int i = 0; i < g.nElems; i++)
			printf(" antes i plane.id %d %d\n", i, g.ptr[i].id);
		
			
		for (int i = 0; i < g.nElems; i++)
			printf("depois i plane.id %d %d\n", i,g.ptr[i].id);
*/

		planes.freeElements();
		return g;

}


	





__device__ float distsq_cuda(Point_cuda *a, Point_cuda *b) { //should be OK
		float res;
		double dx, dy, dz;
//		printf("OLAAA\n");
		dx = a->x - b->x;
		dy = a->y - b->y;
		dz = a->z - b->z;  
		res = dx*dx + dy*dy + dz*dz;
//		printf("res %lf\n", res);
		return res;
}



__device__ void transformPoint_cuda( Point_cuda *p, double *transMat, Point_cuda *np ) //OK
{
   np->x = transMat[  0 ] * p->x + transMat[  1 ] * p->y + transMat[  2 ] * p->z + transMat[  3 ];
   np->y = transMat[  4 ] * p->x + transMat[  5 ] * p->y + transMat[  6 ] * p->z + transMat[  7 ];
   np->z = transMat[  8 ] * p->x + transMat[  9 ] * p->y + transMat[ 10 ] * p->z + transMat[ 11 ];       
}



__device__ bool pointsWithinRange_cuda(Point_cuda *q, double *delta, double *TRANSLATE, double *DIM, int *rangeCount, grid_cuda *g)
{

	//hist[threadIdx.x] += 1;
//		printf("INSIDE distanceCutoff TRANSLATE DIM rangeCount %lf %lf %lf %d\n", *delta, *TRANSLATE, *DIM, *rangeCount);
  Point_cuda p;
	p.x = q->x;
	p.y = q->y;
	p.z = q->z;

/*printf("threadIdx.x %d\n",threadIdx.x);
if (threadIdx.x == 217){
        printf("A p.x %lf\n", p.x);  
        printf("A p.y %lf\n", p.y);
        printf("A p.z %lf\n", p.z);
}*/
  p.x += *TRANSLATE;
  p.y += *TRANSLATE;
  p.z += *TRANSLATE;
/*if (threadIdx.x == 217){
        printf("D p.x %lf\n", p.x);  
        printf("D p.y %lf\n", p.y);
        printf("D p.z %lf\n", p.z);
}*/
  int l = (int)((p.z - (*delta)) / (*DIM));
  int h = (int)((p.z + (*delta)) / (*DIM));



  int size,size1,size2;
  
  dynArray_cuda<line_cuda>  temp1;
  dynArray_cuda<gridcell_cuda> temp2;


  //int tts;
//	printf("CUDA l, h %d %d\n", l, h);
	dynArray_cuda<plane_cuda> S2;
	S2 = g->report_d(l,h);

  (*rangeCount)++;

  if(S2.empty())
  {
		//printf("CUDA s2 empty \n");
    return false;
  }
  size2 = S2.getSize();

//printf("CUDA  l h size2 %d %d %d \n",l, h, size2);


//	printf("CUDA size2 %d\n", size2);

  for(int i = 0; i < size2; i++) 
  {
    (*rangeCount)++;
   
    l = (int)((p.y - (*delta))/ (*DIM));
    h = (int)((p.y + (*delta))/ (*DIM));
    temp1 = (S2.get(i)).report_d(l,h);
    
    size1 = temp1.getSize();



//    printf("CUDA size1 %d i %d\n", size1, i);
//		hist[threadIdx.x] += size1;
    for(int j = 0; j < size1; j++) 
    {
      (*rangeCount)++;

      l = (int)((p.x - (*delta))/ (*DIM));
      h = (int)((p.x + (*delta))/ (*DIM));


      if(!(temp1.get(j).empty()))
      {
				temp2 = (temp1.get(j)).report_d(l,h);
						  	
				size = temp2.getSize();
//				int pbcount;
			//	printf("CUDA  l h size %d %d %d \n",l, h, size);

				Point_cuda *oa = (Point_cuda*) malloc (sizeof (Point_cuda));
		 	  //printf("CUDA S2.get(i).id, temp1.get(j).id l h size i %d %d %d %d %d\n", S2.get(i).id, temp1.get(j).id, l, h, size, i);

				double delsq = (*delta)*(*delta);
				for(int mn=0;mn<size;mn++)
				{
					int atomsInCell = temp2.get(mn).getSize();
	//				printf("atomsInCell %d\n", atomsInCell);
					for(int pq=0;pq<atomsInCell;pq++)
					{
		//				printf("OLA111111111111111111111\n");
						//printf("temp2.get(mn).ptr[pq].x %lf\n", temp2.get(mn).ptr[pq].x);
						oa->x = temp2.get(mn).ptr[pq].x;
						oa->y = temp2.get(mn).ptr[pq].y;
						oa->z = temp2.get(mn).ptr[pq].z;
		//				printf("OLA222222222222222222222\n");
//						printf("delsq %lf\n", delsq);
						float dist = (distsq_cuda(oa, q));
		
						
						if(distsq_cuda(oa, q) <= delsq){ 
				//			printf("CUDA thread dist delsq %d %f %lf\n", threadIdx.x, dist, delsq);
							//printf("TRUE\n");
							return true;
			//				printf("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
						}
					}
				}
				free(oa);
	
      	temp2.freeElements();
	
      }
    }
    temp1.freeElements();
  }
	S2.freeElements();

  return false;
}







//__device__ void pseudoGsol_static(int *offset, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes, grid_cuda *movingPG, double *distanceCutoff, double *TRANSLATE, double *DIM, int *rangeCount)

__device__ void pseudoGsolStatic(grid_cuda *movingPG, int *offset, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATE, double *DIM, int *rangeCount, double *distanceCutoff)//, double *pGsol_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)
{


	int i = threadIdx.x;

	if (i >= 0 && i < *numStaticQPointsOctreeNodes)
	{
		pGsolStatic[i] = 0;
		pGsolHStaticPos[i] = 0;
		pGsolHStaticNeg[i] = 0;
		
		if ( staticQPointsOctreeFlags[ *offset + i ])
		{
			for (int j = staticQPointsOctree[ i ].qPtsStartID; j <= staticQPointsOctree[ i ].qPtsEndID; j++ )
			{	

				Point_cuda p, q;
               
        p.x = staticQPoints[ j ].x;
        p.y = staticQPoints[ j ].y;
        p.z = staticQPoints[ j ].z;

   	    transformPoint_cuda( &p, transI, &q );

    	   if ( pointsWithinRange_cuda( &q, distanceCutoff, TRANSLATE, DIM, rangeCount, movingPG) )
				{
        	( pGsolStatic[i] ) += staticQPoints[ j ].w;

          if ( staticQPoints[ j ].h > 0 ) 
						( pGsolHStaticPos[i] ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
          else 
						( pGsolHStaticNeg[i] ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
       }  
      }            
		}
	}



}


__device__ void pseudoGsolMoving(int *offset, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes, grid_cuda *staticPG, double *distanceCutoff, double *TRANSLATE, double *DIM, int *rangeCount)
{



	int i = threadIdx.x;


	if (i >= 0 && i < *numMovingQPointsOctreeNodes)
	{

		pGsolMoving[i] = 0;
		pGsolHMovingPos[i] = 0;
		pGsolHMovingNeg[i] = 0;

		if ( movingQPointsOctreeFlags[ *offset + i ])
		{
			for (int j = movingQPointsOctree[ i ].qPtsStartID;  j <= movingQPointsOctree[ i ].qPtsEndID; j++ )
			{
				Point_cuda p, q;
               
        p.x = movingQPoints[ j ].x;
        p.y = movingQPoints[ j ].y;
        p.z = movingQPoints[ j ].z;
               
        transformPoint_cuda( &p, trans, &q );
           
        if ( pointsWithinRange_cuda( &q, distanceCutoff, TRANSLATE, DIM, rangeCount, staticPG ) )
        {
        	( pGsolMoving[i] ) += movingQPoints[ j ].w;

          if ( movingQPoints[ j ].h > 0 ) 
						( pGsolHMovingPos[i] ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
          else 
						( pGsolHMovingNeg[i] ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
        }  
     }            
	
		}

	}

}





/*

__device__ void block_sum(double *input, double *results, int n)
{
//	printf("OLAAAAAA11111\n");
	extern __shared__ double sdata[];
	int i = threadIdx.x, tx = threadIdx.x;

	double x = 0;

	if (i < n)
		x = input[i];

	sdata[tx] = x;
	//printf("sdata[tx] %lf %d\n", sdata[tx], tx);
	__syncthreads();

	for (int offset = blockDim.x / 2; offset > 0; offset >>= 1)
	{
		if (tx < offset)
			sdata[tx] += sdata[tx+offset];

		__syncthreads;
	}

//	__shared__ double *b;
	if (threadIdx.x == 0){
		//*results = sdata[0];
		printf("OLAAAAAA %lf\n", sdata[0]);//, *b);
	}
//	*results = sdata[0];
//	printf("OLAAAAAA %lf\n", *results);
	//results[blockIdx.x] = sdata[0];

}
*/


__global__ void kern1(/*static stuff*/ball_cuda *ballsStatic, int *totBallsStatic, int *offsetStatic, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATEstatic, double *DIMstatic, int *rangeCountStatic, double *distanceCutoff /*moving stuff*/, ball_cuda *ballsMoving, int *totBallsMoving, int *offsetMoving, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes,  double *TRANSLATEmoving, double *DIMmoving, int *rangeCountMoving /*others*/)//, double *pGsol_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)
//__global__ void  kern(struct static_params *statP)
{
		//printf("ola\n");



			__shared__ grid_cuda movingGrid;

			if (threadIdx.x == 0)
				movingGrid = copyGrid(ballsStatic, totBallsStatic);

			__syncthreads();

			pseudoGsolStatic(&movingGrid, offsetStatic, transI, pGsolStatic, pGsolHStaticPos, pGsolHStaticNeg, staticQPointsOctreeFlags, staticQPointsOctree, staticQPoints, numStaticQPointsOctreeNodes, TRANSLATEstatic, DIMstatic, rangeCountStatic, distanceCutoff);




		}


__global__ void kern2(/*static stuff*/ball_cuda *ballsStatic, int *totBallsStatic, int *offsetStatic, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATEstatic, double *DIMstatic, int *rangeCountStatic, double *distanceCutoff /*moving stuff*/, ball_cuda *ballsMoving, int *totBallsMoving, int *offsetMoving, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes,  double *TRANSLATEmoving, double *DIMmoving, int *rangeCountMoving /*others*/)//, double *pGsol_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)
//__global__ void  kern(struct static_params *statP)
{
		//printf("ola\n");



			__shared__ grid_cuda staticGrid;

			if (threadIdx.x == 0)
				staticGrid = copyGrid(ballsMoving, totBallsMoving);

			__syncthreads();

			pseudoGsolMoving(offsetMoving, trans, pGsolMoving, pGsolHMovingPos, pGsolHMovingNeg, movingQPointsOctreeFlags, movingQPointsOctree, movingQPoints, numMovingQPointsOctreeNodes, &staticGrid, distanceCutoff, TRANSLATEmoving, DIMmoving, rangeCountMoving);

}

/*




double sum(double *d_input, int n, int block_size, int num_blocks)
{
	double *d_sums = 0;
	cudaMalloc((void**)&d_sums, sizeof(double));// * (num_Blocks+1));
	
	//int smem_sz = block_size* sizeof(double));

//	block_sum<<<num_blocks, block_size>>>(d_input);//, d_sums, n);
	block_sum<<<1,1>>>();//d_input);//, d_sums, n);
//	block_sum<<<1, block_size, smem_sz>>>(d_sums, dsums+num_blocks, num_blocks);

	double result = 0;
	cudaMemcpy(&result, d_sums, sizeof(double),cudaMemcpyDeviceToHost);
	
	return result;
}


*/
int pseudoGsol::call_to_kern(int offsetStatic, double *transI, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE *staticQPointsOctree, QPOINT *staticQPoints, int numStaticQPoints, int numStaticQPointsOctreeNodes, PG *movingPG, double TRANSLATEstatic, double DIMstatic, int rangeCountStatic, PARAMS_IN params, int offsetMoving, double *trans, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE *movingQPointsOctree, QPOINT *movingQPoints, int numMovingQPoints, int numMovingQPointsOctreeNodes, PG *staticPG, double TRANSLATEmoving, double DIMmoving, int rangeCountMoving)
//int pseudoGsol::call_to_kern( PG *movingPG)
{

// convert movingGrid

		vector<ball_cuda> ballsVecMoving_h = convertGrid(movingPG);

		int totBallsMoving = ballsVecMoving_h.size();
		int *totBallsMoving_d;

		int totSizeMoving = totBallsMoving * (sizeof(int)*3+sizeof(Point));

		ball_cuda *arrayMoving_h = (ball_cuda*) malloc(totSizeMoving);
		std::copy(ballsVecMoving_h.begin(), ballsVecMoving_h.end(), arrayMoving_h);
		ball_cuda *arrayMoving_d;


// convert staticGrid

		vector<ball_cuda> ballsVecStatic_h = convertGrid(staticPG);

		int totBallsStatic = ballsVecStatic_h.size();
		int *totBallsStatic_d;

		int totSizeStatic = totBallsStatic * (sizeof(int)*3+sizeof(Point));

		ball_cuda *arrayStatic_h = (ball_cuda*) malloc(totSizeStatic);
		std::copy(ballsVecStatic_h.begin(), ballsVecStatic_h.end(), arrayStatic_h);
		ball_cuda *arrayStatic_d;



// passing structs


		struct QPOINT_CUDA *staticQPoints_h = (QPOINT_CUDA*) malloc(numStaticQPoints * sizeof(QPOINT_CUDA));
		struct QPOINT_CUDA *staticQPoints_d;

		for (int i = 0; i < numStaticQPoints; i++)
		{
			staticQPoints_h[i].x = staticQPoints[i].x;
			staticQPoints_h[i].y = staticQPoints[i].y;
			staticQPoints_h[i].z = staticQPoints[i].z;
			staticQPoints_h[i].w = staticQPoints[i].w;
			staticQPoints_h[i].h = staticQPoints[i].h;
		}			
		



		struct QPOINT_CUDA *movingQPoints_h = (QPOINT_CUDA*) malloc(numMovingQPoints * sizeof(QPOINT_CUDA));
		struct QPOINT_CUDA *movingQPoints_d;

		for (int i = 0; i < numMovingQPoints; i++)
		{
			movingQPoints_h[i].x = movingQPoints[i].x;
			movingQPoints_h[i].y = movingQPoints[i].y;
			movingQPoints_h[i].z = movingQPoints[i].z;
			movingQPoints_h[i].w = movingQPoints[i].w;
			movingQPoints_h[i].h = movingQPoints[i].h;
		}
		

// -------------------------

		struct QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree_h = (QPOINTS_OCTREE_NODE_CUDA*) malloc(numStaticQPointsOctreeNodes * sizeof(QPOINTS_OCTREE_NODE_CUDA));
		struct QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree_d;

		for (int i = 0; i < numStaticQPointsOctreeNodes; i++)
		{
			staticQPointsOctree_h[i].qPtsStartID = staticQPointsOctree[i].qPtsStartID;
			staticQPointsOctree_h[i].qPtsEndID = staticQPointsOctree[i].qPtsEndID;
		}



		struct QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree_h = (QPOINTS_OCTREE_NODE_CUDA*) malloc(numMovingQPointsOctreeNodes * sizeof(QPOINTS_OCTREE_NODE_CUDA));
		struct QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree_d;

		for (int i = 0; i < numMovingQPointsOctreeNodes; i++)
		{
			movingQPointsOctree_h[i].qPtsStartID = movingQPointsOctree[i].qPtsStartID;
			movingQPointsOctree_h[i].qPtsEndID = movingQPointsOctree[i].qPtsEndID;
		}
		

/*
cudaError_t cudaMemcpyAsync 	( 	void *  	dst,
		const void *  	src,
		size_t  	count,
		enum cudaMemcpyKind  	kind,
		cudaStream_t  	stream = 0	 
	) 
*/


		bool *staticQPointsOctreeFlags_d, *movingQPointsOctreeFlags_d;
		double *trans_d, *transI_d;

		double *pGsolStatic_h, *pGsolStatic_d, *pGsolMoving_h, *pGsolMoving_d;
		pGsolStatic_h = (double*) malloc(sizeof(double) * numStaticQPointsOctreeNodes);
		pGsolMoving_h = (double*) malloc(sizeof(double) * numMovingQPointsOctreeNodes);

		double *pGsolHStaticPos_h, *pGsolHStaticPos_d, *pGsolHStaticNeg_h, *pGsolHStaticNeg_d;
		pGsolHStaticPos_h = (double*) malloc(sizeof(double) * numStaticQPointsOctreeNodes);
		pGsolHStaticNeg_h = (double*) malloc(sizeof(double) * numStaticQPointsOctreeNodes);

		double *pGsolHMovingPos_h, *pGsolHMovingPos_d, *pGsolHMovingNeg_h, *pGsolHMovingNeg_d;
		pGsolHMovingPos_h = (double*) malloc(sizeof(double) * numMovingQPointsOctreeNodes);
		pGsolHMovingNeg_h = (double*) malloc(sizeof(double) * numMovingQPointsOctreeNodes);


		struct static_params statS_h, *statS_d;
		
		statS_h.offset = offsetStatic;
		statS_h.numStaticQPointsOctreeNodes = numStaticQPointsOctreeNodes;
		statS_h.distanceCutoff = params.distanceCutoff;
		statS_h.TRANSLATE = TRANSLATEstatic;
		statS_h.DIM = DIMstatic;
		statS_h.rangeCount = rangeCountStatic;


		struct moving_params statM_h, *statM_d;
		
		statM_h.offset = offsetMoving;
		statM_h.numMovingQPointsOctreeNodes = numMovingQPointsOctreeNodes;
		statM_h.distanceCutoff = params.distanceCutoff;
		statM_h.TRANSLATE = TRANSLATEmoving;
		statM_h.DIM = DIMmoving;
		statM_h.rangeCount = rangeCountMoving;

// pG values-----------------------

		//*pGsol = *pGsolHStaticPos = *pGsolHStaticNeg = *pGsolHMovingPos = *pGsolHMovingNeg = 0;

	//	double *pGsol_scalar_d, *pGsolHStaticPos_scalar_d, *pGsolHStaticNeg_scalar_d, *pGsolHMovingPos_scalar_d, *pGsolHMovingNeg_scalar_d;




// cudaMalloc ---------------------------

		cudaDeviceSetLimit (cudaLimitMallocHeapSize, 512*1024*1024);

		cudaMalloc( (void **)&arrayMoving_d, totSizeMoving );
		cudaMalloc( (void **)&totBallsMoving_d, sizeof(int) );

		cudaMalloc( (void **)&arrayStatic_d, totSizeStatic );
		cudaMalloc( (void **)&totBallsStatic_d, sizeof(int) );


		cudaMalloc( (void **)&staticQPoints_d, sizeof(QPOINT_CUDA)*numStaticQPoints );
		cudaMalloc( (void **)&movingQPoints_d, sizeof(QPOINT_CUDA)*numMovingQPoints );

		cudaMalloc( (void **)&staticQPointsOctree_d, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numStaticQPointsOctreeNodes );
		cudaMalloc( (void **)&movingQPointsOctree_d, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numMovingQPointsOctreeNodes );

		cudaMalloc( (void **)&staticQPointsOctreeFlags_d, sizeof(bool)*numStaticQPointsOctreeNodes );
		cudaMalloc( (void **)&movingQPointsOctreeFlags_d, sizeof(bool)*numMovingQPointsOctreeNodes );

		cudaMalloc( (void **)&trans_d, sizeof(double)*12 );
		cudaMalloc( (void **)&transI_d, sizeof(double)*12 );

		cudaMalloc( (void **)&statS_d, 3*(sizeof(int)+sizeof(double)) );//sizeof(static_params) );
		cudaMalloc( (void **)&statM_d, 3*(sizeof(int)+sizeof(double)) );

		
		cudaMalloc( (void **)&pGsolStatic_d, sizeof(double) * numStaticQPointsOctreeNodes);
		cudaMalloc( (void **)&pGsolMoving_d, sizeof(double) * numMovingQPointsOctreeNodes);

		cudaMalloc( (void **)&pGsolHStaticPos_d, sizeof(double) * numStaticQPointsOctreeNodes);
		cudaMalloc( (void **)&pGsolHStaticNeg_d, sizeof(double) * numStaticQPointsOctreeNodes);

		cudaMalloc( (void **)&pGsolHMovingPos_d, sizeof(double) * numMovingQPointsOctreeNodes);
		cudaMalloc( (void **)&pGsolHMovingNeg_d, sizeof(double) * numMovingQPointsOctreeNodes);


	//	cudaMalloc( (void **)&pGsol_scalar_d, sizeof(double));

	//	cudaMalloc( (void **)&pGsolHStaticPos_scalar_d, sizeof(double));
	//	cudaMalloc( (void **)&pGsolHStaticNeg_scalar_d, sizeof(double));

	//	cudaMalloc( (void **)&pGsolHMovingPos_scalar_d, sizeof(double));
	//	cudaMalloc( (void **)&pGsolHMovingNeg_scalar_d, sizeof(double));


		//int *hist_h = (int*) malloc(sizeof(int)*numStaticQPointsOctreeNodes);
//		int *hist_d;
		//cudaMalloc( (void **)&hist_d, sizeof(int) * numStaticQPointsOctreeNodes);

// cudaMemcpy-------------------

		cudaMemcpy( arrayMoving_d, arrayMoving_h, totSizeMoving, cudaMemcpyHostToDevice );
		cudaMemcpy( totBallsMoving_d, &totBallsMoving, sizeof(int), cudaMemcpyHostToDevice );

		cudaMemcpy( arrayStatic_d, arrayStatic_h, totSizeStatic, cudaMemcpyHostToDevice );
		cudaMemcpy( totBallsStatic_d, &totBallsStatic, sizeof(int), cudaMemcpyHostToDevice );

		cudaMemcpy( staticQPoints_d, staticQPoints_h, sizeof(QPOINT_CUDA)*numStaticQPoints, cudaMemcpyHostToDevice );
		cudaMemcpy( movingQPoints_d, movingQPoints_h, sizeof(QPOINT_CUDA)*numMovingQPoints, cudaMemcpyHostToDevice );

		cudaMemcpy( staticQPointsOctree_d, staticQPointsOctree_h, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numStaticQPointsOctreeNodes, cudaMemcpyHostToDevice );
		cudaMemcpy( movingQPointsOctree_d, movingQPointsOctree_h, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numMovingQPointsOctreeNodes, cudaMemcpyHostToDevice );

	

		cudaMemcpy( staticQPointsOctreeFlags_d, staticQPointsOctreeFlags, sizeof(bool)*numStaticQPointsOctreeNodes, cudaMemcpyHostToDevice );
		cudaMemcpy( movingQPointsOctreeFlags_d, movingQPointsOctreeFlags, sizeof(bool)*numMovingQPointsOctreeNodes, cudaMemcpyHostToDevice );


		cudaMemcpy( trans_d, trans, sizeof(double)*12, cudaMemcpyHostToDevice );
		cudaMemcpy( transI_d, transI, sizeof(double)*12, cudaMemcpyHostToDevice );


		cudaMemcpy( statS_d, &statS_h, 3*(sizeof(int)+sizeof(double)), cudaMemcpyHostToDevice );
		cudaMemcpy( statM_d, &statM_h, 3*(sizeof(int)+sizeof(double)), cudaMemcpyHostToDevice );


// call Kernel ---------------------

		int numBlocks = 1;
		int threadsPerBlock;

		threadsPerBlock = numStaticQPointsOctreeNodes;

//		int threadsPerBlock = 1;

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error1: %s\n", cudaGetErrorString(err));

		printf("calling... %d\n", threadsPerBlock);

		kern1<<< numBlocks,threadsPerBlock >>>(arrayMoving_d, totBallsMoving_d, &(statS_d->offset), transI_d, pGsolStatic_d, pGsolHStaticPos_d, pGsolHStaticNeg_d, staticQPointsOctreeFlags_d, staticQPointsOctree_d, staticQPoints_d, &(statS_d->numStaticQPointsOctreeNodes), &(statS_d->TRANSLATE), &(statS_d->DIM), &(statS_d->rangeCount), &(statS_d->distanceCutoff), arrayStatic_d, totBallsStatic_d, &(statM_d->offset), trans_d, pGsolMoving_d, pGsolHMovingPos_d, pGsolHMovingNeg_d, movingQPointsOctreeFlags_d, movingQPointsOctree_d, movingQPoints_d, &(statM_d->numMovingQPointsOctreeNodes), &(statM_d->TRANSLATE), &(statM_d->DIM), &(statM_d->rangeCount));

	threadsPerBlock = numMovingQPointsOctreeNodes;

	kern2<<< numBlocks,threadsPerBlock >>>(arrayMoving_d, totBallsMoving_d, &(statS_d->offset), transI_d, pGsolStatic_d, pGsolHStaticPos_d, pGsolHStaticNeg_d, staticQPointsOctreeFlags_d, staticQPointsOctree_d, staticQPoints_d, &(statS_d->numStaticQPointsOctreeNodes), &(statS_d->TRANSLATE), &(statS_d->DIM), &(statS_d->rangeCount), &(statS_d->distanceCutoff), arrayStatic_d, totBallsStatic_d, &(statM_d->offset), trans_d, pGsolMoving_d, pGsolHMovingPos_d, pGsolHMovingNeg_d, movingQPointsOctreeFlags_d, movingQPointsOctree_d, movingQPoints_d, &(statM_d->numMovingQPointsOctreeNodes), &(statM_d->TRANSLATE), &(statM_d->DIM), &(statM_d->rangeCount));

//(/*static stuff*/ball_cuda *ballsStatic, int *totBallsStatic, int *offsetStatic, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATEstatic, double *DIMstatic, int *rangeCountStatic, double *distanceCutoff /*moving stuff*/, ball_cuda *ballsMoving, int *totBallsMoving, int *offsetMoving, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes,  double *TRANSLATEmoving, double *DIMmoving, int *rangeCountMoving /*others*/, double *pGsol_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)
//__global__ void  kern(struct static_params *statP)


//, pGsol_scalar_d, pGsolHStaticPos_scalar_d, pGsolHStaticNeg_scalar_d, pGsolHMovingPos_scalar_d, pGsolHMovingNeg_scalar_d);



		printf("end...\n");

		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error2: %s\n", cudaGetErrorString(err));


		for (int i = 0; i < numStaticQPointsOctreeNodes; i++){
			pGsolStatic_h[i] = 0;
			pGsolHStaticPos_h[i] = 0;
			pGsolHStaticNeg_h[i] = 0;
		}


		for (int i = 0; i < numMovingQPointsOctreeNodes; i++){
			pGsolMoving_h[i] = 0;
			pGsolHMovingPos_h[i] = 0;
			pGsolHMovingNeg_h[i] = 0;
		}


		cudaMemcpy( pGsolStatic_h, pGsolStatic_d, sizeof(double)* numStaticQPointsOctreeNodes, cudaMemcpyDeviceToHost );
		cudaMemcpy( pGsolHStaticPos_h, pGsolHStaticPos_d, sizeof(double)* numStaticQPointsOctreeNodes, cudaMemcpyDeviceToHost );
		cudaMemcpy( pGsolHStaticNeg_h, pGsolHStaticNeg_d, sizeof(double)* numStaticQPointsOctreeNodes, cudaMemcpyDeviceToHost );	


		cudaMemcpy( pGsolMoving_h, pGsolMoving_d, sizeof(double)* numMovingQPointsOctreeNodes, cudaMemcpyDeviceToHost );
		cudaMemcpy( pGsolHMovingPos_h, pGsolHMovingPos_d, sizeof(double)* numMovingQPointsOctreeNodes, cudaMemcpyDeviceToHost );
		cudaMemcpy( pGsolHMovingNeg_h, pGsolHMovingNeg_d, sizeof(double)* numMovingQPointsOctreeNodes, cudaMemcpyDeviceToHost );	


		//cudaMemcpy( hist_h, hist_d, sizeof(int)*numStaticQPointsOctreeNodes, cudaMemcpyDeviceToHost );


		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error3: %s\n", cudaGetErrorString(err));


// do reductions 
		double sPgSt = 0;
		double sPosSt = 0;
		double sNegSt = 0;

		for (int i = 0; i < numStaticQPointsOctreeNodes; i++){
			sPgSt += pGsolStatic_h[i];
			sPosSt += pGsolHStaticPos_h[i];
			sNegSt += pGsolHStaticNeg_h[i];
		}

		double sPgMov = 0;
		double sPosMov = 0;
		double sNegMov = 0;

		for (int i = 0; i < numMovingQPointsOctreeNodes; i++){
			sPgMov += pGsolMoving_h[i];
			sPosMov += pGsolHMovingPos_h[i];
			sNegMov += pGsolHMovingNeg_h[i];
		}



printf("CUDA *pGsol, *pGsolHStaticPos, *pGsolHStaticNeg, *pGsolHMovingPos, *pGsolHMovingNeg %lf %lf %lf %lf %lf\n", sPgSt+sPgMov, sPosSt, sNegSt, sPosMov, sNegMov);




		cudaFree(arrayMoving_d);
		cudaFree(totBallsMoving_d);
		cudaFree(arrayStatic_d);
		cudaFree(totBallsStatic_d);

		cudaFree(staticQPoints_d);
		cudaFree(movingQPoints_d);

		cudaFree(staticQPointsOctree_d);
		cudaFree(movingQPointsOctree_d);

		cudaFree(staticQPointsOctreeFlags_d);
		cudaFree(movingQPointsOctreeFlags_d);

		cudaFree(trans_d);
		cudaFree(transI_d);

		cudaFree(statS_d);
		cudaFree(statM_d);

		cudaFree(pGsolStatic_d);
		cudaFree(pGsolMoving_d);

		cudaFree(pGsolHStaticPos_d);
		cudaFree(pGsolHStaticNeg_d);

		cudaFree(pGsolHMovingPos_d);
		cudaFree(pGsolHMovingNeg_d);


		free(arrayMoving_h);
		free(arrayStatic_h);

		free(pGsolStatic_h);
		free(pGsolMoving_h);

		free(pGsolHStaticPos_h);
		free(pGsolHStaticNeg_h);
		free(pGsolHMovingPos_h);
		free(pGsolHMovingNeg_h);

		free(staticQPointsOctree_h);
		free(staticQPoints_h);

		free(movingQPointsOctree_h);
		free(movingQPoints_h);

		cudaDeviceReset();

	return 0;
}




vector<ball_cuda> convertGrid(PG *grid_)
{

		int nPlanes;
		int nLines;
		int nCells;
		int nBalls;

		int grid_init, grid_end;	
		int plane_init, plane_end;
		int line_init, line_end;
		
		vector<tuple<plane*> > grid1;

		tuple<plane*> plane1;
		vector<tuple<line*> > plane2;

		tuple<line*> line1;
		vector<tuple<gridcell*> > line2;
		
		tuple<gridcell*> cell1;
		
		vector<ball_cuda> ballsVec_h;
		ball_cuda tempBall_h;

		grid g1 = grid_->getGrid();

		nPlanes = g1.RR.getn();


		grid_init = g1.RR.getOverallMin();
		grid_end = g1.RR.getOverallMax();	
		grid1 = g1.RR.report(grid_init, grid_end);

		
		// goes through all planes
		for (int i = 0; i < nPlanes; i++)
		{
			plane1 = grid1.at(i);
	//		printf("plane ID %d\n", plane1.id);

			nLines = (plane1.ptr)->RR.getn();
//			printf(" nLines %d\n", nLines);
			
			plane_init =(plane1.ptr)->RR.getOverallMin(); 
			plane_end = (plane1.ptr)->RR.getOverallMax();

			plane2 = (plane1.ptr)->RR.report(plane_init, plane_end);

			// goes through all lines
			for (int j = 0; j < nLines; j++)
			{
				line1 = plane2.at(j);
	//			printf("	line ID %d\n", line1.id);

				nCells = (line1.ptr)->RR.getn();
		//		printf("nCells  %d\n", nCells);

				line_init = (line1.ptr)->RR.getOverallMin(); 
				line_end = (line1.ptr)->RR.getOverallMax(); 
				line2 = (line1.ptr)->RR.report(line_init, line_end);

				// goes through all cells
				for (int k = 0; k < nCells; k++)
				{
					cell1 = line2.at(k);
	//				printf("		cell ID %d\n", cell1.id);
					

					nBalls = cell1.ptr->balls.size();
	//				printf("			nBalls %d\n", nBalls);

					// create vector with total number of balls and plane, line, and cell ID
					//goes through all balls
					for (int kk = 0 ; kk < nBalls; kk++)
					{
					
						tempBall_h.cell_id = cell1.id;
						tempBall_h.parentLine = line1.id;
						tempBall_h.parentPlane = plane1.id;
						tempBall_h.ball = Point (cell1.ptr->balls.at(kk)->x, cell1.ptr->balls.at(kk)->y, cell1.ptr->balls.at(kk)->z);

						ballsVec_h.push_back(tempBall_h);

					} //end balls	
				} // end cells
			} // end lines
		} // end planes
		
		return ballsVec_h;

}
