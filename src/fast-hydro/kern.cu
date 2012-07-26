#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

//#include "dynArray_cpu.h"

#include "../../inc/pseudoGsol.h"
#include "pseudoGsol_cuda.h"



__device__ void printGrid(grid_cuda g)
{

		plane_cuda plane;
		line_cuda line;
		gridcell_cuda cell;
		Point_cuda point;

		plane_cuda *planes2;
		line_cuda *lines2;
		gridcell_cuda *cells2;


		printf("CUDA PGrid nPlanes %d\n", g.nElems);
		planes2 = g.report_d(g.ptr[g.nElems-1].id, g.ptr[0].id);
		/*
		for (int i  = 0; i < g.nElems; i++)
			printf("CUDA plane ID %d\n", g.ptr[i].id);
		*/

		for (int i = 0; i < g.nElems; i++)
		{	
			plane = plane_cuda(g.ptr[i].id, g.ptr[i].nElems, g.ptr[i].ptr);
			lines2 = plane.report_d(plane.ptr[plane.nElems-1].id, plane.ptr[0].id); 
			printf("CUDA plane ID %d\n", plane.id);
			printf("CUDA nLines %d\n", plane.nElems);
				
			for (int j = 0; j < plane.nElems; j++)
			{
				line =  line_cuda(plane.ptr[j].id, plane.ptr[j].nElems, plane.ptr[j].ptr);
				cells2 = line.report_d(line.ptr[line.nElems-1].id, line.ptr[0].id);
				printf("CUDA line ID %d\n", line.id);
				printf("CUDA nCells %d\n", line.nElems);				
								
				for (int k = 0; k < line.nElems; k++)
				{
					cell = gridcell_cuda(line.ptr[k].id, line.ptr[k].nElems, line.ptr[k].ptr);
					printf("CUDA cell ID %d\n", cell.id);
					printf("CUDA nBalls %d\n", cell.nElems);
					
					for (int l = 0; l < cell.nElems; l++)
					{
						point = Point_cuda(cell.ptr[l].x, cell.ptr[l].y, cell.ptr[l].z);
//						printf("CUDA Point xx %f\n", point.x);
//						printf("CUDA Point yy %f\n", point.y);
//						printf("CUDA Point zz %f\n", point.z);
					}
//					cell.freePtr();
				}
//				line.freePtr();
			}
//			plane.freePtr();
		}	


//		plane.freePtr();
//		line.freePtr();
//		cell.freePtr();

		free(planes2);
		free(lines2);
		free(cells2);
		
}


__device__ grid_cuda copyGrid(ball_cuda *balls, int *totBalls) //OK
{

//		printf("ola2\n");
		int curCellID, curLineID, curPlaneID;

		// to store temporary lists of planes, lines, an cells
		dynArray_cuda<plane_cuda> planes = dynArray_cuda<plane_cuda>(10); 	
		dynArray_cuda<line_cuda> lines;//  = dynArray_cuda<line_cuda>(); 	
		dynArray_cuda<gridcell_cuda> cells;// = dynArray_cuda<gridcell_cuda>(); 		
		dynArray_cuda<Point> particles;// = dynArray_cuda<Point>(); 		

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

			lines = dynArray_cuda<line_cuda>(10);
			while( balls[i].parentPlane == curPlaneID )		
			{

				cells = dynArray_cuda<gridcell_cuda>(10);
				while ( balls[i].parentLine == curLineID && balls[i].parentPlane == curPlaneID )
				{

					particles = dynArray_cuda<Point>(10);		
					while ( balls[i].cell_id == curCellID && balls[i].parentLine == curLineID && balls[i].parentPlane == curPlaneID )
					{
						particles.push_end(balls[i].ball);
						i++;	
					}	// end cells		

//					printf("CUDA cell id %d\n", curCellID);
//					printf("CUDA nBalls %d\n", particles.getSize());
					
					cell = gridcell_cuda(curCellID, particles.getSize(), particles.getAll());
					cells.push_end(cell);
//					particles.freeElements();
//					cell.freePtr();
					curCellID =  balls[i].cell_id;
				} // end lines

//				printf("CUDA line id %d\n", curLineID);
//				printf("CUDA nCells %d\n", cells.getSize());

				line = line_cuda(curLineID, cells.getSize(), cells.getAll());
				lines.push_end(line);
//				cells.freeElements();
//				line.freePtr();
				curLineID =  balls[i].parentLine;
			} // end planes

//			printf("CUDA plane id %d\n", curPlaneID);
//			printf("CUDA nLines %d\n", lines.getSize());

			plane = plane_cuda(curPlaneID, lines.getSize(), lines.getAll());
			planes.push_end(plane);		
//			lines.freeElements();
//			plane.freePtr();
			curPlaneID =  balls[i].parentPlane;
		} // end all
	
		grid_cuda g = grid_cuda(planes.getSize(), planes.getAll());

		printf("CUDA nPlanes %d\n", planes.getSize());

//		planes.freeElements();
			

		return g;

}

/*


__device__ float distsq_cuda(Point a, Point b) {
		float res;
		double dx, dy, dz;

		dx = a.x - b.x;
		dy = a.y - b.y;
		dz = a.z - b.z;  
		res = dx*dx + dy*dy + dz*dz;

		return res;
}



__device__ void transformPoint_cuda( Point p, double *transMat, Point *np )
{
   np->x = transMat[  0 ] * p.x + transMat[  1 ] * p.y + transMat[  2 ] * p.z + transMat[  3 ];
   np->y = transMat[  4 ] * p.x + transMat[  5 ] * p.y + transMat[  6 ] * p.z + transMat[  7 ];
   np->z = transMat[  8 ] * p.x + transMat[  9 ] * p.y + transMat[ 10 ] * p.z + transMat[ 11 ];       
}



__device__ bool pointsWithinRange_cuda(Point *q, double delta, double TRANSLATE, double DIM, int rangeCount, grid g)
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

						if(oa->distsq_cuda(*oa, *q) <= delsq) 
							return true;
					}
				}
	
      	temp2.clear();
	
      }
    }
    temp1.clear();
  }

  return false;
}







__device__ pseudoGsol_static(int offset, double *transI, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE *staticQPointsOctree, QPOINT *staticQPoints, PG *movingPG, PARAMS_IN params)
{

	int nTot = blockDim * gridDim;

	int i = blockIdx.x;
	int j = threadIdx.x;	
	int ij = i* blockDim.x + j;

	if (ij >= 0 && ij < nTot)
	{
		if ( staticQPointsOctreeFlags[ offset + i ])
		{
			if ( ( blockIdx.x == i &&  j => staticQPointsOctree[ i ].qPtsStartID) && (blockIdx.x == i && j <= staticQPointsOctree[ i ].qPtsEndID)
			{
				
				Point p, q;
               
        p.x = staticQPoints[ j ].x;
        p.y = staticQPoints[ j ].y;
        p.z = staticQPoints[ j ].z;
               
        transformPoint_cuda( p, transI, &q );
               
        if ( movingPG->pointsWithinRange_cuda( &q, params.distanceCutoff ) )
        {
        	( *pGsol ) += staticQPoints[ j ].w;

          if ( staticQPoints[ j ].h > 0 ) 
						( *pGsolHStaticPos ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
          else 
						( *pGsolHStaticNeg ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
        }  
      }            
	
		}

	}

}


__device__ pseudoGsol_moving(int offset, double *trans, double *pGsol, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE *movingQPointsOctree, QPOINT *movingQPoints, PG *staticPG, PARAMS_IN params)
{

	int nTot = blockDim * gridDim;

	int i = blockIdx.x * blockDim.x;
	int j = threadIdx.x;	
	int ij = i*j;

	if (ij >= 0 && ij < nTot)
	{
		if ( movingQPointsOctreeFlags[ offset + i ])
		{
			if ( ( blockIdx.x == i &&  j => movingQPointsOctree[ i ].qPtsStartID) && (blockIdx.x == i && j <= movingQPointsOctree[ i ].qPtsEndID)
			{
			{
				Point p, q;
               
        p.x = movingQPoints[ j ].x;
        p.y = movingQPoints[ j ].y;
        p.z = movingQPoints[ j ].z;
               
        transformPoint_cuda( p, trans, &q );
               
        if ( staticPG->pointsWithinRange_cuda( &q, params.distanceCutoff ) )
        {
        	( *pGsol ) += movingQPoints[ j ].w;

          if ( movingQPoints[ j ].h > 0 ) 
						( *pGsolHMovingPos ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
          else 
						( *pGsolHMovingNeg ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
        }  
      }            
	
		}

	}

}



*/


__global__ void kern(ball_cuda *balls, int *totBalls)
{
		//printf("ola\n");

		grid_cuda g = copyGrid(balls, totBalls);
		

		
		printGrid(g);		

		g.freePtr();
/*
		printf("CUDA totBalls %d\n", *totBalls);



		printf("CUDA cell0 id %d\n", balls[0].cell_id);
		printf("CUDA line0 id %d\n", balls[0].parentLine);
		printf("CUDA plane0 id %d\n", balls[0].parentPlane);

		printf("CUDA ball0 x %f\n", balls[0].ball.x);
		printf("CUDA ball0 y %f\n", balls[0].ball.y);
		printf("CUDA ball0 z %f\n", balls[0].ball.z);

		printf("CUDA cell id %d\n", balls[(*totBalls)-1].cell_id);
		printf("CUDA line id %d\n", balls[(*totBalls)-1].parentLine);
		printf("CUDA plane id %d\n", balls[(*totBalls)-1].parentPlane);

		printf("CUDA ball x %f\n", balls[(*totBalls)-1].ball.x);
		printf("CUDA ball y %f\n", balls[(*totBalls)-1].ball.y);
		printf("CUDA ball z %f\n", balls[(*totBalls)-1].ball.z);

*/
/*
		//printf("ball000 %f\n", balls[9].ball.x);
		dynArray_cuda<Point> arr = dynArray_cuda<Point>();
		for (int i = 0; i < 100; i++){
			arr.push_end(balls[i].ball);
		}
		for (int i = 0; i < 100; i++){
			printf("arr %f\n", arr.get(i).x);
			printf("ball %f\n", balls[i].ball.x);
		}	
		
*/

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


		nPlanes = g1.RR.getRepMap().size();
//		printf("nPlanes %d\n", nPlanes);
			
		
		grid_init = g1.RR.getRepMap().at(0);
		grid_end = g1.RR.getRepMap().at(nPlanes-1);	
		grid1 = g1.RR.report(grid_init, grid_end);
		

		// count the total number of balls
		
		// goes through all planes
		for (int i = 0; i < nPlanes; i++)
		{
			plane1 = grid1.at(i);
//			printf(" plane ID %d\n", plane1.id);

			
			nLines = (plane1.ptr)->RR.getRepMap().size();
//			printf(" nLines %d\n", nLines);
			
			plane_init =(plane1.ptr)->RR.getRepMap().at(0); 
			plane_end = (plane1.ptr)->RR.getRepMap().at(nLines-1);
			plane2 = (plane1.ptr)->RR.report(plane_init, plane_end);
			
			// goes through all lines
			for (int j = 0; j < nLines; j++)
			{
				line1 = plane2.at(j);
//				printf("line ID %d\n", line1.id);
		
				nCells = (line1.ptr)->RR.getRepMap().size();
//				printf("nCells %d\n", nCells);
				
				line_init = (line1.ptr)->RR.getRepMap().at(0); 
				line_end = (line1.ptr)->RR.getRepMap().at(nCells-1); 
				line2 = (line1.ptr)->RR.report(line_init, line_end);
				
				
				// goes through all cells
				for (int k = 0; k < nCells; k++)
				{
					cell1 = line2.at(k);
//					printf("cell ID %d\n", cell1.id);
					

					nBalls = cell1.ptr->balls.size();
//					printf("nBalls %d\n", nBalls);


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




/*
int call_to_kern(int offset, double *transI, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE *staticQPointsOctree, QPOINT *staticQPoints, PG *movingPG, PARAMS_IN params,
double *trans, double *pGsol, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE *movingQPointsOctree, QPOINT *movingQPoints, PG *staticPG,
double TRANSLATE, double DIM, int rangeCount)*/
int call_to_kern( PG *movingPG)
{

		vector<ball_cuda> ballsVec_h = convertGrid(movingPG);

		int totBalls = ballsVec_h.size();
		int *totBalls_d;



		int totSize = totBalls * (sizeof(int)*3+sizeof(Point));

//		printf("tot size %d\n", totSize);

		ball_cuda *arrayBalls_h = (ball_cuda*) malloc(totSize);
		std::copy(ballsVec_h.begin(), ballsVec_h.end(), arrayBalls_h);

		ball_cuda *arrayBalls_d;

/*

		printf("cell0 id %d\n", arrayBalls_h[0].cell_id);
		printf("line0 id %d\n", arrayBalls_h[0].parentLine);
		printf("plane0 id %d\n", arrayBalls_h[0].parentPlane);

		printf("ball0 x %f\n", arrayBalls_h[0].ball.x);
		printf("ball0 y %f\n", arrayBalls_h[0].ball.y);
		printf("ball0 z %f\n", arrayBalls_h[0].ball.z);

		printf("cell id %d\n", arrayBalls_h[totBalls-1].cell_id);
		printf("line id %d\n", arrayBalls_h[totBalls-1].parentLine);
		printf("plane id %d\n", arrayBalls_h[totBalls-1].parentPlane);

		printf("ball x %f\n", arrayBalls_h[totBalls-1].ball.x);
		printf("ball y %f\n", arrayBalls_h[totBalls-1].ball.y);
		printf("ball z %f\n", arrayBalls_h[totBalls-1].ball.z);
		
		printf("totballs %d\n", totBalls);
*/




		cudaDeviceSetLimit (cudaLimitMallocHeapSize, 134217728);
		cudaMalloc( (void **)&arrayBalls_d, totSize );
		cudaMalloc( (void **)&totBalls_d, sizeof(int) );

		cudaMemcpy( arrayBalls_d, arrayBalls_h, totSize, cudaMemcpyHostToDevice );
		cudaMemcpy( totBalls_d, &totBalls, sizeof(int), cudaMemcpyHostToDevice );
		printf("calling...\n");
		kern<<< 1,1 >>>(arrayBalls_d, totBalls_d);
		printf("end...\n");

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error: %s\n", cudaGetErrorString(err));

		cudaFree(arrayBalls_d);
		cudaFree(totBalls_d);

		free(arrayBalls_h);

		cudaDeviceReset();



	return 0;
}




