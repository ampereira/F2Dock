#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#include "pseudoGsol_cuda.h"



// function called inside collectPseudoGsol

int pseudoGsol::call_to_kern(int offsetStatic, double *transI, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE *staticQPointsOctree, QPOINT *staticQPoints, int numStaticQPoints, int numStaticQPointsOctreeNodes, PG *movingPG, double TRANSLATEstatic, double DIMstatic, int rangeCountStatic, PARAMS_IN params, int offsetMoving, double *trans, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE *movingQPointsOctree, QPOINT *movingQPoints, int numMovingQPoints, int numMovingQPointsOctreeNodes, PG *staticPG, double TRANSLATEmoving, double DIMmoving, int rangeCountMoving)
{

// convert movingGrid to array

		vector<ball_cuda> ballsVecMoving_h = convertGrid(movingPG);

		int totBallsMoving = ballsVecMoving_h.size();
		int *totBallsMoving_d;

		int totSizeMoving = totBallsMoving * (sizeof(int)*3+sizeof(Point));

		ball_cuda *arrayMoving_h;
		cudaHostAlloc( (void**)&arrayMoving_h, totSizeMoving, cudaHostAllocDefault);\

		std::copy(ballsVecMoving_h.begin(), ballsVecMoving_h.end(), arrayMoving_h);
		ball_cuda *arrayMoving_d;	


// convert staticGrid to array

		vector<ball_cuda> ballsVecStatic_h = convertGrid(staticPG);

		int totBallsStatic = ballsVecStatic_h.size();
		int *totBallsStatic_d;

		int totSizeStatic = totBallsStatic * (sizeof(int)*3+sizeof(Point));

		ball_cuda *arrayStatic_h;
		cudaHostAlloc( (void**)&arrayStatic_h, totSizeStatic, cudaHostAllocDefault);

		std::copy(ballsVecStatic_h.begin(), ballsVecStatic_h.end(), arrayStatic_h);
		ball_cuda *arrayStatic_d;



// build grids --------------------------------------------

		cudaStream_t streamS, streamM;
		cudaStreamCreate(&streamS);
		cudaStreamCreate(&streamM);

		cudaDeviceSetLimit (cudaLimitMallocHeapSize, 512*1024*1024);

		cudaMalloc( (void **)&arrayMoving_d, totSizeMoving );
		cudaMalloc( (void **)&totBallsMoving_d, sizeof(int) );

		cudaMalloc( (void **)&arrayStatic_d, totSizeStatic );
		cudaMalloc( (void **)&totBallsStatic_d, sizeof(int) );

		cudaMemcpyAsync( arrayMoving_d, arrayMoving_h, totSizeMoving, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( totBallsMoving_d, &totBallsMoving, sizeof(int), cudaMemcpyHostToDevice, streamS );

		cudaMemcpyAsync( arrayStatic_d, arrayStatic_h, totSizeStatic, cudaMemcpyHostToDevice, streamM );
		cudaMemcpyAsync( totBallsStatic_d, &totBallsStatic, sizeof(int), cudaMemcpyHostToDevice, streamM );
	
		

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error1: %s\n", cudaGetErrorString(err));

		
		copyMovingGridtoGPU<<< 1, 1, 0, streamS >>>(arrayMoving_d, totBallsMoving_d);

 		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error2: %s\n", cudaGetErrorString(err));


		copyStaticGridtoGPU<<< 1, 1, 0, streamM >>>(arrayStatic_d, totBallsStatic_d);



		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error3: %s\n", cudaGetErrorString(err));


		cudaFree(arrayMoving_d);
		cudaFree(totBallsMoving_d);
		cudaFree(arrayStatic_d);
		cudaFree(totBallsStatic_d);

		cudaFreeHost(arrayMoving_h);
		cudaFreeHost(arrayStatic_h);


// end building grids---------------------------------------------------



// passing structs

// QPOINT

		struct QPOINT_CUDA *staticQPoints_h, *staticQPoints_d;
		cudaHostAlloc( (void**)&staticQPoints_h, numStaticQPoints * sizeof(QPOINT_CUDA), cudaHostAllocDefault);

		for (int i = 0; i < numStaticQPoints; i++)
		{
			staticQPoints_h[i].x = staticQPoints[i].x;
			staticQPoints_h[i].y = staticQPoints[i].y;
			staticQPoints_h[i].z = staticQPoints[i].z;
			staticQPoints_h[i].w = staticQPoints[i].w;
			staticQPoints_h[i].h = staticQPoints[i].h;
		}			
		

		struct QPOINT_CUDA *movingQPoints_h, *movingQPoints_d;
		cudaHostAlloc( (void**)&movingQPoints_h, numMovingQPoints * sizeof(QPOINT_CUDA), cudaHostAllocDefault);

		for (int i = 0; i < numMovingQPoints; i++)
		{
			movingQPoints_h[i].x = movingQPoints[i].x;
			movingQPoints_h[i].y = movingQPoints[i].y;
			movingQPoints_h[i].z = movingQPoints[i].z;
			movingQPoints_h[i].w = movingQPoints[i].w;
			movingQPoints_h[i].h = movingQPoints[i].h;
		}
		

// QPOINTS_OCTREE_NODE

		struct QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree_h, *staticQPointsOctree_d;
		cudaHostAlloc( (void**)&staticQPointsOctree_h, numStaticQPointsOctreeNodes * sizeof(QPOINTS_OCTREE_NODE_CUDA), cudaHostAllocDefault);


		for (int i = 0; i < numStaticQPointsOctreeNodes; i++)
		{
			staticQPointsOctree_h[i].qPtsStartID = staticQPointsOctree[i].qPtsStartID;
			staticQPointsOctree_h[i].qPtsEndID = staticQPointsOctree[i].qPtsEndID;
		}

		struct QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree_h, *movingQPointsOctree_d;
		cudaHostAlloc( (void**)&movingQPointsOctree_h, numMovingQPointsOctreeNodes * sizeof(QPOINTS_OCTREE_NODE_CUDA), cudaHostAllocDefault);


		for (int i = 0; i < numMovingQPointsOctreeNodes; i++)
		{
			movingQPointsOctree_h[i].qPtsStartID = movingQPointsOctree[i].qPtsStartID;
			movingQPointsOctree_h[i].qPtsEndID = movingQPointsOctree[i].qPtsEndID;
		}

// others

		struct static_params *statS_h, *statS_d;
		cudaHostAlloc( (void**)&statS_h, 3*(sizeof(int)+sizeof(double)), cudaHostAllocDefault);

		statS_h->offset = offsetStatic;
		statS_h->numStaticQPointsOctreeNodes = numStaticQPointsOctreeNodes;
		statS_h->distanceCutoff = params.distanceCutoff;
		statS_h->TRANSLATE = TRANSLATEstatic;
		statS_h->DIM = DIMstatic;
		statS_h->rangeCount = rangeCountStatic;


		struct moving_params *statM_h, *statM_d;
		cudaHostAlloc( (void**)&statM_h, 3*(sizeof(int)+sizeof(double)), cudaHostAllocDefault);

		statM_h->offset = offsetMoving;
		statM_h->numMovingQPointsOctreeNodes = numMovingQPointsOctreeNodes;
		statM_h->distanceCutoff = params.distanceCutoff;
		statM_h->TRANSLATE = TRANSLATEmoving;
		statM_h->DIM = DIMmoving;
		statM_h->rangeCount = rangeCountMoving;


		bool *staticQPointsOctreeFlags_d, *movingQPointsOctreeFlags_d;
		double *trans_d, *transI_d;

		double *pGsolStatic_d, *pGsolMoving_d, *pGsolStatic_h, *pGsolMoving_h;
		double *pGsolHStaticPos_d, *pGsolHStaticNeg_d, *pGsolHStaticPos_h, *pGsolHStaticNeg_h;
		double *pGsolHMovingPos_d, *pGsolHMovingNeg_d, *pGsolHMovingPos_h, *pGsolHMovingNeg_h;




// allocate arrays to return the results

		cudaHostAlloc( (void**)&pGsolStatic_h, sizeof(double) * 4*128, cudaHostAllocDefault);
		cudaHostAlloc( (void**)&pGsolMoving_h, sizeof(double) * 4*128, cudaHostAllocDefault);

		cudaHostAlloc( (void**)&pGsolHStaticPos_h, sizeof(double) * 4*128, cudaHostAllocDefault);
		cudaHostAlloc( (void**)&pGsolHStaticNeg_h, sizeof(double) * 4*128, cudaHostAllocDefault);

		cudaHostAlloc( (void**)&pGsolHMovingPos_h, sizeof(double) * 4*128, cudaHostAllocDefault);
		cudaHostAlloc( (void**)&pGsolHMovingNeg_h, sizeof(double) * 4*128, cudaHostAllocDefault);



// allocate vairables to return the results, in case the reduction is done inside the GPU

		*pGsol = *pGsolHStaticPos = *pGsolHStaticNeg = *pGsolHMovingPos = *pGsolHMovingNeg = 0;

		/*double *pGsolStatic_scalar_d, *pGsolMoving_scalar_d, *pGsolHStaticPos_scalar_d, *pGsolHStaticNeg_scalar_d, *pGsolHMovingPos_scalar_d, *pGsolHMovingNeg_scalar_d;
		
		double *pGsolStatic_scalar_h, *pGsolMoving_scalar_h;

		cudaHostAlloc( (void**)&pGsolStatic_scalar_h, sizeof(double), cudaHostAllocDefault);
		cudaHostAlloc( (void**)&pGsolMoving_scalar_h, sizeof(double), cudaHostAllocDefault);
*/

		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error1: %s\n", cudaGetErrorString(err));



// cudaMalloc ---------------------------


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


		cudaMalloc( (void **)&pGsolStatic_d, sizeof(double) * 4*128);
		cudaMalloc( (void **)&pGsolMoving_d, sizeof(double) * 4*128);

		cudaMalloc( (void **)&pGsolHStaticPos_d, sizeof(double) * 4*128);
		cudaMalloc( (void **)&pGsolHStaticNeg_d, sizeof(double) * 4*128);

		cudaMalloc( (void **)&pGsolHMovingPos_d, sizeof(double) * 4*128);
		cudaMalloc( (void **)&pGsolHMovingNeg_d, sizeof(double) * 4*128);



/*
		cudaMalloc( (void **)&pGsolStatic_scalar_d, sizeof(double));
		cudaMalloc( (void **)&pGsolMoving_scalar_d, sizeof(double));

		cudaMalloc( (void **)&pGsolHStaticPos_scalar_d, sizeof(double));
		cudaMalloc( (void **)&pGsolHStaticNeg_scalar_d, sizeof(double));

		cudaMalloc( (void **)&pGsolHMovingPos_scalar_d, sizeof(double));
		cudaMalloc( (void **)&pGsolHMovingNeg_scalar_d, sizeof(double));
*/

	
// cudaMemcpy-------------------

/*
		cudaMemcpyAsync( arrayMoving_d, arrayMoving_h, totSizeMoving, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( totBallsMoving_d, &totBallsMoving, sizeof(int), cudaMemcpyHostToDevice, streamS );


		cudaMemcpyAsync( arrayStatic_d, arrayStatic_h, totSizeStatic, cudaMemcpyHostToDevice, streamM );
		cudaMemcpyAsync( totBallsStatic_d, &totBallsStatic, sizeof(int), cudaMemcpyHostToDevice, streamM );
*/

		cudaMemcpyAsync( staticQPoints_d, staticQPoints_h, sizeof(QPOINT_CUDA)*numStaticQPoints, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( movingQPoints_d, movingQPoints_h, sizeof(QPOINT_CUDA)*numMovingQPoints, cudaMemcpyHostToDevice, streamM );


		cudaMemcpyAsync( staticQPointsOctree_d, staticQPointsOctree_h, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numStaticQPointsOctreeNodes, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( movingQPointsOctree_d, movingQPointsOctree_h, sizeof(QPOINTS_OCTREE_NODE_CUDA)*numMovingQPointsOctreeNodes, cudaMemcpyHostToDevice, streamM );



		cudaMemcpyAsync( staticQPointsOctreeFlags_d, staticQPointsOctreeFlags, sizeof(bool)*numStaticQPointsOctreeNodes, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( movingQPointsOctreeFlags_d, movingQPointsOctreeFlags, sizeof(bool)*numMovingQPointsOctreeNodes, cudaMemcpyHostToDevice, streamM );



		cudaMemcpyAsync( trans_d, trans, sizeof(double)*12, cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( transI_d, transI, sizeof(double)*12, cudaMemcpyHostToDevice, streamM );


		cudaMemcpyAsync( statS_d, statS_h, 3*(sizeof(int)+sizeof(double)), cudaMemcpyHostToDevice, streamS );
		cudaMemcpyAsync( statM_d, statM_h, 3*(sizeof(int)+sizeof(double)), cudaMemcpyHostToDevice, streamM );




// call Kernel ---------------------


		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error108: %s\n", cudaGetErrorString(err));
);

		kern1<<<4,128, 0, streamS >>>(&(statS_d->offset), transI_d, pGsolStatic_d, pGsolHStaticPos_d, pGsolHStaticNeg_d, staticQPointsOctreeFlags_d, staticQPointsOctree_d, staticQPoints_d, &(statS_d->numStaticQPointsOctreeNodes), &(statS_d->TRANSLATE), &(statS_d->DIM), &(statS_d->rangeCount), &(statS_d->distanceCutoff));//, pGsolStatic_scalar_d, pGsolHStaticPos_scalar_d, pGsolHStaticNeg_scalar_d);

		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error2: %s\n", cudaGetErrorString(err));

	kern2<<< 4,128, 0, streamM >>>(&(statS_d->distanceCutoff), &(statM_d->offset), trans_d, pGsolMoving_d, pGsolHMovingPos_d, pGsolHMovingNeg_d, movingQPointsOctreeFlags_d, movingQPointsOctree_d, movingQPoints_d, &(statM_d->numMovingQPointsOctreeNodes), &(statM_d->TRANSLATE), &(statM_d->DIM), &(statM_d->rangeCount));//, pGsolMoving_scalar_d, pGsolHMovingPos_scalar_d, pGsolHMovingNeg_scalar_d);



		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error3: %s\n", cudaGetErrorString(err));


// copy back the results, either in arrays or in simple variables

/*
	cudaMemcpyAsync( pGsolStatic_scalar_h, pGsolStatic_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamS );	
	cudaMemcpyAsync( pGsolHStaticPos, pGsolHStaticPos_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamS );	
	cudaMemcpyAsync( pGsolHStaticNeg, pGsolHStaticNeg_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamS );	

	cudaMemcpyAsync( pGsolMoving_scalar_h, pGsolMoving_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamM );	
	cudaMemcpyAsync( pGsolHMovingPos, pGsolHMovingPos_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamM );	
	cudaMemcpyAsync( pGsolHMovingNeg, pGsolHMovingNeg_scalar_d, sizeof(double), cudaMemcpyDeviceToHost, streamM );	
*/


	cudaMemcpyAsync( pGsolStatic_h, pGsolStatic_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamS );	
	cudaMemcpyAsync( pGsolHStaticPos_h, pGsolHStaticPos_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamS );	
	cudaMemcpyAsync( pGsolHStaticNeg_h, pGsolHStaticNeg_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamS );	

	cudaMemcpyAsync( pGsolMoving_h, pGsolMoving_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamM );	
	cudaMemcpyAsync( pGsolHMovingPos_h, pGsolHMovingPos_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamM );	
	cudaMemcpyAsync( pGsolHMovingNeg_h, pGsolHMovingNeg_d, sizeof(double) * 4*128, cudaMemcpyDeviceToHost, streamM );	


// synchronize

		cudaStreamSynchronize(streamS);
		cudaStreamSynchronize(streamM);

		err = cudaGetLastError();
		if (err != cudaSuccess) 
		  printf("Error3: %s\n", cudaGetErrorString(err));

		for (int i = 0; i < (4*128); i++)
		{
			//printf("pGsolStatic_h[i] %lf\n", pGsolStatic_h[i]);
			*pGsol += pGsolStatic_h[i];
			*pGsolHStaticPos += pGsolHStaticPos_h[i];
			*pGsolHStaticNeg += pGsolHStaticNeg_h[i];
		}
	
		for (int i = 0; i < (4*128); i++)
		{
			*pGsol += pGsolMoving_h[i];
			*pGsolHMovingPos += pGsolHMovingPos_h[i];
			*pGsolHMovingNeg += pGsolHMovingNeg_h[i];
		}	


	//	printf("CUDA *pGsol, *pGsolHStaticPos, *pGsolHStaticNeg, *pGsolHMovingPos, *pGsolHMovingNeg %lf %lf %lf %lf %lf\n", *pGsol, *pGsolHStaticPos, *pGsolHStaticNeg, *pGsolHMovingPos, *pGsolHMovingNeg);

	
// free memory
		
		cudaStreamDestroy(streamS);
		cudaStreamDestroy(streamM);

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

		cudaFree(pGsolStatic_d);
		cudaFree(pGsolMoving_d);

		cudaFree(pGsolHStaticPos_d);
		cudaFree(pGsolHStaticNeg_d);

		cudaFree(pGsolHMovingPos_d);
		cudaFree(pGsolHMovingNeg_d);
		
		cudaFreeHost(staticQPointsOctree_h);
		cudaFreeHost(staticQPoints_h);

		cudaFreeHost(movingQPointsOctree_h);
		cudaFreeHost(movingQPoints_h);

		cudaFreeHost(statS_h);
		cudaFreeHost(statM_h);

		cudaFreeHost(pGsolStatic_h);
		cudaFreeHost(pGsolMoving_h);

		cudaFreeHost(pGsolHStaticPos_h);
		cudaFreeHost(pGsolHStaticNeg_h);

		cudaFreeHost(pGsolHMovingPos_h);
		cudaFreeHost(pGsolHMovingNeg_h);

		cudaDeviceReset();

	return 0;
}



// prints the grids copied to GPU
__device__ void printGrid(grid_cuda g) 
{

		plane_cuda plane;
		line_cuda line;
		gridcell_cuda cell;
		Point_cuda point;

	//	printf("CUDA PGrid nPlanes %d\n", g.nElems);


		for (int i = 0; i < g.nElems; i++) // goes through all the planes
		{	
			plane = plane_cuda(g.ptr[i].id, g.ptr[i].nElems, g.ptr[i].ptr);
			printf("CUDA plane ID %d\n", plane.id);
		//	printf("CUDA nLines %d\n", plane.nElems);
				
			for (int j = 0; j < plane.nElems; j++) // goes through all the lines
			{
				line =  line_cuda(plane.ptr[j].id, plane.ptr[j].nElems, plane.ptr[j].ptr);
				printf("	CUDA line ID %d\n", line.id);
		//		printf("CUDA nCells %d\n", line.nElems);
				
				for (int k = 0; k < line.nElems; k++) // goes through all the gridcells
				{
					cell = gridcell_cuda(line.ptr[k].id, line.ptr[k].nElems, line.ptr[k].ptr);
					printf("		CUDA cell ID %d\n", cell.id);
					printf("			CUDA nBalls %d\n", cell.nElems);
					
					for (int l = 0; l < cell.nElems; l++) // goes through all the Points
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


// constructs the staticGrid on GPU, is invoked from the host (CPU)

__global__ void copyStaticGridtoGPU(ball_cuda *balls, int *totBalls) //OK
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

		staticGrid = grid_cuda(planes.getSize(), planes.getAll());
		planes.freeElements();

		//printGrid(staticGrid);

}



// constructs the movingGrid on GPU, is invoked from the host (CPU)

__global__ void copyMovingGridtoGPU(ball_cuda *balls, int *totBalls) 
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

		movingGrid = grid_cuda(planes.getSize(), planes.getAll());
		planes.freeElements();

		//printGrid(movingGrid);

}



// is invoked fromt he device (GPU) and constructs a grid

__device__ grid_cuda copyGrid(ball_cuda *balls, int *totBalls) 
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


		planes.freeElements();
		return g;

}




// implements the CPU function distsq on the GPU

__device__ float distsq_cuda(Point_cuda *a, Point_cuda *b) 
{ 
		float res;
		double dx, dy, dz;

		dx = a->x - b->x;
		dy = a->y - b->y;
		dz = a->z - b->z;  
		res = dx*dx + dy*dy + dz*dz;

		return res;
}


// implements the CPU function transformPoint in the GPU
__device__ void transformPoint_cuda( Point_cuda *p, double *transMat, Point_cuda *np ) 
{
   np->x = transMat[  0 ] * p->x + transMat[  1 ] * p->y + transMat[  2 ] * p->z + transMat[  3 ];
   np->y = transMat[  4 ] * p->x + transMat[  5 ] * p->y + transMat[  6 ] * p->z + transMat[  7 ];
   np->z = transMat[  8 ] * p->x + transMat[  9 ] * p->y + transMat[ 10 ] * p->z + transMat[ 11 ];       
}


//implements the CPU function pointsWithinRange in the GPU
__device__ bool pointsWithinRange_cuda(Point_cuda *q, double *delta, double *TRANSLATE, double *DIM, int *rangeCount, grid_cuda *g)
{

  Point_cuda p;
	p.x = q->x;
	p.y = q->y;
	p.z = q->z;

  p.x += *TRANSLATE;
  p.y += *TRANSLATE;
  p.z += *TRANSLATE;

  int l = (int)((p.z - (*delta)) / (*DIM));
  int h = (int)((p.z + (*delta)) / (*DIM));


  int size,size1,size2;
  
  dynArray_cuda<line_cuda>  temp1;
  dynArray_cuda<gridcell_cuda> temp2;

	dynArray_cuda<plane_cuda> S2;
	S2 = g->report_d(l,h);

  (*rangeCount)++;

  if(S2.empty())
    return false;

  size2 = S2.getSize();


  for(int i = 0; i < size2; i++) 
  {
    (*rangeCount)++;
   
    l = (int)((p.y - (*delta))/ (*DIM));
    h = (int)((p.y + (*delta))/ (*DIM));
    temp1 = (S2.get(i)).report_d(l,h);
    
    size1 = temp1.getSize();


    for(int j = 0; j < size1; j++) 
    {
      (*rangeCount)++;

      l = (int)((p.x - (*delta))/ (*DIM));
      h = (int)((p.x + (*delta))/ (*DIM));


      if(!(temp1.get(j).empty()))
      {
				temp2 = (temp1.get(j)).report_d(l,h);
						  	
				size = temp2.getSize();

				Point_cuda *oa = (Point_cuda*) malloc (sizeof (Point_cuda));

				double delsq = (*delta)*(*delta);
				for(int mn=0;mn<size;mn++)
				{
					int atomsInCell = temp2.get(mn).getSize();
					for(int pq=0;pq<atomsInCell;pq++)
					{
						oa->x = temp2.get(mn).ptr[pq].x;
						oa->y = temp2.get(mn).ptr[pq].y;
						oa->z = temp2.get(mn).ptr[pq].z;;
						float dist = (distsq_cuda(oa, q));
		
						
						if(distsq_cuda(oa, q) <= delsq)
							return true;
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




// implements the first part of collectPseudoGsol in GPU
__device__ void pseudoGsolStatic(int *offset, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATE, double *DIM, int *rangeCount, double *distanceCutoff)
{


	int i = blockDim.x*blockIdx.x+threadIdx.x;

	pGsolStatic[i] = 0;
	pGsolHStaticPos[i] = 0;
	pGsolHStaticNeg[i] = 0;

	if (i >= 0 && i < *numStaticQPointsOctreeNodes)
	{
		if ( staticQPointsOctreeFlags[ *offset + i ])
		{

		for( int j = staticQPointsOctree[ i ].qPtsStartID; j <= staticQPointsOctree[ i ].qPtsEndID; j++ )
			{	
				
				Point_cuda p, q;
               
        p.x = staticQPoints[ j ].x;
        p.y = staticQPoints[ j ].y;
        p.z = staticQPoints[ j ].z;

   	    transformPoint_cuda( &p, transI, &q );

    	  if ( pointsWithinRange_cuda( &q, distanceCutoff, TRANSLATE, DIM, rangeCount, &movingGrid) )
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



// implements the second part of collectPseudoGsol in GPU
__device__ void pseudoGsolMoving(int *offset, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes, double *distanceCutoff, double *TRANSLATE, double *DIM, int *rangeCount)
{

	int i = blockDim.x*blockIdx.x+threadIdx.x;


	pGsolMoving[i] = 0;
	pGsolHMovingPos[i] = 0;
	pGsolHMovingNeg[i] = 0;

	if (i >= 0 && i < *numMovingQPointsOctreeNodes)
	{
		if ( movingQPointsOctreeFlags[ *offset + i ])
		{
			for (int j = movingQPointsOctree[ i ].qPtsStartID;  j <= movingQPointsOctree[ i ].qPtsEndID; j++ )
			{
				Point_cuda p, q;
               
        p.x = movingQPoints[ j ].x;
        p.y = movingQPoints[ j ].y;
        p.z = movingQPoints[ j ].z;
               
        transformPoint_cuda( &p, trans, &q );
           
        if ( pointsWithinRange_cuda( &q, distanceCutoff, TRANSLATE, DIM, rangeCount, &staticGrid ) )
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



// both of the next 2 functions are used to do a reduction in CUDA, presently they only do it for threads in the same block, to do the reduction among ALL the threads they have to be modified

__device__ void warpReduce(volatile double *sdata, unsigned int tid) 
{
		int blockSize = 512;
		if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
		if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
		if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
		if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
		if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
		if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}

//template <unsigned int blockSize>
__device__ void reduce4(double *g_idata, double *g_odata, int n) 
{
		__shared__ double sdata[512];
		unsigned int tid = threadIdx.x;

		// each thread loads one element from global to shared mem

		if ( tid < n)
			sdata[tid] = g_idata[tid];	
		else
				sdata[tid] = 0;

		__syncthreads();

		// do reduction in shared mem
		for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
			if (tid < s ){
				sdata[tid] += sdata[tid + s];
			}
			__syncthreads();
		}


		if (tid == 0){
			*g_odata = sdata[0];
		}

}




// the only reason why both of the next 2 function were created instead of calling the pseudoGsol functions directly was so that some pre/pos-processing could be done without making the code "messy"

__global__ void kern1(/*static stuff*/int *offsetStatic, double *transI, double *pGsolStatic, double *pGsolHStaticPos, double *pGsolHStaticNeg, bool *staticQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *staticQPointsOctree, QPOINT_CUDA *staticQPoints, int *numStaticQPointsOctreeNodes,  double *TRANSLATEstatic, double *DIMstatic, int *rangeCountStatic, double *distanceCutoff)// /*others*/, double *pGsolStatic_scalar_d, double *pGsolHStaticPos_scalar_d, double *pGsolHStaticNeg_scalar_d)
{
/*
		__shared__ grid_cuda movingGrid;

		if (threadIdx.x == 0)
			movingGrid = copyGrid(ballsStatic, totBallsStatic);

			__syncthreads();
*/
			pseudoGsolStatic(offsetStatic, transI, pGsolStatic, pGsolHStaticPos, pGsolHStaticNeg, staticQPointsOctreeFlags, staticQPointsOctree, staticQPoints, numStaticQPointsOctreeNodes, TRANSLATEstatic, DIMstatic, rangeCountStatic, distanceCutoff);

/*
		reduce4(pGsolStatic, pGsolStatic_scalar_d, *numStaticQPointsOctreeNodes);
		reduce4(pGsolHStaticPos, pGsolHStaticPos_scalar_d, *numStaticQPointsOctreeNodes);
		reduce4(pGsolHStaticNeg, pGsolHStaticNeg_scalar_d, *numStaticQPointsOctreeNodes);
*/
}


__global__ void kern2(double *distanceCutoff /*moving stuff*/,int *offsetMoving, double *trans, double *pGsolMoving, double *pGsolHMovingPos, double *pGsolHMovingNeg, bool *movingQPointsOctreeFlags, QPOINTS_OCTREE_NODE_CUDA *movingQPointsOctree, QPOINT_CUDA *movingQPoints, int *numMovingQPointsOctreeNodes,  double *TRANSLATEmoving, double *DIMmoving, int *rangeCountMoving)// /*others*/, double *pGsolMoving_scalar_d, double *pGsolHMovingPos_scalar_d, double *pGsolHMovingNeg_scalar_d)
{

	/*
		__shared__	grid_cuda staticGrid;

		if (threadIdx.x == 0)
			staticGrid = copyGrid(ballsMoving, totBallsMoving);

			__syncthreads();
*/
			pseudoGsolMoving(offsetMoving, trans, pGsolMoving, pGsolHMovingPos, pGsolHMovingNeg, movingQPointsOctreeFlags, movingQPointsOctree, movingQPoints, numMovingQPointsOctreeNodes, distanceCutoff, TRANSLATEmoving, DIMmoving, rangeCountMoving);
/*
		reduce4(pGsolMoving, pGsolMoving_scalar_d, *numMovingQPointsOctreeNodes);
		reduce4(pGsolHMovingPos, pGsolHMovingPos_scalar_d, *numMovingQPointsOctreeNodes);
		reduce4(pGsolHMovingNeg, pGsolHMovingNeg_scalar_d, *numMovingQPointsOctreeNodes);
*/
}




// converts the CPU grid to a vector which contains all the points with parents plane, line, and gridcell IDs 'attached'

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

		// create vector with total number of balls and plane, line, and cell ID
		
		for (int i = 0; i < nPlanes; i++)  // goes through all planes
		{
			plane1 = grid1.at(i);
	//		printf("plane ID %d\n", plane1.id);

			nLines = (plane1.ptr)->RR.getn();
//			printf(" nLines %d\n", nLines);
			
			plane_init =(plane1.ptr)->RR.getOverallMin(); 
			plane_end = (plane1.ptr)->RR.getOverallMax();

			plane2 = (plane1.ptr)->RR.report(plane_init, plane_end);

			
			for (int j = 0; j < nLines; j++)  // goes through all lines
			{
				line1 = plane2.at(j);
	//			printf("	line ID %d\n", line1.id);

				nCells = (line1.ptr)->RR.getn();
		//		printf("nCells  %d\n", nCells);

				line_init = (line1.ptr)->RR.getOverallMin(); 
				line_end = (line1.ptr)->RR.getOverallMax(); 
				line2 = (line1.ptr)->RR.report(line_init, line_end);

				
				for (int k = 0; k < nCells; k++)  // goes through all cells
				{
					cell1 = line2.at(k);
	//				printf("		cell ID %d\n", cell1.id);
					

					nBalls = cell1.ptr->balls.size();
	//				printf("			nBalls %d\n", nBalls);
	
					for (int kk = 0 ; kk < nBalls; kk++)  //goes through all balls
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
