#include <mpiModule.h>

namespace MpiModule
{
	int np;		// Number of processes
	int rank;	// Process rank (id)

	int getRank () {
		return rank;
	}

	int getProcs () {
		return np;
	}

	int countLines(char *filename)
	{
		FILE* fp = fopen( filename, "r" );
		int i=0;
		char s[ 2000 ];

		while ( fgets( s, 1999, fp ) != NULL )
		{
			int j = skipWhiteSpaces( s, 0 );
			if ( s[ j ] != '#' ) 
				i++;
		}

		fclose(fp);
		return i;
	}

	bool readRotations(char *rotationFile, PARAMS_IN *p) {
		// number of rotations (value 0 included)
		int numberOfRotations = p->numberOfRotations;
		int nmax = countLines(rotationFile);

		if(nmax < numberOfRotations)
			numberOfRotations = nmax;

		int *sendcount, *displs, chunk, excess;

		sendcount = new int [np];	// amount of rotations of each process (0 index not included)
		displs    = new int [np];	// begin index of the process rotations

		chunk = (numberOfRotations) / np;
		excess = (numberOfRotations) - chunk * np;

		// Calculates the sendcount and displacements arrays for the rotations
		// that each process uses
		for(int j = 0; j < np; ++j){
			if(excess){
				sendcount[j] = chunk + 1;
				--excess;
			}else{
				sendcount[j] = chunk;
			}

			sendcount[j] *= 9;

			if(j == 0)
				displs[j] = 0;
			else
				displs[j] = displs[j - 1] + sendcount[j - 1];
		}
		
		float *rotations = new float [ ( numberOfRotations + 1 ) * 9 ];

		// Root process reads the data and scatters each subset
		if (!rank) {
			float buffer[9];
			ifstream fpRot (rotationFile);

			// Opens the file
			if ( !fpRot.is_open() )
			{
				printf( "Error: Failed to open parameter file %s!\n", rotationFile);
				return false;
			}

			// Reads the rotations
			for (int i = 0; i <= numberOfRotations && !fpRot.eof(); ++i) {
				int off = i*9;

				for (int j = 0; j < 9; ++j)
					fpRot >> rotations[off + j];
					
			}
			fpRot.close();
		}

		float *rotations2 = new float [ sendcount[rank] ];

		// Scatters the rotations that each process has to execute
		MPI_Scatterv(rotations, sendcount, displs, MPI_FLOAT, rotations2, sendcount[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);

		p->rotations = rotations2;
		p->numberOfRotations = sendcount[rank] / 9;

		delete [] rotations;

		return true;
	}

	// Data marshaling to pass between processes (did not manage to do it by void* array)
	void marshalTopValues(vector<ValuePosition3D> localSols, double *arr1, int *arr2){

		for (int i = 0; i < localSols.size(); ++i) {
			int offset1 = i * 23;
			int offset2 = i * 6;

			// Copy doubles
			arr1[offset1     ] = localSols[i].m_Value;
			arr1[offset1 + 1 ] = localSols[i].m_SkinSkinRealValue;
			arr1[offset1 + 2 ] = localSols[i].m_CoreCoreRealValue;
			arr1[offset1 + 3 ] = localSols[i].m_SkinCoreRealValue;
			arr1[offset1 + 4 ] = localSols[i].m_SkinSkinImaginaryValue;
			arr1[offset1 + 5 ] = localSols[i].m_CoreCoreImaginaryValue;
			arr1[offset1 + 6 ] = localSols[i].m_SkinCoreImaginaryValue;
			arr1[offset1 + 7 ] = localSols[i].m_RealValue;
			arr1[offset1 + 8 ] = localSols[i].m_ImaginaryValue;
			arr1[offset1 + 9 ] = localSols[i].m_elecValue;
			arr1[offset1 + 10] = localSols[i].m_hbondValue;
			arr1[offset1 + 11] = localSols[i].m_hydrophobicityValue;
			arr1[offset1 + 12] = localSols[i].m_vdWPotential;
			arr1[offset1 + 13] = localSols[i].m_simpComp;
			arr1[offset1 + 14] = localSols[i].m_pGsol;
			arr1[offset1 + 15] = localSols[i].m_pGsolH;
			arr1[offset1 + 16] = localSols[i].m_delDispE;
			arr1[offset1 + 17] = localSols[i].m_origScore;
			arr1[offset1 + 18] = localSols[i].m_rerankerScore;
			arr1[offset1 + 19] = localSols[i].m_clusterPenalty;
			arr1[offset1 + 20] = localSols[i].m_Translation[0];
			arr1[offset1 + 21] = localSols[i].m_Translation[1];
			arr1[offset1 + 22] = localSols[i].m_Translation[2];

			// Copy ints
			arr2[offset2    ] = localSols[i].m_origRank;
			arr2[offset2 + 1] = localSols[i].m_rerankerRank;
			arr2[offset2 + 2] = localSols[i].m_nClashes;
			arr2[offset2 + 3] = localSols[i].m_RotationIndex;
			arr2[offset2 + 4] = localSols[i].m_FineRotationIndex;
			arr2[offset2 + 5] = localSols[i].m_ConformationIndex;
		}
	}

	// Data marshaling to pass between processes (did not manage to do it by void* array)
	void marshalTopValues(TopValues *localSols, double *arr1, int *arr2){
		int n = localSols->getCurrentNumberOfPositions();
		int i = 0;

		while( n-- ) {
			ValuePosition3D sol;
			localSols->extractMin( sol );

			int offset1 = i * 23;
			int offset2 = i * 6;

			// Copy doubles
			arr1[offset1     ] = sol.m_Value;
			arr1[offset1 + 1 ] = sol.m_SkinSkinRealValue;
			arr1[offset1 + 2 ] = sol.m_CoreCoreRealValue;
			arr1[offset1 + 3 ] = sol.m_SkinCoreRealValue;
			arr1[offset1 + 4 ] = sol.m_SkinSkinImaginaryValue;
			arr1[offset1 + 5 ] = sol.m_CoreCoreImaginaryValue;
			arr1[offset1 + 6 ] = sol.m_SkinCoreImaginaryValue;
			arr1[offset1 + 7 ] = sol.m_RealValue;
			arr1[offset1 + 8 ] = sol.m_ImaginaryValue;
			arr1[offset1 + 9 ] = sol.m_elecValue;
			arr1[offset1 + 10] = sol.m_hbondValue;
			arr1[offset1 + 11] = sol.m_hydrophobicityValue;
			arr1[offset1 + 12] = sol.m_vdWPotential;
			arr1[offset1 + 13] = sol.m_simpComp;
			arr1[offset1 + 14] = sol.m_pGsol;
			arr1[offset1 + 15] = sol.m_pGsolH;
			arr1[offset1 + 16] = sol.m_delDispE;
			arr1[offset1 + 17] = sol.m_origScore;
			arr1[offset1 + 18] = sol.m_rerankerScore;
			arr1[offset1 + 19] = sol.m_clusterPenalty;
			arr1[offset1 + 20] = sol.m_Translation[0];
			arr1[offset1 + 21] = sol.m_Translation[1];
			arr1[offset1 + 22] = sol.m_Translation[2];

			// Copy ints
			arr2[offset2    ] = sol.m_origRank;
			arr2[offset2 + 1] = sol.m_rerankerRank;
			arr2[offset2 + 2] = sol.m_nClashes;
			arr2[offset2 + 3] = sol.m_RotationIndex;
			arr2[offset2 + 4] = sol.m_FineRotationIndex;
			arr2[offset2 + 5] = sol.m_ConformationIndex;

			++i;
		}
	}

	// Data unmarshaling process
	vector<ValuePosition3D> unmarshalTopValues(double *arr1, int *arr2, int size){
		vector<ValuePosition3D> localSols (size);

		for(int i = 0; i < size; ++i){
			int offset = i * 23;

			localSols[i].m_Value 				  = arr1[offset     ];
			localSols[i].m_SkinSkinRealValue 	  = arr1[offset + 1 ];
			localSols[i].m_CoreCoreRealValue	  = arr1[offset + 2 ];
			localSols[i].m_SkinCoreRealValue	  = arr1[offset + 3 ];
			localSols[i].m_SkinSkinImaginaryValue = arr1[offset + 4 ];
			localSols[i].m_CoreCoreImaginaryValue = arr1[offset + 5 ];
			localSols[i].m_SkinCoreImaginaryValue = arr1[offset + 6 ];
			localSols[i].m_RealValue 			  = arr1[offset + 7 ];
			localSols[i].m_ImaginaryValue 		  = arr1[offset + 8 ];
			localSols[i].m_elecValue 			  = arr1[offset + 9 ];
			localSols[i].m_hbondValue 			  = arr1[offset + 10];
			localSols[i].m_hydrophobicityValue    = arr1[offset + 11];
			localSols[i].m_vdWPotential 		  = arr1[offset + 12];
			localSols[i].m_simpComp 			  = arr1[offset + 13];
			localSols[i].m_pGsol 				  = arr1[offset + 14];
			localSols[i].m_pGsolH 				  = arr1[offset + 15];
			localSols[i].m_delDispE 			  = arr1[offset + 16];
			localSols[i].m_origScore 			  = arr1[offset + 17];
			localSols[i].m_rerankerScore 		  = arr1[offset + 18];
			localSols[i].m_clusterPenalty 		  = arr1[offset + 19];
			localSols[i].m_Translation[0] 		  = arr1[offset + 20];
			localSols[i].m_Translation[1] 		  = arr1[offset + 21];
			localSols[i].m_Translation[2] 		  = arr1[offset + 22];

		}

		for(int i = 0; i < size; ++i){
			int offset = i * 6;

			localSols[i].m_origRank 		 = arr2[offset    ];
			localSols[i].m_rerankerRank 	 = arr2[offset + 1];
			localSols[i].m_nClashes 		 = arr2[offset + 2];
			localSols[i].m_RotationIndex 	 = arr2[offset + 3];
			localSols[i].m_FineRotationIndex = arr2[offset + 4];
			localSols[i].m_ConformationIndex = arr2[offset + 5];
		}
		
		return localSols;
	}

	// Merges and broadcasts the TopValues
	TopValues* mergeTopValues(TopValues *inputTopValues, int topValSize, int numFreq){
		int sum = 0, displs1[np], displs2[np], nFinal1[np], nFinal2[np], nFinal[np];
		int nTotal = inputTopValues->getCurrentNumberOfPositions();

		TopValues *outputTopValues = new TopValues (topValSize, numFreq);

		// Gather all the sizes of the TopValues from the processes
		MPI_Allgather(&nTotal, 1, MPI_INT, nFinal, 1, MPI_INT, MPI_COMM_WORLD);

		// Root calculates the indexes
		//if (!rank){
			for ( int i = 0; i < np; ++i ){
				sum += nFinal[i];

				nFinal1[i] = nFinal[i] * 23;
				nFinal2[i] = nFinal[i] * 6;

				if(i){
					displs1[i] = nFinal1[i - 1] + displs1[i - 1];
					displs2[i] = nFinal2[i - 1] + displs2[i - 1];
				}else{
					displs1[i] = 0;
					displs2[i] = 0;
				}
			}
		//}
		//MPI_Bcast(&sum, 1, MPI_INT, 0, MPI_COMM_WORLD);

		double arr1[23 * nTotal], rec1[23 * sum];
		int arr2[6 * nTotal], rec2[6 * sum];
		
		// Marshaling of the data to be transferred
		marshalTopValues(inputTopValues, arr1, arr2);

		MPI_Status st;

		// Gathers all the local solutions and broadcasts to all processes
		MPI_Allgatherv(arr1, nTotal * 23, MPI_DOUBLE, rec1, nFinal1, displs1, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(arr2, nTotal * 6, MPI_INT, rec2, nFinal2, displs2, MPI_INT, MPI_COMM_WORLD);

		vector<ValuePosition3D> inputTopValues2 = unmarshalTopValues(rec1, rec2, sum);

		// Reconstructs the TopValues structure
		for(int i = 0; i < inputTopValues2.size(); ++i)
			outputTopValues->updateTopValues( inputTopValues2[i] );

		return outputTopValues;
	}

	// Merges and scatter the TopValues
	TopValues* scatterTopValues(TopValues *inputTopValues, int topValSize, int numFreq){

		// Scattering of the data
		int chunk, excess, sendcount1[np], sendcount2[np], displs1[np], displs2[np];
		int size = inputTopValues->getCurrentNumberOfPositions();

		TopValues *outputTopValues = new TopValues (topValSize, numFreq);

		double arr1[size * 23], *rec1;
		int arr2[size * 6], *rec2;

		if(!rank){
			chunk = size / np;
			excess = size - chunk * np;

			for(int i = 0; i < np; ++i){
				if(excess){
					sendcount1[i] = chunk * 23 + 23;
					sendcount2[i] = chunk * 6 + 6;

					--excess;
				}else{
					sendcount1[i] = chunk * 23;
					sendcount2[i] = chunk * 6;
				}
				if(i == 0){
					displs1[i] = 0;
					displs2[i] = 0;
				}else{
					displs1[i] = displs1[i - 1] + sendcount1[i - 1];
					displs2[i] = displs2[i - 1] + sendcount2[i - 1];
				}
			}
		}

		marshalTopValues(inputTopValues, arr1, arr2);

		// MPI scattering
		// TODO: tentar eliminar estes bcasts
		MPI_Bcast(sendcount1, np, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(sendcount2, np, MPI_INT, 0, MPI_COMM_WORLD);

		rec1 = new double [sendcount1[rank]];
		rec2 = new int [sendcount2[rank]];

		MPI_Scatterv(arr1, sendcount1, displs1, MPI_DOUBLE, rec1, sendcount1[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(arr2, sendcount2, displs2, MPI_INT, rec2, sendcount2[rank], MPI_INT, 0, MPI_COMM_WORLD);

		vector<ValuePosition3D> inputTopValues2 = unmarshalTopValues(rec1, rec2, sendcount1[rank]/23);

		// Reconstructs the TopValues structure
		for(int i = 0; i < inputTopValues2.size(); ++i)
			outputTopValues->updateTopValues( inputTopValues2[i] );

		return outputTopValues;
	}

	// Merges the rotations from all processes and broadcasts the merged list
	void gatherRotations(PARAMS *pr){

		int total, recvcount[np], displs[np];
		float *rots;

		MPI_Allreduce(&pr->pri->numberOfRotations, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		rots = new float [total * 9];

		int chunk = total / np;
		int excess = total - chunk * np;

		// Calculates the sendcount and displacements arrays for the rotations
		// that each process uses
		for(int j = 0; j < np; ++j){
			if(excess){
				recvcount[j] = chunk + 1;
				--excess;
			}else{
				recvcount[j] = chunk;
			}

			recvcount[j] *= 9;

			if(j == 0)
				displs[j] = 0;
			else
				displs[j] = displs[j - 1] + recvcount[j - 1];
		}

		MPI_Allgatherv(pr->pri->rotations, pr->pri->numberOfRotations, MPI_FLOAT, rots, recvcount, displs, MPI_FLOAT, MPI_COMM_WORLD);

		pr->pri->rotations = rots;
		pr->pri->numberOfRotations = total;
	}

	void Init(int *argc, char*** argv){
		MPI_Init(argc, argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &np);
	}

	void Finalize (){
		MPI_Finalize();
	}

	bool isRoot () {
		if(!rank)
			return true;
		else
			return false;
	}

}