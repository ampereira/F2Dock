#ifndef MPI_MODULE
#define MPI_MODULE

#include <mpi.h>
#include <Docking.h>
#include <TopValues.h>
#include <fstream>

namespace MpiModule
{
	int getRank (void);
	int getProcs (void);

	bool readRotations (char*, PARAMS_IN*);
	void marshalTopValues (vector<ValuePosition3D>, double*, int*);
	void marshalTopValues (TopValues*, double*, int*);
	vector<ValuePosition3D> unmarshalTopValues (double*, int*, int);
	TopValues* mergeTopValues (TopValues*, int, int);
	TopValues* scatterTopValues (TopValues*, int, int);
	void gatherRotations (PARAMS*);
	void Init (int*, char***);
	void Finalize (void);
	bool isRoot (void);
}

#endif