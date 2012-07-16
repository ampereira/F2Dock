/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <string>

using namespace std;

enum {  LINUX, 
	PRISM2
};

enum {  RECEPTOR_PQR, 
	LIGAND_PQR,
	RECEPTOR_F2D, 
	LIGAND_F2D,
	RECEPTOR_RAWN, 
	LIGAND_RAWN,
	RECEPTOR_QUAD, 
	LIGAND_QUAD,
	DOCKING_OUT,
	RERANKING_OUT
};

enum { 	INIT,
	SUBMITTED, 
	RUNNING, 
	COMPLETED, 
	AVAILABLE
};

enum { 	DOCKING, 
	RERANKING, 
	F2DGEN, 
	QUADGEN
};


typedef struct
{
	int jobId;
	int jobType;

	bool isReceptorPDBAvailable;
	string receptorPDBName;

	bool isReceptorPQRAvailable;
	string receptorPQRName;

	bool isReceptorF2dAvailable;
	string receptorF2dName;

	bool isReceptorRAWNAvailable;
	string receptorRAWNName;

	bool isReceptorQuadAvailable;
	string receptorQuadName;

	bool isLigandPDBAvailable;
	string ligandPDBName;

	bool isLigandPQRAvailable;
	string ligandPQRName;

	bool isLigandF2dAvailable;
	string ligandF2dName;

	bool isLigandRAWNAvailable;
	string ligandRAWNName;

	bool isLigandQuadAvailable;
	string ligandQuadName;

	bool isRMSDAvailable;	
	string rmsdFileName;

	string receptorXYZName;
	string ligandXYZName;

	string receptorRAWName;
	string ligandRAWName;

	string receptorIRAWName;
	string ligandIRAWName;

	bool isDockingInputAvailable;
	string dockingInputFileName;

	bool isRerankingInputAvailable;
	string rerankingInputFileName;

	bool isF2dGenInputAvailable;
	string f2dGenInputFileName;

	bool isQuadGenInputAvailable;
	string quadGenInputFileName;

	bool isDockingOutputAvailable;
	string dockingOutputFileName;

	bool isRerankingOutputAvailable;
	string rerankingOutputFileName;

	bool performRerank;
	bool storeIntermediateFiles;

	string scriptFileName;

	bool isLocal;
	int platform;

	double f2dDim;
	double quadDim;
}DockingJobParams;
