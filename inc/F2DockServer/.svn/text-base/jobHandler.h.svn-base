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


#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include "scriptHandler.h"

using namespace std;

class JobHandler
{
	DockingJobParams *dockJob;
	ScriptHandler *sh; 

	void processDockingInput();
	void processRerankingInput();
	void processF2dGenInput();
	void processQuadGenInput();

	bool submitDockingJob(ofstream *script);
	bool submitRerankingJob(ofstream *script);
	bool submitPQRJob(ofstream *script, bool receptor);
	bool submitF2dGenJob(ofstream *script);
	bool submitQuadGenJob(ofstream *script, bool receptor);

	void saveJob();

	public:
	JobHandler(DockingJobParams *dj);
	~JobHandler();

	bool submitJob();
	
	string getFileName(int index);
};


JobHandler::JobHandler(DockingJobParams *dj)
{
	dockJob = dj;	

	sh = new ScriptHandler(dockJob);

	char pidc[20];
	sprintf(pidc, "%d", dockJob->jobId);	

	string pid(pidc);

	dockJob->receptorPDBName = sh->getDataPath() + "/" + pid + "_r.pdb";
	dockJob->ligandPDBName = sh->getDataPath() + "/" + pid + "_l.pdb";

	dockJob->receptorPQRName = sh->getDataPath() + "/" + pid + "_r.pqr";
	dockJob->ligandPQRName = sh->getDataPath() + "/" + pid + "_l.pqr";

	dockJob->receptorRAWName = sh->getDataPath() + "/" + pid + "_r.raw";
	dockJob->ligandRAWName = sh->getDataPath() + "/" + pid + "_l.raw";

	dockJob->receptorIRAWName = sh->getDataPath() + "/" + pid + "_r_imp.raw";
	dockJob->ligandIRAWName = sh->getDataPath() + "/" + pid + "_l_imp.raw";

	dockJob->receptorXYZName = sh->getDataPath() + "/" + pid + "_r.xyz";
	dockJob->ligandXYZName = sh->getDataPath() + "/" + pid + "_l.xyz";

	dockJob->receptorRAWNName = sh->getDataPath() + "/" + pid + "_r.rawn";
	dockJob->ligandRAWNName = sh->getDataPath() + "/" + pid + "_l.rawn";

	dockJob->receptorQuadName = sh->getDataPath() + "/" + pid + "_r.quad";
	dockJob->ligandQuadName = sh->getDataPath() + "/" + pid + "_l.quad";

	dockJob->receptorF2dName = sh->getDataPath() + "/" + pid + "_r.f2d";
	dockJob->ligandF2dName = sh->getDataPath() + "/" + pid + "_l.f2d";

	dockJob->rmsdFileName = sh->getDataPath() + "/" + pid + "_dock.rmsd";
	dockJob->dockingInputFileName = sh->getDataPath() + "/" + pid + "_dock.inp";
	dockJob->dockingOutputFileName = sh->getDataPath() + "/" + pid + "_dock.out";

	dockJob->rerankingInputFileName = sh->getDataPath() + "/" + pid + "_rerank.inp";
	dockJob->rerankingOutputFileName = sh->getDataPath() + "/" + pid + "_rerank.out";

	dockJob->f2dGenInputFileName = sh->getDataPath() + "/" + pid + "_f2dgen.inp";
	dockJob->quadGenInputFileName = sh->getDataPath() + "/" + pid + "_quadgen.inp";

	if(dockJob->platform == PRISM2)
	{
		dockJob->scriptFileName = "/home/docking/F2DockServer/data/Jobs/job" + pid + ".scr";
	}
	else
	{
		dockJob->scriptFileName = sh->getDataPath() + "/" + "job" + pid + ".scr";
	}

	dockJob->isDockingOutputAvailable = false;

	cout<<" Created job handler "<<endl;
}


JobHandler::~JobHandler()
{
	saveJob();
}


void JobHandler::saveJob()
{
	char pidc[20];
	sprintf(pidc, "%d", dockJob->jobId);	

	string pid(pidc);

	string jobDescriptionFileName = sh->getDataPath() + "/" + "job"+pid+ ".txt";
	ofstream job(jobDescriptionFileName.c_str());

	job<<dockJob->jobId<<endl;
	job<<dockJob->jobType<<endl;

	job<<dockJob->receptorPDBName<<endl;
	job<<dockJob->ligandPDBName<<endl;

	job<<dockJob->receptorPQRName<<endl;
	job<<dockJob->ligandPQRName<<endl;

	job<<dockJob->receptorF2dName<<endl;
	job<<dockJob->ligandF2dName<<endl;

	job<<dockJob->receptorRAWNName<<endl;
	job<<dockJob->ligandRAWNName<<endl;

	job<<dockJob->receptorQuadName<<endl;
	job<<dockJob->ligandQuadName<<endl;

	job<<dockJob->dockingOutputFileName<<endl;
	job<<dockJob->rerankingOutputFileName<<endl;

	job<<dockJob->dockingInputFileName<<endl;
	job<<dockJob->rerankingOutputFileName<<endl;

	job.close();
}


string JobHandler::getFileName(int index)
{
	if(index == RECEPTOR_PQR) return dockJob->receptorPQRName;
	else if(index == LIGAND_PQR) return dockJob->ligandPQRName;

	else if(index == RECEPTOR_F2D) return dockJob->receptorF2dName;
	else if(index == LIGAND_F2D) return dockJob->ligandF2dName;

	else if(index == RECEPTOR_RAWN) return dockJob->receptorRAWNName;
	else if(index == LIGAND_RAWN) return dockJob->ligandRAWNName;

	else if(index == RECEPTOR_QUAD) return dockJob->receptorQuadName;
	else if(index == LIGAND_QUAD) return dockJob->ligandQuadName;

	else if(index == DOCKING_OUT) return dockJob->dockingOutputFileName;
	else if(index == RERANKING_OUT) return dockJob->rerankingOutputFileName;

	return NULL;
}


bool JobHandler::submitJob()
{
	int type = dockJob->jobType;

	cout << "Script file name: " << dockJob->scriptFileName.c_str() << endl;

	ofstream jobscript(dockJob->scriptFileName.c_str());	//initializing the job script
	sh->initScript(&jobscript);

	cout<<" script initialized "<<endl; 

	if(type == DOCKING)
	{
		if(!submitDockingJob(&jobscript)) return false;
	}
	else if(type == RERANKING)
	{
		if(!submitRerankingJob(&jobscript)) return false;
	}
	else if(type == F2DGEN)
	{
		if(!submitF2dGenJob(&jobscript)) return false;
	}
	else if(type == QUADGEN)
	{
		if(!submitQuadGenJob(&jobscript,true)) return false;
		if(!submitQuadGenJob(&jobscript,false)) return false;
	}

	cout<<" script prepared "<<endl;

	sh->finalizeScript(&jobscript);			//script has been prepared. Time to execute

	cout<<" script finalized "<<endl;

//	sh->submitScript();

	return true;
}


bool JobHandler::submitDockingJob(ofstream *script)	//Check if sufficient data is available. Prepare the script for docking job. Make sure to include the generation of intermediate files
{
	bool rpdb = dockJob->isReceptorPDBAvailable;
	bool rpqr = dockJob->isReceptorPQRAvailable;
	bool rf2d = dockJob->isReceptorF2dAvailable;
	bool rraw = dockJob->isReceptorRAWNAvailable;
	bool rqua = dockJob->isReceptorQuadAvailable;

	bool lpdb = dockJob->isLigandPDBAvailable;
	bool lpqr = dockJob->isLigandPQRAvailable;
	bool lf2d = dockJob->isLigandF2dAvailable;
	bool lraw = dockJob->isLigandRAWNAvailable;
	bool lqua = dockJob->isLigandQuadAvailable;

	if (!(rpdb || rpqr && (rraw || rqua))) return false;
	if (!(lpdb || lpqr && (lraw || lqua))) return false;
	
	if(!rpqr)
	{
		submitPQRJob(script, true);
	}
	if(!lpqr)
	{
		submitPQRJob(script, false);
	}

	if(!rf2d || !lf2d)
	{
		submitF2dGenJob(script);
	}

	if(!rqua)
	{
		submitQuadGenJob(script,true);
	}
		
	if(!lqua)
	{
		submitQuadGenJob(script,false);
	}

	processDockingInput();				//docking
	sh->prepareDockingScript(script);
	dockJob->isDockingOutputAvailable = true;

	if(dockJob->performRerank)			//reranking if needed
	{
		submitRerankingJob(script);
	}

	return true;
}


bool JobHandler::submitRerankingJob(ofstream *script)
{
        bool rpdb = dockJob->isReceptorPDBAvailable;
        bool rpqr = dockJob->isReceptorPQRAvailable;
        bool rf2d = dockJob->isReceptorF2dAvailable;
        bool rraw = dockJob->isReceptorRAWNAvailable;
        bool rqua = dockJob->isReceptorQuadAvailable;

        bool lpdb = dockJob->isLigandPDBAvailable;
        bool lpqr = dockJob->isLigandPQRAvailable;
        bool lf2d = dockJob->isLigandF2dAvailable;
        bool lraw = dockJob->isLigandRAWNAvailable;
        bool lqua = dockJob->isLigandQuadAvailable;

        if (!(rpdb || rpqr && (rraw || rqua))) return false;
        if (!(lpdb || lpqr && (lraw || lqua))) return false;


	if(!dockJob->isDockingOutputAvailable) return false;

        if(!rpqr)
        {
                submitPQRJob(script, true);
        }
        if(!lpqr)
        {
                submitPQRJob(script, false);
        }

        if(!rqua)
        {
                submitQuadGenJob(script,true);
        }

        if(!lqua)
        {
                submitQuadGenJob(script,false);
        }

	processRerankingInput();
	sh->prepareRerankingScript(script);

	return true;
}


bool JobHandler::submitPQRJob(ofstream *script, bool receptor)
{
	if(receptor)
	{
		if(!dockJob->isReceptorPDBAvailable) return false;

		sh->preparePQRGenScript(script,true);
		dockJob->isReceptorPQRAvailable = true;
	}
	else
	{
		if(!dockJob->isLigandPDBAvailable) return false;
	
		sh->preparePQRGenScript(script,false);
		dockJob->isLigandPQRAvailable = true;
	}

	return true;	
}


bool JobHandler::submitF2dGenJob(ofstream *script)
{
	if(!dockJob->isReceptorPDBAvailable || !dockJob->isLigandPDBAvailable) return false;

	processF2dGenInput();

	if(!dockJob->isReceptorPQRAvailable)
	{
		sh->preparePQRGenScript(script,true);
		dockJob->isReceptorPQRAvailable = true;
	}
	if(!dockJob->isLigandPQRAvailable)
	{
		sh->preparePQRGenScript(script,false);
		dockJob->isLigandPQRAvailable = true;
	}
	
	if(!dockJob->isReceptorRAWNAvailable)
	{
		sh->prepareRAWNGenScript(script,true);
		dockJob->isReceptorRAWNAvailable = true;
	}
	
	if(!dockJob->isLigandRAWNAvailable)
	{
		sh->prepareRAWNGenScript(script,false);
		dockJob->isLigandRAWNAvailable = true;
	}
	
	sh->prepareF2dGenScript(script);
	
	dockJob->isReceptorF2dAvailable = true;
	dockJob->isLigandF2dAvailable = true;

	return true;	
}


bool JobHandler::submitQuadGenJob(ofstream *script, bool receptor)
{
	processQuadGenInput();

	if(receptor)
	{
		if(!dockJob->isReceptorPDBAvailable) return false;

		if(!dockJob->isReceptorPQRAvailable)
		{
			sh->preparePQRGenScript(script,true);
			dockJob->isReceptorPQRAvailable = true;
		}

		if(!dockJob->isReceptorRAWNAvailable)
		{
			sh->prepareRAWNGenScript(script,true);
			dockJob->isReceptorRAWNAvailable = true;
		}

		sh->prepareQuadGenScript(script,true);

		dockJob->isReceptorQuadAvailable = true;
	}
	else
	{
		if(!dockJob->isLigandPDBAvailable) return false;

		if(!dockJob->isLigandPQRAvailable)
		{
			sh->preparePQRGenScript(script,false);
			dockJob->isLigandPQRAvailable = true;
		}

		if(!dockJob->isLigandRAWNAvailable)
		{
			sh->prepareRAWNGenScript(script,false);
			dockJob->isLigandRAWNAvailable = true;
		}

		sh->prepareQuadGenScript(script,false);

		dockJob->isLigandQuadAvailable = true;
	}

	return true;	
}



void JobHandler::processDockingInput()
{
	ofstream input(dockJob->dockingInputFileName.c_str(), ios_base::app);	

	string inputstring = "staticMolecule " + dockJob->receptorF2dName + "\n";
	inputstring += "movingMolecule " + dockJob->ligandF2dName + "\n"; 
	inputstring += "outFile " + dockJob->dockingOutputFileName + "\n";
	inputstring += "staticMoleculePQR " + dockJob->receptorPQRName + "\n";
	inputstring += "movingMoleculePQR " + dockJob->ligandPQRName + "\n";
	inputstring += "staticMoleculeQUAD " + dockJob->receptorQuadName + "\n";
	inputstring += "movingMoleculeQUAD " + dockJob->ligandQuadName + "\n";

	if(dockJob->isRMSDAvailable)
	{
		inputstring += "rmsdAtoms " + dockJob->rmsdFileName + "\n";
	}

	input<<inputstring;

	input.close();				

	return;
}


void JobHandler::processRerankingInput()
{
	ofstream input(dockJob->rerankingInputFileName.c_str(), ios_base::app);	

	string inputstring = "staticMoleculePQR " + dockJob->receptorPQRName + "\n";
	inputstring += "movingMoleculePQR " + dockJob->ligandPQRName + "\n";
	inputstring += "staticMoleculeQUAD " + dockJob->receptorQuadName + "\n";
	inputstring += "movingMoleculeQUAD " + dockJob->ligandQuadName + "\n";
	inputstring += "F2DockOutputFile " + dockJob->dockingOutputFileName + "\n";
	inputstring += "rerankedOutputFile " + dockJob->rerankingOutputFileName + "\n";

	input<< inputstring;

	input.close();				

	return;
}


void JobHandler::processF2dGenInput()
{
	dockJob->f2dDim = 128;
	dockJob->quadDim = 128;

	if(dockJob->isF2dGenInputAvailable)
	{
		ifstream input(dockJob->f2dGenInputFileName.c_str());	

		if(!input) return;

		input >> dockJob->f2dDim;
		input.close();
	}
	return;
}


void JobHandler::processQuadGenInput()
{
	dockJob->f2dDim = 128;
	dockJob->quadDim = 128;

	if(dockJob->isQuadGenInputAvailable)
	{
		ifstream input(dockJob->quadGenInputFileName.c_str());	

		if(!input) return;

		input >> dockJob->quadDim;
		input.close();
	}
}
