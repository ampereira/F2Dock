
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
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cstring>

#include "dockingJobParams.h"

using namespace std;

class ScriptHandler
{
	DockingJobParams *dockJob;

	string MolSurfPath;
	string F2DockPath;
	string GBRerankPath;
	string PDB2PQRPath;
	string DataPath;

	public:
	ScriptHandler(DockingJobParams *dj);
	~ScriptHandler();

	void initScript(ofstream *script);
	void finalizeScript(ofstream *script);
	void submitScript(void);

	void prepareDockingScript(ofstream *script);
	void prepareRerankingScript(ofstream *script);

	void prepareF2dGenScript(ofstream *script);
	void preparePQRGenScript(ofstream *script, bool receptor);
	void prepareRAWNGenScript(ofstream *script, bool receptor);
	void prepareQuadGenScript(ofstream *script, bool receptor);

	void setupStartScript(int id);
	void setupFinishScript(string id);

	string getDataPath(){return DataPath;}
};

ScriptHandler::ScriptHandler(DockingJobParams *dj)
{
	dockJob = dj;

	char s[ 2000 ], line[ 2000 ];  // buffer to read lines
	char *key, *val;	
	char sep[] = " ";

	FILE *fp = fopen("paths.ini","rt");

	if (  fp == NULL )
	{
		printf( "\n\nError: Failed to open paths file-- paths.ini\n\n");
	}

	while ( fgets( s, 1999, fp ) != NULL )
	{
  		if (strlen(s)<3) continue;

		strcpy( line, s );

		key = strtok(s, sep);
		val = strtok(NULL, sep);

		printf("key val = %s %s \n", key, val);

		if (val[strlen(val)-1]=='\n')  // remove unix new line
			val[strlen(val)-1] = '\0';

		if (val[strlen(val)-1]=='\r') // rmove windows new line
			val[strlen(val)-1] = '\0';

		if (strcasecmp(key, "MolSurf")==0) 
		{
			MolSurfPath = val;
		}
		else if (strcasecmp(key, "F2Dock")==0) 
		{
			F2DockPath = val;
		}
		else if (strcasecmp(key, "GBRerank")==0) 
		{
			GBRerankPath = val;
		}
		else if (strcasecmp(key, "PDB2PQR")==0) 
		{
			PDB2PQRPath = val;
		}
		else if (strcasecmp(key, "Data")==0) 
		{
			DataPath = val;
		}

	}
	fclose(fp);
}


void ScriptHandler::setupStartScript(int id)
{
	ofstream stOutFP("gen-start-files.pl");
	if(!stOutFP)
	{
		cout<<"gen-start-files.pl script could not be created."<<endl;
		return;
	}
	else
	{
		stOutFP << "#!/lusr/bin/perl" << endl;
		stOutFP << "print \"\n\nStarting script\n\n\";" << endl;
		stOutFP << "$dummyFile = qq("<<DataPath<<"/"<<id<<"_start.txt);" << endl;
		stOutFP << "open OUT, \"> $dummyFile\" or die \"\nError: Cannot create dummy file!\n\n\"; " << endl;
		stOutFP << "print OUT \"I am done.\n\"; " << endl;
		stOutFP << "close OUT; " << endl;
		stOutFP << endl;
		stOutFP.close();
	}
}

void ScriptHandler::setupFinishScript(string id)
{
	ofstream stOutFP("gen-fin-files.pl");
	if(!stOutFP)
	{
		cout<<"gen-fin-files.pl script could not be created."<<endl;
		return;
	}
	else
	{
		stOutFP << "#!/lusr/bin/perl" << endl;
		stOutFP << "print \"\n\nStarting script\n\n\";" << endl;
		stOutFP << "$dummyFile = qq("<<DataPath<<"/"<<id<<"_done.txt);" << endl;
		stOutFP << "open OUT, \"> $dummyFile\" or die \"\nError: Cannot create dummy file!\n\n\"; " << endl;
		stOutFP << "print OUT \"I am done.\n\"; " << endl;
		stOutFP << "close OUT; " << endl;
		stOutFP << endl;
		stOutFP.close();
	}

}


void ScriptHandler::initScript(ofstream *script)
{
	if(dockJob->platform == PRISM2)
	{
		(*script) << "#!/bin/bash\n";
		(*script) << "#$ -V\n";
		(*script) << "#$ -S /bin/bash\n";
		(*script) << "#$ -N job"<<dockJob->jobId<<"\n";
		(*script) << "#$ -cwd\n";
		(*script) << "#$ -e ./data/Runs/\n";
		(*script) << "#$ -o ./data/Runs/\n";
		(*script) << "export LD_LIBRARY_PATH=/home/docking/F2DockServer:$LD_LIBRARY_PATH\n";
	}
	(*script) << "echo \"Starting...\""<<"\n";
	setupStartScript(dockJob->jobId);
	(*script) << "perl ./gen-start-files.pl "<<"\n";
}


void ScriptHandler::finalizeScript(ofstream *script)
{
	char pidc[20];
	sprintf(pidc, "%d", dockJob->jobId);	

	string pid(pidc);

	(*script) << "echo \"Finished...\""<<"\n";
	setupFinishScript(pid);
	(*script) << "perl ./gen-fin-files.pl "<<pid<<"\n";


	//ofstream jobscript(dockJob->scriptFileName.c_str());	
	//jobscript << (*script);
	//jobscript.close();
}


void ScriptHandler::submitScript(void)
{
	if(dockJob->platform == LINUX)
	{
		string the_args = dockJob->scriptFileName;
		string command ("chmod 700 ");	
		command += the_args;

		cout<<"Executing "<<the_args<<endl;
		std::system(command.c_str()); 
		std::system(the_args.c_str()); 				
	}	
	else if(dockJob->platform == PRISM2)
	{
		string the_args = "qsub " + dockJob->scriptFileName;
		cout<<"Executing "<<the_args<<endl;
		std::system(the_args.c_str()); 				
		cout<<"Submission complete"<<endl;
	}
}


void ScriptHandler::prepareDockingScript(ofstream *script)
{
	(*script) << "echo \"Starting Docking...\""<<"\n";
	(*script) << F2DockPath <<" "<<dockJob->dockingInputFileName<<"\n";
	(*script) << "echo \"Finished Docking...\""<<"\n";
	
}


void ScriptHandler::prepareRerankingScript(ofstream *script)
{
	(*script) << "echo \"Starting Reranking...\""<<"\n";
	(*script) << GBRerankPath <<" "<<dockJob->rerankingInputFileName<<"\n";
	(*script) << "echo \"Finished Reranking...\""<<"\n";
}


void ScriptHandler::preparePQRGenScript(ofstream *script, bool receptor)
{
	(*script) << "echo \"Starting PQRGen...\""<<"\n";

	cout << "python "<<PDB2PQRPath<<" --ff=amber "<<dockJob->receptorPDBName<<" "<<dockJob->receptorPQRName<<"\n";

	if(receptor) 
		(*script) << "python "<<PDB2PQRPath<<" --ff=amber "<<dockJob->receptorPDBName<<" "<<dockJob->receptorPQRName<<"\n";
	else 
		(*script) << "python "<<PDB2PQRPath<<" --ff=amber "<<dockJob->ligandPDBName<<" "<<dockJob->ligandPQRName<<"\n";
		
	(*script) << "echo \"Finished PQRGen...\""<<"\n";
}


void ScriptHandler::prepareF2dGenScript(ofstream *script)
{
	(*script) << "echo \"Starting F2dGen...\""<<"\n";
	//(*script) << "../Python/python ./f2dGenerator.py "<< dockJob->ligandPQRName << " " << dockJob->receptorPQRName << " " << dockJob->f2dDim <<"\n";	
	(*script) << MolSurfPath <<" -surfaceAtoms "<<dockJob->ligandPQRName<<" "<<dockJob->ligandXYZName<<"\n";
	(*script) << MolSurfPath <<" -generateF2d "<< dockJob->ligandPQRName<<" "<<dockJob->ligandXYZName<<" "<<dockJob->ligandF2dName<<" 0 0"<<"\n";
	(*script) << MolSurfPath <<" -populateSASUsingMesh "<<dockJob->receptorPQRName<<" "<<dockJob->receptorRAWNName<<" "<<dockJob->receptorXYZName<<" 1.1 0.6 1.1 0"<<"\n";
	(*script) << MolSurfPath <<" -generateF2d "<< dockJob->receptorPQRName<<" "<<dockJob->receptorXYZName<<" "<<dockJob->receptorF2dName<<" 1 0"<<"\n";	
	(*script) << "echo \"Finished F2dGen...\""<<"\n";	
}


void ScriptHandler::prepareRAWNGenScript(ofstream *script, bool receptor)
{
	(*script) << "echo \"Starting RAWNGen...\""<<"\n";

	if(receptor)
	{
		(*script) << MolSurfPath <<" -surfaceUsingAdaptiveGrid "<<dockJob->receptorPQRName<<" "<<dockJob->receptorRAWName<<" "<< dockJob->quadDim <<"\n";
		(*script) << MolSurfPath <<" -qualityImprove "<<dockJob->receptorRAWName<<" "<<dockJob->receptorIRAWName<<"\n";
		(*script) << MolSurfPath <<" -normals -average "<<dockJob->receptorIRAWName<<" "<<dockJob->receptorRAWNName <<"\n";
	}
	else 
	{
		(*script) << MolSurfPath <<" -surfaceUsingAdaptiveGrid "<<dockJob->ligandPQRName<<" "<<dockJob->ligandRAWName<<" "<< dockJob->quadDim <<"\n";
		(*script) << MolSurfPath <<" -qualityImprove "<<dockJob->ligandRAWName<<" "<<dockJob->ligandIRAWName<<"\n";
		(*script) << MolSurfPath <<" -normals -average "<<dockJob->ligandIRAWName<<" "<<dockJob->ligandRAWNName <<"\n";
	}

	(*script) << "echo \"Finished RAWNGen...\""<<"\n";
}


void ScriptHandler::prepareQuadGenScript(ofstream *script, bool receptor)
{
	(*script) << "echo \"Starting QuadGen...\""<<"\n";

	if(receptor)
		(*script) << MolSurfPath <<" -aSpline -quad "<<dockJob->receptorRAWNName<<" gaussian 1 "<<dockJob->receptorQuadName<<"\n";
	else
		(*script) << MolSurfPath <<" -aSpline -quad "<<dockJob->ligandRAWNName<<" gaussian 1 "<<dockJob->ligandQuadName<<"\n";

	(*script) << "echo \"Finished QuadGen...\""<<"\n";
}
