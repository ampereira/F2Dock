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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/


#include "XmlRPC/XmlRpc.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>

#include "F2DockServer/jobHandler.h"

using namespace XmlRpc;
using namespace std;


// The server
XmlRpcServer s;



//Used to submit jobs

class Submit : public XmlRpcServerMethod
{
public:
  Submit(XmlRpcServer* s) : XmlRpcServerMethod("Submit", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  { 
	result = -1;
	int id;

	// Retrieve and Update the process Id

	ifstream idFP("idfile");

	if(!idFP)
	{
		cout<<"id file not found. Creating a new one."<<endl;
		ofstream idOutFP("idfile");
		if(!idOutFP)
		{
			cout<<"id file could not be created."<<endl;
			return;
		}
		else
		{
			idOutFP << 0;
			idOutFP << endl;
			idOutFP.close();
			result = 0;
			id = 0;
		}
	}
	else
	{
		idFP >> id;
		idFP.close();

		ofstream idOutFP("idfile");

		if(!idOutFP)
		{
			cout<<"id file not updated"<<endl;
			return;
		}
		else
		{		
			idOutFP << (id+1);
			idOutFP << endl;
			idOutFP.close();

			result = id;	
		}
	}


	cout << "Job Id: "<< id << endl;
	
	// Retrieve and Update the process Id

	ifstream platformFP("platform.ini");
	int platform;
	platformFP >> platform;
	platformFP.close();

	//Initializing the paramemters

	DockingJobParams *dockJob = new DockingJobParams();
	dockJob->jobId = id;
	dockJob->platform = platform;
	JobHandler jh(dockJob);



	// reading from the XML

//	ofstream input(dockJob->dockingInputFileName.c_str());	
//	input << string(params[0]);
//	input.close();


	if(string(params[0]).compare("false"))		//receptor PDB
	{
		ofstream file(dockJob->receptorPDBName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->receptorPDBName.c_str()<<endl;

		else
		{
			file << string(params[1]);
			file.close();
			dockJob->isReceptorPDBAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isReceptorPDBAvailable = false;
	}

	if(string(params[2]).compare("false"))		//ligand PDB
	{
		ofstream file(dockJob->ligandPDBName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->ligandPDBName.c_str()<<endl;
		else
		{
			file << string(params[3]);
			file.close();
			dockJob->isLigandPDBAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isLigandPDBAvailable = false;
	}

	if(string(params[4]).compare("false"))		//receptor PQR
	{
		ofstream file(dockJob->receptorPQRName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->receptorPQRName.c_str()<<endl;
		else
		{
			file << string(params[5]);
			file.close();
			dockJob->isReceptorPQRAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isReceptorPQRAvailable = false;
	}

	if(string(params[6]).compare("false"))		//ligand PQR
	{
		ofstream file(dockJob->ligandPQRName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->ligandPQRName.c_str()<<endl;
		else
		{
			file << string(params[7]);
			file.close();
			dockJob->isLigandPQRAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isLigandPQRAvailable = false;
	}

	if(string(params[8]).compare("false"))		//receptor F2D
	{
		ofstream file(dockJob->receptorF2dName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->receptorF2dName.c_str()<<endl;
		else
		{
			file << string(params[9]);
			file.close();
			dockJob->isReceptorF2dAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isReceptorF2dAvailable = false;
	}

	if(string(params[10]).compare("false"))		//ligand F2D
	{
		ofstream file(dockJob->ligandF2dName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->ligandF2dName.c_str()<<endl;
		else
		{
			file << string(params[11]);
			file.close();
			dockJob->isLigandF2dAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isLigandF2dAvailable = false;
	}

	if(string(params[12]).compare("false"))		//receptor RAWN
	{
		ofstream file(dockJob->receptorRAWNName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->receptorRAWNName.c_str()<<endl;
		else
		{
			file << string(params[13]);
			file.close();
			dockJob->isReceptorRAWNAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isReceptorRAWNAvailable = false;
	}

	if(string(params[14]).compare("false"))		//ligand RAWN
	{
		ofstream file(dockJob->ligandRAWNName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->ligandRAWNName.c_str()<<endl;
		else
		{
			file << string(params[15]);
			file.close();
			dockJob->isLigandRAWNAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isLigandRAWNAvailable = false;
	}

	if(string(params[16]).compare("false"))		//receptor QUAD
	{
		ofstream file(dockJob->receptorQuadName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->receptorQuadName.c_str()<<endl;
		else
		{
			file << string(params[17]);
			file.close();
			dockJob->isReceptorQuadAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isReceptorQuadAvailable = false;
	}

	if(string(params[18]).compare("false"))		//ligand QUAD
	{
		ofstream file(dockJob->ligandQuadName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->ligandQuadName.c_str()<<endl;
		else
		{
			file << string(params[19]);
			file.close();
			dockJob->isLigandQuadAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isLigandQuadAvailable = false;
	}

	if(string(params[20]).compare("false"))		//Docking Input
	{
		ofstream file(dockJob->dockingInputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->dockingInputFileName.c_str()<<endl;
		else
		{
			file << string(params[21]);
			file.close();
			dockJob->isDockingInputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isDockingInputAvailable = false;
	}

	if(string(params[22]).compare("false"))		//Reranking Input
	{
		ofstream file(dockJob->rerankingInputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->rerankingInputFileName.c_str()<<endl;
		else
		{
			file << string(params[23]);
			file.close();
			dockJob->isRerankingInputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isRerankingInputAvailable = false;
	}

	if(string(params[24]).compare("false"))		//F2dGen Input
	{
		ofstream file(dockJob->f2dGenInputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->f2dGenInputFileName.c_str()<<endl;
		else
		{
			file << string(params[25]);
			file.close();
			dockJob->isF2dGenInputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isF2dGenInputAvailable = false;
	}

	if(string(params[26]).compare("false"))		//QuadGen Input
	{
		ofstream file(dockJob->quadGenInputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->quadGenInputFileName.c_str()<<endl;
		else
		{
			file << string(params[27]);
			file.close();
			dockJob->isQuadGenInputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isQuadGenInputAvailable = false;
	}


	if(string(params[28]).compare("false"))	//RMSD
	{
		ofstream file(dockJob->rmsdFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->rmsdFileName.c_str()<<endl;
		else
		{
			file << string(params[29]);
			file.close();
			dockJob->isRMSDAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isRMSDAvailable = false;
	}

	if(string(params[30]).compare("false"))		//Docking Output
	{
		ofstream file(dockJob->dockingOutputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->dockingOutputFileName.c_str()<<endl;
		else
		{
			file << string(params[31]);
			file.close();
			dockJob->isDockingOutputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isDockingOutputAvailable = false;
	}

	if(string(params[32]).compare("false"))		//Reranking Output
	{
		ofstream file(dockJob->rerankingOutputFileName.c_str());

		if(!file)
			cout<<"Failed to create file: "<<dockJob->rerankingOutputFileName.c_str()<<endl;
		else
		{
			file << string(params[33]);
			file.close();
			dockJob->isRerankingOutputAvailable = true;
			cout<<"read"<<endl;
		}
	}
	else
	{
		cout<<"not read"<<endl;
		dockJob->isRerankingOutputAvailable = false;
	}

	if(string(params[34]).compare("false"))
		dockJob->performRerank = true;
	else 
		dockJob->performRerank = false;	

	if(string(params[35]).compare("false"))
		dockJob->storeIntermediateFiles = true;
	else 
		dockJob->storeIntermediateFiles = false;	

	dockJob->jobType = int(params[36]);

	cout<<" ***** Read XML values **** "<<endl;

	// handling job
	jh.submitJob();

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

} submit(&s);    // This constructor registers the method with the server




//Used to retrieve the result after the job is complete

class GetOutput : public XmlRpcServerMethod
{
public:
  GetOutput(XmlRpcServer* s) : XmlRpcServerMethod("GetOutput", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {

	cout<<"************** Recieved get output request"<<endl;
	//check if output is ready

	int jobId = int(params[0]);

	char id[40];
	sprintf(id, "%d", jobId);
	string ID(id);

	string doneFileName("./data/Jobs/");
	doneFileName += ID;
	doneFileName += "_done.txt";

	cout<< "checking if"<< doneFileName.c_str() << " is present"<<endl;

	ifstream doneFile(doneFileName.c_str());

	if(!doneFile)
	{
		cout<<"************** Output not ready"<<endl;
		result = "INCOMPLETE";
		return;
	}
	else
	{
		cout<<"************** Output ready"<<endl;

		DockingJobParams *dockJob = new DockingJobParams();

		dockJob->jobId = jobId;

		JobHandler jh(dockJob);

		int requestFor = int(params[1]);
				
		string file_to_open = jh.getFileName(requestFor);

		ifstream output(file_to_open.c_str());
		string line;
		string resultBuffer;

		if (!output) 
		{
			cerr << "Unable to open output File "<<file_to_open.c_str();
			result = "INCOMPLETE";
			return;
		}

		cout << "Opened file " << file_to_open.c_str() << endl;

		while (getline(output, line)) 
		{
			resultBuffer += line + "\n";
		}
		output.close();
		doneFile.close();
	
		result = resultBuffer;

		cout<<"************** Returning Output"<<endl;
		return;
	}

  }


} getOutput(&s);



//Used to poll the server to find out if the job is complete

class getJobStatus : public XmlRpcServerMethod
{
public:
  getJobStatus(XmlRpcServer* s) : XmlRpcServerMethod("getJobStatus", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
	result = SUBMITTED;

	string startFileName = "./data/Jobs/" + string(params[0]) + "_start.txt";
	ifstream startFile(startFileName.c_str());
	if(!startFile)
	{
		return;
	}
	else
	{
		result = RUNNING;
		startFile.close();
	}

	string doneFileName = "./data/Jobs/" + string(params[0]) + "_done.txt";
	ifstream doneFile(doneFileName.c_str());
	if(!doneFile)
	{
		return;
	}
	else
	{
		result = COMPLETED;
		doneFile.close();
		return;
	}
  }

} getJobStatus(&s);



int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: F2DockServer serverPort platform(0 = standalone linux, 1 = prism2, ...)\n";
    return -1;
  }
  int port = atoi(argv[1]);
  int platform = atoi(argv[2]);

  ofstream settings("platform.ini");
  settings << platform << endl;
  settings.close();

  XmlRpc::setVerbosity(5);

  // Create the server socket on the specified port
  s.bindAndListen(port);

  // Enable introspection
  s.enableIntrospection(true);

  // Wait for requests indefinitely
  s.work(-1.0);

  return 0;
}

