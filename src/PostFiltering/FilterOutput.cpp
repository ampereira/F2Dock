#include "fast-clash/clashFilter.h"
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::vector<double> bVector;

using namespace std;

void readPQR(const char* filename, int& numMovingAtoms, double*& movingAtoms)
{

  FILE* pqr = fopen(filename,"rt");

  if(!pqr)
    {
      cout<<"Could not open pqr file"<<endl;
      assert(false);
    }


  char line[ 500 ];
  std::vector<char*> lines;
  
  while (fgets(line, 300, pqr) != NULL)
    {
      if (strncmp(line, "ATOM", 4)==0 || strncmp(line, "HETATM", 6)==0)
	{
	 
	  lines.push_back(line);

	}
    }

  double x, y, z, q, r;
  numMovingAtoms = lines.size();

  movingAtoms = new double[5 * numMovingAtoms];
  for(int i=0, j=0; i<lines.size(); i++) 
    {
      if (sscanf(lines.at(i) + 30, "%lf %lf %lf %lf %lf", &x, &y, &z, &q, &r) == 5)
    {

      movingAtoms[ j++ ] = x;
      movingAtoms[ j++ ] = y;
      movingAtoms[ j++ ] = z;
      movingAtoms[ j++ ] = q;
      movingAtoms[ j++ ] = r;

    }
      
      else {

	std::cout<<"Could not read line "<<i<<" from PQR file.\n";
	assert(false);

      }
    }
  
  fclose(pqr);

}

void readRaw(const char* forbiddenVolFileName, int& numStaticAtoms, double*& staticAtoms)
{

    std::cout<<"Forbidden volume bounded by surface "<<forbiddenVolFileName<<std::endl;

    double x,y,z;
    double charge = 0.0;
    double radius = 10.0;
    FILE* rawfp = fopen(forbiddenVolFileName, "r");

    //  int i = 0;
    if(rawfp) 
      {
	int dum;
	
	if(!(fscanf(rawfp, "%d %d", &numStaticAtoms, &dum)==2))
	  std::cerr<<"Cannot read first line\n";
      }

    else { 

      std::cerr<<"Cannot read forbidden volume raw file. Forbidden volume filter will not be applied.\n";
      assert(false);

    }

    staticAtoms = new double [5 * numStaticAtoms];

    if ( ( staticAtoms == NULL ) ) { 

      std::cerr<<"Cannot allocate memory. Forbidden volume filter will not be applied.\n";
      assert(false);

    }

    for ( int i = 0, j = 0; i < numStaticAtoms; i++ )
      {
	if((fscanf(rawfp, "%lf %lf %lf", &x, &y, &z)==3)) {
	  
	  staticAtoms[ j++ ] = x;
	  staticAtoms[ j++ ] = y;
	  staticAtoms[ j++ ] = z;
	  staticAtoms[ j++ ] = charge;
	  staticAtoms[ j++ ] = radius;

	}
      }
    fclose(rawfp);
}

void readTransformation(const char* filename, std::vector<Matrix>& transformations)
{

  std::ifstream ifile(filename);
  Matrix T(4,4);
  std::string dummy;
  ifile>>dummy;
  while(ifile>>T(0,0)>>T(0,1)>>T(0,2)>>T(0,3)>>T(1,0)>>T(1,1)>>T(1,2)>>T(1,3)>>T(2,0)>>T(2,1)>>T(2,2)>>T(2,3)>>T(3,0)>>T(3,1)>>T(3,2)>>T(3,3)) 
    transformations.push_back(T);
  ifile.close();

}

bool initForbiddenVolumeFilter(const char* forbiddenVolFileName, const char* movingAtomsFileName,
			       clashFilter **cFilter )
{

  std::cout<<"Entering initForbiddenVolumeFilter().\n";
  int numStaticAtoms = 0, numMovingAtoms = 0;//pr->numCentersB;

  double* staticAtoms; double* movingAtoms;

  readRaw(forbiddenVolFileName, numStaticAtoms, staticAtoms);

  readPQR(movingAtomsFileName, numMovingAtoms,  movingAtoms);

  ( *cFilter ) = new clashFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, false, 1.0 );

  delete [] staticAtoms;
  delete [] movingAtoms;
  std::cout<<"Leaving initForbiddenVolumeFilter().\n";

  return true;
}

void applyFilter(clashFilter* fFilter, const std::vector<Matrix>& transformations, const char* transformationFileName)
{

  FILE * fp = fopen(transformationFileName, "w");
  int count = 0;
  int count_2000 = 0;
  std::cout<<"Processing transformations\n";
  for(int i=0; i<transformations.size(); i++) {

    double intVal;
    int nclashes = 0, nclose = 0;
    Matrix m = transformations.at(i);
    CCVOpenGLMath::Matrix transM(m(0,0), m(0,1), m(0,2), m(0,3),
				 m(1,0), m(1,1), m(1,2), m(1,3),
				 m(2,0), m(2,1), m(2,2), m(2,3),
				 m(3,0), m(3,1), m(3,2), m(3,3));
    
    
    fFilter->computeInteractions( transM, &nclashes, &nclose, &intVal );
    
    if ( nclashes == 0 ) {

      fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d\n", 
	      m(0,0), m(0,1), m(0,2), m(0,3),
	      m(1,0), m(1,1), m(1,2), m(1,3),
	      m(2,0), m(2,1), m(2,2), m(2,3),
	      m(3,0), m(3,1), m(3,2), m(3,3),
	      0, 0);
      count++;
      if(i <=2000)
	count_2000++;

    }
  }
  fclose(fp);

  std::cout<<"Number of accepted transformations = "<<count<<std::endl;
  std::cout<<"Number of accepted transformations in the top 2000 = "<<count_2000<<std::endl;
  std::cout<<"Total number of transformations = "<<transformations.size()<<std::endl;

}

int main(int argc, char** argv)
{

  if (argc < 5) {

    std::cout<<"Usage:./FilterOutput <Ligand PQR> <Forbidden Volume RAW> <XForms File> <OutputFile>.\n";
    return -1;

  }

  const char* ligandPQR = argv[1];
  const char* fVolRaw = argv[2];
  const char* XFormsFile = argv[3];
  const char* output = argv[4];
  clashFilter* cFilter;
  std::vector<Matrix> transformations;
  
  initForbiddenVolumeFilter(fVolRaw, ligandPQR,&cFilter);
  readTransformation(XFormsFile,  transformations);
  applyFilter(cFilter, transformations, output);  

}
