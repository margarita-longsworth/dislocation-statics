/*
 * localMD.h
 *
 */

#ifndef MC___LOCALMD_H_
#define MC___LOCALMD_H_

#include "montecarlo_imd.h"
#include "Particle.h"
#include "Cell.h"
#include "mpi.h"

using namespace std;

// Abstract base class
class localMD{

    public:
	    virtual void loadInputFile(char* inputFile)=0;
	    virtual int runLocalMD(vector<Particle>& inConfiguration,
	    		vector<Particle>& outConfiguration, MPI_Comm subComm)=0;
	    virtual int runLocalMD_finiteTemp(vector<Particle>& inConfiguration,
	    		vector<Particle>& outConfiguration,double& eDiff,int rToVID, int vToRID, MPI_Comm subComm)=0;
	    virtual ~localMD();

    protected:
	    localMD();
};

// concrete class
class interfaceIMD:public localMD{

      public:
	      interfaceIMD();
	      virtual ~interfaceIMD();
	      virtual void loadInputFile(char* inputFile);
	      virtual int runLocalMD(vector<Particle>& inConfiguration,
	    		  vector<Particle>& outConfiguration, MPI_Comm subComm);
		  virtual int runLocalMD_finiteTemp(vector<Particle>& inConfiguration,
		    	  vector<Particle>& outConfiguration,double& eDiff,int rToVID, int vToRID, MPI_Comm subComm);
};

#endif /* MC___LOCALMD_H_ */
