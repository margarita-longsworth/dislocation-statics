/*
 * Slave.h
 *
 */

#ifndef SLAVE_H_
#define SLAVE_H_

#include "Domain.h"
#include "localMD.h"
#include "acceptanceCheck.h"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

class Slave{

   public:
	  Slave(int MasterID,string fileName,MPI_Comm subComm);
	  ~Slave();
      void runAsSlave();
      void writeSphereFile(string fileName,vector<Particle>& sphereParticles);
      void getConfigurationFromFile(string fileName,vector<Particle>& configParticles);

   private:
      int         MasterID_;
      string      fileName_;
	  MPI_Status  status_;
	  MPI_Comm    subComm_;
};


#endif /* SLAVE_H_ */
