/*
 * localMD.c++
 *
 */

#include "localMD.h"
#include <iostream>

localMD::localMD(void){

}
localMD::~localMD(void){

}
interfaceIMD::interfaceIMD(void){
}
interfaceIMD::~interfaceIMD(void){

}
void interfaceIMD::loadInputFile(char* inputFile){
	// initialize IMD parameter file
	libimd_init(inputFile);
}

int interfaceIMD::runLocalMD(vector<Particle>& inConfiguration,
		                      vector<Particle>& outConfiguration, MPI_Comm subGroup){

	Particle sphereParticle;
	long nParticles = inConfiguration.size();

	libimd_atom* ptrLibimd_atom;
    ptrLibimd_atom = new libimd_atom[nParticles];

    // load sphere particle attributes
    for(auto i=0; i<nParticles; i++){

    	 sphereParticle = inConfiguration[i];

         //ptrLibimd_atom[i].num  = sphereParticle.getNumber();
         ptrLibimd_atom[i].num  = sphereParticle.getmcID();
         ptrLibimd_atom[i].sort = sphereParticle.getType();
         ptrLibimd_atom[i].mass = sphereParticle.getMass();
         ptrLibimd_atom[i].x    = sphereParticle.getPosition().x;
         ptrLibimd_atom[i].y    = sphereParticle.getPosition().y;
         ptrLibimd_atom[i].z    = sphereParticle.getPosition().z;
         ptrLibimd_atom[i].epot = sphereParticle.getEpot();
    }


    //MPI_Barrier(subGroup);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    //if (myRank ==1)
    	//cout  << "Started IMD  " << endl;
    //MPI_Barrier(subGroup);
    int statusFlag = libimd_run(ptrLibimd_atom,nParticles);
    //MPI_Barrier(subGroup);

    //if (myRank ==1) cout  << "End of IMD  " << endl;
    //MPI_Barrier(subGroup);

    // get updated trajectories back
    for(auto i=0; i<nParticles; i++){

    	//sphereParticle.setNumber(ptrLibimd_atom[i].num);/
    	sphereParticle.setmcID(ptrLibimd_atom[i].num);
    	sphereParticle.setType(ptrLibimd_atom[i].sort);
    	sphereParticle.setMass(ptrLibimd_atom[i].mass);
    	sphereParticle.setPosition(ptrLibimd_atom[i].x,ptrLibimd_atom[i].y,ptrLibimd_atom[i].z);
    	sphereParticle.setEpot(ptrLibimd_atom[i].epot);

    	outConfiguration.push_back(sphereParticle);
    }
    delete[] ptrLibimd_atom;
    return statusFlag;
    //cout << "LOCAL MD: IMD Library run successful! "<<endl;
}

int interfaceIMD::runLocalMD_finiteTemp(vector<Particle>& inConfiguration,
		    	  vector<Particle>& outConfiguration,double& energyDiff,int rToVID, int vToRID, MPI_Comm subComm){

	Particle sphereParticle;
	long nParticles = inConfiguration.size();

	libimd_atom* ptrLibimd_atom;
    ptrLibimd_atom = new libimd_atom[nParticles];

    real eDiff=0.0;

    real* eDiffPtr = &eDiff;
    //eDiff = new real[1];

    // load sphere particle attributes
    for(auto i=0; i<nParticles; i++){
    	 sphereParticle = inConfiguration[i];

         ptrLibimd_atom[i].num  = sphereParticle.getmcID();
         ptrLibimd_atom[i].sort = sphereParticle.getType();
         ptrLibimd_atom[i].mass = sphereParticle.getMass();
         ptrLibimd_atom[i].x    = sphereParticle.getPosition().x;
         ptrLibimd_atom[i].y    = sphereParticle.getPosition().y;
         ptrLibimd_atom[i].z    = sphereParticle.getPosition().z;
         ptrLibimd_atom[i].epot = sphereParticle.getEpot();
    }

    //cout<<"Inside local MD before finite temp call, eDiff :"<<*eDiffPtr<<endl;
    //int statusFlag = libimd_run_finiteTemp(ptrLibimd_atom,nParticles,
    //		eDiffPtr, rToVID, vToRID);
    int statusFlag = 0;

    energyDiff = *eDiffPtr;
    //cout<<"Inside local MD after finite temp call, eDiff :"<<*eDiffPtr<<endl;

    // get updated trajectories back after post equilibration
    for(auto i=0; i<nParticles; i++){
    	sphereParticle.setmcID(ptrLibimd_atom[i].num);
    	sphereParticle.setType(ptrLibimd_atom[i].sort);
    	sphereParticle.setMass(ptrLibimd_atom[i].mass);
    	sphereParticle.setPosition(ptrLibimd_atom[i].x,ptrLibimd_atom[i].y,ptrLibimd_atom[i].z);
    	sphereParticle.setEpot(ptrLibimd_atom[i].epot);

    	outConfiguration.push_back(sphereParticle);
    }
    delete[] ptrLibimd_atom;
    return statusFlag;
}
