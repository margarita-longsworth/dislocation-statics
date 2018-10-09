
/*
 * MasterSlavePostTest.c++
 *
 *
 * Source file implementation of Master-Slave approach for Parallel MonteCarlo(MC) Sampling.
 * Configuration and MonteCarlo parameter files are read. Master process manages the Data
 * structure allocation, random selection, sphere construction and job allocation procedures. Constructed sphere domains
 * are distributed across slave processes, which invoke IMD Library interface for equilibration/relaxation
 * of sphere domains.Each slave process test the acceptance criterion and return its decision to Master process
 * upon completion. Upon acceptance part of the domain involved in sphere gets updated by Master process.
 * A post screening is done to avoid overlapping sphere domains before getting updated to Master domain.
 */

#include "mpi.h"
#include "constructSampleZone.h"
#include "localMD.h"
#include "acceptanceCheck.h"
#include "TrialMove.h"
#include "MonteCarloGlobals.h"
#include "Master.h"
#include "Slave.h"
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <math.h>
#include <time.h>

using namespace std;


int main(int argc,char* argv[]){

    string paramFile      = "paramFile.param";      // MonteCarlo Parameter file
    string sphereParam    = "localMD.param";

	MPI_Init(&argc,&argv);
	int nProcess,pRank;

	MPI_Group MainGroup, slaveGroup;
    MPI_Comm  subComm;

	int deleteRank = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    MPI_Comm_group(MPI_COMM_WORLD, &MainGroup);
    // Sub group
    MPI_Group_excl(MainGroup,1,&deleteRank,&slaveGroup);

    //MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)

    MPI_Comm_create(MPI_COMM_WORLD, slaveGroup, &subComm);

    // Create Slave Process instance
    Slave slave(0,sphereParam,subComm);

    // ----------------------  Master Process ---------------------- //
    if(pRank == 0){ // Rank-0 is designated as Master Process

    	//readMCParamFile(paramFile);
    	//exit(1);

    	// ----------------------------------- Paramfile reading part --------------------
    	// input Parts (to be modified with suitable Param file reader)
    	// (1 1 1) for Master Slave approach

        ifstream paramin(paramFile,std::ios_base::in);
        string keyword;

        // flags for parameters
        int dimFlag=0,stepFlag=0,seedFlag=0,tTypesFlag=0,rTypesFlag=0,boxX_Flag=0,boxY_Flag=0,boxZ_Flag=0;
        int zoneFlag=0, pbcFlag=0, tempFlag=0, radFlag=0,thickFlag=0,cohesFlag=0,flushFlag=0,masterFlag=0,confFlag=0,paramFlag=0;
        int outFlag=0, statFlag=0, trialFlag=0, simIDFlag=0;

        int nConfFlag=0,nFinFlag=0,delEpotFlag=0,biasedPerformed=0,localPerformed=0; // flags for statistical parameters
        int nSwapFlag=0,nAccSwap=0,nRejSwap=0;     // for SWAP type trial move
        int nAddFlag=0,nDelFlag=0,nAddAcc=0,nDelAcc=0,nAddRej=0,nDelRej=0; // for SINGLE type trial move
	int samplingMode=0,stepsInLocalBlock=0,stepsInBiasedBlock=0,sphereRadiusTarget=0,localMovesRadius=0, coveringTimes=0,cylinderRadius=0; //flags for BIAS sampling

        while (!paramin.eof()) {

            if (paramin.get() != '#'){

                 paramin >> keyword;

                 if(keyword == "mc_Cpudim"){
                         paramin >> mccpuDimension.x;
                         paramin >> mccpuDimension.y;
                         paramin >> mccpuDimension.z;
                         dimFlag = 1;
                 }
                 else if(keyword == "mc_Steps"){
                        paramin >>  mcnSteps;     stepFlag = 1;
                 }
                 else if(keyword == "mc_Seed"){// seed for random generators used in MC routine
                        paramin >>  mcSeed;       seedFlag = 1;
                 }
                 else if(keyword == "mc_totalTypes"){
                        paramin >>  mcTotalTypes; tTypesFlag = 1;
                 }
                 else if(keyword == "mc_realTypes"){
                        paramin >>  mcRealTypes;  rTypesFlag = 1;
                 }
                 // Box dimensions
                 else if(keyword == "mc_SimboxX"){
                        paramin >>  mcSimboxX.x;
                        paramin >>  mcSimboxX.y;
                        paramin >>  mcSimboxX.z; boxX_Flag = 1;
                 }

                 else if(keyword == "mc_SimboxY"){
                        paramin >>  mcSimboxY.x;
                        paramin >>  mcSimboxY.y;
                        paramin >>  mcSimboxY.z; boxY_Flag = 1;
                 }

                 else if(keyword == "mc_SimboxZ"){
                        paramin >>  mcSimboxZ.x;
                        paramin >>  mcSimboxZ.y;
                        paramin >>  mcSimboxZ.z; boxZ_Flag = 1;
                 }
                 // Zone dimension (not utilized for MC/localMD approach)
                 else if(keyword == "mc_zoneSize"){
                        paramin >>  mcZoneSize.x;
                        paramin >>  mcZoneSize.y;
                        paramin >>  mcZoneSize.z; zoneFlag = 1;
                 }
                 // Periodic boundary conditions
                 else if(keyword == "mc_pbc"){
                         paramin >> mcPBC.x;
                         paramin >> mcPBC.y;
                         paramin >> mcPBC.z; pbcFlag = 1;
                 }
                 // simulation temperature
                 else if(keyword == "mc_temperature"){
                         paramin >> mcTemperature; tempFlag = 1;
                 }
                 // radius of spherical zone1
                 else if(keyword == "mc_sphereRadius"){
                         paramin >> mcSphereRadius; radFlag = 1;
                 }
                 else if(keyword == "mc_cylinderRadius"){
                         paramin >> mcCylinderRadius; cylinderRadius = 1;
                 }

                 // wall thickness of zone2 and zone3
                 else if(keyword == "mc_sphereWall"){
                         paramin >> mcSphereWallThickness; thickFlag = 1;
                 }
                 // cohesive energy of a single solute
                 // valid for SWAP trial move
                 // otherwise local solute concentration has to be considered for
                 // cohesive energy computation
                 else if(keyword == "mc_cohesiveEnergy"){
                         paramin >> mcCohesive; cohesFlag = 1;
                 }
                 // interval for writing output configuration
                 else if(keyword == "mc_flushInterval"){
                         paramin >> mcFlushInterval; flushFlag = 1;
                 }
                 // flag for Master Slave approach
                 else if(keyword == "mc_MasterFlag"){
               	      paramin >> mcMasterFlag; masterFlag = 1;
                 }
                 // input file configurations for MC simulation
                 else if(keyword == "mc_inputFile"){
                         paramin >> mcInputFile; confFlag = 1;
                 }
                 // final output file configurations after MC simulation
                 else if(keyword == "mc_outputFile"){
               	      paramin >> mcOutputFile; outFlag = 1;
                 }
                 // File containing the evolution of simulation
                 else if(keyword == "mc_statFile"){
               	      paramin >> mcStatFile; statFlag = 1;
                 }
                 // IMD parameter file for localized energy computation
                 else if(keyword == "mc_sphereParam"){
                         paramin >> mcSphereFile; paramFlag = 1;
                 }
                 // 1- RANDOM SINGLE TRIAL MOVE, 2- SWAP TRIAL MOVE 3 CLUSTER TRIAL MOVE
                 else if(keyword == "mc_trialMoveType"){
               	      paramin >> mcTrialMoveType; trialFlag = 1;
               	      if(mcTrialMoveType==3){
               	    	  cerr<<"ERROR: CLUSTER type trial move not implemented."<<endl;
               	    	  exit(1);
               	      }
                 }
                 // simulation ID : 1 - start New ,Otherwise start from intermediate configuration
                 else if(keyword == "mc_simulationID"){
               	      paramin >> mcSimulationID;
               	      simIDFlag = 1;
                 }
                 // statistical parameters for simulation starting from intermediate configurations
                 else if(keyword == "mc_nConflicts"){
                	 paramin >> mcNConflicts; nConfFlag=1;
                 }
                 else if(keyword == "mc_nFinished"){
                	 paramin >> mcNFinished; nFinFlag=1;
                 }
                 else if(keyword == "mc_nSwaps"){
                	 paramin >> mcNSwaps; nSwapFlag=1;
                 }
                 else if(keyword == "mc_nAccepSwaps"){
                	 paramin >> mcNAccepSwaps; nAccSwap=1;
                 }
                 else if(keyword == "mc_nRejeSwaps"){
                	 paramin >> mcNRejeSwaps; nRejSwap=1;
                 }
                 else if(keyword == "mc_totDeltaEpot"){
                	 paramin >> mcTotDeltaEpot; delEpotFlag=1;
                 }
                 else if(keyword == "mc_biasedMovesPerformed"){
                	 paramin >> mcBiasedMovesPerformed; biasedPerformed=1;
                 }
                 else if(keyword == "mc_localMovesPerformed"){
                	 paramin >> mcLocalMovesPerformed; localPerformed=1;
                 }
                 else if(keyword == "mc_nAddition"){
                	 paramin >> mcNAddition; nAddFlag=1;
                 }
                 else if(keyword == "mc_nDeletion"){
                	 paramin >> mcNDeletion; nDelFlag=1;
                 }
                 else if(keyword == "mc_nAddAccepted"){
                	 paramin >> mcNaddAccepted; nAddAcc=1;
                 }
                 else if(keyword == "mc_nDelAccepted"){
                	 paramin >> mcNdelAccepted; nDelAcc=1;
                 }
                 else if(keyword == "mc_nAddRejected"){
                	 paramin >> mcNaddRejected; nAddRej=1;
                 }
                 else if(keyword == "mc_nDelRejected"){
                	 paramin >> mcNdelRejected; nDelRej=1;
                 }
		else if(keyword == "mc_samplingMode"){				
                        paramin >> mcSamplingMode; samplingMode=1;
		}
		else if(keyword == "mc_stepsInLocalBlock"){				
                        paramin >> mcStepsInLocalBlock; stepsInLocalBlock=1;
		}
		else if(keyword == "mc_stepsInBiasedBlock"){				
                        paramin >> mcStepsInBiasedBlock; stepsInBiasedBlock=1;
		}
		else if(keyword == "mc_sphereRadiusTarget"){				
                        paramin >> mcSphereRadiusTarget; sphereRadiusTarget=1;
		}
		else if(keyword == "mc_localMovesRadius"){				
                        paramin >> mcLocalMovesRadius; localMovesRadius=1;
		}
		else if(keyword == "mc_coveringTimes"){				
                        paramin >> mcCoveringTimes; coveringTimes=1;
		}

               }
        }

	    // Starting from new configuration
	    if(mcSimulationID == 1){
	    	mcNConflicts=0;mcNFinished=0;mcNSwaps=0;mcNAccepSwaps=0;mcNRejeSwaps=0;mcBiasedMovesPerformed=0;mcLocalMovesPerformed=0;
	    	mcNAddition=0;mcNDeletion=0;mcNaddAccepted=0;mcNdelAccepted=0;mcNaddRejected=0;mcNdelRejected=0;
	    	mcTotDeltaEpot = 0.0;
	    	// Updating Flags
	        nConfFlag=1;nFinFlag=1;delEpotFlag=1;
	        nSwapFlag=1;nAccSwap=1;nRejSwap=1;biasedPerformed=1;localPerformed=1;
	        nAddFlag=1; nDelFlag=1; nAddAcc=1; nDelAcc=1; nAddRej=1; nDelRej=1; // for SINGLE type trial move
	    }
	    else if(mcSimulationID != 1){ // starting from intermediate configuration

	    	// based on trial Moves
	    	// SINGLE INSERTION/DELETION type trial move
	    	if(mcTrialMoveType==1){
		        nSwapFlag=1;nAccSwap=1;nRejSwap=1; // activate SWAP type trial flags as these parameters were not read
	    	}
	    	// SWAP type trial move
	    	else if(mcTrialMoveType==2){
	    		// activate SINGLE type trial flags
	    		nAddFlag=1; nDelFlag=1; nAddAcc=1; nDelAcc=1; nAddRej=1; nDelRej=1;
	    	}
	    }

        // Throw error in case of missing parameter
        if(dimFlag==0||stepFlag==0||seedFlag==0||tTypesFlag==0||rTypesFlag==0||boxX_Flag==0||boxY_Flag==0||boxZ_Flag==0||
           zoneFlag==0||pbcFlag==0||tempFlag==0||radFlag==0||thickFlag==0||cohesFlag==0||flushFlag==0||masterFlag==0||
		   confFlag==0||paramFlag==0||outFlag==0||statFlag==0||trialFlag==0||simIDFlag==0||
		   nConfFlag==0||nFinFlag==0||delEpotFlag==0|| nSwapFlag==0||nAccSwap==0||nRejSwap==0||
		   nAddFlag==0 ||nDelFlag==0|| nAddAcc==0|| nDelAcc==0|| nAddRej==0|| nDelRej==0 || biasedPerformed==0 || localPerformed==0|| samplingMode==0 || stepsInBiasedBlock==0 || stepsInLocalBlock==0 || sphereRadiusTarget==0 || localMovesRadius==0 || coveringTimes == 0 || cylinderRadius==0){
           if(dimFlag==0){
        	   cerr<<"ERROR: mc_Cpudim not defined "<<endl;
           }
           else if(stepFlag==0){
        	   cerr<<"ERROR: mc_Steps not defined "<<endl;
           }
           else if(seedFlag==0){
        	   cerr<<"ERROR: mc_Seed not defined "<<endl;
           }
           else if(tTypesFlag==0){
        	   cerr<<"ERROR: mc_totalTypes not defined "<<endl;
           }
           else if(rTypesFlag==0){
        	   cerr<<"ERROR: mc_realTypes not defined "<<endl;
           }
           else if(boxX_Flag==0){
        	   cerr<<"ERROR: mc_SimboxX not defined "<<endl;
           }
           else if(boxY_Flag==0){
        	   cerr<<"ERROR: mc_SimboxY not defined "<<endl;
           }
           else if(boxZ_Flag==0){
        	   cerr<<"ERROR: mc_SimboxZ not defined "<<endl;
           }
           else if(zoneFlag==0){
        	   cerr<<"ERROR: mc_zoneSize not defined "<<endl;
           }
           else if(pbcFlag==0){
        	   cerr<<"ERROR: mc_pbc not defined "<<endl;
           }
           else if(tempFlag==0){
        	   cerr<<"ERROR: mc_temperature not defined "<<endl;
           }
           else if(radFlag==0){
        	   cerr<<"ERROR: mc_sphereRadius not defined "<<endl;
           }
           else if(cylinderRadius==0){
        	   cerr<<"ERROR: mc_cylinderRadius not defined "<<endl;
           }
           else if(thickFlag==0){
        	   cerr<<"ERROR: mc_sphereWall not defined "<<endl;
           }
           else if(cohesFlag==0){
        	   cerr<<"ERROR: mc_cohesiveEnergy not defined "<<endl;
           }
           else if(flushFlag==0){
        	   cerr<<"ERROR: mc_flushInterval not defined "<<endl;
           }
           else if(masterFlag==0){
        	   cerr<<"ERROR: mc_MasterFlag not defined "<<endl;
           }
           else if(confFlag==0){
        	   cerr<<"ERROR: mc_inputFile not defined "<<endl;
           }
           else if(outFlag==0){
        	   cerr<<"ERROR: mc_outputFile not defined "<<endl;
           }
           else if(trialFlag==0){
        	   cerr<<"ERROR: mc_trialMoveType not defined "<<endl;
           }
           else if(simIDFlag==0){
        	   cerr<<"ERROR: mc_simulationID not defined "<<endl;
           }
           else if(paramFlag==0){
        	   cerr<<"ERROR: mc_sphereParam not defined "<<endl;
           }
           else if(nConfFlag==0){
        	   cerr<<"ERROR: mc_nConflicts not defined "<<endl;
           }
           else if(nFinFlag==0){
        	   cerr<<"ERROR: mc_nFinished not defined "<<endl;
           }
           else if(delEpotFlag==0){
        	   cerr<<"ERROR: mc_totDeltaEpot not defined "<<endl;
           }
           else if(nSwapFlag==0){
        	   cerr<<"ERROR: mc_nSwaps not defined "<<endl;
           }
           else if(nAccSwap==0){
        	   cerr<<"ERROR: mc_nAccepSwaps not defined "<<endl;
           }
           else if(nRejSwap==0){
        	   cerr<<"ERROR: mc_nRejeSwaps not defined "<<endl;
           }
           else if(biasedPerformed==0){
        	   cerr<<"ERROR: mc_biasedMovesPerformed not defined "<<endl;
           }
           else if(localPerformed==0){
        	   cerr<<"ERROR: mc_localMovesPerformed not defined "<<endl;
           }
           else if(nAddFlag==0){
        	   cerr<<"ERROR: mc_nAddition not defined "<<endl;
           }
           else if(nDelFlag==0){
        	   cerr<<"ERROR: mc_nDeletion not defined "<<endl;
           }
           else if(nAddAcc==0){
        	   cerr<<"ERROR: mc_nAddAccepted not defined "<<endl;
           }
           else if(nDelAcc==0){
        	   cerr<<"ERROR: mc_nDelAccepted not defined "<<endl;
           }
           else if(nAddRej==0){
        	   cerr<<"ERROR: mc_nAddRejected not defined "<<endl;
           }
           else if(nDelRej==0){
        	   cerr<<"ERROR: mc_nDelRejected not defined "<<endl;
           }
           else if(samplingMode==0){
        	   cerr<<"ERROR: mc_samplingMode not defined "<<endl;	
           }
           else if(stepsInBiasedBlock==0){
        	   cerr<<"ERROR: mc_stepsInBiasedBlock not defined "<<endl;	
           } 
           else if(stepsInLocalBlock==0){
        	   cerr<<"ERROR: mc_stepsInLocalBlock not defined "<<endl;	
           }
           else if(sphereRadiusTarget==0){
        	   cerr<<"ERROR: mc_sphereRadiusTarget not defined "<<endl;	
           }
           else if(localMovesRadius==0){
        	   cerr<<"ERROR: mc_localMovesRadius not defined "<<endl;	
           }
           else if(coveringTimes==0){
        	   cerr<<"ERROR: mc_coveringTimes not defined "<<endl;	
           }
           exit(1);
        }

        cout << " Monte Carlo Parameter File reading Successful" << endl;

        cout << " input total Steps   : "  << mcnSteps << endl;
        cout << " input mcSeed        : "  << mcSeed << endl;
        cout << " input CPU Dimension : [ "<< mccpuDimension.x <<" "<<mccpuDimension.y<<" "<<mccpuDimension.z<<" ]"<< endl;
        cout << " Simbox X            : [ "<< mcSimboxX.x <<" "<< mcSimboxX.y<<" "<<mcSimboxX.z<<" ]"<< endl;
        cout << " Simbox Y            : [ "<< mcSimboxY.x <<" "<< mcSimboxY.y<<" "<<mcSimboxY.z<<" ]"<< endl;
        cout << " Simbox Z            : [ "<< mcSimboxZ.x <<" "<< mcSimboxZ.y<<" "<<mcSimboxZ.z<<" ]"<< endl;
        cout << " Zone Size           : [ "<< mcZoneSize.x <<" "<< mcZoneSize.y<<" "<<mcZoneSize.z<<" ]"<< endl;
        cout << " Periodic Boundary   : [ "<< mcPBC.x <<" "<< mcPBC.y<<" "<<mcPBC.z<<" ]"<< endl;
        cout << " input temperature   : "  << mcTemperature << endl;
        cout << " Sphere radius       : "  << mcSphereRadius << endl;
        cout << " Sphere Wall         : "  << mcSphereWallThickness << endl;
        cout << " Cohesive Energy     : "  << mcCohesive << endl;
        cout << " Flush Interval      : "  << mcFlushInterval << endl;
        cout << " Master Flag         : "  << mcMasterFlag    << endl;
        cout << " mcInput file        : "  << mcInputFile     << endl;
        cout << " mcOutput File       : "  << mcOutputFile    << endl;
        cout << " mcStatistics file   : "  << mcStatFile      << endl;
        cout << " sphere Param file   : "  << mcSphereFile    << endl;
        cout << " trial Move type     : "  << mcTrialMoveType << endl;
        cout << " mc Simulation ID    : "  << mcSimulationID  <<endl;
        cout << " Sampling mode       : "  << mcSamplingMode << endl;
        cout << " No of bias moves    : "  << mcStepsInBiasedBlock << endl;
        cout << " No of local moves   : "  << mcStepsInLocalBlock << endl;
        cout << " sphereRadiusTarget  : "  << mcSphereRadiusTarget  <<endl;
        cout << " localMovesRadius    : "  << mcLocalMovesRadius  <<endl;
        cout << " coveringTimes	      : "  << mcCoveringTimes  <<endl;
        cout << " Cylinder Radius     : "  << mcCylinderRadius << endl;
        cout << "=========================================== "<<endl;

        //---------------- End of param file reading part --------------

        //==============================================================
        // Perform dimension check (Sphere Construction)
        if(mcSimboxX.x <= (2*(mcSphereRadius+2*mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: X dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }
        if(mcSimboxY.y <= (2*(mcSphereRadius+2*mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: Y Dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }
        if(mcSimboxZ.z <= (2*(mcSphereRadius+2*mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: Z Dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }

        long prenParticles;

    	// Master Slave Part
    	cout << " Master Process Invoked with rank "<< pRank<< endl;

    	// configuration file reader
    	const int LIMIT=30000;
        cout<< " Input Configuration File : " << mcInputFile << endl;

        ifstream fin(mcInputFile,std::ios_base::in);
        if (!fin.good()){
        	cout<< " Input file does not exist!" << endl;
        	exit(-1);
        }

        // local value holders from file
        long number; int type;
	double mass,epot;
	int accCnt,rejCnt;
	Vector3d position;

        char headerline[LIMIT];
        cout << " File header check " << endl;

        // crunch header part
        while (!fin.eof()){
              fin.getline(headerline,LIMIT,'\n');
              cout << headerline << endl;
              if(headerline[1]=='E'){
                break;
              }
        }

        // read and fill MC vector container
        while((!fin.eof())){ // end of file check
        	   // reading            // appending
               fin>>number;        mcNumbers.push_back(number);
               fin>>type;          mcTypes.push_back(type);
               fin>>mass;          mcMasses.push_back(mass);
               fin>>position.x;
               if (mcPBC.x){
            	   while (position.x < 0.) position.x += mcSimboxX.x;
               	   while (position.x >= mcSimboxX.x) position.x -= mcSimboxX.x;
               }
               mcPositions.push_back(position.x);
               fin>>position.y;
               if (mcPBC.y){
            	   while (position.y < 0.) position.y += mcSimboxY.y;
            	   while (position.y >= mcSimboxY.y) position.y -= mcSimboxY.y;
               }
               mcPositions.push_back(position.y);
               fin>>position.z;
               if (mcPBC.z){
            	   while (position.z < 0.) position.z += mcSimboxZ.z;
            	   while (position.z >= mcSimboxZ.z) position.z -= mcSimboxZ.z;
               }
               mcPositions.push_back(position.z);
               fin>>epot;          mcPotentialEnergies.push_back(epot);
 	       fin>>accCnt;	   mcAcc.push_back(accCnt);
	       fin>>rejCnt;	   mcRej.push_back(rejCnt);
               fin.get(); //Linebreak
               fin.peek(); //Check for eof
        }

        fin.clear();      // resetting bit states
        fin.close();




        cout << "  ======= Post container Check  =========== " << endl;
        cout << " mcNumbers           size : " << mcNumbers.size() << endl;
        cout << " mcTypes             size : " << mcTypes.size() << endl;
        cout << " mcMasses            size : " << mcMasses.size() << endl;
        cout << " mcPositions         size : " << mcPositions.size()/3 << endl;
        cout << " mcPotentialEnergies size : " << mcPotentialEnergies.size() << endl;

        // MonteCarlo cell size
        mcCellSize.x = (mcSphereRadius + (2 * mcSphereWallThickness));
        mcCellSize.y = (mcSphereRadius + (2 * mcSphereWallThickness));
        mcCellSize.z = (mcSphereRadius + (2 * mcSphereWallThickness));

        Vector3d radiusSet;
        radiusSet.x = mcSphereRadius;
        radiusSet.y = mcSphereWallThickness;
        radiusSet.z = mcSphereWallThickness;

        // creating domain and cell construction
	    ptrDomain = Domain::getDomainInstance(pRank,mccpuDimension,mcSimboxX,mcSimboxY,mcSimboxZ,
	   mcPBC,mcRealTypes,mcTotalTypes,radiusSet,mcZoneSize,mcCommunicator,mcCohesive,mcTemperature,mcnSteps,mcFlushInterval,mcMasterFlag,mcSeed,mcSamplingMode,mcStepsInBiasedBlock,mcStepsInLocalBlock,mcSphereRadiusTarget,mcLocalMovesRadius,mcCoveringTimes,mcSphereFile,mcCylinderRadius);

	    mcnParticles = mcNumbers.size();

	    // particle construction

	    ptrDomain->constructParticles(mcnParticles,mcNumbers,mcTypes,mcMasses,mcPositions,mcPotentialEnergies,mcAcc,mcRej);

	    prenParticles = ptrDomain->getnParticlesDomain();

        // evaluate total no of eligible Sampling sites and reference to eligible
        // particles are stored in the reference_wrapper based container sampleSites

	// Estimate eligible sampling sites
        ptrDomain->filterSamplingSites();

	//equilibrationStep 
        ptrDomain->calculateEquilibrationStep();	

        Particle         randomParticle;
        int     slavesCount              = (nProcess-1);
        vector<Particle> sphereHolder;
        vector<Particle> updatedSphere;

	int eqStep = ptrDomain->getEquilibrationStep();
	//int eqStep = 2;

	if ( ptrDomain->getSamplingMode()==2 || ptrDomain->getSamplingMode()==4 || ptrDomain->getSamplingMode()==6 ){

		if ( mcBiasedMovesPerformed == ptrDomain->getStepsInBiasedBlock() ){
			ptrDomain->setBiasedBlock(0);
			ptrDomain->setLocalBlock(1);

//cout << "Simulation starts at MC step " << mcNFinished << endl;
//cout << "The next step will start in LOCAL mode" << endl;
		}

		else if (mcBiasedMovesPerformed!=0 && mcLocalMovesPerformed==0){
			ptrDomain->setBiasedBlock(1);
			ptrDomain->setLocalBlock(0);
//cout << "Simulation starts at MC step " << mcNFinished << endl;
//cout << "The next step will start in BIASED mode" << endl;
		}

		if ( mcLocalMovesPerformed == ptrDomain->getStepsInLocalBlock() || mcNFinished == eqStep ){
			ptrDomain->setBiasedBlock(1);
			ptrDomain->setLocalBlock(0);
//cout << "Simulation starts at MC step " << mcNFinished << endl;
//cout << "The next step will start in BIASED mode" << endl;
		}

		else if (mcBiasedMovesPerformed==0 && mcLocalMovesPerformed!=0){
			ptrDomain->setBiasedBlock(0);
			ptrDomain->setLocalBlock(1);
//cout << "Simulation starts at MC step " << mcNFinished << endl;
//cout << "The next step will start in LOCAL mode" << endl;
		}

		if ( mcNFinished < eqStep && mcBiasedMovesPerformed==0 && mcLocalMovesPerformed==0){
			ptrDomain->setBiasedBlock(0);
			ptrDomain->setLocalBlock(0);
//cout << "Simulation starts at MC step " << mcNFinished << endl;
//cout << "The next step will start in UNBIASED mode" << endl;
		}

	}

        // Initialize Master

        Master master(mcTrialMoveType,mcSeed,mcNConflicts,mcNFinished,
        	   mcNSwaps,mcNAccepSwaps,mcNRejeSwaps,mcNAddition,mcNDeletion,
        	   mcNaddAccepted,mcNdelAccepted,mcNaddRejected,mcNdelRejected,
			   slavesCount,mcTotDeltaEpot,mcBiasedMovesPerformed,mcLocalMovesPerformed,ptrDomain,mcStatFile);

        double master_Start = clock();

        // execute Master process
        master.runAsMaster();

        auto master_End = clock();

        //todo Measuring Master simulation Wall time
        cout<<"Master Wall time : "<<float(master_End-master_Start)/CLOCKS_PER_SEC<<"  secs"<<endl;


	int nTrials               = master.getNTrials();
	long nConflicts           = master.getNConflicts();
	long nFinished            = master.getNFinished();
	long nSwaps               = master.getNSwaps();
	long nAccepSwaps          = master.getNAccepSwaps();
	long nRejeSwaps           = master.getNRejeSwaps();
	double totDeltaEpot       = master.getTotDeltaEpot();
	long biasedMovesPerformed = master.getBiasedMovesPerformed();
	long localMovesPerformed  = master.getLocalMovesPerformed();

   	string parameterFileName  = "paramFile_"+to_string(nFinished)+"_MCSteps.param";


		// invoking File Flush
		ptrDomain->writeCompleteConfiguration(mcOutputFile);
          	ptrDomain->writeParameterFile(parameterFileName,nTrials,nConflicts,nFinished,nSwaps,nAccepSwaps,nRejeSwaps,totDeltaEpot,biasedMovesPerformed,localMovesPerformed);

	    cout << "Particle Container Check postnParticles with MOD:" <<ptrDomain->getnParticlesDomain() << endl;
	    cout << " SIMULATION EXECUTION SUCCESSFULL"<<endl;
	    //------------------------------------------

	    if(prenParticles != ptrDomain->getnParticlesDomain()){
	    	cerr << " ERROR: Inconsistent no of particles after MonteCarlo sampling!! " << endl;
	    }
    }

    // --------------   Slave Part ----------------------
    else{
    	cout << " Slave Process Invoked with rank "<< pRank<< endl;
    	slave.runAsSlave();
    }

	MPI_Finalize();
	return 0;
}







