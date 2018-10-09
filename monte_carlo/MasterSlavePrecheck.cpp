/*
 * MasterSlavePrecheck.c++
 *
 */
/*
 * MasterSlave.c++
 *
 * Source file implementation of Master-Slave approach for Parallel MonteCarlo(MC) Sampling.
 * Configuration and MonteCarlo parameter files are read. Master process manages the Data
 * structure allocation, random selection, sphere construction and job allocation procedures. Constructed sphere domains
 * are distributed across slave processes, which invoke IMD Library interface for equilibration/relaxation
 * of sphere domains.Each slave process test the acceptance criterion and return its decision to Master process
 * upon completion. Upon acceptance part of the domain involved in sphere gets updated by Master process.
 * A prior screening is done to avoid overlapping sphere domains before getting allocated to slave processes.
 */

#include "mpi.h"
#include "constructSampleZone.h"
#include "localMD.h"
#include "acceptanceCheck.h"
#include "TrialMove.h"
#include "MonteCarloGlobals.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>

using namespace std;

int main(int argc,char* argv[]){

	//string inputFile = "inputFileOldlast.chkpt"; // simulation box configuration
    //string outFileName = "outputFileOldlast.chkpt";


	string inputFile   = "DiscSlab_screw_equil.00004.chkpt"; // simulation box configuration
	string outFileName = "DiscSlab_Output.chkpt";

    //string paramFile = "paramFile.chkpt";      // MonteCarlo Parameter file

	MPI_Init(&argc,&argv);

	int nProcess,pRank;
	bool finish = false;
	const int TAG_RUNJOB=0, TAG_TERMINATE=1, TAG_FREE=1, TAG_INITIALIZE=1;

    MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    // ----------------------  Master Process ---------------------- //
    if(pRank == 0){ // Rank-0 is designated as Master Process

    	int flag, source, freeFlag, freeSource;
    	MPI_Status status, freeStatus;
    	MPI_Request request;

    	//==============================================================
    	// input Parts (to be modified with suitable Param file reader)
    	// (1 1 1) for Master Slave approach
	    mccpuDimension.x = 1;
	    mccpuDimension.y = 1;
	    mccpuDimension.z = 1;

	    mcTotalTypes  = 9;
	    mcRealTypes   = 3;

//#X 	3.0082645700000000e+02 0.0000000000000000e+00 0.0000000000000000e+00
//#Y 	0.0000000000000000e+00 2.9889501100000001e+02 0.0000000000000000e+00
//#Z 	0.0000000000000000e+00 0.0000000000000000e+00 5.9362724000000000e+01

	    mcSimboxX.x = 3.00e+02; mcSimboxX.y = 0.0; mcSimboxX.z = 0.0;
	    mcSimboxY.x = 0.0; mcSimboxY.y = 3.00e+02; mcSimboxY.z = 0.0;
	    mcSimboxZ.x = 0.0; mcSimboxZ.y = 0.0; mcSimboxZ.z = 5.9362724000000000e+01;

//	    mcSimboxX.x = 9.7943498000000005e+01; mcSimboxX.y = 0.0; mcSimboxX.z = 0.0;
//	    mcSimboxY.x = 0.0; mcSimboxY.y = 1.0097804400000000e+02; mcSimboxY.z = 0.0;
//	    mcSimboxZ.x = 0.0; mcSimboxZ.y = 0.0; mcSimboxZ.z = 9.8937872999999996e+01;


        // Periodic Boundary Conditions
	    mcPBC.x = 0; mcPBC.y = 0; mcPBC.z = 1;

	    // window size
	    //mcZoneSize.x = 100; mcZoneSize.y = 100; mcZoneSize.z = 100;

	    // Appropriate parameter for Domain decomposition approach
	    mcZoneSize.x = 40; mcZoneSize.y = 40; mcZoneSize.z = 20;

	    // MC simulation and sampling data
	    //mcTemperature          =  15;
//	    mcSphereRadius         =  8.00000;
//	    mcSphereWallThickness  =  1.00000;

	    mcTemperature          =  0.02583;
	    mcSphereRadius         =  20.00000;
	    mcSphereWallThickness  =  5.00000;
        mcnSteps               =  7;
        //==============================================================

    	// Master Slave Part
    	cout << " Master Process Invoked with rank "<< pRank<< endl;

    	// configuration file reader
    	const int LIMIT=30000;
        cout<< " Input Configuration File : " << inputFile << endl;

        ifstream fin(inputFile,std::ios_base::in);

        // local value holders from file
        long number; int type;
        double mass,epot,eamRho;
        Vector3d position,velocity;
        long prenParticles, postnParticles;

        char headerline[LIMIT];
        cout << " File header check " << endl;

        // crunch header part
        while (fin.get() == '#'){
              fin.getline(headerline,LIMIT,'\n');
              cout << headerline << endl;
              continue;
        }

        // read and fill MC vector container

        while((!fin.eof())){ // end of file check
        	   // reading            // pushing
               fin>>number;        mcNumbers.push_back(number);
               fin>>type;          mcTypes.push_back(type);
               fin>>mass;          mcMasses.push_back(mass);
               fin>>position.x;    mcPositions.push_back(position.x);
               fin>>position.y;    mcPositions.push_back(position.y);
               fin>>position.z;    mcPositions.push_back(position.z);
               fin>>velocity.x; fin>>velocity.y; fin>>velocity.z;
               fin>>epot;          mcPotentialEnergies.push_back(epot);
               fin>>eamRho;
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
//        mcCellSize.x = 2*(mcSphereRadius + (2 * mcSphereWallThickness));
//        mcCellSize.y = 2*(mcSphereRadius + (2 * mcSphereWallThickness));
//        mcCellSize.z = 2*(mcSphereRadius + (2 * mcSphereWallThickness));

        mcCellSize.x = (mcSphereRadius + (2 * mcSphereWallThickness));
        mcCellSize.y = (mcSphereRadius + (2 * mcSphereWallThickness));
        mcCellSize.z = (mcSphereRadius + (2 * mcSphereWallThickness));

        Vector3d radiusSet;
        radiusSet.x = mcSphereRadius;
        radiusSet.y = mcSphereWallThickness;
        radiusSet.z = mcSphereWallThickness;
        int masterFlag = 1; // set to one for Master Slave approach !!

        // creating domain and cell construction
	    ptrDomain = Domain::getDomainInstance(pRank,mccpuDimension,mcSimboxX,mcSimboxY,mcSimboxZ,mcCellSize,
	    		mcPBC,mcRealTypes,radiusSet,mcZoneSize,mcCommunicator,masterFlag);

	    mcnParticles = mcNumbers.size();

	    // particle construction
	    ptrDomain->constructParticles(mcnParticles,mcNumbers,mcTypes,mcMasses,mcPositions,mcPotentialEnergies);

	    prenParticles = ptrDomain->getnParticlesDomain();

	    //---------------------------------------------------------
	    cout << "Particle Container Check prenParticles :" <<prenParticles << endl;

	    int Type;
	    long feCount=0, cCount=0, plCount=0;

	    vector<Particle>& currList = ptrDomain->getParticles();
	    // check part species type consistency
	    for(auto i=0; i<prenParticles; i++){

             Type = currList[i].getType();

             if(Type%mcRealTypes == 0){
	    			 feCount++; }

	    	 if(Type%mcRealTypes == 1){
	    			 cCount++; }

	         if(Type%mcRealTypes == 2){
	    			 plCount++; }

	    }

	    cout << " BEFORE UPDATE SPECIES COUNTER CHECK " <<endl;
	    cout << "  feCount  : " <<feCount << endl;
        cout << "  cCount   : " <<cCount << endl;
        cout << "  plCount  : " <<plCount<<endl;

	    vector<long> sampleSites;

	    long nCells,nParticles;
        int typeHold;
        nCells     = ptrDomain->getnCells();

        // evaluate total no of eligible Sampling sites and reference to eligible
        // particles are stored in the reference_wrapper based container sampleSites
        sampleSites = ptrDomain->filterSamplingSites();

        //---------------------------------------------------------
        cout << " *******************************************************" << endl;
        cout << " MAIN CHECK  sample Sites : " << sampleSites.size()<< endl;
        cout << "=======================================================" << endl;



        // random non-overlapping particle list
        vector<Particle> randomList;
        Particle         randomParticle;

        int     listSize              = (nProcess-1);
        double checkDistance          = (mcSphereRadius+2*mcSphereWallThickness);
        double checkDistanceSquare    = checkDistance*checkDistance;

	    // Get the results back from the slaves in no particular order
	    int    nFinished = 0;
        int    freeRank;
        bool   recDecision;
        bool   finished =false;

        int    receivedCount;
        int    currSlaves;
        int    count;
        int  particleCounter;
        long   nSimulation=0;

        vector<Particle> sphereHolder;
        vector<Particle> updatedSphere;
        Vector3d posHolder;

        int recbufind;
        long bufCounter, recBufSize;
        double tempId,tempType,tempMass,tempEpot;
        Vector3d tempPos;
        Particle particleHolder;
        double refZone1Distance;
        Vector3d tempPosition;
    	refZone1Distance  = (mcSphereRadius+mcSphereWallThickness+mcSphereWallThickness);
    	long preCounter, postCounter;

        // computation of masterSteps with available slave processes, total no of available sampling sites
        // and given MC steps.
        int masterSteps;
        masterSteps = (sampleSites.size()/(nProcess-1)) + (sampleSites.size()%(nProcess-1));

        cout << " Computed MasterSteps count :  " << masterSteps << endl;
        currSlaves    = listSize;

        double initData[2] = {static_cast<double> (mcRealTypes), mcTemperature};
        //
        //* Initialization part
        //*
        for(auto i=1; i<=(nProcess-1); i++){
        	// send initialize data
		    MPI_Isend(initData, 2, MPI_DOUBLE, i, TAG_INITIALIZE, MPI_COMM_WORLD,&request);
        }

        // Decision Variables
        int flipType,targetType, resType;

        // Decision Statistics
        long nAddition=0, nDeletion=0, naddAccepted=0, ndelAccepted=0, naddRejected=0, ndelRejected=0;
        vector<int> decisionList;

        //********************************************************************
        //                          Main LOOP
        //                   MonteCarlo Sampling part
        //********************************************************************
        // to be replaced with loop count of Master steps
        while(!finished){


        	// if last step-- corresponding steps has to be included

	        double** jobBuffer;
	        jobBuffer = new double* [currSlaves];
	        long* bufferSizeList;
	        bufferSizeList = new long [currSlaves];

	        if(randomList.size()>0)
	        	  randomList.clear();

	        randomList = ptrDomain->createRandomSampleList(currSlaves,checkDistanceSquare);

		    // sphere Construction
	        for(auto i=0; i<currSlaves; i++){

	        	bufCounter               = 0;
	        	randomParticle           = randomList[i];
	        	posHolder                = randomParticle.getPosition();

	        	// clear existing container
	        	if(sphereHolder.size()>0){
	        		sphereHolder.clear();
	        	}

	    	    // minimum boundary - Particle shift value
	    	    double refXmin = randomParticle.getPosition().x - (refZone1Distance);
	    	    double refYmin = randomParticle.getPosition().y - (refZone1Distance);
	    	    double refZmin = randomParticle.getPosition().z - (refZone1Distance);

	    		// Sphere Construction
	        	ptrDomain->createSphere(posHolder,randomParticle.getmcID(), sphereHolder);

        	    // Appropriate decision list!
	    		// add selected particle into list
	        	targetType = randomParticle.getType();
	        	// Switching placeholder into carbon
	        	if(targetType == 2){
	        		flipType = 1;
	        		decisionList.push_back(1);
	        		nAddition++;
	        	}
	        	// Switching carbon into placeholder
	        	if(targetType == 1){
	        		flipType = 2;
	        		decisionList.push_back(0);
	        		nDeletion++;
	        	}
	    		randomParticle.setType(flipType);

	    		// shift coordinates
	            tempPosition.x = randomParticle.getPosition().x - refXmin;
	            tempPosition.y = randomParticle.getPosition().y - refYmin;
	            tempPosition.z = randomParticle.getPosition().z - refZmin;

	            // adding chosen Particle
	            randomParticle.setPosition(tempPosition.x,tempPosition.y,tempPosition.z);
	    		sphereHolder.push_back(randomParticle);

	    		// Test Part
	    		// check sphere configuration
	    		// opening file stream object
	    		string filename_local = "sphereCheck"+to_string(mcPrank)+".chkpt";

	    		//ofstream fout(filename, ios_base::out);
	    		ofstream fout(filename_local, ios_base::out);
	    		long       fpTotal = sphereHolder.size();
	    		long fileNumber;
	    		int fileType;
	    		double fileMass, fileEpot;
	    		Vector3d filePos;

	    		// print IMD header -- check
	    		fout <<"#F A 1 1 1 3 0 0"<<endl;
	    		fout<<"#C number type mass x y z epot"<<endl;
	    		fout<<"#X "<<2*mcCellSize.x<<" 0"<<" 0 "<<endl;
	    		fout<<"#Y "<<"0 "<<2*mcCellSize.y<<" 0 "<<endl;
	    		fout<<"#Z "<<"0 "<<"0 "<<2*mcCellSize.z<<endl;
	    		fout<<"#E "<<endl;

	    		Particle fileParticle;

	    		for(long fp=0; fp<fpTotal; fp++){

	    			    	fileParticle = sphereHolder[fp];

	    			    	fileNumber   = fileParticle.getmcID();
	    			    	fileType     = fileParticle.getType();
	    			    	fileMass     = fileParticle.getMass();
	    			    	filePos      = fileParticle.getPosition();
	    			    	fileEpot     = fileParticle.getEpot();

	    			    	// file flush
	    			    	fout <<fileNumber
	    			    	            << "  " << fileType
	    			    	            << "  " << setprecision(6) << fileMass
	    			    	            << "  " << setprecision(6) << filePos.x
	    			    	            << "  " << setprecision(6) << filePos.y
	    								<< "  " << setprecision(6) << filePos.z
	    			    	            << "  " << setprecision(6) << fileEpot
	    			    	            << endl;
	    		}
	    		fout.close(); // closing outfile connection


	        	jobBuffer[i]             = new double [sphereHolder.size()*7];
	            bufferSizeList[i]        = sphereHolder.size()*7;

	        	// fill Sphere Particles into Job buffer!!
	        	for(auto j=0; j<sphereHolder.size(); j++){

	        		particleHolder         = sphereHolder[j];
	        		tempId                 = static_cast<double> (particleHolder.getmcID());
	                tempType               = static_cast<double> (particleHolder.getType());
	                tempMass               = particleHolder.getMass();
	        		tempPos                = particleHolder.getPosition();
	        		tempEpot               = particleHolder.getEpot();

	                *(jobBuffer[i]+bufCounter++)  = tempId;
	                *(jobBuffer[i]+bufCounter++)  = tempType;
	                *(jobBuffer[i]+bufCounter++)  = tempMass;
	                *(jobBuffer[i]+bufCounter++)  = tempPos.x;
	                *(jobBuffer[i]+bufCounter++)  = tempPos.y;
	                *(jobBuffer[i]+bufCounter++)  = tempPos.z;
	                *(jobBuffer[i]+bufCounter++)  = tempEpot;
	        	}

	        }

            //***********************************************************
	    	//*                  Communication Part
	        //*
	        //***********************************************************
	        count =0;
	        for(auto rank=1; rank<=(nProcess-1); rank++){

			    // send buffer size
			    MPI_Isend(&bufferSizeList[count], 1, MPI_LONG, rank, TAG_RUNJOB, MPI_COMM_WORLD,&request);

			    // send buffer content
			    MPI_Isend(jobBuffer[count], (bufferSizeList[count]), MPI_DOUBLE, rank, TAG_RUNJOB, MPI_COMM_WORLD,&request);

                count++;
	        }
	        receivedCount = 0;

	        //*************************************************************************************
	        //             Receiving results part
	        //*************************************************************************************

	        while(receivedCount<currSlaves){

	 	        MPI_Iprobe(MPI_ANY_SOURCE, TAG_RUNJOB, MPI_COMM_WORLD, &flag, &status);

	 	        if (flag){

	 	  	         source = status.MPI_SOURCE;
	 	  	         MPI_Recv(&recDecision, 1, MPI_C_BOOL, source,TAG_RUNJOB, MPI_COMM_WORLD, &status);

	 	  	         if(recDecision){

	 	  	        	MPI_Recv(&recBufSize, 1, MPI_LONG, source,TAG_RUNJOB, MPI_COMM_WORLD, &status);

	 	  	        	// allocate receive buffer
	    			    double* receiveBuffer;
	    			    receiveBuffer = new double [recBufSize];

	    			    //import buffer
	    			    MPI_Recv(receiveBuffer,recBufSize, MPI_DOUBLE, source,TAG_RUNJOB, MPI_COMM_WORLD, &status);

	    			    recbufind = 0;
	    			    particleCounter = recBufSize/7;

	    			    // clear existing container
	    			    if(updatedSphere.size()>0){
	    			    	updatedSphere.clear();
	    			    }
	    			    // fill back in particle container
	    			    for(auto i=0; i<particleCounter; i++){
	    			    	particleHolder.setmcID(receiveBuffer[recbufind++]);
	    			    	particleHolder.setType(receiveBuffer[recbufind++]);
                            particleHolder.setMass(receiveBuffer[recbufind++]);
	    			    	particleHolder.setPosition(receiveBuffer[recbufind++],receiveBuffer[recbufind++],
	    			    			                   receiveBuffer[recbufind++]);
	    			    	particleHolder.setEpot(receiveBuffer[recbufind++]);
	    			    	updatedSphere.push_back(particleHolder);
	    			    }

	    			    // identify the centre particle, shift its coordinates delete core sphere from domain and
	    			    // add the inner core list to Domain

	    			    // update part
	    			    randomParticle = randomList[source-1];

	    	        	resType = randomParticle.getType();

	    	        	// Acceptance of carbon addition
	    	        	if(resType == 2){
	    	        		naddAccepted++;
	    	        	}

	    	        	// Acceptance of carbon deletion
	    	        	if(resType == 1){
	    	        		ndelAccepted++;
	    	        	}


                        // Shifting sphere coordinates back to the natural coordinate system
	       		        ptrDomain->shiftSphere(updatedSphere,mcSimboxX,
	       		        			      mcSimboxY,mcSimboxZ,randomParticle.getPosition());

	       		        // precounter Part
	       		        preCounter = ptrDomain->getnParticlesDomain();
//	       		        cout <<" PRECOUNTER  :" << preCounter << endl;


	       			    long feCount=0, cCount=0, plCount=0;
	       			    // check part species type consistency

	       			    vector<Particle>& preParticles = ptrDomain->getParticles();

	    	    		// Upon acceptance update the configuration in Domain
	       		        ptrDomain->updateDomain(updatedSphere,randomParticle.getmcID() );

	    	    		//CP-------------------------------------------------------
	       			    feCount=0; cCount=0; plCount=0;
	       			    // check part species type consistency

	       			    vector<Particle>& postParticles = ptrDomain->getParticles();

	       		        postCounter = ptrDomain->getnParticlesDomain();

	    			    delete receiveBuffer;
	 	  	         }
	 	  	         else{
		 	  	        //cout << "Received REJECTION "<<recDecision<<" from Process :" << source << endl;

		    			    randomParticle = randomList[source-1];
		    	        	resType = randomParticle.getType();

		    	        	// Rejection of carbon addition
		    	        	if(resType == 2){
		    	        		naddRejected++;
		    	        	}

		    	        	// Rejection of carbon deletion
		    	        	if(resType == 1){
		    	        		ndelRejected++;
		    	        	}
	 	  	         }
	 	  	         receivedCount++;
	 	  	         nSimulation++;
	 	  	     }
	        }

		      //deleting internal buffers
		    for(auto i=0; i<currSlaves; i++){
		    	delete[] jobBuffer[i];
		    }
		    delete[] bufferSizeList;
		    delete[] jobBuffer;

  	      	if (nSimulation == mcnSteps){
  	      	 	 finished=true;
  	      	 	 cout << "============================================" << endl;
  	      	 	 cout << " Montecarlo Sampling Successful !! " << endl;

  	      	 	 cout << " Total MonteCarlo Steps : " << mcnSteps  << endl;
  	      	 	 cout << " No of addition         : " << nAddition << endl;
  	      	 	 cout << " No of deletion         : " << nDeletion << endl;
  	      	 	 cout << " No of accepted addition: " << naddAccepted << endl;
  	      	 	 cout << " No of rejected addition: " << naddRejected << endl;
  	      	     cout << " No of accepted deletion: " << ndelAccepted << endl;
  	      	     cout << " No of rejected deletion: " << ndelRejected << endl;
  	      	     cout << "============================================" << endl;

  	      	}
	    } // loop over masterSteps


	    //Termination of slave processes
	    for (auto i=1 ; i<nProcess ; ++i){
     	     MPI_Isend(NULL, 0, MPI_LONG, i, TAG_TERMINATE, MPI_COMM_WORLD,&request);
	    }

	    postnParticles = ptrDomain->getnParticlesDomain();
		feCount=0; cCount=0; plCount=0;
		vector<Particle>& finalParticles = ptrDomain->getParticles();

	    // new file writer part!!
	    ofstream fout(outFileName, ios_base::out);

	    // print IMD header -- check

	    fout <<"#F A 1 1 1 3 0 0"<< endl;
	    fout <<"#C number type mass x y z Epot"<< endl;
	    fout <<"#X "<<mcSimboxX.x<<" "<<mcSimboxX.y<<" "<<mcSimboxX.z<<" "<< endl;
	    fout <<"#Y "<<mcSimboxY.x<<" "<<mcSimboxY.y<<" "<<mcSimboxY.z<<" "<< endl;
	    fout <<"#Z "<<mcSimboxZ.x<<" "<<mcSimboxZ.y<<" "<<mcSimboxZ.z<<" "<< endl;
	    fout <<"#E "<< endl;

	    vector<Particle>& globalParticles = ptrDomain->getParticles();

	    // loop over particles in Domain
	    for(auto i=0; i<ptrDomain->getnParticlesDomain(); i++){

           number    = globalParticles[i].getNumber();
           type      = globalParticles[i].getType();
           mass      = globalParticles[i].getMass();
           position  = globalParticles[i].getPosition();
	       epot      = globalParticles[i].getEpot();

	       fout <<number
	       << "  " << type
	       << "  " << setprecision(6) << mass
	       << "  " << setprecision(6) << position.x
	       << "  " << setprecision(6) << position.y
	       << "  " << setprecision(6) << position.z
	       << "  " << setprecision(6) << epot
	       << endl;

	    }

	    fout.close(); // closing outfile connection


	    //Decision history during file writer part!!
	    string decFileName = "decisionHistory.txt";
	    ofstream dout(decFileName, ios_base::out);

	    dout <<"#   MC-Step               Decision "<< endl;

	    // loop over particles in Domain
	    for(auto i=0; i<decisionList.size(); i++){

	       dout <<"  "<<(i+1)
	       << "                           " << decisionList[i]
	       << endl;
	    }

	    dout.close(); // closing decision file connection

        //CP
	    cout << "Particle Container Check postnParticles with MOD:" <<postnParticles << endl;
	    //------------------------------------------

	    if(prenParticles != postnParticles){
	    	cerr << " ERROR: Inconsistent no of particles after MonteCarlo sampling!! " << endl;
	    }
    }

    // --------------   Slave Part ----------------------
    else{

    	cout << " Slave Process Invoked with rank "<< pRank<< endl;

    	int Master=0;        // Master Process
    	MPI_Status status;
    	MPI_Request request;

    	bool slaveDecision;
        bool finish;
    	bool slaveFree=true;
        bool initFlag=true;
        long bufferLength, buffInd, exportLength,expInd;
        int  particleCounter;
        int realTypes;
        double initData[2];
        double temperature;

        if(initFlag){
	      MPI_Recv(&initData, 2, MPI_DOUBLE, Master, TAG_INITIALIZE,
		       MPI_COMM_WORLD, &status);
	      initFlag = false;
        }

	    realTypes        = static_cast<int> (initData[0]);
        temperature      = initData[1];


        interfaceIMD        mdObject;
        nvtIMD              acceptObject(realTypes,temperature);
        Particle            particleObject;
    	vector<Particle>    inputSphere;
   	    vector<Particle>    updatedSphere;

   	    string fileName    = "sphereEquil_"+to_string(pRank)+".param";
   	    string sysCommand  = "cp sphereEquil.param "+fileName ;
   	    string remCommand  = "rm "+fileName;

   	    char *cpCommand    =  new char[sysCommand.length() + 1];
   	    char *rmCommand    =  new char[remCommand.length() + 1];

   	    sysCommand.copy(cpCommand, sysCommand.length(),0);
   	    remCommand.copy(rmCommand, remCommand.length(),0);

   	    system(cpCommand);

   	    // getting char* from string
   	    char *paramFile = new char[fileName.length() + 1];
   	    fileName.copy(paramFile, fileName.length(),0);

   	    mdObject.loadInputFile(paramFile);

    	while(slaveFree){

    		particleCounter=0;

    		// Receive buffer size
    	    MPI_Recv(&bufferLength, 1, MPI_LONG, Master, MPI_ANY_TAG,
    		       MPI_COMM_WORLD, &status);

    	    switch (status.MPI_TAG)
    		{
    		     case TAG_RUNJOB:

    			     // Allocate buffers
    			     double* localBuffer;
    			     localBuffer = new double [bufferLength];

    		    	 // Receive Buffers
    		    	 MPI_Recv(localBuffer, bufferLength, MPI_DOUBLE, Master, TAG_RUNJOB,
    		    		       MPI_COMM_WORLD, &status);

    		    	 particleCounter = bufferLength/7; // each particle has 7 attributes
    		    	 buffInd         = 0;

    		    	 // Fill particle container from buffer
    		    	 for(auto i=0; i<particleCounter; i++){
    		    		 //particleObject.setNumber(static_cast<long> (localBuffer[buffInd++]));
    		    		 particleObject.setmcID(static_cast<long> (localBuffer[buffInd++]));
                         particleObject.setType(static_cast<int> (localBuffer[buffInd++]));
                         particleObject.setMass(localBuffer[buffInd++]);
                         particleObject.setPosition(localBuffer[buffInd++],
                        		 localBuffer[buffInd++],localBuffer[buffInd++]);
                         particleObject.setEpot(localBuffer[buffInd++]);
                         inputSphere.push_back(particleObject);
    		    	 }

    		    	 // invoke Local MD call
    		    	 mdObject.runLocalMD(inputSphere,updatedSphere);

    		    	 // acceptance Call
    		    	 slaveDecision = acceptObject.isAccepted(inputSphere,updatedSphere);

    		    	 if(slaveDecision){
    		    	     exportLength  = updatedSphere.size()*7;
        		    	 // collect Results
                         double* exportBuffer;
                         exportBuffer = new double[exportLength];
                         expInd       = 0;

                         for(auto i=0; i<updatedSphere.size(); i++){
                             particleObject         = updatedSphere[i];
                        	 //exportBuffer[expInd++] = static_cast<double> (particleObject.getNumber());
                        	 exportBuffer[expInd++] = static_cast<double> (particleObject.getmcID());
                        	 exportBuffer[expInd++] = static_cast<double> (particleObject.getType());
                        	 exportBuffer[expInd++] = particleObject.getMass();
                        	 exportBuffer[expInd++] = particleObject.getPosition().x;
                        	 exportBuffer[expInd++] = particleObject.getPosition().y;
                        	 exportBuffer[expInd++] = particleObject.getPosition().z;
                        	 exportBuffer[expInd++] = particleObject.getEpot();
                         }

                         // send decision
                         MPI_Send(&slaveDecision, 1, MPI_C_BOOL, Master, TAG_RUNJOB, MPI_COMM_WORLD);

                         // send export buffer size
                         MPI_Send(&exportLength,1, MPI_LONG, Master, TAG_RUNJOB, MPI_COMM_WORLD);

                         // export buffer
                         MPI_Send(exportBuffer,exportLength, MPI_DOUBLE, Master, TAG_RUNJOB, MPI_COMM_WORLD);

                         delete exportBuffer;
    		    	 }

    		    	 else{
    		    		 exportLength           = 1;
                         expInd                 = 0;

                         // send decision
                         MPI_Send(&slaveDecision, 1, MPI_C_BOOL, Master, TAG_RUNJOB, MPI_COMM_WORLD);
    		    	 }

    		    	 delete localBuffer;
    		    	 inputSphere.clear();
    		    	 updatedSphere.clear();
    		         break;

    		     // Message tag to terminate Slave process
    		     case TAG_TERMINATE:
        		     finish    = true;
        		     slaveFree = false;
        		     break;
    		}
    	}

    	system(rmCommand);
    	delete[] cpCommand;
    	delete[] paramFile;
    }

	MPI_Finalize();
	return 0;
}




