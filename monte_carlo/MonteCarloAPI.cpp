/*
 * MonteCarloAPI.c++
 *
 */


#include "MonteCarloHeaders.h"
#include <iostream>
#include   <string>
#include <fstream>
#include <iomanip>

using namespace std;


extern "C" void initializeMonteCarlo(int* mdcpuDimension, int mdTotalTypes,
		        int mdRealTypes,int* mdRestriction, double* mdSimbox,
				double mdTemperature,double mdSphereRadius,
				double mdSphereWallThickness,int mdTotalMCSteps,
				int* mdPBC,double* mdZone,MPI_Comm mdCommunicator, int debug){

	    debugFlag  = debug;

	    mccpuDimension.x = *(mdcpuDimension++);
	    mccpuDimension.y = *(mdcpuDimension++);
	    mccpuDimension.z = *(mdcpuDimension++);

	    mcTotalTypes  = mdTotalTypes;
	    mcRealTypes   = mdRealTypes;

	    mcSimboxX.x = *(mdSimbox++); mcSimboxX.y = *(mdSimbox++); mcSimboxX.z = *(mdSimbox++);
	    mcSimboxY.x = *(mdSimbox++); mcSimboxY.y = *(mdSimbox++); mcSimboxY.z = *(mdSimbox++);
	    mcSimboxZ.x = *(mdSimbox++); mcSimboxZ.y = *(mdSimbox++); mcSimboxZ.z = *(mdSimbox++);

	    mcRestriction = new IntVector3d[mcTotalTypes];

	    // restriction vector for all element types
	    int count=0;
	    for(int i=0; i<mcTotalTypes; i++){
	    	mcRestriction[i].x = *(mdRestriction+count); count++;
	    	mcRestriction[i].y = *(mdRestriction+count); count++;
	    	mcRestriction[i].z = *(mdRestriction+count); count++;
	    }
        // Periodic Boundary Conditions
	    mcPBC.x = mdPBC[0]; mcPBC.y = mdPBC[1]; mcPBC.z = mdPBC[2];

	    // window size
	    mcZoneSize.x = mdZone[0]; mcZoneSize.y = mdZone[1];	mcZoneSize.z = mdZone[2];

	    // MC simulation and sampling data
	    mcTemperature          =  mdTemperature;
	    mcSphereRadius         =  mdSphereRadius;
	    mcSphereWallThickness  =  mdSphereWallThickness;
        mcnSteps               =  mdTotalMCSteps;
        mcCommunicator         =  mdCommunicator;

        // Perform dimension check (Sphere Construction)
        if(mcSimboxX.x <= (2*(mcSphereRadius+mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: X dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }
        if(mcSimboxY.y <= (2*(mcSphereRadius+mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: Y Dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }
        if(mcSimboxZ.z <= (2*(mcSphereRadius+mcSphereWallThickness)) ){
        	cerr << " MonteCarlo: Z Dimension is too small for Sphere Construction."<<endl;
                exit(1);
        }

        // MonteCarlo cell size
        mcCellSize.x = mcSphereRadius + mcSphereWallThickness;
        mcCellSize.y = mcSphereRadius + mcSphereWallThickness;
        mcCellSize.z = mcSphereRadius + mcSphereWallThickness;

	    MPI_Comm_size(mcCommunicator, &mcnProcess);
	    MPI_Comm_rank(mcCommunicator, &mcPrank);

        Vector3d radiusSet;
        radiusSet.x = mcSphereRadius;
        radiusSet.y = mcSphereWallThickness;
        radiusSet.z = mcSphereWallThickness;

        int masterFlag = 0;

	    ptrDomain = Domain::getDomainInstance(mcPrank,mccpuDimension,mcSimboxX,mcSimboxY,mcSimboxZ,mcCellSize,
	    		mcPBC,mcRealTypes,radiusSet,mcZoneSize,mcCommunicator,masterFlag);

        if(debugFlag){
        	cout<<"*****************************************************"<<endl;
        	cout<<"*****************************************************"<<endl;
        	cout<<"Parameters passing check -- initializeMonteCarlo"     << endl;
			cout<<"*****************************************************"<<endl;
        	cout<<"mccpuDimension [ :"<< mccpuDimension.x<<" "
        			                  << mccpuDimension.y<<" "
									  << mccpuDimension.z<<" "
        			                  <<"]"<<endl;
        	cout<<"mcSimboxX size : [ " << mcSimboxX.x << " "
        			                    << mcSimboxX.y << " "
										<< mcSimboxX.z << " "
        			                << "]"<< endl;
        	cout<<"mcSimboxY size : [ " << mcSimboxY.x << " "
        			                    << mcSimboxY.y << " "
										<< mcSimboxY.z << " "
        			                << "]"<< endl;
        	cout<<"mcSimboxZ size : [ " << mcSimboxZ.x << " "
        			                    << mcSimboxZ.y << " "
										<< mcSimboxZ.z << " "
        			                << "]"<< endl;

        	cout<<"mcZonesize  : [ " << mcZoneSize.x << " "
        			                    << mcZoneSize.y << " "
										<< mcZoneSize.z << " "
        			                << "]"<< endl;


        	for( auto i=0; i<mcTotalTypes; i++){

            	cout<<"Restriction ["<< i << "] : "<< mcRestriction[i].x <<" "
            			                           << mcRestriction[i].y <<" "
												   << mcRestriction[i].z <<" "
            			                           << endl;

        	}

        	cout<< "mcTemperature  :"            << mcTemperature  << endl;
        	cout<< "mcSphereRadius :"            << mcSphereRadius << endl;
        	cout<< "mcSphereWallThickness :"     << mcSphereWallThickness<< endl;
        	cout<< "mcnSteps : "                 <<mcnSteps << endl;

        	cout<<"===============================================" << endl;
        	cout<<"           sample neighbor check               " << endl;
        	cout<<"===============================================" << endl;

        	for(auto i=0; i<ptrDomain->getCell(0).getSampleNBLSize(); i++){
        		cout <<"Neighbor ["<<i<<" ]  Glob coordinate : "<<ptrDomain->getCell(0).getSampleNeighbor(i).x<<" "
        				                                        <<ptrDomain->getCell(0).getSampleNeighbor(i).y<<" "
																<<ptrDomain->getCell(0).getSampleNeighbor(i).z<<endl;
        	}

        	cout<<"===============================================" << endl;

        	cout<<"===============================================" << endl;
        	cout<<"           sphere neighbor check               " << endl;
        	cout<<"===============================================" << endl;

        	for(auto i=0; i<ptrDomain->getCell(0).getSphereNBLSize(); i++){
        		cout <<"Sphere Neighbor ["<<i<<" ]  Glob coordinate : "<<ptrDomain->getCell(0).getSphereNeighbor(i).x<<" "
        				                                        <<ptrDomain->getCell(0).getSphereNeighbor(i).y<<" "
																<<ptrDomain->getCell(0).getSphereNeighbor(i).z<<endl;
        	}

        	cout<<"===============================================" << endl;

        }

}

extern "C" void packConfigurationToMonteCarlo(long mdmcnParticles,long *mdmcNumbers,int *mdmcTypes,
		double  *mdmcMasses,double  *mdmcPositions , double *mdmcPotentialEnergies){


	mcnParticles = mdmcnParticles; // get total particles from MD
	long m=0;
	for(long i=0; i<mcnParticles; i++){

		// filling data container with corresponding values
		mcNumbers.push_back(mdmcNumbers[i]);
		mcTypes.push_back(mdmcTypes[i]);
		mcMasses.push_back(mdmcMasses[i]);
 		mcPositions.push_back(mdmcPositions[m++]); // holds x
 		mcPositions.push_back(mdmcPositions[m++]); // holds y
 		mcPositions.push_back(mdmcPositions[m++]); // holds z
		mcPotentialEnergies.push_back(mdmcPotentialEnergies[i]);
	}

	// check for unallocated vector
	if(mcNumbers.empty()||mcTypes.empty()||mcMasses.empty()||mcPositions.empty()||mcPotentialEnergies.empty()){
		std::cerr<<"MONTECARLO: problem allocating resources for MD data"<<endl;
	}
	else if((mcNumbers.size()==mcnParticles)&&(mcTypes.size()==mcnParticles)&&(mcMasses.size()==mcnParticles)
			&&((mcPositions.size()/3)==mcnParticles) && (mcPotentialEnergies.size()==mcnParticles)){
		std::cout<<" ------------------------------------------- "<<endl;
		std::cout<<" Process id : "<<mcPrank<<endl;
		std::cout<<" Total particles received : " <<mcnParticles <<endl;
		std::cout<<"MONTECARLO: Resource allocation successful " <<endl;
		std::cout<<" ------------------------------------------- "<<endl;
    }
	else {
		std::cerr<<"MONTECARLO: unknown error with allocation" <<endl;
	}


}

extern "C" void performMonteCarlo(int mdPrank,long *mdnParticles,long **mdNumbers,int **mdTypes,
		double **mdMasses,double **mdPositions){

	cout<<" MONTECARLO:Performing MonteCarlo simulation "<<endl;

    ptrDomain->constructParticles(mcnParticles,mcNumbers,mcTypes,mcMasses,mcPositions,mcPotentialEnergies);

    Particle randomParticle;
    Decision decision;
    bool isAcceptance;
    int nTrialmoves = 2;
    long nSites=0;
    int seed = 38;
    double refZone1Distance;
    Vector3d tempPosition;
	refZone1Distance  = (mcSphereRadius+mcSphereWallThickness+mcSphereWallThickness);

    vector<Particle>    inputSphere;
    vector<Particle>    updatedSphere;

    addDeleteCarbon     trialMove(nTrialmoves);
    randomSampleZone    randomZone(mcZoneSize,seed);
    interfaceIMD        mdObject;
    nvtIMD               acceptObject(mcRealTypes,mcTemperature);

    string fileName = "sphereEquil_"+to_string(mcPrank)+".param";

    // getting char* from string

    char *paramFile = new char[fileName.length() + 1];
    fileName.copy(paramFile, fileName.length(),0);

    mdObject.loadInputFile(paramFile);

    // Monte Carlo sampling Loop

    for(auto i=0; i<mcnSteps; i++){

    	//decision       = trialMove.getDecision();

    	// decision hard coded part
    	decision.flipType   = 1;
    	decision.targetType = 2;

    	cout <<"////////////////////////////////////////////////"<< endl;

    	cout << " Decision Check  "<<endl;
    	cout << " decision.flipType " <<decision.flipType << endl;
    	cout << " decision.targetType" << decision.targetType << endl;

    	nSites = randomZone.estimateSamplingSites(decision.targetType,ptrDomain);
    	cout << " Random no of sites " << nSites << endl;
    	cout <<"////////////////////////////////////////////////"<< endl;

    	if(nSites!=0){
    		randomParticle = randomZone.findRandomParticle(nSites,decision.targetType,ptrDomain);

    	    // minimum boundary - Particle shift value
    	    double refXmin = randomParticle.getPosition().x - (refZone1Distance);
    	    double refYmin = randomParticle.getPosition().y - (refZone1Distance);
    	    double refZmin = randomParticle.getPosition().z - (refZone1Distance);

    		cout << " selected random particle info " << endl;
    		cout << " randomParticle.getNumber() "<<randomParticle.getNumber()<<endl;
    		cout << " randomParticle.getType() " <<randomParticle.getType()<<endl;
    		cout << " randomParticle.getMass() "<< randomParticle.getMass()<<endl;
    		cout << " randomParticle.getPosition() [ " << randomParticle.getPosition().x<<" "
    		    			                                   << randomParticle.getPosition().y<<" "
    														   << randomParticle.getPosition().z<<" ]"<<endl;
    		cout << " randomParticle.getEpot() "<<randomParticle.getEpot()<< endl;
    		cout <<"********************************************"<<endl;

    		//inputSphere  = ptrDomain->createSphere(randomParticle.getPosition());

    		inputSphere  = ptrDomain->createSphere(randomParticle.getPosition(),randomParticle.getNumber());

    		// add selected particle into list
    		randomParticle.setType(decision.flipType);
    		// shift coordinates
            tempPosition.x = randomParticle.getPosition().x - refXmin;
            tempPosition.y = randomParticle.getPosition().y - refYmin;
            tempPosition.z = randomParticle.getPosition().z - refZmin;

            randomParticle.setPosition(tempPosition.x,tempPosition.y,tempPosition.z);
    		inputSphere.push_back(randomParticle);


    		//# also update mass if necessary

    		cout << " Sphere construction Check" << endl;
    		cout << " Total particles in Sphere : "<< inputSphere.size()<< endl;

    		// check sphere configuration
    	    // opening file stream object
    	    string filename_local = "MC_sphere_config_p_"+to_string(mcPrank)+".chkpt";

    	    //ofstream fout(filename, ios_base::out);
    	    ofstream fout(filename_local, ios_base::out);
    	    long fpTotal = inputSphere.size();

    	    long fileNumber;
    	    int fileType;
    	    double fileMass, fileEpot;
    	    Vector3d filePos;

    	    // print IMD header -- check
    	    fout <<"#F A 1 1 1 3 0 0"<<endl;
    	    fout<<"#C number type mass x y z epot"<<endl;
    	    fout<<"#X 24.009552 0 0 "<<endl;
    	    fout<<"#Y 0 24.009552 0 "<<endl;
    	    fout<<"#Z 0 0 24.009552 "<<endl;
    	    fout<<"#E "<<endl;

    	    Particle fileParticle;

    	    for(long fp=0; fp<fpTotal; fp++){

    	    	    fileParticle = inputSphere[fp];
    	        	fileNumber   = fileParticle.getNumber();
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

   		    mdObject.runLocalMD(inputSphere,updatedSphere);

   		    // copy CALL
   		    //mdObject.runLocalMD(inputSphere,updatedSphere);

    		isAcceptance  = acceptObject.isAccepted(inputSphere,updatedSphere);

    		if(isAcceptance){
    		     ptrDomain->shiftSphere(updatedSphere,mcSimboxX,
    		        			      mcSimboxY,mcSimboxZ,randomParticle.getPosition());
    		     // Old implementation
    		     //ptrDomain->updateDomain(updatedSphere);
    		     ptrDomain->updateDomain(updatedSphere,randomParticle);
    		}
  	        inputSphere.clear();
		    updatedSphere.clear();
    	}

    	else{
    		cout<<"Sampling sites for chosen particle type not exist"<< endl;
    	}

    }

    delete paramFile;

    // clear existing MC container
    mcNumbers.clear();
    mcTypes.clear();
    mcMasses.clear();
    mcPositions.clear();
    mcPotentialEnergies.clear();

    // fill back MC container
    Particle particleObj;
    Cell     cellObj;

    for(auto i=0; i<ptrDomain->getnCells(); i++){
        cellObj = ptrDomain->getCell(i);

    	for(auto j=0; j<cellObj.getnParticlesCell(); j++){

    		particleObj = cellObj.getParticle(j);

    		mcNumbers.push_back(particleObj.getNumber());
    		mcTypes.push_back(particleObj.getType());
    		mcMasses.push_back(particleObj.getMass());
    		mcPositions.push_back(particleObj.getPosition().x);
    		mcPositions.push_back(particleObj.getPosition().y);
    		mcPositions.push_back(particleObj.getPosition().z);
    		mcPotentialEnergies.push_back(particleObj.getEpot());
    	}

    }

    long updateTotal = mcNumbers.size(); // updated no of particles

    cout<< "MC CONTAINER -- PARTICLES COUNT : " <<updateTotal << endl;

	mdPrank         = mcPrank;               // process rank
	mdnParticles    = &updateTotal;          // current total particles
	*mdNumbers      = mcNumbers.data();      // atom id
    *mdTypes        = mcTypes.data();        // atom types
	*mdMasses       = mcMasses.data();       // atom mass
    *mdPositions    = mcPositions.data();    // atom positions

    if(debugFlag){
//    	cout<<"========================================"<<endl;
//    	cout<<" Post Particles Check"<< endl;
//    	cout<<"========================================"<<endl;
//
//    	cout << " Particle contribution check " << endl;
//
//    	for( auto i=0; i<ptrDomain->getnCells();  i++){
//    		cout << "Cell No   : "<<i<<" and Particles Count : "<<ptrDomain->getCell(i).getnParticlesCell()<<endl;
//    	}
//    	cout<<"========================================"<<endl;

//    	for( auto i=0; i<ptrDomain->getCell(0).getnParticlesCell();  i++){
//            cout << "Particle No   : "<<ptrDomain->getCell(0).getParticle(i).getNumber()<<endl;
//            cout << "Particle Type : "<<ptrDomain->getCell(0).getParticle(i).getType()<<endl;
//            cout << "Particle mass : "<<ptrDomain->getCell(0).getParticle(i).getMass()<<endl;
//            cout << "Particle position : "<<ptrDomain->getCell(0).getParticle(i).getPosition().x<<" "
//            		                      <<ptrDomain->getCell(0).getParticle(i).getPosition().y<<" "
//										  <<ptrDomain->getCell(0).getParticle(i).getPosition().z<<" "
//            		                              <<endl;
//            cout << "Particle Epot : "<<ptrDomain->getCell(0).getParticle(i).getEpot()<<endl;
//
//            cout<<"**********************************************"<<endl;
//    	}

    	cout<<"========================================"<<endl;
    }

    cout<<" MONTECARLO: MonteCarlo simulation completed !"<<endl;
}


extern "C" void importPotentialEnergy(double **mdPotentialEnergies){

}

extern "C" void cleanMonteCarlo(){

	cout<<" ------------------------------------------- "<<endl;
	cout<<"MONTECARLO : Cleaning Data structures "<<endl;

	//ptrDomain->~Domain();

	mcNumbers.clear();
	mcTypes.clear();
	mcMasses.clear();
	mcPositions.clear();
	mcVelocities.clear();
	mcPotentialEnergies.clear();

	cout<<  " MONTECARLO : Post Clean Check      " <<endl;
	cout<<  " size : mcNumbers             " << static_cast<int> (mcNumbers.size())<<endl;
	cout<<  " size : mcTypes               " << static_cast<int> (mcTypes.size())  << endl;
	cout<<  " size : mcMasses              " << static_cast<int> (mcMasses.size()) << endl;
	cout<<  " size : mcPositions           " << static_cast<int> (mcPositions.size())<< endl;
	cout<<  " size : mcVelocities          " << static_cast<int> (mcVelocities.size())<< endl;
	cout<<  " size : mcPotentialEnergies   " << static_cast<int> (mcPotentialEnergies.size())<< endl;
	cout<<" ------------------------------------------- "<<endl;

}

extern "C" void shutDownMonteCarlo(){
    cout<<"MONTECARLO : shutting Down MonteCarlo"<<endl;
}
