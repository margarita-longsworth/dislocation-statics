/*
 * constructSampleZone.c++
 *
 */
#include "constructSampleZone.h"
#include <iostream>

constructSampleZone::constructSampleZone(Vector3d zoneSize,int seed){
    zoneSize_            = zoneSize;
    distributionXSample_ = std::uniform_real_distribution<double>(0.0,zoneSize_.x);
    distributionYSample_ = std::uniform_real_distribution<double>(0.0,zoneSize_.y);
    distributionZSample_ = std::uniform_real_distribution<double>(0.0,zoneSize_.z);
    generatorSample_     = std::mt19937(seed);
}
constructSampleZone::~constructSampleZone(void){

}

Vector3d constructSampleZone::getRandomPosition(){
	Vector3d randomPosition;
	randomPosition.x = distributionXSample_(generatorSample_);
	randomPosition.y = distributionYSample_(generatorSample_);
	randomPosition.z = distributionZSample_(generatorSample_);
	return randomPosition;
}

randomSampleZone::~randomSampleZone(void){

}
Vector3d constructSampleZone::getZoneSize(void){
	return zoneSize_;
}


long randomSampleZone::estimateSamplingSites(int speciesType,Domain* ptrDomain){

	Cell targetCell;
    int randomFactor[3]={0,0,0};
    long nSites=0, locSites=0;

    Vector3d zoneSize = getZoneSize();
    int pRank=0;
    pRank= ptrDomain->getDomainNo();

    cout << "pRank check " << pRank << endl;
    if(pRank==0){
        Vector3d randomPosition = getRandomPosition();

		randomFactor[0] = static_cast<int> (randomPosition.x);
		randomFactor[1] = static_cast<int> (randomPosition.y);
		randomFactor[2] = static_cast<int> (randomPosition.z);
		cout <<"randomPosition check "<< endl;
		cout <<" randomPosition.x  " << randomPosition.x << endl;
		cout <<" randomPosition.y  " << randomPosition.y << endl;
		cout <<" randomPosition.z  " << randomPosition.z << endl;

		cout <<" randomFactor[0]   " << randomFactor[0] << endl;
		cout <<" randomFactor[1]   " << randomFactor[1] << endl;
		cout <<" randomFactor[2]   " << randomFactor[2] << endl;

	}

	//MPI_Bcast(randomFactor,3,MPI_INT,0,mcCommunicator);

	VectorRange domainRange = ptrDomain->getDomainRange();

    cout << " Domain range check " << endl;
    cout << "Domain range.xmin   " <<domainRange.xMin << endl;
    cout << "Domain range.xmax   " <<domainRange.xMax << endl;
    cout << "Domain range.ymin   " <<domainRange.yMin << endl;
    cout << "Domain range.ymax   " <<domainRange.yMax << endl;
    cout << "Domain range.zmin   " <<domainRange.zMin << endl;
    cout << "Domain range.zmax   " <<domainRange.zMax << endl;

	// Random Centre of Mass
	centreOfMass_.x = (domainRange.xMax-domainRange.xMin)/randomFactor[0];
	centreOfMass_.y = (domainRange.yMax-domainRange.yMin)/randomFactor[1];
	centreOfMass_.z = (domainRange.zMax-domainRange.zMin)/randomFactor[2];

	cout <<" centreOfMass_.x   " << centreOfMass_.x << endl;
	cout <<" centreOfMass_.y   " << centreOfMass_.y << endl;
	cout <<" centreOfMass_.z   " << centreOfMass_.z << endl;

	IntVector3d cellGlobalcoordinate = ptrDomain->findParticleGlobalCoordinate(centreOfMass_);
    cout << " Cell Global coordinate " << endl;
    cout << "cellGlobalcoordinate.x  " <<cellGlobalcoordinate.x << endl;
    cout << "cellGlobalcoordinate.y  " <<cellGlobalcoordinate.y << endl;
    cout << "cellGlobalcoordinate.z  " <<cellGlobalcoordinate.z << endl;

	IntVector3d cellLocalcoordinate  = ptrDomain->computeCellLocalCoordinate(cellGlobalcoordinate);

	IntVector3d nblGlobal,nblLocal;
	Cell nbCell;
	Particle nbParticle;
    cout << " Cell local coordinate " << endl;
    cout << "cellLocalcoordinate.x  " <<cellLocalcoordinate.x << endl;
    cout << "cellLocalcoordinate.y  " <<cellLocalcoordinate.y << endl;
    cout << "cellLocalcoordinate.z  " <<cellLocalcoordinate.z << endl;
	targetCell = ptrDomain->findCellwithLocalCoordinate(cellLocalcoordinate);

	int rTypes = ptrDomain->getRealTypes();

	// loop over sample neighbor cell
	for(int i=0; i<targetCell.getSampleNBLSize(); i++){

		locSites=0;
        nblGlobal = targetCell.getSampleNeighbor(i);

        // for MPI communication make respective calls here
        // find CPU no from Global Coordinate

        nblLocal  = ptrDomain->computeCellLocalCoordinate(nblGlobal);
		nbCell    = ptrDomain->findCellwithLocalCoordinate(nblLocal);

		// cell level estimation of species
		nbCell.computeSpeciesCount(rTypes);
		locSites =nbCell.getSpeciesCount(speciesType);
		sitesAtCell_.push_back(locSites);

		nSites += locSites;
	}

	return nSites;
}

Particle randomSampleZone::findRandomParticle(int nSites,int particleType,Domain* ptrDomain){

	int type=0,randomIndex=0;
	Particle randomParticle;
	Cell targetCell;

	int residue=0,jumpFlag=0, currSites=0, placeFlag=0,particleFlag=0;

	distributionParticle_ = std::uniform_int_distribution<int>(0,nSites-1);

	randomIndex= distributionParticle_(generatorParticle_);
	int counter = 0;

	IntVector3d cellGlobalcoordinate = ptrDomain->findParticleGlobalCoordinate(centreOfMass_);
	IntVector3d cellLocalcoordinate  = ptrDomain->computeCellLocalCoordinate(cellGlobalcoordinate);

	IntVector3d nblGlobal,nblLocal;
	Cell nbCell;
	Particle nbParticle;
	targetCell = ptrDomain->findCellwithLocalCoordinate(cellLocalcoordinate);

	cout <<"CHOSEN PARTICLE TYPE "<<particleType<< endl;
	cout <<"Random Index : "<< randomIndex <<endl;

	int skipFlag=0;
	int i=0;

	// loop over sample neighbor cell
	while ((i<targetCell.getSampleNBLSize()) && (skipFlag!=1) ){
		   currSites += sitesAtCell_[i];
           jumpFlag = randomIndex/currSites;

    	   placeFlag    = currSites - sitesAtCell_[i];
           particleFlag = randomIndex - placeFlag;

//           cout << " OUTSIDE Check   :" << endl;
//           cout << " jumpFlag   :" << jumpFlag <<endl;
//           cout << " placeFlag  :" << placeFlag<< endl;
//           cout << " particleFlag   :" <<particleFlag<< endl;
//           cout << "particleType" << particleType << endl;
//           cout << " Random index : "<< randomIndex << endl;
//           cout << " counter  : "<< counter << endl;

           // cell location is identified // Check zero indexing issue
           if(jumpFlag==0){

        	   nblGlobal = targetCell.getSampleNeighbor(i);

               nblLocal  = ptrDomain->computeCellLocalCoordinate(nblGlobal);
        	   nbCell    = ptrDomain->findCellwithLocalCoordinate(nblLocal);

        	   counter=0; // new inclusion
               for(auto j=0; j<nbCell.getnParticlesCell(); j++){

                       nbParticle = nbCell.getParticle(j);

                       if((counter == particleFlag) && (nbParticle.getType()== particleType)){
                    	   randomParticle = nbParticle;
                           cout << " *** findRandomParticle **** "<< endl;
                           cout << " Preliminary Check   :" << endl;
                           cout << " jumpFlag   :" << jumpFlag <<endl;
                           cout << " placeFlag  :" << placeFlag<< endl;
                           cout << " particleFlag   :" <<particleFlag<< endl;
                           cout << "particleType" << particleType << endl;
                           cout << " counter  : "<< counter << endl;

                           cout << " Random Particle Identified " << endl;
                           cout << "nbParticle.getNumber() " <<nbParticle.getNumber()<< endl;
                           cout << "nbParticle.getType() " <<nbParticle.getType()<< endl;
                           cout << "nbParticle.getMass() " <<nbParticle.getMass()<< endl;
                           counter++;
                           skipFlag = 1;
                           break;
                       }
                       if(nbParticle.getType()== particleType){
                           counter++;
                       }
               }
           }
           i++;
	}
	return randomParticle;
}
