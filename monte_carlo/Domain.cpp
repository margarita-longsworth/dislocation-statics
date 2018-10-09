/*
 * Domain.c++
 *
 */
#include <iostream>
#include "Domain.h"

using namespace std;

Domain* Domain::getDomainInstance(int domainNo,IntVector3d cpuDimension,Vector3d boxSizeX,
        Vector3d boxSizeY, Vector3d boxSizeZ, IntVector3d PBC,
		int realTypes,int totalTypes,Vector3d  sphereLayerRadii,Vector3d zoneSize,MPI_Comm communicator,
		double cohesiveEnergy, double temperature,long nMCSteps,long writeInterval,int masterFlag,int userSeed, int SamplingMode, long StepsInBiasedBlock, long StepsInLocalBlock,double sphereRadiusTarget, int localMovesRadius,int coveringTimes,string mcSphereFile,double cylinderRadius){

	static Domain objDomain(domainNo,cpuDimension,boxSizeX,
            boxSizeY,boxSizeZ,PBC,realTypes,totalTypes,sphereLayerRadii,zoneSize,communicator,
			cohesiveEnergy,temperature,nMCSteps,writeInterval,masterFlag,userSeed,SamplingMode,StepsInBiasedBlock,StepsInLocalBlock,sphereRadiusTarget,localMovesRadius,coveringTimes,mcSphereFile,cylinderRadius);

	return &objDomain;
}

Domain::Domain(int domainNo,IntVector3d cpuDimension,Vector3d boxSizeX,
        Vector3d boxSizeY, Vector3d boxSizeZ, IntVector3d PBC,
		int realTypes,int totalTypes,Vector3d  sphereLayerRadii,Vector3d zoneSize,MPI_Comm communicator,
		double cohesiveEnergy, double temperature,long nMCSteps,long writeInterval,int masterFlag,int userSeed, int SamplingMode, long StepsInBiasedBlock, long StepsInLocalBlock, double sphereRadiusTarget, int localMovesRadius, int coveringTimes, string mcSphereFile, double cylinderRadius){

	cout << "Domain Constructor Invoked " << endl;

    // Seeding number generator
    std::random_device rd;

    // use random Seed if user seed is not provided
    if(userSeed==0){
    	userSeed_ = rd();
    	generator_.seed(rd());
    }
    else{
    	userSeed_ = userSeed;
    	generator_.seed(userSeed);
    }

	domainNo_         = domainNo;
	cpuDimension_     = cpuDimension;
	sphereLayerRadii_ = sphereLayerRadii;
    masterFlag_       = masterFlag;

    totalSizeX_       = boxSizeX;
    totalSizeY_       = boxSizeY;
    totalSizeZ_       = boxSizeZ;

	domainSizeX_.x = boxSizeX.x / cpuDimension.x;
	domainSizeX_.y = boxSizeX.y / cpuDimension.x;
	domainSizeX_.z = boxSizeX.z / cpuDimension.x;

	domainSizeY_.x = boxSizeY.x / cpuDimension.y;
	domainSizeY_.y = boxSizeY.y / cpuDimension.y;
	domainSizeY_.z = boxSizeY.z / cpuDimension.y;

	domainSizeZ_.x = boxSizeZ.x / cpuDimension.z;
	domainSizeZ_.y = boxSizeZ.y / cpuDimension.z;
	domainSizeZ_.z = boxSizeZ.z / cpuDimension.z;

	double minCellSize = sphereLayerRadii_.x+sphereLayerRadii_.y+sphereLayerRadii_.z;

	double globX = boxSizeX.x / minCellSize;
	double globY = boxSizeY.y / minCellSize;
	double globZ = boxSizeZ.z / minCellSize;

	globalCellArray_.x = static_cast<int> (globX);
	globalCellArray_.y = static_cast<int> (globY);
	globalCellArray_.z = static_cast<int> (globZ);

	// Update Cell Size
    cellSize_.x = boxSizeX.x/globalCellArray_.x;
    cellSize_.y = boxSizeY.y/globalCellArray_.y;
    cellSize_.z = boxSizeZ.z/globalCellArray_.z;

	domainCellArray_.x = globalCellArray_.x/cpuDimension.x ;
	domainCellArray_.y = globalCellArray_.y/cpuDimension.y ;
	domainCellArray_.z = globalCellArray_.z/cpuDimension.z;

	if( domainCellArray_.x != (globalCellArray_.x/cpuDimension.x))
		std::cerr<<" Incompatible global and domain cell array dimension "<<endl;
	if( domainCellArray_.y != (globalCellArray_.y/cpuDimension.y))
		std::cerr<<" Incompatible global and domain cell array dimension "<<endl;
	if( domainCellArray_.z != (globalCellArray_.z/cpuDimension.z))
		std::cerr<<" Incompatible global and domain cell array dimension "<<endl;

	int grid_coord[3];

	if(masterFlag_){
		grid_coord[0]=0; grid_coord[1]=0; grid_coord[2]=0;
	}
	else{
		MPI_Cart_coords(communicator,domainNo,3,grid_coord);
	}

	domainPosition_.x = grid_coord[0];
	domainPosition_.y = grid_coord[1];
	domainPosition_.z = grid_coord[2];

	domainRange_.xMin = (boxSizeX.x/cpuDimension.x) * domainPosition_.x;
	domainRange_.xMax = (boxSizeX.x/cpuDimension.x) * (domainPosition_.x+1);
	domainRange_.yMin = (boxSizeY.y/cpuDimension.y) * domainPosition_.y;
	domainRange_.yMax = (boxSizeY.y/cpuDimension.y) * (domainPosition_.y+1);
	domainRange_.zMin = (boxSizeZ.z/cpuDimension.z) * domainPosition_.z;
	domainRange_.zMax = (boxSizeZ.z/cpuDimension.z) * (domainPosition_.z+1);

    cout << " cellSize_ check " << endl;
    cout << " cellSize_.x    :  " <<  cellSize_.x << endl;
    cout << " cellSize_.y    :  " <<  cellSize_.y << endl;
    cout << " cellSize_.z    :  " <<  cellSize_.z << endl;
    cout << " ==================================="<<endl;

    cout << " domainSize check " << endl;
    cout << " domainSizeX_    :  [ " <<  domainSizeX_.x <<" "<<domainSizeX_.y<<" "<<domainSizeX_.z<<" "<<"]"<< endl;
    cout << " domainSizeY_    :  [ " <<  domainSizeY_.x <<" "<<domainSizeY_.y<<" "<<domainSizeY_.z<<" "<<"]"<< endl;
    cout << " domainSizeZ_    :  [ " <<  domainSizeZ_.x <<" "<<domainSizeZ_.y<<" "<<domainSizeZ_.z<<" "<<"]"<< endl;
    cout << " ==================================="<<endl;
    cout << " ==================================="<<endl;

    cout << " globalCellArray_ check " << endl;
    cout << " globalCellArray_.x    :  " <<  globalCellArray_.x << endl;
    cout << " globalCellArray_.y    :  " <<  globalCellArray_.y << endl;
    cout << " globalCellArray_.z    :  " <<  globalCellArray_.z << endl;
    cout << " ==================================="<<endl;

    cout << " domainCellArray_ check " << endl;
    cout << " domainCellArray_.x    :  " <<  domainCellArray_.x << endl;
    cout << " domainCellArray_.y    :  " <<  domainCellArray_.y << endl;
    cout << " domainCellArray_.z    :  " <<  domainCellArray_.z << endl;
    cout << " ==================================="<<endl;

    cout << " domainPosition_ check " << endl;
    cout << " domainPosition_.x    :  " <<  domainPosition_.x << endl;
    cout << " domainPosition_.y    :  " <<  domainPosition_.y << endl;
    cout << " domainPosition_.z    :  " <<  domainPosition_.z << endl;
    cout << " ==================================="<<endl;

    cout << " domainRange_ check " << endl;
    cout << " domainRange_.xMin    :  " <<  domainRange_.xMin << endl;
    cout << " domainRange_.xMax    :  " <<  domainRange_.xMax << endl;
    cout << " domainRange_.yMin    :  " <<  domainRange_.yMin << endl;
    cout << " domainRange_.yMax    :  " <<  domainRange_.yMax << endl;
    cout << " domainRange_.zMin    :  " <<  domainRange_.zMin << endl;
    cout << " domainRange_.zMax    :  " <<  domainRange_.zMax << endl;

    cout << " ==================================="<<endl;

    PBC_              = PBC;
    realTypes_        = realTypes;
    totalTypes_       = totalTypes;
    zoneSize_         = zoneSize;
    communicator_     = communicator;
    cohesiveEnergy_   = cohesiveEnergy;
    temperature_      = temperature;
    mcSphereFile_     = mcSphereFile;

    totalMCSteps_     = nMCSteps;
    writeMCinterval_  = writeInterval;

    constructCells();
    createSphereNeighbors();

    samplingMode_	= SamplingMode;
    stepsInBiasedBlock_ = StepsInBiasedBlock;
    stepsInLocalBlock_  = StepsInLocalBlock;
    Kappa_              = 0; 
    KappaR_             = 0; 
    sphereRadiusTarget_ = sphereRadiusTarget; 
    localMovesRadius_ = localMovesRadius; 
    cylinderRadius_ = cylinderRadius; 
    coveringTimes_      = coveringTimes;

    // Domain decomposition approach
    if(masterFlag_ == 0)
            createSampleNeighbors();

    // Otherwise MasterSlave approach
}

Domain::~Domain(void){
	cout << "Domain Destructor Invoked " << endl;
	cells_.clear();
}

void Domain::addParticle(Particle p){
	    particles_.push_back(p);
}
void Domain::setParticle(Particle p,long index){
	      particles_.at(index)=p;
}
void Domain::deleteParticle(long index){

	      // assign last particle at index to be deleted
	      particles_[index] = particles_.back();

	      // delete last element
	      particles_.pop_back();
}

const vector<Particle>&   Domain::getParticles(void) const{
        return particles_;
}

void Domain::computeCellSpeciesCount(int cellNo){

	int* counter;
	counter = new int[realTypes_];

	int value=0;
    long particleID;
    int  refType;

	// initialization
	for(auto j=0; j<realTypes_; j++){
		counter[j]=0;
	}

	long total                 = cells_[cellNo].getnParticlesCell();

	for(auto i=0; i<total; i++){
		 particleID        = cells_[cellNo].getParticleID(i);
		 refType           = particles_[particleID].getType();
		 value             = refType%realTypes_;
         counter[value]++;
	}

    // Assigning Species Count
	cellSpecies_[cellNo].x = counter[0];
	cellSpecies_[cellNo].y = counter[1];
	cellSpecies_[cellNo].z = counter[2];

	delete[] counter;
}

void Domain::setCell(Cell& cell,int index){
	cells_.at(index) = cell;
}
void Domain::deleteCell(int index){
	cells_.erase(cells_.begin() + index);
}
int Domain::getnCells(void) const{
    return cells_.size();
}
int Domain::getDomainNo(){
	return domainNo_;
}
IntVector3d Domain::findParticleGlobalCoordinate(Vector3d position) {

	IntVector3d cellGlobalCoordinate;
	// PBC Check
    if(PBC_.x==1 || PBC_.y==1 || PBC_.z==1){
    	if(position.x >= domainSizeX_.x)
    		position.x-=domainSizeX_.x*PBC_.x;

    	else if(position.x < 0.0)
    		position.x+=domainSizeX_.x*PBC_.x;

	    if(position.y >= domainSizeY_.y)
	    	position.y-=domainSizeY_.y*PBC_.y;

	    else if(position.y < 0.0)
	    	position.y+=domainSizeY_.y*PBC_.y;

	    if(position.z >= domainSizeZ_.z)
	    	position.z-=domainSizeZ_.z*PBC_.z;

	    else if(position.z < 0.0)
	    	position.z+=domainSizeZ_.z*PBC_.z;

    }

	cellGlobalCoordinate.x = static_cast<int> (floor(position.x/cellSize_.x));
	cellGlobalCoordinate.y = static_cast<int> (floor(position.y/cellSize_.y));
	cellGlobalCoordinate.z = static_cast<int> (floor(position.z/cellSize_.z));



	// Index correction
	if(cellGlobalCoordinate.x >= globalCellArray_.x){
		if (PBC_.x)
			cellGlobalCoordinate.x -= globalCellArray_.x;
		else cellGlobalCoordinate.x = globalCellArray_.x-1;
	}
	if(cellGlobalCoordinate.x < 0){
		if (PBC_.x)
			cellGlobalCoordinate.x += globalCellArray_.x;
		else cellGlobalCoordinate.x = 0;
	}

	if(cellGlobalCoordinate.y >= globalCellArray_.y){
		if (PBC_.y)
			cellGlobalCoordinate.y -= globalCellArray_.y;
		else cellGlobalCoordinate.y = globalCellArray_.y-1;
	}
	if(cellGlobalCoordinate.y < 0){
		if (PBC_.y)
			cellGlobalCoordinate.y += globalCellArray_.y;
		else cellGlobalCoordinate.y = 0;
	}
	if(cellGlobalCoordinate.z >= globalCellArray_.z){
		if (PBC_.z)
			cellGlobalCoordinate.z -= globalCellArray_.z;
		else cellGlobalCoordinate.z = globalCellArray_.z-1;
	}
	if(cellGlobalCoordinate.z < 0){
		if (PBC_.z)
			cellGlobalCoordinate.z += globalCellArray_.z;
		else cellGlobalCoordinate.z = 0;
	}
	return cellGlobalCoordinate;
}

IntVector3d Domain::computeCellLocalCoordinate(IntVector3d cellGlobalcoordinate){

	IntVector3d cellLocalcoordinate;

	cellLocalcoordinate.x  = cellGlobalcoordinate.x - (domainPosition_.x * domainCellArray_.x);
	cellLocalcoordinate.y  = cellGlobalcoordinate.y - (domainPosition_.y * domainCellArray_.y);
	cellLocalcoordinate.z  = cellGlobalcoordinate.z - (domainPosition_.z * domainCellArray_.z);

	return cellLocalcoordinate;
}

Cell & Domain::getCell(int index){
	Cell& refCell = cells_[index];
	return refCell;
}

Cell & Domain::findCellwithLocalCoordinate(IntVector3d cellLocalcoordinate) {

	int cellIndex = (cellLocalcoordinate.x * (domainCellArray_.y*domainCellArray_.z) +
			  cellLocalcoordinate.y * (domainCellArray_.z) +
			  cellLocalcoordinate.z);

	Cell& targetCell = cells_[cellIndex];
	return targetCell;
}

void Domain::constructCells(){


	// Cell construction
	Cell cellObject;
	int cellNo=0,cellCounter=0;
	IntVector3d cellGlobalcoordinate,cellLocalcoordinate;

	for (int i=0; i<domainCellArray_.x; i++){

	  for(int j=0; j<domainCellArray_.y; j++){

		for(int k=0; k<domainCellArray_.z; k++){

			cellNo = cellCounter;

        	cellGlobalcoordinate.x = i + (i * domainPosition_.x);
        	cellGlobalcoordinate.y = j + (j * domainPosition_.y);
        	cellGlobalcoordinate.z = k + (k * domainPosition_.z);

        	cellLocalcoordinate.x  = cellGlobalcoordinate.x - (domainPosition_.x * domainCellArray_.x);
        	cellLocalcoordinate.y  = cellGlobalcoordinate.y - (domainPosition_.y * domainCellArray_.y);
        	cellLocalcoordinate.z  = cellGlobalcoordinate.z - (domainPosition_.z * domainCellArray_.z);

            cellObject.setCellNo(cellNo);
       		cellObject.setGlobalCoordinate(cellGlobalcoordinate);
            cellObject.setLocalCoordinate(cellLocalcoordinate);

            cells_.push_back(cellObject);
            cellCounter++;

		}

	  }

	}
	cout << "MONTECARLO: Cell construction successful " <<endl;
	cout << " No of Cells : " << static_cast<int> (cells_.size()) << endl;

}

void Domain::constructParticles(long nparticles, vector<long> number,
		                  vector<int> type,vector<double> mass,
						  vector<double> position,vector<double> epot, vector<int> accCounter, vector<int> rejCounter){

	 int cellIndex;
	 Particle particle;
	 Vector3d positionHolder;
	 IntVector3d cellGlobalcoordinate, cellLocalcoordinate;

	 for(auto i=0;i<nparticles;i++){

		 particle.setNumber(number.at(i));
		 // assigning unique MC-ID
		 particle.setmcID(i);
		 particle.setType(type.at(i));
		 particle.setMass(mass.at(i));
		 particle.setSuccessFactor(accCounter.at(i));
		 particle.setFailureFactor(rejCounter.at(i));

		 positionHolder.x = position.at((i*3)+0);
		 positionHolder.y = position.at((i*3)+1);
		 positionHolder.z = position.at((i*3)+2);

		 particle.setPosition(positionHolder.x,positionHolder.y,positionHolder.z);
		 particle.setEpot(epot.at(i));

		 // append Particle into list
		 particles_.push_back(particle);

		 cellGlobalcoordinate   = findParticleGlobalCoordinate(positionHolder);

		 cellLocalcoordinate.x  = cellGlobalcoordinate.x - (domainPosition_.x * domainCellArray_.x);
		 cellLocalcoordinate.y  = cellGlobalcoordinate.y - (domainPosition_.y * domainCellArray_.y);
		 cellLocalcoordinate.z  = cellGlobalcoordinate.z - (domainPosition_.z * domainCellArray_.z);

		 cellIndex = (cellLocalcoordinate.x * (domainCellArray_.y*domainCellArray_.z) +
				  cellLocalcoordinate.y * (domainCellArray_.z) +
				  cellLocalcoordinate.z);

		 Cell& cell = getCell(cellIndex);
		 cell.addParticleID(particle.getmcID());

		 // Check if required
		 setCell(cell,cellIndex);
	}
	cout << "MONTECARLO: Particle construction successful " <<endl;

}


void Domain::createSphereNeighbors(){

	Cell cellObject;

	IntVector3d sphereNeighbor;
	IntVector3d cellGlobalcoordinate;
    IntVector3d shiftVector;

    bool boundaryCase;
    bool xAvail,yAvail,zAvail;

	for( decltype(cells_.size()) cellNo=0; cellNo<cells_.size(); cellNo++){

		boundaryCase = false;
		cellObject = getCell(cellNo);
		cellGlobalcoordinate = cellObject.getGlobalCoordinate();

		// Cell Boundary Check
		if( (cellGlobalcoordinate.x==0) || (cellGlobalcoordinate.x==globalCellArray_.x-1) ||
			(cellGlobalcoordinate.y==0) || (cellGlobalcoordinate.y==globalCellArray_.y-1) ||
			(cellGlobalcoordinate.z==0) || (cellGlobalcoordinate.z==globalCellArray_.z-1)){

			boundaryCase = true;
		}

		// Non boundary Cell
		if(!boundaryCase){

			for(int xIndex=-1; xIndex<=+1; xIndex++){
				for(int yIndex=-1; yIndex<=+1; yIndex++){
					for(int zIndex=-1; zIndex<=+1; zIndex++){

					  sphereNeighbor.x = cellGlobalcoordinate.x + xIndex;
					  sphereNeighbor.y = cellGlobalcoordinate.y + yIndex;
					  sphereNeighbor.z = cellGlobalcoordinate.z + zIndex;

					  // shift Vector is always zero
					  shiftVector.x =0; shiftVector.y =0; shiftVector.z =0;

					  cellObject.addSphereNeighbor(sphereNeighbor);
					  cellObject.addShiftIndex(shiftVector);

					}
				}
			}
		}

		// Boundary Cell
		else if(boundaryCase){

			for(int xIndex=-1; xIndex<=+1; xIndex++){
				for(int yIndex=-1; yIndex<=+1; yIndex++){
					for(int zIndex=-1; zIndex<=+1; zIndex++){

					  sphereNeighbor.x = cellGlobalcoordinate.x + xIndex;
					  sphereNeighbor.y = cellGlobalcoordinate.y + yIndex;
					  sphereNeighbor.z = cellGlobalcoordinate.z + zIndex;

					  xAvail = true; yAvail = true; zAvail = true;

					  // Initialize shift Vector
					  shiftVector.x =0; shiftVector.y =0; shiftVector.z =0;

	                  // for Periodic boundary conditions
					  if(sphereNeighbor.x < 0){
						  sphereNeighbor.x = (globalCellArray_.x-1);
		            	  shiftVector.x    = -1;
					      if(!PBC_.x)
					    	  xAvail = false;
					  }

					  if(sphereNeighbor.y < 0){
						  sphereNeighbor.y = (globalCellArray_.y-1);
						  shiftVector.y    = -1;
					      if(!PBC_.y)
					    	  yAvail = false;
					  }

					  if(sphereNeighbor.z < 0){
						  sphereNeighbor.z = (globalCellArray_.z-1);
						  shiftVector.z    = -1;
					      if(!PBC_.z)
					    	  zAvail = false;
					  }

					  if(sphereNeighbor.x > (globalCellArray_.x-1)){
						  sphereNeighbor.x = 0;
		        	      shiftVector.x = +1;
					      if(!PBC_.x)
					    	  xAvail = false;
					  }
					  if(sphereNeighbor.y > (globalCellArray_.y-1)){
						  sphereNeighbor.y = 0;
						  shiftVector.y = +1;
					      if(!PBC_.y)
					    	  yAvail = false;
					  }
					  if(sphereNeighbor.z > (globalCellArray_.z-1)){
						  sphereNeighbor.z = 0;
						  shiftVector.z = +1;
					      if(!PBC_.z)
					    	  zAvail = false;
					  }

                      if( xAvail && yAvail && zAvail){
					     cellObject.addSphereNeighbor(sphereNeighbor);
					     cellObject.addShiftIndex(shiftVector);
                      }

					}
				}
			}
		}
 		setCell(cellObject,cellNo);
	}
}

void Domain::createSampleNeighbors(){

	Cell cellObject;

	IntVector3d sampleNeighbor;
	IntVector3d cellGlobalcoordinate;

    int zoneX,zoneY,zoneZ; // span of sample zone(in multiples of cell)

    zoneX = static_cast<int> (zoneSize_.x/cellSize_.x);
    zoneY = static_cast<int> (zoneSize_.y/cellSize_.y);
    zoneZ = static_cast<int> (zoneSize_.z/cellSize_.z);

    // zone maximum size check
    if(zoneX>(domainCellArray_.x - 2)){
           std::cerr<< "MONTECARLO: sample zone x size is beyond threshold "<< endl;
           exit(1);
    }
    if(zoneY>(domainCellArray_.y - 2)){
           std::cerr<< "MONTECARLO: sample zone y size is beyond threshold "<< endl;
           exit(1);
    }
    if(zoneZ>(domainCellArray_.z - 2)){
           std::cerr<< "MONTECARLO: sample zone z size is beyond threshold "<< endl;
           exit(1);
    }

	int startX,startY,startZ;
	int endX,endY,endZ;

	for( decltype(cells_.size()) cellNo=0; cellNo<cells_.size(); cellNo++){

		cellObject = getCell(cellNo);
		cellGlobalcoordinate = cellObject.getGlobalCoordinate();

		startX = (zoneX/2) + 1 - (zoneX%2);
		startY = (zoneY/2) + 1 - (zoneY%2);
		startZ = (zoneZ/2) + 1 - (zoneZ%2);

		endX = startX; endY = startY; endZ = startZ;

		if(((cellGlobalcoordinate.x == globalCellArray_.x-1) ||
			(cellGlobalcoordinate.x == 0)) && (PBC_.x==0)){
			endX   = startX - (startX * (cellGlobalcoordinate.x/(globalCellArray_.x-1)));
			startX = startX * (cellGlobalcoordinate.x/(globalCellArray_.x-1));

		}
		if(((cellGlobalcoordinate.y == globalCellArray_.y-1) ||
			(cellGlobalcoordinate.y == 0)) && (PBC_.y==0)){
			endY   = startY - (startY * (cellGlobalcoordinate.y/(globalCellArray_.y-1)));
			startY = startY * (cellGlobalcoordinate.y/(globalCellArray_.y-1));
		}
		if(((cellGlobalcoordinate.z == globalCellArray_.z-1) ||
			(cellGlobalcoordinate.z == 0)) && (PBC_.z==0)){
			endZ   = startZ - (startZ * (cellGlobalcoordinate.z/(globalCellArray_.z-1)));
			startZ = startZ * (cellGlobalcoordinate.z/(globalCellArray_.z-1));
		}


		for(int xIndex=(-1*startX); xIndex<=endX; xIndex++){
			for(int yIndex=(-1*startY); yIndex<=endY; yIndex++){
				for(int zIndex=(-1*startZ); zIndex<=endZ; zIndex++){

				  sampleNeighbor.x = cellGlobalcoordinate.x + xIndex;
				  sampleNeighbor.y = cellGlobalcoordinate.y + yIndex;
				  sampleNeighbor.z = cellGlobalcoordinate.z + zIndex;

                  // for Periodic boundary conditions
				  if(sampleNeighbor.x < 0)
					  sampleNeighbor.x = globalCellArray_.x - (1 * -xIndex);
				  if(sampleNeighbor.y < 0)
					  sampleNeighbor.y = globalCellArray_.y - (1 * -yIndex);
				  if(sampleNeighbor.z < 0)
					  sampleNeighbor.z = globalCellArray_.z - (1 * -zIndex);

				  if(sampleNeighbor.x == globalCellArray_.x)
					  sampleNeighbor.x = xIndex-1;
				  if(sampleNeighbor.y == globalCellArray_.y)
					  sampleNeighbor.y = yIndex-1;
				  if(sampleNeighbor.z == globalCellArray_.z)
					  sampleNeighbor.z = zIndex-1;

				  cellObject.addSampleNeighbor(sampleNeighbor);
				}
			}
		}
 		setCell(cellObject,cellNo);

	}

}

long Domain::getnParticlesDomain(){
	return particles_.size();
}

VectorRange Domain::getDomainRange(){
	VectorRange domainRange;
	cout << " Domain::getDomainRange " << endl;
	cout << "domainRange_.xMin "<<domainRange_.xMin<<endl;
    cout << "domainRange_.xMax "<<domainRange_.xMax<<endl;
	cout << "domainRange_.yMin "<<domainRange_.yMin<<endl;
	cout << "domainRange_.yMax "<<domainRange_.yMax<<endl;
	cout << "domainRange_.zMin "<<domainRange_.zMin<<endl;
	cout << "domainRange_.zMax "<<domainRange_.zMax<<endl;

	domainRange = domainRange_;
	return domainRange;
}

void Domain::setRealTypes(int realTypes){
	realTypes_ = realTypes;
}
int Domain::getRealTypes(){
	return realTypes_;
}

void Domain::setTotalTypes(int totalTypes){
	totalTypes_ = totalTypes;
}

int Domain::getTotalTypes(){
	return totalTypes_;
}
void Domain::setPBC(IntVector3d  PBC){
	PBC_ = PBC;
}
void Domain::setSphereLayer(Vector3d  sphereLayerRadii){
	sphereLayerRadii_ = sphereLayerRadii;
}

void Domain::createConfiguration(const vector<Swap>& swapList,vector<Particle>& sphereParticles){
	if (swapList.size() > 2) {
		//TODO Implement and TEST for >2 particles
		cerr << "ERROR: NOT IMPLEMENTED FOR MORE THAN TWO PARTICLES" << endl;
		exit (EXIT_FAILURE);
	}
	if(swapList.size()==1){
		createSphere(swapList[0].ID,sphereParticles);
		return;
	}
	else { // many particles(not currently imppleme)
		const vector<Particle>& refParticles = getParticles();

		const Vector3d& pos1 = refParticles[swapList[0].ID].getPosition();
		const Vector3d& pos2 = refParticles[swapList[1].ID].getPosition();

		if (!doSpheresOverlap(pos1, pos2)) {
			createSphere(swapList[0].ID, sphereParticles);
			createSphere(swapList[1].ID, sphereParticles);
			return;
		}

		else {
			vector<Particle> tmpList1;
			vector<Particle> tmpList2;
			createSphere(swapList[0].ID, tmpList1);
			createSphere(swapList[1].ID, tmpList2);

			//sorting two sphere List
			std::sort(tmpList1.begin(), tmpList1.end(),
					[](const Particle &left,const Particle &right) {
						return left.getmcID() < right.getmcID();
					});
			std::sort(tmpList2.begin(), tmpList2.end(),
					[](const Particle &left,const Particle &right) {
						return left.getmcID() < right.getmcID();
					});

			size_t index1 = 0, index2 = 0;
			// find the union of two sphere list
			while ((index1 < tmpList1.size()) && (index2 < tmpList2.size())) {
				if (tmpList1[index1].getmcID() == tmpList2[index2].getmcID()) {
					tmpList1[index1].setType(min(tmpList1[index1].getType(), tmpList2[index2].getType()));
					if(tmpList1[index1].getmcID() != swapList[0].ID && tmpList1[index1].getmcID() != swapList[1].ID)
					    sphereParticles.push_back(tmpList1[index1]);
					index1++;
					index2++;
				} else if (tmpList1[index1].getmcID()< tmpList2[index2].getmcID()) {
					if(tmpList1[index1].getmcID() != swapList[0].ID && tmpList1[index1].getmcID() != swapList[1].ID)
					    sphereParticles.push_back(tmpList1[index1]);
					index1++;
				} else {
					if(tmpList2[index2].getmcID() != swapList[0].ID && tmpList2[index2].getmcID() != swapList[1].ID)
						sphereParticles.push_back(tmpList2[index2]);
					index2++;
				}
			}

			while ((index1 < tmpList1.size())) {
				sphereParticles.push_back(tmpList1[index1]);
				index1++;
			}
			while ((index2 < tmpList2.size())) {
				sphereParticles.push_back(tmpList2[index2]);
				index2++;
			}
		}
	}// end of many particles part

}


void  Domain::createSphere(long particleNo, vector<Particle>& sphereParticles){

	Cell targetCell;
	Cell nbCell;

	// sphere construction routine
	Particle    sphereParticle;
	IntVector3d cellGlobalCoordinate,nbGlobalCoordinate;
	IntVector3d cellLocalCoordinate, nbLocalCoordinate;

	double refZone1Distance,refZone2Distance, distanceSquare;
	double refZone1DistanceSquare, refZone2DistanceSquare, refSphereDistanceSquare;

    Vector3d tempPosition;
    //bool boundaryIssue=false;
    int xFactor=0, yFactor=0, zFactor=0;

	const vector<Particle>& refParticles = getParticles();
	Vector3d position = refParticles[particleNo].getPosition();

	cellGlobalCoordinate   = findParticleGlobalCoordinate(position);
	cellLocalCoordinate    = computeCellLocalCoordinate(cellGlobalCoordinate);
	targetCell             = findCellwithLocalCoordinate(cellLocalCoordinate);

	double sphereRadii    = sphereLayerRadii_.x;
	double zone1Radii     = sphereLayerRadii_.y;
	double zone2Radii     = sphereLayerRadii_.z;

	refZone2Distance       = ( sphereRadii + zone1Radii + zone2Radii );
	refZone2DistanceSquare = refZone2Distance * refZone2Distance;
	refZone1Distance        =  sphereRadii + zone1Radii;
    refZone1DistanceSquare  =  refZone1Distance * refZone1Distance;
    refSphereDistanceSquare = sphereRadii * sphereRadii;

    for(decltype(targetCell.getSphereNBLSize()) i=0; i<targetCell.getSphereNBLSize(); i++){

    		xFactor=0; yFactor=0; zFactor=0; // Reinitialize Factors
            nbGlobalCoordinate = targetCell.getSphereNeighbor(i);

            // for MPI communication make respective calls here
            // find CPU no from Global Coordinate
            nbLocalCoordinate = computeCellLocalCoordinate(nbGlobalCoordinate);
            nbCell            = findCellwithLocalCoordinate(nbLocalCoordinate);

            IntVector3d shiftVector = targetCell.getShiftIndex(i); // Wrap around corrections if any
            xFactor = shiftVector.x; yFactor = shiftVector.y; zFactor = shiftVector.z;

            for(auto j=0; j< nbCell.getnParticlesCell(); j++){

            	long particleID     = nbCell.getParticleID(j);

                // access particle from particles_
            	sphereParticle = particles_[particleID];
            	tempPosition   = sphereParticle.getPosition();

            	// Wrap around correction
            	tempPosition.x = tempPosition.x + (xFactor * totalSizeX_.x);
            	tempPosition.y = tempPosition.y + (yFactor * totalSizeY_.y);
            	tempPosition.z = tempPosition.z + (zFactor * totalSizeZ_.z);

            	distanceSquare = computeDistanceSquare(position,tempPosition);

            	// ZONE0 - (CORE SPHERE)
            	if ((distanceSquare <=refSphereDistanceSquare) && (sphereParticle.getmcID() != particleNo)){
                    sphereParticles.push_back(sphereParticle);
            	}
                // Zone1 (MIDDLE LAYER)
            	else if( (distanceSquare <= refZone1DistanceSquare) && (distanceSquare>refSphereDistanceSquare) ){
                    sphereParticle.setType(sphereParticle.getType()+totalTypes_);
                    sphereParticles.push_back(sphereParticle);
            	}
            	// Zone2 (OUTERMOST LAYER)
            	else if((distanceSquare <= refZone2DistanceSquare) && (distanceSquare>refZone1DistanceSquare)){
                    if(sphereParticle.getType()%realTypes_!=2){ // ignoring placeholder particles
            		   sphereParticle.setType(sphereParticle.getType()+(totalTypes_*2));
            		   sphereParticles.push_back(sphereParticle);
                    }
            	}

            }// loop over particles
    }// loop over cells
}

void  Domain::createSphereTarget(long particleNo, vector<Particle>& sphereParticles){

	Cell targetCell;
	Cell nbCell;

	// sphere construction routine
	Particle    sphereParticle;
	IntVector3d cellGlobalCoordinate,nbGlobalCoordinate;
	IntVector3d cellLocalCoordinate, nbLocalCoordinate;

//	double refZone1Distance,refZone2Distance, distanceSquare;
//	double refZone1DistanceSquare, refZone2DistanceSquare, refSphereDistanceSquare;
	double refSphereDistanceSquare, distanceSquare;

    Vector3d tempPosition;
    //bool boundaryIssue=false;
    int xFactor=0, yFactor=0, zFactor=0;

	const vector<Particle>& refParticles = getParticles();
	Vector3d position = refParticles[particleNo].getPosition();

	cellGlobalCoordinate   = findParticleGlobalCoordinate(position);
	cellLocalCoordinate    = computeCellLocalCoordinate(cellGlobalCoordinate);
	targetCell             = findCellwithLocalCoordinate(cellLocalCoordinate);

	double sphereRadii    = getSphereRadiusTarget();
//	double zone1Radii     = sphereLayerRadii_.y;
//	double zone2Radii     = sphereLayerRadii_.z;

//	refZone2Distance       = ( sphereRadii + zone1Radii + zone2Radii );
//	refZone2DistanceSquare = refZone2Distance * refZone2Distance;
//	refZone1Distance        =  sphereRadii + zone1Radii;
//    refZone1DistanceSquare  =  refZone1Distance * refZone1Distance;
    refSphereDistanceSquare = sphereRadii * sphereRadii;

    for(decltype(targetCell.getSphereNBLSize()) i=0; i<targetCell.getSphereNBLSize(); i++){

    		xFactor=0; yFactor=0; zFactor=0; // Reinitialize Factors
            nbGlobalCoordinate = targetCell.getSphereNeighbor(i);

            // for MPI communication make respective calls here
            // find CPU no from Global Coordinate
            nbLocalCoordinate = computeCellLocalCoordinate(nbGlobalCoordinate);
            nbCell            = findCellwithLocalCoordinate(nbLocalCoordinate);

            IntVector3d shiftVector = targetCell.getShiftIndex(i); // Wrap around corrections if any
            xFactor = shiftVector.x; yFactor = shiftVector.y; zFactor = shiftVector.z;

            for(auto j=0; j< nbCell.getnParticlesCell(); j++){

            	long particleID     = nbCell.getParticleID(j);

                // access particle from particles_
            	sphereParticle = particles_[particleID];
            	tempPosition   = sphereParticle.getPosition();

            	// Wrap around correction
            	tempPosition.x = tempPosition.x + (xFactor * totalSizeX_.x);
            	tempPosition.y = tempPosition.y + (yFactor * totalSizeY_.y);
            	tempPosition.z = tempPosition.z + (zFactor * totalSizeZ_.z);

            	distanceSquare = computeDistanceSquare(position,tempPosition);

            	// ZONE0 - (CORE SPHERE)
            	if ((distanceSquare <=refSphereDistanceSquare) && (sphereParticle.getmcID() != particleNo)){
                    sphereParticles.push_back(sphereParticle);
            	}
                // Zone1 (MIDDLE LAYER)
//            	else if( (distanceSquare <= refZone1DistanceSquare) && (distanceSquare>refSphereDistanceSquare) ){
//                    sphereParticle.setType(sphereParticle.getType()+totalTypes_);
//                    sphereParticles.push_back(sphereParticle);
//            	}
            	// Zone2 (OUTERMOST LAYER)
//            	else if((distanceSquare <= refZone2DistanceSquare) && (distanceSquare>refZone1DistanceSquare)){
//                    if(sphereParticle.getType()%realTypes_!=2){ // ignoring placeholder particles
//            		   sphereParticle.setType(sphereParticle.getType()+(totalTypes_*2));
//            		   sphereParticles.push_back(sphereParticle);
//                    }
//            	}

            }// loop over particles
    }// loop over cells
}


void Domain::shiftSphere(vector<Particle>& sphereParticles,Vector3d position){

	double     refZone1Distance;
	Vector3d    tempPosition;

	double sphereRadii = sphereLayerRadii_.x;
	double zone1Radii  = sphereLayerRadii_.y;
	double zone2Radii  = sphereLayerRadii_.z;

	refZone1Distance       = (sphereRadii+zone1Radii+zone2Radii);

    // minimum boundary - Particle shift value
    double refXmin = position.x - (refZone1Distance);
    double refYmin = position.y - (refZone1Distance);
    double refZmin = position.z - (refZone1Distance);

    for (auto& sphereParticle : sphereParticles){
    		tempPosition = sphereParticle.getPosition();

    		// shift coordinates back to native coordinate system
            tempPosition.x +=  refXmin;
            tempPosition.y += refYmin;
            tempPosition.z += refZmin;

            // PBC effect correction
            if(tempPosition.x>domainSizeX_.x)
            	tempPosition.x = tempPosition.x - (domainSizeX_.x*PBC_.x);
            if(tempPosition.y>domainSizeY_.y)
            	tempPosition.y = tempPosition.y - (domainSizeY_.y*PBC_.y);
            if(tempPosition.z>domainSizeZ_.z)
            	tempPosition.z = tempPosition.z - (domainSizeZ_.z*PBC_.z);
            if(tempPosition.x<0.0)
            	tempPosition.x = tempPosition.x + (domainSizeX_.x*PBC_.x);
            if(tempPosition.y<0.0)
            	tempPosition.y = tempPosition.y + (domainSizeY_.y*PBC_.y);
            if(tempPosition.z<0.0)
            	tempPosition.z = tempPosition.z + (domainSizeZ_.z*PBC_.z);

            sphereParticle.setPosition(tempPosition.x,tempPosition.y,tempPosition.z);

    }
}


void Domain::updateDomain(vector<Particle>& config , vector<Swap>& swapList){

	 for (const auto& updateObj : config) {
		long  mcID = updateObj.getmcID();

	    // check for sphere core/ zone1
	    if(updateObj.getType() < totalTypes_){

	    	// Get Particle older information from Domain data base
			Vector3d oldPosition                = particles_[mcID].getPosition();
			IntVector3d oldCellGlobalcoordinate = findParticleGlobalCoordinate(oldPosition);
			IntVector3d oldCellLocalcoordinate;
		    oldCellLocalcoordinate.x            = oldCellGlobalcoordinate.x - (domainPosition_.x * domainCellArray_.x);
		    oldCellLocalcoordinate.y            = oldCellGlobalcoordinate.y - (domainPosition_.y * domainCellArray_.y);
		    oldCellLocalcoordinate.z            = oldCellGlobalcoordinate.z - (domainPosition_.z * domainCellArray_.z);
		    int oldCellNo                       = findCellwithLocalCoordinate(oldCellLocalcoordinate).getCellNo();

		    // Get Particle current information from Domain data base
			Vector3d newPosition                = updateObj.getPosition();
			IntVector3d newCellGlobalcoordinate = findParticleGlobalCoordinate(newPosition);
			IntVector3d newCellLocalcoordinate;
		    newCellLocalcoordinate.x            = newCellGlobalcoordinate.x - (domainPosition_.x * domainCellArray_.x);
		    newCellLocalcoordinate.y            = newCellGlobalcoordinate.y - (domainPosition_.y * domainCellArray_.y);
		    newCellLocalcoordinate.z            = newCellGlobalcoordinate.z - (domainPosition_.z * domainCellArray_.z);
		    int newCellNo                       = findCellwithLocalCoordinate(newCellLocalcoordinate).getCellNo();

            // Update type in Domain for the swap list
		    for(decltype(swapList.size()) i=0; i<swapList.size();i++){
                long ID = swapList[i].ID;

	            if(mcID == ID){
	         	   particles_[mcID].setType(updateObj.getType());
	         	   // Update targetSites list
	         	   if(updateObj.getType()==TARGET){
	         		   addTargetSite(mcID);
	         		   deleteSampleSite(mcID);
	         	   }
	         	   else if(updateObj.getType()==SAMPLE){
	         		   deleteTargetSite(mcID);
	         		   addSampleSite(mcID);
	         	   }
	            }
		    }

            // Particle exchange between Cells
		    if(oldCellNo != newCellNo){
		    	// deleting id from Old cell
		    	Cell& oldCell                    = getCell(oldCellNo);
		    	oldCell.deleteParticleID(updateObj.getmcID());

		    	// adding id to New Cell
		    	Cell& newCell                    = getCell(newCellNo);
		    	newCell.addParticleID(updateObj.getmcID());

		    }
			particles_[mcID].setPosition(newPosition.x,newPosition.y,newPosition.z); // updating position
			particles_[mcID].setEpot(updateObj.getEpot()); // Updating Epot
	    }
	    // redefine in terms of distance!!
	    else if((updateObj.getType()>=totalTypes_) && (updateObj.getType()<2*totalTypes_)){// else zone 2
	    	particles_[mcID].setEpot(updateObj.getEpot()); // Update Epot
	    }

	 }// loop over config
}

void Domain::filterSamplingSites(void){

    // loop over particles
    for(const auto& p : getParticles()){
    	int typeHold = p.getType();

    	// filter Sample Type
    	if(typeHold == 2){
    		samplingSites_.push_back(p.getmcID());
    	}
    	// filter Target type
    	else if(typeHold == 1){
    		targetSites_.push_back(p.getmcID());
    	}
    }

    cout << " Initial Configuration Sample Sites : " << samplingSites_.size()<< endl;
    cout << " Initial Configuration Target Sites : " << targetSites_.size()<< endl;

}


long   Domain::getRandomId(int targetType){
	auto size = targetType==1?targetSites_.size()-1 : samplingSites_.size()-1;
    std::uniform_int_distribution<int> distribution(0,size);
	if (targetType == 1)
		return targetSites_[distribution(generator_)];
	else
		return samplingSites_[distribution(generator_)];
}
bool Domain::doEnergyZoneOverlap(Vector3d v1,Vector3d v2){
	double radius = getSphereRadius() + getWallThickness();
	double distSquare = getMinimumImageSquaredDistance(v1,v2);
	return (distSquare <= 4*radius*radius);
}

bool Domain::doSpheresOverlap(Vector3d v1,Vector3d v2){
	double radius = getSphereRadius() + 2*getWallThickness();
	double distSquare = getMinimumImageSquaredDistance(v1,v2);
	return (distSquare <= 4*radius*radius); // Check Distances for Spatial conflicting cases
}

void  Domain::writeSphereConfiguration(string filename,vector<Particle>& sphereParticles){

    // new file writer part!!
    ofstream fout(filename, ios_base::out);

    //  Sphere Core radii + SphereWall Thickness + SphereWall Thickness;
    //double boxSize = 2*(sphereLayerRadii_.x + sphereLayerRadii_.y + sphereLayerRadii_.z);

    const vector<Particle>& refParticles = getParticles();

    // print IMD header -- check
    fout <<"#F A 1 1 1 3 0 1"<< endl;
    fout <<"#C number type mass x y z Epot"<< endl;
    fout <<"#X "<<totalSizeX_.x<<" "<<totalSizeX_.y<<" "<<totalSizeX_.z<<" "<< endl;
    fout <<"#Y "<<totalSizeY_.x<<" "<<totalSizeY_.y<<" "<<totalSizeY_.z<<" "<< endl;
    fout <<"#Z "<<totalSizeZ_.x<<" "<<totalSizeZ_.y<<" "<<totalSizeZ_.z<<" "<< endl;
    fout <<"#E "<< endl;


    // loop over particles in the container
    for(const auto& p : sphereParticles){
       fout << refParticles[p.getmcID()].getNumber()
       << "  " << p.getType()
       << "  " << setprecision(6) << p.getMass()
       << "  " << setprecision(6) << p.getPosition().x
       << "  " << setprecision(6) << p.getPosition().y
       << "  " << setprecision(6) << p.getPosition().z
       << "  " << setprecision(6) << p.getEpot()
       << endl;

    }

    fout.close(); // closing outfile connection

}

void  Domain::writeConfiguration(string filename){
                       
                   cout <<" WRITE CONFIG CHECK"<<endl;

		    // new file writer part!!
		    ofstream fout(filename, ios_base::out);

		    // print IMD header -- check
		    fout <<"#F A 1 1 1 3 0 2"<< endl;
		    fout <<"#C number type mass x y z Epot"<< endl;
		    fout <<"#X "<<totalSizeX_.x<<" "<<totalSizeX_.y<<" "<<totalSizeX_.z<<" "<< endl;
		    fout <<"#Y "<<totalSizeY_.x<<" "<<totalSizeY_.y<<" "<<totalSizeY_.z<<" "<< endl;
		    fout <<"#Z "<<totalSizeZ_.x<<" "<<totalSizeZ_.y<<" "<<totalSizeZ_.z<<" "<< endl;
		    fout <<"##PBC "<<PBC_.x<<" "<<PBC_.y<<" "<<PBC_.z<<" "<< endl;
		    fout <<"#E "<< endl;

		    // loop over particles in Domain
		    for(const auto& p : getParticles()){

                    //Quick hack for filtering Carbon and boundary Fe atoms
                    // NOTE: NOT TO BE INCLUDED IN STANDARD VERSION

                      if(p.getType()==0 || p.getType()==1 || p.getType()==3){
                           //Marking sampled sites
		       fout << p.getNumber()
		       << "  " << p.getType()
		       << "  " << setprecision(6) << p.getMass()
		       << "  " << setprecision(6) << p.getPosition().x
		       << "  " << setprecision(6) << p.getPosition().y
		       << "  " << setprecision(6) << p.getPosition().z
                       << "  " << setprecision(6) << p.getEpot()
		       << endl;
		      }
 
                    }   

		    fout.close(); // closing outfile connection

}

void  Domain::writeCompleteConfiguration(string filename){
                       
                   cout <<" WRITE CONFIG CHECK"<<endl;

		    // new file writer part!!
		    ofstream fout(filename, ios_base::out);



		    // print IMD header -- check
		    fout <<"#F A 1 1 1 3 0 2"<< endl;
		    fout <<"#C number type mass x y z Epot AccCounter RejCounter"<< endl;
		    fout <<"#X "<<totalSizeX_.x<<" "<<totalSizeX_.y<<" "<<totalSizeX_.z<<" "<< endl;
		    fout <<"#Y "<<totalSizeY_.x<<" "<<totalSizeY_.y<<" "<<totalSizeY_.z<<" "<< endl;
		    fout <<"#Z "<<totalSizeZ_.x<<" "<<totalSizeZ_.y<<" "<<totalSizeZ_.z<<" "<< endl;
		    fout <<"##PBC "<<PBC_.x<<" "<<PBC_.y<<" "<<PBC_.z<<" "<< endl;
		    fout <<"#E "<< endl;


		    // loop over particles in Domain
		    for(const auto& p : getParticles()){

		       fout << p.getNumber()
		       << "  " << p.getType()
		       << "  " << setprecision(6) << p.getMass()
		       << "  " << setprecision(6) << p.getPosition().x
		       << "  " << setprecision(6) << p.getPosition().y
		       << "  " << setprecision(6) << p.getPosition().z
		       << "  " << setprecision(6) << p.getEpot()
		       << "  " << p.getSuccessFactor()
		       << "  " << p.getFailureFactor()
		       << endl;
		  
                      }   

		    fout.close(); // closing outfile connection

}

void Domain::writeParameterFile(string filename, int trialType, long conflicts,long finished,long swaps,long accSwaps,long rejSwaps,double totDeltaEpot,long biasedMovesPerformed,long localMovesPerformed){

	string inputFile = "Simu_MC_output_"+to_string(finished)+".chkpt";

	ofstream fout(filename, ios_base::out);

	fout <<"#param file for Monte Carlo simulation" << endl;
	fout <<"#The units are angstrom, eV and amu for the mass " << endl;
	fout << "mc_Cpudim"     << "		"<< cpuDimension_.x << " " << cpuDimension_.y << " " << cpuDimension_.z << endl;
	fout << "mc_Steps"      << "		"<< "MCSTEPS" << endl;
	fout << "mc_Seed"       << "		"<< userSeed_ << endl;
	fout << "mc_totalTypes" << "		"<< totalTypes_ << endl;
	fout << "mc_realTypes"     << "		"<< realTypes_ << endl;

	fout << " " << endl;

	fout <<"mc_SimboxX "       << "		"<<totalSizeX_.x<<" "<<totalSizeX_.y<<" "<<totalSizeX_.z<<" "<< endl;
	fout <<"mc_SimboxY "       << "		"<<totalSizeY_.x<<" "<<totalSizeY_.y<<" "<<totalSizeY_.z<<" "<< endl;
	fout <<"mc_SimboxZ "       << "		"<<totalSizeZ_.x<<" "<<totalSizeZ_.y<<" "<<totalSizeZ_.z<<" "<< endl;
	fout <<"mc_zoneSize"       << "		"<< zoneSize_.x << " " << zoneSize_.y <<" " << zoneSize_.z << " " << endl;
	fout <<"mc_pbc "           << "		"<<PBC_.x<<" "<<PBC_.y<<" "<<PBC_.z<<" "<< endl;
	fout <<"mc_temperature"    << "		"<< temperature_ << endl;
	fout <<"mc_sphereRadius"   << "		"<< sphereLayerRadii_.x << endl;
	fout <<"mc_sphereWall"     << "		"<< sphereLayerRadii_.y << endl;
	fout <<"mc_cohesiveEnergy" << "		"<< cohesiveEnergy_ << endl;
	fout <<"mc_flushInterval"  << "		"<< writeMCinterval_ << endl;
	fout <<"mc_MasterFlag"     << "		"<< masterFlag_ << endl;
	fout <<"mc_inputFile"      << "		"<< inputFile << endl;
	fout <<"mc_outputFile"     << "		"<< "Simu_MC_output_MCSTEPS.chkpt" << endl;
	fout <<"mc_statFile"       << "		"<< "simStat_new_MCSTEPS.txt" << endl;
	fout <<"mc_sphereParam"    << "		"<< mcSphereFile_ << endl;

	fout << " " << endl;

	fout << "# 1- RANDOM SINGLE TRIAL MOVE, 2- SWAP TRIAL MOVE 3 CLUSTER TRIAL MOVE" << endl;
	fout << "mc_trialMoveType" << "		"<< trialType << endl;

	fout << " " << endl;

	fout << "# 1 - start New ,Otherwise start from intermediate configuration" << endl;              
	fout <<"mc_simulationID"   << "		"<< "0" << endl; //BY DEFAULT

	fout << " " << endl;

	fout << "# Statistical counters for SWAP type trial move (not required if mc_simulationID=1 or mc_trialMoveType=1)" << endl;
	fout << "mc_nConflicts"           << "		"<< conflicts << endl;
	fout << "mc_nFinished"            << "		"<< finished << endl;
	fout << "mc_nSwaps"               << "		"<< swaps << endl;
	fout << "mc_nAccepSwaps"          << "		"<< accSwaps << endl;
	fout << "mc_nRejeSwaps"           << "		"<< rejSwaps << endl;
	fout << "mc_totDeltaEpot"         << "		"<< totDeltaEpot << endl;
	fout << "mc_biasedMovesPerformed" << "		"<< biasedMovesPerformed << endl;
	fout << "mc_localMovesPerformed"  << "		"<< localMovesPerformed << endl;

	fout << " " << endl;

	fout << "# Biased sampling " << endl;

	fout <<"mc_samplingMode"       << "		"<< samplingMode_ << endl;
	fout <<"mc_stepsInBiasedBlock" << "		"<< stepsInBiasedBlock_ << endl;
	fout <<"mc_stepsInLocalBlock"  << "		"<< stepsInLocalBlock_ << endl;
	fout <<"mc_sphereRadiusTarget" << "		"<< sphereRadiusTarget_ << endl;
	fout <<"mc_localMovesRadius"   << "		"<< localMovesRadius_ << endl;
	fout <<"mc_cylinderRadius"   << "		"<< cylinderRadius_ << "	" <<"#Just used in mono -> mc_samplingMode==8"<< endl;
	fout <<"mc_coveringTimes"      << "		"<< coveringTimes_ << endl;

	fout.close(); 

}

double Domain::getMinimumImageSquaredDistance(Vector3d refPos,Vector3d currPos){
	 Vector3d dir;
	 dir.x = refPos.x - currPos.x;
	 dir.y = refPos.y - currPos.y;
	 dir.z = refPos.z - currPos.z;

	 if (PBC_.x){
		if (dir.x > totalSizeX_.x * 0.5) dir.x -= totalSizeX_.x;
		else if (dir.x < -totalSizeX_.x * 0.5) dir.x += totalSizeX_.x;
	}
	if (PBC_.y){
		if (dir.y > totalSizeY_.y * 0.5) dir.y -= totalSizeY_.y;
		else if (dir.y < -totalSizeY_.y * 0.5) dir.y += totalSizeY_.y;
	}
	if (PBC_.z){
		if (dir.z > totalSizeZ_.z * 0.5f) dir.z -= totalSizeZ_.z;
		else if (dir.z < -totalSizeZ_.z * 0.5f) dir.z += totalSizeZ_.z;
	}
	return dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
}

void Domain::updateSampleFrequency(vector<Swap>& swapList){
	for(decltype(swapList.size()) i=0; i<swapList.size();i++){
	    long mcID = swapList[i].ID;
		particles_[mcID].updateFrequency(1);
	}
}

double  Domain::getTemperature(void){
      return temperature_;
}
double  Domain::getCohesiveEnergy(void){
      return cohesiveEnergy_;
}
double  Domain::getSphereRadius(void){
	// sphereLayerRadii_.x contains the core sphere radius
	return sphereLayerRadii_.x;
}
double  Domain::getWallThickness(void){
	// sphereLayerRadii_.y contains the core sphere Wall thickness
	return sphereLayerRadii_.y;
}
Vector3d  Domain::getBoxSizeX(){
	return totalSizeX_;
}
Vector3d  Domain::getBoxSizeY(){
	return totalSizeY_;
}
Vector3d  Domain::getBoxSizeZ(){
	return totalSizeZ_;
}
long   Domain::getMCWriteInterval(){
     	return writeMCinterval_;
}
long   Domain::getMCTotalSteps(){
        return totalMCSteps_;
}

void  Domain::addSampleSite(long mcID){
       samplingSites_.push_back(mcID);
}

void  Domain::deleteSampleSite(long mcID){
	long ID;

	for( decltype(samplingSites_.size()) i=0; i < samplingSites_.size(); i++){
		ID = samplingSites_[i];

		if(ID == mcID){
		      // assign last element from mcIDList_ at index to be deleted
			  samplingSites_[i] = samplingSites_.back();
		      // delete last element
			  samplingSites_.pop_back();
		      return;
		}
	}
	cerr<<"ERROR: SAMPLE PARTICLE WITH ID "<< mcID<<" NOT FOUND FOR DELETION"<<endl;
	exit(1);
}

void Domain::addTargetSite(long mcID){
     	targetSites_.push_back(mcID);
}
void Domain::deleteTargetSite(long mcID){
	long ID;

	for( decltype(targetSites_.size()) i=0; i < targetSites_.size(); i++){
		ID = targetSites_[i];

		if(ID == mcID){
		      // assign last element from mcIDList_ at index to be deleted
			  targetSites_[i] = targetSites_.back();
		      // delete last element
			  targetSites_.pop_back();
		      return;
		}
	}
	cerr<<"ERROR: TARGET PARTICLE WITH ID "<< mcID<<" NOT FOUND FOR DELETION"<<endl;
	exit(1);
}

int Domain::getUserSeed(){
	return userSeed_;
}

long Domain::getTargetListSize(void) const{
	return targetSites_.size();
}

long Domain::getSampleListSize(void) const{
	return samplingSites_.size();
}
double Domain::getTotalEnergy(void) const{
	int newType;
	KahanSum energyGlobal;

    for(const auto& newParticle : particles_){
        newType     = newParticle.getType();
        // ignoring placeholders and boundary
        if( (newType < realTypes_) && (newType%realTypes_ != 2) ){
    		energyGlobal.add(newParticle.getEpot());
    	}
    }
    return energyGlobal.getSum();
}

void Domain::increaseSuccessOrFailureFactor(vector<Particle>& Config, bool Decision){	
	
	int configID;

	if(Decision){
		for(auto& particle : Config){
			if( particle.getType() == 1 || particle.getType() == 2 ){		
				configID = particle.getmcID();
				particles_[configID].increaseSuccessFactor();
			}
		}		
	}
	else{
		for(auto& particle : Config){
			if( particle.getType() == 1 || particle.getType() == 2 ){
				configID = particle.getmcID();
				particles_[configID].increaseFailureFactor();
			}
		}
	}

}

long Domain::getBiasedId(int targetType){		
	
	if (targetType == 1){
		auto size = targetSites_.size()-1; 
		std::uniform_int_distribution<int> distribution(0,size);
		return targetSites_[distribution(generator_)];
	}

	else{
		long randomBiasedID = biased_distribution();	
		return randomBiasedID;			
	}
}

long Domain::biased_distribution(void){

	if ( getSamplingMode()==5 || getSamplingMode()==6 ){

		const vector<Particle>& globalList = getParticles();
		vector<Particle> configHolder = globalList;	
		int bellCenter = getBellCenter();
		int samplingSitesSize = samplingSites_.size();
		double *samplingSitesProbability = new double[samplingSitesSize];	
		double *shuffledVacanciesIDs = new double[samplingSitesSize];
		double totalProbability = 0;
		int index = 0;
		int typeHold = 0;	

		for ( const auto& p : configHolder ){
			typeHold = p.getType();
			if (typeHold == 2){	
				totalProbability += exp( -(p.getSuccessFactor() - bellCenter )*(p.getSuccessFactor() - bellCenter ) );	
				samplingSitesProbability[index] = exp( -(p.getSuccessFactor() - bellCenter )*(p.getSuccessFactor() - bellCenter ) );	
				shuffledVacanciesIDs[index] = p.getmcID();
				index++;
			}
		}

		long biasedIndex = 0;					
		double sum0 = 0;					
		double sum1 = sum0 + samplingSitesProbability[0];

		std::uniform_real_distribution<> dis(0, totalProbability);		
		double randomNumber = dis(generator_);		
		long biasedID;	

		while(!( sum0 <= randomNumber && randomNumber <= sum1 )){			
			biasedIndex++;
			sum0 = sum1;
			sum1 = sum1 + samplingSitesProbability[biasedIndex];	
		}

		biasedID = shuffledVacanciesIDs[biasedIndex];

		delete[] shuffledVacanciesIDs;
	    	delete[] samplingSitesProbability;
		return biasedID;

	}

	else{

		const vector<Particle>& globalList = getParticles();
		vector<Particle> configHolder = globalList;
		int typeHold = 0;	
	 
		double Kappa = getKappa();
		double KappaR = getKappaR();

		int samplingSitesSize = samplingSites_.size();
		double *samplingSitesProbability = new double[samplingSitesSize];	
		double *shuffledVacanciesIDs = new double[samplingSitesSize];
		double totalProbability = 0;
		int index = 0;

		for ( const auto& p : configHolder ){
			typeHold = p.getType();
			if (typeHold == 2){	
				totalProbability += exp( Kappa*  p.getSuccessFactor() - KappaR* p.getFailureFactor() );	
				samplingSitesProbability[index] = exp( Kappa *  p.getSuccessFactor() - KappaR* p.getFailureFactor()  );	
				shuffledVacanciesIDs[index] = p.getmcID();
				index++;
			}
		}

		long biasedIndex = 0;					
		double sum0 = 0;					
		double sum1 = sum0 + samplingSitesProbability[0];

		std::uniform_real_distribution<> dis(0, totalProbability);		
		double randomNumber = dis(generator_);		
		long biasedID;	

		while(!( sum0 <= randomNumber && randomNumber <= sum1 )){			
			biasedIndex++;
			sum0 = sum1;
			sum1 = sum1 + samplingSitesProbability[biasedIndex];	
		}

		biasedID = shuffledVacanciesIDs[biasedIndex];

		delete[] shuffledVacanciesIDs;
	    	delete[] samplingSitesProbability;
		return biasedID;

	}

}


void Domain::calculateEquilibrationStep(void){

	int totalVacancies = samplingSites_.size();
	double sphereVolume = 4.0/3 * 3.1416 *sphereRadiusTarget_*sphereRadiusTarget_*sphereRadiusTarget_;
	double cellVolume = 3*3*3;
	int cellsInSphere = int (sphereVolume / cellVolume );
	int vacanciesInSphere = cellsInSphere*6;
	int trialsNumber = totalVacancies / vacanciesInSphere * coveringTimes_  ;
	equilibrationStep_ = trialsNumber; 

}

long Domain::getLocalId(long targetID){

	double maxRadiusDistance = getLocalMovesRadius(); 
	double dx, dy, dz, distanceToSample;
	long sampleID;
	const vector<Particle>& refParticles = getParticles();
	const Vector3d& targetPos = refParticles[targetID].getPosition();

	for ( decltype(samplingSites_.size()) i=0; i < samplingSites_.size(); i++){
		sampleID = getRandomId(SAMPLE);
		const Vector3d& samplePos = refParticles[sampleID].getPosition();
		dx = samplePos.x - targetPos.x;	
		dy = samplePos.y - targetPos.y;
		dz = samplePos.z - targetPos.z;
		distanceToSample = sqrt ( dx*dx + dy*dy + dz*dz); 
		if ( distanceToSample <= maxRadiusDistance ){	
			return sampleID;					
		}
	}
	
	return 0; 
}

long Domain::getMonoId(long targetID){

	double dxSample, dySample, dxTarget, dyTarget;
	double distanceToTheCenterSample, distanceToTheCenterTarget;
	long sampleID;
	const vector<Particle>& refParticles = getParticles();
	double centerCoords = getCylinderRadius();

	const Vector3d& targetPos = refParticles[targetID].getPosition();
	dxTarget = targetPos.x - centerCoords;	
	dyTarget = targetPos.y - centerCoords;
	distanceToTheCenterTarget = sqrt ( dxTarget*dxTarget + dyTarget*dyTarget); 

	for ( decltype(samplingSites_.size()) i=0; i < samplingSites_.size(); i++){
		sampleID = getRandomId(SAMPLE);
		const Vector3d& samplePos = refParticles[sampleID].getPosition();
		dxSample = samplePos.x - centerCoords;	
		dySample = samplePos.y - centerCoords;
		distanceToTheCenterSample = sqrt ( dxSample*dxSample + dySample*dySample); 

		if ( distanceToTheCenterSample <= distanceToTheCenterTarget ){
			return sampleID;					
		}
	}
	
	return 0;
}


