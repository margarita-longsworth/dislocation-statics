/*
 * Cell.c++
 *
 */
#include "Cell.h"
#include <iostream>
using namespace std;

Cell::Cell(void){
	      cellNo_                = 0;
	      globalCoordinate_      = {0,0,0};
	      localCoordinate_       = {0,0,0};
}
Cell::~Cell(void){
	  sphereNeighbors_.clear();
	  sampleNeighbors_.clear();
}

void Cell::addParticleID(long mcID){
	mcIDList_.push_back(mcID);
}

long Cell::getParticleID(long index){
    return mcIDList_[index];
}
void Cell::deleteParticleID(long mcID){
	long ID;

	for( decltype(mcIDList_.size()) i=0; i < mcIDList_.size(); i++){
		ID = mcIDList_[i];

		if(ID == mcID){
		      // assign last element from mcIDList_ at index to be deleted
			  mcIDList_[i] = mcIDList_.back();
		      // delete last element
			  mcIDList_.pop_back();
		      break;
		}
	}
}

void Cell::setCellNo(int cellNo){
	      cellNo_=cellNo;
}
void Cell::setGlobalCoordinate(IntVector3d cellGlobalCoordinate){
	      globalCoordinate_=cellGlobalCoordinate;
}
void Cell::setLocalCoordinate(IntVector3d  cellLocalCoordinate){
          localCoordinate_=cellLocalCoordinate;
}
void Cell::addSphereNeighbor(IntVector3d sphereNeighbor){
	sphereNeighbors_.push_back(sphereNeighbor);
}
void Cell::addSampleNeighbor(IntVector3d sampleNeighbor){
	sampleNeighbors_.push_back(sampleNeighbor);
}

void Cell::addShiftIndex(IntVector3d shiftIndex){
	shiftIndex_.push_back(shiftIndex);
}

int   Cell::getCellNo(){
	    return cellNo_;
}
long  Cell::getnParticlesCell(){
	    return mcIDList_.size();
}

int Cell::getSphereNBLSize(void){
	return sphereNeighbors_.size();
}
int Cell::getSampleNBLSize(void){
	return sampleNeighbors_.size();
}
IntVector3d  Cell::getGlobalCoordinate(){
	    return globalCoordinate_;
}
IntVector3d  Cell::getLocalCoordinate(){
        return localCoordinate_;
}
IntVector3d  Cell::getSphereNeighbor(int index){
	return sphereNeighbors_[index];
}
IntVector3d Cell::getSampleNeighbor(int index){
	return sampleNeighbors_[index];
}
IntVector3d   Cell::getShiftIndex(int index){
	return shiftIndex_[index];
}
