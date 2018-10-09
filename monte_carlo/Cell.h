/*
 * Cell.h
 *
 */

#ifndef MC___CELL_H_
#define MC___CELL_H_
#include <vector>
#include "IntVector3d.h"
#include "VectorRange.h"

using namespace std;

class Cell{

    public:

	    Cell();
	    ~Cell();

	    void addParticleID(long mcID);
	    void deleteParticleID(long mcID);
	    void setCellNo(int cellNo);
        void setGlobalCoordinate(IntVector3d cellGlobalCoordinate);
        void setLocalCoordinate(IntVector3d  cellLocalCoordinate);
        void addSphereNeighbor(IntVector3d sphereNeighbor);
        void addSampleNeighbor(IntVector3d sampleNeighbor);

        void addShiftIndex(IntVector3d shiftIndex);

	    long                   getParticleID(long index);
	    int                    getCellNo();
	    long                   getnParticlesCell();
        int                    getSphereNBLSize();
        int                    getSampleNBLSize();
	    IntVector3d            getGlobalCoordinate();
	    IntVector3d            getLocalCoordinate();
        IntVector3d            getSphereNeighbor(int index);
        IntVector3d            getSampleNeighbor(int index);

        IntVector3d            getShiftIndex(int index);

    private:
		int                   cellNo_;
	    IntVector3d           globalCoordinate_;
		IntVector3d           localCoordinate_;
	    vector<long>          mcIDList_ ;
		vector<IntVector3d>   sampleNeighbors_;
		vector<IntVector3d>   sphereNeighbors_;
		vector<IntVector3d>   shiftIndex_;
};

#endif /* MC___CELL_H_ */
