/*
 * UtilityMethods.h
 *
 */

#ifndef MC___UTILITYMETHODS_H_
#define MC___UTILITYMETHODS_H_

#include "Vector3d.h"
#include "Particle.h"
#include <vector>
using namespace std;

// Machine States definition
// RUNNING     - Job that currently running on any worker/slave processor
// FINISHED    - Job that is completed by any worker processor
// DISCARDED   - Manager/Master invalidated that job irrespective of its decision due to conflict
// READY       - indicate the status that the job is ready for assignment to any worker processor
// UNAVAILABLE - indicate that job at this particular index is INVALID or EXPIRED for processing

enum State {RUNNING,FINISHED,DISCARDED,READY,UNAVAILABLE}; // The possible machine states
enum RepeatMode {STATE,SPATIAL,NONE};
enum Move {SINGLE,SWAP,CLUSTER};
enum Type {TARGET=1,SAMPLE=2};

class Swap{
    public:
	    long ID;
	    int  type;
};

class Sites{
     public:
	     long nAssumedSamples;
	     long nAssmumedTargets;
};
class KahanSum{

  private:
	 double sum;
	 double correction;
  public:
	 KahanSum(void);
	 void add(double d);
	 void add(vector<double> d);
	 double getSum();

};
class JobResult {
    public:
		  vector<Particle> particles;
		  bool     result;
		  bool     isConverged;
	      double   deltaEnergy;
};

class statisticsResult {
    public: 
	double	totalBoundedProbability;
	double	totalProbability;
	long	totalBoundedVacancies;
	long	totalBoundedAcceptances;
	long	totalBoundedRejections;


};

double computeDistanceSquare(Vector3d position1,Vector3d position2);

#endif /* MC___UTILITYMETHODS_H_ */
