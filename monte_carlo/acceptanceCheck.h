/*
 * acceptanceCheck.h
 *
 */

#ifndef MC___ACCEPTANCECHECK_H_
#define MC___ACCEPTANCECHECK_H_

#include <cmath>
#include <random>
#include "UtilityMethods.h"
#include "Domain.h"

using namespace std;

class acceptanceCheck{
    public:
	    virtual JobResult isAccepted(Domain& dom, vector<Particle>& currentConfiguration,vector<Swap>& swapList, long MCStep) = 0;
	    virtual JobResult isAccepted_finiteTemp(double energyDiff) = 0;
	    virtual ~acceptanceCheck();
		std::mt19937 generatorRandom;

    protected:
	    acceptanceCheck();
};

class nvtIMD:public acceptanceCheck{
    public:
	    nvtIMD(int inSeed,int realTypes,int totalTypes,double temperature, double cohesiveEnergy);
	    ~nvtIMD();
	    virtual JobResult isAccepted(Domain& dom, vector<Particle>& currentConfiguration,vector<Swap>& swapList, long MCStep);
	    virtual JobResult isAccepted_finiteTemp(double energyDiff);

    private:
	    int       realTypes_;
	    int       totalTypes_;
	    double   temperature_;
	    double   cohesiveEnergy_;
};

#endif /* MC___ACCEPTANCECHECK_H_ */
