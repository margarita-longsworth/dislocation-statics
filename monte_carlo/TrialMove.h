/*
 * TrialMove.h
 *
 */

#ifndef MC___TRIALMOVE_H_
#define MC___TRIALMOVE_H_

#include <random>
#include <vector>
#include <iostream>
#include "UtilityMethods.h"
#include "Domain.h"

class TrialMove{	//Abstract class TrialMove
     public:
	    virtual vector<Swap> getSwapList(Domain& dom, long currentMCStep) = 0;//Because of this method IS abstract class
	    void setNtrialmove();	//Never defined but compiler does not complain
	    int getNtrialmove();
	    int getRandomNo();
	    virtual ~TrialMove();

     protected:
	    TrialMove(int nTrialmoves,int inSeed);
     private:
	    int nTrialmoves_;
	    std::mt19937 generatorTrial;
	    std::uniform_int_distribution<int> distributionTrial;
};

//addDeleteParticle is concrete class of TrialMove
//addDeleteParticle has to define getSwapList. This is for trial move =2, but for other moves, 
//other implementations of getSwapList should be defined. That's why you have this abstract class. 
class addDeleteParticle:public TrialMove{//addDeleteParticle is concrete class inherited from TrialMove
     public:
	    addDeleteParticle(int nTrialmoves,int inSeed):
	    	TrialMove(nTrialmoves,inSeed){};	//The constructor of addDeleteParticle is the same as in TrialMove, with those arguments. First you have to initialize TrialMove before you can call addDeleteParticle then. 
						

	    virtual vector<Swap> getSwapList(Domain& dom,long currentMCStep);	//Here is getSwapList defined
        virtual ~addDeleteParticle();
};

#endif /* MC___TRIALMOVE_H_ */
