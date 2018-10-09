/*
 * acceptanceCheck.c++
 *
 */
#include "Particle.h"
#include "acceptanceCheck.h"
#include <iostream>

acceptanceCheck::acceptanceCheck(void){
}
acceptanceCheck::~acceptanceCheck(void){

}
nvtIMD::nvtIMD(int inSeed,int realTypes,int totalTypes,double temperature, double cohesiveEnergy){
	realTypes_      = realTypes;
	totalTypes_     = totalTypes;
	temperature_    = temperature;
	cohesiveEnergy_ = cohesiveEnergy;
	generatorRandom.seed(inSeed*4);
}
nvtIMD::~nvtIMD(void){

}
JobResult nvtIMD::isAccepted(Domain& dom, vector<Particle>& currentConfiguration,vector<Swap>& swapList, long MCStep){

    JobResult jobResult_;
    bool decision=false;
    double totalEnergyOld,totalEnergyNew;
    double energyDifference, val;
    double transProb=0.0, localRandom = 0.0;
    int newType, switchFlag=0;
    KahanSum energyOld, energyNew;
    double k = dom.getKappa();
    double kR = dom.getKappaR();
    long insertID = swapList[0].ID;
    long deleteID = swapList[1].ID;
    long nacc_v = dom.getParticles()[insertID].getSuccessFactor();
    long nrej_v = dom.getParticles()[insertID].getFailureFactor();
    long nacc_c = dom.getParticles()[deleteID].getSuccessFactor();
    long nrej_c = dom.getParticles()[deleteID].getFailureFactor();
    double wt = exp( k*nacc_c - kR*nrej_c )*1.0 / exp( k*nacc_v - kR*nrej_v );

    int nAdded=0, nDeleted=0;

    for(decltype(swapList.size()) i=0; i<swapList.size(); i++){
    	if(swapList[i].type == TARGET)
    	     nAdded++;
    	else
    		nDeleted++;
    }
    // Cohesive Flag estimation (for Carbon Addition/Deletion )
    switchFlag = nAdded - nDeleted;
    // ( 1Joule = 6.2418076274889206e+18 eV)
    //boltz_factor = 1.3806504e-23 * 6.2418076274889206e+18
    double boltz  = 1.3806504e-23 * 6.2418076274889206e+18; // in eV/K
    //double tempK  = temperature_ * 11605.50; // temperature in K

    std::uniform_real_distribution<double> distribution(0.0,1.0);
    const vector<Particle>& domParticle = dom.getParticles();

    for(const auto& newParticle : currentConfiguration){
        newType     = newParticle.getType();
    	// Ignore Placeholders and Zone3 for Energy Computation
    	if( newType <(2*totalTypes_) ){
    		if (newType%realTypes_ != 2 ) energyNew.add(newParticle.getEpot()); //TODO replace 2 by SAMPLE and test
    		const Particle& oldParticle = domParticle[newParticle.getmcID()];
    		if(oldParticle.getType()%realTypes_ !=2)
    		   energyOld.add(oldParticle.getEpot());
    	}
    }

    totalEnergyOld = energyOld.getSum();
    totalEnergyNew = energyNew.getSum();

    // For Energy Computation
    // totalEnergyNew - contribution of real particles in Zone 1 & 2 with trial move
    // totalEnergyOld - contribution of real particles in Zone 1 & 2 with out trial move

    energyDifference = totalEnergyNew - (totalEnergyOld + switchFlag*cohesiveEnergy_);

    if(energyDifference < 0.0){
	dom.setDetailedBalanceA(0.0);
	dom.setDetailedBalanceB(0.0);
	dom.setDetailedBalanceAB(0.0);
    	decision = true;
    }

    else if(energyDifference >= 0.0){	
    	if(temperature_!=0.0){	
	   if ( dom.getSamplingMode()==0 || dom.getSamplingMode()==5 ){ 
              val = (-1 * energyDifference)/(temperature_*boltz);
	      transProb = exp(val); 
	      dom.setDetailedBalanceA(1.0);
	      dom.setDetailedBalanceB( exp(val) );
	      dom.setDetailedBalanceAB( transProb );
           }
           if ( dom.getSamplingMode()==1 || dom.getSamplingMode()==3 ){ //BSx or fixedBS
	      if ( MCStep >= dom.getEquilibrationStep() ){
   	         val = (-1 * energyDifference)/(temperature_*boltz);
	         transProb = wt*exp(val); 
                 dom.setDetailedBalanceA(wt);
	         dom.setDetailedBalanceB( exp(val) );
                 dom.setDetailedBalanceAB( transProb );
	      }
	      else{
                 val = (-1 * energyDifference)/(temperature_*boltz);
	         transProb = exp(val); 
	         dom.setDetailedBalanceA(1.0);
	         dom.setDetailedBalanceB( exp(val) );
	         dom.setDetailedBalanceAB( transProb );
	      }
           }
           if ( dom.getSamplingMode()==2 || dom.getSamplingMode()==4 ){ //BSx or fixed in tandem
		   if ( dom.getLocalBlock()==1 && dom.getBiasedBlock()==0 ){
	    	      val = (-1 * energyDifference)/(temperature_*boltz);
	    	      transProb = exp(val); 
		      dom.setDetailedBalanceA(1.0);
		      dom.setDetailedBalanceB(exp(val));
		      dom.setDetailedBalanceAB(transProb);
		   }
		   if ( dom.getLocalBlock()==0 && dom.getBiasedBlock()==1 ){
	    	      val = (-1 * energyDifference)/(temperature_*boltz);
	    	      transProb = wt*exp(val); 
		      dom.setDetailedBalanceA(wt);
		      dom.setDetailedBalanceB( exp(val) );
		      dom.setDetailedBalanceAB( transProb );
		   }	  
		   if ( dom.getLocalBlock()==0 && dom.getBiasedBlock()==0 ){
		      val = (-1 * energyDifference)/(temperature_*boltz);
		      transProb = exp(val); 
		      dom.setDetailedBalanceA(1.0);
		      dom.setDetailedBalanceB( exp(val) );
		      dom.setDetailedBalanceAB( transProb );
		   }
           }   	  
    	}
    	else if(temperature_ == 0.0){
    		// Molecular Statics case
            // transition probability
     	   transProb = 0.0;
    	}

		localRandom = distribution(generatorRandom);

    	if(transProb > localRandom){
    		 decision = true;
    	}
    }
    // Update JobResult
    jobResult_.result      = decision;
    jobResult_.deltaEnergy = energyDifference;

	return jobResult_;


}

JobResult nvtIMD::isAccepted_finiteTemp(double energyDiff){

	JobResult jobResult_;
	bool decision=false;
    double val,transProb=0.0, localRandom = 0.0;

    std::uniform_real_distribution<double> distribution(0.0,1.0);

    // Energy Difference is computed and readily available for finite temperature
    // cases from the worker, which could be used for evaluation

    double boltz  = 1.3806504e-23 * 6.2418076274889206e+18; // in eV/K

    if(energyDiff < 0.0){
    	decision = true;
    }
    else if(energyDiff >= 0.0){
    	if(temperature_!=0.0){
    	   val = (-1 * energyDiff)/(temperature_*boltz);
           // transition probability
    	   transProb = exp(val);
    	}
    	else{
            // transition probability
     	   transProb = 0.0;
    	}
		localRandom = distribution(generatorRandom);


    	if(transProb > localRandom){
    		 decision = true;
    	}
    }

    // Update JobResult
    jobResult_.result      = decision;
    jobResult_.deltaEnergy = energyDiff;

    return jobResult_;
}

