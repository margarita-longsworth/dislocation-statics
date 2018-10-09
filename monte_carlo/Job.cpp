/*
 * Job.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: ganeshfx
 */

#include "Job.h"
#include <iostream>
//#include "tr"
Job::Job() {}

Job::~Job() {}

Job::Job(long stepNo,Master& master,Domain& domain, long currentMCStep){
	reInitJob(0, master, domain, currentMCStep);
	StepNO   = stepNo;
}
void Job::modifyJobForSpatialConflict(Master& master,Domain& domain, long currentMCStep){
	if( repeatMode != SPATIAL || state != DISCARDED){
		std::cerr << "ERROR: JOB WAS NOT IN DISCARDED SPATIAL CONFLICT" << std::endl;
		exit(1);
	}
	reInitJob(versionCounter+1, master, domain, currentMCStep);
}
int Job::discardDependencies(deque<Job>& jobQueue){
	if (jobQueue.empty()) return 0;
	long headStep = jobQueue.front().StepNO;

	int counter=0;
	for(auto &p: dependancyList){
		Job& j = jobQueue[p.first-headStep];
		if(j.StepNO != p.first){
			cerr<<"ERROR: HEAD STEP NO DOES NOT MATCH "<<endl;
		    exit(1);
		}
		// Just to avoid marking a UNAVAILABLE job from repeated discarding
		// in case of conflict
		// State which is marked as "UNAVAILABLE" should not be over-written!
		if(j.versionCounter == p.second && j.state!= UNAVAILABLE){
			if (j.state!= DISCARDED ) counter++;
			j.state = DISCARDED;
			j.repeatMode = SPATIAL;
		}
//		//todo Just for testing
//		cout<<"##################################"<<endl;
//		cout<<" CALL: Job::discardDependencies " << endl;
//		cout<<" ACCESSED JOB STEP NO :"<<j.StepNO<<endl;
//		cout<<" ACCESSED JOB STATE   :"<<j.state <<endl;
//		cout<<"##################################"<<endl;
//		//todo Just for testing
	}
	return counter;
}

void Job::createDependencies(deque<Job>& jobQueue,Domain& domain){

	for(auto& mySwap : swapList){ // loop over the swap list of my Job
	    Vector3d refPosition = domain.getParticles()[mySwap.ID].getPosition();
	    auto mystate = std::make_pair(StepNO, versionCounter);

	    for(auto &job : jobQueue){ // loop over job queue

	    	//exclude my own Job from comparison and also needless to check overlapping with an UNAVAILABLE or DISCARDED job index in jobQueue.

			if(job.StepNO ==StepNO || job.state == UNAVAILABLE || job.state == DISCARDED)
				continue;

			    for(auto& otherSwap : job.swapList){ // loop over particle chosen for swapping for each Job in queue
					Vector3d  currPosition = domain.getParticles()[otherSwap.ID].getPosition();
					bool spatialConflict = domain.doEnergyZoneOverlap(refPosition, currPosition);

					if(spatialConflict){
						if( job.StepNO < StepNO){        // Jobs in early state
				    		// dependency list contains future jobs who are spatially conflicting
				    		job.dependancyList.push_back(mystate);
				    		break;
				        }
						else {   // Jobs in future state
//							//todo Just for testing
//							cout<<"##################################"<<endl;
//							cout<<" CALL: Job::createDependencies " << endl;
//							cout<<" ACCESSED JOB STEP NO :"<<job.StepNO<<endl;
//							cout<<" ACCESSED JOB STATE   :"<<job.state <<endl;
//							cout<<" CREATED DEPENDENCY FOR STEP NO :"<<StepNO<<endl;
//							cout<<"##################################"<<endl;
//							//todo Just for testing
							// construct New dependency list considering current state
							auto targetPair =  std::make_pair(job.StepNO,job.versionCounter);
				    		// add these future Jobs into my dependency list which are spatially conflicting
				            dependancyList.push_back(targetPair);
				            break;
					    }
					}
			    }

		} // End of Job queue loop

	}// End of swap list
}

void Job::clearInvalidatedJobData(void){
	clearLists(); //free memory for this invalidated job
}

void Job::reInitJob(int versionNo, Master& master,Domain& domain, long currentMCStep){
	clearLists();

	swapList               = master.getTrialMove()->getSwapList(domain, currentMCStep);

	// Random seed order check for different configuration for chosen carbon and virtual atom
	const vector<Particle>& refParticles = domain.getParticles();

	//################

    exertedMove            = master.getMoveType();
    myResult.deltaEnergy   = 0.0;
    state                  = READY;
    sites.nAssumedSamples  = domain.getSampleListSize();
    sites.nAssmumedTargets = domain.getTargetListSize();
    repeatMode             = NONE;
    versionCounter         = versionNo;
	myResult.result = false;
	myResult.particles.clear();
}

bool Job::testForStateConflict(Master& master,Domain& dom){

	if(exertedMove == SWAP || exertedMove == CLUSTER){
		return false;
	}
	long    nCurrent;
	long    nAssumed = 0;

	if(!(addedTargets.empty() && deletedTargets.empty())){

		std::uniform_real_distribution<double>   distribution(0.0,1.0);

		std::sort(addedTargets.begin(),addedTargets.end());
		std::sort(deletedTargets.begin(),deletedTargets.end());

		vector<long> addedListDiff;   // addedList   - deletedList
		vector<long> deletedListDiff; // deletedList - addedList

		// Perform Set operation to eliminate common entries
		std::set_difference(addedTargets.begin(), addedTargets.end(),
				deletedTargets.begin(), deletedTargets.end(), std::inserter(addedListDiff, addedListDiff.begin()));

	    std::set_difference(deletedTargets.begin(), deletedTargets.end(),
	    		addedTargets.begin(), addedTargets.end(),std::inserter(deletedListDiff, deletedListDiff.begin()));

	    addedTargets = addedListDiff;
	    deletedTargets = deletedListDiff;

	    // Check 1 Addition/Deletion
		if(swapList[0].type == TARGET){
			// Addition of Target particles, so assumed and current Samples have to be considered
			nAssumed   = sites.nAssumedSamples;
			nCurrent   = dom.getSampleListSize();
		}else{
			// Deletion of Target particles, so assumed and current Targets have to be considered
	        nAssumed  = sites.nAssmumedTargets;
	        nCurrent  = dom.getTargetListSize();
		}

		if( nAssumed > nCurrent ){        // Case 1
	         return false;
		}

		double probOfCorrectChoice = (static_cast<double>(nAssumed))/nCurrent;
		double p = distribution(master.generatorRand_);

	    if(p > probOfCorrectChoice){
	    	cout <<"nAssumed: " << nAssumed<<endl;
	    	cout <<"nCurrent: " << nCurrent<<endl;
	    	cout <<"random No:" << p <<endl;
	    	cout<<"probOfCorrectChoice"<<probOfCorrectChoice<<endl;
	        state             =  DISCARDED;
	    	repeatMode        =  STATE;      // Mark the repeat mode as State conflict
	    	return true;
	    }
	}
    return false;
}
void Job::modifyJobForStateConflict(Master& master,Domain& domain){

	if( repeatMode != STATE || state != DISCARDED){
		std::cerr << "ERROR: JOB WAS NOT IN DISCARDED STATE CONFLICT" << std::endl;
		exit(1);
	}
	if(exertedMove == SWAP || exertedMove == CLUSTER){
		cerr<<"ERROR: STATE CONFLICT MODIFICATION NOT APPLICABLE FOR SWAP/CLUSTER TRIAL MOVE"<<endl;
	    exit(1);
	}
	auto& changedList = (swapList[0].type == TARGET) ? deletedTargets : addedTargets;
	auto intDistro    =  std::uniform_int_distribution<int>(0,changedList.size()-1);
    int  index        =  intDistro(master.generatorRand_);
	swapList[0].ID    =  changedList[index];
	cout<<" Particle ID       : " <<swapList[0].ID<<endl;
	state             = READY;
	repeatMode        = NONE;
	versionCounter++;

	clearLists();
}
int Job::getVersionNo(){
	return versionCounter;
}

void Job::clearLists(){
	//TODO comment
	swapList.clear();
	dependancyList.clear();
	addedTargets.clear();
	deletedTargets.clear();
}

