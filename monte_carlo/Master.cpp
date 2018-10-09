/*
 * Master.c++
 *
 */

#include "Master.h"
#include <iostream>

using namespace std;

//trialMove is an abstract class (its virtual member is vector<Swap> getSwapList(Domain& dom) = 0)
//addDeleteParticle is concrete class of it. Okay. 
//Master can create and object called trialMove by calling the constructor of addDeleteParticle. 

Master::Master(int nTrials,int seed,long mc_nConflicts,long mc_nFinished,
		long mc_nSwaps,long mc_nAccepSwaps,long mc_nRejeSwaps,long mc_nAddition,long mc_nDeletion,
		long mc_naddAccepted, long mc_ndelAccepted, long mc_naddRejected, long mc_ndelRejected,
		int nSlaves,double mc_totDeltaEpot,long mc_biasedMovesPerformed,long mc_localMovesPerformed,Domain* globDomain,string statFileName): trialMove(nTrials, seed){

    nTrials_		= nTrials;
    nSlaves_            = nSlaves;
    masterDomain_       = globDomain;
    nConflicts          = mc_nConflicts;
    nFinished           = mc_nFinished;
    startStep           = mc_nFinished;
    biasedMovesPerformed= mc_biasedMovesPerformed;
    localMovesPerformed = mc_localMovesPerformed;
    nAssigned           = 0;
    //todo temporary global variables (nPeformed, nInvalidated might not be available in future versions)
    nPerformed          = 0;
    nInvalidated        = 0;
    //todo remove above global variables if necessary
    nSwaps              = mc_nSwaps;
    nAccepSwaps         = mc_nAccepSwaps;
    nRejeSwaps          = mc_nRejeSwaps;
    totDeltaEpot        = mc_totDeltaEpot;
    TAG_RUNJOB = 0; TAG_TERMINATE = 1;  TAG_INITIALIZE = 1; TAG_FREE=2;
    nAddition=mc_nAddition; nDeletion=mc_nDeletion; naddAccepted=mc_naddAccepted;
    ndelAccepted=mc_ndelAccepted; naddRejected=mc_naddRejected; ndelRejected=mc_ndelRejected;
    trialMove = addDeleteParticle(nTrials,seed);	//trialMove is an object type addDeleteParticle

    // new file writer part!!
    resultOut.open(statFileName, ios_base::out);

    if(nTrials==1){		//Move move gets assigned values
		move = SINGLE;}
    else if(nTrials==2){
		move = SWAP;}
    else{
		move = CLUSTER;}

    generatorRand_.seed(seed);
 }
Master::~Master(void){
     //delete trialMove;
}
Move Master::getMoveType(void){
	return move;
}
TrialMove* Master::getTrialMove(){
	return &trialMove;
}

Job* Master::findDiscardedJob(void){
	for (Job& job : jobQueue_){
		if(job.state == DISCARDED){
			return &job;
		}
	}
	return nullptr;
}

Job* Master::findDiscardedJobAndInvalidate(void){
	for (Job& job : jobQueue_){
		if(job.state == DISCARDED){   // find existing discarded Job in jobQueue
//			//todo Just for testing
//			cout<<"##################################"<<endl;
//			cout<<" CALL: Job::findDiscardedJobAndInvalidate " << endl;
//			cout<<" ACCESSED JOB STEP NO :"<<job.StepNO<<endl;
//			cout<<" ACCESSED JOB STATE   :"<<job.state <<endl;
//			//todo Just for testing
			job.state = UNAVAILABLE; // indicate that job at this particular index is INVALID or EXPIRED for further processing
//			cout<<" MODIFIED JOB STATE   :"<<job.state <<endl;
//			cout<<"##################################"<<endl;
//			//todo
			nInvalidated++;             // temporary global counter for cross-checking nAssigned
			job.clearInvalidatedJobData(); // No future Jobs should depend on UNAVAILABLE job
			return &job;
		}
	}
	return nullptr;
}

bool Master::isQueueHeadReady(){
	if (jobQueue_.empty()) return false;
//	//todo Just for testing
//	cout<<"##################################"<<endl;
//	cout<<" CALL: Job::isQueueHeadReady " << endl;
//	cout<<" ACCESSED JOB STEP NO :"<<jobQueue_.front().StepNO<<endl;
//	cout<<" ACCESSED JOB STATE   :"<<jobQueue_.front().state <<endl;
//	//todo Just for testing
	// return true if the head job is 'FINISHED' or 'UNAVAILABLE'
	return (jobQueue_.front().state == FINISHED || jobQueue_.front().state == UNAVAILABLE	);
}

Job* Master::getNextJob(){
	long  totalMCSteps = masterDomain_->getMCTotalSteps();

	long n = 0; //NO of RUNNING or FINISHED jobs in jobQueue
	for(const auto& job : jobQueue_){
		 // ignore DISCARDED and UNAVAILABLE jobs for counting
		 if(job.state!=DISCARDED && job.state!=UNAVAILABLE) n++;
	}

	if((n+nFinished) < totalMCSteps ){
		//todo in future stable version findDiscardedJob() might not be available
		//Job*  jobPtr      = findDiscardedJob();
		Job*  jobPtr      = findDiscardedJobAndInvalidate(); // return null pointer for FINISHED and UNAVAILABLE jobs
	    bool  isDiscardedAndInvalidated = (jobPtr != nullptr);

	    if(isDiscardedAndInvalidated){
//	    	//todo Just for testing
//	    	cout<<"##################################"<<endl;
//	    	cout<<" CALL: Master::getNextJob() -- isDiscardedAndInvalidated" << endl;
//	    	cout<<" ACCESSED JOB STEP NO :"<<jobPtr->StepNO<<endl;
//	    	cout<<" ACCESSED JOB STATE   :"<<jobPtr->state <<endl;
//	    	cout<<"##################################"<<endl;
//	    	//todo Just for testing

	        // Re-insertion of Job at same location in job queue is not required any more
            // therefore, update the statistical counters, create a new job with incremented counter and
	    	// append it to the end of job queue
	        nConflicts++;
	    }
	    //  Assign new job and push it to job queue
	    Job job_new = Job((nAssigned+1),*this,*masterDomain_,nFinished);
		jobQueue_.push_back(job_new);
		Job* jobPtr_new  = &(jobQueue_.back());
		//todo - nPerformed might not be available in future versions
		nPerformed++; // temporary global counter for cross-checking nAssigned
		nAssigned++;

		jobPtr_new->createDependencies(jobQueue_,*masterDomain_);
		return jobPtr_new;
	}
	else{
		return nullptr;
	}
}

void Master::sendJobToWorker(Job& inJob,int processID){

	const vector<Particle>& globalList   = masterDomain_->getParticles();
    string fileName = "twoSphere_config.chkpt";
	vector<Particle>   configHolder;

    masterDomain_->createConfiguration(inJob.swapList, configHolder);

    for (const auto& s : inJob.swapList){ // adding chosen Particle into configuration
    	Particle  particle  = globalList[s.ID];

/* ALWAYS SWITCH THE TYPES RIGHT AWAY: NO FINITE TEMPERATURE LIBRARY IS USED*/ 

    	particle.setType(s.type);


/*
    	// enable switch
    	// change types of chosen particle before sending it to worker (only for T=0K)
    	// for FINITE TEMP , chosen particle types will be changed during transmutation inside imd-library
    	if(masterDomain_->getTemperature() == 0.0){
    		particle.setType(s.type);
    	}
*/

    	//cout<<" Particle ID : "<<particle.getNumber()<<" and type "<<particle.getType()<<endl;
    	configHolder.push_back(particle);
    }

    // update particle freq for
    // visualizing particles which involved in sphere construction 
//    for (auto& p : configHolder){
//    	globalList[p.getmcID()].updateFrequency(1);
//    }

    //masterDomain_->writeSphereConfiguration(fileName,configHolder);
    long bufferSize   = (configHolder.size()*7); // Contain buffer Size and flip Type

	// ----  Fill buffer with entries from constructed sphere
	double* jobBuffer;
	jobBuffer            = new double [bufferSize];
 	long* bufferSizeList = new long [5];

    bufferSizeList[0]   = bufferSize;
    bufferSizeList[1]   = inJob.StepNO;
    bufferSizeList[2]   = inJob.getVersionNo();

    long rtoVID =-1, vtoRID=-1;

    // for T=0 K cases, make rtoVID and vtoRID as -1
    // otherwise get ID values from swapList

 /* NEVER CALL THE FINITE LIBRARY I.E. JUST LEAVE rtoVID and vtoRID as -1  */

/*

    double Temp = masterDomain_->getTemperature();

    int count=1;

    if(Temp != 0.0){ // finiteTemp
    	for(const auto& s: inJob.swapList){
    		//cout << " I come here for time : "<<count<<endl;

    		if(globalList[s.ID].getType() == TARGET){
//    			cout <<"inside sample loop"<<endl;
    			//cout << " rToVID : "<<s.ID<<endl;
    			rtoVID = s.ID;
    		}
    		else if(globalList[s.ID].getType() == SAMPLE){
//    			cout <<"inside target loop"<<endl;
    			//cout << " vToRID : "<<s.ID<<endl;
    			vtoRID = s.ID;
    		}
    		else{
    			cerr<<"ERROR: Unknown swap type encountered for finiteTemp case with type :"<<globalList[s.ID].getType() <<endl;
    			exit(1);
    		}
    		count++;
    	}
    }

*/

	bufferSizeList[3] = rtoVID;
	bufferSizeList[4] = vtoRID;

    long bufCounter=0;

    for(const auto& p : configHolder){ // Pack configuration contents into buffer
           jobBuffer[bufCounter++]   = static_cast<double> (p.getmcID());
           jobBuffer[bufCounter++]   = static_cast<double> (p.getType());
           jobBuffer[bufCounter++]   = p.getMass();
           jobBuffer[bufCounter++]   = p.getPosition().x;
           jobBuffer[bufCounter++]   = p.getPosition().y;
           jobBuffer[bufCounter++]   = p.getPosition().z;
           jobBuffer[bufCounter++]   = p.getEpot();
	}
    MPI_Send(bufferSizeList, 5, MPI_LONG, processID, TAG_RUNJOB, MPI_COMM_WORLD);
    MPI_Send(jobBuffer, bufferSize, MPI_DOUBLE, processID, TAG_RUNJOB, MPI_COMM_WORLD);

    delete[] jobBuffer;
    delete[] bufferSizeList;

    inJob.state      = RUNNING;
    inJob.nParticles = configHolder.size();
    inJob.workerID   = processID;

//    //todo Testing part
//    cout<< "***************************************************"<<endl;
//    cout <<"CALL: SEND JOB TO WORKER"<<endl;
//    cout << "Master assigned Job to process rank " << processID << " for step no "
//         << inJob.StepNO << " with " << inJob.nParticles << " particles" << endl;
//    cout <<"Chosen ID :"<<inJob.swapList[0].ID<<" and its Type :"<<inJob.swapList[0].type<<endl;
//    cout <<"Chosen ID :"<<inJob.swapList[1].ID<<" and its Type :"<<inJob.swapList[1].type<<endl;
//    cout<< "***************************************************"<<endl;
//    //todo Testing part
}

void Master::deleteJobQueueHead(void){
//	   //todo
//		cout<<"##################################"<<endl;
//		cout<<" CALL: Master::deleteJobQueueHead " << endl;
//		cout<<" current jobQueue_ size : "<< jobQueue_.size()<<endl;
//		cout<<" TO BE DELETED JOB ID :"<<jobQueue_.front().StepNO<<endl;
//
//	   //todo
	    jobQueue_.pop_front();
//	    cout<<"POST CHECK jobQueue_ size : "<< jobQueue_.size()<<endl;
//		cout<<"##################################"<<endl;
	   //todo
}

void Master::fetchResult(int slaveRank){

	// This method obtain the results from the given Process Rank and store the
	// configurations at corresponding index in jobQueue via member result
    int source = slaveRank;
    MPI_Status status;
    double recResult[5];
    // Receive results from worker process
    MPI_Recv(recResult, 5, MPI_DOUBLE, source,TAG_RUNJOB, MPI_COMM_WORLD, &status);

    long targetStepNo = static_cast<long> (recResult[1]);
    long targetVerNo  = static_cast<long> (recResult[2]);
    long recBufSize   = static_cast<long> (recResult[3]);

/*

    double* eDiffPtr;
    eDiffPtr=nullptr;

     double eDiff;

    // for finiteTemp case energyDiff is computed by worker, therefore store the collected result
    if(masterDomain_->getTemperature() != 0.0){
         eDiff    = recResult[4];
         //cout <<"Master fetched eDiff :"<< eDiff<<endl;
         eDiffPtr = &eDiff;
    }

*/

    bool isRequired = false;
	Job& headJob = jobQueue_.front();

	if(targetStepNo >= headJob.StepNO){ // otherwise Job is no longer exist in queue

		int index = targetStepNo - headJob.StepNO;
		Job& job = jobQueue_[index];

		if (job.StepNO != targetStepNo){
			cerr << "ERROR: Wrong StepNO" << endl;
			exit(1);
		}

		// Receive result configurations from Worker, only if the Job state is NOT marked as DISCARDED and UNAVAILABLE
		if( (job.getVersionNo() == targetVerNo) && (job.state != DISCARDED) && (job.state != UNAVAILABLE) ){
            isRequired = true;
		    //send Requirement Flag
			MPI_Send(&isRequired,1, MPI_C_BOOL, source, TAG_RUNJOB, MPI_COMM_WORLD);

			// allocate receive buffer
			double* receiveBuffer = new double [recBufSize];

			//import buffer
			MPI_Recv(receiveBuffer,recBufSize, MPI_DOUBLE, source,TAG_RUNJOB, MPI_COMM_WORLD, &status);

			jobQueue_[index].myResult.particles.clear();
            jobQueue_[index].myResult.isConverged = (static_cast<int>(recResult[0])==1) ;

/*
            // for finiteTemp case
            if(eDiffPtr != nullptr){
            	jobQueue_[index].myResult.deltaEnergy = *eDiffPtr;
            	//cout<<"Master assigned eDiff to jobqueue :"<<jobQueue_[index].myResult.deltaEnergy<<endl;
            }
*/

		    long recbufind = 0l;
			long particleCounter = recBufSize/7l;

            if(jobQueue_[index].nParticles != particleCounter){
                cerr<<" total sent     :" <<jobQueue_[index].nParticles<<" to worker "<<jobQueue_[index].workerID << endl;
                cerr<<" total received :" <<particleCounter <<" from worker "<<slaveRank<<endl;
            	cerr<<"ERROR: occurred with different number of send and receive particles for step NO"<<jobQueue_[index].StepNO<<endl;
            	exit(1);
            }

			// fill back buffer into particle container
			for(auto i=0; i<particleCounter; i++){
			    Particle particleHolder;
			    particleHolder.setmcID(receiveBuffer[recbufind++]);
			    particleHolder.setType(receiveBuffer[recbufind++]);
		        particleHolder.setMass(receiveBuffer[recbufind++]);

		        double recPosX = receiveBuffer[recbufind++];
		        double recPosY = receiveBuffer[recbufind++];
		        double recPosZ = receiveBuffer[recbufind++];

			    particleHolder.setPosition(recPosX,recPosY,recPosZ);
			    particleHolder.setEpot(receiveBuffer[recbufind++]);
			    jobQueue_[index].myResult.particles.push_back(particleHolder);
			}

			// Mark job as finished if the machine state is running (otherwise retain its state)
			if(job.state == RUNNING)
				 job.state = FINISHED;

			// deallocate buffer
			delete[] receiveBuffer;
			return;
		}
	}
	//send Requirement Flag
	MPI_Send(&isRequired,1, MPI_C_BOOL, source, TAG_RUNJOB, MPI_COMM_WORLD);
}

void Master::evaluateDecision(Job& inJob){

	 //todo KEEP THE acceptObject initialization inside Master constructor
     int       inSeed        = masterDomain_->getUserSeed();
     int       RealTypes     = masterDomain_->getRealTypes();
     int       TotalTypes    = masterDomain_->getTotalTypes();
     double    Temperature   = masterDomain_->getTemperature();
     double    Cohesive      = masterDomain_->getCohesiveEnergy();

     static nvtIMD    acceptObject(inSeed,RealTypes,TotalTypes,Temperature,Cohesive); //TODO MUST be updated !
     JobResult jobResult;

	/* ALWAYS USE THE T=0K LIBRARY FOR RELAXATION EVEN FOR HIGHER TEMPERATURE SIMULATIONS */ 

     jobResult = acceptObject.isAccepted(*masterDomain_,inJob.myResult.particles,inJob.swapList,nFinished );
     inJob.myResult.deltaEnergy = jobResult.deltaEnergy; // Delta energy computed inside isAccepted()

/* 

     if(Temperature==0.0){ // evaluate acceptance criterion and get decision (for T=0K case!)
     	  jobResult = acceptObject.isAccepted(*masterDomain_,
   			     inJob.myResult.particles,inJob.swapList);
     	  inJob.myResult.deltaEnergy = jobResult.deltaEnergy; // Delta energy computed inside isAccepted()
     }
     else{
    	 //cout<<"Master evaluate Decision for finiteTemp case"<<endl;
    	 // Delta energy for finiteTemp case computed by worker and send to Manager
    	 jobResult = acceptObject.isAccepted_finiteTemp(inJob.myResult.deltaEnergy);
    	 string RESULT=(jobResult.result)?"ACCEPTED":"REJECTED";

    	 //cout<<" And the evaluated decision : "<< RESULT<<endl;
     }

*/

	 // consider the convergence criterion for T=0K case (for T>0K isCoverged always true)
	 inJob.myResult.result      = ((inJob.myResult.isConverged) && (jobResult.result)); //CRUCIAL STEP



}

bool Master::processJob(Job& inJob){

   vector<Particle> sphereAroundSelectedSamplingSite;

   bool  isStateConflict = inJob.testForStateConflict(*this,*masterDomain_);

   long nDiscarded = 0;

   if (!isStateConflict){ 

      if(inJob.state == FINISHED){

	if ( masterDomain_->getSamplingMode()==1 ){ //BSx

		if ( nFinished >= masterDomain_->getEquilibrationStep() ){

			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList;
			vector<Particle> activeParticles; 
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double accRatio = 0;
			long maxAcc = 0;
			long maxRej = 0;
			long totalSuccess = 0;
			long totalFailure = 0;
			int activeVacancies = 0;
			int inactiveVacancies = 0;
			int totalSuccessInActiveRegion = 0;
			int vacanciesSameMaximum = 0;
			double avgSuccessInActiveRegion = 0;
			double avgFailureInActiveRegion = 0;
			double avgSuccessInFixedRegion = 0;
			double avgFailureInFixedRegion = 0;
			int totalFailureInActiveRegion = 0;
			double ratio = 0;
			double largestRatio = 0;
			long totalVacancies = masterDomain_->getSampleListSize();
			double biasStrength = 0;	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalSuccess += p.getSuccessFactor();
				}
			}	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalFailure += p.getFailureFactor();
				}
			}

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){	
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
						maxSuccessFactor = configHolder[i].getSuccessFactor(); 
					}
				}
			}	

			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}

			for ( const auto& p : configHolder ){							
				typeHold = p.getType();
				if (typeHold == 2){	
					if ( p.getSuccessFactor() == maxSuccessFactor ){
						activeVacancies += 1;	
						vacanciesSameMaximum += 1;	
						totalSuccessInActiveRegion += p.getSuccessFactor(); 
						totalFailureInActiveRegion += p.getFailureFactor(); 
						accRatio += p.getSuccessFactor()*1.0/totalSuccess;
						if ( accRatio >= 0.10 ){
							break;
						} 		
					}
				}
			}
	
			while ( accRatio <= 0.10 ) {

				vacanciesSameMaximum = 0;
				maxSuccessFactor -= 1.0; 

				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == maxSuccessFactor ){
							activeVacancies += 1;	
							vacanciesSameMaximum += 1;
							totalSuccessInActiveRegion += p.getSuccessFactor(); 	
							totalFailureInActiveRegion += p.getFailureFactor(); 
							accRatio += p.getSuccessFactor()*1.0/totalSuccess;
							if ( accRatio >= 0.10 ){
								break;
							} 		
						}
					}
				}

			}	
	
			inactiveVacancies = totalVacancies - activeVacancies;

			avgSuccessInActiveRegion = totalSuccessInActiveRegion*1.0/activeVacancies;
			avgFailureInActiveRegion = totalFailureInActiveRegion*1.0/activeVacancies;
			avgSuccessInFixedRegion = (totalSuccess - totalSuccessInActiveRegion)*1.0/inactiveVacancies;
			avgFailureInFixedRegion = (totalFailure - totalFailureInActiveRegion)*1.0/inactiveVacancies;

			masterDomain_-> setActiveVacancies ( activeVacancies );
			masterDomain_-> setInactiveActiveRatio ( inactiveVacancies*1.0/activeVacancies );
			masterDomain_-> setActiveRegionSuccess ( avgSuccessInActiveRegion );
			masterDomain_-> setActiveRegionFailure ( avgFailureInActiveRegion );
			masterDomain_-> setInactiveRegionSuccess ( avgSuccessInFixedRegion );
			masterDomain_-> setInactiveRegionFailure ( avgFailureInFixedRegion );

			biasStrength = maxAcc* ( log (0.3*inactiveVacancies*1.0/ activeVacancies ) + ( avgFailureInActiveRegion - avgFailureInFixedRegion )*1.0/maxRej )*1.0/ ( avgSuccessInActiveRegion -  avgSuccessInFixedRegion  ) ;
		 	
			masterDomain_->setBiasStrength( biasStrength );
			masterDomain_->setMaxSuccess( maxAcc ); 
			masterDomain_->setMaxFailure( maxRej ); 

			masterDomain_->setKappa( biasStrength / maxAcc );
			masterDomain_->setKappaR ( 1.0 / maxRej );
		
		}
	}	


	if ( masterDomain_->getSamplingMode()==2 ){ //BSx_tandem

		if ( masterDomain_->getBiasedBlock()==1 && masterDomain_->getLocalBlock()==0 ){

			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList;
			vector<Particle> activeParticles; 
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double accRatio = 0;
			long maxAcc = 0;
			long maxRej = 0;
			long totalSuccess = 0;
			long totalFailure = 0;
			int activeVacancies = 0;
			int inactiveVacancies = 0;
			int totalSuccessInActiveRegion = 0;
			int vacanciesSameMaximum = 0;
			double avgSuccessInActiveRegion = 0;
			double avgFailureInActiveRegion = 0;
			double avgSuccessInFixedRegion = 0;
			double avgFailureInFixedRegion = 0;
			int totalFailureInActiveRegion = 0;
			double ratio = 0;
			double largestRatio = 0;
			long totalVacancies = masterDomain_->getSampleListSize();
			double biasStrength = 0;	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalSuccess += p.getSuccessFactor();
				}
			}	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalFailure += p.getFailureFactor();
				}
			}

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){	
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
						maxSuccessFactor = configHolder[i].getSuccessFactor(); 
					}
				}
			}	

			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}

			for ( const auto& p : configHolder ){							
				typeHold = p.getType();
				if (typeHold == 2){	
					if ( p.getSuccessFactor() == maxSuccessFactor ){
						activeVacancies += 1;	
						vacanciesSameMaximum += 1;	
						totalSuccessInActiveRegion += p.getSuccessFactor(); 
						totalFailureInActiveRegion += p.getFailureFactor(); 
						accRatio += p.getSuccessFactor()*1.0/totalSuccess;
						if ( accRatio >= 0.10 ){
							break;
						} 		
					}
				}
			}
	
			while ( accRatio <= 0.10 ) {

				vacanciesSameMaximum = 0;
				maxSuccessFactor -= 1.0; 

				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == maxSuccessFactor ){
							activeVacancies += 1;	
							vacanciesSameMaximum += 1;
							totalSuccessInActiveRegion += p.getSuccessFactor(); 	
							totalFailureInActiveRegion += p.getFailureFactor(); 
							accRatio += p.getSuccessFactor()*1.0/totalSuccess;
							if ( accRatio >= 0.10 ){
								break;
							} 		
						}
					}
				}

			}	
	
			inactiveVacancies = totalVacancies - activeVacancies;

			avgSuccessInActiveRegion = totalSuccessInActiveRegion*1.0/activeVacancies;
			avgFailureInActiveRegion = totalFailureInActiveRegion*1.0/activeVacancies;
			avgSuccessInFixedRegion = (totalSuccess - totalSuccessInActiveRegion)*1.0/inactiveVacancies;
			avgFailureInFixedRegion = (totalFailure - totalFailureInActiveRegion)*1.0/inactiveVacancies;

			masterDomain_-> setActiveVacancies ( activeVacancies );
			masterDomain_-> setInactiveActiveRatio ( inactiveVacancies*1.0/activeVacancies );
			masterDomain_-> setActiveRegionSuccess ( avgSuccessInActiveRegion );
			masterDomain_-> setActiveRegionFailure ( avgFailureInActiveRegion );
			masterDomain_-> setInactiveRegionSuccess ( avgSuccessInFixedRegion );
			masterDomain_-> setInactiveRegionFailure ( avgFailureInFixedRegion );

			biasStrength = maxAcc* ( log (0.3*inactiveVacancies*1.0/ activeVacancies ) + ( avgFailureInActiveRegion - avgFailureInFixedRegion )*1.0/maxRej )*1.0/ ( avgSuccessInActiveRegion -  avgSuccessInFixedRegion  ) ;
		 	
			masterDomain_->setBiasStrength( biasStrength );
			masterDomain_->setMaxSuccess( maxAcc ); 
			masterDomain_->setMaxFailure( maxRej ); 

			masterDomain_->setKappa( biasStrength / maxAcc );
			masterDomain_->setKappaR ( 1.0 / maxRej );

		}
	}	
	
	if ( masterDomain_->getSamplingMode()==3 ){ //fixedBS

		if ( nFinished >= masterDomain_->getEquilibrationStep() ){

			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList; 
			std::random_shuffle ( configHolder.begin(), configHolder.end() );
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double largestRatio = 0; 
			double ratio = 0;
			long maxAcc = 0;
			long maxRej = 0;

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){
				if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
					maxSuccessFactor = configHolder[i].getSuccessFactor(); 
				}
			}	
			
			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}	

			masterDomain_->setKappa(8.0 / maxAcc );
			masterDomain_->setKappaR ( 1.0 / maxRej );
				
		}
	}	

	if ( masterDomain_->getSamplingMode()==4 ){ //fixedBS_tandem

		if ( masterDomain_->getBiasedBlock()==1 && masterDomain_->getLocalBlock()==0 ){

			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList; 
			std::random_shuffle ( configHolder.begin(), configHolder.end() );
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double largestRatio = 0; 
			double ratio = 0;
			long maxAcc = 0;
			long maxRej = 0;

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){
				if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
					maxSuccessFactor = configHolder[i].getSuccessFactor(); 
				}
			}	
			
			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}	

			masterDomain_->setKappa(8.0 / maxAcc );
			masterDomain_->setKappaR ( 1.0 / maxRej );
				
		}
	}	

	if ( masterDomain_->getSamplingMode()==5 ){ //Bell algorithm

		if ( nFinished >= masterDomain_->getEquilibrationStep() ){
			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList;
			vector<Particle> activeParticles; 
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double accRatio = 0;
			long maxAcc = 0;
			long maxRej = 0;
			long totalSuccess = 0;
			long totalFailure = 0;
			int activeVacancies = 0;
			int totalSuccessInActiveRegion = 0;
			int vacanciesSameMaximum = 0;
			double avgSuccessInActiveRegion = 0;
			double avgFailureInActiveRegion = 0;
			int totalFailureInActiveRegion = 0;
			double ratio = 0;
			double largestRatio = 0;

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalSuccess += p.getSuccessFactor();
				}
			}	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalFailure += p.getFailureFactor();
				}
			}

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){	
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
						maxSuccessFactor = configHolder[i].getSuccessFactor(); 
					}
				}
			}	

			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}

			for ( const auto& p : configHolder ){							
				typeHold = p.getType();
				if (typeHold == 2){	
					if ( p.getSuccessFactor() == maxSuccessFactor ){
						activeVacancies += 1;	
						vacanciesSameMaximum += 1;	
						totalSuccessInActiveRegion += p.getSuccessFactor(); 
						totalFailureInActiveRegion += p.getFailureFactor(); 
						accRatio += p.getSuccessFactor()*1.0/totalSuccess;
						if ( accRatio >= 0.15 ){
							break;
						} 		
					}
				}
			}

			while ( accRatio <= 0.15 ) {

				vacanciesSameMaximum = 0;
				maxSuccessFactor -= 1.0; 

				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == maxSuccessFactor ){
							activeVacancies += 1;	
							vacanciesSameMaximum += 1;
							totalSuccessInActiveRegion += p.getSuccessFactor(); 	
							totalFailureInActiveRegion += p.getFailureFactor(); 
							accRatio += p.getSuccessFactor()*1.0/totalSuccess;
							if ( accRatio >= 0.15 ){
								break;
							} 		
						}
					}
				}

			}	
	
			avgSuccessInActiveRegion = totalSuccessInActiveRegion*1.0/activeVacancies;
			avgFailureInActiveRegion = totalFailureInActiveRegion*1.0/activeVacancies;

			int vacanciesInBellCenter = 0;
			double totalRejectionsUnderBell = 0;
			double avgRejectionsUnderBell = 0;
			int bellCenter = masterDomain_->getBellCenter();
			int instantSuccess = masterDomain_->getInstantAvgAcc();
			int instantFailure = masterDomain_->getInstantAvgRej();
			bool IDNotFound = true;

			if ( bellCenter == 0 ){
				masterDomain_->setBellCenter( maxAcc );
				bellCenter = masterDomain_->getBellCenter(); 
				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == bellCenter ){
								vacanciesInBellCenter += 1;	
								totalRejectionsUnderBell += p.getFailureFactor(); 
						}
					}
				}
				avgRejectionsUnderBell = totalRejectionsUnderBell*1.0/vacanciesInBellCenter;
				masterDomain_->setAvgRejectionsUnderBell( avgRejectionsUnderBell );	
				masterDomain_->setInstantAvgAcc ( floor(avgSuccessInActiveRegion) );
				masterDomain_->setInstantAvgRej ( maxRej );
				instantSuccess = masterDomain_->getInstantAvgAcc();
				instantFailure = masterDomain_->getInstantAvgRej();
				if ( (avgRejectionsUnderBell >= instantFailure) && (bellCenter == floor(instantSuccess)) ){
					masterDomain_->setBellCenter ( maxAcc );
				}
				else if ( avgRejectionsUnderBell >= instantFailure ){
					masterDomain_->setBellCenter( bellCenter-1 );
				}
			}	
			else{
				while(IDNotFound){
					bellCenter = masterDomain_->getBellCenter(); 
					for ( const auto& p : configHolder ){							
						typeHold = p.getType();
						if (typeHold == 2){	
							if ( p.getSuccessFactor() == bellCenter ){
									vacanciesInBellCenter += 1;	
									totalRejectionsUnderBell += p.getFailureFactor(); 
									IDNotFound = false;
							}
						}
					}
					if (IDNotFound){
						bellCenter = bellCenter - 1;
						masterDomain_->setBellCenter ( bellCenter );
					}			
				}

				avgRejectionsUnderBell = floor(totalRejectionsUnderBell*1.0/vacanciesInBellCenter);	
				masterDomain_->setAvgRejectionsUnderBell( avgRejectionsUnderBell );
				if ( (avgRejectionsUnderBell >= instantFailure) && (bellCenter == floor(instantSuccess)) ){
					masterDomain_->setBellCenter ( maxAcc );
					masterDomain_->setInstantAvgAcc ( floor(avgSuccessInActiveRegion) );
					masterDomain_->setInstantAvgRej ( maxRej );
				}
				else if ( avgRejectionsUnderBell >= instantFailure ){
					masterDomain_->setBellCenter( bellCenter-1 );
				}
			}
			
			masterDomain_->setMaxSuccess( maxAcc ); 
			masterDomain_->setMaxFailure( maxRej ); 
			masterDomain_->setAvgSuccess(avgSuccessInActiveRegion ); 
		       	masterDomain_->setAvgFailure( avgFailureInActiveRegion );

		}
	}

	if ( masterDomain_->getSamplingMode()==6 ){ //Bell-tandem

		if ( masterDomain_->getBiasedBlock()==1 && masterDomain_->getLocalBlock()==0 ){

			const vector<Particle>& globalList = masterDomain_->getParticles();
			vector<Particle> configHolder = globalList;
			vector<Particle> activeParticles; 
			long maxSuccessFactor = 0;
			int typeHold = 0;	
			int maxIndex = 0; 
			double accRatio = 0;
			long maxAcc = 0;
			long maxRej = 0;
			long totalSuccess = 0;
			long totalFailure = 0;
			int activeVacancies = 0;
			int totalSuccessInActiveRegion = 0;
			int vacanciesSameMaximum = 0;
			double avgSuccessInActiveRegion = 0;
			double avgFailureInActiveRegion = 0;
			int totalFailureInActiveRegion = 0;
			double ratio = 0;
			double largestRatio = 0;

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalSuccess += p.getSuccessFactor();
				}
			}	

			for ( auto& p: configHolder ){	
				typeHold = p.getType();
				if (typeHold == 2){	
					totalFailure += p.getFailureFactor();
				}
			}

			for ( decltype(configHolder.size()) i = 0 ; i < configHolder.size(); i++ ){
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					maxSuccessFactor = configHolder[i].getSuccessFactor();
					maxIndex = i;
					break; 
				}
			}

			for ( decltype(configHolder.size()) i = maxIndex ; i < configHolder.size(); i++ ){	
				typeHold = configHolder[i].getType();
				if (typeHold == 2){	
					if ( maxSuccessFactor < configHolder[i].getSuccessFactor() ){
						maxSuccessFactor = configHolder[i].getSuccessFactor(); 
					}
				}
			}	

			vector<Particle> mostImportantVacancies;

			for ( const auto& p : configHolder ){
				if ( p.getSuccessFactor() == maxSuccessFactor ){
					mostImportantVacancies.push_back(p);
				}
			}

			for ( auto& p : mostImportantVacancies ){
				ratio = p.getSuccessFactor()*1.0 / p.getFailureFactor();
				if ( largestRatio < ratio ){
					largestRatio = ratio;
					maxAcc = p.getSuccessFactor();
					maxRej = p.getFailureFactor();
				}
			
			}

			for ( const auto& p : configHolder ){							
				typeHold = p.getType();
				if (typeHold == 2){	
					if ( p.getSuccessFactor() == maxSuccessFactor ){
						activeVacancies += 1;	
						vacanciesSameMaximum += 1;	
						totalSuccessInActiveRegion += p.getSuccessFactor(); 
						totalFailureInActiveRegion += p.getFailureFactor(); 
						accRatio += p.getSuccessFactor()*1.0/totalSuccess;
						if ( accRatio >= 0.15 ){
							break;
						} 		
					}
				}
			}

			while ( accRatio <= 0.15 ) {

				vacanciesSameMaximum = 0;
				maxSuccessFactor -= 1.0; 

				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == maxSuccessFactor ){
							activeVacancies += 1;	
							vacanciesSameMaximum += 1;
							totalSuccessInActiveRegion += p.getSuccessFactor(); 	
							totalFailureInActiveRegion += p.getFailureFactor(); 
							accRatio += p.getSuccessFactor()*1.0/totalSuccess;
							if ( accRatio >= 0.15 ){
								break;
							} 		
						}
					}
				}

			}	
	
			avgSuccessInActiveRegion = totalSuccessInActiveRegion*1.0/activeVacancies;
			avgFailureInActiveRegion = totalFailureInActiveRegion*1.0/activeVacancies;

			int vacanciesInBellCenter = 0;
			double totalRejectionsUnderBell = 0;
			double avgRejectionsUnderBell = 0;
			int bellCenter = masterDomain_->getBellCenter();
			int instantSuccess = masterDomain_->getInstantAvgAcc();
			int instantFailure = masterDomain_->getInstantAvgRej();
			bool IDNotFound = true;

			if ( bellCenter == 0 ){
				masterDomain_->setBellCenter( maxAcc );
				bellCenter = masterDomain_->getBellCenter(); 
				for ( const auto& p : configHolder ){							
					typeHold = p.getType();
					if (typeHold == 2){	
						if ( p.getSuccessFactor() == bellCenter ){
								vacanciesInBellCenter += 1;	
								totalRejectionsUnderBell += p.getFailureFactor(); 
						}
					}
				}
				avgRejectionsUnderBell = totalRejectionsUnderBell*1.0/vacanciesInBellCenter;
				masterDomain_->setAvgRejectionsUnderBell( avgRejectionsUnderBell );	
				masterDomain_->setInstantAvgAcc ( floor(avgSuccessInActiveRegion) );
				masterDomain_->setInstantAvgRej ( maxRej );
				instantSuccess = masterDomain_->getInstantAvgAcc();
				instantFailure = masterDomain_->getInstantAvgRej();
				if ( (avgRejectionsUnderBell >= instantFailure) && (bellCenter == floor(instantSuccess)) ){
					masterDomain_->setBellCenter ( maxAcc );
				}
				else if ( avgRejectionsUnderBell >= instantFailure ){
					masterDomain_->setBellCenter( bellCenter-1 );
				}
			}	
			else{
				while(IDNotFound){
					bellCenter = masterDomain_->getBellCenter(); 
					for ( const auto& p : configHolder ){							
						typeHold = p.getType();
						if (typeHold == 2){	
							if ( p.getSuccessFactor() == bellCenter ){
									vacanciesInBellCenter += 1;	
									totalRejectionsUnderBell += p.getFailureFactor(); 
									IDNotFound = false;
							}
						}
					}
					if (IDNotFound){
						bellCenter = bellCenter - 1;
						masterDomain_->setBellCenter ( bellCenter );
					}			
				}

				avgRejectionsUnderBell = floor(totalRejectionsUnderBell*1.0/vacanciesInBellCenter);	
				masterDomain_->setAvgRejectionsUnderBell( avgRejectionsUnderBell );
				if ( (avgRejectionsUnderBell >= instantFailure) && (bellCenter == floor(instantSuccess)) ){
					masterDomain_->setBellCenter ( maxAcc );
					masterDomain_->setInstantAvgAcc ( floor(avgSuccessInActiveRegion) );
					masterDomain_->setInstantAvgRej ( maxRej );
				}
				else if ( avgRejectionsUnderBell >= instantFailure ){
					masterDomain_->setBellCenter( bellCenter-1 );
				}
			}
			
			masterDomain_->setMaxSuccess( maxAcc ); 
			masterDomain_->setMaxFailure( maxRej ); 
			masterDomain_->setAvgSuccess(avgSuccessInActiveRegion ); 
		       	masterDomain_->setAvgFailure( avgFailureInActiveRegion );

		}
	}

         evaluateDecision(inJob);
         int jobRunning = 0, jobFinished = 0;
         int jobUnavailable = 0;

         for(const auto& j: jobQueue_){
            if(j.state == RUNNING)          jobRunning++;
            else if(j.state == FINISHED)    jobFinished++;
            else if(j.state == UNAVAILABLE) jobUnavailable++;
         }

         if(inJob.myResult.result==true){

            masterDomain_->updateDomain(inJob.myResult.particles,inJob.swapList);

            nDiscarded = inJob.discardDependencies(jobQueue_);

            for(decltype(inJob.swapList.size()) i=0; i<inJob.swapList.size(); i++){

               long particleID = inJob.swapList[i].ID;

       	       for(auto& job : jobQueue_){

       	          if( job.state == FINISHED || job.state == RUNNING ){ 
           	     if (inJob.swapList[i].type == TARGET )
           	        job.addedTargets.push_back(particleID);
           	     else
           	        job.deletedTargets.push_back(particleID);
       		  }

       	       }

            }

            totDeltaEpot += inJob.myResult.deltaEnergy;

         }

        updateStatCounters(inJob);
       	nFinished++;

	int eqStep = masterDomain_->getEquilibrationStep();

	if ( masterDomain_->getSamplingMode()==2 || masterDomain_->getSamplingMode()==4 || masterDomain_->getSamplingMode()==6 ){

		if ( masterDomain_->getBiasedBlock()==0 && masterDomain_->getLocalBlock()==0 ){
//cout << "The trial move performed at MC step " << nFinished << " was UNBIASED" << endl;
		}

		if ( masterDomain_->getBiasedBlock()==0 && masterDomain_->getLocalBlock()==1 ){
			biasedMovesPerformed = 0;
			localMovesPerformed +=1;
//cout << "The trial move performed at MC step " << nFinished << " was LOCAL" << endl;
//cout << "biasedMovesPerformed: " << biasedMovesPerformed << endl;
//cout << "localMovesPerformed: " << localMovesPerformed << endl;
		}

		if ( masterDomain_->getBiasedBlock()==1 && masterDomain_->getLocalBlock()==0 ){
			biasedMovesPerformed += 1;
			localMovesPerformed = 0;
//cout << "The trial move performed at MC step " << nFinished << " was BIASED" << endl;
//cout << "biasedMovesPerformed: " << biasedMovesPerformed << endl;
//cout << "localMovesPerformed: " << localMovesPerformed << endl;
		}

		if ( biasedMovesPerformed == masterDomain_->getStepsInBiasedBlock() ){
			masterDomain_->setBiasedBlock(0);
			masterDomain_->setLocalBlock(1);
//cout << "At MC step " << nFinished << ":" << endl; 
//cout << "LOCAL MODE IS SET UP FOR THE NEXT TRIAL MOVE" << endl;
		}

		if ( localMovesPerformed == masterDomain_->getStepsInLocalBlock() ||  nFinished== eqStep ){
			masterDomain_->setBiasedBlock(1);
			masterDomain_->setLocalBlock(0);
//cout << "At MC step " << nFinished << ":" << endl; 
//cout << "BIASED MODE IS SET UP FOR THE NEXT TRIAL MOVE" << endl;
		}

	}

	if ( masterDomain_->getSamplingMode()!=0 && masterDomain_->getSamplingMode()!=7 && masterDomain_->getSamplingMode()!=8 ){	
		masterDomain_->createSphereTarget(inJob.swapList[0].ID, sphereAroundSelectedSamplingSite); 
		sphereAroundSelectedSamplingSite.push_back( masterDomain_->getParticles()[inJob.swapList[0].ID] );	
		masterDomain_->increaseSuccessOrFailureFactor(sphereAroundSelectedSamplingSite, inJob.myResult.result);
	}
	
	long insertID = inJob.swapList[0].ID;
	Vector3d add_p = masterDomain_->getParticles()[insertID].getPosition();
	long deleteID = inJob.swapList[1].ID;
	Vector3d del_p = masterDomain_->getParticles()[deleteID].getPosition();

	       	 resultOut<<"   "<<(nFinished)<<"   "<<inJob.myResult.result<<"  "<<setprecision(8)<<inJob.myResult.deltaEnergy<<"   "<<setprecision(8)<< totDeltaEpot <<"   " 
		 << insertID<< "   "<< add_p.x<<"   " <<  add_p.y<<"   " <<  add_p.z << "   " << deleteID<<"   "<< del_p.x<<"   " <<  del_p.y<<"   " <<  del_p.z <<"   " 
		 << setprecision(5) <<(static_cast<float> (nAccepSwaps)/nFinished)*100<< "   " <<setprecision(5) <<(static_cast<float> (nRejeSwaps)/nFinished)*100        
		 << endl;
	
       	 string status = (inJob.myResult.result) ? "ACCEPTED" : "REJECTED";
       	 long       writeInterval= masterDomain_->getMCWriteInterval();

       	 if(nFinished==masterDomain_->getMCTotalSteps())

       	    cout<<"Wall time required to perform   "<<nFinished<< "  steps "<<"with time: "<<(clock()-manager_Start_time)/CLOCKS_PER_SEC <<"  (secs)"<<endl;

   	 if((nFinished%writeInterval == 0) && (nFinished>0)){

   	    cout <<"Flushing file information for Step No: "<< inJob.StepNO << endl;
   	    string interFile = "configuration_MC_"+to_string(nFinished)+".chkpt";
   	    masterDomain_->writeConfiguration(interFile);

         }

         masterDomain_->updateSampleFrequency(inJob.swapList);
         deleteJobQueueHead();
         return true;

       } 

       else if(inJob.state==UNAVAILABLE){

          masterDomain_->updateSampleFrequency(inJob.swapList);
    	  deleteJobQueueHead();
    	  return true;

       }

   }    

   if(inJob.exertedMove == SWAP || inJob.exertedMove == CLUSTER ){
      cerr<<"ERROR: STATE CONFLICT NOT SUPPORTED FOR SWAP TRIAL MOVE "<<endl;
      exit(1);
   }

   return false;

}


void Master::writeStatistics(void){

	if(getMoveType() == SINGLE){
        cout << " No of particle insertion     : " << nAddition << endl;
 	 	cout << " No of particle deletion      : " << nDeletion << endl;
 	 	cout << " No of accepted insertion     : " << naddAccepted << endl;
 	 	cout << " No of rejected insertion     : " << naddRejected << endl;
 	    cout << " No of accepted deletion      : " << ndelAccepted << endl;
 	    cout << " No of rejected deletion      : " << ndelRejected << endl;
	}
	else if(getMoveType() == SWAP){
		cout<<" No of attempted Swap moves : "<<nSwaps<<endl;
		cout<<" No of Accepted  Swap moves : "<<nAccepSwaps<<endl;
		cout<<" No of Rejected  Swap moves : "<<nRejeSwaps<<endl;
		cout<<" No performed biased moves  : "<<biasedMovesPerformed<<endl;
		cout<<" No performed local moves   : "<<localMovesPerformed<<endl;
	}
	cout << " No of encountered conflicts  : " << nConflicts << endl;
}

void Master::updateStatCounters(Job& inJob){

	for(decltype(inJob.swapList.size()) i=0;i<inJob.swapList.size();i++){
		if((inJob.exertedMove != SWAP) && (inJob.swapList[i].type ==TARGET)){
			//cout <<" Trial Move : INSERTION"<<endl;
			nAddition++;
			(inJob.myResult.result) ? naddAccepted++ : naddRejected++;
		}else if((inJob.exertedMove != SWAP) && (inJob.swapList[i].type ==SAMPLE)){
			//cout <<" Trial Move : DELETION"<<endl;
			nDeletion++;
			(inJob.myResult.result) ? ndelAccepted++ : ndelRejected++;
		}
		else if((inJob.exertedMove == SWAP) && (inJob.state != DISCARDED)){
	        //cout<<"Trial Move : SWAP "<<endl;
	        nSwaps++;
	        (inJob.myResult.result) ? nAccepSwaps++ : nRejeSwaps++;
	        return;
		}
	}
}

bool Master::updateJobIfStateConflict(Job& inJob){

	if(inJob.exertedMove == SWAP || inJob.exertedMove == CLUSTER){
		cerr<<"ERROR: UPDATE STATE CONFLICT NOT POSSIBLE"<<endl;
		exit(1);
	}

	long         nCurrent =0, nAssumed = 0;
	bool         stateConflict  = false;
	bool         isChanged      = false;
	double       stateRandom    = 0.0;
	vector<long> newList;

	std::uniform_real_distribution<double>   distribution(0.0,1.0);
	std::uniform_int_distribution<int>       intDistro;

    vector<long> addedList   = inJob.addedTargets;
	vector<long> deletedList = inJob.deletedTargets;

	// special case shall be handled even if the list is empty
	std::sort(addedList.begin(),addedList.end());
	std::sort(deletedList.begin(),deletedList.end());

	vector<long> addedListDiff;   // addedList   - deletedList
	vector<long> deletedListDiff; // deletedList - addedList

	// Perform Set operation to eliminate common entries
	std::set_difference(addedList.begin(), addedList.end(),
			deletedList.begin(), deletedList.end(), std::inserter(addedListDiff, addedListDiff.begin()));

    std::set_difference(deletedList.begin(), deletedList.end(),
    		addedList.begin(), addedList.end(),std::inserter(deletedListDiff, deletedListDiff.begin()));

    // Check 1 Addition/Deletion
	if(inJob.swapList[0].type == TARGET){

		// Addition of Target particles, so assumed and current Samples have to be considered
		nAssumed   = inJob.sites.nAssumedSamples;
		nCurrent   = masterDomain_->getSampleListSize();

		// Addition implies change from Sample to Target
		// to check whether the handled particle exist currently
		isChanged = (std::find(addedListDiff.begin(), addedListDiff.end(), inJob.swapList[0].ID)
		            != addedListDiff.end()) ? true : false;
	}else{

		// Deletion of Target particles, so assumed and current Targets have to be considered
        nAssumed  = inJob.sites.nAssmumedTargets;
        nCurrent  = masterDomain_->getTargetListSize();

        // Deletion implies change from Target to Sample
        // to check whether the handled particle exist currently
		isChanged = (std::find(deletedListDiff.begin(),deletedListDiff.end(), inJob.swapList[0].ID)
		             != deletedListDiff.end()) ? true : false;
	}

	if(nCurrent > nAssumed){        // Case 1
		 if(inJob.swapList[0].type == TARGET){
			  newList  = inJob.deletedTargets; // (targets turned into samples)
		 }else{
			  newList  = inJob.addedTargets;   // (samples turned into targets)
		 }
	}
	else if(nCurrent <= nAssumed){ // Case 2
		 if(inJob.swapList[0].type == TARGET){
			  newList  = inJob.addedTargets;   // (samples turned into targets)
		 }else{
			  newList  = inJob.deletedTargets; // (targets turned into samples)
		 }
	}

	double additionalProb = (static_cast<double> (newList.size())/nCurrent);
	stateRandom  = distribution(generatorRand_);

	if( (stateRandom <= additionalProb) || isChanged ){
		 stateConflict     =  true;
		 intDistro         =  std::uniform_int_distribution<int>(0,newList.size()-1);
         //int  index        =  intDistro(generatorRand_);
		 //long particleID   =  newList[index];
		 // Update Job
		 //inJob.particleID  =  particleID;
		 inJob.state       =  DISCARDED;
		 inJob.repeatMode  =  STATE;      // Mark the repeat mode as State conflict
	}
	else{
		 stateConflict = false;
	}
	return stateConflict;
}

void Master::processJobQueue() {
	while (isQueueHeadReady()) {
		processJob(jobQueue_.front());
	}
}

void Master::runAsMaster(void){

	//todo performance measurement
    manager_Start_time = clock();

    long       totalMCSteps = masterDomain_->getMCTotalSteps();
    bool       finished =false;
    MPI_Status  status;
    
    if(totalMCSteps>0 && nFinished>0){
    	if(totalMCSteps==nFinished){
    	    cerr<<"ERROR: mc_Steps and mc_nFinished are same"<<endl;
    	    exit(1);
    	}
    }

    resultOut<<"# StepNo    Res   deltaEpot   totEpot   VID   X   Y   Z   CID   X   Y   Z   accepRate   rejecRate"<<endl;

    // Perform this main loop as long as the required number of
    // MonteCarlo (MC) steps are completed.
    while((!finished)){

    	    // Check for free/Completed Slave process
    	    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
 	        int source    = status.MPI_SOURCE;
 	        //cout << "Master probed from Source"<< source<<" with tag "<<status.MPI_TAG<<endl;

 	        // worker process finished with results
 	        if (status.MPI_TAG == TAG_RUNJOB){
 	        	    fetchResult(source);

 			        processJobQueue();
 	 	            Job* jobPtr = getNextJob();
 	 	            if(jobPtr!=nullptr){
 	 	                  sendJobToWorker(*jobPtr,source);}
 	 	            else{
 				    	// terminate the slave process
 				    	MPI_Send(NULL, 0, MPI_LONG, source, TAG_TERMINATE, MPI_COMM_WORLD);
 				    	nSlaves_--;
 				    	cout <<" nFinished    : "<<nFinished<<endl;
 				    	cout <<" totalMCSteps : "<<totalMCSteps<<endl;
 				    	cout <<" nPerformed   : "<<nPerformed<<endl;
 				    	cout <<" nInvalidated : "<<nInvalidated<<endl;
 				    	cout <<"jobQueue_.size()"<<jobQueue_.size()<<endl;
 				    	cout <<" No of available workers :"<<nSlaves_;
 	 	            }
 	        }

 	        // Free worker process without results
 	        else if(status.MPI_TAG == TAG_FREE){
 	            long freeFlag;
 	        	MPI_Recv(&freeFlag, 1, MPI_LONG, source,TAG_FREE, MPI_COMM_WORLD, &status);
 	        	//cout <<" Master received info about Free processor "<< source <<endl;
	 	        Job* jobPtr = getNextJob();  //TODO replace if possible
	 	        if(jobPtr!=nullptr)
	 	             sendJobToWorker(*jobPtr,source);
			else{
			  //cout << "Processors being killed because termination problem ID: " << source << endl;  
			  //cout << "Their status is: " << status.MPI_TAG << endl;
			  //The prescribed number of MC steps has been completed. Nevertheless, some processors are still longing to do more.
			  //Those slaves must be terminated. 
			  MPI_Send(NULL, 0, MPI_LONG, source, TAG_TERMINATE, MPI_COMM_WORLD);
			  nSlaves_--;
			  cout <<" No of available workers :"<<nSlaves_;
			}
 	        }
		
  	      	if ((nFinished == totalMCSteps) && (jobQueue_.empty()) ){
	      	 	   finished=true;
	      	 	   cout << "============================================" << endl;
	      	 	   cout << " MonteCarlo Sampling Successful !! " << endl;
	      	 	   cout << " Total Required  MonteCarlo Steps       : " << nFinished  << endl;
	      	 	   cout << " Total Performed MonteCarlo Steps       : " << nPerformed << endl;
	      	 	   cout << " Total Invalidated MonteCarlo Steps       : " <<nInvalidated << endl;
	      	 	   writeStatistics();
	      	       cout << "============================================" << endl;
	      	}

		
    } // loop over masterSteps
    cout <<" Master Process Successfully terminated." << endl;

    resultOut.close(); // closing outfile connection

}
