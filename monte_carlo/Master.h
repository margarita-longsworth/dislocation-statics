/*
 * Master.h
 *
 */

#ifndef MASTER_H_
#define MASTER_H_

#include "TrialMove.h"
#include "acceptanceCheck.h"
#include "Domain.h"
#include "Job.h"
#include "algorithm"
#include <math.h>
#include <time.h>


class Job;

class Master{

     public:

	     Master(int nTrials,int seed,long mc_nConflicts,long mc_nFinished,
	    		long mc_nSwaps,long mc_nAccepSwaps,long mc_nRejeSwaps,long mc_nAddition,long mc_nDeletion,
				long mc_naddAccepted, long mc_ndelAccepted, long mc_naddRejected, long mc_ndelRejected,
				int nSlaves,double mc_totDeltaEpot,long mc_biasedMovesPerformed,long mc_localMovesPerformed,Domain* globDomain,string statFileName);

	    ~Master();
	     TrialMove* getTrialMove();


	     int getNTrials(void){return nTrials_;}
	     long getNConflicts(void){return nConflicts;}
	     long getNFinished(void){return nFinished;}
	     long getNSwaps(void){return nSwaps;}
	     long getNAccepSwaps(void){return nAccepSwaps;}
	     long getNRejeSwaps(void){return nRejeSwaps;}
	     double getTotDeltaEpot(void){return totDeltaEpot;}
	     long getBiasedMovesPerformed(void){return biasedMovesPerformed;}
	     long getLocalMovesPerformed(void){return localMovesPerformed;}

	     Move       getMoveType(void);
	     void              runAsMaster();
	     std::mt19937       generatorRand_; // TODO get method

     private:
	     Domain*            masterDomain_; //From master constructor in MasterSlavePostTest
	     int  TAG_RUNJOB, TAG_TERMINATE, TAG_INITIALIZE, TAG_FREE;
	     long nAddition,  nDeletion, naddAccepted, ndelAccepted, naddRejected, ndelRejected; //From master constructor in MasterSlavePostTest
	     long nSwaps,nAccepSwaps,nRejeSwaps; //From master constructor in MasterSlavePostTest
	     long nConflicts, nFinished, nAssigned,startStep, nPerformed, nInvalidated, biasedMovesPerformed, localMovesPerformed; //From master constructor in MasterSlavePostTest
	     int  nSlaves_; //From master constructor in MasterSlavePostTest
	     std::deque<Job>    jobQueue_;
             addDeleteParticle   trialMove; //trialMove gets assigned vals. calling cstr addDeleteParticle
             Move                move;
             ofstream            resultOut;
             double              totDeltaEpot;
	     int		 nTrials_;

         //todo different timer variables
             double              manager_Start_time;

	     void                processJobQueue();
             bool                isQueueHeadReady();
             Job*                findDiscardedJob();
             Job*                findDiscardedJobAndInvalidate(); // special method for modified jobQueue processing approach
             bool                updateJobIfStateConflict(Job& inJobStatPtr);
             Job*                getNextJob(void);
             void                sendJobToWorker(Job& inJob,int processID);
             void                deleteJobQueueHead();
             void                evaluateDecision(Job& inJob);
             void                fetchResult(int slaveRank);
             bool                processJob(Job& inJob);
             void                writeStatistics();
             void                updateStatCounters(Job& inJob);
};


#endif /* MASTER_H_ */
