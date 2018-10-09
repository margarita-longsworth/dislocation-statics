/*
 * Job.h
 *
 */

#ifndef JOB_H_
#define JOB_H_

#include <vector>
#include "UtilityMethods.h"
#include "Master.h"
using std::vector;
using std::pair;

class Master;

class Job {
public:
	Job();
	Job(long stepNo, Master& master,Domain& domain,long currentMCStep);
	virtual ~Job();
	void     modifyJobForSpatialConflict(Master& master,Domain& domain, long currentMCStep);
	bool     testForStateConflict(Master& master, Domain& dom);
	void     modifyJobForStateConflict(Master& master,Domain& domain);
	int      discardDependencies(deque<Job>& jobQueue);
	void     createDependencies(deque<Job>& jobQueue,Domain& domain);
	//todo   this method might not be available in future versions (if Seg fault arises)
	void     clearInvalidatedJobData(void);

public:       // TODO to be private
		     long               StepNO;
		     Sites               sites;
		     Move                exertedMove;
		     vector<Swap>        swapList;
       	             vector<long>        addedTargets;
		     vector<long>        deletedTargets;
	 	     JobResult           myResult;
	 	     State               state;
	 	     RepeatMode          repeatMode;
	 	     long                nParticles;
	 	     int                 workerID;
	 	     int                 getVersionNo();

private:
	 	    int                 versionCounter;
	 	    vector<pair<long, int> >  dependancyList;
			void reInitJob(int versionNo,Master& master,Domain& domain,long currentMCStep);
	 	    void clearLists();
};

#endif /* JOB_H_ */
