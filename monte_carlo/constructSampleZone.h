/*
 * constructSampleZone.h
 *
 */

#ifndef MC___CONSTRUCTSAMPLEZONE_H_
#define MC___CONSTRUCTSAMPLEZONE_H_

#include <random>
#include "Domain.h"

using namespace std;

// Abstract base class
class constructSampleZone{

  public:
    virtual long estimateSamplingSites(int speciesType,Domain* ptrDomain) = 0;
    virtual Particle findRandomParticle(int nSites,int particleType,Domain* ptrDomain) = 0;
	virtual ~constructSampleZone();
	Vector3d getZoneSize();
	Vector3d getRandomPosition();
  protected:
    constructSampleZone(Vector3d zoneSize,int seed);

  private:
    Vector3d      zoneSize_;
	std::mt19937  generatorSample_;
    std::uniform_real_distribution<double>  distributionXSample_;
    std::uniform_real_distribution<double>  distributionYSample_;
    std::uniform_real_distribution<double>  distributionZSample_;

};

// Concrete Class
class randomSampleZone: public constructSampleZone{

  public:
	randomSampleZone(Vector3d zoneSize,int seed):
		constructSampleZone(zoneSize,seed){};
	virtual long estimateSamplingSites(int speciesType,Domain* ptrDomain);
	virtual Particle findRandomParticle(int nSites,int particleType,Domain* ptrDomain);
	virtual ~randomSampleZone();

  private:
	Vector3d      centreOfMass_;
	vector<long>  sitesAtCell_;
//	std::mt19937  generatorParticle_;
	std::uniform_int_distribution<int> 	distributionParticle_;

};

#endif /* MC___CONSTRUCTSAMPLEZONE_H_ */
