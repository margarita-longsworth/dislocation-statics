/*
 * Particle.h
 *
 */

#ifndef MC___PARTICLE_H_
#define MC___PARTICLE_H_

#include "Vector3d.h"

class Particle{


    public:
	   Particle();
       ~Particle();
       void setType (int inType);
	   void setNumber (long inNumber);
	   void setmcID (long mcID);
       void setMass (double inMass);
	   void setPosition (double inPosx,double inPosy,double inPosz);
	   void setEpot (double inEpot);
	   void updateFrequency (int frequency);

	   int          getType() const;
	   long         getNumber() const;
	   long         getmcID() const;
	   double       getMass() const;
	   double       getEpot() const;
	   Vector3d	getPosition() const;
	   int		getFrequency() const;

           void		increaseSuccessFactor();
           void		increaseFailureFactor();
	   void 	setSuccessFactor( int successFactor);
           void 	setFailureFactor( int failureFactor);
           long		getSuccessFactor() const;
           long		getFailureFactor() const;

    private:
	   int          type_;
	   int          sampleFreq_;
	   long         number_;
	   long         mcID_;
	   double       mass_;
	   double       epot_;
	   Vector3d	position_;

           long		successFactor_;
           long		failureFactor_;

};

#endif /* MC___PARTICLE_H_ */
