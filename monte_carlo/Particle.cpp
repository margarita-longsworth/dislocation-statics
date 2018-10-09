/*
 * Particle.c++
 *
 */
#include "Particle.h"
using namespace std;

Particle::Particle(void){
	type_       = 0;
	number_     = 0;
	mcID_       = 0;
	mass_       = 0.0;
	epot_       = 0.0;
	position_   = {0.0,0.0,0.0};
	sampleFreq_ = 0;

	successFactor_ = 0;
	failureFactor_ = 0;
}

Particle::~Particle(void){

}

void Particle::setType (int inType)        { type_   = inType;   }
void Particle::setNumber (long inNumber)   { number_ = inNumber; }
void Particle::setmcID (long mcID)         { mcID_   = mcID;     }
void Particle::setMass (double inMass)     { mass_   = inMass;   }
void Particle::setEpot (double inEpot)     { epot_   = inEpot;   }
void Particle::setPosition (double inPosx,double inPosy,double inPosz){
		   position_.x = inPosx;
		   position_.y = inPosy;
		   position_.z = inPosz;
}
void Particle::updateFrequency(int frequency){
	     sampleFreq_ += frequency;
}

int      Particle::getType()     const   { return type_;      }
long     Particle::getNumber()   const  { return number_;    }
long     Particle::getmcID()     const  { return mcID_;      }
double   Particle::getMass()     const  { return mass_;      }
double   Particle::getEpot()     const  { return epot_;      }
Vector3d Particle::getPosition() const  { return position_;  }

int      Particle::getFrequency()   const  { return sampleFreq_;}

void Particle::increaseSuccessFactor() { successFactor_ += 1;}
void Particle::increaseFailureFactor() { failureFactor_ += 1;}

long Particle::getSuccessFactor() const { return successFactor_; }
long Particle::getFailureFactor() const { return failureFactor_; }

void Particle::setSuccessFactor( int successFactor ) { successFactor_ = successFactor; }
void Particle::setFailureFactor( int failureFactor ) { failureFactor_ = failureFactor; }






