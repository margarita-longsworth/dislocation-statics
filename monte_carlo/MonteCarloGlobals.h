/*
 * MonteCarloGlobals.h
 *
 */

#ifndef MC___MONTECARLOGLOBALS_H_
#define MC___MONTECARLOGLOBALS_H_

#include "Domain.h"

using namespace std;

IntVector3d mccpuDimension;
int mcTotalTypes, mcRealTypes;
int debugFlag;

// Standard Template Library

vector<int>     mcTypes;
vector<long>    mcNumbers;
vector<double>  mcMasses;
vector<double>  mcPositions;
vector<double>  mcVelocities;
vector<double>  mcPotentialEnergies;
vector<int>	mcAcc;
vector<int>	mcRej;
vector<int>	mcNewAcc;
vector<int>	mcNewRej;

Vector3d mcSimboxX,mcSimboxY,mcSimboxZ;

IntVector3d* mcRestriction = nullptr;
IntVector3d mcPBC;
Vector3d mcZoneSize;

double mcTemperature;
double mcSphereRadius;
double mcSphereWallThickness;
double mcCohesive;
long   mcFlushInterval;
Vector3d mcCellSize;

long   mcnParticles;
int    mcnSteps;
int    mcPrank, mcnProcess;
int    mcSeed;
int    mcMasterFlag; // set to one for Master Slave approach !!
int    mcTrialMoveType;

long mcNConflicts,mcNFinished,mcNSwaps,mcNAccepSwaps,mcNRejeSwaps,mcNAddition,mcNDeletion;
long mcNaddAccepted,mcNdelAccepted,mcNaddRejected, mcNdelRejected, mcBiasedMovesPerformed, mcLocalMovesPerformed;
double mcTotDeltaEpot;

string    mcSphereFile;
string    mcInputFile;
string    mcOutputFile, mcStatFile;

int mcSimulationID;
Domain* ptrDomain = nullptr;

MPI_Comm mcCommunicator;

//Biased sampling

int mcSamplingMode;
long mcStepsInLocalBlock;
long mcStepsInBiasedBlock;
double mcSphereRadiusTarget;
double mcLocalMovesRadius;
double mcCylinderRadius;
int mcCoveringTimes;

#endif /* MC___MONTECARLOGLOBALS_H_ */
