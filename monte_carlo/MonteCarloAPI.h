/*
 * MonteCarloAPI.h
 *
 */

#ifndef MC___MONTECARLOAPI_H_
#define MC___MONTECARLOAPI_H_

#include "mpi.h"
#ifdef __cplusplus
extern "C" {
#endif

    // prepare MonteCarlo simulation setup from MD
    void initializeMonteCarlo(int* mdcpuDimension, int mdTotalTypes,int mdRealTypes,
    		int* mdRestriction,double* mdSimbox,double mdTemperature,
			double mdSphereRadius,double mdSphereWallThickness,int mdTotalMCSteps,
			int* mdPBC,double* mdZone,MPI_Comm mdCommunicator,int debug);

    // export required data to MonteCarlo
    void packConfigurationToMonteCarlo(long mdmcnParticles,long *mdmcNumbers,int *mdmcTypes,
		double  *mdmcMasses,double  *mdmcPositions , double *mdmcPotentialEnergies);

    // perform monteCarlo sampling and import updated data
    void performMonteCarlo(int mdPrank,long *mdnParticles,long **mdNumbers,int **mdTypes,
    		double **mdMasses,double **mdPositions);

    // import Potential energy to MD
    void importPotentialEnergy(double **mdPotentialEnergies);

    // clear all data structures created by MonteCarlo
    void cleanMonteCarlo();

    // clear all MonteCarlo simulation setup
    void shutDownMonteCarlo();

#ifdef __cplusplus
}
#endif

#endif /* MC___MONTECARLOAPI_H_ */
