/*
 * imd_montecarlo.h
 *
 *  Header file for Montecarlo interface
 *  Created on: Sep 3, 2014
 *      Author: ganeshfx
 */

#include "montecarlo_types.h"

#ifndef IMD_MONTECARLO_H_
#define IMD_MONTECARLO_H_

/*  function prototypes */

void packconfig_to_montecarlo(void);

void do_montecarlo(int myid,long *montecarlo_tatoms_cpu,long *montecarlo_atomnumber,long *montecarlo_atomtypes,
		double  *montecarlo_atommass,double  *montecarlo_positions);

void export_config(int myid,long *montecarlo_tatoms_cpu,long *montecarlo_atomnumber,long *montecarlo_atomtypes,
		double  *montecarlo_atommass,double  *montecarlo_positions);

void unpackconfig_from_montecarlo(void);

#endif /* IMD_MONTECARLO_H_ */
