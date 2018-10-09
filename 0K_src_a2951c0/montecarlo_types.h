/***********************************************************************************************
 *                                    montecarlo_types.h
 *
 *                      Data types definitions for Montecarlo_API interface
 *
 *  Created on: Sep 9, 2014
 /**********************************************************************************************/

#ifndef MONTECARLO_TYPES_H_
#define MONTECARLO_TYPES_H_

/* data types part */

long rand_ref;                                 // dummy variable

int     montecarlo_ncells;                     //  total no of cells
int montecarlo_celldim;                        //  local cell array dimension [per cpu including buffer]
int montecarlo_globdim;                        //  global cell array dimension

long    *montecarlo_atomtypes INIT(0);         //  atom type number
long    montecarlo_natoms INIT(0);             //  total number of real atoms
long    montecarlo_patoms INIT(0);             //  total number of placeholder atoms
long    montecarlo_tatoms_cpu INIT(0);         //  sum of real & placeholder atoms [per cpu]
long    *montecarlo_atomnumber INIT(0);        //  atom id
double  *montecarlo_atommass INIT(NULL);         //  atom mass
double  *montecarlo_positions INIT(NULL);        //  atom positions
double  *montecarlo_velocities INIT(NULL);       //  atom velocities
double  r_sample INIT(0.0);                      //  radius for sphere to be sampled

#endif /* MONTECARLO_TYPES_H_ */
