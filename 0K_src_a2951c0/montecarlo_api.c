/*****************************************************************************************************
 *
 *
 *                                 Generic Montecarlo API
 *  specifications to be included
 *  Created on: Sep 4, 2014
 *      Author: ganeshfx
 /****************************************************************************************************/

/*  MC routine */

#include<stdio.h>
#include<stdlib.h>


void do_montecarlo(int myid,long *montecarlo_tatoms_cpu,long *montecarlo_atomnumber,long *montecarlo_atomtypes,
		double  *montecarlo_atommass,double  *montecarlo_positions){

	int i,j, pid;
	long rand_no, tatoms;


	tatoms = *montecarlo_tatoms_cpu;
	pid = myid; // process id



	printf(" ------------------------------------------------- \n ");
	printf(" From process id: %d --- Hello from Montecarlo !!  \n",pid);


	/* selecting a random no */

	srand(time(NULL));
	rand_no=(long) rand()%tatoms;

	printf(" length of data array : %ld \n ",tatoms);
	printf(" generated random number : %ld \n", rand_no);

	while(montecarlo_atomtypes[rand_no]!=2){
		printf(" Chosen particle is not a placeholder and \n ");
		rand_no=(long) rand()%tatoms;
		printf(" process id: %d, atom type : %ld and atom id : %ld \n ",pid,montecarlo_atomtypes[rand_no],montecarlo_atomnumber[rand_no]);
	}
	printf(" ------------------------------------------------- \n ");



	/*  switching placeholder into real atom  (for FeC system, into C) */
	montecarlo_atomtypes[rand_no]= 1; // later include option in parameter file

	/* post check */

	printf(" corresponding atomic particle id from data: %ld \n", montecarlo_atomnumber[rand_no]);
	printf(" corresponding atomic type  from data: %ld \n", montecarlo_atomtypes[rand_no]);
	printf(" corresponding atomic mass from data: %f \n", montecarlo_atommass[rand_no]);
	printf(" corresponding atomic position_X from data: %f \n", montecarlo_positions[rand_no*3]);
	printf(" corresponding atomic position_Y from data: %f \n", montecarlo_positions[rand_no*3+1]);
	printf(" corresponding atomic position_Z from data: %f \n", montecarlo_positions[rand_no*3+2]);

}

void export_config(int myid,long *montecarlo_tatoms_cpu,long *montecarlo_atomnumber,long *montecarlo_atomtypes,
		double  *montecarlo_atommass,double  *montecarlo_positions){

	int pid=myid;
	long tot_particles = *montecarlo_tatoms_cpu;

	printf(" configuration is being packed up from Montecarlo \n");
	printf(" total particles count : %ld from process id: %d \n ",tot_particles,pid);
	printf(" End of packing form Montecarlo \n ");
}
