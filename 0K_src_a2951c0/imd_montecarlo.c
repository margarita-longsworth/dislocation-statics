/***************************************************************************************************
*
* imd_montecarlo.c -- perform Montecarlo simulation for the given configuration
*
***************************************************************************************************
*
***************************************************************************************************
*
// * $version:1.0 $
// * $Date   :31/07/2014 $
// *
// **************************************************************************************************/

#include "imd.h"
#include "/home/users/ganeshfx/git/MC++/MonteCarloAPI.h"

// interface(MD-MC) parameters

long       md_mc_tatoms_cpu;   // total particles per cpu
int        md_mc_pid;          // processor id

// MD related parameters
int        md_ncells;
int        md_pid;            // processor id
long      md_tatoms_cpu;      // total particles per cpu
int*      md_atomtypes;
long*     md_atomnumber;
double*   md_atommass;
double*   md_positions;
double*   md_epot;

// check lines
int*      md_restriction;
int*      md_cpu_dim;         //  CPU dimension for MC
double*   md_simbox;          //  Simbox physical dimensions
int       md_tot_types;       //   no of total element types (CHECK for correctness??)
int       md_real_types;      //   no of real element types

void mc_blackbox(){


	// check globals.h

	md_tot_types      = vtypes;                                  // no of total element types (CHECK for correctness??)
	md_real_types     = ntypes;                                  // no of real element types

	md_restriction        = malloc(md_tot_types*3 * sizeof *md_restriction);
	md_cpu_dim            = malloc(3 * sizeof *md_cpu_dim);
	md_simbox             = malloc(9 * sizeof *md_simbox);

	int ind;
	int count = 0;

	for( ind=0; ind<md_tot_types; ind++){
		md_restriction[count++]  = (restrictions+ind)->x;
		md_restriction[count++]  = (restrictions+ind)->y;
		md_restriction[count++]  = (restrictions+ind)->z;
	}

	md_simbox[0] = box_x.x;
	md_simbox[1] = box_x.y;
	md_simbox[2] = box_x.z;

	md_simbox[3] = box_y.x;
	md_simbox[4] = box_y.y;
	md_simbox[5] = box_y.z;

	md_simbox[6] = box_z.x;
	md_simbox[7] = box_z.y;
	md_simbox[8] = box_z.z;

	md_cpu_dim[0] = cpu_dim.x;
	md_cpu_dim[1] = cpu_dim.y;
	md_cpu_dim[2] = cpu_dim.z;


	// later to be included via param file

	double md_temperature = 15;          // system temperature
	//double md_temperature = temperature; // system temperature
	double rsweep        = 8.000000;          // sphere radius
	double wall_thick    = 2.000000;          // sphere wall thickness
	int nMCsteps          = 10;          // no of MC steps
	MPI_Comm comm         =  cpugrid;
	int debug_Flag        = 1;

	int mdPBC[3] = {pbc_dirs.x,pbc_dirs.y,pbc_dirs.z};
	double mdZone[3] = {40.0,40.0,40.0};

	// calling IMD method
	fix_cells();

    // local parameters
	int i,j,k,m,n;

	md_mc_pid = myid;
	md_ncells = ncells;

	/* loop over cell array*/
    // get total atoms in CPU
	for(i=0; i<md_ncells ; ++i ){
		md_mc_tatoms_cpu+=CELLPTR(i)->n;
	}

	printf("-------------------------------------------------\n");
	printf("              IMD :some sanity check                   \n");
	printf("     Total Particles : %ld \n",md_mc_tatoms_cpu);
	printf("     Total Cells    : %d  \n",md_ncells);
	printf("     Process ID     : %d  \n",md_mc_pid);

	printf("-------------------------------------------------\n");

	/* critical declaration -1 */

	md_tatoms_cpu = md_mc_tatoms_cpu;
	md_pid    = myid;

	/* Memory allocations */

	int *md_mc_atomtypes     = malloc(md_mc_tatoms_cpu * sizeof *md_mc_atomtypes);
	long *md_mc_atomnumber    = malloc(md_mc_tatoms_cpu * sizeof *md_mc_atomnumber);
	double *md_mc_atommass    = malloc(md_mc_tatoms_cpu * sizeof *md_mc_atommass);
	double *md_mc_positions   = malloc(md_mc_tatoms_cpu * DIM * sizeof *md_mc_positions);
	double *md_mc_velocities  = malloc(md_mc_tatoms_cpu * DIM * sizeof *md_mc_velocities); // check if really necessary??
    double *md_mc_epot =       malloc(md_mc_tatoms_cpu * sizeof *md_mc_epot);

	/*  critical declaration 2 */
	md_atomtypes = NULL;
	md_atomnumber= NULL;
	md_atommass =  NULL;
	md_positions = NULL;
	md_epot      = NULL;

	/* Null pointer check */

	if (md_mc_atomtypes == NULL || md_mc_atomnumber == NULL || \
		md_mc_atommass== NULL || md_mc_positions == NULL ||\
		md_mc_velocities == NULL || md_mc_epot == NULL) \
					error("Could not allocate memory for montecarlo blackbox ");


    /* Phase -1 : packing data to Montecarlo */

	m=0;
	n=0;
	k=0;

    for(j=0; j<md_ncells; j++){
    	cell *p=CELLPTR(j);
    	for(i=0; i < p->n; i++ ){

    	   md_mc_atomtypes[m]      =  SORTE(p,i);
    	   md_mc_atomnumber[m]     =  NUMMER(p,i);
    	   md_mc_atommass[m]       =  (double)MASSE(p,i);

    	   md_mc_positions[k++]    =  (double)ORT(p,i,X);
    	   md_mc_positions[k++]    =  (double)ORT(p,i,Y);
    	   md_mc_positions[k++]    =  (double)ORT(p,i,Z);

    	   md_mc_velocities[n++]   =  (double)((IMPULS(p,i,X) / MASSE(p,i)));
    	   md_mc_velocities[n++]   =  (double)((IMPULS(p,i,Y) / MASSE(p,i)));
    	   md_mc_velocities[n++]   =  (double)((IMPULS(p,i,Z) / MASSE(p,i)));

    	   md_mc_epot[m]           =  POTENG(p,i);

		   m++;
        }
    }

    // clear all contents

    for(j=0; j<md_ncells; j++){
    	cell *p=CELLPTR(j);
    	for(i=0; i < p->n; i++ ){

    	   SORTE(p,i)  = 0;
    	   NUMMER(p,i) = 0;
    	   MASSE(p,i)  = 0.0;

    	   ORT(p,i,X) = 0.0;
    	   ORT(p,i,Y) = 0.0;
    	   ORT(p,i,Z) = 0.0;

    	   IMPULS(p,i,X) = 0.0;
    	   IMPULS(p,i,Y) = 0.0;
    	   IMPULS(p,i,Z) = 0.0;

    	   POTENG(p,i)  = 0.0;

        }
    }

    // Empty all cells //

    for (k=0; k<nallcells; k++) cell_array[k].n = 0;

    // initializing MonteCarlo
    initializeMonteCarlo(md_cpu_dim,md_tot_types,md_real_types,md_restriction,md_simbox,
    		md_temperature,rsweep,wall_thick,nMCsteps,mdPBC,mdZone,comm,debug_Flag);


 	// invoking Pack method to Montecarlo (MC-Interface) //
    packConfigurationToMonteCarlo(md_mc_tatoms_cpu,md_mc_atomnumber,md_mc_atomtypes,md_mc_atommass,md_mc_positions,md_mc_epot);

    free(md_mc_atomnumber);
    free(md_mc_atomtypes);
    free(md_mc_atommass);
    free(md_mc_positions);
    free(md_mc_velocities);

    long*    md_atomnumber;
    int*     md_atomtypes;
    double*  md_atommass;
    double*  md_positions;

    // Perform Montecarlo simulation (mc-Interface) //
    performMonteCarlo(md_pid,&md_tatoms_cpu,&md_atomnumber,&md_atomtypes,&md_atommass,&md_positions);

     // get values back
    long total_particle = *(&md_tatoms_cpu);

    // Unpacking data from Montecarlo Interface
    long i_up, j_up,m_up,k_up,s_up;   //_up -local variable for unpacking
    m_up=0,k_up=0;

    int to_cpu, count_atom;
    long nparticles, tmp;
    long addnumber = 0;

     // IMD data types

     minicell* to;
     ivektor cellc;
     real pos_x,pos_y,pos_z;
     cell* input;

     // allocate a new cell //
     input =  (cell *) malloc(sizeof(cell));
     if (0==input) error("Cannot allocate input cell in unpack config.");
     input->n_max=0;

     alloc_cell(input, 1); // allocate memory for new cell

     // initialization of parameters

     natoms  = 0;
     patoms  = 0;
     nactive = 0;
     nactive_placeholders = 0;

     for (i_up=0; i_up<ntypes; i_up++)
    	 num_sort[i_up]=0;

     for (i_up=0; i_up<vtypes; i_up++)
    	 num_vsort[i_up]=0;

     // loop over size of local data array - identical to loop over each line in file

     for ( j_up=0;j_up< total_particle; j_up++){

 	// updating values

     	input->n = 1;
         s_up = md_atomtypes[m_up];
         NUMMER(input,0) = md_atomnumber[m_up];
         SORTE (input,0) = MOD(s_up,ntypes);

         if (s_up>=vtypes) error("atom type must not exceed total_types");

         VSORTE(input,0) = s_up;
         MASSE (input,0) = md_atommass[m_up];

         pos_x = md_positions[k_up++];
         pos_y = md_positions[k_up++];
         pos_z = md_positions[k_up++];

         ORT(input,0,X) = pos_x;
         ORT(input,0,Y) = pos_y;

 #ifndef TWOD
         ORT(input,0,Z) = pos_z;
 #endif

         IMPULS(input,0,X) = 0;
         IMPULS(input,0,Y) = 0;

 #ifndef TWOD
         IMPULS(input,0,Z) = 0;
 #endif

         KRAFT(input,0,X) = 0;
         KRAFT(input,0,Y) = 0;
 #ifndef TWOD
         KRAFT(input,0,Z) = 0;
 #endif

         cellc = cell_coord(pos_x,pos_y,pos_z); // give the global cell coordinate

         to_cpu = cpu_coord(cellc);
         count_atom = 0;

         if(*(&md_pid) != myid) error(" process mismatch during unpacking ");

         cellc = local_cell_coord(cellc);
         to = PTR_VV(cell_array,cellc,cell_dim);

         MOVE_ATOM(to,input,0);
         count_atom=1;

         m_up++;

         if(count_atom){

 #ifdef VIRTUAL_ATOMS

                int q_isReal= real_atom[VSORTE(input,0)]; // 1-real ; 0-virtual
                natoms += q_isReal;
                patoms += 1-q_isReal;

 #else
     	       natoms++;
 #endif

 #ifdef VIRTUAL_ATOMS
                nactive += (long) (restrictions+s_up)->x * q_isReal ;
                nactive += (long) (restrictions+s_up)->y * q_isReal;

                nactive_placeholders +=(long) (restrictions+s_up)->x - q_isReal ;
                nactive_placeholders +=(long) (restrictions+s_up)->y - q_isReal;
 #ifndef TWOD
                nactive += (long) (restrictions+s_up)->z * q_isReal;
                nactive_placeholders +=(long) (restrictions+s_up)->z - q_isReal;
 #endif

 #else
                nactive += (long) (restrictions+s_up)->x;
                nactive += (long) (restrictions+s_up)->y;
 #ifndef TWOD
                nactive += (long) (restrictions+s_up)->z;
 #endif
 #endif
                num_sort [ SORTE(input,0)]++;
                num_vsort[VSORTE(input,0)]++;
     	    }


     } // Loop over all particles in CPU

     // calling Montecarlo interface - clear memory created during montecarlo routine
     cleanMonteCarlo();

     // shutdown MC
     shutDownMonteCarlo();

     	  /* Add the number of atoms read (and kept) by each CPU */
 #ifdef MPI
     MPI_Allreduce( &natoms,  &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     natoms = tmp;

     MPI_Allreduce( &nactive, &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     nactive = tmp;

     for (i_up=0; i_up<ntypes; i_up++) {
     	 MPI_Allreduce( &num_sort[i_up], &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     	 num_sort[i_up] = tmp;
     }
     for (i_up=0; i_up<vtypes; i_up++) {
     	 MPI_Allreduce( &num_vsort[i_up], &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     	 num_vsort[i_up] = tmp;
     }

 #ifdef VIRTUAL_ATOMS
     MPI_Allreduce( &nactive_placeholders, &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     nactive_placeholders = tmp;
     MPI_Allreduce( &patoms, &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
     patoms = tmp;
 #endif

 #endif // MPI


}// End of mc_blackbox
