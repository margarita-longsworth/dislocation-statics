#include "imd.h"

void init_bopfox(){
	int i;
	double boxsize[9];
	int pbc[3];
	int cpuDim[3];
	int error = 0;

	boxsize[0] = box_x.x; boxsize[1] = box_y.x; boxsize[2] = box_z.x;
	boxsize[3] = box_x.y; boxsize[4] = box_y.y; boxsize[5] = box_z.y;
	boxsize[6] = box_x.z; boxsize[7] = box_y.z; boxsize[8] = box_z.z;

	pbc[0] = pbc_dirs.x; pbc[1] = pbc_dirs.y; pbc[2] = pbc_dirs.z;
	cpuDim[0] = cpu_dim.x; cpuDim[1] = cpu_dim.y; cpuDim[2] = cpu_dim.z;

	boplib_init_parameters_(boxsize, &bopfox_rcut, &bopfox_r2cut,
			&bopfox_rskin, &bopfox_rthickskin, cpuDim, pbc, &cpugrid, &error);

	for (i=0; i<ntypes; i++){
		size_t l = strlen(bopfox_typename[i]);
		boplib_set_type_(&i, bopfox_typename[i], &l);
	}

	if (error != 0){
		error("Calling bopfox returned an error");
	}
}

void exit_bopfox(){
	boplib_exit_();
}

void calc_forces_bopfox(int steps){
	if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();

	int error = 0;
	fix_cells();
	pack_bopfox();

#ifdef DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Starting Bopfox calculation\n");
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	boplib_energy_force_(&bopfox_nAtoms, bopfoxInputPositions,
			bopfoxInputTypes, bopfoxResultEnergy, bopfoxResultForce, &error);

	if (error != 0){
		error("Calling bopfox returned an error");
	}

#ifdef DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Bopfox calculation finished\n");
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	unpack_bopfox();
}

void pack_bopfox(void){
	int k, n, m, i;

	/* (re-)allocate data structures if necessary */
	bopfox_nAtoms = 0;
	for (k = 0; k < NCELLS; ++k)
		bopfox_nAtoms += CELLPTR(k)->n;

	if (bopfox_nAtomsMax < bopfox_nAtoms) {
		bopfox_nAtomsMax = (int)(1.1 * bopfox_nAtoms);
		free(bopfoxInputPositions);
		free(bopfoxInputTypes);
		free(bopfoxResultEnergy);
		free(bopfoxResultForce);
		bopfoxInputPositions = malloc(bopfox_nAtomsMax * sizeof(double) * DIM);
		bopfoxInputTypes = malloc(bopfox_nAtomsMax * sizeof(int));
		bopfoxResultEnergy = malloc(bopfox_nAtomsMax * sizeof(double));
		bopfoxResultForce = malloc(bopfox_nAtomsMax * sizeof(double) * DIM);
		if (bopfoxInputPositions == NULL || bopfoxInputTypes == NULL ||
				bopfoxResultEnergy == NULL || bopfoxResultForce == NULL) error("Could not allocate bopfox data");
	}

	/* clear forces and energies */
	clear_forces();

	/* collect data from cell array */
	n = 0;
	m = 0;
	for (k = 0; k < NCELLS; ++k) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; ++i) {
			bopfoxInputPositions[n++] = ORT(p,i,X);
			bopfoxInputPositions[n++] = ORT(p,i,Y);
			bopfoxInputPositions[n++] = ORT(p,i,Z);
			bopfoxInputTypes[m++] = SORTE(p,i);
		}
	}

	printf("Send %i atoms to bopfox cpu%i\n",bopfox_nAtoms, myid);
}

void unpack_bopfox(void){
	real e, pot2, pot1 = 0.0;
	int i, k;
	int n = 0;
	int m = 0;
	for (k = 0; k < NCELLS; ++k) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; ++i) {
			KRAFT(p,i,X) = bopfoxResultForce[n++];
			KRAFT(p,i,Y) = bopfoxResultForce[n++];
			KRAFT(p,i,Z) = bopfoxResultForce[n++];
			e = bopfoxResultEnergy[m++];
			POTENG(p,i) = e;
			pot1 += e;

			if (isnan(KRAFT(p,i,X)) || isnan(KRAFT(p,i,X) || isnan(KRAFT(p,i,X)))){
				printf("Force on atom %i is NAN\n", NUMMER(p,i));
			}
			if (isnan(POTENG(p,i))){
				printf("Warning: Potential energy on atom %i is NAN\n", NUMMER(p,i));
			}
			if (isinf(KRAFT(p,i,X)) || isinf(KRAFT(p,i,X) || isinf(KRAFT(p,i,X)))){
				printf("Force on atom %i is infinite\n", NUMMER(p,i));
			}
			if (isinf(POTENG(p,i))){
				printf("Warning: Potential energy on atom %i is infinite\n", NUMMER(p,i));
			}
			if (POTENG(p,i)>1000.){
				printf("Warning: Potential energy on atom %i is very large: %f\n", NUMMER(p,i),POTENG(p,i));
			}
		}
	}

	/* sum up potential energy */
#ifdef MPI
	MPI_Allreduce(&pot1, &pot2, 1, MPI_DOUBLE, MPI_SUM, cpugrid);
	tot_pot_energy += pot2;
#else
	tot_pot_energy += pot1;
#endif
}

void clear_forces(void){
	int k, i;
	tot_pot_energy = 0.0;
	virial = vir_xx = vir_yy = vir_xy = vir_zz = vir_yz = vir_zx = 0.0;
	for (k = 0; k < nallcells; k++) {
		cell *p = cell_array + k;
		for (i = 0; i < p->n; i++) {
			KRAFT(p,i,X) = 0.0;
			KRAFT(p,i,Y) = 0.0;
			KRAFT(p,i,Z) = 0.0;
#if defined (STRESS_TENS)
			PRESSTENS(p,i,xx) = 0.0;
			PRESSTENS(p,i,yy) = 0.0;
			PRESSTENS(p,i,xy) = 0.0;
			PRESSTENS(p,i,zz) = 0.0;
			PRESSTENS(p,i,yz) = 0.0;
			PRESSTENS(p,i,zx) = 0.0;
#endif
			POTENG(p,i) = 0.0;
		}
	}
}
