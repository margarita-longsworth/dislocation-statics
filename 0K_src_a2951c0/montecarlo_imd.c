#define MAIN
#include "imd.h"
#include "montecarlo_imd.h"

//Local routines
void insertAtoms(libimd_atom* atoms, int n);
#ifdef RELAX
int isRelaxed(void);
#endif

void LIB_PUBLIC libimd_init(char* parameterFile){
	printf("******Init IMD as library ******\n");

	read_parameters(parameterFile, 0);
	setup_potentials();

	if(box_from_header)
		error("Box size must be given! \"box_from_header\" must be 0");
	if(box_x.x == 0. || box_y.y == 0. || box_z.z == 0.)
		error("Box size is not initialized!");

	make_box();

	printf("******Init IMD as library finished ******\n");
}

int LIB_PUBLIC libimd_run(libimd_atom* atoms, int n){
	int i,j;
	int success = 1;

	imd_timer timer;
	imd_start_timer(&timer);

	insertAtoms(atoms,n);

	compute_energy_between_fixed_atoms = 1;
	have_valid_nbl = 0;
	calc_forces(0);
	compute_energy_between_fixed_atoms = 0;
	have_valid_nbl = 0;

#ifdef NVT
	if (ensemble == ENS_NVT)
		eta = 0.;
#endif
#ifdef GLOK
	if (glok_start <= steps_min) glok_start=steps_min;
	success = 0;
#endif

	maxwell(temperature);

	for (steps = 0; steps < steps_max; steps++) {
		calc_forces(steps);

#ifdef GLOK
		if (ensemble == ENS_GLOK) update_glok();
#endif

		move_atoms();

#ifdef NBLIST
		check_nblist();
#else
    	fix_cells();
#endif

#ifdef GLOK
#ifdef VIRTUAL_ATOMS
    	if(steps%update_step == 0){	//Only check relaxed state at timesteps
    								//where placeholders interactions are computed
#endif
    		int relaxed = isRelaxed();
    		if (ensemble == ENS_GLOK && relaxed) {
    			success = 1;
    			break;
    		}
#endif
    	}
	}


#ifdef NBLIST
	have_valid_nbl = 0;
#endif
	calc_forces(0);

	imd_stop_timer(&timer);
	printf("Runtime: %f, Steps: %i, fnorm: %e\n",timer.total, steps, SQRT( fnorm / nactive));

	n=0;
	for (i=0; i<ncells;++i){
		cell *c = CELLPTR(i);
		for (j = 0; j < c->n; j++) {
			atoms[n].x = ORT(c,j,X);
			atoms[n].y = ORT(c,j,Y);
			atoms[n].z = ORT(c,j,Z);
			atoms[n].num = NUMMER(c,j);
			atoms[n].epot = POTENG(c,j);
			atoms[n].mass = MASSE(c,j);
			atoms[n].sort = VSORTE(c,j);
			n++;
		}
	}
	return success;
}

void insertAtoms(libimd_atom* atoms, int n){
	int i;

	//delete all existing atoms
	//reset state variables
	int cells = cell_dim.x * cell_dim.y * cell_dim.z;
	for (i=0; i< cells ;++i){
		cell *c = cell_array+i;
		c->n=0;
	}

	nactive = 0;
	natoms  = 0;
#ifdef VIRTUAL_ATOMS
	nactive_placeholders = 0;
	patoms  = 0;
#endif
#ifdef NBLIST
	deallocate_nblist();
#endif

	// allocate input cell for one atom
	cell *incell = malloc(sizeof *incell);
	if (0 == incell) error("Cannot allocate input cell.");
	incell->n_max = 0;
	alloc_cell(incell, 1);

	for (i = 0; i < n; ++i) {

		ivektor cellc;
		cell* to;
		shortint s = atoms[i].sort;

		incell->ort[0] = (real)atoms[i].x;
		incell->ort[1] = (real)atoms[i].y;
		incell->ort[2] = (real)atoms[i].z;
		incell->nummer[0] = atoms[i].num;
		incell->masse[0]  = atoms[i].mass;
		incell->vsorte[0] = s;
		incell->sorte[0]  = s%ntypes;

		incell->impuls[0] = 0.0;
		incell->impuls[1] = 0.0;
		incell->impuls[2] = 0.0;

		incell->kraft[0] = 0.0;
		incell->kraft[1] = 0.0;
		incell->kraft[2] = 0.0;

		cellc = cell_coord(incell->ort[0], incell->ort[1], incell->ort[2]);
		cellc = local_cell_coord(cellc);
		to = PTR_VV(cell_array,cellc,cell_dim);
		incell->n = 1;
		INSERT_ATOM(to, incell, 0);

#ifdef VIRTUAL_ATOMS
		int is_real = real_atom[s];
		nactive += ((restrictions+s)->x+(restrictions+s)->y+(restrictions+s)->z) * is_real;
		nactive_placeholders += ((restrictions+s)->x+(restrictions+s)->y+(restrictions+s)->z) * (1-is_real);
		natoms  += is_real;
		patoms  += (1-is_real);
#else
		nactive += (restrictions+s)->x+(restrictions+s)->y+(restrictions+s)->z;
		natoms  ++;
#endif
	}

	alloc_cell(incell, 0);
	free(incell);
}


#ifdef RELAX
int isRelaxed(void) {
    is_relaxed = 0;
    if ((ensemble==ENS_MIK) || (ensemble==ENS_GLOK)) {

        int stop = 0;
        real fnorm2, ekin, epot, delta_epot;

        fnorm2 = SQRT( fnorm / nactive );
        ekin   = 2 * tot_kin_energy / nactive;
        epot   = tot_pot_energy / natoms;
        delta_epot = old_epot - epot;
        if (delta_epot < 0) delta_epot = -delta_epot;

        if ((ekin  <  ekin_threshold) || (fnorm2 < fnorm_threshold) ||
            (delta_epot < delta_epot_threshold)) return 1;

        old_epot = epot;
    }
    return 0;
}
#endif
