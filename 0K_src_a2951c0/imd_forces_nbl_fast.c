/******************************************************************************
 *
 * IMD -- The ITAP Molecular Dynamics Program
 *
 * Copyright 1996-2011 Institute for Theoretical and Applied Physics,
 * University of Stuttgart, D-70550 Stuttgart
 *
 ******************************************************************************/

/******************************************************************************
 *
 * imd_forces_nbl_fast.c -- force loop with neighbor lists
 * optimized for EAM potentials.
 * During computing the potential energy and the electron density,
 * the gradient of each (!) pair interaction is cached in an additional array
 * Later, the forces of both parts are summed together.
 * This has the advantage that electron density and its gradient can be
 * computed at the same time and the force summation loop needs only to be
 * executed once.
 *
 ******************************************************************************/


#include "imd.h"
#include "potaccess.h"

#if   defined(FOURPOINT)
#define   PAIR_INT_FAST   PAIR_INT4_FAST
#elif defined(SPLINE)
#define   PAIR_INT_FAST   PAIR_INT_SP_FAST
#else
#define   PAIR_INT_FAST   PAIR_INT3_FAST
#endif


#define PAIR_INT_SP_FAST(pot, grd, pt, col, r2)					             \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  step  = (pt).step[col];                                                    \
  r2a   = r2 * istep;                                                        \
  k     = POS_TRUNC(r2a);                                                    \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k    += (col)*(pt).maxsteps;                                               \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k++;                                                                       \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1.;                                                        \
  b2    = b * b - 1.;                                                        \
  st6   = step / 6.;                                                         \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = a * p1 + b * p2 + (a * a2 * d21 + b * b2 * d22) * st6 * step;       \
  grd = 2.*((p2 - p1) * istep + ((3.*b2 + 2.) * d22 - (3.*a2 + 2.) * d21) * st6);\
}

#define PAIR_INT4_FAST(pot, grd, pt, col, r2)                     			 \
{                                                                            \
  real r2a, istep, chi, *t;									                 \
  real fac[4], dfac[4];	      												 \
  int  k;                                                                    \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation */                                        \
  fac[0] = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                         \
  fac[1] =        0.5 * (chi*chi-1.0) * (chi-2.0);                           \
  fac[2] =       -0.5 * chi * (chi+1.0) * (chi-2.0);                         \
  fac[3] =  (1.0/6.0) * chi * (chi*chi-1.0);                                 \
                                                                             \
  /* factors for the interpolation of the derivative */                      \
  dfac[0] = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                            \
  dfac[1] =        0.5 * ((3.0*chi-4.0)*chi-1.0);                            \
  dfac[2] =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                            \
  dfac[3] =    1.0/6.0 * (3.0*chi*chi-1.0);                                  \
                                                                             \
  /* intermediate values */                                                  \
  t = PTR_2D((pt).table, (col), k-1, 0, (pt).maxsteps);                      \
  /* potential energy */                                                     \
  pot = fac[0]*t[0] + fac[1]*t[1] + fac[2]*t[2] + fac[3]*t[3];               \
  /* twice the derivative */                                                 \
  grd = 2. * istep * (dfac[0]*t[0] + dfac[1]*t[1] + dfac[2]*t[2] + dfac[3]*t[3]);\
}

#define PAIR_INT3_FAST(pot, grd, pt, col, r2)                                \
{                                                                            \
  real r2a, istep, chi, dv, d2v, *t;                                         \
  int k;                                                                     \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  k    = (int) (r2a);                                                        \
  chi   = r2a - k;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  t = PTR_2D((pt).table, (col),k, 0,(pt).maxsteps);                          \
  dv  = t[1] - t[0];                                                         \
  d2v = t[2] - 2. * t[1] + t[0];                                             \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot = t[0] + chi * dv + 0.5 * chi * (chi - 1.) * d2v;                      \
  grd = 2. * istep * (dv + (chi - 0.5) * d2v);                               \
}

#define NBLMINLEN 10000

struct ppair {
	int p[4];
};

struct ppair** restrict pairLists = NULL;
int* restrict pairsListLengths = NULL;
int* restrict pairsListMaxLengths = NULL;
real* restrict cutoffRadii = NULL;
int initialized = 0;

//Temporary storages for computed potential gradients
real* restrict grad = NULL;
#ifdef EAM2
real* restrict rho_grad1 = NULL;
real* restrict rho_grad2 = NULL;
#endif
int gradListSize = 0;

void init(void){
	int i,j,m;
	int n = ntypes*ntypes;
	//Allocate the list of the number of pairs
	pairsListLengths = malloc(n* sizeof *pairsListLengths);
	if (pairsListLengths == NULL) error("Cannot allocate memory for pairsListLengths");
	pairsListMaxLengths = malloc(n* sizeof *pairsListMaxLengths);
	if (pairsListMaxLengths == NULL) error("Cannot allocate memory for pairsListMaxLengths");
	for (i=0; i<n; i++){
		pairsListLengths[i] = 0;
		pairsListMaxLengths[i] = NBLMINLEN;
	}

	//Compute/read the individual cut-off radii for the different pairs
	//if different cut-off radii are used for individual steps, use the maximal value
	cutoffRadii = malloc(n * sizeof *cutoffRadii);
	for (i=0; i<ntypes;i++){
		for (j=0; j<ntypes;j++){
			m = i*ntypes + j;
			cutoffRadii[m] = pair_pot.end[m]+nbl_margin;
#ifdef EAM2
			cutoffRadii[m] = MAX(cutoffRadii[m],rho_h_tab.end[m]+nbl_margin);
#endif
		}
	}

	pairLists = malloc(n * sizeof *pairLists);
	if (pairLists == NULL) error("Cannot allocate memory for pairLists");

	for (i=0; i<n; i++){
		pairLists[i] = malloc(pairsListMaxLengths[i] * sizeof *pairLists[i]);
		if (pairLists[i] == NULL)
			error("Cannot allocate memory for pairLists[]");
	}

	initialized = 1;
}

void deallocate_nblist(void){
	if (!initialized) return;
	int i;
	for (i=0; i<ntypes * ntypes;i++)
		free(pairLists[i]);

	free(pairLists); pairLists = NULL;
	free(pairsListLengths); pairsListLengths = NULL;
	free(pairsListMaxLengths); pairsListMaxLengths = NULL;
	free(cutoffRadii);	cutoffRadii = NULL;

	if(gradListSize != 0){
		free(grad);	grad = NULL;
#ifdef EAM2
		free(rho_grad1);	rho_grad1 = NULL;
		free(rho_grad2);	rho_grad2 = NULL;
#endif
		gradListSize = 0;
	}


	have_valid_nbl = 0;
	initialized = 0;
}

void make_nblist(void){
	int i,j, k, n;

	int* fullyFixedTypes = malloc(vtypes * sizeof *fullyFixedTypes);
	for (i=0; i<vtypes; i++)
		fullyFixedTypes[i] = (((restrictions + i)->x + (restrictions + i)->y + (restrictions + i)->z)==0 ? 1 : 0);

	/* update reference positions */
#ifdef OMP
#pragma omp parallel for
#endif
	for (k = 0; k < ncells; k++) {
		cell *p = cell_array + cnbrs[k].np;
		const int n = p->n;
#ifdef INTEL_SIMD
#pragma ivdep
#endif
		for (i = 0; i < n; i++) {
			NBL_POS(p, i, X) = ORT(p, i, X);
			NBL_POS(p, i, Y) = ORT(p, i, Y);
			NBL_POS(p, i, Z) = ORT(p, i, Z);
		}
	}

	n = ntypes*ntypes;

	for (j = 0; j<n; j++)
		pairsListLengths[j] = 0;

	/* for all cells */
	int c;
	for (c = 0; c < ncells; c++) {
		int i, c1 = cnbrs[c].np;
		cell *p = cell_array + c1;

		/* for each atom in cell */
		for (i = 0; i < p->n; i++) {
			int m;
			vektor d1;

			d1.x = ORT(p, i, X);
			d1.y = ORT(p, i, Y);
			d1.z = ORT(p, i, Z);
			int is = SORTE(p,i);
			int ivs = VSORTE(p,i);

			/* for each neighboring atom */
			for (m = 0; m < NNBCELL; m++) {

				int c2, jstart, j;
				real r2;
				cell *q;

				c2 = cnbrs[c].nq[m];
				if (c2 < 0) continue;
				if (c2 == c1) jstart = i + 1;
				else jstart = 0;

				q = cell_array + c2;

				for (j = jstart; j < q->n; j++) {
					vektor d;
					d.x = ORT(q,j,X) - d1.x;
					d.y = ORT(q,j,Y) - d1.y;
					d.z = ORT(q,j,Z) - d1.z;

					int js = SORTE(q,j);

					r2 = SPROD(d, d);
					n = is*ntypes + js;


						int jvs = VSORTE(q,j);
						if ( (fullyFixedTypes[ivs] && fullyFixedTypes[jvs]) ^ compute_energy_between_fixed_atoms) continue;


					//Test if this pair of atoms is inside the cutoff radius
					if (r2 <= cutoffRadii[n]) {
						k = pairsListLengths[n]++;
						pairLists[n][k].p[0] = c1;
						pairLists[n][k].p[1] = c2;
						pairLists[n][k].p[2] = i;
						pairLists[n][k].p[3] = j;
						//Double the size of the table, if running out of space
						if(pairsListLengths[n] == pairsListMaxLengths[n]){
							pairsListMaxLengths[n] *= 2;
							pairLists[n] = realloc(pairLists[n], pairsListMaxLengths[n]*sizeof *pairLists[n]);
							if (pairLists[n] == NULL)
								error("Cannot reallocate memory for pairLists[]");
						}
					}
				}
			}
		}
	}

	free(fullyFixedTypes);
	have_valid_nbl = 1;
	nbl_count++;
}

/******************************************************************************
 *
 *  calc_forces
 *
 ******************************************************************************/

void calc_forces(int steps){
	int anyPbcEnabled = SPROD(pbc_dirs, pbc_dirs);

	int i,j,k,n;
	const int nPairs = ntypes*ntypes;

	if (!initialized){
		init();
		initialized = 1;
	}

	if (0 == have_valid_nbl) {
#ifdef MPI
		/* check message buffer size */
		if (0 == nbl_count % BUFSTEP) setup_buffers();
#endif
		/* update cell decomposition */
		fix_cells();
	}

	/* fill the buffer cells */
	if (anyPbcEnabled)
		send_cells(copy_cell, pack_cell, unpack_cell);

	/* make new neighbor lists */
	if (0 == have_valid_nbl) make_nblist();


	/* clear global accumulation variables */
	tot_pot_energy = 0.0;

	virial = 0.0;
	vir_xx = 0.0;
	vir_yy = 0.0;
	vir_xy = 0.0;
	vir_zz = 0.0;
	vir_yz = 0.0;
	vir_zx = 0.0;
	nfc++;

	/* clear per atom accumulation variables, also in buffer cells */
#ifdef OMP
#pragma omp parallel for
#endif
	for (k = 0; k < nallcells; k++) {
		cell *p = cell_array + k;
		const int n = p->n;
#ifdef INTEL_SIMD
#pragma ivdep
#endif
		for (i = 0; i < n; i++) {
			KRAFT(p,i,X) = 0.0;
			KRAFT(p,i,Y) = 0.0;
			KRAFT(p,i,Z) = 0.0;
			POTENG(p,i)  = 0.0;
#ifdef EAM2
			EAM_RHO(p,i) = 0.0;
#endif
#if defined(STRESS_TENS)
			PRESSTENS(p,i,xx) = 0.0;
			PRESSTENS(p,i,yy) = 0.0;
			PRESSTENS(p,i,xy) = 0.0;
			PRESSTENS(p,i,zz) = 0.0;
			PRESSTENS(p,i,yz) = 0.0;
			PRESSTENS(p,i,zx) = 0.0;
#endif
		}
	}

	/* clear total forces */
#ifdef RIGID
	if ( nsuperatoms>0 )
	for(i=0; i<nsuperatoms; i++) {
		superforce[i].x = 0.0;
		superforce[i].y = 0.0;
		superforce[i].z = 0.0;
	}
#endif


	int sumList = 0;
	for (i = 0; i<nPairs; i++){
		sumList += pairsListLengths[i];
	}

	if(sumList>gradListSize){
		gradListSize = (int)(nbl_size*sumList);
		if(grad) free(grad);
		grad = malloc(gradListSize * sizeof *grad);
		if (grad == NULL) error("Cannot allocate memory for grad");
#ifdef EAM2
		if(rho_grad1) free(rho_grad1);
		if(rho_grad2) free(rho_grad2);
		rho_grad1 = malloc(gradListSize * sizeof *rho_grad1);
		if (rho_grad1 == NULL) error("Cannot allocate memory for rho_grad1");
		rho_grad2 = malloc(gradListSize * sizeof *rho_grad2);
		if (rho_grad2 == NULL) error("Cannot allocate memory for rho_grad2");
#endif
	}

	int cellpairOffset = 0;
	for (n = 0; n<nPairs; n++){
		const struct ppair* restrict pair = pairLists[n];

		const int m = pairsListLengths[n];

		const real potBegin = pair_pot.begin[n];
		const real potEnd = pair_pot.end[n];

		const int type1 = n / ntypes;
		const int type2 = n % ntypes;
		const int col1 = n;
		const int col2 = type2 * ntypes + type1;

#ifdef VIRTUAL_ATOMS
		const int firstIsReal = real_atom[type1];
		const int secondIsReal = real_atom[type2];

		const int first_addResults = (secondIsReal || !firstIsReal) ? 1:0;
		const int second_addResults = (firstIsReal || !secondIsReal) ? 1:0;

		//Skip interaction of virtual atoms in certain steps
		int ignore_non_real_step=steps%update_step;
		if(ignore_non_real_step && (!firstIsReal || !secondIsReal)) continue;
#else
		const int first_addResults = 1;
		const int second_addResults = 1;
#endif

#ifdef EAM2
		const real rhoBeginCol1 = rho_h_tab.begin[col1];
		const real rhoEndCol1 = rho_h_tab.end[col1];
		const real rhoBeginCol2 = rho_h_tab.begin[col2];
		const real rhoEndCol2 = rho_h_tab.end[col2];
#endif

		//Precompute squared distances and the pair-potential part
		for (i=0; i<m; i++){
			vektor v;
			cell *p = cell_array+pair[i].p[0];
			cell *q = cell_array+pair[i].p[1];
			int n_i = pair[i].p[2];
			int n_j = pair[i].p[3];
			//Compute squared distance
			v.x = ORT(q, n_j, X) - ORT(p, n_i, X);
			v.y = ORT(q, n_j, Y) - ORT(p, n_i, Y);
			v.z = ORT(q, n_j, Z) - ORT(p, n_i, Z);
			real rd = SPROD(v,v);

			if (rd < pair_pot.end[n]){
				//Clamp distance into the valid range of the potential table
				real r = MAX( 0.0, rd - potBegin);
				//Compute pair potential and store epot and its gradient
				real epot;
				PAIR_INT_FAST(epot, grad[i+cellpairOffset], pair_pot, col1, r);
				if(first_addResults)  POTENG(p, n_i) += epot * 0.5;
				if(second_addResults) POTENG(q, n_j) += epot * 0.5;
			}

#ifdef EAM2
			//Compute the electron density and the gradients
			if(type1 == type2){
				//Both atoms are of the same type. Electron density and gradient are
				//identical for both
				if (rd < rhoEndCol1){
					//Clamp distance into the valid range of the potential table
					real r = MAX( 0.0, rd - rhoBeginCol1);
					//Compute electron density rho and the gradient, store the gradient for later
					real rho;
					PAIR_INT_FAST(rho, rho_grad1[i+cellpairOffset], rho_h_tab, n, r);
					EAM_RHO(p, n_i) += rho;
					EAM_RHO(q, n_j) += rho;
				}
			} else {
				//Different atomic types in the interaction
				//Compute individual electron density
				if (rd < rhoEndCol1){
					real r = MAX( 0.0, rd-rhoBeginCol1);
					real rho;
					PAIR_INT_FAST(rho, rho_grad1[i+cellpairOffset], rho_h_tab, col1, r);
					EAM_RHO(p, n_i) += rho;
				} else {
					rho_grad1[i+cellpairOffset] = 0.;
				}
				if (rd < rhoEndCol2){
					real s = MAX( 0.0, rd-rhoBeginCol2);
					real rho;
					PAIR_INT_FAST(rho, rho_grad2[i+cellpairOffset], rho_h_tab, col2, s);
					EAM_RHO(q, n_j) += rho;
				} else {
					rho_grad2[i+cellpairOffset] = 0.;
				}

			}
#endif

		}

		cellpairOffset+=m;	//Increase the offset for gradient arrays
	} // pairs n

#ifdef EAM2

	//Sum up rho over buffer cells
	if (anyPbcEnabled)
		send_forces(add_rho,pack_rho,unpack_add_rho);


	if (!compute_energy_between_fixed_atoms)
	/* compute embedding energy and its derivative */
	for (k=0; k<ncells; k++) {
		cell *p = CELLPTR(k);
		real pot, tmp, tr;
		const int n=p->n;
#ifdef INTEL_SIMD
#pragma ivdep
#endif
#ifdef OMP
#pragma omp parallel for
#endif
		for (i=0; i<n; i++) {
			int sorte = SORTE(p,i);
			EAM_RHO(p,i) += REF_RHO(p,i);
			real r = MAX( 0.0, EAM_RHO(p,i) - embed_pot.begin[sorte]);
			PAIR_INT_FAST(pot, EAM_DF(p,i), embed_pot, sorte, r);
			POTENG(p,i) += pot;
		}
	}

	/* distribute derivative of embedding energy */
	if (anyPbcEnabled)
		send_cells(copy_dF,pack_dF,unpack_dF);
#endif

	//Reduce all results
	cellpairOffset = 0;
	if (!compute_energy_between_fixed_atoms)
	for (n = 0; n < nPairs; n++) {
		const struct ppair* restrict pair = pairLists[n];

		const int m = pairsListLengths[n];

		const int type1 = n / ntypes;
		const int type2 = n % ntypes;
		const int col1 = n;
		const int col2 = type2 * ntypes + type1;

#ifdef VIRTUAL_ATOMS
		const int firstIsReal = real_atom[type1];
		const int secondIsReal = real_atom[type2];
		const int bothReal = (firstIsReal && secondIsReal) ? 1 : 0;

		const int first_addResults = (firstIsReal || !secondIsReal) ? 1:0;
		const int second_addResults = (secondIsReal || !firstIsReal) ? 1:0;

		//Skip interaction of virtual atoms in certain steps
		int ignore_non_real_step=steps%update_step;
		if(ignore_non_real_step && (!firstIsReal || !secondIsReal)) continue;
#else
		const int first_addResults = 1;
		const int second_addResults = 1;
		const int bothReal = 1;
#endif

#ifdef EAM2
		real rhoCut = MAX(rho_h_tab.end[col1], rho_h_tab.end[col2]);
#endif

		for (i=0; i<m; i++) {
			vektor v, force;

			cell *p = cell_array+pair[i].p[0];
			cell *q = cell_array+pair[i].p[1];
			int n_i = pair[i].p[2];
			int n_j = pair[i].p[3];

			v.x = ORT(q, n_j, X) - ORT(p, n_i, X);
			v.y = ORT(q, n_j, Y) - ORT(p, n_i, Y);
			v.z = ORT(q, n_j, Z) - ORT(p, n_i, Z);
			real r2 = SPROD(v,v);

			real g = 0.;
			if (r2 < pair_pot.end[n]){
				g = grad[i+cellpairOffset];
			}

#ifdef EAM2
			//Add the gradient of the electron density if needed
			if (r2 < rhoCut) {
				if (type1==type2)
					g += 0.5 * (EAM_DF(p,n_i) + EAM_DF(q,n_j)) * rho_grad1[i+cellpairOffset];
				else
					g += 0.5 * (EAM_DF(p,n_i) * rho_grad1[i+cellpairOffset] + EAM_DF(q,n_j) * rho_grad2[i+cellpairOffset]);
			}
#endif

			if(g!=0){
				force.x = v.x * g;
				force.y = v.y * g;
				force.z = v.z * g;

				if(first_addResults){
					KRAFT(q, n_j,X) -= force.x;
					KRAFT(q, n_j,Y) -= force.y;
					KRAFT(q, n_j,Z) -= force.z;
				}

				if(second_addResults){
					KRAFT(p, n_i,X) += force.x;
					KRAFT(p, n_i,Y) += force.y;
					KRAFT(p, n_i,Z) += force.z;
				}

				if(bothReal){
#ifdef P_AXIAL
				vir_xx -= v.x * force.x;
				vir_yy -= v.y * force.y;
				vir_zz -= v.z * force.z;
#else
				virial       -= SPROD(v,force);
#endif
				}
#ifdef STRESS_TENS
				if (do_press_calc) {
					/* avoid double counting of the virial */
					force.x *= 0.5;
					force.y *= 0.5;
					force.z *= 0.5;

					PRESSTENS(p, n_i,xx) -= v.x * force.x * first_addResults;
					PRESSTENS(p, n_i,yy) -= v.y * force.y * first_addResults;
					PRESSTENS(p, n_i,xy) -= v.x * force.y * first_addResults;
					PRESSTENS(p, n_i,zz) -= v.z * force.z * first_addResults;
					PRESSTENS(p, n_i,yz) -= v.y * force.z * first_addResults;
					PRESSTENS(p, n_i,zx) -= v.z * force.x * first_addResults;
					PRESSTENS(q, n_j,zx) -= v.z * force.x * second_addResults;
					PRESSTENS(q, n_j,xx) -= v.x * force.x * second_addResults;
					PRESSTENS(q, n_j,yy) -= v.y * force.y * second_addResults;
					PRESSTENS(q, n_j,zz) -= v.z * force.z * second_addResults;
					PRESSTENS(q, n_j,xy) -= v.x * force.y * second_addResults;
					PRESSTENS(q, n_j,yz) -= v.y * force.z * second_addResults;
				}
#endif
			}
		}
		cellpairOffset+=m;
	}//n pairs

#ifdef VIRTUAL_ATOMS
	tot_pot_energy_placeholders = 0.0;
#endif

	//Sum total potential energy
	for (k=0; k<nallcells; k++) {
		cell *p = cell_array+k;
		int i;
		const int n = p->n;
		for (i=0; i<n; i++){
#ifdef VIRTUAL_ATOMS
          if (!real_atom[SORTE(p,i)]) tot_pot_energy_placeholders+=POTENG(p,i);
          else
#endif
      	  tot_pot_energy+=POTENG(p,i);
		}
	}

#ifdef MPI
	real tmpvec1[8], tmpvec2[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	/* sum up results of different CPUs */
	tmpvec1[0] = tot_pot_energy;
	tmpvec1[1] = virial;
	tmpvec1[2] = vir_xx;
	tmpvec1[3] = vir_yy;
	tmpvec1[4] = vir_zz;
	tmpvec1[5] = vir_xy;
	tmpvec1[6] = vir_yz;
	tmpvec1[7] = vir_zx;
	MPI_Allreduce( tmpvec1, tmpvec2, 8, REAL, MPI_SUM, cpugrid);
	tot_pot_energy = tmpvec2[0];
	virial = tmpvec2[1];
	vir_xx = tmpvec2[2];
	vir_yy = tmpvec2[3];
	vir_zz = tmpvec2[4];
	vir_xy = tmpvec2[5];
	vir_yz = tmpvec2[6];
	vir_zx = tmpvec2[7];

#ifdef VIRTUAL_ATOMS
	if (!(steps%update_step)) {
		/* sum up results of different CPUs */
		MPI_Allreduce(&tot_pot_energy_placeholders, tmpvec2, 1, REAL, MPI_SUM, cpugrid);
		tot_pot_energy_placeholders = tmpvec2[0];
	}
#endif
#endif

	/* add forces back to original cells/cpus */
	if (anyPbcEnabled)
		send_forces(add_forces, pack_forces, unpack_forces);

	if (compute_energy_between_fixed_atoms){
		for (k = 0; k < nallcells; k++) {
			cell *p = cell_array + k;
			const int n = p->n;
			for (i = 0; i < n; i++) {
				REF_EPOT(p,i) = POTENG(p,i) ;
#ifdef EAM2
				REF_RHO(p,i)  = EAM_RHO(p,i);
#endif
			}
		}
	}
	
	for (k = 0; k < NCELLS; k++) {
		cell *p = CELLPTR(k);
		const int n = p->n;
		for (i = 0; i < n; i++) {
			POTENG(p,i) += REF_EPOT(p, i);
		}
	}
}

/******************************************************************************
 *
 *  check_nblist
 *
 ******************************************************************************/

void check_nblist(){
	real r2, max1 = 0.0, max2;
	vektor d;
	int k;

	/* compare with reference positions */
	for (k = 0; k < NCELLS; k++) {
		int i;
		cell *p = CELLPTR(k);
		const int n = p->n;
		for (i = 0; i < n; i++) {
			d.x = ORT(p,i,X) - NBL_POS(p, i, X);
			d.y = ORT(p,i,Y) - NBL_POS(p, i, Y);
			d.z = ORT(p,i,Z) - NBL_POS(p, i, Z);

			r2 = SPROD(d, d);
			if (r2 > max1) max1 = r2;
		}
	}

#ifdef MPI
	MPI_Allreduce( &max1, &max2, 1, REAL, MPI_MAX, cpugrid);
#else
	max2 = max1;
#endif
	if (max2 > SQR(0.5 * nbl_margin)) have_valid_nbl = 0;
}

