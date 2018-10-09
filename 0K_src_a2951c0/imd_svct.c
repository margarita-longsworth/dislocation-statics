/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_svct.c -- Routines for Slipvector-analysis
*
******************************************************************************/

/******************************************************************************
* $Revision: 1.1 $
* $Date: 2009/22/04 11:54:51 $
******************************************************************************/

#include "imd.h"

void init_svct(void) {
	int k;

	/* update neighbor table cutoff */
	if (ada_nbr_r2cut == 0.) {
		error("ada_nbr_r2cut is 0., cannot compute SlipVector");
	}

	for (k = 0; k < NCELLS; k++) {
		int i;
		cell* p;
		p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			SLIPREFERENCE(p, i, X) = ORT(p,i,X);
			SLIPREFERENCE(p, i, Y) = ORT(p,i,Y);
			SLIPREFERENCE(p, i, Z) = ORT(p,i,Z);
			SLIPVECTOR(p, i, X) = 0.;
			SLIPVECTOR(p, i, Y) = 0.;
			SLIPVECTOR(p, i, Z) = 0.;
		}
	}

}

void do_svct(void) {
	int j, k, nneigh, i, nslip;
	vektor sv, d;
	/*real hx = box_x.x * (cell_dim.x - 2) / global_cell_dim.x;
	real hy = box_y.y * (cell_dim.y - 2) / global_cell_dim.y;
	real hz = box_z.z * (cell_dim.z - 2) / global_cell_dim.z;*/

	do_neightab_complete();

	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {

			/*if (pic_ur.x != (real) 0)*/ /* if pic_ur.x still 0, write everything */

				/*if ((ORT(p,i,X) < pic_ll.x) || (ORT(p,i,X) > pic_ur.x)
						|| (ORT(p,i,Z) < pic_ll.z) || (ORT(p,i,Z) > pic_ur.z)
						||
						(ORT(p,i,Y) < pic_ll.y) || (ORT(p,i,Y) > pic_ur.y))
					continue; */

			nslip = 0;
			nneigh = NEIGH(p,i)->n;

			sv.x = 0.;
			sv.y = 0.;
			sv.z = 0.;

			for (j = 0; j < nneigh; j++) {
				cell *q = NEIGH(p,i)->cl[j];
				int ii = NEIGH(p,i)->num[j];

				d.x = SLIPREFERENCE(p,i,X) - SLIPREFERENCE(q,ii,X);
				d.y = SLIPREFERENCE(p,i,Y) - SLIPREFERENCE(q,ii,Y);
				d.z = SLIPREFERENCE(p,i,Z) - SLIPREFERENCE(q,ii,Z);

				/*d.x -= NEIGH(p,i)->dist[3* j ];
				d.y -= NEIGH(p,i)->dist[3* j + 1];
				d.z -= NEIGH(p,i)->dist[3* j + 2];*/
				d.x -= ORT(p,i,X) - ORT(q,ii,X);
				d.y -= ORT(p,i,Y) - ORT(q,ii,Y);
				d.z -= ORT(p,i,Z) - ORT(q,ii,Z);

/*
				reduce_displacement(&d);

				while (d.x >= hx - 0.1)
					d.x -= hx;
				while (d.x <= -hx + 0.1)
					d.x += hx;
				while (d.y >= hy - 0.1)
					d.y -= hy;
				while (d.y <= -hy + 0.1)
					d.y += hy;
				while (d.z >= hz - 0.1)
					d.z -= hz;
				while (d.z <= -hz + 0.1)
					d.z += hz;*/

				/*
				 * Ignore slipped neighbors with relative
				 * slip less than 0.4A (0.4^2 = 0.16)
			 	 */
				if (SPROD(d,d) > 0.16) {
					sv.x += d.x;
					sv.y += d.y;
					sv.z += d.z;
					nslip++;
				}
			}

			if (nslip > 0) {
				sv.x /= nslip;
				sv.y /= nslip;
				sv.z /= nslip;
			}

			SLIPVECTOR(p,i,X) += sv.x;
			SLIPVECTOR(p,i,Y) += sv.y;
			SLIPVECTOR(p,i,Z) += sv.z;
		}
	}
}

void update_svct_ref(void) {
	int k;
	for (k = 0; k < NCELLS; k++) {
		int i;
		cell* p;
		p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			SLIPREFERENCE(p, i, X) = ORT(p,i,X);
			SLIPREFERENCE(p, i, Y) = ORT(p,i,Y);
			SLIPREFERENCE(p, i, Z) = ORT(p,i,Z);
		}
	}
}
