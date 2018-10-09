/******************************************************************************
 *
 * IMD -- The ITAP Molecular Dynamics Program
 *
 * Copyright 1996-2006 Institute for Theoretical and Applied Physics,
 * University of Stuttgart, D-70550 Stuttgart
 *
 ******************************************************************************/

/******************************************************************************
 *
 * imd_stress -- calculates stress field
 *
 * Voronoi analysis according to Allen, Tildesley
 * Parameter reading and cell decomposition from imd_pair
 *
 ******************************************************************************/

/******************************************************************************
 * $Revision: 1.9 $
 * $Date: 2006/07/12 09:08:49 $
 ******************************************************************************/

#include "imd.h"

void write_header_stress(FILE *out) {
	char c;
	time_t now;

	/* format line */
	if (binary_output)
		c = is_big_endian ? 'b' : 'l';
	else
		c = 'A';


	/* contents line */
#ifndef TWOD
	fprintf(out, "#F %c 1 3 7\n", c);
	fprintf(out, "#C number x y z s_xx s_yy s_zz s_yz s_zx s_xy s_vM\n");
#else
	fprintf(out, "#F %c 1 2 4\n", c);
	fprintf(out, "#C number x y s_xx s_yy s_xy s_vM\n");
#endif

	/* box lines */
	fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x, box_x.y, box_x.z);
	fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x, box_y.y, box_y.z);
	fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x, box_z.y, box_z.z);

	/* generation date and endheader line */
	time(&now);
	fprintf(out, "## Generated on %s", ctime(&now));
	fprintf(out, "## by %s (version of %s)\n", progname, DATE);
	fprintf(out, "#E\n");
}

void write_atoms_stress(FILE *out) {
	int i, k, n, valid, len = 0;
	float tmp;
	vektor cellSize;
	ivektor cellCoords;

	i_or_f *data;
	sym_tensor average;

	cellSize.x = box_x.x /  global_cell_dim.x;
	cellSize.y = box_y.y /  global_cell_dim.y;
	cellSize.z = box_z.z /  global_cell_dim.z;

	if (fullAtomicResolution){
		/* Output temporal averaged stress per atom */
		for (k = 0; k < NCELLS; k++) {
			cell *p = CELLPTR(k);
			for (i = 0; i < p->n; i++) {
				if (p->vol[i] > 0.) tmp = 1. / (p->vol[i]*stressAverage_stepsAveraged);
				else tmp = 0.;

				if (binary_output) {
					n = 0;
					data = (i_or_f *) (outbuf + len);
					data[n++].i = (int) NUMMER(p,i);
					data[n++].f = (float) ORT (p,i,X);
					data[n++].f = (float) ORT (p,i,Y);
#ifndef TWOD
					data[n++].f = (float) ORT (p,i,Z);
#endif
					data[n++].f = (float) PRESSTENSSUM(p,i,xx) * tmp;
					data[n++].f = (float) PRESSTENSSUM(p,i,yy) * tmp;
#ifndef TWOD
					data[n++].f = (float) PRESSTENSSUM(p,i,zz) * tmp;
					data[n++].f = (float) PRESSTENSSUM(p,i,yz) * tmp;
					data[n++].f = (float) PRESSTENSSUM(p,i,zx) * tmp;
#endif
					data[n++].f = (float) PRESSTENSSUM(p,i,xy) * tmp;
					data[n++].f = (float) getVonMisesStress(p,i);
					len += n * sizeof(i_or_f);
				}
				/* ASCII output */
				else {
#ifdef TWOD
					len += sprintf(outbuf + len, "%d %12f %12f %12f %12f %12f %12f\n",
							NUMMER(p,i), ORT(p,i,X), ORT(p,i,Y),
							PRESSTENSSUM(p,i,xx) * tmp, PRESSTENSSUM(p,i,yy) * tmp,
							PRESSTENSSUM(p,i,xy) * tmp), getVonMisesStress(p,i);
#else
					len += sprintf(outbuf + len, "%d %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f\n",
						NUMMER(p,i), ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z),
						PRESSTENSSUM(p,i,xx) * tmp, PRESSTENSSUM(p,i,yy) * tmp,
						PRESSTENSSUM(p,i,zz) * tmp, PRESSTENSSUM(p,i,yz) * tmp,
						PRESSTENSSUM(p,i,zx) * tmp, PRESSTENSSUM(p,i,xy) * tmp,
						getVonMisesStress(p,i));
#endif
				}
				/* flush or send outbuf if it is full */
				if (len > outbuf_size - 256)
					flush_outbuf(out, &len, OUTBUF_TAG);
			}
		}
	} else {
		/* Output average of temporal averaged stress in whole cell */
		for (k = 0; k < NCELLS; k++) {
			cell *p = CELLPTR(k);

			average.xx = 0.;
			average.yy = 0.;
			average.xy = 0.;
#ifndef TWOD
			average.zz = 0.;
			average.yz = 0.;
			average.zx = 0.;
#endif

			valid = 0;
			for (i = 0; i < p->n; i++) {
				if (p->vol[i] > 0.) {
					tmp = 1. / p->vol[i];
					valid++;

					average.xx += PRESSTENSSUM(p, i, xx) * tmp;
					average.yy += PRESSTENSSUM(p, i, yy) * tmp;
#ifndef TWOD
					average.zz += PRESSTENSSUM(p, i, zz) * tmp;
					average.yz += PRESSTENSSUM(p, i, yz) * tmp;
					average.zx += PRESSTENSSUM(p, i, zx) * tmp;
#endif
					average.xy += PRESSTENSSUM(p, i, xy) * tmp;
				}
			}

			if (valid > 0) {
				/* Calculate the global cell coordinate from the first atom in a cell*/
#ifdef TWOD
				cellCoords = cell_coord(ORT(p,0,X), ORT(p,0,Y));
#else
				cellCoords = cell_coord(ORT(p,0,X), ORT(p,0,Y), ORT(p,0,Z));
#endif

				tmp = 1. / (valid * stressAverage_delta);

				if (binary_output) {
					n = 0;

					data = (i_or_f *) (outbuf + len);
					data[n++].i = (int) myid * NCELLS + k;
					data[n++].f = (float) cellSize.x * (cellCoords.x + 0.5);
					data[n++].f = (float) cellSize.y * (cellCoords.y + 0.5);
#ifndef TWOD
					data[n++].f = (float) cellSize.z * (cellCoords.z + 0.5);
#endif
					data[n++].f = (float) (average.xx * tmp);
					data[n++].f = (float) (average.yy * tmp);
#ifndef TWOD
					data[n++].f = (float) (average.zz * tmp);
					data[n++].f = (float) (average.yz * tmp);
					data[n++].f = (float) (average.zx * tmp);
#endif
					data[n++].f = (float) (average.xy * tmp);

					len += n * sizeof(i_or_f);
				}
				/* ASCII output */
				else {

#ifdef TWOD
					len += sprintf(outbuf + len, "%d %12f %12f %12f %12f %12f \n",
							myid*NCELLS+k, cellSize.x * (cellCoords.x + 0.5),
							cellSize.y * (cellCoords.y + 0.5),
							average.xx * tmp, average.yy * tmp, average.xy * tmp);
#else
					len += sprintf(outbuf + len, "%d %12f %12f %12f %12f %12f %12f %12f %12f %12f \n",
								myid * NCELLS + k, cellSize.x * (cellCoords.x + 0.5),
								cellSize.y * (cellCoords.y + 0.5), cellSize.z * (cellCoords.z + 0.5),
								average.xx * tmp, average.yy * tmp,
								average.zz * tmp, average.yz * tmp,
								average.zx * tmp, average.xy * tmp);
#endif
				}
			}
			/* flush or send outbuf if it is full */
			if (len > outbuf_size - 256)
				flush_outbuf(out, &len, OUTBUF_TAG);
		}
	}
	flush_outbuf(out, &len, OUTBUF_TAG + 1);
}

void initPressTensorSum(void){
	if (ada_nbr_r2cut == 0.) {
		error("ada_nbr_r2cut is 0., cannot compute stress averages");
	}
	resetPressTensorSum();
}

void resetPressTensorSum(void) {
	int i,k;
	cell* p;

	if (!stressAverage_addedThisStep) return;
	for (k = 0; k < NCELLS; k++) {
		p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			PRESSTENSSUM(p,i,xx) = 0.;
			PRESSTENSSUM(p,i,yy) = 0.;
			PRESSTENSSUM(p,i,xy) = 0.;
#ifndef TWOD
			PRESSTENSSUM(p,i,yz) = 0.;
			PRESSTENSSUM(p,i,zx) = 0.;
			PRESSTENSSUM(p,i,zz) = 0.;
#endif
			VOLUME(p,i) = 0.;
		}
	}
	stressAverage_stepsAveraged = 0;
}

real getVonMisesStress(cell *p, int i){
	real tmp;

	if (p->vol[i] > 0.) tmp = 1. / (p->vol[i]*stressAverage_stepsAveraged);
	else tmp = 0.;

	real s_xx = PRESSTENSSUM(p,i,xx) * tmp;
	real s_yy = PRESSTENSSUM(p,i,yy) * tmp;
	real s_xy = PRESSTENSSUM(p,i,xy) * tmp;

#ifndef TWOD
	real s_yz = PRESSTENSSUM(p,i,yz) * tmp;
	real s_zx = PRESSTENSSUM(p,i,zx) * tmp;
	real s_zz = PRESSTENSSUM(p,i,zz) * tmp;

	real vM = SQR(s_xx - s_yy) + SQR(s_yy - s_zz) + SQR(s_zz - s_xx);
	vM += 6*(SQR(s_yz) + SQR(s_zx) + SQR(s_xy));
	vM *= 0.5;
#else
	real vM = SQR(s_xx) + SQR(s_yy) - s_xx * s_yy + 3*SQR(s_xy);
#endif
	return SQRT(vM);
}

void addPressTensor(void) {
	int i,k;
	cell* p;

	if (stressAverage_addedThisStep) return;
	for (k = 0; k < NCELLS; k++) {
		p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			PRESSTENSSUM(p,i,xx) += PRESSTENS(p,i,xx);
			PRESSTENSSUM(p,i,yy) += PRESSTENS(p,i,yy);
			PRESSTENSSUM(p,i,xy) += PRESSTENS(p,i,xy);
#ifndef TWOD
			PRESSTENSSUM(p,i,zz) += PRESSTENS(p,i,zz);
			PRESSTENSSUM(p,i,yz) += PRESSTENS(p,i,yz);
			PRESSTENSSUM(p,i,zx) += PRESSTENS(p,i,zx);
#endif
		}
	}
	stressAverage_stepsAveraged++;
	stressAverage_addedThisStep = 1;
}

/******************************************************************************
 *
 *  calculateVolumePerAtom -- calculates neighbouring points for Voronoi construction
 *             and calls calculation of volume/area
 *
 ******************************************************************************/

void calculateVolumePerAtom(void) {
	cell *p, *q;
	int i, j, k, ii;
	int nneigh;
	real atomic_volume = 0.;

	do_neightab_complete();

	vektor *candcoord = malloc(neigh_len * sizeof(vektor));
	real *canddist2 = malloc(neigh_len * sizeof(real));
	if (canddist2 == NULL || candcoord == NULL)
		error("Cannot allocate arrays in calculateVolumePerAtom");

	/* for each cell */
	for (k = 0; k < NCELLS; k++) {
		p = CELLPTR(k);

		/* for each atom in first cell */
		for (i = 0; i < p->n; ++i) {
			nneigh = NEIGH(p, i)->n;

			/* If there are less than four (three) points, a polyhedron (polygon) cannot
			 be constructed */
#ifndef TWOD
			if (nneigh < 4)
				atomic_volume = 0.0;
#else
			if (nneigh < 3 ) atomic_volume = 0.0;
#endif
			else {
				for (j = 0; j < nneigh; j++) {
					q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];

					candcoord[j].x = ORT(p,i,X) - ORT(q,ii,X);
					candcoord[j].y = ORT(p,i,Y) - ORT(q,ii,Y);
#ifndef TWOD
					candcoord[j].z = ORT(p,i,Z) - ORT(q,ii,Z);
#endif
					/*reduce_displacement(&candcoord[j]);*/
					canddist2[j] = SPROD(candcoord[j], candcoord[j]);
				}

				/* Sort candidates in ascending order of distance */
				sort(nneigh, candcoord, canddist2);

				/* Perform Voronoi analysis */
#ifdef TWOD
				atomic_volume = do_voronoi_2d(nneigh, candcoord, canddist2);
#else
				atomic_volume = do_voronoi_3d(nneigh, candcoord, canddist2);
#endif
			}

			VOLUME(p,i) = atomic_volume;
		}
	}

	free(candcoord);
	free(canddist2);
}

/******************************************************************************
 *
 *  sort -- Sorts candidates for neighbour atoms in increasing order of distance
 *
 ******************************************************************************/

void sort(int neighnum, vektor *candcoord, real *canddist2) {
	int i, j;
	vektor tmp;
	real tmpdist2;

	for (i = (neighnum - 1); i > 0; --i)
		for (j = 0; j < i; ++j)
			if (canddist2[j] > canddist2[j + 1]) {
				tmp.x = candcoord[j].x;
				tmp.y = candcoord[j].y;
#ifndef TWOD
				tmp.z = candcoord[j].z;
#endif
				tmpdist2 = canddist2[j];

				candcoord[j].x = candcoord[j + 1].x;
				candcoord[j].y = candcoord[j + 1].y;
#ifndef TWOD
				candcoord[j].z = candcoord[j + 1].z;
#endif
				canddist2[j] = canddist2[j + 1];

				candcoord[j + 1].x = tmp.x;
				candcoord[j + 1].y = tmp.y;
#ifndef TWOD
				candcoord[j + 1].z = tmp.z;
#endif
				canddist2[j + 1] = tmpdist2;
			}

}

#ifndef TWOD

/******************************************************************************
 *
 *  do_voronoi_3d -- Calculates Voronoi cells and volume, 3d version
 *
 ******************************************************************************/

real do_voronoi_3d(int neighnum, vektor *candcoord, real *canddist2)

{
    real const TOL = 1.0e-6;
	real const TOL2 = 1.0e-10;
    real const TOL_VERT2 = 1.0e-6;
    real const TOL_DIST2 = 2.0e-11;
    int const NUM = 200;

    typedef vektor vektorstr[NUM];

	int i, j, k, l, n;
	real ab, bc, ca, da, db, dc, det, detinv, tmp;
	vektor icoord, jcoord, kcoord;
	real idist2, jdist2, kdist2, atomic_volume;
	vektor tmpvek, tmpvertex, vertex[NUM];
	int ok, vertexcount, vertexnum, facesnum, edgesnum;
	int *vertexnumi;
	real area_i, height;
	real sin, cos, maxcos;
	int mink, ord[NUM], index[NUM], surfind[NUM];
	vektor *coord, center;
	vektorstr *vertexloc;

	/* Allocate memory for vertices */
	vertexnumi = (int *) malloc(neighnum * sizeof(int));
	coord = (vektor *) malloc(neighnum * sizeof(vektor));
	vertexloc = (vektorstr *) malloc(neighnum * sizeof(vektorstr));

	if (vertexloc == NULL || coord == NULL || vertexnumi == NULL)
		error("Cannot allocate memory for vertices!\n");

	mink = 0;
	vertexcount = 0;
	atomic_volume = 0.0;

	/* Each possible vertex defined by the intersection of 3 planes is examined */
	for (i = 0; i < neighnum - 2; ++i) {
		icoord.x = candcoord[i].x;
		icoord.y = candcoord[i].y;
		icoord.z = candcoord[i].z;
		idist2 = -canddist2[i];

		for (j = i + 1; j < neighnum - 1; ++j) {
			jcoord.x = candcoord[j].x;
			jcoord.y = candcoord[j].y;
			jcoord.z = candcoord[j].z;
			jdist2 = -canddist2[j];

			ab = icoord.x * jcoord.y - jcoord.x * icoord.y;
			bc = icoord.y * jcoord.z - jcoord.y * icoord.z;
			ca = icoord.z * jcoord.x - jcoord.z * icoord.x;
			da = idist2 * jcoord.x - jdist2 * icoord.x;
			db = idist2 * jcoord.y - jdist2 * icoord.y;
			dc = idist2 * jcoord.z - jdist2 * icoord.z;

			for (k = j + 1; k < neighnum; ++k) {
				kcoord.x = candcoord[k].x;
				kcoord.y = candcoord[k].y;
				kcoord.z = candcoord[k].z;
				kdist2 = -canddist2[k];

				det = kcoord.x * bc + kcoord.y * ca + kcoord.z * ab;

				/* Check whether planes intersect */
				if (SQR(det) > TOL2) {
					detinv = 1.0 / det;

					tmpvertex.x
							= (-kdist2 * bc + kcoord.y * dc - kcoord.z * db)
									* detinv;
					tmpvertex.y
							= (-kcoord.x * dc - kdist2 * ca + kcoord.z * da)
									* detinv;
					tmpvertex.z = (kcoord.x * db - kcoord.y * da - kdist2 * ab)
							* detinv;

					/* Check whether vertex belongs to the Voronoi cell */
					l = 0;
					ok = 1;

					do {
						if (l != i && l != j && l != k)
							ok = (SPROD( candcoord[l], tmpvertex )
									<= canddist2[l] + TOL_VERT2);

						++l;

					} while (ok && (l < neighnum));

					if (ok) {
						vertex[vertexcount].x = 0.5 * tmpvertex.x;
						vertex[vertexcount].y = 0.5 * tmpvertex.y;
						vertex[vertexcount].z = 0.5 * tmpvertex.z;

						++vertexcount;
					}

				} /* Planes intersect */

			} /* k */
		} /* j */
	} /* i */

	vertexnum = vertexcount;

	/* Check whether some vertices coincide */
	for (i = 0; i < vertexnum; ++i) {
		index[i] = 1;
		for (j = i + 1; j < vertexnum; ++j) {
			tmpvek.x = vertex[j].x - vertex[i].x;
			tmpvek.y = vertex[j].y - vertex[i].y;
			tmpvek.z = vertex[j].z - vertex[i].z;

			if (SPROD( tmpvek, tmpvek) < TOL2)
				index[i] = 0;
		}
	}

	/* Remove coincident vertices */
	j = 0;
	for (i = 0; i < vertexnum; ++i)
		if (index[i] != 0) {
			vertex[j].x = vertex[i].x;
			vertex[j].y = vertex[i].y;
			vertex[j].z = vertex[i].z;

			++j;
		}

	vertexnum = j;

	/* Number of vertices of Voronoi cell must be greater than 3 */
	if (vertexnum > 3) {
		/* Check whether faces exist */
		facesnum = 0;

		/* Each neighbour atom i corresponds to at most one surface *
		 * Sum over all surfaces */
		for (i = 0; i < neighnum; ++i) {
			/* Coordinates of center of surface i */
			coord[i].x = 0.5 * candcoord[i].x;
			coord[i].y = 0.5 * candcoord[i].y;
			coord[i].z = 0.5 * candcoord[i].z;

			/* Look for vertices that belong to surface i */
			vertexnumi[i] = 0;
			for (j = 0; j < vertexnum; ++j) {
				surfind[j] = 0;

				vertexloc[i][j].x = vertex[j].x - coord[i].x;
				vertexloc[i][j].y = vertex[j].y - coord[i].y;
				vertexloc[i][j].z = vertex[j].z - coord[i].z;

				tmp = SPROD(coord[i],vertexloc[i][j]);

				if (SQR(tmp) < TOL_DIST2) {
					/* vertex j belongs to surface i */
					surfind[j] = 1;
					++vertexnumi[i];
				}
			}

			/* Surface i exists */
			if (vertexnumi[i] > 2) {
				++facesnum;

				/* Compute coordinates of vertices belonging to surface i */
				k = 0;
				for (j = 0; j < vertexnum; ++j)
					if (surfind[j] == 1) {
						vertexloc[i][k].x = vertexloc[i][j].x;
						vertexloc[i][k].y = vertexloc[i][j].y;
						vertexloc[i][k].z = vertexloc[i][j].z;

						++k;
					}
			}
			/* Transform into center of mass system */
			center.x = 0.0;
			center.y = 0.0;
			center.z = 0.0;

			if (vertexnumi[i] > 2) {
				for (j = 0; j < vertexnumi[i]; ++j) {
					center.x += vertexloc[i][j].x;
					center.y += vertexloc[i][j].y;
					center.z += vertexloc[i][j].z;
				}

				tmp = 1.0 / vertexnumi[i];
				center.x *= tmp;
				center.y *= tmp;
				center.z *= tmp;

				for (j = 0; j < vertexnumi[i]; ++j) {
					vertexloc[i][j].x -= center.x;
					vertexloc[i][j].y -= center.y;
					vertexloc[i][j].z -= center.z;
				}

			}

		} /* i */

		/* Number of edges of Voronoi cell */
		edgesnum = 0;

		for (n = 0; n < neighnum; ++n)
			if (vertexnumi[n] > 2)
				edgesnum += vertexnumi[n];

		edgesnum /= 2;

		/* Check whether Euler relation holds */
		if ((vertexnum - edgesnum + facesnum) == 2) {
			/* Compute volume of Voronoi cell */

			/* For all potential faces */
			for (i = 0; i < neighnum; ++i)
				/* Surface i exists */
				if (vertexnumi[i] > 2) {
					/* Sort vertices of face i */
					ord[0] = 0;
					for (j = 0; j < vertexnumi[i] - 1; ++j) {
						maxcos = -1.0;
						for (k = 0; k < vertexnumi[i]; ++k) {
							tmpvek.x = vertexloc[i][k].y
									* vertexloc[i][ord[j]].z
									- vertexloc[i][k].z
											* vertexloc[i][ord[j]].y;
							tmpvek.y = vertexloc[i][k].z
									* vertexloc[i][ord[j]].x
									- vertexloc[i][k].x
											* vertexloc[i][ord[j]].z;
							tmpvek.z = vertexloc[i][k].x
									* vertexloc[i][ord[j]].y
									- vertexloc[i][k].y
											* vertexloc[i][ord[j]].x;

							sin = SPROD( tmpvek, coord[i]);

							if (sin > TOL) {
								cos
										= SPROD( vertexloc[i][k], vertexloc[i][ord[j]] )
												/ sqrt(
														SPROD(vertexloc[i][k],vertexloc[i][k]))
												/ sqrt(
														SPROD(vertexloc[i][ord[j]],vertexloc[i][ord[j]]));
								if (cos > maxcos) {
									maxcos = cos;
									mink = k;
								}
							}
						}

						ord[j + 1] = mink;

					}

					/* Compute area of surface i */
					area_i = 0.0;
					height = sqrt(SPROD( coord[i], coord[i] ));
					tmp = 1.0 / height;

					for (j = 0; j < vertexnumi[i] - 1; ++j) {
						tmpvek.x = vertexloc[i][ord[j + 1]].y
								* vertexloc[i][ord[j]].z - vertexloc[i][ord[j
								+ 1]].z * vertexloc[i][ord[j]].y;
						tmpvek.y = vertexloc[i][ord[j + 1]].z
								* vertexloc[i][ord[j]].x - vertexloc[i][ord[j
								+ 1]].x * vertexloc[i][ord[j]].z;
						tmpvek.z = vertexloc[i][ord[j + 1]].x
								* vertexloc[i][ord[j]].y - vertexloc[i][ord[j
								+ 1]].y * vertexloc[i][ord[j]].x;

						area_i += 0.5 * SPROD( tmpvek, coord[i] ) * tmp;

					}
					tmpvek.x = vertexloc[i][ord[0]].y
							* vertexloc[i][ord[vertexnumi[i] - 1]].z
							- vertexloc[i][ord[0]].z
									* vertexloc[i][ord[vertexnumi[i] - 1]].y;
					tmpvek.y = vertexloc[i][ord[0]].z
							* vertexloc[i][ord[vertexnumi[i] - 1]].x
							- vertexloc[i][ord[0]].x
									* vertexloc[i][ord[vertexnumi[i] - 1]].z;
					tmpvek.z = vertexloc[i][ord[0]].x
							* vertexloc[i][ord[vertexnumi[i] - 1]].y
							- vertexloc[i][ord[0]].y
									* vertexloc[i][ord[vertexnumi[i] - 1]].x;

					area_i += 0.5 * SPROD( tmpvek, coord[i] ) * tmp;

					/* Volume of Voronoi cell */
					atomic_volume += area_i * height / 3.0;

				} /* vertexnum[i] > 2 */

		} /* Euler relation holds */

	} /* Number of vertices > 3 */

	free(vertexloc);
	free(coord);
	free(vertexnumi);

	return atomic_volume;
}

#else

/******************************************************************************
 *
 *  do_voronoi_2d -- Calculates Voronoi cells and volume, 2d version
 *
 ******************************************************************************/

real do_voronoi_2d(int neighnum, vektor *candcoord, real *canddist2)

{
	int i, j, l, n;
	real det, detinv;
	vektor icoord, jcoord, tmpvertex, vertex[NUM];
	int ok, index[NUM];
	int vertexcount, vertexnum, edgesnum, atomic_area;
	int *edges;
	ivektor ivertex[NUM];
	real idist2, jdist2;
	real sin, cos, maxcos;
	int minj, ord[NUM];

	/* Allocate memory for data of edges */
	edges = (int *) malloc( neighnum * sizeof(int) );

	if( edges == NULL )
	error("Cannot allocate memory for vertices");

	vertexcount = 0;
	atomic_area = 0.0;

	/* Each possible vertex defined by the intersection of 2 lines is examined */
	for (i=0; i<(neighnum-1); ++i)
	{
		icoord.x = candcoord[i].x;
		icoord.y = candcoord[i].y;
		idist2 = -canddist2[i];

		for (j=i+1; j<neighnum; ++j)
		{
			jcoord.x = candcoord[j].x;
			jcoord.y = candcoord[j].y;
			jdist2 = -canddist2[j];

			det = icoord.x * jcoord.y - icoord.y * jcoord.x;

			/* check whether edges intersect */
			if ( SQR(det) > TOL2)
			{
				detinv = 1.0 / det;

				tmpvertex.x = ( icoord.y * jdist2 - jcoord.y * idist2 ) * detinv;
				tmpvertex.y = ( jcoord.x * idist2 - icoord.x * jdist2 ) * detinv;

				/* Check whether vertex belongs to voronoi cell */
				l = 0;
				ok = 1;

				do {
					if( l!=i && l!=j )
					ok = ( SPROD(candcoord[l] , tmpvertex ) <= canddist2[l] );

					++l;

				}while ( (ok==1) && (l<neighnum) );

				if( ok==1 )
				{
					ivertex[vertexcount].x = i;
					ivertex[vertexcount].y = j;

					vertex[vertexcount].x = 0.5 * tmpvertex.x;
					vertex[vertexcount].y = 0.5 * tmpvertex.y;

					++vertexcount;
				}
			}
		}
	}

	vertexnum = vertexcount;

	/* Check whether some vertices coincide */
	for ( i=0; i<vertexnum; ++i )
	{
		index[i] = 1;
		for (j=i+1; j<vertexnum; ++j )
		if ( (SQR(vertex[j].x-vertex[i].x)<TOL2) && (SQR(vertex[j].y-vertex[i].y)<TOL2) )
		index[i] = 0;
	}

	/* Remove coincident vertices */
	j = 0;
	for ( i=0; i<vertexnum; ++i )
	if ( index[i] != 0 )
	{
		ivertex[j].x = ivertex[i].x;
		ivertex[j].y = ivertex[i].y;

		vertex[j].x = vertex[i].x;
		vertex[j].y = vertex[i].y;

		++j;
	}

	vertexnum = j;

	/* Number of vertices of Voronoi cell must be greater than 2 */
	if ( vertexnum < 3 ) atomic_area = 0.0;
	else
	{
		/* Initialization */
		for (n=0; n<neighnum; ++n)
		edges[n] = 0;

		for (n=0; n<vertexnum; ++n)
		{
			++edges[ivertex[n].x];
			++edges[ivertex[n].y];
		}

		/* Number of edges of Voronoi cell */
		edgesnum = 0;
		for (n=0; n<neighnum; ++n)
		{
			edgesnum += edges[n];
		}

		edgesnum /= 2;

		/* Check whether number of vertices equals number of edges */
		if ( edgesnum == vertexnum )
		{

			/* Statistics */
			if( neighnum > maxneigh ) maxneigh = neighnum;
			sumneigh += neighnum;
			if( vertexnum > maxvert ) maxvert = vertexnum;
			sumvert += vertexnum;
			if( edgesnum > maxedges ) maxedges = edgesnum;
			sumedges += edgesnum;

			++atomcount;

			/* Order vertices */
			ord[0] = 0;
			for ( i=0; i<vertexnum-1; ++i)
			{
				maxcos = -1.0;
				for (j=0; j<vertexnum; ++j)
				{
					sin = vertex[j].x * vertex[ord[i]].y - vertex[j].y * vertex[ord[i]].x;
					if ( sin > TOL )
					{
						cos = SPROD( vertex[j], vertex[ord[i]] )/ sqrt(SPROD(vertex[j],vertex[j]))/ sqrt(SPROD(vertex[ord[i]],vertex[ord[i]]));
						if ( cos > maxcos )
						{
							maxcos = cos;
							minj = j;
						}
					}
				}

				ord[i+1] = minj;

			}

			/* Compute area of voronoi cell */
			for (i=0; i<vertexnum-1; ++i)
				atomic_area += 0.5 * (vertex[ord[i+1]].x * vertex[ord[i]].y - vertex[ord[i+1]].y * vertex[ord[i]].x);

			atomic_area += 0.5 * (vertex[ord[0]].x * vertex[ord[vertexnum-1]].y - vertex[ord[0]].y * vertex[ord[vertexnum-1]].x);

		} /* number of edges == number of vertices */

		else atomic_area = 0.0;

	} /* vertesnum < 3 */

	free(edges);

	return atomic_area;
}

#endif

