#define MAIN

#define LIB_PUBLIC __attribute__ ((visibility ("default")))

#ifdef __cplusplus
#define CEXT extern "C"
#else
#define CEXT
#endif



/* data structure to store a potential table or a function table */
typedef struct {
  double x,y,z;
  double mass, epot;
  int sort;
  int num;
} libimd_atom;

/**
 * Initializes and configures IMD to be used in an MC cycle.
 * Must be called exactly once before the very first time a MC cycle is performed
 * Afterwards the existing configuration is reused
 */
CEXT void LIB_PUBLIC libimd_init(char* parameterFile);

/**
 * Runs the simulation for an MC cycle.
 * The pointer atoms must hold n atoms that are to be inserted into the
 * MD domains. Position, mass and sort and num of each atoms must be properly defined.
 * Epot is ignored during import.
 *
 * After the particles are inserted, an simulation of the given numbers of steps in
 * the parameter file is performed.
 * The final atomic configurations is copied back into the given array, including updated
 * values of Epot and positions.
 *
 * Supported ensembles are NVE, NVT and GLOK. Multistep simulations are not supported.
 *
 * WARNING: The order of atoms in this list is not maintained
 */
CEXT int LIB_PUBLIC libimd_run(libimd_atom* atoms, int n);
CEXT int LIB_PUBLIC libimd_run_finiteTemp(libimd_atom* atoms, int n,
		real* energyDifference, int realToVirtualID, int virtualToRealID);
