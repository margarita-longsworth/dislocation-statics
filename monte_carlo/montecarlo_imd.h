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

typedef double real;

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

/**
 * Runs the simulation for an MC cycle for finite temperatures.
 * The pointer atoms must hold n atoms that are to be inserted into the
 * MD domains. Position, mass and sort and num of each atoms must be properly defined.
 * Epot is ignored during import.
 *
 * With a given ID one particle can be transformed from a virtual to a real particle and
 * another one from real to virtual. If the ID does not exist, not transformation is performed
 * This simulation run consists internally of multiple stpes.
 *
 * 1: Run several steps to equilibrate the system after the temperature is initialized
 * 2: Gather average potential energies of all real atoms over a number of steps
 * 3: Transmute selected particles between real<->virtual
 * 4: Run further equilibration steps
 * 5: Gather average potential energies of all real atoms over a number of steps
 *
 * The difference in energy computed in steps 2 and 5 is copied in the argument energyDifference
 *
 * The final atomic configurations is copied back into the given array, including updated
 * values of positions and the potential energy at the last step of the simulation.
 *
 * Supported ensembles are NVE, NVT. Multistep simulations are not supported.
 *
 * WARNING: The order of atoms in this list is not maintained
 */
CEXT int LIB_PUBLIC libimd_run_finiteTemp(libimd_atom* atoms, int n,
		real* energyDifference, int realToVirtualID, int virtualToRealID);
