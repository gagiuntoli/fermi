/* frees all the memory */

#include "funct.h"

int ferfini(void){
    // In case that the process was spawned by a parent, he is waiting for fermi to finish with an MPI_Barrier
    // With this order sended, the parent can read the output
    MPI_Comm parent;
    // Obtain an intercommunicator to the parent MPI job
    MPI_Comm_get_parent(&parent);
    // Check if this process is a spawned one and if so enter the barrier
    if (parent != MPI_COMM_NULL)
      MPI_Barrier(parent);

    SlepcFinalize();
    
    return 0;
}
