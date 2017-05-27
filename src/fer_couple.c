/* 

Communications routines for Fermi throught MPI 

Authors: Federico Caccia, 
Guido Giuntoli   (guido.giuntoli@bsc.es)

-> In this implementation Fermi can be connected with only one code

*/

#include <mpi.h>
#include "fermi.h"

#define  NTRIES_MAX 5

int fer_couple(int order)
{

    switch(order){
         
	 /****************************************************************************************************/
	 case COUPLE_INIT:

             fer_coinit(
	     break;

	 /****************************************************************************************************/
	 case COUPLE_RECV:

	     break;
	 
	 /****************************************************************************************************/
	 case COUPLE_SEND:

	     break;

	 /****************************************************************************************************/
	 case COUPLE_ENDS:

	     break;

	 /****************************************************************************************************/
    }

    return 0;
}

int fer_coinit(MPI_Comm * couple_comm, char * server_n)
{

    int  ntry, ierr;
    char port_n[MPI_MAX_PORT_NAME];

    // server_n is the string corresponding to the port we want to communicate
    // this name comes from input
    PetscPrintf(world," Stablishing MPI connection with %s:\n", server_n);

    // we look to the port name we want to communicate with
    ntry = 0;
    while(ntry < NTRIES_MAX){
	ierr = MPI_Lookup_name (server_n, MPI_INFO_NULL, port_n);
	if (ierr == 0){
	    PetscPrintf(world,"ntry:%d -> Can't find service, trying again\n", ntry);
	}
	ntry++;
    }
    if(ierr == 0){
	PetscPrintf(world,"Service found at port: %s after %d tries.\n", port_n, ntry);
	ntry = 0;
    }
    else{
	PetscPrintf(world,"fer_coinit: problem searching port name\n");
	return 1;
    }

    // connecting server's port 
    // doubt : MPI_COMM_SELF ?
    while(ntry < 5){
	PetscPrintf(world,"Connecting...\n");
	ierr = MPI_Comm_connect(port_n, MPI_INFO_NULL, 0, MPI_COMM_SELF, couple_comm);
	if(ierr == 0){
	    PetscPrintf(world,"Can't connect to service, re-trying.\n");
	}
	ntry = ntry + 1;
    }
    if (ierr != 0){
	PetscPrintf(world,"fer_coinit: problem connecting to port.\n");
	return 1;
    }
    PetscPrintf(world,"Connected:0K\n");


//    PetscPrintf(world,"Receiving general parameters.\n");
//
//    int tag = 100;
//    int status;
//  
//    ierr = MPI_Recv (&\textcolor{OliveGreen}{t_0}, 1, MPI_DOUBLE_PRECISION, 0, tag, Comp_Comm, &status);
//    ierr = MPI_Recv (&\textcolor{OliveGreen}{N_t}, 1, MPI_INTEGER, 0, tag, Comp_Comm, &status);
//    ierr = MPI_Recv (&\textcolor{OliveGreen}{N_input_var}, 1, MPI_INTEGER, 0, tag, Comp_Comm, &status);
//    ierr = MPI_Recv (&\textcolor{OliveGreen}{N_output_var}, 1, MPI_INTEGER, 0, tag, Comp_Comm, &status);

//    Check consistency in data
//
//    \textcolor{Gray}{Receiving control instruction:}
//    \textcolor{Gray}{0: restart step / 1: continue / 2: abort}
//
//
//    ierr = MPI_Recv (&\textcolor{OliveGreen}{order}, 1, MPI_DOUBLE_PRECISION, 0, tag, Comp_Comm, &status)

    return 0;
}
