#include <iostream>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "nsgaii_options.h"
#include "nsgaii_core.h"


int main(int argc, char* argv[])
{

    boost::mpi::environment env;
    boost::mpi::communicator world;

    std::string value;
      if (world.rank() == 0) {
        value = "Hello, World!";
      }

      boost::mpi::broadcast(world, value, 0);

      std::cout << "Process #" << world.rank() << " says " << value
                << std::endl;

      return 0;

    NSGAII::options opt;

    if (world.rank() == 0){
        NSGAII::readInputParameters(argc, argv,opt);
        NSGAII::Population pop(opt);
        pop.initializePopulation();

    }
    world.barrier();

    while (true){
        if (world.rank() == 0){
            //boost::mpi::broadcast(world, opt.PopulationSize,0);
        }



        break;
    }

    std::cout << "Process " << world.rank() << " says " << opt.PopulationSize << std::endl;







/*
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    std::cout << "Process " << world_rank << " of " << world_size << std::endl;


    // Finalize the MPI environment.
    MPI_Finalize();
    */

    return 0;
}
