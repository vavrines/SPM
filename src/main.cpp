/***
 * Finite Volume Method for Linear Boltzmann Equation
 * (c) 2021 Tianbai Xiao
 *
 ***/

#include "common/config.hpp"
#include "common/io.hpp"
#include "solvers/solverbase.hpp"
#include <mpi.h>
#include <string>

int main( int argc, char** argv ) {
    // preprocess
    MPI_Init( &argc, &argv );
    std::string filename = ParseArguments( argc, argv );
    Config* config       = new Config( filename );

    // log
    PrintLogHeader( filename );

    // solver
    SolverBase* solver = SolverBase::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    // postprocess
    delete solver;
    delete config;
    MPI_Finalize();

    return EXIT_SUCCESS;
}
