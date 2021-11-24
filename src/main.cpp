/* Finite Volume Method for Linear Boltzmann Equation */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL KITRT_ARRAY_API
#include <mpi.h>
#include <string>

#include "common/config.hpp"
#include "common/io.hpp"
#include "solvers/solverbase.hpp"

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );
    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Print input file and run info to file
    PrintLogHeader( filename );

    // Build solver
    SolverBase* solver = SolverBase::Create( config );

    // Run solver and export
    solver->Solve();
    solver->PrintVolumeOutput();
    delete solver;

    delete config;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
