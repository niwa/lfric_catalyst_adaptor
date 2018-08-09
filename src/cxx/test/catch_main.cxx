#define CATCH_CONFIG_RUNNER
#include "Catch2/catch.hpp"
#include "coprocessor.h"
#include <mpi.h>

int main( int argc, char * argv[] ) {

  MPI_Init(&argc, &argv);

  MPI_Fint mpi_fortran_comm = MPI_Comm_c2f(MPI_COMM_WORLD);

  // Only test C++ pipeline
  int usePythonPipeline = 0;
  char pythonScript[4] = "";

  // Set up C++ pipeline
  int visualisationFrequency = 2;
  char outputFileName[4] = "abc";

  coprocessor_initialize(visualisationFrequency, outputFileName,
                         mpi_fortran_comm, usePythonPipeline,
                         pythonScript);

  int result = Catch::Session().run( argc, argv );

  coprocessor_finalize();

  MPI_Finalize();

  return result;

}
