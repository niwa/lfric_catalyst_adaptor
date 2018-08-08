#ifndef COPROCESSOR_H
#define COPROCESSOR_H

#include <mpi.h>

#include <vtkCPProcessor.h>
#include <vtkCPDataDescription.h>

//
// Coprocessor interface for a Fortran code
//

// We define our own C-style API rather than using Catalyst's "CPAdaptorAPI"
// layer, to get access to features such as custom MPI communicators.
//
// The API is used for talking to Catalyst, register our visualisation pipeline etc.
// Users can choose between a C++ pipeline and a scripted Python pipeline.

// We need to keep track of the coprocessor and simulation data states
extern vtkCPProcessor * Processor;
extern vtkCPDataDescription * dataDescription;

extern "C" {

  void coprocessor_initialize(const int visualisationFrequency, const char * outputFileName,
                              const MPI_Fint mpi_fortran_comm, const int usePythonPipeline,
                              const char * pythonScript);

  void coprocessor_requestdatadescription(int * timeStep, double * time,
                                          int * coprocessThisTimeStep);

  void coprocessor_needtocreategrid(int * needGrid);

  void coprocessor_coprocess();

  void coprocessor_finalize();

}

#endif
