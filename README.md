# Catalyst Adaptor
ParaView Catalyst adaptor example for a Fortran code

This package builds a library for visualising simulation data with a simple VTK visualisation pipeline. The pipeline can be defined either in C++ or using a Python script.

## Building the adaptor

To build this code, you will need to build and install the "Catalyst-v5.4.1-Base-Enable-Python-Essentials-Extras-Rendering-Base" source code package on your system.

Once this is done, build the code using CMake as follows:
```
mkdir build
cd build
cmake .. -DParaView_DIR=/path/to/catalyst/install/directory/lib/cmake/paraview-5.4
make
```
If you want to build a debug version of the code, add ```-DCMAKE_BUILD_TYPE=Debug``` to the CMake configuration, or use the ```ccmake``` configuration tool. You can add additional compiler flags using the ```-DCMAKE_CXX_FLAGS=``` option.

On a Cray XC50 system, the following build setup should work:
```
cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_EXE_LINKER_FLAGS=-dynamic -DParaView_DIR=/path/to/catalyst/install/directory/lib/cmake/paraview-5.4
```
Note that dynamic linking simplifies the linking process of the Fortran application significantly.

## Running a simulation with the adaptor

The Catalyst adaptor and libraries are usually dynamically linked. If the build system of your code does not hardcode shared library paths, you will need to set (possibly adapting ParaView version)
```
export LD_LIBRARY_PATH=/path/to/catalyst/installation/lib/paraview-5.4:$LD_LIBRARY_PATH
```
If you want to use the Python pipeline, set ```PYTHONPATH``` to something like
```
export PYTHONPATH=/path/to/catalyst/installation/lib/paraview-5.4/site-packages:/path/to/catalyst/installation/lib/paraview-5.4/site-packages/vtk:$PYTHONPATH
```

## Python scripts

The repository includes a number of Python scripts which define visualisation pipelines or provide some post-processing functionality.

### full_output.py
Simple Python pipeline for writing the model grid and data field to a VTK file.

### spherical_slice.py

Simple Python pipeline for creating spherical slices of model grid with a preset radius, which are written into a VTK polydata file. Full output of the model grid and data field can also be produced by setting the corresponding flag in the pipeline script.

### map_project.py

This Python program expects a spherical slice (as produced by the ```spherical_slice.py``` visualisation pipeline) in VTK polydata format as input and produces a VTK polydata file with a map projection as output. The program can handle partitioned datasets, but computing map projections for multiple timesteps is not supported yet.

Running
```
./map_project.py input.vtp output.vtp
```
computes a Mollweide map projection. Use flag ```--list-projections``` to get a list of projections and their short names (projections are provide by the PROJ library). Short names can be used to set another projection with the ```--projname``` flag, e.g., ```--projname=gall```.
