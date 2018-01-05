# catalyst_adaptor
ParaView Catalyst adaptor example for a Fortran code

To build this code, you will need to build and install the "Catalyst-v5.4.1-Base-Essentials-Extras-Rendering-Base" source code package on your system.

Once this is done, build the code using CMake as follows:
```
mkdir build
cd build
cmake .. -DParaView_DIR=/path/to/catalyst/install/directory/lib/cmake/paraview-5.4
make
```
On a Cray XC50 system, the following build setup should work:
```
cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_EXE_LINKER_FLAGS=-dynamic -DParaView_DIR=/path/to/catalyst/install/directory/lib/cmake/paraview-5.4
```
Note that dynamic linking simplifies the linking process of the Fortran application significantly.
