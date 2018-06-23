# Coarrays example
Example library utilising OpenCoarrays module.

## Compilation & running the program
To run the single-image version of the example program:
```
cmake CMakeLists.txt
make all
./Coarrays N MODE
```
This will start the program for NxN matrix. The MODE parameter specifies which algorithm should be used (0 for matrix multiplication, 1 for gaussian elimination).
This should work with both [gfortran](https://gcc.gnu.org/wiki/GFortran) or [ifort](https://software.intel.com/en-us/fortran-compilers) set as the defualt CMake Fortran compiler (preferably the former).<br />
</br>
To generate a python module compile the library and run the following:
```
F77=gfortran #fortran compiler of your choice
F90=gfortran
CC=gcc
f2py -c -m src.matrix_lib src/matrix_lib.F90
```
</br>
Doxygen config file included.

## Multiple images
```
caf $FORTRAN_FLAGS src/matrix_lib_coarr.F90 -o Matrix_Coarr
cafrun -n THREADS_NUMBER ./Matrx_Coarr N MODE
```
The multi-image version of the library (**src/matrix_lib.F90**) should be MPI compatible, albeit not tested.

## Results
For the detailed results check **results** directory.</br>
**Matrix multiplication**
For N ∈ <100, 2500> with step of 100, with -02 compilation flag.<br />
![mm](https://github.com/kasprzyckit/fortran-ca/blob/master/results/mm.png)
<br />
On the image above:
* black - **Single-image version**
* magneta - **Multi-image, one image**
* cyan - **Multi-image, two images**
* red - **Multi-image, three images**

**Gaussian elimination**
For N ∈ <100, 2500> with step of 100, with -02 compilation flag.
![ge](https://github.com/kasprzyckit/fortran-ca/blob/master/results/ge.png)
* black - **Single-image version**
* magneta - **Multi-image, one image**
* cyan - **Multi-image, two images**
* red - **Multi-image, three images**
