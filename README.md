# hmmvp-okada
Construct an H-matrix and compute matrix-vector products of the form B*x, B(rs,:)*x, and B(rs,cs)*x(cs).
***

## Installation

This isn't a pretty way to do this but it's what I have at the moment

Just copy and paste these two blocks into a shell

```
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Hd.cpp -o src/Hd.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Compress.cpp -o src/Compress.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Hmat.cpp -o src/Hmat.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/HmatIo.cpp -o src/HmatIo.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/KeyValueFile.cpp -o src/KeyValueFile.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/CodeAnalysis.cpp -o src/CodeAnalysis.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Mpi.cpp -o src/Mpi.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/CHmat.cpp -o src/CHmat.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/SFHmat.cpp -o src/SFHmat.o
ar rucs lib/libhmmvp_omp.a src/Hd.o src/Compress.o src/Hmat.o src/HmatIo.o src/KeyValueFile.o src/CodeAnalysis.o src/Mpi.o src/CHmat.o src/SFHmat.o
```
and 
```
gfortran src/hmmvpbuild.cpp external/dc3omp.o src/Hd.o src/Compress.o src/Hmat.o src/HmatIo.o src/KeyValueFile.o src/CodeAnalysis.o src/Mpi.o src/CHmat.o src/SFHmat.o -lstdc++ -fopenmp -llapack -lblas -o bin/hmmvpbuild_omp
```
If you want to get rid of those pesky .o files, `make clean` should do the trick
