# RWPT_unsat
The folders "/1D_homogeneous", "/1D_stratified" and "/3D_heterogeneous" provide the inputs and output files required to reproduce the corresponding RWPT simulations in unsaturated conditions. 
For each transport simulation the file "rw3d.inp" set the input parameters up. The file "rw3d.nam" set the name of input and output files up. 

The folder "/source" provides the source files of the new implementation of the random-walk particle tracking code RW3D. 
To run it, the code has to be re-compiled on any PC, preferentially using intel fotran compilers. The Lapack library has to be linked. 
