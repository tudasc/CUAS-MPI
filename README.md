# CUAS-MPI &middot; [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
## What is CUAS-MPI?
CUAS-MPI is an MPI parallel implementation of the Confined-Unconfined Aquifer System model ([CUAS18](#ref-CUAS-2018)) for subglacial hydrology.
The implementation relies on the parallel data structures and linear system solvers provided by the Portable, Extensible Toolkit for Scientific Computation ([PETSc](https://petsc.org/)).
The model uses a one-layer equivalent porous medium (EPM) approach for efficient (channel-like) and inefficient (cavity-like) subglacial water transport.
A two-dimensional Darcy-type groundwater flow equation with spatially and temporarily varying hydraulic transmissivity is solved considering confined and unconfined aquifer conditions.

## How to use it?

One of the integration tests can be used to generate a simple setup to explore the modelling choices and command line options.
The example below needs `ncks` and `ncap2` from the [NCO toolkit](https://nco.sourceforge.net/) to manipulate the NetCDF files.

```
# modifiy according to your installation
CUAS_BUILD_DIR=$CUAS_ROOT/CUAS-MPI/cmake-build-debug/

# number of MPI processes
NN=4

#
# generate simple input file from example integration test
#
exact=$CUAS_BUILD_DIR/test/cuascore/integration/test-exactSteadySolution
mpirun -n $NN $exact 1000.0 101 101 31 86400 out0.nc
# convert output to CUAS-MPI input and forcing file format
ncks -O -d time,-1 -v topg,bnd_mask,thk out0.nc input.nc
ncap2 -O -s "bmelt=watersource * 365*24*3600" out0.nc forcing.nc

# run a simple experiment 
#
cuas=$CUAS_BUILD_DIR/tools/cuas.exe

# set-up the solver
TOL="-ksp_rtol 1e-7 -ksp_atol 1e-15 -ksp_max_it 10000 -ksp_converged_use_min_initial_residual_norm"
export PETSC_OPTIONS="-options_left -ksp_initial_guess_nonzero -pc_type bjacobi -ksp_type gmres $TOL"

# make use of many options for this example
mpirun -n $NN  $cuas --totaltime '15 days' --dt '1 hour' --saveEvery 1 --verbose --outputSize large \
       --doChannels --Tmax 100  --Tmin 1.e-08 --initialHead Nzero  $opts \
       --conductivity 10 --layerThickness 0.1 \
       --flowConstant 3.4e-24 --cavityBeta 5e-4 --basalVelocityIce 1e-6 --supplyMultiplier 1.0 \
       --forcingFile forcing.nc  \
       input.nc output.nc
```

## How to install?

### Requirements

- c++ compiler
  - at least support for c++ 17
  - the code is tested using gcc/10.2.0
- MPI
  - we use OpenMPI
    - https://www.open-mpi.org/
    - version 4.0.x and 4.1.x
- PETSc
  - https://petsc.org/
  - we tested version 3.14.6
  - with MPI support
- netcdf-c
  - https://www.unidata.ucar.edu/software/netcdf/
  - we use version 4.7.4
  - with MPI and hdf5 1.8.22

### Build

We expect that *PETSC_ROOT* and *LIBNETCDF_ROOT* are set, being the path to the toplevel diretory of your PETSc and NetCDF installation.
And *CUAS_INSTALL_PREFIX* is the path which you want to use for the CUAS-MPI installation.
Fullfilling these requirements, you should be able to run:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DPETSC_DIR=$PETSC_ROOT -DNETCDF_DIR=$LIBNETCDF_ROOT
cmake --build build
cmake --install build --prefix CUAS_INSTALL_PREFIX
```

Beside the default cmake options like *CMAKE_BUILD_TYPE* we use the following options:
| Option | Default | Description                                                                                                  |
| --- | :---: |--------------------------------------------------------------------------------------------------------------|
| `BUILD_TESTS` | `OFF` | Enables build of unit and integration tests |

## References

<table style="border:0px">
<tr>
    <td valign="top"><a name="ref-CUAS-2018"></a>[CUAS18]</td>
    <td>Beyer, Sebastian and Kleiner, Thomas and Aizinger, Vadym and Rückamp, Martin and Humbert, Angelika
    <a href=https://doi.org/10.5194/tc-12-3931-2018>
    A confined–unconfined aquifer model for subglacial hydrology and its application to the Northeast Greenland Ice Stream</a>.
    In <i>The Cryosphere</i>, pages 3931–3947, 2018.</td>
</tr>
</table>
