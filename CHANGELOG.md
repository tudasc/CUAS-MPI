# CUAS-MPI Change Log

All notable changes to this project will be documented in this file.

## 0.2.0

### Changes in Model Behavior

- implements a consistent flux discretization scheme and a test (Langtangen&Linge,2017, doi: 10.1007/978-3-319-55456-3) with discontinuous transmissivity, but continuous flux (!222)
- adds more options for the initial head: 'Nzero', 'Nopc', 'low', 'mid', 'high', 'topg' or a valid number (!221)
- adds a steady state solution (Galluzzo2011, doi:10.17077/etd.kbfbts1p) with discontinuous transmissivity (!176)
- enables semi-implicit time stepping via command-line option 'timeSteppingTheta' (!196)
- accounts for the layer thickness consistently (pressure to head and head to pressure). Most notable in the output files, where now the effective pressure is zero when CUAS-MPI is initialized with --initialHead Nzero (!195).
- removes unintended transmissivity update for Dirichlet points (!175)
- renames ConstantForcing to SteadyForcing and TimeForcing to TimeDependentForcing (!194)
- adds new types of Forcing: ScalarTimeDependentForcing, BufferedForcing and MultiForcing, adds a new setup routine for Forcing (!178,!187,!212)
- refactors the interface of Forcing to be clearer (!210,!218)
- adds new class TimeIntegrator, wich handles the time steps in a well structured and extendable way (!204,!205)
- introduces interface WaterSource, which is used to define water in the model by NetCDF files and is used for coupling in the future (!215,!216)
- improves the output handling, which now distinguish between mutable and constant fields (!208)
- adds an option to provide a coordinates file which is also copied to the output NetCDF (!184)

### Softare Engineering

- applies our default class structure on e.g. CUASSolver and SolutionHandler, splits functions into subroutines and improves the function interfaces (!179,!190,!199,!202,!203,!220)
- verbose solver is now handled by the caller receiving a convergence information from solver (!183)
- removes month as a unit and improves parsing of time (!193)
- minor code improvements: fix clang-tidy issues, uses row, col instead of i,j, fix typos (!177,!191,!200)
- clarifies "How to install" and adds a list of projects using CUAS-MPI in README (!182)
- restructures and parallelizes the CI pipeline (!209)
- renames BUILD_TESTS to CUAS_ENABLE_TESTS (!226)
- adds two options to build docs using Doxygn or DoxygenSphinxBreathe. This is work in progress. (!113,!227)

### bug fix

- fix the last time step in cases of time overshooting (!197)
- fix the start time has not to be January 1st anymore (!211)
- fix the spacing of axis to ensure even axis spacing (!181)
- uses PETSC_NULLPTR for compatibility with future PETSc versions (!186)
- fix cmake deprecation warnings to ensure compatibility with future versions (!214)

## 0.1.0

  - first public version on [GitHub](https://github.com/tudasc/CUAS-MPI)
