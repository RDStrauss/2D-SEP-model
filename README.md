# SEP_propagator-2D

This is the code for the spatially 2D solar energetic particle (SEP) transport model of *Strauss & Fichtner [2015]*

## Using SEP_propagator

There are a few parameters in the Fortran90 code which can be changed:
- `lambda`: The parallel radial mean free path; taken to be independent of radial distance.
- `lambda_perp`: The perpendicular mean free path; usually a constant factor of parallel
- `energy`: The particle's kinetic energy in MeV.
- `species`: 1 for electrons or 2 for protons.  
- `totaltime`: The time the simulation must be computed for in hours.
- For the adopted Reid-Axford injection profile, there is also an `acceleration_time` and `escape_time` in hours to control its shape.
- `injection_broadness`: The Gaussian width of the injected source.
- `phi_source`: The position of the peak injection (peak of Gaussian).

To compile the Fortran90 code, use either the Intel Fortran compiler with
```
gfortran -fopenmp -O2 Version_1.f90
```
and run the executable `./a.out`. In most cases the stack size should first be increased (or risk getting a segmentation fault),

```
ulimit -s unlimited
```

## Test 3 - SDE benchmark run

This folder contains the model comparison with the SDE-based results of *Droge et al. [2010]* for 4 MeV protons using a delta-like injection. Also contains results from the 1D SEP transport model. This is a special model set-up specifically to mimic a 1D model solution: (1) Perpendicular diffusion is set to zero, (2) a circumsolar injection in azimuth is hard-coded, (3) the radial dependence of the parallel mean-free-path is taken from *Droge et al. [2010]*, and (4) a delta-like injection in time is hardcoded. This is probably not the version to use when looking at any data.

## Example

Add example....

## Disclosure and Notice

The  code  is published under the Creative Commons license, but is not intended to be used for commercial applications. We ask anyone using this model to reference Strauss & Fichtner [2015] in all research outputs and to contact the authors when used extensively.

## References

[Heita, P.K.N. 2019. Numerical investigation of solar energetic particle transport between the Sun, Earth, and Mars. *MSc thesis*. North-West University, South Africa.](https://dspace.nwu.ac.za/handle/10394/33865)

[Strauss, R.D. and Fichnter, H. 2015. On aspects pertaining to the perpendicular diffusion of solar energetic particles. *The Astrophysical Journal*, 801: 29.](https://ui.adsabs.harvard.edu/abs/2015ApJ...801...29S/abstract)

van den Berg, J.P., Strauss, R.D., and Effenberger, F. 2020. A primer on focused solar energetic particle transport. *Space Science Reviews*, https://link.springer.com/article/10.1007/s11214-020-00771-x

