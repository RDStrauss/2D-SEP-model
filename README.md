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
- By default we've added placeholders for different spacecraft positions with a comments on the coordinate system used.
- The number of compute cores is hardcoded by `omp_set_num_threads(4)`.

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

## 7_feb_2010_benchmark

We benchmark directly with the 3D SDE results of *Droge et al. [2014]* for a specific event (7 Feb 2010 event) using the same input parameters. This comparison is also summarised in *Steyn et al. [2020]*. Two scenarios are shown: One with the spacecraft positions set to their actual position for the event, and a second where the magnetic field footpoint locations of each spacecraft is used. We suspect Droge et al. used the footpoint locations in their simulations to account for different solar wind speeds, and hence, different magnetic connectivity between the spacecraft and source when a longitudinally average solar wind speed is used. 

## 14_aug_2010_benchmark

We benchmark directly with the 3D SDE results of *Droge et al. [2016]* for a specific event (14 Aug 2010 event) using the same input parameters. Two scenarios are shown: One with the spacecraft positions set to their actual position for the event, and a second where the magnetic field footpoint locations of each spacecraft is used. We suspect Droge et al. used the footpoint locations in their simulations to account for different solar wind speeds, and hence, different magnetic connectivity between the spacecraft and source when a longitudinally average solar wind speed is used. 

## Example

The folder contains a quick multi-spacecraft example, including plotting routines for time profiles and contour plots.

## Quick convergence study

Still on the ToDO list, sadly....

## Disclosure and Notice

The  code  is published under the Creative Commons license, but is not intended to be used for commercial applications. We ask anyone using this model to reference Strauss & Fichtner [2015] in all research outputs and to contact the authors when used extensively.

## References

[Heita, P.K.N. 2019. Numerical investigation of solar energetic particle transport between the Sun, Earth, and Mars. *MSc thesis*. North-West University, South Africa.](https://dspace.nwu.ac.za/handle/10394/33865)

[Strauss, R.D. and Fichnter, H. 2015. On aspects pertaining to the perpendicular diffusion of solar energetic particles. *The Astrophysical Journal*, 801: 29.](https://ui.adsabs.harvard.edu/abs/2015ApJ...801...29S/abstract)

van den Berg, J.P., Strauss, R.D., and Effenberger, F. 2020. A primer on focused solar energetic particle transport. *Space Science Reviews*, (https://link.springer.com/article/10.1007/s11214-020-00771-x)

[W. Dröge, Y. Y. Kartavykh, N. Dresing, B. Heber, A. Klassen. 2014. Wide longitudinal distribution of interplanetary electrons following the 7 February 2010 solar event: Observations and transport modeling. *J. Geophys. Res.*, 119, 6074–6094](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014JA019933)

[W. Dröge, Y. Y. Kartavykh, N. Dresing, A. Klassen. 2016. MULTI-SPACECRAFT OBSERVATIONS AND TRANSPORT MODELING OF ENERGETIC ELECTRONS FOR A
SERIES OF SOLAR PARTICLE EVENTS IN AUGUST 2010. *ApJ*, 826, 134](https://iopscience.iop.org/article/10.3847/0004-637X/826/2/134)

[Ruhann Steyn, Du Toit Strauss, Frederic Effenberger and Daniel Pacheco4. 2020. The soft X-ray Neupert effect as a proxy for solar energetic particle injection: A proof-of-concept physics-based forecasting model. *JSWSC*, 10, 64](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014JA019933](https://www.swsc-journal.org/articles/swsc/full_html/2020/01/swsc200079/swsc200079.html)
