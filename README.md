# :fire: InstaHeat :fire: Gravitational **Insta**bilities during Metric Pre**heat**ing

InstaHeat is a 3D pseudo-spectral preheating code including metric scalar perturbations to to analyze gravitational instabilities.

It evolves a single scalar field in the FLRW background metric including scalar metric perturbations and backreactions through a phase of metric preheating after the end of chaotic large field inflation. Due to parametric resonance the evolution becomes non-linear and the produced tensor perturbations are used to compute the power spectrum of the generated gravitational waves.

>:exclamation: **Disclaimer:** InstaHeat is scientific code under currently not under active development. Expect adventures. InstaHeat is written with speed in mind and, perhaps, not enough consideration for readability. However, we put quite some effort into a detailed [documentation](#build-documentation) and explanation of the [parameters and options](#parameters-and-options).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Build Process](#build-process)
  - [Dependencies](#dependencies)
  - [Build InstaHeat](#build-instaheat)
  - [Build Documentation](#build-documentation)
- [Physics & References](#physics-&-references)
  - [A Short Intro](#a-short-intro)
- [Numerics](#numerics)
- [Output](#output)
- [Parameters and Options](#parameters-and-options)
- [Previous Codes](#previous-codes)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Build Process

### Dependencies

The code has been tested under MacOS and Linux. It uses the _Fastest Fourier Transform in the West_ [FFTW3](http://www.fftw.org/), the GNU Scientific Library [GSL](https://www.gnu.org/software/gsl/), the [HDF5](https://www.hdfgroup.org/HDF5/) data model as well as the [OpenMP](http://openmp.org/wp/) API.

On a Mac probably the most convenient way to install current versions of FFTW3, GSL and HDF5 is via [homebrew](http://brew.sh/) (`h5utils` is not necessary, but provides useful tools):

```sh
brew install gsl hdf5 h5utils
brew install fftw --with-openmp
```

Your C compiler should implement OpenMP. If in doubt, install the GNU Compiler Collection [GCC](https://gcc.gnu.org/). Again, on a Mac simply type

```sh
brew install gcc
```

### Build InstaHeat

Navigate to the `src/` directory and simply try

```sh
make
```

If you get any errors you might have to modify the Library and Include search paths in the `Makefile`. We implemented a very simple check for the OS (only Linux and MacOS)

```make
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBS = -L/lib -L/share/sw/free/gsl/2.1/lib -g -shlib -lfftw3_threads -lfftw3 -lgsl -lgslcblas -lm -fopenmp
    CFLAGS += -I/share -I/share/sw/free/gsl/2.1/include
endif
ifeq ($(UNAME_S),Darwin)
    LIBS = -L/usr/local/lib -g -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -fopenmp
    CFLAGS += -march=native -I/usr/local/include
    export HDF5_CC = gcc-6
    export HDF5_CLINKER = gcc-6
endif
```

The `LIBS` and `CFLAGS` flags are simply what we used on your systems and you might want to change them. On MacOS `gcc` actually points to Apple's LLVM clang frontend, thus we have to tell `h5cc` explicitly to use the GNU compiler `gcc-6`. You might want to use another compiler.

Now `make` should do the job and produce an executable called `run` in the `src/` directory.

>:exclamation: You have to specify the parameters **before** the build process! Read the [Parameters and Options](#parameters-and-options) section before running InstaHeat.

### Build Documentation

InstaHeat is extensively annotated with [Doxygen](http://www.stack.nl/~dimitri/doxygen/) comments. A detailed documentation of the code can be exported to `html` and `latex`. To build the documentation for both, simply type

```sh
make doc
```

The documentation then resides in the created `doc` folder on the same level as the repository folder.

## Physics & References

### A Short Intro

InstaHeat starts at the end of chaotic large field inflation. After the inflaton rolled down to the bottom of the potential, it oscillates around the minimum with decaying amplitude. Due to cosmic inflation the vacuum fluctuations have been stretched out such that the inflaton field is almost homogeneous and can be approximated by a classical field with small perturbations on all scales. Thus typical initial conditions for InstaHeat are given by the Bunch-Davies vacuum.

We then evolve the inflaton field \f$\phi\f$ together with a scalar metric perturbation \f$\psi\f$ and the scale factor \f$a\f$ in the perturbed FLRW metric \f$ds^2 = - (1 + 2 \psi) dt^2 + a^2 \delta_{ij} (1 - 2 \psi) dx^i dx^j\f$. The initial perturbations \f$\psi\f$ are chosen to be compliant with the constraint from the Einstein equations.

Throughout the evolution we can compute the generated tensorial metric perturbations in the transverse traceless gauge alongside the the inflaton field and the metric perturbation. However, we only include backreactions from the scalar perturbations, not the tensorial part. From the tensor perturbations we compute the power spectru of the generated gravitational waves.

## Numerics

InstaHeat is written in C99. It is a pseudo-spectral code, i.e. we are constantly going back and forth between real space and Fourier space. Those Fourier transforms are expensive and take up roughly 80% to 85% of the runtime. On the other hand, it brings the advantage of exponentially small errors in spatial derivatives and thereby allows us to reach high precision and accuracy on comparatively small grids.

For the time integration of the equations of motion one can choose between several explicit Runge Kutta steppers. The available options are

* Runge Kutta 4 (classic RK4 with fixed time steps)
* Runge Kutta Fehlberg 4(5) (adaptive time steps)
* Runge Kutta Cash-Carp 4(5) (adaptive time steps)
* Runge Kutta Dormand-Prince 8(9) (adaptive time steps)
* **Runge Kutta Dormand-Prince 8(5,3), default** (adaptive time steps)

Details about the numerics can be found in the extensive documentation of the code under (tbd).

## Output

InstaHeat as flexible and extensive output options. All output of a run is bundled into a single `.h5` file. Additionally to basically all the parameters of the run, the scale factor \f$a\f$ and the physical time \f$t\f$ one can choose any combination of the following:

The full field (on part) of the simulation grid, the mean, the variance, the maximum, the minimum, the power spectrum for

* \f$\phi\f$ - The scalar inflaton field.
* \f$\dot{\phi}\f$ - The temporal derivative of the scalar inflaton field.
* \f$\psi\f$ - The scalar metric perturbation.
* \f$\dot{\psi}\f$ - The temporal derivative of the scalar metric perturbation.
* \f$\rho\f$ - The energy density.
* \f$p\f$ - The pressure.

Detailed information about the output options can be found in section [Parameters and Options](#parameters-and-options).

## Parameters and Options

The simulation parameters and output options are explained in great detail in the documentation of the parameters. This is also important for the build process. You can find the docs

* [here](doc_parameters.md) if you are viewing this on GitHub, or
* [here](@ref docparameters) if you are viewing this on Doxygen.

## Previous Codes

Previous codes include but are not limited to:

* LatticeEasy: [paper](http://arxiv.org/abs/hep-ph/0011159), [code](http://felderbooks.com/latticeeasy/index)
* ClusterEasy: [paper](http://arxiv.org/abs/0712.0813), [code](http://www.felderbooks.com/latticeeasy/index)
* CudaEasy: [paper](http://arxiv.org/abs/0911.5692)
* DEFROST: [paper](http://arxiv.org/abs/0809.4904), [code](http://www.sfu.ca/physics/cosmology/defrost/)
* PSpectRe: [paper](http://arxiv.org/abs/1005.1921), [code](http://cosmology.auckland.ac.nz/2011/10/16/pspectre/)
* HLattice: [paper](http://arxiv.org/abs/1102.0227), [code](http://www.cita.utoronto.ca/~zqhuang/hlat/)
* GABE: [paper](http://arxiv.org/abs/1305.0561), [code](http://cosmo.kenyon.edu/gabe.html)
* PyCool: [paper](http://arxiv.org/abs/1201.5029), [code](https://github.com/jtksai/PyCOOL)
