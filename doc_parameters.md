# Parameter file documentation (`parameters.sh`)

This is an explanation of the parameters in the file `parameters.sh`. We want to keep the parameter file itself short and concise to allow for quick changes. Therefore we do not clutter it by comments, but outsource the documentation to this file.

__IMPORTANT:__ Unlike parameters that are passed to a program at runtime, our parameters are read even _before_ compilation. The build process initiated by `make` (see `Makefile`) goes as follows:

1. From within the `Makefile` the executable script `configure` is called. It performs the following steps:
    * It loads the values specified in `parameters.sh`.
    * It makes a copy of `main_template.h` -- the template main header -- and saves it as `main.h`.
    * It replaces placeholders in `main.h` by the values found in `parameters.sh`. Hence, almost all parameters are provided as preprocessor defines in `main.h` before compilation.

    This step is not yet optimal and will probably be enhanced in the future.

2. The source files, i.e. the `.c` and `.h` files (including the newly generated `main.h`) are compiled and linked, see `Makefile` for details.

3. The final executable is called `run` and can be run with `./run`. It does not take any arguments anymore. All parameters have been determined during the build process.

This procedure has several advantages in terms of optimization and also usability:

* Everything is controlled from a single file `parameters.sh` that clearly lists all options and parameters. It contains solely the definition of the parameters and nothing else.
* For a given parameter file, only the necessary parts of the source are compiled. Thereby we can optimize memory usage and avoid superfluous allocations as well as unused or dead code.
* The runtime is optimized due to elimination of conditionals that determine program flow during runtime. Most of the program flow is already determined at compile time.
* Most parameters are available as fixed constant numbers before compilation starts. This information gives the compiler a lot of room for optimizations and all sorts of compiler magic.

## How to handle the `parameters.sh` file

* Do not delete or comment any parameters in `parameters.sh`. All of them are required and have to be specified.
* Do not change the names of any parameters.
* All parameter values are given as strings, i.e. enclosed by `" "`.
* Boolean values are encoded by `"1"` (true) and `"0"` (false).
* Numerical values can be given in any format that is known to C99. Additionally one can use `"PI"` and `"-PI"`.
* You can add comments to the file, starting with #. (It is a normal shell script.)
* The order of the parameters does not matter. We grouped parameters into the following categories by default:
    - __common__: These are changed frequently in different simulation runs. Feel free to play around with them.
    - __uncommon__: These are not changed as often. However, you might want to try some special things. Some of them are useful for debugging.
    - __rare__: There is hardly any reason to change these parameters. Only perform changes if you know exactly what you are doing and why you are doing it.
    - __output__: Here we specify the output of the simulation.

    These categories are only rough guidelines. In various cases you might well want to change _uncommon_ or rare_ values frequently.

In this documentation, we choose the following grouping of the parameters:

## Simulation volume

* `GRIDPOINTS_X`: (integer, >0) The number of spatial gridpoints in the x-direction
* `GRIDPOINTS_Y`: (integer, >0) The number of spatial gridpoints in the y-direction
* `GRIDPOINTS_Z`: (integer, >0) The number of spatial gridpoints in the z-direction
* `SPATIAL_LOWER_BOUND_X`: (double) Lower bound of the simulation box in the x-direction
* `SPATIAL_UPPER_BOUND_X`: (double) Upper bound of the simulation box in the x-direction
* `SPATIAL_LOWER_BOUND_Y`: (double) Lower bound of the simulation box in the y-direction
* `SPATIAL_UPPER_BOUND_Y`: (double) Upper bound of the simulation box in the y-direction
* `SPATIAL_LOWER_BOUND_Z`: (double) Lower bound of the simulation box in the z-direction
* `SPATIAL_UPPER_BOUND_Z`: (double) Upper bound of the simulation box in the z-direction

__Remarks__:

* The code can handle 1, 2 and 3 dimensional simulations. The dimension is implicitly deduced from the `GRIDPOINTS` in X, Y, and Z direction. See examples for further information. If there is only one gridpoint in a certain direction, the `SPATIAL_LOWER_BOUND` is used as a value for this point.

* The simulation always runs in an interval, rectangle or cuboid (depending on the dimension) with periodic boundary conditions in each direction.

* Integer powers of 2 (generally numbers that factor into small primes only) for the number of gridpoints result in the fastest code, due to the discrete Fourier transforms performed by FFTW3. In general any combination of numbers are supported. (One exception is when `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"` in which case the number of gridpoints has to be equal in each direction and the box length hast to be 10 in each direction. See [initial conditions](#initial-conditions) for details.) However, we only extensively tested equal numbers of gridpoints in each direction, all being integer powers of 2.

* The number of spatial gridpoints is also the number of Fourier modes used in the spectral parts of the code. (There is no padding in momentum space, hence one has to choose sufficiently many gridpoints to avoid aliasing. See also the section on filtering for more information.)

* The `SPATIAL_UPPER_BOUND` must be larger than `SPATIAL_LOWER_BOUND` for each direction.

* For the simulation usually only the box length in each direction (i.e. `SPATIAL_UPPER_BOUND - SPATIAL_LOWER_BOUND`) is relevant. However, the actual spatial gridpoint coordinates are required, if the initial conditions are constructed from internal functions. See [initial conditions](#initial-conditions) for the construction of initial conditions.

__Examples__:

* For a one-dimensional simulation with `N` gridpoints, one has to specify
    ```C
    GRIDPOINTS_X="N"
    GRIDPOINTS_Y="1"
    GRIDPOINTS_Z="1"
    ```

* For a two-dimensional simulation with `N` times `M` gridpoints, one has to specify
    ```C
    GRIDPOINTS_X="N"
    GRIDPOINTS_Y="M"
    GRIDPOINTS_Z="1"
    ```

## Time evolution

* `INITIAL_TIME`: (double) The starting time of the simulation. Since the evolution equations do not explicitly depend on time, this can be chosen arbitrarily. We recommend leaving it at 0.0.

* `FINAL_TIME`: (double, >`INITIAL_TIME`) The final time of the simulation. This is obviously one of the key drivers for the overall runtime of the simulation.

* `DELTA_T`: (double, >0) The initial value for the time step. Its meaning depends on the chosen integration routine [program flow](#program-flow).
    - If the evolution is performed by the fixed time step RK4 routine (see `rk4.c`) this is the fixed time step used throughout the simulation (see remarks for one exception).
    - If the evolution is performed by the adaptive time step Dormand Prince 8(5,3) or the Runge Kutta Felberg 4(5) routine (see`dopri853.c` and `rkf45.c` respectively),the time step is adjusted after each step. Thus `DELTA_T` is just the initial try for the first step.

* `MINIMAL_DELTA_T`: (double, >0) This is only relevant if the chosen integration routine has adaptive step sizes (i.e. `dopri853` or `rkf45`), see [program flow](#program-flow). It gives a lower bound on the step size. Once the integration routine would have to reduce the step size below this limit to satisfy the tolerances, the integration is terminated.

* `MAX_STEPS`: (integer, >0) The maximal number of steps performed by the integration routine. If (`FINAL_TIME` - `INITIAL_TIME`) / `DELTA_T` is larger than `MAX_STEPS`, the integration will not start.

* `MAX_DT_HUBBLE_FRACTION`: (double, >0) This is only relevant if the chosen integration routine is the Dormand Prince 8(5,3) stepper (adaptive stepsize) [program flow](#program-flow). To avoid large timesteps in the beginning of the evolution that could lead to instabilities, we limit the time step from above by `MAX_DT_HUBBLE_FRACTION` times the Hubble time $$1/H$$. Typical values are on the order of $$10^{-3}$$ to $$10^{-2}$$.

__Remarks__:

* Even for the fixed time step RK4 method, the time step might be different for the very last step. To ensure that the simulation always ends exactly at the specified `FINAL_TIME`, the time step might be adjusted for the very last step.

* When using one of the adaptive time step integrators (`dopri853` or `rkf45`), we recommend a very small `DELTA_T`. The routines will quickly increase the time steps, if they are sufficiently small initially, such that the code does not have to do much extra work. However, if the initial time step is too large, small errors in the first steps can lead to a seriously flawed evolution later on.

* Technically, `MAX_STEPS` is only a sensible strict bound on the number of time steps for the fixed time step RK4 (`rk4`) routine. However, it is also used for the adaptive step size routines as a proxy whether the parameters make sense. In case one uses a very small initial `DELTA_T` which is expected to grow by a few orders of magnitude, `MAX_STEPS` can easily abort the integration and should thus be chosen very large.

## Initial conditions

* `A_INITIAL`: (double, >0) The initial value for the scale factor a, see TODO(link equation)

* `INITIAL_CONDITIONS`: This parameter determines how to obtain the initial conditions for the fields $$\phi$$, $$\dot{\phi}$$ and in certain cases also $$\psi$$ and $$\dot{\psi}$$. The valid options are:
    - `"IC_FROM_INTERNAL_FUNCTION"`: The initial $$\phi$$ and $$\dot{\phi}$$ are given as functions `phi_init` and `dphi_init` in `setup.c`. These are evaluated on the spatial grid specified by [simulation volume](#simulation-volume). The initial values for $$\psi$$ are computed via `mk_psi` in `toolbox.c`.
    - `"IC_FROM_H5_FILE"`: One can load initial conditions from a previous simulation. The following conditions have to be satisfied:
          + The `INITIAL_DATAPATH` parameter in section [initial conditions](#initial-conditions) points to the output file of a previous simulation from which a certain timeslice should be used as initial conditions for the current run.
          + The `INITIAL_TIME` parameter in section [initial conditions](#initial-conditions) has to lie in the simulation region of the previous output file, i.e. `INITIAL_TIME` has to be larger than `INITIAL_TIME` of the previous simulation and smaller than `FINAL_TIME` of the previous simulation.
          + The number of gridpoints (`GRIDPOINTS_{X,Y,Z}` parameters in section [simulation volume](#simulation-volume)) have to coincide with the number of gridpoints from the previous simulation.

      Then the program will find the time slice _closest_ to the specified `INITIAL_TIME` in the previous output file and use all necessary values ($$a$$, $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$) as initial conditions for the new simulation. Note that the simulation might therefore not start _precisely_ at the specified `INITIAL_TIME`.

    - `"IC_FROM_DAT_FILE_WITH_PSI"` and `"IC_FROM_DAT_FILE_WITHOUT_PSI"`: One can load initial conditions from a `.dat` file. This feature was implemented during development of the code when initial conditions were provided from a collaborator in a `.dat` file. There is an option to also import $$\psi$$ and $$\dot{\psi}$$ from this `.dat` file, or to compute them from the provided $$\phi$$ and $$\dot{\phi}$$ via `mk_psi` in `toolbox.c`. If you want to read initial conditions from a separately generated file, you might change the method `read_initial_data` in `io.c` accordingly. Also make sure to read and understand the memory layout of the field variables TODO(link, where?).
    - `"IC_FROM_BUNCH_DAVIES"`: The initial conditions are given by the Bunch Davies vacuum. We followed the implementation of DEFROST TODO(link, ref) to compute initial values for $$\phi$$ and $$\dot{\phi}$$. For detailed information see TODO(link paper, coderef). Subsequently we compute $$\psi$$ and $$\dot{psi}$$ via `mk_psi` in `toolbox.c`.

* `INITIAL_DATAPATH`: See section [file IO](#file-io).

* `INITIAL_TIME`: See section [simulation volume](#simulation-volume).

## Potential parameters

* `MASS`: (double, >0) The mass parameter used in the code. Note that this value can be easily rescaled, hence has little physical relevance.

* `INFLATON_MASS`: (double, >0) If `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"`, this value is the inflaton mass in units of the Planck mass. While the actual mass parameter used in the code (for example when computing the potential) can be rescaled to arbitrary values, `INFLATON_MASS` is a physical quantity and only enters once when setting the amplitude of initial fluctuations in the Bunch Davies vacuum. Note that in the code $$8 \pi G = 1$$, hence we set the reduced Planck mass to one. `INFLATON_MASS` is the only mass related parameter in the code that carries physical meaning.

* `MASS_KARSTEN`: (double, >0) Obsolete, intended for developmental use only. This parameter was usde to scale certain amplitudes for comparison with Karsten Jedamzik's code.

## File IO

* `DATAPATH`: (string) The full path to the output file. All output of a single run is bundled into one `.h5` file. The `DATAPATH` can be relative to the working directory. It also needs to include the filename and `.h5` extension. The file does not have to exist. It will be created. If it does already exist, it will be __overwritten__. __Example__: `DATAPATH="../data/myrun.h5"`.

* `INITIAL_DATAPATH`: (string) The full path to the file that contains initial conditions for the simulation.
    - If `INITIAL_CONDITIONS="IC_FROM_DAT_FILE_WITH_PSI"` or `INITIAL_CONDITIONS="IC_FROM_DAT_FILE_WITHOUT_PSI"`, this is the full path (including filename and `.dat` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.
    - If `INITIAL_CONDITIONS="IC_FROM_H5_FILE"`, this is the full path (including filename and `.h5` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.

* `VERSION_CONTROL`: If git or mercurial is used as a version control system, one can add a hash of the current revision to the output file. This helps identify the data later. This option was most useful during development. __Warning__: Since the program has to create a pipe, invoke a shell and call git or mercurial as well as receive the output, this is not guaranteed to work on different systems. Choose the `VERSION_CONTROL_NONE` option when in doubt. Valid options are
  - `"VERSION_CONTROL_HG"`: Try to get the hash of the current mercurial revision and include it in the output.
  - `"VERSION_CONTROL_GIT"`: Try to get the hash of the current git revision and include it in the output.
  - `"VERSION_CONTROL_NONE"`: Do not try to get the hash of the current revision. Nothing is added to the output file. __Note__: This is the recommended choice, because it is quite possible that the program does not compile or fail at runtime on different systems.

* `WRITE_OUT_BUFFER_NUMBER`: (integer, >0) The program uses an internal buffer where output is accumulated before it is written to disk. Since rare and large writeouts are usually faster than frequent small ones, this gives a speedup despite the memory overhead and copy instructions to the buffers. `WRITE_OUT_BUFFER_NUMBER` specifies how many timeslices of the output should be buffered before disk access. If one outputs full 3D fields ([see output](#output)), memory constraints might not allow for a large buffer or slow down the simulation.

* `POWER_SPECTRUM_BINS`: (integer, >0) The number of bins used to compute the power spectra. This is relevant for the gravitational wave spectrum as well as for any of the outputs (see [output](#output)) ending in `POWER_SPECTRUM`. In the computation of the power spectrum, modes are binned into this number of bins. See TODO(link thesis) for details.

* `TIME_STEP_SKIPS`: (integer, >0) The number of time steps skipped between outputs. To avoid large output files, one can skip `TIME_STEP_SKIPS` many time slices, before writing a time slice to disk again.

* `STRIDE_X`: (integer, >0) The stride in the x-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than are used for the computation internally.

* `STRIDE_Y`: (integer, >0) The stride in the y-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than are used for the computation internally.

* `STRIDE_Z`: (integer, >0) The stride in the z-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than are used for the computation internally.

__Remarks__:

* Given the number of gridpoints $$n$$ in the x-direction by `GRIDPOINTS_X` and the stride $$s$$ in the x-direction by `STRIDE_X`, the number of gridpoints in the output is computed by $$(n+s-1)/s$$ rounded down towards the next integer. The output gridpoints are evenly spaced within the number of gridpoints used for internal computations.

__Examples__:

* If we specify

```C
GRIDPOINTS_X="256"
GRIDPOINTS_Y="200"
GRIDPOINTS_Z="59"

STRIDE_X="4"
STRIDE_Y="1"
STRIDE_Z="6"
```

then the computation would be done on a $$256 \times 200 \times 59$$ grid. However the output of the fields would be on a $$64 \times 200 \times 10$$.

## Performance parameters

* `THREAD_NUMBER`: (integer, >=0) The number of threads used in OpenMP and also for the discrete Fourier transforms performed by FFTW3. __Important__: If `THREAD_NUMBER="0"`, the number of threads is determined automatically by `omp_get_max_threads()`. The default value is 0.

* `FFTW_DEFAULT_FLAG`: The planning flag used in the discrete Fourier transforms performed by FFTW3. Take a look at the FFTW3 documentation TODO(link) for more. The valid options are:
    - `"FFTW_ESTIMATE"`
    - `"FFTW_MEASURE"`
    - `"FFTW_PATIENT"`
    - `"FFTW_EXHAUSTIVE"`

## Special paramaters for the Dormand Prince 8(5,3) integrator

The parameters in this section are only relevant if one of the adaptive time step size integration routine is used, i.e. if `INTEGRATION_METHOD="DOPRI853"` or `INTEGRATION_METHOD="RKF45"`.

* `RELATIVE_TOLERANCE`: The relative tolerance in the integration routine. See TODO(link numerical recipes, paper?) for more information.

* `ABSOLUTE_TOLERANCE`: The absolute tolerance in the integration routine. See TODO(link numerical recipes, paper?) for more information.

* `SMALLEST_SCALING`: The smallest possible rescaling of the step size in any step.

* `LARGEST_SCALING`: The largest possible rescaling of the step size in any step.

* `BETA`: This is an internal parameter of the adaptive step size control for PI control of the step size following the rescaling method of Lund. Only change this if you know exactly what you are doing. See TODO(link numerical recipes, Gustaffson paper) for more information.

* `SAFE`: This is an internal safety parameter of the adaptive step size control, which makes the next step size more likely to be accepted. Only change this if you know exactly what you are doing. See TODO(link numerical recipes, hairer book) for more information.

## Program flow

* `INTEGRATION_METHOD`: There are three integration routines available:
    - `"RK4"`: The standard fourth order Runge Kutte stepper with fixed time step size.
    - `"RKF45"`: A fourth order Runge Kutta method with 5th order error estimation for adaptive time stepping. We are using the GSL implementation TODO(link to gsl).
    - `"DOPRI853"`: A more sophisticated adaptive Dormand Prince stepper of 8th order with 5th and 3rd order error estimation for adaptive time stepping. A detailed description can be found in TODO(link numerical recipes).

* `ENABLE_FFT_FILTER`: Switch for a spectral filter for the fields. If switched on, at each time step the highest modes of $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$ are cut off to avoid aliasing. A more detailed description can be found in TODO(link to thesis).

## Miscellaneous

* `SEED`: The seed for the random number generation used for example to create the Bunch Davies vacuum if `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"` or also for `INITIAL_CONDITIONS="IC_FROM_INTERNAL_FUNCTION"` (depending on what the functions do). We use the GSL implementation of the Mersenne Twister.

* `FFTW_SIMD_STRIDE`: To understand what this parameter is doing, we need to go into some implementation details. As a general rule, don't change it unless you get a runtime error telling you to double the parameter. In this case, double the value and try again. _Detailed explanation_: The field that is evolved by the integration routines bundles all relevant physical quantities in one large array. Thereby we can provide a simple interface to `mk_rhs()` in `toolbox.c` and do not need unnecessarily many similar statements in the update steps and temporary arrays for the integration routines. This means that we have to index into those large arrays to access the individual physical quantities and especially also perform Fourier transforms on parts of the large array starting at various indices. For speed and simplicity, we only set up one plan for Fourier transforms and one for inverse Fourier transforms and reuse it for different arrays or index positions within arrays. For FFTW3 to be able to take advantage of SIMD instructions, we need to ensure, that all memory which we did not explicitly create a plan for is aligned properly. Depending on `PSI_METHOD` and the grid dimensions, proper alignment might be broken. Then `FFTW_SIMD_STRIDE` allows us to squeeze in some irrelevant memory dragged along in the evolution to restore proper alignment.

## Output

Most of the simulation parameters are always present in the output (and cannot be switched off). The optional values are

* `PHI`: The scalar field $$\phi$$.

* `DPHI`: The temporal derivative of the scalar field $$\dot{\phi}$$.

* `PSI`: The scalar metric perturbation $$\psi$$.

* `DPSI`: The scalar metric perturbation $$\dot{\psi}$$.

* `RHO`: The energy density $$\rho$$.

* `PHI_SUMMARY`: The mean value, the variance, the minimum and the maximum of the scalar field $$\phi$$ on each timeslice (in this order).

* `DPHI_SUMMARY`: The mean value, the variance, the minimum and the maximum of the temporal derivative of the scalar field $$\dot{\phi}$$ on each timeslice (in this order).

* `PSI_SUMMARY`: The mean value, the variance, the minimum and the maximum of the scalar metric perturbation $$\psi$$ on each timeslice (in this order).

* `DPSI_SUMMARY`: The mean value, the variance, the minimum and the maximum of the scalar metric perturbation $$\dot{\psi}$$ on each timeslice (in this order).

* `RHO_SUMMARY`: The mean value, the variance, the minimum and the maximum of the energy density $$\rho$$ on each timeslice (in this order).

* `PHI_POWER_SPECTRUM`: The power spectrum of the field $$\phi$$.

* `DPHI_POWER_SPECTRUM`: The power spectrum of the field $$\dot{\phi}$$.

* `PSI_POWER_SPECTRUM`: The power spectrum of the field $$\psi$$.

__Remarks__:

* `PHI`, `DPHI`, `PSI`, `DPSI`, `RHO` result in `xout*yout*zout` double values where {x,y,z}out are determined from `GRIDPOINTS_{X,Y,Z}` and `STRIDE_{X,Y,Z}` (see [file IO](#file-io)) on __each__ time slice of the output (per field).

* All values ending in `SUMMARY` result in four double values on each time slice of the output (per field).

* All values ending in `POWER_SPECTRUM` result in `POWER_SPECTRUM_BINS` double values on each time slice of the output (per field).
