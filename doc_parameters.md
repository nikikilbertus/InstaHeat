# Parameter file documentation (`parameters.sh`)

This is an explanation of the parameters in the file `parameters.sh`. We want to keep the parameter file itself short and concise to allow for quick changes. Therefore we do not clutter it by comments, but outsource the documentation to this file.

__IMPORTANT:__ Unlike parameters that are passed to a program at runtime, our parameters are read even _before_ compilation. The build process initiated by `make` (see `Makefile`) goes as follows:

1. From within the `Makefile` the executable script `configure` is called. It performs the following steps:
    * It loads the values specified in `parameters.sh`.
    * It makes a copy of `main_template.h` - the template main header - and saves it as `main.h`.
    * It replaces placeholders in `main.h` by the values found in `parameters.sh`. Hence, almost all parameters are provided as preprocessor defines in `main.h` before compilation.

>This step is not yet optimal and will probably be enhanced in the future.

2. The source files, i.e. the `.c` and `.h` files (including the newly generated `main.h`) are compiled and linked, see `Makefile` for details.

3. The final executable is called `run` by default and can be run with `./run`. It does not take any arguments anymore. All parameters have been determined during the build process.

This procedure has several advantages in terms of optimization and usability:

* All settings are controlled from a single file `parameters.sh` that clearly lists all options and parameters. It contains solely the definition of the parameters and is thus short and concise.
* For a given parameter file, only the necessary parts of the source are compiled. Thereby we can optimize memory usage and avoid superfluous allocations as well as unused or dead code.
* The runtime is optimized due to elimination of conditionals that determine program flow during runtime. A large portion of the program flow is already determined at compile time.
* Most parameters are available as fixed constant numbers during compilation. This information gives the compiler a lot of room for optimizations and all sorts of compiler magic.

## How to handle the `parameters.sh` file

* Do not delete or comment any parameters in `parameters.sh`. All of them are required and have to be specified.
* Do not change the names of any parameters.
* All parameter values are given as strings, i.e. enclosed by `" "`.
* Boolean values are encoded by `"1"` (true) and `"0"` (false).
* Numerical values can be given in any format that is known to C99. Additionally one can use `"PI"` and `"-PI"`.
* You can add comments to the file, starting with `#`. (It is handled as a regular shell script.)
* The order of the parameters does not matter. We grouped parameters into the following categories based on our usage during debugging and running our simulations:
    - __common__: These are changed frequently in different simulation runs. Feel free to play around with them.
    - __uncommon__: These are not changed as often. However, you might want to try some special scenarios. Some of them are useful for debugging.
    - __rare__: There is hardly any reason to change these parameters. Only perform changes if you know exactly what you are doing and why you are doing it.
    - __output__: All output specifications are bundled in a separate section.

    These categories are only rough guidelines. In various cases you might want to change _uncommon_ or _rare_ values also.

In this documentation, we choose another categorization of the parameters:

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

* The code can handle 1, 2 and 3 dimensional simulations. The dimension is implicitly deduced from the `GRIDPOINTS` in X, Y, and Z direction. See examples for further information. If there is only one gridpoint in a certain direction, the `SPATIAL_LOWER_BOUND` is used as a value for this point. Some parts of the code require `GRIDPOINTS_X >= GRIDPOINTS_Y >= GRIDPOINTS_Z` and we recommend to always stick to this convention.

* The simulation always runs in an interval (1D), rectangle (2D) or cuboid (3D) with periodic boundary conditions in each direction.

* Integer powers of 2 (generally numbers that factor into small primes only) for the number of gridpoints result in the fastest code, due to the discrete Fourier transforms performed by FFTW3. In general any combination of numbers is supported. However, we only extensively tested equal numbers of gridpoints in each direction, all factoring to powers of 2 and 3 only.

* The number of spatial gridpoints is also the number of Fourier modes used in the spectral parts of the code. (There is no padding in momentum space, hence one has to choose sufficiently many gridpoints to avoid aliasing. See also the section on filtering for more information.)

* The `SPATIAL_UPPER_BOUND` must be larger than `SPATIAL_LOWER_BOUND` for each direction.

* For most parts of the simulation only the box length in each direction (i.e. `SPATIAL_UPPER_BOUND - SPATIAL_LOWER_BOUND`) is relevant. However, the actual spatial gridpoint coordinates are of importance, if the initial conditions are constructed from internal functions. See [initial conditions](#initial-conditions) for the construction of initial conditions.

__Examples__:

* For a one-dimensional simulation with `N` gridpoints, one has to specify
    ```C
    GRIDPOINTS_X="N"
    GRIDPOINTS_Y="1"
    GRIDPOINTS_Z="1"
    ```

* For a two-dimensional simulation with `N` times `M` gridpoints, one has to specify (with `N >= M`)
    ```C
    GRIDPOINTS_X="N"
    GRIDPOINTS_Y="M"
    GRIDPOINTS_Z="1"
    ```

## Time evolution

* `INITIAL_TIME`: (double) The starting time of the simulation. Since the evolution equations do not explicitly depend on time, this can be chosen arbitrarily. We recommend leaving it at 0.0.

* `FINAL_TIME`: (double, >`INITIAL_TIME`) The final time of the simulation. This is obviously one of the key drivers for the overall runtime of the simulation. A reasonable final time typically depends on various other parameter values and is found by trial and error.

* `DELTA_T`: (double, >0) The initial value for the time step. Its meaning depends on the chosen integration routine [program flow](#program-flow).
    - If the evolution is performed by the fixed time step RK4 routine (see `rk4.c`) this is the fixed time step used throughout the simulation (see remarks for one exception).
    - If the evolution is performed by one of the other adaptive time step routines `RKF45`, `RKCK`, `DOPRI89`, `DOPRI853` (see `dopri853.c` and `gsl_stepper.c` respectively),the time step is adjusted after each step. Thus `DELTA_T` is just the initial try for the first step.

* `MINIMAL_DELTA_T`: (double, >0) This is only relevant if the chosen integration routine has adaptive step sizes, i.e. `RKF45`, `RKCK`, `DOPRI89` or `DOPRI853`, see [program flow](#program-flow). It gives a lower bound on the step size. Once the integration routine would have to reduce the step size below this limit to satisfy the tolerances, the integration is terminated.

* `MAX_STEPS`: (integer, >0) The maximal number of steps performed by the integration routine.

* `MAX_RUNTIME`: (double, >0) The maximal wall clock runtime of the program in seconds. This is a handy termination criterion when running InstaHeat on a cluster with limited run time per job. For example, setting `MAX_RUNTIME="172500"` will automatically wrap up the program and flush all buffers to disk shortly before the simulation has run 48 hours. If a run is cancelled by the scheduler due to a time limit, the output file might be corrupt and unreadable, because it has not been closed properly. You can set it to a huge value if you do not want to restrict the simulation by a certain runtime.

* `MAX_DT_HUBBLE_FRACTION`: (double, >0) This is only relevant if the chosen integration routine is the Dormand Prince 8(5,3) stepper (`DOPRI853`), see [program flow](#program-flow). To avoid large timesteps in the beginning of the evolution that could lead to instabilities, we limit the time step from above by `MAX_DT_HUBBLE_FRACTION` times the Hubble time $$1/H$$. Typical values are on the order of $$10^{-3}$$ to $$10^{-2}$$.

__Remarks__:

* Even for the fixed time step `RK4` method, the time step might be different for the very last step. To ensure that the simulation always ends exactly at the specified `FINAL_TIME`, the time step might be adjusted for the very last step.

* When using one of the adaptive time step integrators (`RKF45`, `RKCK`, `DOPRI89` or `DOPRI853`), we recommend a very small `DELTA_T`, e.g. $$10^{-6}$$. The routines will quickly increase the time steps, if they are sufficiently small initially such that the code does not have to do much extra work. However, if the initial time step is too large, small errors in the first steps can lead to a seriously flawed evolution later on.

## Initial conditions

* `A_INITIAL`: (double, >0) The initial value for the scale factor $$a$$.

* `INITIAL_CONDITIONS`: This parameter determines how to obtain the initial conditions for the fields $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$ and $$a$$. The valid options are:
    - `"IC_FROM_INTERNAL_FUNCTION"`: The initial $$\phi$$ and $$\dot{\phi}$$ are given as functions `phi_init` and `dphi_init` in `setup.c`. These are evaluated on the spatial grid specified by [simulation volume](#simulation-volume). The initial values for $$\psi$$ are computed via `mk_initial_psi` in `setup.c`.
    - `"IC_FROM_H5_FILE"`: One can load initial conditions from a previous simulation. The following conditions have to be satisfied:
          + The `INITIAL_DATAPATH` parameter in section [initial conditions](#initial-conditions) points to the output file of a previous simulation from which a certain timeslice should be used as initial conditions for the current run.
          + The `INITIAL_TIME` parameter in section [initial conditions](#initial-conditions) has to lie in the simulation region of the previous output file, i.e. `INITIAL_TIME` has to be larger than `INITIAL_TIME` of the previous simulation and smaller than `FINAL_TIME` of the previous simulation.
          + The number of gridpoints (`GRIDPOINTS_{X,Y,Z}` parameters in section [simulation volume](#simulation-volume)) have to coincide with the number of gridpoints from the previous simulation. In general, the spatial volume simulation parameters have to compatible.

      Then the program will find the time slice _closest_ to the specified `INITIAL_TIME` in the previous output file and use all necessary values ($$a$$, $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$, $$t$$) as initial conditions for the new simulation. Note that the simulation might therefore not start _precisely_ at the specified `INITIAL_TIME`.

    - `"IC_FROM_DAT_FILE_WITH_PSI"` and `"IC_FROM_DAT_FILE_WITHOUT_PSI"`: One can load initial conditions from a `.dat` file. This feature was implemented during development of the code when initial conditions were provided from a collaborator in a `.dat` file. There is an option to also import $$\psi$$ and $$\dot{\psi}$$ from this `.dat` file, or to compute them from the provided $$\phi$$ and $$\dot{\phi}$$ via `mk_initial_psi` in `setup.c`. If you want to read initial conditions from a separately generated file, you might change the method `read_initial_data` in `io.c` accordingly. Also make sure to read and understand the memory layout of the field variables TODO(link, where?).
    - `"IC_FROM_BUNCH_DAVIES"`: The initial conditions are given by the Bunch Davies vacuum. Our implementation is similar in nature to the one used in DEFROST TODO(link, ref). We compute initial values for $$\phi$$ and $$\dot{\phi}$$ following the Bunch Davies spectrum. For detailed information see TODO(link paper, coderef). Subsequently we compute $$\psi$$ and $$\dot{psi}$$ via `mk_initial_psi` in `setup.c`. See also `BUNCH_DAVIES_CUTOFF` in this section.

* `BUNCH_DAVIES_CUTOFF`: (integer, >=0) The Bunch Davies vacumm has an ultra violet divergence that we circumvent by an exponential cutoff. If `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"`, this parameter sets the value of the cutoff as the number of gridpoints. It has to be smaller than half the number of gridpoints (i.e. we cut off the spectrum below the Nyquist frequency). For example, if we have a grid with $$64$$ gridpoints in each direction, `BUNCH_DAVIES_CUTOFF` has to be smaller or equal than $$32$$. If `BUNCH_DAVIES_CUTOFF="0"`, the cutoff is adjusted to the current grid automatically and chosen as large as possible. This parameter also controls the smoothness of the initial conditions for the Bunch Davies vacuum, i.e. a small value corresponds to few modes, i.e. very smooth and well behaved initial conditions.

* `INITIAL_DATAPATH`: See section [file IO](#file-io).

* `INITIAL_TIME`: See section [Time evolution](#time-evolution).

## Potential parameters

* `MASS`: (double, >0) The mass parameter used in the code. Note that this value can be easily rescaled, hence has no physical relevance.

* `INFLATON_MASS`: (double, >0) If `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"`, this value is the inflaton mass in units of the Planck mass. While the actual mass parameter used in the code (for example when computing the potential) can be rescaled to arbitrary values, `INFLATON_MASS` is a physical quantity and only enters once when setting the amplitude of initial fluctuations in the Bunch Davies vacuum. Note that in the code $$8 \pi G = 1$$, hence we set the reduced Planck mass to one. `INFLATON_MASS` is the only mass related parameter in the code that carries physical meaning.

* `MASS_KARSTEN`: (double, >0) Deprecated! Intended for developmental use only. This parameter was usde to scale certain amplitudes for comparison with Karsten Jedamzik's code.

## File IO

* `DATAPATH`: (string) The full path to the output file. All output of a single run is bundled into one `.h5` file. The `DATAPATH` can be relative to the working directory or absolute. It also needs to include the filename and `.h5` extension. The file does not have to exist. It will be created. If it does already exist, it will be __overwritten__. __Example__: `DATAPATH="../data/testrun.h5"`.

* `INITIAL_DATAPATH`: (string) The full path to the file that contains initial conditions for the simulation.
    - If `INITIAL_CONDITIONS="IC_FROM_DAT_FILE_WITH_PSI"` or `INITIAL_CONDITIONS="IC_FROM_DAT_FILE_WITHOUT_PSI"`, this is the full path (including filename and `.dat` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.
    - If `INITIAL_CONDITIONS="IC_FROM_H5_FILE"`, this is the full path (including filename and `.h5` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.

* `VERSION_CONTROL`: If git or mercurial is used as a version control system, one can add a hash of the current revision to the output file. This helps identify the data later. This option was most useful during development. __Warning__: Since the program has to create a pipe, invoke a shell and call git or mercurial as well as receive the output, this is not guaranteed to work on different systems. It was only tested on MacOS 10.11. Choose the `VERSION_CONTROL_NONE` option when in doubt. Valid options are
  - `"VERSION_CONTROL_HG"`: Try to get the hash of the current mercurial revision and include it in the output.
  - `"VERSION_CONTROL_GIT"`: Try to get the hash of the current git revision and include it in the output.
  - `"VERSION_CONTROL_NONE"`: Do not try to get the hash of the current revision. The corresponding dataset will be missing in the output file. __Note__: This is the recommended choice, because it is quite possible that the program does fail at runtime on different systems.

* `WRITE_OUT_BUFFER_NUMBER`: (integer, >0) InstaHeat uses an internal buffer to accumulate output before writing to disk. Since rare and large writeouts are usually faster than frequent small ones, this gives a speedup despite the memory overhead and copy instructions to the buffers. `WRITE_OUT_BUFFER_NUMBER` specifies how many timeslices of the output should be buffered before disk access. If one outputs full 3D fields ([see output](#output)), memory constraints might not allow for a large buffer or slow down the simulation. When in doubt we recommend a value of 1, i.e. keeping only the current time slice in the buffer.

* `POWER_SPECTRUM_BINS`: (integer, >0) The number of bins used to compute the power spectra. This is relevant for the gravitational wave spectrum as well as for any of the outputs (see [output](#output)) ending in `POWER_SPECTRUM`. In the computation of the power spectrum, modes are binned into this number of bins. See TODO(link thesis) for details.

* `TIME_STEP_SKIPS`: (integer, >0) The number of time steps skipped between outputs. To avoid large output files, one can skip `TIME_STEP_SKIPS` many time slices, before writing a time slice to disk again.

* `GSL_OUTPUT_NUMBER`: (integer, >0) If `INTEGRATION_METHOD` is set to one of `"RKF45"`, `"RKCK"` or `"DOPRI89"`, which are all implemented via GSL, we do not have direct access to every timeslice. Instead, we discretize the integration interval from `INITIAL_TIME` to `FINAL_TIME` by `GSL_OUTPUT_NUMBER` many equally spaced points. On each of those points in time we write out the desired data. __Note__: If `TIME_STEP_SKIPS` is set to a value greater than 1 the output will contain less than `GSL_OUTPUT_NUMBER` many timeslices.

* `STRIDE_X`: (integer, >0) The stride in the x-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than the one used for the computation internally.

* `STRIDE_Y`: (integer, >0) The stride in the y-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than the one used for the computation internally.

* `STRIDE_Z`: (integer, >0) The stride in the z-direction of the output of fields. To avoid large ouput files, one can output the fields on a smaller grid than the one used for the computation internally.

__Remarks__:

* Given the number of gridpoints $$n$$ in the x-direction by `GRIDPOINTS_X` and the stride $$s$$ in the x-direction by `STRIDE_X`, the number of gridpoints in the output is computed by $$(n+s-1)/s$$ rounded down to the next integer. The output gridpoints are evenly spaced within the number of gridpoints used for internal computations.

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

then the computation would be done on a $$256 \times 200 \times 59$$ grid. However the output of the fields would be on a $$64 \times 200 \times 10$$ grid.

## Performance parameters

* `THREAD_NUMBER`: (integer, >=0) The number of threads used in OpenMP and also for the discrete Fourier transforms performed by FFTW3. __Important__: If `THREAD_NUMBER="0"`, the number of threads is determined automatically by `omp_get_max_threads()`. While this can be helpful when the target architecture is unknown, we strongly recommend running a short analysis to determine the optimal number of threads. It might differ significantly from the number of processors/cores on the target machine.

* `FFTW_DEFAULT_FLAG`: The planning flag used in the discrete Fourier transforms performed by FFTW3. Take a look at the FFTW3 documentation TODO(link) for more information. The valid options are:
    - `"FFTW_ESTIMATE"`
    - `"FFTW_MEASURE"`
    - `"FFTW_PATIENT"`
    - `"FFTW_EXHAUSTIVE"`
    For short runs we recommend `"FFTW_ESTIMATE"` and for long runs we recommend `"FFTW_PATIENT"`.

## Special paramaters for adaptive time step integration routines

The parameters in this section are only relevant if one of the adaptive time step integration routine is used, i.e. if `INTEGRATION_METHOD` is set to `"RKF45"`, `"RKCK"`, `"DOPRI89"`, or `"DOPRI853"`.

* `RELATIVE_TOLERANCE`: (double, >=0) The relative tolerance in the integration routine. See TODO(link numerical recipes, paper?) for more information.

* `ABSOLUTE_TOLERANCE`: (double, >=0) The absolute tolerance in the integration routine. See TODO(link numerical recipes, paper?) for more information.

* `SMALLEST_SCALING`: (double, >0, <=1) Only relevant for `DOPRI853`. The smallest possible rescaling of the step size in any step.

* `LARGEST_SCALING`: (double, >=1) Only relevant for `DOPRI853`. The largest possible rescaling of the step size in any step.

* `BETA`: (double) Only relevant for `DOPRI853`. This is an internal parameter of the adaptive step size control for PI control of the step size following the rescaling method of Lund. Only change this if you know exactly what you are doing. See TODO(link numerical recipes, Gustaffson paper) for more information.

* `SAFE`: (double, >0, <=1) Only relevant for `DOPRI853`. This is an internal safety parameter of the adaptive step size control, which makes the next step size more likely to be accepted. Only change this if you know exactly what you are doing. See TODO(link numerical recipes, hairer book) for more information. The default value is 0.9.

## Program flow

* `INTEGRATION_METHOD`: There are three integration routines available:
    - `"RK4"`: The standard 4th order Runge Kutta stepper with fixed time step size.
    - `"RKF45"`: A 4th order Runge Kutta Felberg method with 5th order error estimation for adaptive time stepping. We use the GSL implementation TODO(link to gsl).
    - `"RKCK"`: A 4th order Runge Kutta Cash-Karp method with 5th order error estimation for adaptive time stepping. We use the GSL implementation TODO(link to gsl).
    - `"DOPRI89"`: A 8th order Runge Kutta Dormand Prince method with 9th order error estimation for adaptive time stepping. We use the GSL implementation TODO(link to gsl).
    - `"DOPRI853"`: A 8th order Runge Kutta Dormand Prince method with 5th and 3rd order error estimation for adaptive time stepping. A detailed description can be found in TODO(link numerical recipes).

* `ENABLE_FFT_FILTER`: (boolean) Switch for a spectral filter for the fields. If switched on, at each time step the highest modes of $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$ are cut off to avoid aliasing. A more detailed description can be found in TODO(link to thesis).

* `ENABLE_GW`: (boolean) Switch for the extraction of gravitational waves. If switched on, at each time step the generated gravitational wave spectrum is computed and added to the output.

* `ENABLE_FOLLOWUP`: (boolean) Switch for the output of all fields on the very last timeslice. If switched on, regardless of the output settings on the very last timeslice all the fields $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$ and $$a$$ are writtten to the output file on the last time slice of the simulation. We can then use these field values as initial conditions for a subsequent run using the `INITIAL_CONDITIONS="IC_FROM_H5_FILE"` option.

## Miscellaneous

* `SEED`: (integer, >0) The seed for the random number generation used for example to create the Bunch Davies vacuum if `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"` or also for `INITIAL_CONDITIONS="IC_FROM_INTERNAL_FUNCTION"` (depending on what the functions do). We use the GSL implementation of the Mersenne Twister. TODO(link)

* `ENABLE_STIFFNESSCHECK`: (boolean) Only relevant for `DOPRI853`. Use a simple method to check whether the evolution of the equations becomes a stiff problem and abort the evolution if so. This was used mostly for debugging. According to this simple test we have never encountered stiffness.

* `ENABLE_TIMING`: (boolean) One can monitor the wall clock execution time of various parts of the program. If `ENABLE_TIMING` is switched on, the output will contain various datasets with the wall clock runtime of certain parts of the program. The overall runtime is always monitored and unaffected by this parameter. Since the system call to obtain wall clock times can be quite time consuming, we recommend only turning this on for debugging or optimization purposes.

## Output

Most of the simulation parameters are always present in the output (and cannot be switched off). Each parameter takes either the value `"0"` (false) or `"1"` (true). The optional output consists of

* `OUTPUT_PHI`: The scalar field $$\phi$$ on the entire output grid on each output timeslice.

* `OUTPUT_DPHI`: The temporal derivative of the scalar field $$\dot{\phi}$$ on the entire output grid on each output timeslice.

* `OUTPUT_PSI`: The scalar metric perturbation $$\psi$$ on the entire output grid.

* `OUTPUT_DPSI`: The scalar metric perturbation $$\dot{\psi}$$ on the entire output grid on each output timeslice.

* `OUTPUT_RHO`: The energy density $$\rho$$ on the entire output grid on each output timeslice.

* `OUTPUT_PHI_SMRY`: The mean value, the variance, the minimum and the maximum of the scalar field $$\phi$$ on each output timeslice (in this order).

* `OUTPUT_DPHI_SMRY`: The mean value, the variance, the minimum and the maximum of the temporal derivative of the scalar field $$\dot{\phi}$$ on each output timeslice (in this order).

* `OUTPUT_PSI_SMRY`: The mean value, the variance, the minimum and the maximum of the scalar metric perturbation $$\psi$$ on each output timeslice (in this order).

* `OUTPUT_DPSI_SMRY`: The mean value, the variance, the minimum and the maximum of the temporal derivative of the scalar metric perturbation $$\dot{\psi}$$ on each output timeslice (in this order).

* `OUTPUT_RHO_SMRY`: The mean value, the variance, the minimum and the maximum of the energy density $$\rho$$ on each output timeslice (in this order).

* `OUTPUT_PRESSURE_SMRY`: The mean value, the variance, the minimum and the maximum of the pressure $$p$$ on each output timeslice (in this order).

* `OUTPUT_PHI_PS`: The power spectrum of the field $$\phi$$ on each output timeslice.

* `OUTPUT_PSI_PS`: The power spectrum of the field $$\psi$$ on each output timeslice.

* `OUTPUT_RHO_PS`: The power spectrum of the field $$\rho$$ on each output timeslice.

* `OUTPUT_CONSTRAINTS`: The Hamiltonian constraint in $$l_2$$ and $$l_{\infty}$$ norm as well as the momentum constraint in $$l_2$$ and $$l_{\infty}$$ norm (in this order), where the various terms of the constraints have been combined in such a way to give 0. This means that the closer those 4 values are to 0, the better the constraints are fulfilled.

__Remarks__:

* `PHI`, `DPHI`, `PSI`, `DPSI`, `RHO` result in `xout*yout*zout` double values where {x,y,z}out are determined from `GRIDPOINTS_{X,Y,Z}` and `STRIDE_{X,Y,Z}` (see [file IO](#file-io)) on __each__ output timeslice (per field).

* All values ending in `SUMMARY` result in four double values on each output timeslice (per field).

* All values ending in `POWER_SPECTRUM` result in `POWER_SPECTRUM_BINS` double values on each output timeslice (per field).

* The number of output timeslices can be adjusted by `TIME_STEP_SKIPS` and by `GSL_OUTPUT_NUMBER` if `INTEGRATION_METHOD` is set to one of `"RKF45"`, `"RKCK"` or `"DOPRI89"`.