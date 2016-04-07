# Parameter file documentation (`parameters.sh`)

This is an explanation of the parameters in the file `parameters.sh`. We want to keep the parameter file itself short and concise to allow for quick changes. Therefore we do not clutter it by comments, but outsource the documentation to this document.

__IMPORTANT:__ Unlike parameters that are passed to the program at runtime, here parameters are read _before_ compilation. In the build process initiated by `make`, first the `parameters.sh` file is read by a shell script and the values are filled directly into the source code. Then the code is compiled. The executable program does not take any arguments anymore. After compilation no parameters can be changed. This allows us to compile only the code needed for a specific choice of parameters in a highly optimized way. Thereby we can optimize program size as well as memory usage and runtime, because we save potentially unnecessary memory allocations and conditionals.

## How to handle this file

* Do not delete any parameters. All of them are needed.
* Do not change the names of any parameters.
* All parameter values are given as strings, i.e. enclosed by `" "`.
* Boolean values are encoded by 1 (true) and 0 (false).
* Numerical values can be given in any format that is known by C99. One can use `"PI"` for the value of pi.
* You can add comments to the file, starting with #.
* The order of the parameters does not matter. We grouped parameters into groups:
    - __common__: These are changed frequently in different simulation. Feel free to play around with them.
    - __uncommon__: These are not changed as often. However, you might want to try some special things or tune some of the parameters. Some of them are useful for debugging.
    - __rare__: There is hardly any reason to change these parameters. Only perform changes if you know exactly what you are doing and why you are doing it.
    - __output__: Here we specify the output of the simulation.

In the documentation, we will group the parameters according to other criteria.

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
* The simulation always runs in an interval, rectangle or cuboid with periodic boundary conditions in each direction.
* The number of spatial gridpoints is also the number of Fourier modes used in the spectral parts of the code. (There is no padding in momentum space, hence one has to choose sufficiently many gridpoints to avoid aliasing. See also the section on filtering for more information.)
* The `SPATIAL_UPPER_BOUND` must be larger than `SPATIAL_LOWER_BOUND` for each direction.
* For the simulation results usually only the box length in each direction (i.e. `SPATIAL_UPPER_BOUND - SPATIAL_LOWER_BOUND`) is relevant. However, the actual spatial gridpoint coordinates are needed, if the initial conditions are constructed from internal functions. See [initial conditions](#initial-conditions) for the construction of initial conditions.

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

* `INITIAL_TIME`: (double) The start time of the simulation. Since the evolution equations do not explicitly depend on time, this can be chosen arbitrarily.
* `FINAL_TIME`: (double, >`INITIAL_TIME`) The end time of the simulation. This is one of the key drivers of the duration of the simulation.
* `DELTA_T`: (double, >0) The initial value for the time step. Its meaning depends on the chosen integration routine [program flow](#program-flow).
    - If the evolution is performed by the fixed time step RK4 routine this is the fixed time step used throughout the simulation (see remarks for one exception).
    - If the evolution is performed by the adaptive time step Dormand Prince 853 routine, the time step is adjusted. Thus `DELTA_T` is just the initial try.
* `MINIMAL_DELTA_T`: (double, >0) This is only relevant if the chosen integration routine is the Dormand Prince 853 stepper (adaptive stepsize) [program flow](#program-flow). It gives a lower bound on the step size. Once the integration routine tries to reduce the step size below this limit, the integration is terminated.
* `MAX_STEPS`: (integer, >0) The maximal number of steps performed by the integration routine.
* `MAX_DT_HUBBLE_FRACTION`: (double, >0) This is only relevant if the chosen integration routine is the Dormand Prince 853 stepper (adaptive stepsize) [program flow](#program-flow). To avoid large timesteps in the beginning of the evolution that could lead to instabilities, we limit the time step from above by `MAX_DT_HUBBLE_FRACTION` times the Hubble time $1/H$. A typical values is on the order of $10^{-3}$ to $10^{-2}$.

__Remarks__:
* Even for the fixed time step RK4 method, the time step might be different for the very last step. To ensure that the simulation always ends exactly at the specified `FINAL_TIME`, the time step might be adjusted for the very last step.

## Initial conditions

* `A_INITIAL`: (double, >0) The initial value for the scale factor a, see TODO (equation links?)
* `INITIAL_CONDITIONS`: This parameter determines how to obtain the initial conditions for the fields $\phi$, $\dot{\phi}$ and in certain cases also $\psi$ and $\dot{\psi}$. The valid options are:
    - `"IC_FROM_INTERNAL_FUNCTION"`: The initial $\phi$ and $\dot{\phi}$ are given as functions `phi_init` and `dphi_init` in `setup.c`. These are evaluated on the spatial grid specified by [simulation volume](#simulation-volume). The initial values for $\psi$ are computed via TODO (how to state this?).
    - `"IC_FROM_H5_FILE"`: One can load initial conditions from a previous simulation. The following conditions have to be satisfied:
          + The `INITIAL_DATAPATH` parameter in section [initial conditions](#initial-conditions) points to the output file of a previous simulation from which a certain timeslice should be used as initial conditions for the current run.
          + The `INITIAL_TIME` parameter in section [initial conditions](#initial-conditions) has to lie in the simulation region of the previous output file, i.e. `INITIAL_TIME` has to be larger than `INITIAL_TIME` of the previous simulation and smaller than `FINAL_TIME` of the previous simulation.
          + The number of gridpoints (`GRIDPOINTS_{X,Y,Z}` parameters in section [simulation volume](#simulation-volume)) have to coincide with the number of gridpoints from the previous simulation.

      Then the program will find the time slice __closest__ to the specified `INITIAL_TIME` in the previous output file and use all necessary values ($a$, $\phi$, $\dot{\phi}$, $\psi$, $\dot{\psi}$) as initial conditions for the new simulation. Note that the simulation might therefore not start exactly at the specified `INITIAL_TIME`.

    - `"IC_FROM_DAT_FILE"`: One can load initial conditions from a `.dat` file. This feature was implemented during development of the code when initial conditions were provided from a collaborator in a `.dat` file. It is obsolete now. If you want to read initial conditions from a separately generated file you might change the method `read_initial_data()` in `filehandling.c` accordingly. Also make sure to read and understand the memory layout of the field variables TODO(link, where?).
    - `"IC_FROM_BUNCH_DAVIES"`: The initial conditions are given by the Bunch Davies vacuum. We followed the implementation of DEFROST TODO(link, ref) to compute initial values for $\phi$ and $\dot{\phi}. For detailed information see TODO(paper, coderef). Subsequently we compute $\psi$ and $\dot{psi}$ from TODO(how to formulate this?).
* `INITIAL_DATAPATH`: See section [file IO parameters](#file-io-parameters).
* `INITIAL_TIME`: See section [simulation volume](#simulation-volume).

## Potential parameters

* `MASS`: (double, >0) The mass parameter used in the code. Note that this value can be easily rescaled, hence has little physical relevance.
* `MASS_PLANCK`: (double, >0) If `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"`, this value is the ratio Planck mass/inflaton mass. This determines the amplitude of the vacuum fluctuations of the Bunch Davies vacuum. Note in the code $8 \pi G = 1$, hence the reduced Planck mass is one. Also the `MASS` is a parameter that does not itself carry physical meaning. This `MASS_PLANCK` is the only value that carries physical meaning.
* `MASS_KARSTEN`: (double, >0) Obsolete. Internal use only.

## File IO

* `DATAPATH`: (string) The full path to the output file. All output of a single run is bundled into one `.h5` file. The `DATAPATH` can be relative to the working directory. It also needs to include the filename and `.h5` extension. The file does not have to exist. It will be created. If it does already exist, it will be __overwritten__. __Example__: `DATAPATH="../data/myrun.h5"`.
* `INITIAL_DATAPATH`: (string) The full path to the file that contains initial conditions for the simulation.
    - If `INITIAL_CONDITIONS="IC_FROM_DAT_FILE"`, this is the full path (including filename and `.dat` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.
    - If `INITIAL_CONDITIONS="IC_FROM_H5_FILE"`, this is the full path (including filename and `.h5` extension) to the file. See [initial conditions](#initial-conditions) for detailed instructions.
* `VERSION_CONTROL`: If git or mercurial is used as a version control system, one can add a hash of the current revision to the output file. This helps identify the data later. This was most useful during development. Since the program has to call shell command and receive their output, this is not guaranteed to work on different systems. Valid options are
  - `"VERSION_CONTROL_HG"`: Try to get the hash of the current mercurial revision and include it in the output.
  - `"VERSION_CONTROL_GIT"`: Try to get the hash of the current git revision and include it in the output.
  - `"VERSION_CONTROL_NONE"`: Do not try to get the hash of the current revision. Nothing is added to the output file.
* `WRITE_OUT_BUFFER_NUMBER`: (integer, >0) TODO
* `POWER_SPECTRUM_BINS`: (integer, >0) TODO
* `TIME_STEP_SKIPS`: (integer, >0) TODO
* `STRIDE_X`: (integer, >0) TODO
* `STRIDE_Y`: (integer, >0) TODO
* `STRIDE_Z`: (integer, >0) TODO

## Performance parameters

* `THREAD_NUMBER`: (integer, >=0) The number of threads used in openMP and also for the discrete Fourier transforms performed by FFTW3. __Important__: If `THREAD_NUMBER="0"`, the number of threads is determined automatically by `omp_get_max_threads()`. This is the default value.
* `FFTW_DEFAULT_FLAG`: The planning flag used in the discrete Fourier transforms performed by FFTW3. Take a look at the FFTW3 documentation TODO (link) for more. The valid options are:
    - `"FFTW_ESTIMATE"`
    - `"FFTW_MEASURE"`
    - `"FFTW_PATIENT"`
    - `"FFTW_EXHAUSTIVE"`

## Special paramaters for the Dormand Prince 853 integrator

* `RELATIVE_TOLERANCE`: TODO
* `ABSOLUTE_TOLERANCE`: TODO
* `SMALLEST_SCALING`: TODO
* `LARGEST_SCALING`: TODO
* `BETA`: TODO
* `SAFE`: TODO

## Program flow

* `INTEGRATION_METHOD`: TODO
* `PSI_METHOD`: TODO
* `ENABLE_FFT_FILTER`: TODO

## Miscellaneous

* `SEED`: The seed for the random number generation used for example to create the Bunch Davies vacuum if `INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"` or also for `INITIAL_CONDITIONS="IC_FROM_INTERNAL_FUNCTION"` (depending on what the functions do).

## Output

TODO
