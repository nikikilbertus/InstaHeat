dsetsAll = {'a', 'bunch_davies_cutoff', 'calls_rhs', 'commit_hash', 'dimension', ...
'dphi_summary', 'dpsi_summay', 'filter', 'followup', 'gravitational_wave_spectrum', ...
'gravitational_waves', 'gridpoints_internal', 'gridpoints_output', 'inflaton_mass', ...
'mass', 'max_dt_hubble_fraction', 'method', 'phi_power_spectrum', 'phi_summary', ...
'pressure_summary', 'psi_power_spectrum', 'psi_summary', 'rho_power_spectrum', ...
'rho_summary', 'runtime_constraints', 'runtime_copy_buffer', 'runtime_elliptic', ...
'runtime_fftw', 'runtime_fftwplan', 'runtime_filter', 'runtime_stepper', 'runtime_stt', ...
'runtime_summaries', 'runtime_total', 'runtime_writeout', 'seed', 'spatial_bounds_x', ...
'spatial_bounds_y', 'spatial_bounds_z', 'steps_bad', 'steps_ok', 'steps_total', ...
'strides_space', 'strides_time', 'time', 'tolerances', 'phi', 'psi', 'dphi', 'dpsi', ...
'rho', 'h1_summary', 'h2_summary', 'constraints'};

dsetsNames = {'a', 'cutoff', 'calls_rhs', 'commit_hash', 'dim', ...
'dphiS', 'dpsiS', 'filter', 'followup', 'gwps', ...
'gwQ', 'Nint', 'N', 'inflaton_mass', ...
'mass', 'max_dt_hubble_fraction', 'method', 'phips', 'phiS', ...
'pressureS', 'psips', 'psiS', 'rhops', ...
'rhoS', 'runtime_constraints', 'runtime_copy_buffer', 'runtime_elliptic', ...
'runtime_fftw', 'runtime_fftwplan', 'runtime_filter', 'runtime_stepper', 'runtime_stt', ...
'runtime_summaries', 'runtime_total', 'runtime_writeout', 'seed', 'spatial_bounds_x', ...
'spatial_bounds_y', 'spatial_bounds_z', 'steps_bad', 'steps_ok', 'steps_total', ...
'strides_space', 'strides_time', 't', 'tols', 'phi', 'psi', 'dphi', 'dpsi', ...
'rho', 'h1S', 'h2S', 'constraints'};

dsetsRuntime = {'runtime_constraints', 'runtime_copy_buffer', 'runtime_elliptic', ...
'runtime_fftw', 'runtime_fftwplan', 'runtime_filter', 'runtime_stepper', 'runtime_stt', ...
'runtime_summaries', 'runtime_total', 'runtime_writeout'};

dsetsSummary = {'phi_summary', 'dphi_summary', 'psi_summary', 'dpsi_summary', ...
'rho_summary', 'pressure_summary', 'h1_summary', 'h2_summary'};

dsetsSteps = {'steps_total', 'steps_ok', 'steps_bad'};

dsetsPars = {' mass', 'inflaton_mass', 'dimension', 'seed', 'strides_time', ...
'strides_space', 'tolerances', 'gridpoints_internal', 'gridpoints_output', ...
'spatial_bounds_x', 'spatial_bounds_y', 'spatial_bounds_z', 'max_dt_hubble_fraction', ...
'filter', 'followup'};

dsetsSpectra = {'phi_power_spectrum', 'psi_power_spectrum', 'rho_power_spectrum'};

dictNames = containers.Map(dsetsAll, dsetsNames);
dictNamesInv = containers.Map(dsetsNames, dsetsAll);