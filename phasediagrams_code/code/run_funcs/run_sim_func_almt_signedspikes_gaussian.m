function [] = run_sim_func_almt_signedspikes_gaussian(N,spfrac,randseed,numangles)
%RUN_SIM_FUNC_ALMT_SIGNEDSPIKES_GAUSSIAN 
% Run single reconstruction. Wrapper of core_sim_func_matrixtype to
% simplify parameter sweep to create ALMT phase diagram reconstructions.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

%% Include necessary packages

addpath ../ext/mosek_wrap/
addpath ../ext/AIRtools_1.0/
addpath ../sim_funcs/

%% Set up the image class and problem_type

image_id = 'signedspikes';

problem_type = 'BP';

%% Set solver

solver = 'mosek_wrap';

solver_pars.precision = 'medium';
solver_pars.verbose   = 3;

lowlev_pars = struct;
lowlev_pars.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';

%% Set up the sampling matrix.

matrixtype = 'gaussian';
angle_interval = nan; 


%% Directory to save results to

resultspath = fullfile(...
    '../../data_raw/almt/gaussian',...
    image_id,...
    solver);
    
%% Run the core simulation which saves results to the directory resultspath

core_sim_func_matrixtype(N,spfrac,randseed,numangles,image_id,...
    solver,solver_pars,lowlev_pars,resultspath, problem_type, ...
    matrixtype, angle_interval)

