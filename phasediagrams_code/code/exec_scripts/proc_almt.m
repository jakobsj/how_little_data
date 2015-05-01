%
% Script to be run for processing the raw ALMT reconstruction data into
% reconstruction errors wrt the original image for saving into single
% mat-file for later display as ALMT phase diagram.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.
%

clear
clc

%% Includes
addpath ../ext/parobj/
addpath ../load_funcs/
addpath ../sim_funcs/

%% Set fixed parameters
N = 64;
spfrac_list = 0.025:0.025:1-0.025;
seed_list = 0:99;
na_list = 1:26;

% Put together as cell array to simplify similar calls
fp = {N,spfrac_list,seed_list,na_list};

%%
image_id = 'signedspikes';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'BP';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'signedspikes';
matrix_type = 'gaussian';
problem_type = 'BP';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'BP_NONNEG';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'gaussian';
problem_type = 'BP_NONNEG';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'TV_EQ_2D_NEUM';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'gaussian';
problem_type = 'TV_EQ_2D_NEUM';
proc_almt_fun(fp{:}, image_id, matrix_type, problem_type)
