%
% Script to be run for processing the raw DT reconstruction data into
% reconstruction errors wrt the original image for saving into single
% mat-file for later display as DT phase diagram.
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
seed_list = 0:99;
na_list = 1:25;
k_div_n_list = (1:32) / 32;

% Put together as cell array to simplify similar calls
fp = {N,seed_list,na_list,k_div_n_list};

%%
image_id = 'signedspikes';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'BP';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'signedspikes';
matrix_type = 'gaussian';
problem_type = 'BP';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'signedspikes';
matrix_type = 'fanbeam_rand';
problem_type = 'BP';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'signedspikes';
matrix_type = 'random_rays';
problem_type = 'BP';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'gaussian';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'fanbeam_rand';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'spikes';
matrix_type = 'random_rays';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'fftpower_2_0';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'fftpower_2_0';
matrix_type = 'gaussian';
problem_type = 'BP_NONNEG';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'fanbeam_equi_offset20';
problem_type = 'TV_EQ_2D_NEUM';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'gaussian';
problem_type = 'TV_EQ_2D_NEUM';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'fanbeam_rand';
problem_type = 'TV_EQ_2D_NEUM';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)

%%
image_id = 'altprojisotv';
matrix_type = 'random_rays';
problem_type = 'TV_EQ_2D_NEUM';
proc_dt_fun(fp{:}, image_id, matrix_type, problem_type)
