function proc_dt_fun(N, seed_list, na_list, k_div_n_list, image_id,...
    matrix_type, problem_type)
%PROC_DT_FUN   Process raw recons. to rel. errors for DT phase diagram.
% proc_dt_fun(N, seed_list, na_list, k_div_n_list, image_id,...
%     matrix_type, problem_type)
% Read reconstructed images from data files, compute reconstruction error
% wrt. original and save for all images into single mat-file for later
% display as DT phase diagram.
%
% Inputs:
%  - N: Problem size.
%  - seed_list: List of random seeds.
%  - na_list: List of number of angles.
%  - k_div_n_list: List of values of k/N, where k is target sparsity.
%  - image_id: String with the type of image.
%  - matrix_type: String with the type of sampling matrix.
%  - problem_type: String with the problem type.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

solver = 'mosek_wrap';

%% Set path to read raw data from
data_raw_path = fullfile(...
    '../../data_raw/dt', ...
    matrix_type, ...
    image_id, ...
    solver, ...
    sprintf('N_%d',N));

%% Set path to write processed data to
data_proc_path = fullfile(...
    '../../data_processed_user/dt', ...
    matrix_type, ...
    image_id);

% Create if non-existing
if ~exist(data_proc_path,'dir')
    mkdir(data_proc_path)
end

%% Set up parameter object

parvals_k_div_n = {...
    N, ...
    seed_list, ...
    na_list, ...
    k_div_n_list, ...
    }

% First, change pararray to k values due to file naming convention.
po = parobj;
po.setValues( parvals_k_div_n );
po.buildArray() 
new_array = po.array;
new_array(:,4) = po.array(:,4).*po.array(:,3)*2*N;
po.setArray(new_array);

po.setNames( {'N','seed','numangles','k'} );
po.setTypes( {'%d', '%d', '%d', '%d'} );
po.setStub( sprintf('res_%s',image_id) );


% File name post-component
post = ['_',problem_type,'.mat']

% Sweep over all parameters and load data and compute relative 2-norm error
% wrt original.
ell2ell1 = cell2mat(...
    po.loadResults('load_single_dt_result', data_raw_path, post));

% Save
proc_filename = ['all_summary',post];
save(fullfile(data_proc_path,proc_filename),'po','ell2ell1')