function proc_almt_fun(N, spfrac_list, seed_list, na_list, image_id,...
    matrix_type, problem_type)
%PROC_ALMT_FUN   Process raw recons. to rel. errors for ALMT phase diagram.
% proc_almt_fun(N, spfrac_list, seed_list, na_list, image_id,...
%     matrix_type, problem_type)
% Read reconstructed images from data files, compute reconstruction error
% wrt. original and save for all images into single mat-file for later
% display as ALMT phase diagram.
%
% Inputs:
%  - N: Problem size.
%  - spfrac_list: List of relative sparsities.
%  - seed_list: List of random seeds.
%  - na_list: List of number of angles.
%  - image_id: String with the type of image.
%  - matrix_type: String with the type of sampling matrix.
%  - problem_type: String with the problem type.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

solver = 'mosek_wrap';

%% Set path to read raw data from
data_raw_path = fullfile(...
    '../../data_raw/almt', ...
    matrix_type, ...
    image_id, ...
    solver, ...
    sprintf('N_%d',N));

%% Set path to write processed data to
data_proc_path = fullfile(...
    '../../data_processed_user/almt', ...
    matrix_type, ...
    image_id);

% Create if non-existing
if ~exist(data_proc_path,'dir')
    mkdir(data_proc_path)
end

%% Set up parameter object

parvals = {...
    N, ...
    spfrac_list, ...
    seed_list, ...
    na_list, ...
    }

po = parobj;

po.setValues( parvals );
po.setNames( {'N','spfrac','seed','numangles'} );
po.setTypes( {'%d', '%1.3e', '%d', '%d'} );
po.setStub( sprintf('res_%s',image_id) );
po.buildArray()

% File name post-component
post = ['_',problem_type,'.mat']

% Sweep over all parameters and load data and compute relative 2-norm error
% wrt original.
ell2ell1 = cell2mat(po.loadResults('load_single_almt_result', ...
                              data_raw_path, ...
                              post));

% Save
proc_filename = ['all_summary',post];
save(fullfile(data_proc_path,proc_filename),'po','ell2ell1')