function [ell2err,ell1rel,info_sol] = load_single_dt_result(savefilename)
%LOAD_SINGLE_DT_RESULT   Load recon. and compute error wrt original image
% [ell2err,ell1rel,info_sol] = load_single_dt_result(savefilename)
%
% Inputs:
%  - savefilename : Full filename (including path) to reconstruction.
%
% Outputs:
%  - ell2err: Relative 2-norm error
%  - ell1rel: Relative difference in 1-norm.
%  - info_sol: 1 if the solution status is 'Solved', otherwise 0.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

% Load the result
load(savefilename)

% Load the original
idx = strfind(savefilename,'BP');
if isempty(idx)
    idx = strfind(savefilename,'TV');
end
savefilename_xorig = [savefilename(1:idx-1),'xorig.mat'];
load(savefilename_xorig)

if exist('x_BP','var')
    x_sol = x_BP;
    info_sol = double(strcmp(info_BP.status,'Solved'));
elseif exist('x_BP_NONNEG','var')
    x_sol = x_BP_NONNEG;
    info_sol = double(strcmp(info_BP_NONNEG.status,'Solved'));
elseif exist('x_ITV','var')
    x_sol = x_ITV;
    info_sol = double(strcmp(info_ITV.status,'Solved'));
end

% Compute relative 2-norm error and 1-norm objective difference.
ell2err = norm(x_sol(:)-x_orig(:),2) / norm(x_orig(:));
ell1rel = abs(norm(x_sol(:),1) - norm(x_orig(:),1)) / norm(x_orig(:),1);
