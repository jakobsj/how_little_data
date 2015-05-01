function [x,f,R] = mskITVslv(A,b,varargin)
%  The call  [x,f,R] = mskITVslv(A,b,D)
%  solves the problem (eq constrained isotropic TV)
%
%     min_x  sum_k || D_k*x ||_2
%      s.t.   A*x = b
%
%  where D_k is a 2 by length(x) matrix extracted from D as
%
%  D_k = D( [k, k+length(x)], : ).
%
%  Fourth argument can be a struct with 
%  parameters. For now, options are
%
%   p.echo = level of printing
%          :  0   (no printing at all) [default]
%          :  1,2 (a little printing)
%          :  3   (full print, lots of info)
%   p.tol  = Tolerance on sol
%          : verylow
%          : low
%          : medium [default]
%          : high
%          : veryhigh
%
%  For full flexibility, the most general 
%  (and low level) MOSEK routine is used.
%  it is called mosekopt. 
%  see MOSEK documentation for more info.
%
% Anders Skajaa, Jakob S. Joergensen (jakj@dtu.dk), 2014.

if nargin > 2
    D = varargin{1};
else
    error('A thrid argument, D, must be given.');
end
if nargin > 3
    p = varargin{2};
else
    p = struct;
end
if nargin > 4
    lowlev_pars = varargin{3};
else
    lowlev_pars = struct;
end

if ~isfield(p,'echo')
    p.echo = 0;
end
if ~isfield(p,'tol')
    p.tol = 'medium';
end
if ~isfield(p,'type')
    p.type = 'direct';
end

A        = sparse(A);
[m,n]    = size(A);

cmd = 'minimize';
cmd = sprintf('%s echo(%i)',cmd,p.echo);

switch p.type
    case 'direct'
        Dcell = cell(n,1);
        for k = 1:n
            Dcell{k} = D([k,k+n],:);
        end
        
        
        prob.a = [sparse(2*n+size(A,1),n),...
            [cell2mat(Dcell);A], ...
            speye(2*n+size(A,1),2*n)];
        prob.buc = [zeros(2*n,1); b];
        prob.blc = prob.buc;
        for k = 1:n
            prob.cones{k}.type = 'MSK_CT_QUAD';
            prob.cones{k}.sub = [k, 2*n + (k-1)*2+[1,2]];
        end
        prob.c = [ones(n,1); zeros(n+2*n,1)];
        
        % Simple bounds on t, x, y, (z, zhat)
        prob.blx = [];
        prob.bux = [];
    case 'svd'
        Dcell = cell(n,1);
        for k = 1:n
            Dcell{k} = D([k,k+n],:);
        end
        
        [U,s,V] = svd(full(A)); 
        s = diag(s);
        [NA,x_part] = nullandleastnormsol(U,s,V,b);
        [sNA1,sNA2] = size(NA);
        
        prob.a = [sparse(2*n+n,n),...
            [cell2mat(Dcell);speye(n,n)], ...
            speye(2*n+n,2*n),...
            [sparse(2*n,sNA2);-NA]];
        prob.buc = [zeros(2*n,1); x_part];
        prob.blc = prob.buc;
        for k = 1:n
            prob.cones{k}.type = 'MSK_CT_QUAD';
            prob.cones{k}.sub = [k, 2*n + (k-1)*2+[1,2]];
        end
        prob.c = [ones(n,1); zeros(n+2*n+sNA2,1)];

        % Simple bounds on t, x, y, (z, zhat)
        prob.blx = [];
        prob.bux = [];
end

NoptTol = 100;
if ischar(p.tol)
    switch p.tol
        case 'veryverylow'
            Ctol = 1e-2;
        case 'verylow'
            Ctol = 1e-3;
        case 'low'
            Ctol = 1e-5;
        case 'medium'
            Ctol = 1e-8;
        case 'high'
            Ctol = 1e-10;
            NoptTol = 50;
        case 'veryhigh'
            Ctol = 1e-14;
            NoptTol = 10;
    end
elseif isfloat(p.tol)
    if p.tol > 0.9
        fprintf('Warning: p.tol meaninglessly large\n')
        fprintf('         using p.tol = ''medium''\n');
        p.tol = 1e-8;
    end
    Ctol = p.tol;
end

%Controls primal feasibility
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = Ctol;

%Controls dual feasibility
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = Ctol;

%Controls relative gap
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = Ctol;

% Near optimal tol:
param.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL = NoptTol;

%Controls when the problem is declared infeasible
param.MSK_DPAR_INTPNT_TOL_INFEAS = Ctol;

%Controls when the complementarity is reduced enough
param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = Ctol;

% this is important:
if m >= n
    param.MSK_IPAR_PRESOLVE_LINDEP_USE = 'MSK_OFF';
end

% Pass options in lowlev_pars into MOSEK by creating the same fields in
% param.
fnames = fieldnames(lowlev_pars);
for k = 1:length(fnames)
    s = fnames{k};
    param.(s) = lowlev_pars.(s);
end

[r,res] = mosekopt(cmd,prob,param);
try
    xx      = res.sol.itr.xx(n +(1:n));
    f      = res.sol.itr.pobjval;
catch excep
    xx = NaN;
    f = NaN;
end

R.res  = res;
R.ef   = r;
x      = xx;
