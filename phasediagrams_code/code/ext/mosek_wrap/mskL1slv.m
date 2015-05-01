function [x,f,R] = mskL1slv(A,b,varargin)
%  The call  [x,f,R] = mskL1slv(A,b)
%  solves the problem
%
%     min_x  || x ||_1
%      s.t.   A*x = b
%
%  Third argument can be a matrix D
%  in that case, [x,f,R] = mskL1slv(A,b,D)
%  solves the problem
%
%     min_x  || D*x ||_1
%      s.t.    A*x = b
%
%  Fourth argument can be a struct with 
%  parameters. For now, options are
%
%   p.echo = level of printing
%          :  0   (no printing at all) [default]
%          :  1,2 (a little printing)
%          :  3   (full print, lots of info)
%
%   p.type = formulation of problem
%          : 'lubound'  [default]
%            upper and lower bounds x by a 
%            new variable. Better sparsity 
%            but more variables than option below
%          : 'abssplit' 
%            splits x in positive/negative parts
%            and solves as LP with simple 
%            lower bounds on each component
%            fewer variables, but worse sparsity
%   
%   p.nonneg = Use nonnegativity constraint
%            : 'off' [default]
%              Do not impose nonnegativity.
%              'on'
%              Restrict x to be nonnegative
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
    D = [];
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

if ~isfield(p,'type')
    p.type = 'lubound';
end
if ~isfield(p,'echo')
    p.echo = 0;
end
if ~isfield(p,'tol')
    p.tol = 'medium';
end
if ~isfield(p,'nonneg')
    p.nonneg = 'off';
end

Disempty = isempty(D);
A        = sparse(A);
[m,n]    = size(A);
In       = speye(n);
en       = ones(n,1);
subtract = [];

if ~Disempty
    [mD,junk] = size(D);
    ImD       = speye(mD);
end

cmd = 'minimize';
cmd = sprintf('%s echo(%i)',cmd,p.echo);

% Nonnegativity only implemented for empty D
if strcmp(p.nonneg, 'on') && Disempty
    prob.a   = A;
    if min(b) >= 0
        prob.buc = b; 
        prob.blc = b; 
    else
        prob.buc = b;
        prob.blc = b;
    end
    
    prob.bux = inf(n,1);
    prob.blx = zeros(n,1);
    prob.c   = en;
else
    switch p.type
        case 'abssplit'
            if Disempty
                prob.a   = [A,-A];
                prob.buc = b;
                prob.blc = b;
                prob.bux = [];
                prob.blx = zeros(2*n,1);
                prob.c   = [en;en];
                subtract = n+1:2*n;
            else
                prob.a   = ...
                    [A,sparse(m,2*mD);
                    D,-ImD,ImD];
                prob.buc = [b;zeros(mD,1)];
                prob.blc = prob.buc;
                prob.bux = [];
                prob.blx = [-inf*ones(n,1);zeros(2*mD,1)];
                prob.c   = [zeros(n,1);ones(2*mD,1)];
                
            end
            
        case 'lubound'
            if Disempty
                prob.a   = ...
                    [A,sparse(m,n);
                    In,-In;
                    -In,-In];
                prob.buc = [b;zeros(2*n,1)];
                prob.blc = [b;-inf*ones(2*n,1)];
                prob.bux = [];
                prob.blx = [];
                prob.c   = [zeros(n,1);en];
            else
                prob.a    = ...
                    [A,sparse(m,2*mD);...
                    D,-ImD,sparse(mD,mD);...
                    sparse(mD,n),ImD,-ImD;...
                    sparse(mD,n),-ImD,-ImD];
                prob.buc = [b;zeros(3*mD,1)];
                prob.blc = [b;zeros(mD,1);-inf*ones(2*mD,1)];
                prob.bux = [];
                prob.blx = [];
                prob.c   = [zeros(n+mD,1);ones(mD,1)];
            end
            
        case 'luboundelim'
            
            if Disempty
                prob.a   = ...
                    [A,sparse(m,n);...
                    -2*speye(n),speye(n)];
                prob.buc = [b;sparse(n,1)];
                prob.blc = [b;-inf*ones(n,1)];
                prob.bux = [];
                prob.blx = [-inf*ones(n,1);zeros(n,1)];
                prob.c   = [en;-en];
            else
                em = ones(mD,1);
                prob.a   = ...
                    [A,sparse(m,mD);...
                    -2*D,speye(mD)];
                prob.buc = [b;sparse(mD,1)];
                prob.blc = [b;-inf*ones(mD,1)];
                prob.bux = [];
                prob.blx = [-inf*ones(n,1);zeros(mD,1)];
                prob.c   = [D'*em;-em];
            end
    end
end


NoptTol = 100;
if ischar(p.tol)
    switch p.tol
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

% If nonnegative, switch off presolve, since some problems are determined
% to be infeasible
if strcmp(p.nonneg, 'on') && Disempty
    param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
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
    xx      = res.sol.itr.xx(1:n);
    if ~isempty(subtract)
        xx  = xx - res.sol.itr.xx(subtract);
    end
    f      = res.sol.itr.pobjval;
catch excep
    xx = NaN;
    f = NaN;
end

R.res  = res;
R.ef   = r;
x      = xx;
