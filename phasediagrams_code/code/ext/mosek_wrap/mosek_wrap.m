function [xsol,info] = mosek_wrap(prob, solver_pars, lowlev_pars)
%MOSEK_WRAP  Wrapper to solve selected problems using MOSEK
% [xsol,info] = mosek_wrap(prob, solver_pars, lowlev_pars)
%
% Inputs:
%  - prob: struct describing the optimization problem, fields:
%          name: String with problem name, see switch-cases below.
%          mask: logical array with true indicating the pixels of the image
%          A: The measurement matrix.
%          b: The right-hand side/data vector.
%  - solver_pars: struct holding high-level options for the solver, fields:
%          verbose: how much to print, 3=most, 2,1=some, 0=none.
%          precision: How precise/accurate a solution to solve for:
%                     veryhigh, high, medium (default), low, verylow.
%
% Outputs:
%  - xsol: The solution to the problem min J(x) st A*x=b, possibly with 
%          elementwise x>=0 for J being either the L1-norm (BP), L1-norm 
%          with nonnegativity constraint (BP_NONNEG) or isotropic TV with 
%          appropriate boundary conditions ('TV_EQ_2D_ZERO',
%          'TV_EQ_2D_NONE', 'TV_EQ_2D_NEUM').
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

if nargin < 2
    solver_pars = struct;
end

if nargin < 3
    lowlev_pars = struct;
end

info = init_info();

% Specify solver parameters
if nargin > 1
    if isfield(solver_pars,'verbose') && ~isempty(solver_pars.verbose)
        p.echo = solver_pars.verbose;
    else % The default
        p.echo = 3;
    end
    if isfield(solver_pars,'precision') && ~isempty(solver_pars.precision)
        p.tol = solver_pars.precision;
    else
        p.tol ='medium';
    end
end
    

switch prob.name
    
    % Basis pursuit.
    case 'BP'
        
        % Allow a mask
        is_masked = isfield(prob,'mask');
        if is_masked
            Norig = size(prob.A,2);
            prob = set_up_masking(prob);
        end
        N = size(prob.A,2);
        
        tic
        [xsol,info.optval,R] = mskL1slv(prob.A,prob.b,[],p,lowlev_pars);
        info.time = toc;
        
        if is_masked
            xsol = post_proc_masking(xsol,prob,Norig);
        end
        
        if strcmp(R.res.sol.itr.solsta,'OPTIMAL') || ...
                strcmp(R.res.sol.itr.solsta,'NEAR_OPTIMAL')
            info.status = 'Solved';
        else
            info.status = 'Failed';
        end
    
    % Nonnegativity-constrained basis pursuit.
    case 'BP_NONNEG'
        
        % Allow a mask
        is_masked = isfield(prob,'mask');
        if is_masked
            Norig = size(prob.A,2);
            prob = set_up_masking(prob);
        end
        N = size(prob.A,2);
        
        %p.type = 'luboundelim';
        
        p.nonneg = 'on';
        
        tic
        [xsol,info.optval,R] = mskL1slv(prob.A,prob.b,[],p,lowlev_pars);
        info.time = toc;
        
        if is_masked
            xsol = post_proc_masking(xsol,prob,Norig);
        end
        
        if strcmp(R.res.sol.itr.solsta,'OPTIMAL') || ...
                strcmp(R.res.sol.itr.solsta,'NEAR_OPTIMAL')
            info.status = 'Solved';
        else
            info.status = 'Failed';
        end
    
    % Isotropic TV in 2D with zero/none/neumann BC
    case {'TV_EQ_2D_ZERO','TV_EQ_2D_NONE', 'TV_EQ_2D_NEUM'}
        
        % Allow a mask
        is_masked = isfield(prob,'mask');
        if is_masked
            N = size(prob.A,2);
            prob = set_up_masking(prob);
        end
        Nmasked = size(prob.A,2);
        sqrtN = round(sqrt(N));
        
        switch prob.name(10:13)
            case 'ZERO'
                L = spdiags([-ones(sqrtN,1),ones(sqrtN,1)],...
                    [0,1],sqrtN,sqrtN);
                I = speye(sqrtN,sqrtN);
                D = [kron(I,L);
                     kron(L,I)];
                if is_masked
                    D = D([prob.mask(:);prob.mask(:)],prob.mask(:));
                end
            case 'NONE'
                L = spdiags(ones(sqrtN,1)*[-1,1],[0,1],sqrtN-1,sqrtN);
                I = speye(sqrtN);
                D = [kron(I,L); 
                     kron(L,I)];
                if is_masked
                    % Remove columns for pixels outside mask
                    D = D(:,prob.mask(:));
                    % Remove rows with only one component left
                    D = D(sum(D~=0,2)>1,:);
                end
            case 'NEUM'
                L = spdiags(ones(sqrtN,1)*[-1,1],[0,1],sqrtN-1,sqrtN);
                L(sqrtN,sqrtN) = 0;
                I = speye(sqrtN);
                D = [kron(I,L); 
                     kron(L,I)];
                if is_masked
                    % Remove columns and rows for pixels outside mask
                    D = D([prob.mask(:);prob.mask(:)],prob.mask);
                    % Set to zero the gradient with components outside mask
                    D(sum(D~=0,2)==1,:) = 0;
                end
        end
        
        % Depending on size of problem, use more efficient formulation.
        [sA1,sA2] = size(prob.A);
        if sA1 < 0.4*sA2
            p.type = 'direct';
        else
            p.type = 'svd';
        end
        
        tic
        [xsol,info.optval,R] = mskITVslv(prob.A,prob.b,D,p,lowlev_pars);
        info.time = toc;
        
        if is_masked
            xsol = post_proc_masking(xsol,prob,N);
        end
        
        if strcmp(R.res.sol.itr.solsta,'OPTIMAL') || ...
                strcmp(R.res.sol.itr.solsta,'NEAR_OPTIMAL')
            info.status = 'Solved';
        else
            info.status = 'Failed';
        end

    otherwise
        error('Unsupported optimization problem given')
end


% Aux functions.


function info = init_info()
info = struct;
info.time   = [];
info.status = [];
info.optval = [];
info.alg_info = [];


function prob = set_up_masking(prob)
prob.A = prob.A(:,prob.mask(:));


function xsol = post_proc_masking(xsol,prob,Norig)
xtemp = xsol;
xsol = zeros(Norig,1);
xsol(prob.mask(:)) = xtemp;
