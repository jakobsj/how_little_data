function [] = core_sim_func_k(N,randseed,numangles,k,image_id,...
    solver,solver_pars,lowlev_pars,resultspath, problem_type, ...
    matrixtype, angle_interval)
%CORE_SIM_FUNC_K Main simulation function for DT phase diagrams.
% core_sim_func_k(N,randseed,numangles,k,image_id,...
%    solver,solver_pars,lowlev_pars,resultspath, problem_type, ...
%    matrixtype, angle_interval)
%
% Inputs:
% - N : Problem size (integer).
% - randseed : Random seed of particular image instance.
% - numangles : The number of tomographic projections.
% - k : Target sparsity of test image.
% - image_id: String with supported image type.
% - solver: String with supported name of solver.
% - solver_pars: High-level paramters to the solver.
% - lowlev_pars: Low-level parameters to the solver.
% - resultspath: Path to directory for writing outputs.
% - problem_type: String with type of problem.
% - matrixtype: String with type of matrix.
% - angle_interval: 2-vector with angular range, degrees, default [0,360].
% 
% Outputs: None. Results are saved to disk.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

if nargin < 11
    matrixtype = 'gaussian';
end
if nargin < 12
    angle_interval = [0,360];
end

% Geometry parameters for fanbeamtomo
p = 2*N; % num rays per view
R = 2;
R2 = R*N;
v1 = [-N/2;-R2];
v2 = [N/2;-R2];
w = acosd(( v1'*v2)/(v1'*v1) );
isDisp = 0;

%% Mask for test image
[prob.mask, N_pix] = build_mask(N);

%%

is_save_nan = k > N_pix;

%%

savefilename = fullfile(resultspath,...
    sprintf('N_%d',N),...
    sprintf('res_%s_N_%d_seed_%d_numangles_%d_k_%d',...
    image_id,N,randseed,numangles,k))

% Create directory if not existing
savefiledir = fileparts(savefilename);
if ~exist(savefiledir,'dir')
    mkdir(savefiledir)
end

%%
if ~is_save_nan
    
    switch image_id
        case 'spikes'
            image_opts.randstate = randseed;
            image_opts.k = k; 
            image_opts.mask = prob.mask;
        case 'signedspikes'
            image_opts.randstate = randseed;
            image_opts.k = k; 
            image_opts.mask = prob.mask;
        case 'fftpower_2_0'
            image_opts.pp = 2.0;
            image_opts.randnstate = randseed;
            image_opts.k = k; %round(spfrac*N_pix);
            image_opts.mask = prob.mask;
        case 'altprojisotv'
            opts.randnstate = randseed;
            opts.mask = prob.mask;
            opts.k = k
            opts.maxiter = 1e6;
            opts.epsilon = 1e-10;
            opts.numtoskip = 1e5;
            opts.maxrepetitions = 1e4;
            D = get_D(prob.mask,'neumann');
        otherwise
            error('unknown image_id given')
    end
    
    switch image_id
        case 'altprojisotv'
            if exist([savefilename,'_xorig.mat'],'file')
                load([savefilename,'_xorig.mat'])
                is_orig_loaded = true;
            else
                [X_orig,v_orig,testim_info] = ...
                    get_test_image_isotv_qr_repeat(image_id,D,opts);
                x_orig = X_orig(:);
                is_orig_loaded = false;
            end
        otherwise
            X_orig = get_test_image(image_id,N,image_opts);
            x_orig = X_orig(:);
    end
    
    
    
    %% Run recs
    
    
    switch matrixtype
        case 'gaussian'
            randn('state',randseed*10)
            prob.A = randn(2*N*numangles,N^2);
        case {'fanbeam_equi','fanbeam_equi_offset'}
            theta = linspace(angle_interval(1),angle_interval(2),numangles+1);
            theta = theta(1:end-1);
            prob.A = fanbeamtomo(N,theta,p,R,w,isDisp);
        case 'fanbeam_rand'
            rand('state',randseed*10)
            theta = rand(1,numangles);
            theta = (angle_interval(2)-angle_interval(1))*theta + ...
                angle_interval(1);
            prob.A = fanbeamtomo(N,theta,p,R,w,isDisp);
        case 'random_rays'
            num_rays = 2*N*numangles;
            rand('state',randseed*10);
            theta = 180*rand(num_rays,1);
            x0 = N*rand(num_rays,1) - N/2;
            prob.A = paralleltomo_randomrays(N,theta,x0,isDisp);
    end

    
    [sA1,sA2] = size(prob.A)    
    
    % The perfect data from the image
    prob.b = prob.A*x_orig;
else
    x_orig = NaN;
end

switch image_id
    case 'altprojisotv'
        if ~is_orig_loaded
            save([savefilename,'_xorig.mat'],'x_orig','v_orig','testim_info')
        end
    otherwise
        save([savefilename,'_xorig.mat'],'x_orig')
end

switch problem_type
    case 'BP'
        prob.name = 'BP';
        if is_save_nan
            x_BP = NaN;
            info_BP = NaN;
        else
            fprintf('###################  %s  ################\n',prob.name)
            [x_BP, info_BP] = feval(solver,prob,solver_pars,lowlev_pars);
            if strcmp(info_BP.status, 'Failed')  && (sA1 > N_pix)
                x_BP(prob.mask(:)) = prob.A(:,prob.mask(:)) \ prob.b;
                info_BP.optval = norm(x_BP,1);
                info_BP.status = 'Solved';
            end
        end
        save([savefilename,'_',prob.name,'.mat'],'x_BP','info_BP')
        
    case 'BP_NONNEG'
        prob.name = 'BP_NONNEG';
        if is_save_nan
            x_BP_NONNEG = NaN;
            info_BP_NONNEG = NaN;
        else
            fprintf('###################  %s  ################\n',prob.name)
            
            [x_BP_NONNEG, info_BP_NONNEG] = feval(solver,prob,solver_pars,lowlev_pars);
            if strcmp(info_BP_NONNEG.status, 'Failed')  && (sA1 > N_pix)
                x_BP_NONNEG(prob.mask(:)) = prob.A(:,prob.mask(:)) \ prob.b;
                info_BP_NONNEG.optval = norm(x_BP_NONNEG,1);
                info_BP_NONNEG.status = 'Solved';
            end
        end
        save([savefilename,'_',prob.name,'.mat'],'x_BP_NONNEG','info_BP_NONNEG')
        
    case 'TV_EQ_2D_NEUM'
        prob.name = 'TV_EQ_2D_NEUM';
        
        % Loop over increasingly low solver precision
        precision_list = {'veryverylow',...
            'verylow',...
            'low',...
            'medium',...
            'high',...
            'veryhigh',...
            }
        precision_number = find(strcmp(solver_pars.precision,precision_list))+1;
        
        % Normalize problem
        rownorms = sqrt(sum(prob.A.^2,2));
        S = spdiags(full(1./rownorms.^2),0,sA1,sA1);
        prob.A = S*prob.A;
        prob.b = S*prob.b;
        
        if is_save_nan
            x_ITV = NaN;
            info_ITV = NaN;
        else
            info_ITV.status = 'No';
            fprintf('###################  %s  ################\n',prob.name)
            while precision_number > 1 && ~strcmp(info_ITV.status,'Solved')
                precision_number = precision_number - 1;
                solver_pars.precision = precision_list{precision_number}
                [x_ITV, info_ITV] = feval(solver,prob,solver_pars,lowlev_pars);
            end
            
            info_ITV.precision = solver_pars.precision;
            if strcmp(info_ITV.status, 'Failed')  && (sA1 > N_pix)
                x_ITV(prob.mask(:)) = prob.A(:,prob.mask(:)) \ prob.b;
                info_ITV.optval = sum(   sqrt(sum(reshape(D*x_ITV(prob.mask(:)),size(D,1)/2,2).^2,2)) );
                info_ITV.status = 'Solved';
            end
        end
        save([savefilename,'_',prob.name,'.mat'],'x_ITV','info_ITV')
end