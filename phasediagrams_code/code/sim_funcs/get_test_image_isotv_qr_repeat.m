function [X,v,info] = get_test_image_isotv_qr_repeat(id,D,opts)
%GET_TEST_IMAGE_ISOTV_QR_REPEAT Find s-sparse element in matrix range.
% [X,v,info] = get_test_image_isotv_qr_repeat(id,D,opts) 
% finds an element of the range of D with card(supp(x)) = s.
%
% Inputs: 
%  - id : test image type, currently only 'altprojisotv' supported
%  - D :  Matrix
%  - opts : Struct of options with fields:
%    opts.k : Sparsity
%    opts.randnstate : State integer for randn
%    opts.epsilon :  Stopping criterion (e.g. 10^-12)
%    opts.maxrepetitions : Maximal number of times to repeat procedure
%    opts.numtoskip : Number of state integers to increase when repeat.
%    opts.mask : logical array holding mask of image.
%
% Outputs: 
% - X : Image for which x=X(mask) satisfies v=D*x.
% - v : s-sparse element in range of D.
% - info : Struct with information about the run of the algorithm.
%
% Christian Kruschel and Jakob S. Joergensen (jakj@dtu.dk), 2014.

if ~isfield(opts,'randnstate')
    opts.randnstate = 0;
end
if ~isfield(opts,'epsilon')
    opts.epsilon = 1e-10;
end
if ~isfield(opts,'maxrepetitions')
    opts.maxrepetitions = 1000;
end
if ~isfield(opts,'numtoskip')
    opts.numtoskip = 100000;
end

switch id
    case 'altprojisotv'
        
        n = size(D,2);
        n2 = size(D,1)/2;
        
        [Q,R] = qr(D);
        r = size(D,2)-1;
        Q1 = Q(:,1:r);
        R11 = R(1:r,1:r);
        
        info.status = 0;
        rep_counter = 0;
        info.repetitions = 0;
        
        while (info.status==0) && rep_counter < opts.maxrepetitions
            
            randn('state', opts.randnstate + rep_counter*opts.numtoskip);
            xstart = randn(n,1);
            v = D*xstart;
            t = sqrt(v(1:n2).^2 + v(n2+1:end).^2);
            rep_counter = rep_counter + 1;
            fprintf('Rep: %d\n',rep_counter)
            
            j = 0;
            b = true;
            l = 0;
        
            while l ~= opts.k && b
                % Increment iteration counter.
                j = j + 1;
                
                % Save old v for comparing with new.
                v2 = v;
                
                % Keep only opts.k elements with largest magnitudes t of v.
                [~,I] = sort(t);
                v([I(1:end-opts.k),n2+I(1:end-opts.k)]) = 0;
                
                % Project onto range of D.
                v = D*[R11\(Q1'*v);0];
                %v = P*v;
                
                % Compute the magnitudes of v.
                t = sqrt(v(1:n2).^2 + v(n2+1:end).^2);
                
                % Compute sparsity.
                l = sum(t>opts.epsilon);
                
                % Update stopping criterion.
                b = j < opts.maxiter && norm(v-v2) > opts.epsilon;
                
                % Print progress information to screen.
                fprintf('iter %d of %d. Sparsity %d of target %d.\n',...
                    j,opts.maxiter, l, opts.k);
            end
            
            X = zeros(size(opts.mask));
            %X(opts.mask) = Dp*v;
            
            % Get the x based on the QR. This x is not minimum norm.
            xqr = [R11\(Q1'*v);0];
            
            % Update to minimum-norm xqrmod = xqr + t*e by minimizing
            % over scalar t: ||xqr + t*e||_2^2, which can be solved
            % analystically
            e = ones(size(D,2),1);
            t = -e'*xqr / (e'*e);
            xqrmod = xqr + t*e;
            X(opts.mask) = xqrmod;
            
            if l == opts.k
                info.status = 1;
            end
            info.repetitions = rep_counter;
            
        end
        
    otherwise
        error('Unknown image_id given')
end