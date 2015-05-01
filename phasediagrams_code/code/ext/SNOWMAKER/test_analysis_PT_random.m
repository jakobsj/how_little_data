%% Generate PT diagram for analysis L1

% Length of signal (ambient dimension)
d = 200;
nTrials = 50;

relTol = 1e-5; % Relative tolerance for success

dictSkip = 50;
measSkip = 5; 
nnzSkip = 5;

% Number of elements in the sparisifying dictionary 
nn = d:dictSkip:round(2*d);

% Number of measurements
mm = 1:measSkip:d;

% Sparsity of Dict*x0
ss = 1:nnzSkip:d; % It's impossible to have more than d zeros in Dict*x0

% Store the results.  NOTE ORDER OF ARGUMENTS
results = zeros(numel(mm),numel(ss),numel(nn));

fprintf('\n');
for ii = 1:numel(nn)
    n = nn(ii); % Size of dictionary
    
    for jj = 1:numel(mm)
        m = mm(jj); % Number of measurments
        
        for kk = 1:numel(ss)
            s = ss(kk); % Sparsity of x0 (number of ZEROS)
            
            % Check for auto-failure conditions
            if (d-s) > m
                continue;
            end
            fprintf('(ii,jj,kk) = (%d,%d,%d) of (%d,%d,%d)\n',...
                ii,jj,kk,numel(nn),numel(mm),numel(ss));
            tic; 
            for trial = 1:nTrials
                fprintf('.');
                
                % Sparsifying dictionary
                Dict = randstiefel(n,d);

                % Location of zero elements
                zeroLocs  = randsample(n,s);

                % Basis for nullspace
                Null = null(Dict(zeroLocs,:));

                % Project a random vector onto this nullspace:
                x0 = Null*(Null'*randn(d,1));

                % Generate a measurement operator
                M = randstiefel(d,m)';

                % Generate the measurement
                b = M*x0;

                % Solve
                cvx_begin quiet
                   
                   variable x(d,1);
                   minimize norm(Dict*x,1);
                   subject to
                      M*x == b;
                cvx_end

                if norm(x-x0) < max(relTol*norm(x0),1e-10) % Success!
                    results(jj,kk,ii) = results(jj,kk,ii)+1; % NOTE ORDER OF SUBSCRIPTS
                end
                    
            end
            toc;
        end
    end
    
    
end
% Save the results

timestamp = datestr(now);
save('test_analysis_PT_random.mat');
save(['test_analysis_PT_random ',timestamp,'.mat']);

%% Plot the results
load('test_analysis_PT_random.mat');

%%

p_sig = @(sig,lambda) exp(-lambda.^2./4./(sig.^2 + lambda/3));
tail_bd = @(sdim) 2 * (d-sdim) .* p_sig(sqrt(d-sdim),d-sdim);
clf;
for ii = 1:numel(nn)
    
    n = nn(ii);
    % ss = Number of ZEROS
    kk = n-ss; % kk = Number of NONZEROS in analysis vector

    
    subplot(2,3,ii);
    cla;
    imagesc(kk,mm,results(:,:,ii));
    axis square;
    colormap bone;
    title(sprintf('n = %d',nn(ii)));
    ylabel('m: # Measurements');
    xlabel('(s): # Non Zeros');
    axis xy
    
    hold on;
    
    % Compute the statdim
    %l1 = nn(ii)*snow('l1',(nn(ii)-(d-ss))./nn(ii)) - (nn(ii)-d);
    %l1 = d*snow('l1',(nn(ii)-ss)./nn(ii));
    
    % Intersection rule with the image
    l1 = n*snow('l1',kk./n)+(d-n); 
    
    b = 4/3*log(2)*(n-kk);
    c = 4*log(2)*(n-kk).*min(l1,n-l1);
    
    num_meas = l1 + (b +sqrt(b.^2 + 4*c))/2;
    
    %num_meas = l1 + sqrt(10/3*log(2)*max(n-kk,0)*d);
    %plot(d-ss,l1,d-ss, l1 - tail_bd(l1));
    plot(kk,num_meas,'LineWidth',2);
end
