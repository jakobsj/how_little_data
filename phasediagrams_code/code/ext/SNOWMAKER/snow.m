function sdim = snow(name,varargin)
% SNOW: Compute statistical dimensions for various cones
% 
% USAGE: SDIM = SNOW("NAME",Param1,Param2,...) 
% 
% INPUTS:
%  -NAME  The name determines the type of cone, e.g., circular, orthant, l1
%         descent cone, etc. The implemented methods appear below.
%
%  -Param Some cones require parameters, such as sparsity levels, angles, 
%         or relative ranks.  See below for more details about the argument 
%         that each name requires.
%
%
% OUTPUTS:
%  -SDIM  The relative statistical dimension, so that
%                  SDIM*D = statistical dimension + O(1),
%         where D is the ambient dimension. The methods we implement are
%         typically accurate up to first order, but check with the
%         reference below for details.
%
% IMPLEMENTED METHODS:
%  
%  L1:    The statistical dimension of the descent cone to the l1 norm at 
%         sparse vectors.  Takes one parameter, the relative sparsity. The
%         parameter can be a vector containing multiple values of relative
%         sparsities (see example below).
%
%  S1:    The statistical dimension of the descent cone to the Schatten
%         1-norm, also known as the "nuclear" norm. Takes two parameters,
%         the relative rank (between zero and one) and the aspect ratio of
%         the matrix (m/n for an m x n matrix). The parameters can be 
%         vectors, in which case the output is an array consiting of the
%         statistical dimension for each aspect ratio/relative rank pair
%         (see the example below).
%
%  circ:  The statistical dimension of a circular cone.  Takes one or two
%         parameters, the base angle and the ambient dimension (optional).
%         Including the ambient dimension results in a small correction
%         term that makes the approximation much more accurate,
%         particularly for small ambient dimensions; the dimension infinity
%         is allowed, and is equivalent to not specifying the dimension.
%         The parameters can be vectors, in which case the output is an
%         array consisting of statistical dimensions for each
%         angle/dimension pair (see the example below).
%         
%
%         
% EXAMPLES:
%
%          % Circular cones
%          angles = linspace(0,pi/2,100);
%          D = [10,50,inf]; % Dimensions
%          sdim = snow('circ',angles,D);
%          
%          subplot(1,3,1);
%          plot(angles,sdim(:,1),angles,sdim(:,2),angles,sdim(:,3));         
%          xlabel('Angle of cone'); ylabel('Relative stat. dim.');
%          title('Stat. dim. vs. angle for circular cone');
%          legend('D=10','D=50','D -> \infty','Location','SE');
%          axis([0,pi/2,0,1])
%          axis square;
% %%      
%          % The descent cone of the l1 and l1+ norm at sparse vectors
%          rel_sparsity = linspace(0,1,100);
%          sdim = snow('l1',rel_sparsity);
%          sdimp = snow('l1+',rel_sparsity);
%          subplot(1,3,2);
%          plot(rel_sparsity,sdim,rel_sparsity,sdimp);
%          xlabel('Relative sparsity'); ylabel('Relative stat. dim.'); 
%          title('Stat. dim. vs. sparsity for l1 norm');
%          legend('L1','L1+','Location','SE');
%          axis square;
%          
% %%         
%          % The descent cone of the Schatten 1 norm
%          aspect_ratio = [0,1/2,1];
%          rel_rank = linspace(0,1,100);
%          sdim = snow('S1',rel_rank,aspect_ratio);
%          subplot(1,3,3);
%          plot(rel_rank,sdim(:,1), rel_rank,sdim(:,2), rel_rank,sdim(:,3));         
%          xlabel('Relative rank'); ylabel('Relative stat. dim.');
%          title('Stat. dim. vs. rank for S1 norm');
%          legend('Rel. rank -> 0','Rel. rank = 1/2','Rel. rank = 1','Location','SE');
%          axis square;
%          
%
%
% KNOWN ISSUES:
%
%  At extreme parameter values---usually near zero or one---the results 
%  may be inaccurate or fail due to numerical quadrature issues.
%


% Copyright (c) 2014, Michael B. McCoy.  
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met: 
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution. 
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 


% VERSION 0.1.1.1

% Thanks to John Bruer for constructive comments

switch lower(name)
    case {'lone','l1','ell 1','1'}
        taus = varargin{1};
        sdim = sdim_lone(taus,false); % "false" means no positivity restriction

    case {'lone+','l1+','ell 1 +','1+'}
        taus = varargin{1};
        sdim = sdim_lone(taus,true); % "true means positivity restriction
        
    case {'orthant','self-dual','sdp'}
        sdim = 1/2;

    case {'s1','nuclear','*','schatten-1','schatten'}
        rhos = varargin{1};
        nus = varargin{2};
        sdim = sdim_sone(rhos,nus);
    
    case {'circular','circ'}
        alphas = varargin{1};
        
        if numel(varargin)>1
            D = varargin{2};
        else
            D = inf;
        end
        sdim = sdim_circ(alphas,D);
        
            
        
    otherwise
        error('SNOWMAKER:BadMethod','Undefined method: %s',name);
end

end % End of main function

%%
function [sdim_norm,opt_loc]=sdim_lone(taus,nonneg)
% SDIM_LONE The statistical dimension of l1 descent cone.
%
% INPUTS (*=optional)
%   taus:    Vector of normalized sparsities (k/d) between zero and one.
%   nonneg*: Compute value for nonnegative sparse vectors (default = FALSE)
% 
% OUTPUTS:
%   sdim_norm: The normalized statistic dimension (sdim_norm = sdim/d).
%   opt_loc*: The optimal integration parameter

if nargin < 2
    nonneg = false;
end


sdim_norm= zeros(size(taus));
opt_loc = zeros(size(taus));

for ii = 1:numel(taus)
    
    tau = taus(ii);
    
    if nonneg == false
        fcn = @(t) tau/(1-tau) - sqrt(2/pi)/t* exp(-t^2/2) + erfc(t/sqrt(2));
    else
        fcn = @(t) tau/(1-tau) - sqrt(1/2/pi)/t* exp(-t^2/2) + (1/2)*erfc(t/sqrt(2));
    end
    
    lowerBd = erfcinv(min(sqrt(2)/(1-tau),1-eps(1)));  % Nearly tight as tau->1
    upperBd = sqrt(2/pi)*(1/tau-1); % Probably not tight as tau->0.
    
    try
        %TODO: Deal with some special known cases, use approximations for
        %out-of-bound numbers, etc.
        if tau < eps(1)
            sdim_norm(ii) = 0;
        elseif tau > 1-eps
            sdim_norm(ii) = 1;
        else % Actually do the computation then
            t = fzero(fcn,[lowerBd,upperBd]);
            if nonneg == false
                sdim_norm(ii) = tau*(1+t^2) + (1-tau) * (1+t^2) * erfc(t/sqrt(2)) ...
                    - (1-tau) * t *sqrt(2/pi)*exp(-t^2/2);
            else
                sdim_norm(ii) = tau*(1+t^2) + (1-tau) * (1+t^2) * (1/2)* erfc(t/sqrt(2)) ...
                    - (1-tau) * t *sqrt(1/2/pi)*exp(-t^2/2);
            end
        end
        
    catch errId
        warning(errId.identifier,errId.message);
        sdim_norm(ii) = NaN;
    end
    
    
    if nargout > 1
        opt_loc(ii) = t;
    end
end
end % End of sdim_lone function

%%
function opt_val=sdim_sone(rhos,nus)
% SDIM_SONE Numerically compute Schatten-1 norm statistic dimension
%
% (*=optional)
%
% Input: 
%   nu:  Aspect ratio of matrix (p/q for a p x q matrix)
%   rho: Normalized rank (can be a vector).  Takes values between 0 and 1
%
% Output: 
%   opt_loc: Location of minimizing parameter


nus = nus(:);
if min(nus(:))<0
    error('SNOWMAKER:SONE:NegativeAspectRatio',...
        'Aspect ratio must be greater than zero');
end


% Make sure that the aspect ratio lies between zero and one.  Results in
% the same answer.
nus = min(nus(:),1./nus(:));

rhos = rhos(:);

if min(nus(nus>0))<1e-2
    warning('SNOWMAKER:ExtremeValue',...
        strcat('One of the aspect ratios is very small,',...
        'which may lead to inaccurate computations.',...
        ' Consider replacing with limiting value nu = 0.'));
end


opt_val = zeros(length(rhos),length(nus));
%opt_loc = zeros(length(rhos),length(nus));

% Define auxiliary functions, to be integrated later
f1 = @(x,t)(sqrt(x./t)-1);
f2 = @(x,t)(sqrt(x)-sqrt(t)).^2;
dens = @(x,a,b)sqrt((b-x).*max((x-a),0))./max(x,eps(1)); % Regularize

for jj = 1:numel(nus)
    nu = nus(jj);
    
    for ii = 1:numel(rhos)
        rho = rhos(ii)*nu;
        if nu < eps(1)
            opt_val(ii,jj) = 2*rhos(ii) - rhos(ii).^2;
        else
            
    
            if rho < eps(1)
                opt_val(ii,jj) = 0;
            elseif rho > nu - eps(1)
                opt_val(ii,jj) = 1;
            else
                % Integration boundaries for Marcenko-Pastur integral
                y = (nu-rho)./(1-rho);
                a = (1-sqrt(y))^2;
                b = (1+sqrt(y))^2;
                rhs = 2*pi*rho./(1-rho);% Right-hand side value (integrand must equal this)
                
                tmpfun = @(x,t)f1(x,t).*dens(x,a,b);
                f = @(t)(rhs-quadgk(@(x)tmpfun(x,t),max(t,a),b));
                
                try
                    z = fzero(@(t)f(t),[eps(1),b]);
                catch errId
                    
                    if strcmp(errId.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign')
                        warning('SNOWMAKER:ExtremeValue',...
                            'Unable to complete root finding for realtive rank. Try a value farther from zero or one');
                        z = NaN;
                    else
                        % Rethrow
                        throw(errId);
                    end
                    
                end
                
                %opt_loc(ii,jj) = z;
                
                if isnan(z)
                    opt_val(ii,jj) = NaN;
                else
                    % Compute upper bound on statdim
                    tmpint = quad(@(x)f2(x,z).*dens(x,a,b),z,b);
                    opt_val(ii,jj) = (rho*(nu+1-rho+(1-rho)*z)+...
                        (1-rho)^2*tmpint/(2*pi))/nu;
                end
                
            end
        end
    end
end
end % End sdim_sone function

%%
function sdim = sdim_circ(alphas,D)
% SDIM_CIRC Calculate the statistical dimension of circular cones
%
% Input:
%   alphas: Angle or angles of cones. Either scaler or vector. 
%   D:      Dimension of cone. Scalar or vector, Inf allowed.  
%
% Output:
%   sdim:   The normalized statistical dimension, arrays of size
%           length(alphas) by length(D).  
%   
alphas = alphas(:);
D = D(:);

sdim = zeros(length(alphas),length(D));

for ii = 1:numel(D)
    d = D(ii);
    
    if d < inf
        %keyboard;
        for jj = 1:numel(alphas)
            ang = alphas(jj);
            
            f = @(phi) 1.*(phi <= ang)+ ...
                cos(phi-ang).^2 .*(ang <= phi).*(phi <= pi/2 + ang);
            integrand = @(phi) f(phi).*sin(phi).^(d/ii-2);

            sdim(jj,ii) =  (1/d)*(...
                ii* exp(log(d/ii/sqrt(pi))+gammaln(d/2/ii)...
                -gammaln((d/ii-1)/2))*quadgk(integrand,0,pi/2+ang));
        end
        
    else % d == infinity
        sdim(:,ii) = sin(alphas).^2;
    end
    
end

end