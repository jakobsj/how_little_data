function [mask, N_pix] = build_mask(N,radius)
%BUILD_MASK    Create logical array of disk-shaped mask
% [mask, N_pix] = build_mask(N,radius)
% 
% Inputs:
%  - N          Side length of square domain.
%  - radius     Relative radius, the maximal value of 1 gives a disk
%               inscribed in the square domain, tangent to the borders
%               (optional, default=1).
% Outputs:
%  - mask       Logical array holding mask.
%  - N_pix      Number of "true" entries, ie pixels inside the mask.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

% Set default radius.
if nargin < 2
    radius = 1;
end

%% Mask for test image
mx = -N/2+0.5:N/2-0.5;
[MX,MY] = meshgrid(mx/(N/2));
mask = MX.^2 + MY.^2 <= radius^2;
N_pix = sum(mask(:));