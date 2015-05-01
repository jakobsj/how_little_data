function A = paralleltomo_randomrays(N,theta,x0,isDisp)
%PARALLELTOMO_RANDOMRAYS Creates a 2D tomography matrix with random rays.
%
%   [A b x theta p d] = paralleltomo_randomrays(N,theta,x0)
%   [A b x theta p d] = paralleltomo_randomrays(N,theta,x0,isDisp)
%
% This function creates a 2D tomography test problem with an N-times-N
% domain. The vectors theta and x0, which must have same length, specify
% angle and distance to the origin of an individual ray.
%
% Input:
%   N           Scalar denoting the number of discretization intervals in
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the angles in degrees.
%   x0          Vector containing the distance to origin of the rays.
%   isDisp      If isDisp is non-zero it specifies the time in seconds
%               to pause in the display of the rays. If zero (the default),
%               no display is shown.
%
% Output:
%   A           Coefficient matrix with N^2 columns and nA*p rows,
%               where nA is the number of angles, i.e., length(theta).
%
% See also: fanbeamtomo, seismictomo.
%
% This function is modified from the function paralleltomo from AIR Tools,
% http://www2.compute.dtu.dk/~pcha/AIRtools/
%
% Jakob Heide J�rgensen, Maria Saxild-Hansen and Per Christian Hansen,
% June 21, 2011, DTU Informatics.
%
% Reference: A. C. Kak and M. Slaney, Principles of Computerized
% Tomographic Imaging, SIAM, Philadelphia, 2001.
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

% Default illustration:
if nargin < 4 || isempty(isDisp)
    isDisp = 0;
end

if length(theta) ~= length(x0)
    error('Number of angular and line coordinates must be the same')
end
num_rays = length(theta);

% The starting values both the x and the y coordinates.
y0 = zeros(num_rays,1);

% The intersection lines.
x = (-N/2:N/2)';
y = x;

% Initialize vectors that contains the row numbers, the column numbers and
% the values for creating the matrix A effiecently.
rows = zeros(2*N*num_rays,1);
cols = rows;
vals = rows;
idxend = 0;

% Prepare for illustration
if isDisp
    AA = rand(N);
    figure
end

% Illustration of the domain
if isDisp
    clf
    pause(isDisp)
    imagesc((-N/2+.5):(N/2-0.5),(-N/2+.5):(N/2-0.5),AA), colormap gray,
    hold on
    axis xy
    axis equal
    axis([-N/2-1 N/2+1 -N/2-1 N/2+1])
end

% Loop over the rays
for i = 1:num_rays
    cur_theta = theta(i);
    cur_x0 = x0(i);
    cur_y0 = y0(i);
    
    % All the starting points for the current angle.
    x0theta = cosd(cur_theta)*cur_x0-sind(cur_theta)*cur_y0;
    y0theta = sind(cur_theta)*cur_x0+cosd(cur_theta)*cur_y0;
    
    % The direction vector for all the rays corresponding to the current
    % angle.
    a = -sind(cur_theta);
    b = cosd(cur_theta);
    
    % Use the parametrisation of line to get the y-coordinates of
    % intersections with x = k, i.e. x constant.
    tx = (x - x0theta)/a;
    yx = b*tx + y0theta;
    
    % Use the parametrisation of line to get the x-coordinates of
    % intersections with y = k, i.e. y constant.
    ty = (y - y0theta)/b;
    xy = a*ty + x0theta;
    
    % Illustration of the rays
    if isDisp
        
        plot(x,yx,'-','color',[220 0 0]/255,'linewidth',1.5)
        plot(xy,y,'-','color',[220 0 0]/255,'linewidth',1.5)
        
        set(gca,'Xticklabel',{})
        set(gca,'Yticklabel',{})
        pause(isDisp)
    end
    
    % Collect the intersection times and coordinates.
    t = [tx; ty];
    xxy = [x; xy];
    yxy = [yx; y];
    
    % Sort the coordinates according to intersection time.
    [t I] = sort(t);
    xxy = xxy(I);
    yxy = yxy(I);
    
    % Skip the points outside the box.
    I = (xxy >= -N/2 & xxy <= N/2 & yxy >= -N/2 & yxy <= N/2);
    xxy = xxy(I);
    yxy = yxy(I);
    
    % Skip double points.
    I = (abs(diff(xxy)) <= 1e-10 & abs(diff(yxy)) <= 1e-10);
    xxy(I) = [];
    yxy(I) = [];
    
    % Calculate the length within cell and determines the number of
    % cells which is hit.
    d = sqrt(diff(xxy).^2 + diff(yxy).^2);
    numvals = numel(d);
    
    % Store the values inside the box.
    if numvals > 0
        
        % If the ray is on the boundary of the box in the top or to the
        % right the ray does not by definition lie with in a valid cell.
        if ~((b == 0 && abs(y0theta - N/2) < 1e-15) || ...
                (a == 0 && abs(x0theta - N/2) < 1e-15)       )
            
            % Calculates the midpoints of the line within the cells.
            xm = 0.5*(xxy(1:end-1)+xxy(2:end)) + N/2;
            ym = 0.5*(yxy(1:end-1)+yxy(2:end)) + N/2;
            
            % Translate the midpoint coordinates to index.
            col = floor(xm)*N + (N - floor(ym));
            
            % Create the indices to store the values to vector for
            % later creation of A matrix.
            idxstart = idxend + 1;
            idxend = idxstart + numvals - 1;
            idx = idxstart:idxend;
            
            % Store row numbers, column numbers and values.
            %rows(idx) = (i-1)*p + j;
            rows(idx) = i;
            cols(idx) = col;
            vals(idx) = d;
            
        end
    end
    
end

% Truncate excess zeros.
rows = rows(1:idxend);
cols = cols(1:idxend);
vals = vals(1:idxend);

% Create sparse matrix A from the stored values.
A = sparse(rows,cols,vals,num_rays,N^2);
