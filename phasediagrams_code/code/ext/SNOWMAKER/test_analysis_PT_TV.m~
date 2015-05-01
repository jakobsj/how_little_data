%% Generate PT diagram for analysis L1

% Length of signal (ambient dimension)
d = 200;

% Number of elements in the sparisifying dictionary 
n= d;


% Sparsifying dictionary = first difference
Dict = diag(ones(d,1),0)+ diag(-ones(d-1,1),1);


% Sparsity (number of nonzeros)
nnz = 38;

% Number of measurements
m = 100;

% Location of zero elements
zeroLocs  = randsample(n,n-nnz);

% Basis for nullspace
Null = null(Dict(zeroLocs,:));

% Project a random vector onto this nullspace:
x = Null*(Null'*randn(d,1));

% Generate a measurement operator
M = randstiefel(d,m)';

% Generate the measurement
b = M*x;

%
cvx_begin
    variable y(d,1);
    minimize norm(Dict*y,1);
    subject to 
        M*y == b;
cvx_end

norm(y-x)/norm(x)