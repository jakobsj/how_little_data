% UTV Tools.
% Version 1.0  9-August-11
%
% Iterative ART Methods.
%   kaczmarz       - Kaczmarz's method (often referred to as ART).
%   randkaczmarz   - Randomized Kaczmarz method.
%   symkaczmarz    - Symmetric Kaczmarz method.
%
% Iterative SIRT Methods.
%   cav            - Component Averaging (CAV) method.
%   cimmino        - Cimmino's projection method.
%   drop           - Diagonally Relaxed Orthogonal Projections (DROP) method.
%   landweber      - The classical Landweber method.
%   sart            - Simultaneous Algebraic Reconstruction Technique (SART) method.
%
% Training Routines.
%   trainDPME       - Training method for the stopping rules DP and ME.
%   trainLambdaART  - Training to determine optimal lambda for ART methods.
%   trainLambdaSIRT - Training to determine optimal lambda for SIRT method.
%
% Test Problems.
%   fanbeamtomo    - Creates a 2D tomography test problem using fan beams.
%   paralleltomo    - Creates a 2D tomography test problem using parallel beams. 
%   seismictomo    - Creates a seismic tomography test problem.
%
% Demonstration Scripts.
%   ARTdemo        - Demonstrates the use of, and the results from, the ART methods.
%   nonnegdemo     - Demonstrates the use of nonnegativity constraints.
%   SIRTdemo       - Demonstrates the use of, and the results from, the SIRT methods.
%   trainingdemo   - Demonstrates the use of the training methods.
%
% Auxiliary Routines.
%   calczeta       - Calculates the roots of a specific polynomial g(z) of degree k.
%   rzr            - Remove zero rows of A and the corresponding elements of b.
