function [w, x] = FGA1d(alpha, vepsExp, finalTime, right_x, initWavefun, potential)
% FGA1D Solver for 1-dim fractional Schrodinger equation with high frequency
% wave using frozen Gaussian approximation (FGA) method.
%    Inputs:
%        alpha       -- order of fractional operator
%        vepsExp     -- exponent of veps (scaled Planck constant), 
%                       veps = 2 ^ vepsExp
%        finalTime   -- final time of evolution
%        right_x     -- right endpoint of domain of x
%                       left endpoint of domian of x is 0
%        initWavefun -- function handle for initial wavefunction
%                       u0 = initWavefun(x, veps)
%        potential   -- function handle of potential 
%                       [V, DV, D2V] = potential(Q)
%    Outputs:
%       x -- samples on [0, right_x] 
%       w -- solution at (x, finalTime)
%
%    See also FGA2d, initial_deomposition_1d, odes_delta_1d, wave_reconstruction_1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


veps = 2 ^ vepsExp;  % scaled Planck constant

% Setup mesh grid
dx = veps;
nx = floor( (right_x - 0) / dx);
dy = dx;  % use the same mesh size for x and y
ny = nx;

% number of y grid included in each stepsize of q, nydq := ny / nq, dq := dy * nqdq
nydq = floor( 2^(-vepsExp/2) / 2 );  
% number of points included in a Gaussian kernel
kernelSize = floor( 2^(-vepsExp/2) ) * 2^3;  

dt = 1e-2;
nt = floor( finalTime / dt + 1e-6);

% Initialization
x = 0 : dx : right_x;  % mesh on axis x, left endpoint is 0
x = x(1 : end-1)';  % shape: (nx, 1)
u0 = initWavefun(x, veps);

odes = @odes_delta_1d;

% Main loop
[A0, S0, Q0, P0, nGB] = initial_decomposition_1d(u0, veps, dy, ny, kernelSize, nydq); 
DzQ0 = ones(size(Q0));
DzP0 = -1i * ones(size(P0));
[A, S, Q, P] = time_evolution(A0, S0, Q0, P0, DzQ0, DzP0, dt, nt, alpha, odes, potential, veps);
w = wave_reconstruction_1d(veps, A, S, Q, P, nGB, x, dx, nx, kernelSize);

end