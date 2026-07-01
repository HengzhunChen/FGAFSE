function [ww, xx] = FGA2d(alpha, vepsExp, finalTime, right_x, initWave, potential_2d, delta, solver)
% FGA2D Solver for 2-dim fractional Schrodinger equation with high frequency
% wave using frozen Gaussian approximation (FGA) method.
%     Inputs:
%         alpha       -- order of fractional operator       
%         vepsExp     -- exponent of veps (scaled Plank constant),
%                        veps = 2 ^ vepsExp
%         finalTime   -- final time of evolution
%         right_x     -- right endpoint of domain of x1(x2)
%         initWavefun -- function handle for initial wavefunction
%                        u0 = initWavefun(X, Y, veps)
%         potential_2d
%                     -- function handle of potential
%                        [V, DV, D2V] = potential_2d(Q1, Q2)
%         delta       -- regularization parameter
%                        default: delta = veps
%         solver      -- ODE solver: 'rk4' (default) or 'symplectic'
%     Outputs:
%         xx -- coordinates of x1 axis on 2 dimension mesh grid(i.e., matrix)
%               coordinates of x2 axis is xx';
%         ww -- solution at (x1, x2, finalTime) 
%
%    See also FGA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


veps = 2 ^ vepsExp;  % veps: varepsilon, the scaled Planck constant 

if nargin < 7 || isempty(delta)
    delta = veps;
end

if nargin < 8 || isempty(solver)
    solver = 'rk4';
end
switch solver
    case 'rk4'
        dt = 1e-3;
    case {'symplectic', 'stormer-verlet'}
        dt = 1e-3;
end

% Setup mesh grid
dx = veps;
nx = floor( (right_x - 0) / dx );  % number of mesh grids of x1, x2
dy = dx;  % mesh size of y1, y2 axis, use same mesh size for x and y
ny = nx;  % number of mesh grids of y1, y2

% number of y grid included in each stepsize of q, nydq := ny / dq, dq := ny * nydq
nydq = max(1, floor( 2^(-vepsExp/2) / 2 ));
% numer of points included in a Gaussian kernel
kernelSize = floor( 2^(-vepsExp/2) ) * 2^3;

% Initialization

% *****************************************************************************
% Sample mesh: x1 axis and x2 axis use the same grid size
%
% [0, 0, 0, ..., 0]          [0, 1, 2, ..., nx]
% [1, 1, 1, ..., 1]     +    [0, 1, 2, ..., nx]  = X + Y
%         ...                         ...
% [nx, nx, ..., nx]          [0, 1, 2, ..., nx]
%
% i.e., from xx and xx' we can get the mesh we want.
% Here, we suppose rows of matrix as x1 axis and columns as x2 axis.
% However, in meshgrid() function, 
% it uses columns of matrix to denote x1 axis and rows to denote x2 axis.
% *****************************************************************************

x = linspace(0, right_x, nx+1);  % shape: (1, nx+1)
x = x(1 : end-1)';  % shape: (nx, 1)
xx = x * ones(1, nx);
u0 = initWave(xx, xx', veps);

% Main loop
[A0, S0, Q0, P0, nGB] = initial_decomposition_2d( ...
    u0, veps, dy, ny, kernelSize, nydq);
DzQ0 = repmat([1, 0, 0, 1], nGB, 1);
DzP0 = repmat( -1i * [1, 0, 0, 1], nGB, 1);

odes = @(Q, P, DzQ, DzP) ...
    odes_delta_2d(Q, P, DzQ, DzP, alpha, delta, potential_2d);
evalKinetic = @(P) kinetic_delta_2d(P, alpha, delta);
evalPotential = @(Q) potential_2d(Q(:, 1), Q(:, 2));

[A, S, Q, P] = time_evolution(A0, S0, Q0, P0, DzQ0, DzP0, ...
    dt, finalTime, odes, solver, evalKinetic, evalPotential);

ww = wave_reconstruction_2d(veps, A, S, Q, P, nGB, dx, nx, kernelSize);

end
