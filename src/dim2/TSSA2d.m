function [uu, xx] = TSSA2d(alpha, vepsExp, finalTime, right_x, dt, initWavefun, potential)
% TSSA2D Solver fo 2-dim Schrodinger equation using time splitting spectral 
% approximation (TSSA) method.
%     Inputs:
%         alpha       -- order of fractional operator
%         vepsExp     -- exponent of veps (scaled Plank constant),
%                        veps = 2 ^ vepsExp
%         finalTime   -- final time of evolution
%         right_x     -- right endpoints of domain of x1 and x2
%                        left endpoints of domain of x1 and x2 is 0
%         dt          -- mesh size of axis t
%         initWavefun -- function handle for initial wavefunction
%                        u0 = initWavefun(X, Y, veps)
%         potential   -- function handle of potential 
%                        V = potential(Q1, Q2)
%     Outputs:
%         xx -- coordinates of x1 axis on 2 dimensional mesh grid(i.e, matrix)
%               coordinates of x2 axis is xx';
%         uu -- solution at (x1, x2, finalTime)
%
%     See also TSSA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


veps = 2 ^ vepsExp;

dx = veps;
nx = floor((right_x - 0) / dx);
nt = floor( finalTime / dt + 1e-6);

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
xx = x * ones(1, nx);  % shape: (nx, nx)

[V, ~, ~] = potential(xx, xx');
uu = initWavefun(xx, xx', veps);


% Main loop
p = [0 : nx/2 - 1, -nx/2 : -1] * 2 * pi / (right_x - 0);  % shape: (nx, 1)
pp = p' * ones(1, nx);  % shape: (nx, nx)
% fLaplace = (abs(pp) .^ alpha + abs(pp') .^ alpha) / alpha;
fLaplace = (sqrt(pp .^ 2 + pp' .^ 2) .^ alpha) / alpha;


% *****************************************************************************
% Option 1: first-order time splitting spectral approximation
% *****************************************************************************
% for i = 1 : nt
%     %-----------------------------------------
%     % step 1. The Schrodinger
%     %-----------------------------------------
%     uu = fft2(uu);
%     uu = exp( -1i * (veps ^ (alpha - 1)) * fLaplace * dt ) .* uu;
%     uu = ifft2(uu);
%     %-----------------------------------------
%     % step 2. The external potential
%     %-----------------------------------------
%     uu = exp( -1i / veps * V * dt) .* uu;
% end

% *****************************************************************************
% Option 2: Strang splitting spectral approximation
% *****************************************************************************
for i = 1 : nt
    %-----------------------------------------
    % step 1. The external potential
    %-----------------------------------------
    uu = exp( -1i / veps * V * dt / 2) .* uu;
    %-----------------------------------------
    % step 2. The Schrodinger
    %-----------------------------------------
    uu = fft2(uu);
    uu = exp( -1i * (veps ^ (alpha - 1)) * fLaplace * dt) .* uu;
    uu = ifft2(uu);
    %-----------------------------------------
    % step 3. The external potential
    %-----------------------------------------
    uu = exp( -1i / veps * V * dt / 2) .* uu;
end

end