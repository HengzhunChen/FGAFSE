function [u, x] = TSSA1d(alpha, vespExp, finalTime, right_x, dt, dx, initWavefun, potential)
% TSSA1D Solver for 1-dim Schrodinger equation using time splitting spectral
% approximation (TSSA) method.
%    Inputs: 
%        alpha       -- order of fractional operator
%        vepsExp     -- exponent of veps (scaled Planck constant), 
%                       veps = 2 ^ vepsExp
%        finalTime   -- final time of evolution
%        right_x     -- right endpoint of domain of x
%                       left endpoint of domain of x is 0
%        dt          -- mesh size of axis t
%        dx          -- mesh size of axis x
%        initWavefun -- function handle for initial wavefunction
%                       u0 = initWavefun(x, veps)
%        potential   -- function handle of potential 
%                       V = potential(Q)
%    Outputs:
%        x -- samples on [0, right_x] 
%        u -- solution at (x, finalTime)
%
%    See also TSSA2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


veps = 2 ^ vespExp;

nx = floor( (right_x - 0) / dx );  % number of mesh grids of x    
nt = floor( finalTime / dt + 1e-6);  % number of mesh grids of t

% Initialization
x = 0 : dx : right_x;  % mesh on axis x, left endpoint is 0
x = x(1 : end-1)';  % shape: (nx, 1)
u = initWavefun(x, veps);
V = potential(x);

% Main loop
k = [0 : nx/2 - 1, -nx/2 : -1]' * 2 * pi / (right_x - 0);
fLaplace = (abs(k) .^ alpha) / alpha;


% ************************************************************
% Option 1: first-order time splitting spectral approximation
% ************************************************************
% for i = 1 : nt
%     %-----------------------------------------
%     % step 1. The Laplace
%     %-----------------------------------------
%     u = fft(u);
%     u = exp( -1i * (veps ^ (alpha - 1)) * fLaplace * dt) .* u;
%     u = ifft(u);
%     %-----------------------------------------
%     % step 2. The external potential
%     %-----------------------------------------
%     u = exp( -1i / veps * V * dt) .* u;
% end


% ***********************************************************
% Option 2: Strang splitting spectral approximation
% ***********************************************************
for i = 1 : nt
    %-----------------------------------------
    % step 1. The external potential
    %-----------------------------------------
    u = exp( -1i / veps * V * dt / 2) .* u;
    %-----------------------------------------
    % step 2. The Laplace
    %-----------------------------------------
    u = fft(u);
    u = exp( -1i * (veps ^ (alpha - 1)) * fLaplace * dt) .* u;
    u = ifft(u);
    %-----------------------------------------
    % step 3. The external potential
    %-----------------------------------------
    u = exp( -1i / veps * V * dt / 2) .* u;
end

end