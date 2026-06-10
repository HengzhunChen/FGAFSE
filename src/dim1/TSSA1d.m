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
%        dt          -- maximum time step
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

% Initialization
x = 0 : dx : right_x;  % mesh on axis x, left endpoint is 0
x = x(1 : end-1)';  % shape: (nx, 1)
u = initWavefun(x, veps);
V = potential(x);

% Main loop
k = [0 : nx/2 - 1, -nx/2 : -1]' * 2 * pi / (right_x - 0);
fLaplace = (abs(k) .^ alpha) / alpha;

numFullSteps = floor(finalTime / dt);
remainingTime = finalTime - numFullSteps * dt;
timeTolerance = 1e-12;
if remainingTime <= timeTolerance
    remainingTime = 0;
end
numSteps = numFullSteps + (remainingTime > 0);

% ************************************************************
% Option 1: first-order time splitting spectral approximation
% ************************************************************

% kineticPhase = exp(-1i * veps^(alpha - 1) * fLaplace * dt);
% potentialPhase = exp(-1i / veps * V * dt);
% for i = 1 : numSteps
%     if i > numFullSteps
%         kineticPhase = exp(-1i * veps^(alpha - 1) * fLaplace * remainingTime);
%         potentialPhase = exp(-1i / veps * V * remainingTime);
%     end
%
%     u = fft(u);
%     u = kineticPhase .* u;
%     u = ifft(u);
%     u = potentialPhase .* u;
% end

% ***********************************************************
% Option 2: Strang splitting spectral approximation
% ***********************************************************

kineticPhase = exp(-1i * veps^(alpha - 1) * fLaplace * dt);
potentialHalfPhase = exp(-1i / veps * V * dt / 2);
for i = 1 : numSteps
    if i > numFullSteps
        kineticPhase = exp(-1i * veps^(alpha - 1) * fLaplace * remainingTime);
        potentialHalfPhase = exp(-1i / veps * V * remainingTime / 2);
    end

    u = potentialHalfPhase .* u;
    u = fft(u);
    u = kineticPhase .* u;
    u = ifft(u);
    u = potentialHalfPhase .* u;
end

end
