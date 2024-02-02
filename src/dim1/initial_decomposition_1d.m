function [A, S, Q, P, nGB] = initial_decomposition_1d(u0, veps, dy, ny, GKsize, nydq)
% INITIAL_DECOMPOSITION_1D function to decompose the initial wavefunction into
% Gaussian beams. 
%    Gaussian beam: (Q, P) such that the absolute value of "amplitude" A greater
%    than a given threshold.
%    Inputs: 
%        u0     -- initial wavefunction over a uniform mesh
%        veps   -- scaled Planck constant 
%        dy     -- stepsize of y axis
%        ny     -- num of mesh grids of y axis
%        GKsize -- size of the Gaussian kernel, i.e.,
%                  num of y grid contained in a Gaussian kernel,
%                  at the same time: num of mesh grids of p
%        nydq   -- number of y grid included in a stepsize of q
%    Outputs:
%        nGB -- num of Gaussian beams
%        A   -- Amplitude array of size nGB, multipling stepsize dp, dq
%        S   -- Action array of size nGB
%        Q   -- Position array of size nGB
%        P   -- Momentum array of size nGB
% 
%    See also FGA1d, wave_reconstruction_1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


% *******************************************************************************
%                     Mesh grid strategy for q, p axes
% Decomposition formula
%
%   A(t,q,p) = \int a*u0(y) * exp(-1i/veps * p*(y-q) + |y-q|^2 / (2*veps)) * dy.
%
% We firstly give a uniform grid for y consistent with the scaled Planck
% constant veps, then sample the q grid points over the y grid uniformly with
% gap dy*nydq.
%
% Note that only those y points within the kernel size of a center will be taked
% into account, number of p grid points is the size of a Gaussain Kernel. 
% Step size of p mesh is designed for converting the integral into standard form 
% of FFT, i.e.,
%
%   1i/veps * p_k * y_j = 1i/veps * k*dp * j*dy = 1i* 2*pi * k * j / N. 
%
% *******************************************************************************

nq = floor( ny / nydq );
dq = dy * nydq;  % stepsize of q
np = GKsize;
dp = veps * 2 * pi / (dy * GKsize);  % mesh size of p axis, stepsize of p

Qmesh = (0 : nq-1)' * dq;
Pmesh = [0 : GKsize/2 - 1, -GKsize/2 : -1]' * dp;

index_Gk = ( -GKsize/2 : GKsize/2 - 1 )';  % index of Gaussian kernel
Gk_filter = exp( -(index_Gk * dy).^2 / (2*veps) );


A_pq = zeros(np, nq);
% Compute integral with axis y using FFT
iq = 0;  % index of q
for i = 0 : nydq : ny-1
    % for each q(i), (il, ir) is a kernel centered at q(i)
    il = max(i - GKsize/2, 0) + 1;  % index of left endpoint
    ir = min(i + GKsize/2 - 1, ny - 1) + 1; % index of right endpoint

    ubox = zeros(GKsize, 1);
    ubox((il : ir) - i + GKsize/2) = u0(il : ir);  % ir-il may not equal to GKsize
    ubox = ubox .* Gk_filter;
    ubox = fft(ubox) * dy;  % use FFT to help simulate the integral

    psi = sqrt(2) * ubox * dq * dp / (2 * pi * veps)^(3/2);  % sqrt(2) comes from a(0,q,p)
    psi = psi .* exp(1i / veps * Pmesh * (GKsize/2) * dy) ;
    % with fft we sum j from 0 to N-1, we transform back by adding N/2 to j

    iq = iq + 1;
    if iq > nq
        break;
    end
    A_pq(:, iq) = psi;
end
Q_pq = ones(np, 1) * Qmesh';
P_pq = Pmesh * ones(1, nq);

% Set a threshold to discard small values
threshold = 1e-8;

AM = max(abs(A_pq), [], 'all') * threshold;
index = abs(A_pq) >= AM;

nGB = sum(index, 'all');  % number of Gaussian beams
A = A_pq(index);
A = A(:);
Q = Q_pq(index);
Q = Q(:);
P = P_pq(index);
P = P(:);
S = zeros(nGB, 1);

% Verify whether initial decomposition is correct
right_x = dy * ny;
y = 0 : dy : right_x;
y = y(1 : end-1)';
w0 = wave_reconstruction_1d(veps, A, S, Q, P, nGB, y, dy, ny, GKsize);
errInit = sqrt( sum( abs(w0 - u0).^2 ) * dy );
fprintf('L2 error of initial decomposition: %e\n', errInit);

end