function [A, S, Q, P, nGB] = initial_decomposition_2d(u0, veps, dy, ny, GKsize, nydq)
% INITIAL_DECOMPOSITION_2D function to decompose the initial wavefunction into
% Gaussian beams. 
%    Gaussian beam: (Q, P) such that the absolute value of "amplitude" A greater
%    than a given threshold.
%    Inputs: 
%        u0     -- initial wavefunction in a uniform mesh
%        veps   -- scaled Planck constant 
%        dy     -- stepsize of y1, y2 axis
%        ny     -- num of mesh grids of y1, y2 axis
%        GKsize -- size of the Gaussian kernel, i.e.,
%                  num of y grid contained in a Gaussian kernel,
%                  at the same time: num of mesh grids of p
%        nydq   -- number of y grid included in a stepsize of q
%    Outputs:
%        nGB -- num of Gaussians beams
%        A   -- Amplitude array of size nGB, multipling steps dp, dq, scalar
%        S   -- Action array of size nGB, scalar
%        Q   -- Position array of size nGB, Q = [Q1, Q2]
%        P   -- Momentum array of size nGB, P = [P1, P2]
%
%    See aslo FGA2d, wave_reconstruction_2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


Gk1 = (-GKsize/2 : GKsize/2 - 1)' * ones(1, GKsize);
Gk2 = Gk1';
Gk_filter = exp( -( (Gk1 * dy).^2 + (Gk2 * dy).^2 ) / (2*veps) );

nq = floor( ny / nydq );  % number of mesh grids of q
dq = dy * nydq;  % stepsize of q
np = GKsize;  % number of mesh grids of p 
dp = veps * 2 * pi / (dy * GKsize);  % mesh size of p axis, stepsize of p

Qmesh1 = (0 : nq-1)' * ones(1, nq) * dq;
Qmesh2 = Qmesh1';
Pmesh1 = [0: GKsize/2 - 1, -GKsize/2 : -1]' * ones(1, GKsize) * dp;
Pmesh2 = Pmesh1';

nGBMax = 2 ^ 30;  % 2^20 == 1048576
ntemp = min(nGBMax, np * np * nq * nq);

% Compute integral with axes y1, y2 using FFT
A_pq = zeros(np, np, nq, nq);
for iq1 = 1 : nq
    i = (iq1 - 1) * nydq;
    for iq2 = 1 : nq
        j = (iq2 - 1) * nydq;

        % for each q(i, j), (il, ir) and (jl, jr) give the local kernel box
        il = max( i - GKsize/2, 0 ) + 1;  % index of left endpoint of y1
        ir = min( i + GKsize/2 - 1, ny - 1 ) + 1;  % index of right endpoint of y1
        jl = max( j - GKsize/2, 0 ) + 1;  % index of left endpoint of y2
        jr = min( j + GKsize/2 - 1, ny - 1 ) + 1;  % index of right endpoint of y2

        ubox = zeros(GKsize, GKsize);
        ubox( (il : ir) - i + GKsize/2, (jl : jr) - j + GKsize/2 ) = u0(il : ir, jl: jr);
        ubox = ubox .* Gk_filter;
        ubox = fft2(ubox) * dy * dy;  % use FFT to help simulate the integral

        psi = 2 * ubox * dq * dq * dp * dp / ((2*pi*veps) ^ 3);  % 2 comes from a(0,q,p)
        psi = psi .* exp( 1i / veps * Pmesh1 * GKsize * dy / 2 ) ...
                  .* exp( 1i / veps * Pmesh2 * GKsize * dy / 2 );
        % with fft we sum j from 0 to N-1, we change back by adding N/2 to j

        A_pq(:, :, iq1, iq2) = psi;
    end
end

% Set a threshold to discard small values
threshold = 1e-5;

AM = max(abs(A_pq), [], 'all') * threshold;
index = abs(A_pq) >= AM;

nGB = sum(index, 'all');  % number of Gaussian beams
if nGB > ntemp
    error("Exceed the maximum beam number ...");
end

A = A_pq(index);
A = A(:);

[ip1, ip2, iq1, iq2] = ind2sub(size(A_pq), find(index));
idP = sub2ind([np, np], ip1, ip2);
Q = [(iq1 - 1) * dq, (iq2 - 1) * dq];
P = [Pmesh1(idP), Pmesh2(idP)];
S = zeros(nGB, 1);


% Verifying whether initial decomposition is correct
w = wave_reconstruction_2d(veps, A, S, Q, P, nGB, dy, ny, GKsize);
errInit = sqrt( sum(sum( abs(w - u0).^2 )) * dy * dy );
fprintf('L2 error of initial decomposition: %e, number of Gaussian beams: %d\n', errInit, nGB);

end
