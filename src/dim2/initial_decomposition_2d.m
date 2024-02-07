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

nGBMax = 2 ^ 20;  % 2^20 == 1048576
ntemp = min(nGBMax, np * np * nq * nq);
A = zeros(ntemp, 1);
Q = zeros(ntemp, 2);  % Q = [Q1, Q2]
P = zeros(ntemp, 2);  % P = [P1, P2]    
    
% Compute integral with axes y1, y2 using FFT
thresh = 1e-4;
num = 0;
for i = 0 : nydq : ny-1
    for j = 0 : nydq : ny-1
        % for each q(i), (il, ir) is kernel center at q(i)
        il = max( i - GKsize/2, 0 ) + 1;  % index of left endpoint of y1
        ir = min( i + GKsize/2 - 1, ny - 1 ) + 1;  % index of right endpoint of y1
        jl = max( j - GKsize/2, 0 ) + 1;  % index of left endpoint of y2
        jr = min( j + GKsize/2 - 1, ny - 1 ) + 1;  % index of right endpoint of y2
        
        ubox = zeros(GKsize, GKsize);
        ubox( (il : ir) - i + GKsize/2, (jl : jr) - j + GKsize/2 ) = u0(il : ir, jl: jr);  
        % ir-il may not equal to GKsize
        ubox1 = ubox .* Gk_filter;

        u0Max = sum(sum( abs(u0(il: ir, jl: jr)).^2 )) * dy * dy;
        if u0Max <= 1e-7
            s0 = 1;
        else
            s0 = u0Max / max(1e-7, sum(sum( abs(ubox1).^2 )) * dq * dq);
        end
        scal = max(s0, 1);
        
        ubox = fft2(ubox1) * dy * dy;  % use FFT to help simulate the integral

        psi = 2 * ubox * dq * dq * dp *dp / ((2*pi*veps) ^ 3);  % 2 comes from a(0,q,p)
        psi = psi .* exp( 1i / veps * Pmesh1 * GKsize * dy / 2 ) ...
                  .* exp( 1i / veps * Pmesh2 * GKsize * dy / 2 );
        % with fft we sum j from 0 to N-1, we change back by adding N/2 to j

        iq1 = i / nydq + 1;  % index of q1
        iq2 = j / nydq + 1;  % index of q2
        if iq1 > nq || iq2 > nq
            break;
        end
        
        psiM = max( max(abs(psi), [], 'all') * thresh, 1e-7 );

        index = ( abs(psi) >= psiM * scal );
        length = sum(sum(index));
        if (num + length) > ntemp
            error("Exceed the maximum beam number ...");
        end
        A(num + 1 : num + length) = psi(index);
        Q(num + 1 : num + length, 1) = Qmesh1(iq1, iq2);
        Q(num + 1 : num + length, 2) = Qmesh2(iq1, iq2);
        P(num + 1 : num + length, 1) = Pmesh1(index);
        P(num + 1 : num + length, 2) = Pmesh2(index);
        num = num + length;
    end
end

nGB = num;
A = A(1 : nGB);
Q = Q(1 : nGB, :);
P = P(1 : nGB, :);
S = zeros(nGB, 1);


% Verifying whether initial decomposition is correct
w = wave_reconstruction_2d(veps, A, S, Q, P, nGB, dy, ny, GKsize);
errInit = sqrt( sum(sum( abs(w - u0).^2 )) * dy * dy );
fprintf('L2 error of initial decomposition: %e, number of Gaussian beams: %d\n', errInit, nGB);

end