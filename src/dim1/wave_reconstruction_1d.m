function w = wave_reconstruction_1d(veps, A, S, Q, P, nGB, x, dx, nx, GKsize)
% WAVE_RECONSTRUCTION_1D function to sum all the Gaussian beams up and return 
% a wavefunction.
%    Inputs:
%        veps   -- scaled Planck constant
%        nGB    -- num of Gaussian beams
%        A, S, Q, P 
%               -- FGA parameter arrays of size nGB
%        x      -- a uniform mesh grids where the wavefunction evaluated
%        dx     -- mesh step size
%        nx     -- num of mesh grids
%        GKsize -- size of the Gaussian kernel, i.e.,
%                  num of y grid contained in a Gaussian kernel,
%    Output:
%        w      -- wavefunction, array of size nx, w = w(x)
%
%    See also FGA1d, initial_decomposition_1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


% For each Q(i), find out all the x grids included in it kernel range
% and sum up into w(x).

right_x = dx * nx;
w = zeros(nx, 1);
for i = 1 : nGB
    if (Q(i) < 0 || Q(i) > right_x)  % Q(i) lies in domain of x
        continue
    end
    Q_id = floor( Q(i) / dx + 1e-6 );  % index of Q(i) in mesh grid of x
    il = max( -GKsize/2, 0 - Q_id ) + 1;  % left part
    ir = min( GKsize/2 - 1, nx - 1 - Q_id ) + 1;  % right part
    id_x = Q_id + (il : ir);  % index set of x whose kernel including Q(i)
    % center at Q(i), left part = right part = kernel / 2
    w(id_x) = w(id_x) + ...
                A(i) * exp( -((x(id_x) - Q(i)) .^2) / (2 * veps) + ...
                1i / veps * (S(i) + P(i) * (x(id_x) - Q(i))) );
end

end