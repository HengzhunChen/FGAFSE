function w = wave_reconstruction_2d(veps, A, S, Q, P, nGB, dx, nx, GKsize)
% WAVE_RECONSTRUCTION_2D function to sum all the Gaussian beams up and return 
% a wavefunction.
%    Inputs:
%        veps   -- scaled Planck constant
%        nGB    -- num of Gaussian beams
%        A, S, Q, P 
%               -- FGA parameter arrays of size nGB
%        dx     -- mesh step size in x1 and x2
%        nx     -- num of mesh grids in x1 and x2
%        GKsize -- size of the Gaussian kernel, i.e.,
%                  num of y grid contained in a Gaussian kernel,
%    Output:
%        w      -- wavefunction, matrix of size nx * nx, w = w(x)
%
%    See also FGA2d, initial_decomposition_2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


% For every [Q1(i), Q2(i)], find out all the x grids included in it kernel range
% and sum up into w(x).

right_x = dx * nx;
w = zeros(nx, nx);

for i = 1 : nGB
    Q1 = Q(i, 1);
    Q2 = Q(i, 2);
    if (Q1 < 0 || Q1 > right_x || Q2 < 0 || Q2 > right_x)  
        % Q(i) should lie in the domain of x
        continue  
    end    

    Q1_id = floor( Q1 / dx + 1e-6 );  % index of Q(i) in mesh grid of x
    Q2_id = floor( Q2 / dx + 1e-6 );
    il = max( -GKsize/2, 0 - Q1_id ) + 1;
    ir = min( GKsize/2 - 1, nx - 1 - Q1_id ) + 1;
    jl = max( -GKsize/2, 0 - Q2_id ) + 1;
    jr = min( GKsize/2 - 1, nx - 1 - Q2_id ) + 1;
    id_x1 = Q1_id + (il : ir);  % index set of x whose kernel including Q(i)
    id_x2 = Q2_id + (jl : jr);

    %!!! NOTE: index is not equal to mesh multiplying step size
    X1 = (id_x1 - 1)' * ones(size(id_x2)) * dx;  % shape: (col(id_x1), col(id_x2))
    X2 = dx * ones(size(id_x1))' * (id_x2 - 1);  % shape: (col(id_x1), col(id_x2))

    w(id_x1, id_x2) = w(id_x1, id_x2) + A(i) * ...
        exp( -((X1 - Q1) .^2) / (2*veps) - ((X2 - Q2) .^2) / (2*veps) + ...
             1i / veps * ( S(i) + P(i, 1) * (X1 - Q1) + P(i, 2) * (X2 - Q2) ) );
end

end