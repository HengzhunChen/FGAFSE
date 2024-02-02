function [A, S, Q, P, DzQ, DzP] = time_evolution(A0, S0, Q0, P0, DzQ0, DzP0, dt, Nt, alpha, Odes, potential, veps)
% TIME_EVOLUTION function to evolve the FGA ODEs using 4-th order 
% Runge-Kutta (RK) method for 1-dimensional problem
%    Inputs:
%        A0, S0, Q0, P0, DzQ0, DzP0 
%                  -- initial FGA parameter arrays of size nGB
%        dt        -- step size of time t
%        Nt        -- num of time steps to evolve
%        alpha     -- order of fractional operator
%        Odes      -- function handle of the FGA ODES
%                     [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = ...
%                             Odes(Q, P, DzQ, DzP, alpha, potential)
%        potential -- function handle of potential and its derivatives
%        veps      -- scaled Planck constant
%    Outputs:
%        A, S, Q, P, DzQ, DzP 
%                  -- FGA parameter arrays of size nGB
%    A, S are scalar while Q, P may be vectors.
%
%    See also FGA1d, odes_1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


% Initialization
Qi = Q0;
Pi = P0;
Si = S0;
logAi = log(A0);
DzQi = DzQ0;
DzPi = DzP0;

% Main loop
for tt = 1 : Nt
    % find slope vector k1
    [kQ1, kP1, kS1, klogA1, kDzQ1, kDzP1] = ...
            Odes(Qi, Pi, DzQi, DzPi, alpha, potential, veps);
    
    % find slope vector k2
    Q = Qi + kQ1 * dt / 2;
    P = Pi + kP1 * dt / 2;
    DzQ = DzQi + kDzQ1 * dt / 2;
    DzP = DzPi + kDzP1 * dt / 2;
    [kQ2, kP2, kS2, klogA2, kDzQ2, kDzP2] = ...
            Odes(Q, P, DzQ, DzP, alpha, potential, veps);
    
    % find slope vector k3       
    Q = Qi + kQ2 * dt / 2;
    P = Pi + kP2 * dt / 2;
    DzQ = DzQi + kDzQ2 * dt / 2;
    DzP = DzPi + kDzP2 * dt / 2;
    [kQ3, kP3, kS3, klogA3, kDzQ3, kDzP3] = ...
            Odes(Q, P, DzQ, DzP, alpha, potential, veps);
    
    % find slope vector k4
    Q = Qi + kQ3 * dt;
    P = Pi + kP3 * dt;
    DzQ = DzQi + kDzQ3 * dt;
    DzP = DzPi + kDzP3 * dt;
    [kQ4, kP4, kS4, klogA4, kDzQ4, kDzP4] = ...
            Odes(Q, P, DzQ, DzP, alpha, potential, veps);

    % evolution
    Qi = Qi + (kQ1 + 2 * kQ2 + 2 * kQ3 + kQ4) * dt / 6;
    Pi = Pi + (kP1 + 2 * kP2 + 2 * kP3 + kP4) * dt / 6;
    Si = Si + (kS1 + 2 * kS2 + 2 * kS3 + kS4) * dt / 6;
    DzQi = DzQi + (kDzQ1 + 2 * kDzQ2 + 2 * kDzQ3 + kDzQ4) * dt / 6;
    DzPi = DzPi + (kDzP1 + 2 * kDzP2 + 2 * kDzP3 + kDzP4) * dt / 6;
    logAi = logAi + (klogA1 + 2 * klogA2 + 2 * klogA3 + klogA4) * dt / 6;        
end

Q = Qi;
P = Pi;
S = Si;
A = exp(logAi);
DzQ = DzQi;
DzP = DzPi;

end