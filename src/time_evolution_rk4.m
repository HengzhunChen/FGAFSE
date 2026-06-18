function [A, S, Q, P, DzQ, DzP] = time_evolution_rk4( ...
        A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, odes)
% TIME_EVOLUTION_RK4 Evolve the FGA ODEs using fourth-order Runge-Kutta.
%    Inputs:
%        A0, S0, Q0, P0, DzQ0, DzP0
%                  -- initial FGA parameter arrays of size nGB
%        dt        -- maximum step size of time t
%        finalTime -- duration of time evolution
%        odes      -- closed function handle of the FGA ODEs
%                     [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = ...
%                             odes(Q, P, DzQ, DzP)
%    Outputs:
%        A, S, Q, P, DzQ, DzP
%                  -- FGA parameter arrays of size nGB
%
%    See also TIME_EVOLUTION, TIME_EVOLUTION_SYMPLECTIC.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


Qn = Q0;
Pn = P0;
Sn = S0;
logAn = zeros(size(A0));
DzQn = DzQ0;
DzPn = DzP0;

numFullSteps = floor(finalTime / dt);
remainingTime = finalTime - numFullSteps * dt;
timeTolerance = 1e-12;
if remainingTime <= timeTolerance
    remainingTime = 0;
end
numSteps = numFullSteps + (remainingTime > 0);

for tt = 1 : numSteps
    stepSize = dt;
    if tt > numFullSteps
        stepSize = remainingTime;
    end

    [kQ1, kP1, kS1, klogA1, kDzQ1, kDzP1] = odes(Qn, Pn, DzQn, DzPn);

    Q = Qn + kQ1 * (stepSize / 2);
    P = Pn + kP1 * (stepSize / 2);
    DzQ = DzQn + kDzQ1 * (stepSize / 2);
    DzP = DzPn + kDzP1 * (stepSize / 2);
    [kQ2, kP2, kS2, klogA2, kDzQ2, kDzP2] = odes(Q, P, DzQ, DzP);

    Q = Qn + kQ2 * (stepSize / 2);
    P = Pn + kP2 * (stepSize / 2);
    DzQ = DzQn + kDzQ2 * (stepSize / 2);
    DzP = DzPn + kDzP2 * (stepSize / 2);
    [kQ3, kP3, kS3, klogA3, kDzQ3, kDzP3] = odes(Q, P, DzQ, DzP);

    Q = Qn + kQ3 * stepSize;
    P = Pn + kP3 * stepSize;
    DzQ = DzQn + kDzQ3 * stepSize;
    DzP = DzPn + kDzP3 * stepSize;
    [kQ4, kP4, kS4, klogA4, kDzQ4, kDzP4] = odes(Q, P, DzQ, DzP);

    Qn = Qn + (kQ1 + 2 * kQ2 + 2 * kQ3 + kQ4) * (stepSize / 6);
    Pn = Pn + (kP1 + 2 * kP2 + 2 * kP3 + kP4) * (stepSize / 6);
    Sn = Sn + (kS1 + 2 * kS2 + 2 * kS3 + kS4) * (stepSize / 6);
    DzQn = DzQn + (kDzQ1 + 2 * kDzQ2 + 2 * kDzQ3 + kDzQ4) * (stepSize / 6);
    DzPn = DzPn + (kDzP1 + 2 * kDzP2 + 2 * kDzP3 + kDzP4) * (stepSize / 6);
    logAn = logAn + (klogA1 + 2 * klogA2 + 2 * klogA3 + klogA4) * (stepSize / 6);
end

Q = Qn;
P = Pn;
S = Sn;
A = A0 .* exp(logAn);
DzQ = DzQn;
DzP = DzPn;

end
