function [A, S, Q, P, DzQ, DzP] = time_evolution_symplectic( ...
        A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, evalKinetic, evalPotential)
% TIME_EVOLUTION_SYMPLECTIC Evolve the FGA ODEs by Stormer-Verlet splitting.
%    Inputs:
%        A0, S0, Q0, P0, DzQ0, DzP0
%                  -- initial FGA parameter arrays of size nGB
%        dt        -- maximum step size of time t
%        finalTime -- duration of time evolution
%        evalKinetic
%                  -- closed function handle of kinetic data
%                     [T, DT, D2T] = evalKinetic(P)
%        evalPotential
%                  -- closed function handle of potential data
%                     [V, DV, D2V] = evalPotential(Q)
%                     It should adapt the dimension-specific user potential
%                     to accept Q as an nGB-by-d array.
%    Outputs:
%        A, S, Q, P, DzQ, DzP
%                  -- FGA parameter arrays of size nGB
%
%    The ray and tangent variables follow a kick-drift-kick map.
%    The action is advanced by the compatible splitting quadrature,
%    and log(A) by a trapezoidal rule.
%
%    See also TIME_EVOLUTION, TIME_EVOLUTION_RK4.

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

    [Qn, Pn, Sn, logAn, DzQn, DzPn] = stormer_verlet_step( ...
        Qn, Pn, Sn, logAn, DzQn, DzPn, stepSize, evalKinetic, evalPotential);
end

Q = Qn;
P = Pn;
S = Sn;
A = A0 .* exp(logAn);
DzQ = DzQn;
DzP = DzPn;

end

function [Q_new, P_new, S_new, logA_new, DzQ_new, DzP_new] = stormer_verlet_step( ...
        Q_old, P_old, S_old, logA_old, DzQ_old, DzP_old, h, evalKinetic, evalPotential)

    [~, ~, D2T_old] = evalKinetic(P_old);
    [V_old, DV_old, D2V_old] = evalPotential(Q_old);
    DlogA_old = log_amplitude_rate(DzQ_old, DzP_old, D2T_old, D2V_old);

    % Kick-drift-kick:
    %   P_{n+1/2} = P_n - h/2 * d_Q V(Q_n)
    %   Q_{n+1}   = Q_n + h * d_P T(P_{n+1/2})
    %   P_{n+1}   = P_{n+1/2} - h/2 * d_Q V(Q_{n+1})
    % with the tangent variables advanced by the corresponding Hessians.
    P_half = P_old - (h / 2) * DV_old;
    DzP_half = DzP_old ...
        - (h / 2) * multiply_flattened_matrices(DzQ_old, D2V_old);

    [T_half, DT_half, D2T_half] = evalKinetic(P_half);
    Q_new = Q_old + h * DT_half;
    DzQ_new = DzQ_old + h * multiply_flattened_matrices(DzP_half, D2T_half);

    [V_new, DV_new, D2V_new] = evalPotential(Q_new);
    P_new = P_half - (h / 2) * DV_new;
    DzP_new = DzP_half ...
        - (h / 2) * multiply_flattened_matrices(DzQ_new, D2V_new);

    [~, ~, D2T_new] = evalKinetic(P_new);
    DlogA_new = log_amplitude_rate(DzQ_new, DzP_new, D2T_new, D2V_new);
    kineticAction = sum(P_half .* DT_half, 2) - T_half;
    S_new = S_old + h * kineticAction - (h / 2) * (V_old + V_new);
    logA_new = logA_old + (h / 2) * (DlogA_old + DlogA_new);

end

function DlogA = log_amplitude_rate(DzQ, DzP, D2T, D2V)
    DtDzQ = multiply_flattened_matrices(DzP, D2T);
    DtDzP = -multiply_flattened_matrices(DzQ, D2V);
    Z = DzQ + 1i * DzP;
    DZ = DtDzQ + 1i * DtDzP;

    DlogA = 0.5 * trace_inverse_product_flattened(Z, DZ);
end

function C = multiply_flattened_matrices(A, B)
    % Rows store row-major flattened matrices; compute C_i = A_i * B_i.
    dimension = sqrt(size(A, 2));
    if dimension == 1
        C = A .* B;
    elseif dimension == 2
        C = zeros(size(A));
        C(:, 1) = A(:, 1) .* B(:, 1) + A(:, 2) .* B(:, 3);
        C(:, 2) = A(:, 1) .* B(:, 2) + A(:, 2) .* B(:, 4);
        C(:, 3) = A(:, 3) .* B(:, 1) + A(:, 4) .* B(:, 3);
        C(:, 4) = A(:, 3) .* B(:, 2) + A(:, 4) .* B(:, 4);
    else
        error('FGAFSE:UnsupportedDimension', ...
            'The symplectic solver currently supports dimensions 1 and 2.');
    end
end

function value = trace_inverse_product_flattened(A, B)
    % Rows store row-major flattened matrices; compute trace(A_i^{-1} B_i).
    dimension = sqrt(size(A, 2));
    zTolerance = 1e-12;
    if dimension == 1
        if any(isnan(A) | isinf(A), 'all') || any(abs(A) < zTolerance, 'all')
            error('FGAFSE:SingularZ', 'Numerical loss of invertibility in Z.');
        end
        value = B ./ A;
    elseif dimension == 2
        detA = A(:, 1) .* A(:, 4) - A(:, 2) .* A(:, 3);
        if any(isnan(detA) | isinf(detA), 'all') || any(abs(detA) < zTolerance, 'all')
            error('FGAFSE:SingularZ', 'Numerical loss of invertibility in Z.');
        end
        value = ( A(:, 4) .* B(:, 1) - A(:, 2) .* B(:, 3) ...
            - A(:, 3) .* B(:, 2) + A(:, 1) .* B(:, 4) ) ./ detA;
    else
        error('FGAFSE:UnsupportedDimension', ...
            'The symplectic solver currently supports dimensions 1 and 2.');
    end
end
