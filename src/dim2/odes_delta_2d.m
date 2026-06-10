function [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = odes_delta_2d(Q, P, DzQ, DzP, alpha, delta, potential)
% ODES_DELTA_2D Compute derivative values of modified FGA odes at Q, P, DzQ, DzP.
%     Inputs:
%         Q, P, DzQ, DzP -- Gaussian parameters in FGA
%         alpha          -- order of fractional operator
%         delta          -- regularization parameter
%         potential      -- function handle of potential and its derivatives
%                           [V, DV, D2V] = potential(Q1, Q2)
%     Outputs: 
%         DQ, DP, DS, DlogA, DtDzQ, DtDzP
%                        -- Derivative of FGA odes at Q, P, DzQ, DzP
%     Note: DS, DlogA are scalars in shape of (nGB, 1),
%           DQ, DP are matrices in shape of (nGB, 2),
%           DzQ, DzP, DtDzQ, DtDzP are matrices in shape of (nGB, 4).
%
%     See also time_evolution, FGA2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


P2 = sum(P.^2, 2);
D2T = zeros(size(DzP));
if alpha == 2
    % Recover the unregularized quadratic kinetic energy exactly.
    Pdelta = P2;
    kineticScale = ones(size(P2));
    D2T(:, 1) = 1;
    D2T(:, 4) = 1;
else
    Pdelta = P2 + delta^2;
    kineticScale = Pdelta.^((alpha-2)/2);
    kineticCorrection = (alpha - 2) * kineticScale ./ Pdelta;
    D2T(:, 1) = kineticScale + kineticCorrection .* P(:, 1).^2;
    D2T(:, 2) = kineticCorrection .* P(:, 1) .* P(:, 2);
    D2T(:, 3) = D2T(:, 2);
    D2T(:, 4) = kineticScale + kineticCorrection .* P(:, 2).^2;
end

[V, DV, D2V] = potential(Q(:, 1), Q(:, 2));

T = Pdelta.^(alpha/2) / alpha;
DQ = kineticScale .* P;
DP = -DV;
DS = sum(P .* DQ, 2) - T - V;

DtDzQ = multiply_2x2(DzP, D2T);
DtDzP = -multiply_2x2(DzQ, D2V);

Z = DzQ + 1i * DzP;
detZ = Z(:, 1) .* Z(:, 4) - Z(:, 2) .* Z(:, 3);
zTolerance = 1e-12;
if any(isnan(detZ) | isinf(detZ), 'all') || any(abs(detZ) < zTolerance, 'all')
    error('FGAFSE:SingularZ', 'Numerical loss of invertibility in Z.');
end

DZ = DtDzQ + 1i * DtDzP;
DlogA = 0.5 * trace_inverse_product_2x2(Z, DZ, detZ);

end

% -------------------------------------------------------------------
% Helper functions
% -------------------------------------------------------------------

function C = multiply_2x2(A, B)
    % Multiply batches of flattened row-major 2-by-2 matrices.

    C = zeros(size(A));
    C(:, 1) = A(:, 1) .* B(:, 1) + A(:, 2) .* B(:, 3);
    C(:, 2) = A(:, 1) .* B(:, 2) + A(:, 2) .* B(:, 4);
    C(:, 3) = A(:, 3) .* B(:, 1) + A(:, 4) .* B(:, 3);
    C(:, 4) = A(:, 3) .* B(:, 2) + A(:, 4) .* B(:, 4);
end

function value = trace_inverse_product_2x2(A, B, detA)
    % Compute trace(A^{-1} B) for flattened row-major 2-by-2 matrices.

    value = ( A(:, 4) .* B(:, 1) - A(:, 2) .* B(:, 3) ...
        - A(:, 3) .* B(:, 2) + A(:, 1) .* B(:, 4) ) ./ detA;
end
