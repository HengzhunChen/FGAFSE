function [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = odes_delta_1d(Q, P, DzQ, DzP, alpha, delta, potential)
% ODES_DELTA_1D Compute derivative values of modified FGA odes at Q, P, DzQ, DzP
%    Inputs:
%        Q, P, DzQ, DzP -- Gaussian parameters in FGA
%        alpha          -- order of fractional operator
%        delta          -- regularization parameter
%        potential      -- function handle of potential and its derivatives
%                          [V, DV, D2V] = potential(Q)
%    Outputs:
%        DQ, DP, DS, DlogA, DtDzQ, DtDzP
%                       -- Derivative of FGA odes at Q, P, DzQ, DzP
%
%    See also time_evolution, FGA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


P2 = P.^2;
if alpha == 2
    % Recover the unregularized quadratic kinetic energy exactly.
    Pdelta = P2;
    kineticScale = ones(size(P));
else
    Pdelta = P2 + delta^2;
    kineticScale = Pdelta.^((alpha-2)/2);
end

[V, DV, D2V] = potential(Q);

T = Pdelta.^(alpha/2) / alpha;
DQ = kineticScale .* P;
DP = -DV;
DS = P .* DQ - T - V;

if alpha == 2
    D2T = ones(size(P));
else
    D2T = kineticScale .* (1 + (alpha - 2) * P2 ./ Pdelta);
end

DtDzQ = DzP .* D2T;
DtDzP = -DzQ .* D2V;

Z = DzQ + 1i * DzP;
zTolerance = 1e-12;
if any(isnan(Z) | isinf(Z), 'all') || any(abs(Z) < zTolerance, 'all')
    error('FGAFSE:SingularZ', 'Numerical loss of invertibility in Z.');
end

DZ = DtDzQ + 1i * DtDzP;
DlogA = 0.5 * DZ ./ Z;

end
