function [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = odes_delta_1d(Q, P, DzQ, DzP, alpha, delta, potential_1d)
% ODES_DELTA_1D Compute derivative values of modified FGA odes at Q, P, DzQ, DzP
%    Inputs:
%        Q, P, DzQ, DzP -- Gaussian parameters in FGA
%        alpha          -- order of fractional operator
%        delta          -- regularization parameter
%        potential_1d   -- function handle of potential and its derivatives
%                          [V, DV, D2V] = potential_1d(Q)
%    Outputs:
%        DQ, DP, DS, DlogA, DtDzQ, DtDzP
%                       -- Derivative of FGA odes at Q, P, DzQ, DzP
%
%    See also time_evolution, FGA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


[T, DQ, D2T] = kinetic_delta_1d(P, alpha, delta);
[V, DV, D2V] = potential_1d(Q);

DP = -DV;
DS = P .* DQ - T - V;

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
