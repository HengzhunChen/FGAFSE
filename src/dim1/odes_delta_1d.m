function [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = odes_delta_1d(Q, P, DzQ, DzP, alpha, potential, veps)
% ODES_DELTA_1D Compute derivative values of modified FGA odes at Q, P, DzQ, DzP
%    Inputs:
%        Q, P, DzQ, DzP -- Gaussian parameters in FGA
%        alpha          -- order of fractional operator
%        potential      -- function handle of potential and its derivatives
%                          [V, DV, D2V] = potential(Q)
%        veps           -- scaled Planck constant
%    Outputs:
%        DQ, DP, DS, DlogA, DtDzQ, DtDzP
%                       -- Derivative of FGA odes at Q, P, DzQ, DzP
%
%    See also time_evolution, FGA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


[V, DV, D2V] = potential(Q);

if alpha == 2
    DQ = P;
    DP = -DV;
    DS = 0.5 * P.^2 - V;
    
    Z = DzQ + 1i * DzP;
    eps = 1e-12;
    Z_inv = 1 ./ (Z - min( sign(abs(Z) - eps), 0 ) * eps);
    DtDzQ = DzP;
    DtDzP = -DzQ .* D2V;
    
    DZ = DtDzQ + 1i * DtDzP;
    DlogA = 0.5 * Z_inv .* DZ;
else
    delta = veps;  % default value
    % delta = veps^(6/11);  % this option comes from current theoretical analysis 

    Pdelta = abs(P).^2 + delta^2;

    DQ = ( Pdelta.^((alpha- 2)/2) ) .* P;
    DP = -DV;
    DS = (1 - 1/alpha) * Pdelta.^(alpha/2)  - V;

    Z = DzQ + 1i * DzP;
    eps = 1e-12;
    Z_inv = 1 ./ (Z - min( sign(abs(Z) - eps), 0 ) * eps);

    DtDzQ = DzP .* ( (alpha - 2) * Pdelta.^((alpha- 4)/2) .* P .* P + Pdelta.^((alpha- 2)/2) );
    DtDzP = -DzQ .* D2V;

    DZ = DtDzQ + 1i * DtDzP;
    DlogA = 0.5 * Z_inv .* DZ;
end

end