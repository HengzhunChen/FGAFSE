function [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = odes_delta_2d(Q, P, DzQ, DzP, alpha, potential, veps)
% ODES_DELTA_2D Compute derivative values of modified FGA odes at Q, P, DzQ, DzP.
%     Inputs:
%         Q, P, DzQ, DzP -- Gaussian parameters in FGA
%         alpha          -- order of fractional operator
%         potential      -- function handle of potential and its derivatives
%                           [V, DV, D2V] = potential(Q1, Q2)
%         veps           -- scaled Planck constant
%     Outputs: 
%         DQ, DP, DS, DlogA, DtDZQ, DtDzP
%                        -- Derivative of FGA odes at Q, P, DzQ, DzP
%     Note: DS, DlogA are scalars in shape of (nGB, 1),
%           DQ, DP are matrices in shape of (nGB, 2),
%           DzQ, DzP, DtDzQ, DtDzP are matrices in shape of (nGB, 4).
%
%     See also time_evolution, FGA2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


[V, DV, D2V] = potential(Q(:, 1), Q(:, 2));

if alpha == 2
    DQ = P;
    DP = -DV;
    DS = 0.5 * (sum(P.^2, 2)) - V;
    
    DtDzQ = DzP;
    DtDzP = zeros(size(DzP));
    % NOTE: DtDzP = -DzQ * D2V;
    DtDzP(:, 1) = -( DzQ(:, 1) .* D2V(:, 1) + DzQ(:, 2) .* D2V(:, 3) );
    DtDzP(:, 2) = -( DzQ(:, 1) .* D2V(:, 2) + DzQ(:, 2) .* D2V(:, 4) );
    DtDzP(:, 3) = -( DzQ(:, 3) .* D2V(:, 1) + DzQ(:, 4) .* D2V(:, 3) );
    DtDzP(:, 4) = -( DzQ(:, 3) .* D2V(:, 2) + DzQ(:, 4) .* D2V(:, 4) );
    
    Z = DzQ + 1i * DzP;
    detZ = Z(:, 1) .* Z(:, 4) - Z(:, 2) .* Z(:, 3);
    Z_inv = zeros(size(Z));
    Z_inv(:, 1) = Z(:, 4);
    Z_inv(:, 2) = -Z(:, 2);
    Z_inv(:, 3) = -Z(:, 3);
    Z_inv(:, 4) = Z(:, 1);
    eps = 1e-12;
    Z_inv = Z_inv ./ (detZ - min( sign(abs(detZ) - eps), 0 ) * eps);
    
    DZ = DtDzQ + 1i * DtDzP;
    % NOTE: DlogA = 0.5 * Z_inv * DZ;
    DlogA = 0.5 * ( Z_inv(:, 1) .* DZ(:, 1) + Z_inv(:, 2) .* DZ(:, 3) +...
                    Z_inv(:, 3) .* DZ(:, 2) + Z_inv(:, 4) .* DZ(:, 4) );
else
    delta = veps;  % default value
    % delta = veps ^ (3/5);  % this option comes from current theorectical analysis
    % delta = veps^2;

    Pdelta = sum(P.^2, 2) + delta^2;  % shape: (nGB, 1)

    DQ = Pdelta .^ ((alpha-2)/2) .* P;
    DP = -DV;
    DS = (1 - 1/alpha) * Pdelta .^ (alpha/2) - V;

    % ---------------------------------------------------------------------------------------
    % NOTE: Be careful with the derivative calculation here to match the dimension.
    %   DtDzQ 
    % = (alpha - 2) * Pdelta^((alpha-4)/2) * DzP * P * P' + Pdelta^((alpha-2)/2) * DzP
    % = Pdelta^((alpha-2)/2) * DzP * ( (alpha-2) * Pdelta^(-1) * P * P' + eye(2) )
    % = scalar * DzP * temp
    % temp = (alpha-2) * Pdelta^(-1) * P * P' + eye(2), scalar = Pdelta^((alpha-2)/2),
    % temp is a symmetric matrix.
    % ---------------------------------------------------------------------------------------
    temp = zeros(size(DzQ));
    scal = (alpha - 2) ./ Pdelta;
    temp(:, 1) = scal .* P(:, 1) .* P(:, 1) + 1.0;
    temp(:, 2) = scal .* P(:, 1) .* P(:, 2);
    temp(:, 3) = temp(:, 2);
    temp(:, 4) = scal .* P(:, 2) .* P(:, 2) + 1.0;
    DtDzQ = zeros(size(DzQ));
    scalar = Pdelta .^ ((alpha-2)/2);
    DtDzQ(:, 1) = scalar .* ( DzP(:, 1) .* temp(:, 1) + DzP(:, 2) .* temp(:, 3) );
    DtDzQ(:, 2) = scalar .* ( DzP(:, 1) .* temp(:, 2) + DzP(:, 2) .* temp(:, 4) );
    DtDzQ(:, 3) = scalar .* ( DzP(:, 3) .* temp(:, 1) + DzP(:, 4) .* temp(:, 2) );
    DtDzQ(:, 4) = scalar .* ( DzP(:, 3) .* temp(:, 2) + DzP(:, 4) .* temp(:, 4) );


    DtDzP = zeros(size(DzP));
    % NOTE: DtDzP = -DzQ * D2V;
    DtDzP(:, 1) = -( DzQ(:, 1) .* D2V(:, 1) + DzQ(:, 2) .* D2V(:, 3) );
    DtDzP(:, 2) = -( DzQ(:, 1) .* D2V(:, 2) + DzQ(:, 2) .* D2V(:, 4) );
    DtDzP(:, 3) = -( DzQ(:, 3) .* D2V(:, 1) + DzQ(:, 4) .* D2V(:, 3) );
    DtDzP(:, 4) = -( DzQ(:, 3) .* D2V(:, 2) + DzQ(:, 4) .* D2V(:, 4) );

    Z = DzQ + 1i * DzP;
    detZ = Z(:, 1) .* Z(:, 4) - Z(:, 2) .* Z(:, 3);
    Z_inv = zeros(size(Z));
    Z_inv(:, 1) = Z(:, 4);
    Z_inv(:, 2) = -Z(:, 2);
    Z_inv(:, 3) = -Z(:, 3);
    Z_inv(:, 4) = Z(:, 1);
    eps = 1e-12;
    Z_inv = Z_inv ./ (detZ - min( sign(abs(detZ) - eps), 0 ) * eps);

    DZ = DtDzQ + 1i * DtDzP;
    % NOTE: DlogA = 0.5 * trace(Z_inv * DZ);
    DlogA = 0.5 * ( Z_inv(:, 1) .* DZ(:, 1) + Z_inv(:, 2) .* DZ(:, 3) + ...
                    Z_inv(:, 3) .* DZ(:, 2) + Z_inv(:, 4) .* DZ(:, 4) );
end

end