function [T, DT, D2T] = kinetic_delta_1d(P, alpha, delta)
% KINETIC_DELTA_1D Evaluate modified 1-D kinetic data.

P2 = P.^2;
if alpha == 2
    Pdelta = P2;
    kineticScale = ones(size(P));
    D2T = ones(size(P));
else
    Pdelta = P2 + delta^2;
    kineticScale = Pdelta.^((alpha - 2) / 2);
    D2T = kineticScale .* (1 + (alpha - 2) * P2 ./ Pdelta);
end

T = Pdelta.^(alpha / 2) / alpha;
DT = kineticScale .* P;

end
