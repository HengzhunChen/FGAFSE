function [T, DT, D2T] = kinetic_delta_2d(P, alpha, delta)
% KINETIC_DELTA_2D Evaluate modified 2-D kinetic data.

P2 = sum(P.^2, 2);
D2T = zeros(size(P, 1), 4);
if alpha == 2
    Pdelta = P2;
    kineticScale = ones(size(P2));
    D2T(:, 1) = 1;
    D2T(:, 4) = 1;
else
    Pdelta = P2 + delta^2;
    kineticScale = Pdelta.^((alpha - 2) / 2);
    kineticCorrection = (alpha - 2) * kineticScale ./ Pdelta;
    D2T(:, 1) = kineticScale + kineticCorrection .* P(:, 1).^2;
    D2T(:, 2) = kineticCorrection .* P(:, 1) .* P(:, 2);
    D2T(:, 3) = D2T(:, 2);
    D2T(:, 4) = kineticScale + kineticCorrection .* P(:, 2).^2;
end

T = Pdelta.^(alpha / 2) / alpha;
DT = kineticScale .* P;

end
