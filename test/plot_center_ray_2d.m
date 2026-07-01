%% Plot the 2D center ray
%
% The center ray is the single Hamiltonian trajectory starting from the
% maximum of the initial amplitude and the central phase gradient.
%
% This is a diagnostic for checking when the central trajectory approaches
% the singular momentum point P = (0, 0).

% Add functions into file path
testPath = fileparts(mfilename('fullpath'));
repoRoot = fullfile(testPath, '..');
cd(repoRoot)
FGAFSE_startup();
cd(testPath)


% ------------------------------------------------------------
% Parameters setting
% ------------------------------------------------------------
right_x = 8;
final_time = 6;
alphaList = [1.1, 1.3, 1.5, 1.7, 1.9];
alphaSelected = 1.1;
vepsExp = -8;
deltaPower = 1;
dt = 1e-3;

veps = 2^vepsExp;
delta = veps^deltaPower;
Q0 = [4, 4];  % intial position center
P0 = [0, 1];  % intial momentum center

figureFolder = fullfile(testPath, 'figures', 'center_ray_2d');
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

fprintf('2D center ray diagnostic\n');
fprintf('Q0 = [%.6f, %.6f], P0 = [%.6f, %.6f]\n', Q0, P0);
fprintf('right_x = %.6f, final_time = %.6f\n', right_x, final_time);
fprintf('veps = 2^(%d) = %.6e, delta = veps^(%.6f) = %.6e\n', ...
        vepsExp, veps, deltaPower, delta);
fprintf('center ray dt = %.6e\n\n', dt);


% ------------------------------------------------------------
% Plot one selected center ray in detail
% ------------------------------------------------------------
[time, Q, P, diagnostics] = evolve_center_ray( ...
    alphaSelected, delta, final_time, dt, Q0, P0);
Pnorm = sqrt(sum(P.^2, 2));
[minPnorm, minIndex] = min(Pnorm);

fprintf('alpha = %.6f: min |P_center(t)| = %.6e at t = %.6f\n', ...
        alphaSelected, minPnorm, time(minIndex));

figure;
subplot(2, 2, 1)
plot(Q(:, 1), Q(:, 2), 'LineWidth', 1.5);
hold on
plot(Q(1, 1), Q(1, 2), 'ko', 'MarkerFaceColor', 'k');
plot(Q(end, 1), Q(end, 2), 'ks', 'MarkerFaceColor', 'w');
hold off
box on
grid on
axis equal
xlim([0, right_x])
ylim([0, right_x])
xlabel('Q_1')
ylabel('Q_2')
title('Center ray in position')

subplot(2, 2, 2)
plot(P(:, 1), P(:, 2), 'LineWidth', 1.5);
hold on
plot(P(1, 1), P(1, 2), 'ko', 'MarkerFaceColor', 'k');
plot(P(end, 1), P(end, 2), 'ks', 'MarkerFaceColor', 'w');
plot(0, 0, 'rx', 'LineWidth', 1.5, 'MarkerSize', 8);
hold off
box on
grid on
axis equal
xlabel('P_1')
ylabel('P_2')
title('Center ray in momentum')

subplot(2, 2, 3)
plot(time, Q(:, 1), 'LineWidth', 1.5);
hold on
plot(time, Q(:, 2), 'LineWidth', 1.5);
hold off
box on
grid on
xlabel('t')
ylabel('Q')
legend('Q_1', 'Q_2', 'Location', 'best')
title('Position components')

subplot(2, 2, 4)
plot(time, P(:, 1), 'LineWidth', 1.5);
hold on
plot(time, P(:, 2), 'LineWidth', 1.5);
plot(time, Pnorm, 'k--', 'LineWidth', 1.5);
plot(time(minIndex), minPnorm, 'ro', 'MarkerFaceColor', 'r');
hold off
box on
grid on
xlabel('t')
ylabel('P')
legend('P_1', 'P_2', '|P|', 'min |P|', 'Location', 'best')
title('Momentum components')

sgtitle(sprintf('2D center ray, alpha = %.2f, delta = eps^{%.3g}', ...
        alphaSelected, deltaPower));
saveas(gcf, fullfile(figureFolder, 'center_ray_alpha_1p1.png'));


% ------------------------------------------------------------
% Plot FGA quantities along the selected center ray
% ------------------------------------------------------------
figure;
subplot(2, 2, 1)
semilogy(time, Pnorm, 'LineWidth', 1.5);
hold on
yline(delta, 'k--', '\delta');
plot(time(minIndex), minPnorm, 'ro', 'MarkerFaceColor', 'r');
hold off
box on
grid on
xlabel('t')
ylabel('|P|')
title('Distance to zero momentum')

subplot(2, 2, 2)
semilogy(time, diagnostics.D2Tnorm, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('|D^2T^\delta(P)|_F')
title('Kinetic Hessian norm')

subplot(2, 2, 3)
semilogy(time, diagnostics.detZabs, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('|det Z|')
title('Variational determinant')

subplot(2, 2, 4)
semilogy(time, diagnostics.ampAbs, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('|a|')
title('FGA amplitude')

sgtitle(sprintf('2D center-ray FGA quantities, alpha = %.2f', alphaSelected));
saveas(gcf, fullfile(figureFolder, 'center_ray_fga_quantities_alpha_1p1.png'));


% ------------------------------------------------------------
% Compare momentum norm for all alpha values
% ------------------------------------------------------------
figure;
hold on
box on
grid on
for alpha = alphaList
    [time, ~, P] = evolve_center_ray(alpha, delta, final_time, dt, Q0, P0);
    Pnorm = sqrt(sum(P.^2, 2));
    [minPnorm, minIndex] = min(Pnorm);
    plot(time, Pnorm, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\alpha = %.1f, min %.2e at t %.3f', ...
                                alpha, minPnorm, time(minIndex)));
    fprintf('alpha = %.6f: min |P_center(t)| = %.6e at t = %.6f\n', ...
            alpha, minPnorm, time(minIndex));
end
yline(delta, 'k--', 'DisplayName', '\delta');
xlabel('t')
ylabel('|P_{center}(t)|')
legend('Location', 'best')
title(sprintf('2D center ray momentum norm, delta = eps^{%.3g}', deltaPower));
hold off
saveas(gcf, fullfile(figureFolder, 'center_ray_momentum_norm.png'));


% ------------------------------------------------------------
% Local functions
% ------------------------------------------------------------

function [time, Q, P, diagnostics] = evolve_center_ray(alpha, delta, finalTime, dt, Q0, P0)
    numFullSteps = floor(finalTime / dt);
    remainingTime = finalTime - numFullSteps * dt;
    timeTolerance = 1e-12;
    if remainingTime <= timeTolerance
        remainingTime = 0;
    end
    numSteps = numFullSteps + (remainingTime > 0);

    time = zeros(numSteps + 1, 1);
    Q = zeros(numSteps + 1, 2);
    P = zeros(numSteps + 1, 2);
    DzQ = zeros(numSteps + 1, 4);
    DzP = zeros(numSteps + 1, 4);
    logAmp = zeros(numSteps + 1, 1);
    diagnostics.D2Tnorm = zeros(numSteps + 1, 1);
    diagnostics.detZabs = zeros(numSteps + 1, 1);
    diagnostics.ampAbs = zeros(numSteps + 1, 1);

    Q(1, :) = Q0;
    P(1, :) = P0;
    DzQ(1, :) = [1, 0, 0, 1];
    DzP(1, :) = -1i * [1, 0, 0, 1];
    logAmp(1) = log(2);
    diagnostics = update_diagnostics( ...
        diagnostics, 1, P(1, :), DzQ(1, :), DzP(1, :), logAmp(1), ...
        alpha, delta);

    for step = 1 : numSteps
        stepSize = dt;
        if step > numFullSteps
            stepSize = remainingTime;
        end

        [Q(step + 1, :), P(step + 1, :), DzQ(step + 1, :), ...
            DzP(step + 1, :), logAmp(step + 1)] ...
        = rk4_center_ray_step( ...
            Q(step, :), P(step, :), DzQ(step, :), DzP(step, :), ...
            logAmp(step), alpha, delta, stepSize);
        time(step + 1) = time(step) + stepSize;
        diagnostics = update_diagnostics( ...
            diagnostics, step + 1, P(step + 1, :), DzQ(step + 1, :), ...
            DzP(step + 1, :), logAmp(step + 1), alpha, delta);
    end
end

function [Qnew, Pnew, DzQnew, DzPnew, logAmpNew] = rk4_center_ray_step( ...
    Q, P, DzQ, DzP, logAmp, alpha, delta, dt)

    [kQ1, kP1, kDzQ1, kDzP1, kLogA1] = center_ray_rhs( ...
        Q, P, DzQ, DzP, alpha, delta);
    [kQ2, kP2, kDzQ2, kDzP2, kLogA2] = center_ray_rhs( ...
        Q + kQ1 * dt / 2, P + kP1 * dt / 2, ...
        DzQ + kDzQ1 * dt / 2, DzP + kDzP1 * dt / 2, alpha, delta);
    [kQ3, kP3, kDzQ3, kDzP3, kLogA3] = center_ray_rhs( ...
        Q + kQ2 * dt / 2, P + kP2 * dt / 2, ...
        DzQ + kDzQ2 * dt / 2, DzP + kDzP2 * dt / 2, alpha, delta);
    [kQ4, kP4, kDzQ4, kDzP4, kLogA4] = center_ray_rhs( ...
        Q + kQ3 * dt, P + kP3 * dt, ...
        DzQ + kDzQ3 * dt, DzP + kDzP3 * dt, alpha, delta);

    Qnew = Q + (kQ1 + 2 * kQ2 + 2 * kQ3 + kQ4) * dt / 6;
    Pnew = P + (kP1 + 2 * kP2 + 2 * kP3 + kP4) * dt / 6;
    DzQnew = DzQ + (kDzQ1 + 2 * kDzQ2 + 2 * kDzQ3 + kDzQ4) * dt / 6;
    DzPnew = DzP + (kDzP1 + 2 * kDzP2 + 2 * kDzP3 + kDzP4) * dt / 6;
    logAmpNew = logAmp + (kLogA1 + 2 * kLogA2 + 2 * kLogA3 + kLogA4) * dt / 6;
end

function [DQ, DP, DtDzQ, DtDzP, DlogAmp] = center_ray_rhs( ...
    Q, P, DzQ, DzP, alpha, delta)

    [~, DQ, D2T] = kinetic_delta_2d(P, alpha, delta);
    [~, DV, D2V] = potential(Q(:, 1), Q(:, 2));
    DP = -DV;
    DtDzQ = multiply_2x2(DzP, D2T);
    DtDzP = -multiply_2x2(DzQ, D2V);

    Z = DzQ + 1i * DzP;
    DZ = DtDzQ + 1i * DtDzP;
    detZ = Z(:, 1) .* Z(:, 4) - Z(:, 2) .* Z(:, 3);
    DlogAmp = 0.5 * trace_inverse_product_2x2(Z, DZ, detZ);
end

function diagnostics = update_diagnostics( ...
    diagnostics, index, P, DzQ, DzP, logAmp, alpha, delta)

    [~, ~, D2T] = kinetic_delta_2d(P, alpha, delta);
    Z = DzQ + 1i * DzP;
    detZ = Z(:, 1) .* Z(:, 4) - Z(:, 2) .* Z(:, 3);
    diagnostics.D2Tnorm(index) = sqrt(sum(abs(D2T).^2, 2));
    diagnostics.detZabs(index) = abs(detZ);
    diagnostics.ampAbs(index) = abs(exp(logAmp));
end

function C = multiply_2x2(A, B)
    C = zeros(size(A));
    C(:, 1) = A(:, 1) .* B(:, 1) + A(:, 2) .* B(:, 3);
    C(:, 2) = A(:, 1) .* B(:, 2) + A(:, 2) .* B(:, 4);
    C(:, 3) = A(:, 3) .* B(:, 1) + A(:, 4) .* B(:, 3);
    C(:, 4) = A(:, 3) .* B(:, 2) + A(:, 4) .* B(:, 4);
end

function value = trace_inverse_product_2x2(A, B, detA)
    value = ( A(:, 4) .* B(:, 1) - A(:, 2) .* B(:, 3) ...
        - A(:, 3) .* B(:, 2) + A(:, 1) .* B(:, 4) ) ./ detA;
end

function u0 = initWave(X, Y, veps)
    r = sqrt( (X - 4).^2 + (Y - 4).^2 );
    u0 = sqrt(128 / pi) * exp(-64 * r.^2) .* exp(1i * (Y - 1) / veps);
end

function [V, DV, D2V] = potential(Q1, Q2)
    a = 0.1;
    b1 = 4;
    b2 = 4;
    V = a * ((Q1 - b1).^2 + (Q2 - b2).^2);
    DV = 2 * a * [Q1 - b1, Q2 - b2];
    D2V = repmat(2 * a * [1, 0, 0, 1], size(Q1));
end
