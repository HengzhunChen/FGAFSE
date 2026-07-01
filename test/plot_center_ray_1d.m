%% Plot the 1D center ray
%
% The center ray is the single Hamiltonian trajectory starting from the
% maximum of the initial amplitude and the central phase gradient.

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
dt = 1e-4;

veps = 2^vepsExp;
delta = veps^deltaPower;
Q0 = 4;
P0 = 1;

figureFolder = fullfile(testPath, 'figures', 'center_ray_1d');
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

fprintf('1D center ray diagnostic\n');
fprintf('Q0 = %.6f, P0 = %.6f\n', Q0, P0);
fprintf('right_x = %.6f, final_time = %.6f\n', right_x, final_time);
fprintf('veps = 2^(%d) = %.6e, delta = veps^(%.6f) = %.6e\n', ...
        vepsExp, veps, deltaPower, delta);
fprintf('center ray dt = %.6e\n\n', dt);


% ------------------------------------------------------------
% Plot one selected center ray in detail
% ------------------------------------------------------------
[time, Q, P, diagnostics] = evolve_center_ray( ...
    alphaSelected, delta, final_time, dt, Q0, P0);
Pabs = abs(P);
[minPabs, minIndex] = min(Pabs);

fprintf('alpha = %.6f: min |P_center(t)| = %.6e at t = %.6f\n', ...
        alphaSelected, minPabs, time(minIndex));

figure;
subplot(1, 2, 1)
plot(time, Q, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('Q')
ylim([0, right_x])
title('Center ray position')

subplot(1, 2, 2)
plot(time, P, 'LineWidth', 1.5);
hold on
plot(time, Pabs, 'k--', 'LineWidth', 1.5);
plot(time(minIndex), minPabs, 'ro', 'MarkerFaceColor', 'r');
hold off
box on
grid on
xlabel('t')
ylabel('P')
legend('P', '|P|', 'min |P|', 'Location', 'best')
title('Center ray momentum')

sgtitle(sprintf('1D center ray, alpha = %.2f, delta = eps^{%.3g}', ...
        alphaSelected, deltaPower));
saveas(gcf, fullfile(figureFolder, 'center_ray_alpha_1p1.png'));


% ------------------------------------------------------------
% Plot FGA quantities along the selected center ray
% ------------------------------------------------------------
figure;
subplot(2, 2, 1)
semilogy(time, Pabs, 'LineWidth', 1.5);
hold on
yline(delta, 'k--', '\delta');
plot(time(minIndex), minPabs, 'ro', 'MarkerFaceColor', 'r');
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
ylabel('|D^2T^\delta(P)|')
title('Kinetic Hessian norm')

subplot(2, 2, 3)
semilogy(time, diagnostics.detZabs, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('|Z|')
title('Variational determinant')

subplot(2, 2, 4)
semilogy(time, diagnostics.ampAbs, 'LineWidth', 1.5);
box on
grid on
xlabel('t')
ylabel('|a|')
title('FGA amplitude')

sgtitle(sprintf('1D center-ray FGA quantities, alpha = %.2f', alphaSelected));
saveas(gcf, fullfile(figureFolder, 'center_ray_fga_quantities_alpha_1p1.png'));


% ------------------------------------------------------------
% Compare momentum magnitude for all alpha values
% ------------------------------------------------------------
figure;
hold on
box on
grid on
for alpha = alphaList
    [time, ~, P] = evolve_center_ray(alpha, delta, final_time, dt, Q0, P0);
    Pabs = abs(P);
    [minPabs, minIndex] = min(Pabs);
    plot(time, Pabs, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\alpha = %.1f, min %.2e at t %.3f', ...
                                alpha, minPabs, time(minIndex)));
    fprintf('alpha = %.6f: min |P_center(t)| = %.6e at t = %.6f\n', ...
            alpha, minPabs, time(minIndex));
end
yline(delta, 'k--', 'DisplayName', '\delta');
xlabel('t')
ylabel('|P_{center}(t)|')
legend('Location', 'best')
title(sprintf('1D center ray momentum magnitude, delta = eps^{%.3g}', deltaPower));
hold off
saveas(gcf, fullfile(figureFolder, 'center_ray_momentum_abs.png'));


% ------------------------------------------------------------
% Local functions
% ------------------------------------------------------------

function [time, Q, P, diagnostics] = evolve_center_ray( ...
    alpha, delta, finalTime, dt, Q0, P0)

    numFullSteps = floor(finalTime / dt);
    remainingTime = finalTime - numFullSteps * dt;
    timeTolerance = 1e-12;
    if remainingTime <= timeTolerance
        remainingTime = 0;
    end
    numSteps = numFullSteps + (remainingTime > 0);

    time = zeros(numSteps + 1, 1);
    Q = zeros(numSteps + 1, 1);
    P = zeros(numSteps + 1, 1);
    DzQ = zeros(numSteps + 1, 1);
    DzP = zeros(numSteps + 1, 1);
    logAmp = zeros(numSteps + 1, 1);
    diagnostics.D2Tnorm = zeros(numSteps + 1, 1);
    diagnostics.detZabs = zeros(numSteps + 1, 1);
    diagnostics.ampAbs = zeros(numSteps + 1, 1);

    Q(1) = Q0;
    P(1) = P0;
    DzQ(1) = 1;
    DzP(1) = -1i;
    logAmp(1) = log(sqrt(2));
    diagnostics = update_diagnostics( ...
        diagnostics, 1, P(1), DzQ(1), DzP(1), logAmp(1), ...
        alpha, delta);

    for step = 1 : numSteps
        stepSize = dt;
        if step > numFullSteps
            stepSize = remainingTime;
        end

        [Q(step + 1), P(step + 1), DzQ(step + 1), ...
            DzP(step + 1), logAmp(step + 1)] = rk4_center_ray_step( ...
            Q(step), P(step), DzQ(step), DzP(step), logAmp(step), ...
            alpha, delta, stepSize);
        time(step + 1) = time(step) + stepSize;
        diagnostics = update_diagnostics( ...
            diagnostics, step + 1, P(step + 1), DzQ(step + 1), ...
            DzP(step + 1), logAmp(step + 1), alpha, delta);
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

    [~, DQ, D2T] = kinetic_delta_1d(P, alpha, delta);
    [~, DV, D2V] = potential(Q);
    DP = -DV;
    DtDzQ = DzP .* D2T;
    DtDzP = -DzQ .* D2V;

    Z = DzQ + 1i * DzP;
    DZ = DtDzQ + 1i * DtDzP;
    DlogAmp = 0.5 * DZ ./ Z;
end

function diagnostics = update_diagnostics( ...
    diagnostics, index, P, DzQ, DzP, logAmp, alpha, delta)

    [~, ~, D2T] = kinetic_delta_1d(P, alpha, delta);
    Z = DzQ + 1i * DzP;
    diagnostics.D2Tnorm(index) = abs(D2T);
    diagnostics.detZabs(index) = abs(Z);
    diagnostics.ampAbs(index) = abs(exp(logAmp));
end

function u0 = initWave(x, veps)
    beta = 1;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - 3).^2) ...
        .* exp(1i / veps * beta * x);
end

function [V, DV, D2V] = potential(Q)
    a = 0.1;
    b = 4;
    V = a * (Q - b).^2;
    DV = 2 * a * (Q - b);
    D2V = 2 * a;
end
