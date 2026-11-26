function aggregate = plot_error_decay_2d_saved(resultFiles, alpha, vepsExp, figName)
%PLOT_ERROR_DECAY_2D_SAVED Plot a study from independent alpha result files.

if numel(resultFiles) ~= numel(alpha)
    error('FGAFSE:ResultCountMismatch', ...
        'There must be one result file for each alpha value.');
end
alpha = alpha(:).';
vepsExp = vepsExp(:).';
nalpha = numel(alpha);
nveps = numel(vepsExp);
err_L2 = zeros(nalpha, nveps);
time_FGA = zeros(nalpha, nveps);
time_TSSA = zeros(nalpha, nveps);

for i = 1:nalpha
    if ~exist(resultFiles{i}, 'file')
        error('FGAFSE:MissingResult', 'Missing result file: %s', resultFiles{i});
    end
    saved = load(resultFiles{i});
    if ~isequal(saved.vepsExp(:).', vepsExp) || ...
            abs(saved.alpha - alpha(i)) > 10 * eps(max(1, abs(alpha(i))))
        error('FGAFSE:ConfigurationMismatch', ...
            'Configuration mismatch in %s.', resultFiles{i});
    end
    if numel(saved.err_L2) ~= nveps || numel(saved.time_FGA) ~= nveps ...
            || numel(saved.time_TSSA) ~= nveps
        error('FGAFSE:ResultShapeMismatch', ...
            'Unexpected result size in %s.', resultFiles{i});
    end
    err_L2(i, :) = saved.err_L2(:).';
    time_FGA(i, :) = saved.time_FGA(:).';
    time_TSSA(i, :) = saved.time_TSSA(:).';
end

[alpha, order] = sort(alpha);
err_L2 = err_L2(order, :);
time_FGA = time_FGA(order, :);
time_TSSA = time_TSSA(order, :);

slopes = zeros(nalpha, 1);
intercepts = zeros(nalpha, 1);
styles = {'-*', '-o', '-^', '-s', '-d'};
fig = figure('Visible', 'off');
hold on;
box on;
grid on;
for i = 1:nalpha
    fit = polyfit(vepsExp, log2(err_L2(i, :)), 1);
    slopes(i) = fit(1);
    intercepts(i) = fit(2);
    plot(-vepsExp, log2(err_L2(i, :)), ...
        styles{1 + mod(i - 1, numel(styles))}, 'LineWidth', 1);
end
xlabel('$-\log_2(\varepsilon)$', 'FontSize', 15, 'Interpreter', 'latex');
ylabel('$\log_2$($L^2$ error)', 'FontSize', 15, 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'Box', 'on', 'TickDir', 'out');
legStr = cell(nalpha, 1);
for i = 1:nalpha
    legStr{i} = sprintf('$\\alpha=%g$, slope: %g', alpha(i), slopes(i));
end
legend(legStr, 'Interpreter', 'latex', 'FontSize', 14, ...
    'Box', 'on', 'Location', 'best');
hold off;

if isempty(figName)
    figName = fullfile(pwd, 'L2_err_veps.png');
end
[figFolder, figBase, figExt] = fileparts(figName);
if isempty(figFolder)
    figFolder = '.';
end
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end
if strcmpi(figExt, '.eps')
    saveas(fig, figName, 'epsc');
    saveas(fig, fullfile(figFolder, [figBase, '.png']), 'png');
else
    if isempty(figExt)
        figName = fullfile(figFolder, [figBase, '.png']);
    end
    saveas(fig, figName, 'png');
end
close(fig);

aggregate = struct('alpha', alpha, 'vepsExp', vepsExp, ...
    'err_L2', err_L2, 'time_FGA', time_FGA, 'time_TSSA', time_TSSA, ...
    'slope', slopes, 'intercept', intercepts);
end
