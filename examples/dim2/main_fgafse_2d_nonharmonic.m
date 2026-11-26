function aggregate = main_fgafse_2d_nonharmonic(mode, deltaPower, alphaValue)
%MAIN_FGAFSE_2D_NONHARMONIC Run or plot one non-harmonic 2D study.
%
% Examples:
%   main_fgafse_2d_nonharmonic('compute', 1, 1.3)
%   main_fgafse_2d_nonharmonic('compute', 7/12, 1.8)
%   main_fgafse_2d_nonharmonic('plot', 1)
%   main_fgafse_2d_nonharmonic('plot', 7/12)
%
% Each compute call saves its own data and error-decay plot after the
% epsilon sweep finishes. Plot mode combines all alpha values for one power.

if nargin < 2
    error('FGAFSE:Usage', ...
        ['Use main_fgafse_2d_nonharmonic(''compute'', deltaPower, alpha) ', ...
         'or main_fgafse_2d_nonharmonic(''plot'', deltaPower).']);
end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('FGAFSE:InvalidMode', 'mode must be ''compute'' or ''plot''.');
end
mode = lower(char(mode));
if ~ismember(mode, {'compute', 'plot'})
    error('FGAFSE:InvalidMode', 'mode must be ''compute'' or ''plot''.');
end
if ~isscalar(deltaPower) || ~isnumeric(deltaPower) || ~isfinite(deltaPower)
    error('FGAFSE:InvalidDeltaPower', 'deltaPower must be a finite scalar.');
end

exampleDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
cleanup = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, '..', '..'));
FGAFSE_startup();
cd(exampleDir);
helperDir = fullfile(exampleDir, 'nonharmonic_helpers');
if ~exist(helperDir, 'dir')
    error('FGAFSE:MissingHelpers', ...
        'Missing nonharmonic helper directory: %s', helperDir);
end
addpath(helperDir);

vepsExp = [-5, -6, -7, -8];
right_x = 8;
final_time = 5;
useCachedTSSA = true;
fgaDt = 1e-3;
quadraticScale = 1/16;
quarticScale = 1/128; % Adjust this to control the anharmonicity.
potentialName = sprintf(['V=%.16g*|y|^2 + %.16g*|y|^4, ' ...
    'y=(x1-4,x2-4)'], quadraticScale, quarticScale);
potential = @(Q1, Q2) nonharmonicPotential(Q1, Q2, ...
    quadraticScale, quarticScale);

study = selectStudy(deltaPower);
timeTag = sprintf('T_%g', final_time);
potentialTag = sprintf('quadratic_%g_quartic_%g', ...
    quadraticScale, quarticScale);
resultRoot = fullfile(exampleDir, 'data', 'nonharmonic_2d', timeTag, ...
    potentialTag, study.tag);
figureRoot = fullfile(exampleDir, 'figures', 'nonharmonic_2d', timeTag, ...
    potentialTag, study.tag);
tssaFolder = fullfile(resultRoot, 'tssa');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end
if ~exist(tssaFolder, 'dir')
    mkdir(tssaFolder);
end

switch mode
    case 'compute'
        if nargin < 3 || ~isscalar(alphaValue) || ~isnumeric(alphaValue) ...
                || ~isfinite(alphaValue)
            error('FGAFSE:InvalidAlpha', ...
                'compute mode requires one finite scalar alpha.');
        end

        result = error_decay_2d_single(alphaValue, vepsExp, right_x, ...
            final_time, @initWave, potential, useCachedTSSA, ...
            study.power, fgaDt, tssaFolder);
        result.deltaPower = study.power;
        result.quadraticScale = quadraticScale;
        result.quarticScale = quarticScale;
        result.potentialName = potentialName;
        resultFile = fullfile(resultRoot, alphaFileName(alphaValue));
        temporaryFile = [tempname(resultRoot), '.mat'];
        save(temporaryFile, '-struct', 'result');
        movefile(temporaryFile, resultFile, 'f');
        fprintf('Saved result: %s\n', resultFile);

        % Plot this alpha immediately; no other alpha jobs are needed.
        [~, alphaTag] = fileparts(resultFile);
        figureName = fullfile(figureRoot, [alphaTag, '.eps']);
        plot_error_decay_2d_saved({resultFile}, alphaValue, vepsExp, figureName);
        fprintf('Saved per-alpha plots: %s (.eps and .png)\n', ...
            fullfile(figureRoot, alphaTag));
        aggregate = result;

    case 'plot'
        if nargin >= 3 && ~isempty(alphaValue)
            error('FGAFSE:InvalidPlotArguments', ...
                'plot mode takes only mode and deltaPower.');
        end
        resultFiles = cell(numel(study.alpha), 1);
        for i = 1:numel(study.alpha)
            resultFiles{i} = fullfile(resultRoot, alphaFileName(study.alpha(i)));
            legacyFile = fullfile(resultRoot, legacyAlphaFileName(study.alpha(i)));
            if ~exist(resultFiles{i}, 'file') && exist(legacyFile, 'file')
                resultFiles{i} = legacyFile;
            end
        end
        figureName = fullfile(figureRoot, study.figure);
        aggregate = plot_error_decay_2d_saved(resultFiles, study.alpha, ...
            vepsExp, figureName);
        aggregate.right_x = right_x;
        aggregate.final_time = final_time;
        aggregate.deltaPower = study.power;
        aggregate.fgaDt = fgaDt;
        aggregate.quadraticScale = quadraticScale;
        aggregate.quarticScale = quarticScale;
        aggregate.potentialName = potentialName;
        save(fullfile(resultRoot, 'aggregate.mat'), '-struct', 'aggregate');
        save(fullfile(resultRoot, ['L2_error_nonharmonic_2d_' study.tag '.mat']), ...
            '-struct', 'aggregate');
        fprintf('Saved aggregate: %s\n', fullfile(resultRoot, 'aggregate.mat'));
end
end

function study = selectStudy(deltaPower)
    if abs(deltaPower - 1) <= 10 * eps(max(1, abs(deltaPower)))
        study.power = 1;
        study.tag = 'delta_eps';
        study.alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
        study.figure = 'L2err_nonharmonic_2d_delta_eps.eps';
    elseif abs(deltaPower - 7/12) <= 10 * eps(max(1, abs(deltaPower)))
        study.power = 7/12;
        study.tag = 'delta_eps_7div12';
        study.alpha = [1.75, 1.8, 1.85, 1.9, 1.95];
        study.figure = 'L2err_nonharmonic_2d_delta_eps_7div12.eps';
    else
        error('FGAFSE:UnknownDeltaPower', 'deltaPower must be 1 or 7/12.');
    end
end

function fileName = alphaFileName(alpha)
    % Avoid exposing binary floating-point tails in human-readable names.
    fileName = sprintf('alpha_%.12g.mat', alpha);
end

function fileName = legacyAlphaFileName(alpha)
    fileName = sprintf('alpha_%0.17g.mat', alpha);
end


function u0 = initWave(X, Y, veps)
    center = [4, 4];
    initialMomentum = [0, 1];
    r2 = (X - center(1)).^2 + (Y - center(2)).^2;
    u0 = sqrt(128 / pi) * exp(-64 * r2) ...
        .* exp(1i / veps * (initialMomentum(1) * (X - center(1)) ...
        + initialMomentum(2) * (Y - center(2))));
end

function [V, DV, D2V] = nonharmonicPotential(Q1, Q2, quadraticScale, quarticScale)
    y1 = Q1 - 4;
    y2 = Q2 - 4;
    radiusSquared = y1.^2 + y2.^2;
    V = quadraticScale * radiusSquared + quarticScale * radiusSquared.^2;
    DV = [(2 * quadraticScale + 4 * quarticScale * radiusSquared) .* y1, ...
          (2 * quadraticScale + 4 * quarticScale * radiusSquared) .* y2];
    h11 = 2 * quadraticScale + 4 * quarticScale * radiusSquared ...
        + 8 * quarticScale * y1.^2;
    h22 = 2 * quadraticScale + 4 * quarticScale * radiusSquared ...
        + 8 * quarticScale * y2.^2;
    h12 = 8 * quarticScale * y1 .* y2;
    D2V = [h11, h12, h12, h22];
end
