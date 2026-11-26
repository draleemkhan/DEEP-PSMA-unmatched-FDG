function FIG_paper3_figures_v1_1
% FIG_paper3_figures_v1_1  PAPER 3 manuscript figures 1-3 (greyscale), MATLAB R2025a.
%
% Version 1.1, 25 November 2026. Changes from v1.0 (layout and export only; data and values unchanged):
%   - TeX interpreter used only for strings that need subscripts (v1.0 TeX text was rendered with
%     widened word spacing and two long Fig. 1 lines were truncated); box text drawn line by line;
%     long subscript lines split; Fig. 3 group headers drawn as bold text, tau as a Unicode character.
%   - TIFF written with print(-r1200) at the full 174 mm figure width, then converted to 8-bit
%     greyscale and saved with LZW compression at 1200 dpi (v1.0 TIFFs were RGB and tightly cropped).
%
% Version 1.0, 25 November 2026. Approved by the first author (M.A. Khan): figures in MATLAB,
% black/white/greyscale only, Fig. 1 drawn in MATLAB, TIFF at 1200 dpi.
%
% Fig. 1  Study workflow with volume accounting (values reported in manuscript draft 2.3).
% Fig. 2  Per-patient UFVF_ref with identification interval, read from the primary patient table
%         06_output\P3_conn26_20261112_001915\P3_patients_20261112_001915.csv
%         (MD5 2f88309e1be6eff7c499c6c1d36b50ed; checked before use). DMTVF_ref in the file = UFVF_ref.
% Fig. 3  Dependence on analysis choice (values reported in manuscript draft 2.3, Tables 2, 3, 5
%         and Online Resource 1 Table S5; each value commented with its source).
%
% Output: new folder 06_output\FIG_<yyyyMMdd_HHmmss>\ with, per figure,
%   FigN.tif  1200 dpi, greyscale (EJNMMI: line art >= 1200 dpi)
%   FigN.eps  vector, fonts embedded as text (EJNMMI preferred vector format)
%   FigN.pdf  vector (for checking)
% and FIG_manifest.json (script, input and output MD5s, MATLAB version, self-test results).
% Figure width 174 mm (EJNMMI single-column text area); Arial; no titles inside the artwork.
% Nothing existing is modified.

root   = fileparts(fileparts(mfilename('fullpath')));        % ...\PAPER3_Discordance
inCSV  = fullfile(root, '06_output', 'P3_conn26_20261112_001915', 'P3_patients_20261112_001915.csv');
expMD5 = '2f88309e1be6eff7c499c6c1d36b50ed';

% ---------- input check ----------
gotMD5 = md5file(inCSV);
assert(strcmp(gotMD5, expMD5), 'Input MD5 mismatch: %s (expected %s)', gotMD5, expMD5);

T = readtable(inCSV, 'TextType', 'string');
tests = selfTests(T);
fprintf('Self-tests passed: %d of %d\n', sum([tests.pass]), numel(tests));
assert(all([tests.pass]), 'Self-test failure; no figures written.');

% ---------- output folder ----------
stamp  = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outDir = fullfile(root, '06_output', ['FIG_' stamp]);
assert(~exist(outDir, 'dir'), 'Output folder already exists.');
mkdir(outDir);

old = setDefaults();
cleanup = onCleanup(@() restoreDefaults(old)); %#ok<NASGU>

f1 = makeFig1();   out1 = exportFig(f1, outDir, 'Fig1'); close(f1);
f2 = makeFig2(T);  out2 = exportFig(f2, outDir, 'Fig2'); close(f2);
f3 = makeFig3();   out3 = exportFig(f3, outDir, 'Fig3'); close(f3);

% ---------- manifest ----------
M.script         = [mfilename '.m'];
M.script_md5     = md5file([mfilename('fullpath') '.m']);
M.matlab_version = version;
M.run_time       = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
M.input          = struct('file', inCSV, 'md5', gotMD5);
M.self_tests     = tests;
M.outputs        = [out1, out2, out3];
fid = fopen(fullfile(outDir, 'FIG_manifest.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(M, 'PrettyPrint', true)); fclose(fid);
fprintf('Figures written to %s\n', outDir);
end

% =====================================================================================
function tests = selfTests(T)
% Reproduce the reported primary results from the patient table before plotting.
u  = 100 * T.DMTVF_ref;  up = 100 * T.DMTVF_upper;
tests = struct('name', {}, 'value', {}, 'expected', {}, 'pass', {});
add = @(n, v, e) struct('name', n, 'value', v, 'expected', e, 'pass', isequal(v, e));
tests(end+1) = add('n_patients',                 height(T),                       100);
tests(end+1) = add('median_UFVF_pct_2dp',        round(median(u), 2),             0.13);
tests(end+1) = add('median_upper_pct_2dp',       round(median(up), 2),            0.27);
tests(end+1) = add('patients_UFVF_gt0',          sum(u > 0),                      58);
tests(end+1) = add('patients_UFVF_gt1pct',       sum(u > 1),                      23);
tests(end+1) = add('patients_UFVF_gt5pct',       sum(u > 5),                      2);
tests(end+1) = add('patients_UFVF_gt10pct',      sum(u > 10),                     1);
tests(end+1) = add('cohort_VU_mL_1dp',           round(sum(T.V_U_mL), 1),         144.9);
tests(end+1) = add('cohort_VFDG_mL_1dp',         round(sum(T.V_FDG_total_mL), 1), 23917.2);
tests(end+1) = add('cohort_VI_registration_1dp', round(sum(T.V_I_registration_mL), 1), 457.8);
end

% =====================================================================================
function f = makeFig1()
f = newFig(160);
ax = axes(f, 'Units', 'normalized', 'Position', [0 0 1 1]);
hold(ax, 'on'); axis(ax, 'off'); xlim(ax, [0 174]); ylim(ax, [0 160]);
G1 = 0.93; G2 = 0.80;                                    % light and mid grey fills
cx = 58;
dot = char(183); tms = char(215); nd = char(8211);

box(ax, 8, 145, 100, 13, {'DEEP-PSMA training set: 100 patients', ...
    'paired PSMA and FDG PET/CT, expert TTB masks', 'supplied rigid registration'}, 1, '-', 0.6);
arrow(ax, cx, 145, cx, 140);
box(ax, 8, 120, 100, 20, {'Tumour-blind registration (CT only)', ...
    'trunk: translation-corrected T_A (97 patients)', 'head and pelvis: regional rigid refinement', ...
    'blinded visual QC: 29 failed regions in 23 patients'}, 1, '-', 0.6);
box(ax, 114, 120, 56, 20, {'Locked protocol and SAP', ['before any PSMA' nd 'FDG'], ...
    'correspondence was computed'}, G1, '--', 0.6);
arrow(ax, cx, 120, cx, 115);
box(ax, 8, 102, 100, 13, {'3,838 FDG components (26-connectivity)', '23,917.2 mL native FDG volume'}, 1, '-', 0.6);
box(ax, 114, 102, 56, 13, {'V_S: 871 components, 56.5 mL', 'below eligibility (denominator only)'}, G1, '-', 0.6);
arrow(ax, 108, 108.5, 114, 108.5);
arrow(ax, cx, 102, cx, 97);
box(ax, 8, 82, 100, 15, {'2,967 eligible components', 'mapped to the PSMA PET grid; coverage check', ...
    ['matching within ' char(964) '_i (voxel diagonal, median 6.10 mm)']}, 1, '-', 0.6);
% branch to the three outcomes
plot(ax, [cx cx], [82 78], 'k-', 'LineWidth', 0.6);
plot(ax, [30 144], [78 78], 'k-', 'LineWidth', 0.6);
for x = [30 87 144], arrow(ax, x, 78, x, 74); end
box(ax, 2, 56, 56, 18, {'V_I correspondence-indeterminate', '584.1 mL', ...
    ['registration 457.8 ' dot ' band 8.8'], 'coverage 117.5'}, G1, '--', 0.6);
box(ax, 59, 56, 56, 18, {'V_M matched', '2,632 components', '23,131.7 mL'}, 1, '-', 0.6);
box(ax, 116, 56, 56, 18, {'V_U unmatched', '217 components', '144.9 mL'}, G2, '-', 1.4);
arrow(ax, 87, 56, 87, 51);
box(ax, 22, 35, 130, 16, {'Primary endpoint per patient', 'UFVF_{ref} = V_U / V_{FDG,total}  (median 0.13%)', ...
    'upper bound (V_U + V_I) / V_{FDG,total}  (median 0.27%)', 'methodological identification interval'}, 1, '-', 0.6);
plot(ax, [87 87], [35 33.5], 'k-', 'LineWidth', 0.6);
plot(ax, [45 132], [33.5 33.5], 'k-', 'LineWidth', 0.6);
for x = [45 132], arrow(ax, x, 33.5, x, 31); end
box(ax, 2, 9, 86, 22, {'Prespecified analyses', ['matching rules ' dot ' voxel level ' dot ' connectivity ' dot ' eligibility'], ...
    ['registration-QC handling ' dot ' visceral organs'], ['geometric stress test (2 levels ' tms ' 50 perturbations)'], ...
    'PSMA uptake and anatomy of unmatched components'}, 1, '-', 0.6);
box(ax, 92, 9, 80, 22, {'Post-hoc analyses', ['stability ' tms ' PSMA uptake cross-tabulation'], ...
    'registration-strategy comparator', '(supplied registration; conditional', 'on the adopted QC decisions)'}, 1, '--', 0.6);
putText(ax, 2, 5.2, 'V_{FDG,total} = V_S + V_I + V_M + V_U in every patient (native FDG volumes)', 6.5, 'left');
putText(ax, 2, 2.2, 'Dashed outline: governance, correspondence-indeterminate or post-hoc', 6.5, 'left');
end

% =====================================================================================
function f = makeFig2(T)
f = newFig(80);
ax = axes(f); hold(ax, 'on');
d  = sortrows(T, {'DMTVF_ref', 'DMTVF_upper'});
x  = (1:height(d))';  lo = 100 * d.DMTVF_ref;  hi = 100 * d.DMTVF_upper;
reg = d.V_I_registration_mL > 0;
CAP = 15;  GL = 0.55 * [1 1 1];
yline(ax, 1, '--', 'Color', 0.45 * [1 1 1], 'LineWidth', 0.6);
yline(ax, 5, '--', 'Color', 0.45 * [1 1 1], 'LineWidth', 0.6);
text(ax, 101.5, 1, '1%', 'FontSize', 7, 'VerticalAlignment', 'middle', 'Interpreter', 'none');
text(ax, 101.5, 5, '5%', 'FontSize', 7, 'VerticalAlignment', 'middle', 'Interpreter', 'none');
for i = 1:numel(x)
    if hi(i) > lo(i) + 1e-9
        plot(ax, [x(i) x(i)], [lo(i) min(hi(i), CAP)], '-', 'Color', GL, 'LineWidth', 1.0);
        if hi(i) > CAP
            plot(ax, x(i), CAP, '^', 'MarkerSize', 4, 'MarkerFaceColor', GL, 'MarkerEdgeColor', GL);
            text(ax, x(i) + 1.2, CAP - 0.2, sprintf('upper bound %.1f%%', hi(i)), 'FontSize', 6.5, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Interpreter', 'none');
        end
    end
end
h1 = plot(ax, x(~reg), lo(~reg), 'o', 'MarkerSize', 3.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h2 = plot(ax, x(reg),  lo(reg),  'd', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
h3 = plot(ax, NaN, NaN, '-', 'Color', GL, 'LineWidth', 1.0);
lg = legend(ax, [h1 h2 h3], {'UFVF_{ref}', 'UFVF_{ref}, patient with registration-indeterminate volume', ...
    'identification interval (to upper bound)'}, 'Location', 'northoutside', 'NumColumns', 3, 'Box', 'off', 'FontSize', 7, ...
    'Interpreter', 'tex');
lg.ItemTokenSize = [12 8];
xlim(ax, [0 103]); ylim(ax, [-0.3 CAP + 0.8]);
xticks(ax, [1 25 50 75 100]); yticks(ax, 0:5:15); ax.TickLabelInterpreter = 'none';
xlabel(ax, 'Patients ordered by UFVF_{ref}', 'Interpreter', 'tex'); ylabel(ax, 'Fraction of FDG tumour volume (%)', 'Interpreter', 'none');
styleAxes(ax);
end

% =====================================================================================
function f = makeFig3()
TOT = 23917.2;                                            % cohort native FDG volume, mL
tms = char(215); deg = char(176);
% type: H header, P primary, S prespecified, R range (stress test), X post-hoc
% columns: label, type, cohort V_U (mL) [median lo hi], patients > 5% [median lo hi]
R = {
 'Primary',                          'H', [],                     []
 'Reference rule (voxel diagonal)',     'P', 144.9,                  2      % Table 2
 'Matching rule',                    'H', [],                     []
 ['Permissive (2 ' tms ' max voxel)'],  'S', 129.0,                  2      % Table 2
 'Strict (max voxel)',                  'S', 185.9,                  5      % Table 2
 ['Overlap only (' char(964) ' = 0)'],             'S', 287.7,                  11     % Table 2
 'Granularity',                      'H', [],                     []
 'Voxel level',                         'S', 968.8,                  31     % Results, voxel level
 'Partition and eligibility',        'H', [],                     []
 '6-connectivity',                      'S', 153.7,                  3      % Table S5
 '18-connectivity',                     'S', 151.3,                  3      % Table S5
 'Eligibility rule B',                  'S', 130.0,                  2      % Table S5
 'Registration-QC handling',         'H', [],                     []
 'S1 grade-1 regions recovered',        'S', 148.2,                  2      % Table S5
 'S2 non-assessable heads released',    'S', 144.9,                  2      % Table S5
 'S3b majority rule',                   'S', 154.2,                  2      % Table S5
 'S3c split components',                'S', 152.3,                  2      % Table S5
 'S4 failed necks released',            'S', 163.3,                  2      % Table S5
 'Geometric stress test',            'H', [],                     []
 ['Lower level (5.16 mm, 2.29' deg ')'],  'R', [241.8 171.4 332.8],  [7.5 3 14]    % Table 3
 ['Higher level (8.50 mm, 4.47' deg ')'], 'R', [709.7 242.8 1094.2], [30 11 44]    % Table 3
 'Registration strategy (post-hoc)', 'H', [],                     []
 'Supplied registration',               'X', 314.7,                  10     % Table 5
 };
n = size(R, 1); y = (n:-1:1)';
f = newFig(120);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axA = nexttile(tl); axB = nexttile(tl);
for ax = [axA axB], hold(ax, 'on'); end
GL = 0.55 * [1 1 1];
for i = 1:n
    t = R{i, 2};
    if t == 'H'
        if i > 1
            for ax = [axA axB], yline(ax, y(i) + 0.5, '-', 'Color', 0.82 * [1 1 1], 'LineWidth', 0.4); end
        end
        continue
    end
    va = R{i, 3} * 100 / TOT;  vb = R{i, 4};
    for k = 1:2
        if k == 1, ax = axA; v = va; else, ax = axB; v = vb; end
        if t == 'R', plot(ax, v(2:3), [y(i) y(i)], '-', 'Color', GL, 'LineWidth', 1.6); end
        switch t
            case 'P', plot(ax, v(1), y(i), 's', 'MarkerSize', 5,   'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
            case {'S', 'R'}, plot(ax, v(1), y(i), 'o', 'MarkerSize', 4.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
            case 'X', plot(ax, v(1), y(i), '^', 'MarkerSize', 5,   'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
        end
    end
end
xline(axA, 144.9 * 100 / TOT, '--', 'Color', 0.45 * [1 1 1], 'LineWidth', 0.6);
xline(axB, 2, '--', 'Color', 0.45 * [1 1 1], 'LineWidth', 0.6);
isH = strcmp(R(:, 2), 'H');
yd = y(~isH); ld = R(~isH, 1);
[yd, o] = sort(yd); ld = ld(o);
for ax = [axA axB], ylim(ax, [0.4 n + 0.6]); styleAxes(ax); ax.TickLabelInterpreter = 'none'; end
set(axA, 'YTick', yd, 'YTickLabel', ld);
set(axB, 'YTick', yd, 'YTickLabel', []);
for i = find(isH)'
    text(axA, -0.035, (y(i) - 0.4) / (n + 0.2), R{i, 1}, 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 7.5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'none', 'Clipping', 'off');
end
xlim(axA, [0 5]); xlim(axB, [0 46]);
xlabel(axA, 'Cohort unmatched FDG volume (% of total)', 'Interpreter', 'none');
xlabel(axB, 'Patients with unmatched fraction > 5%', 'Interpreter', 'none');
title(axA, 'a', 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none'); axA.TitleHorizontalAlignment = 'left';
title(axB, 'b', 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none'); axB.TitleHorizontalAlignment = 'left';
h = gobjects(1, 4);
h(1) = plot(axB, NaN, NaN, 's', 'MarkerSize', 5,   'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h(2) = plot(axB, NaN, NaN, 'o', 'MarkerSize', 4.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h(3) = plot(axB, NaN, NaN, '-', 'Color', GL, 'LineWidth', 1.6);
h(4) = plot(axB, NaN, NaN, '^', 'MarkerSize', 5,   'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
lg = legend(axB, h, {'Primary analysis', 'Prespecified analysis', 'Range over 50 perturbations (marker: median)', ...
    'Post-hoc analysis'}, 'NumColumns', 2, 'Box', 'off', 'FontSize', 7, 'Interpreter', 'none');
lg.Layout.Tile = 'south'; lg.ItemTokenSize = [12 8];
end

% =====================================================================================
function f = newFig(h_mm)
f = figure('Units', 'centimeters', 'Position', [2 2 17.4 h_mm / 10], 'Color', 'w', ...
    'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto');
end

function styleAxes(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, 'FontSize', 7.5, 'XColor', 'k', 'YColor', 'k', ...
    'TickLabelInterpreter', 'none');
end

function box(ax, x, y, w, h, str, g, ls, lw)
rectangle(ax, 'Position', [x y w h], 'Curvature', [0.04 0.12], 'FaceColor', g * [1 1 1], ...
    'EdgeColor', 'k', 'LineStyle', ls, 'LineWidth', lw);
if ischar(str), str = {str}; end
dy = 7 * 0.3528 * 1.35;                                   % line pitch in mm for 7-pt text
n = numel(str); y0 = y + h / 2 + (n - 1) / 2 * dy;
for k = 1:n
    putText(ax, x + w / 2, y0 - (k - 1) * dy, str{k}, 7, 'center');
end
end

function putText(ax, x, y, s, fs, ha)
% TeX only where needed (subscripts/superscripts/Greek), plain text otherwise.
if any(s == '_') || any(s == '^') || any(s == '\'), interp = 'tex'; else, interp = 'none'; end
text(ax, x, y, s, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle', 'FontSize', fs, ...
    'Interpreter', interp, 'Clipping', 'off');
end

function arrow(ax, x1, y1, x2, y2)
L = 1.8; Wd = 1.3;                                        % arrowhead length and width, mm
v = [x2 - x1, y2 - y1]; v = v / norm(v); p = [-v(2), v(1)];
b = [x2 y2] - L * v;
plot(ax, [x1 b(1)], [y1 b(2)], 'k-', 'LineWidth', 0.6);
patch(ax, [x2, b(1) + Wd / 2 * p(1), b(1) - Wd / 2 * p(1)], [y2, b(2) + Wd / 2 * p(2), b(2) - Wd / 2 * p(2)], ...
    'k', 'EdgeColor', 'none');
end

function out = exportFig(f, outDir, name)
drawnow;
tif = fullfile(outDir, [name '.tif']); epsf = fullfile(outDir, [name '.eps']); pdf = fullfile(outDir, [name '.pdf']);
set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
print(f, tif, '-dtiff', '-r1200');
img = imread(tif);
if size(img, 3) == 3, img = rgb2gray(img); end
imwrite(img, tif, 'tif', 'Resolution', 1200, 'Compression', 'lzw');
exportgraphics(f, epsf, 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(f, pdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
info = imfinfo(tif);
if isempty(info.XResolution), info.XResolution = NaN; end
out = struct('figure', name, ...
    'tif', [name '.tif'], 'tif_md5', md5file(tif), 'tif_px', [info.Width info.Height], ...
    'tif_dpi', info.XResolution, 'tif_colortype', info.ColorType, 'tif_width_mm', 25.4 * info.Width / info.XResolution, ...
    'eps', [name '.eps'], 'eps_md5', md5file(epsf), 'pdf', [name '.pdf'], 'pdf_md5', md5file(pdf));
fprintf('%s: %d x %d px at %g dpi (%.1f mm wide)\n', name, info.Width, info.Height, info.XResolution, out.tif_width_mm);
end

function h = md5file(p)
fid = fopen(p, 'r'); assert(fid > 0, 'Cannot open %s', p);
b = fread(fid, inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(b(:), 'int8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function old = setDefaults()
names = {'defaultAxesFontName', 'defaultTextFontName', 'defaultLegendFontName', 'defaultAxesFontSize'};
old = cellfun(@(n) get(groot, n), names, 'UniformOutput', false);
set(groot, 'defaultAxesFontName', 'Arial', 'defaultTextFontName', 'Arial', ...
    'defaultLegendFontName', 'Arial', 'defaultAxesFontSize', 7.5);
old = [names; old];
end

function restoreDefaults(old)
for k = 1:size(old, 2), set(groot, old{1, k}, old{2, k}); end
end
