function FIG4_anatomy_v1_0
% FIG4_anatomy_v1_0  Manuscript Fig. 4: anatomical distribution of the 217 primary unmatched FDG
% components (greyscale), MATLAB R2025a. Base MATLAB only.
%
% Version 1.0, 25 November 2026. Derived from FIG_S1_anatomy_v1_0.m (MD5 5857c2cdb0abe5b521ed19302da6cb10)
% after the first author decided to move the figure from Online Resource 1 (Fig. S1) to the main text
% (Fig. 4). Changes: output names (Fig4, FIG4_<timestamp>, FIG4_manifest.json) and a wider left margin
% so that the longest group heading is not clipped. Data, self-tests and style unchanged.
%
% FIG_S1_anatomy_v1_0 history: version 1.0, 25 November 2026. Approved by the first author (M.A. Khan): supplementary figure in
% MATLAB, black/white/greyscale, TIFF at 1200 dpi; same style and export as FIG_paper3_figures_v1_1.m.
%
% Input : 06_output\AD_20261112_125432\AD_categories_20261112_125432.csv (MD5 checked)
%         (output of AD_anatomical_distribution_v1_0.py; values reported in Online Resource 1, Table S6)
% Output: new folder 06_output\FIG4_<yyyyMMdd_HHmmss>\ with Fig4.tif (1200 dpi, 8-bit greyscale,
%         LZW), Fig4.eps and Fig4.pdf (vector), FIG4_manifest.json. Nothing existing is modified.

root   = fileparts(fileparts(mfilename('fullpath')));
inCSV  = fullfile(root, '06_output', 'AD_20261112_125432', 'AD_categories_20261112_125432.csv');
expMD5 = '93ff205be4d376dac760210015dbbaf5';
gotMD5 = md5file(inCSV);
assert(strcmp(gotMD5, expMD5), 'Input MD5 mismatch: %s (expected %s)', gotMD5, expMD5);
T = readtable(inCSV, 'TextType', 'string', 'Delimiter', ',');
T.subcategory(ismissing(T.subcategory)) = "";

% ---------- rows: group header, then (category, subcategory, display label) ----------
G = {
 'Spine', {
   'bone', 'spinal: cervical',  'Cervical'
   'bone', 'spinal: thoracic',  'Thoracic'
   'bone', 'spinal: lumbar',    'Lumbar'
   'bone', 'spinal: sacral',    'Sacral'}
 'Other axial skeleton', {
   'bone', 'pelvic bones',               'Pelvic bones'
   'bone', 'ribs and costal cartilages', 'Ribs and costal cartilages'
   'bone', 'sternum',                    'Sternum'
   'bone', 'skull',                      'Skull'}
 'Shoulder girdle and proximal limbs', {
   'bone', 'scapulae',  'Scapulae'
   'bone', 'clavicles', 'Clavicles'
   'bone', 'humeri',    'Humeri'
   'bone', 'femurs',    'Femurs'}
 'Soft tissue (location only)', {
   'soft tissue', 'supradiaphragmatic location',    'Above T12/L1 (diaphragm proxy)'
   'soft tissue', 'retroperitoneal level',          'Retroperitoneal level'
   'soft tissue', 'iliac zone (region unassigned)', 'Iliac zone'}
 'Visceral organs', {
   'visceral organ', 'liver',       'Liver'
   'visceral organ', 'colon',       'Colon'
   'visceral organ', 'brain',       'Brain'
   'visceral organ', 'small bowel', 'Small bowel'
   'visceral organ', 'spleen',      'Spleen'
   'visceral organ', 'heart',       'Heart'}
 'Site-indeterminate', {
   'site-indeterminate', 'anatomy unavailable',    'Anatomy unavailable'
   'site-indeterminate', 'anatomically ambiguous', 'Anatomically ambiguous'}
 };

% ---------- assemble and self-test ----------
rows = struct('label', {}, 'header', {}, 'mL', {}, 'n', {}, 'pts', {});
allML = 0; allN = 0;
for g = 1:size(G, 1)
    items = G{g, 2}; gML = 0; gN = 0; sub = struct('label', {}, 'header', {}, 'mL', {}, 'n', {}, 'pts', {});
    for k = 1:size(items, 1)
        r = T(T.category == items{k, 1} & T.subcategory == items{k, 2}, :);
        assert(height(r) == 1, 'Category not found: %s / %s', items{k, 1}, items{k, 2});
        sub(end+1) = struct('label', items{k, 3}, 'header', false, 'mL', r.native_mL, ...
            'n', r.components, 'pts', r.patients); %#ok<AGROW>
        gML = gML + r.native_mL; gN = gN + r.components;
    end
    rows(end+1) = struct('label', sprintf('%s (%.1f mL)', G{g, 1}, gML), 'header', true, ...
        'mL', gML, 'n', gN, 'pts', NaN); %#ok<AGROW>
    rows = [rows, sub]; %#ok<AGROW>
    allML = allML + gML; allN = allN + gN;
end
tot = T(T.category == "total", :);
st = struct('name', {}, 'value', {}, 'expected', {}, 'pass', {});
st(end+1) = chk('components, all rows',  allN, 217);
st(end+1) = chk('components, total row', tot.components, 217);
st(end+1) = chk('volume 1dp, all rows',  round(allML, 1), 144.9);
st(end+1) = chk('volume 1dp, total row', round(tot.native_mL, 1), 144.9);
st(end+1) = chk('bone 1dp (three skeletal groups)', round(rows(1).mL + rows(6).mL + rows(11).mL, 1), 115.5);
st(end+1) = chk('bone components',       rows(1).n + rows(6).n + rows(11).n, 153);
st(end+1) = chk('spine 1dp',             round(rows(1).mL, 1), 28.2);
st(end+1) = chk('other axial 1dp',       round(rows(6).mL, 1), 37.3);
st(end+1) = chk('shoulder/limbs 1dp',    round(rows(11).mL, 1), 49.9);
st(end+1) = chk('soft tissue 1dp',       round(rows(16).mL, 1), 25.0);
st(end+1) = chk('visceral 1dp',          round(rows(20).mL, 1), 3.9);
fprintf('Self-tests passed: %d of %d\n', sum([st.pass]), numel(st));
assert(all([st.pass]), 'Self-test failure; no figure written.');

% ---------- figure ----------
old = setDefaults(); cleanup = onCleanup(@() restoreDefaults(old)); %#ok<NASGU>
n = numel(rows); y = (n:-1:1)';
f = figure('Units', 'centimeters', 'Position', [2 2 17.4 13], 'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperPositionMode', 'auto');
ax = axes(f, 'Position', [0.37 0.09 0.55 0.88]); hold(ax, 'on');
GL = 0.45 * [1 1 1];
for i = 1:n
    if rows(i).header
        if i > 1, yline(ax, y(i) + 0.5, '-', 'Color', 0.82 * [1 1 1], 'LineWidth', 0.4); end
        continue
    end
    barh(ax, y(i), rows(i).mL, 0.62, 'FaceColor', GL, 'EdgeColor', 'none');
    text(ax, rows(i).mL + 0.35, y(i), sprintf('%d / %d', rows(i).n, rows(i).pts), 'FontSize', 6.5, ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none');
end
isH = [rows.header];
yd = y(~isH); ld = {rows(~isH).label};
[yd, o] = sort(yd); ld = ld(o);
set(ax, 'YTick', yd, 'YTickLabel', ld, 'TickLabelInterpreter', 'none', 'Box', 'off', 'TickDir', 'out', ...
    'LineWidth', 0.6, 'FontSize', 7.5, 'XColor', 'k', 'YColor', 'k');
ylim(ax, [0.4 n + 0.6]); xlim(ax, [0 25]);
for i = find(isH)
    text(ax, -0.02, (y(i) - 0.4) / (n + 0.2), rows(i).label, 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 7.5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end
xlabel(ax, 'Native FDG volume of unmatched components (mL)', 'Interpreter', 'none');
text(ax, 1, 0.012, 'Labels: components / patients', 'Units', 'normalized', 'FontSize', 6.5, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');

% ---------- export ----------
outDir = fullfile(root, '06_output', ['FIG4_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
assert(~exist(outDir, 'dir')); mkdir(outDir);
out = exportFig(f, outDir, 'Fig4'); close(f);
M.script = [mfilename '.m']; M.script_md5 = md5file([mfilename('fullpath') '.m']);
M.matlab_version = version; M.run_time = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
M.input = struct('file', inCSV, 'md5', gotMD5); M.self_tests = st; M.outputs = out;
fid = fopen(fullfile(outDir, 'FIG4_manifest.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(M, 'PrettyPrint', true)); fclose(fid);
fprintf('Figure written to %s\n', outDir);
end

% =====================================================================================
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
out = struct('figure', name, 'tif', [name '.tif'], 'tif_md5', md5file(tif), 'tif_px', [info.Width info.Height], ...
    'tif_dpi', info.XResolution, 'tif_colortype', info.ColorType, 'tif_width_mm', 25.4 * info.Width / info.XResolution, ...
    'eps', [name '.eps'], 'eps_md5', md5file(epsf), 'pdf', [name '.pdf'], 'pdf_md5', md5file(pdf));
fprintf('%s: %d x %d px at %g dpi (%.1f mm wide, %s)\n', name, info.Width, info.Height, info.XResolution, ...
    out.tif_width_mm, info.ColorType);
end

function s = chk(name, v, e)
s = struct('name', name, 'value', v, 'expected', e, 'pass', isequal(v, e) || abs(v - e) < 1e-9);
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
