function A3_intrareader_kappa_v1_0
% A3_intrareader_kappa_v1_0  Intra-reader agreement of the head/pelvis visual registration QC
% (hidden repeats VR101-VR120), MATLAB R2025a. Base MATLAB only (no toolbox functions).
%
% Version 1.0, 25 November 2026. Approved by the first author (M.A. Khan) to re-run, in MATLAB,
% the calculation first done in memory on 25 Nov 2026 (A3_intrareader_kappa_note_20261125.md).
% Follows A3_visualQC_analysis_plan_20261108.md, section 6:
%   primary   - exact and within-one-grade agreement for the A2 grade, the 3a grade and the preference,
%               separately for head (n = 3) and pelvis (n = 17);
%   secondary - pelvis only: weighted kappa for grades, unweighted kappa for preference.
% Kappa weights were not specified in the plan; linear weights are the main value (decision of
% 25 Nov 2026), quadratic weights are reported for completeness. Grade categories 0-3;
% preference categories A2 / equivalent / A3a.
%
% Input : 06_output\A3_visualQC_unblinded_20261108\A3_repeat_pairs.csv (MD5 checked)
% Output: new folder 06_output\A3_kappa_<yyyyMMdd_HHmmss>\ with
%         A3_intrareader_agreement.csv, A3_intrareader_contingency.txt, A3_kappa_manifest.json
% Nothing existing is modified.

root   = fileparts(fileparts(mfilename('fullpath')));
inCSV  = fullfile(root, '06_output', 'A3_visualQC_unblinded_20261108', 'A3_repeat_pairs.csv');
expMD5 = 'b10920ee2a62af34bd7e8795fdcd8c81';
gotMD5 = md5file(inCSV);
assert(strcmp(gotMD5, expMD5), 'Input MD5 mismatch: %s (expected %s)', gotMD5, expMD5);

% ---------- synthetic self-tests of the kappa function ----------
st = struct('name', {}, 'value', {}, 'expected', {}, 'pass', {});
st(end+1) = chk('perfect agreement, unweighted', kappaTab([5 0; 0 7], 'none'), 1);
st(end+1) = chk('textbook 2x2 [20 5; 10 15], unweighted', kappaTab([20 5; 10 15], 'none'), 0.4);
st(end+1) = chk('3x3 linear equals unweighted for 2 used categories', ...
    kappaTab([20 5 0; 10 15 0; 0 0 0], 'linear'), 0.4);
st(end+1) = chk('constant rater gives 0', kappaTab([0 1; 0 16], 'linear'), 0);

% ---------- data ----------
T = readtable(inCSV, 'TextType', 'string');
assert(height(T) == 20, 'Expected 20 repeat pairs.');
grades = 0:3;  prefs = ["A2" "equivalent" "A3a"];
R = table();  txt = strings(0, 1);
for region = ["head" "pelvis"]
    S = T(T.region == region, :);  n = height(S);
    for item = ["A2" "A3a" "pref"]
        if item == "pref"
            a = S.pref_first;  b = S.pref_repeat;  cats = prefs;
            exact = sum(a == b);  within1 = NaN;
        else
            a = S.(item + "_first");  b = S.(item + "_repeat");  cats = grades;
            exact = sum(a == b);  within1 = sum(abs(a - b) <= 1);
        end
        O = tabulateTab(a, b, cats);
        kl = NaN; kq = NaN; ku = NaN;
        if region == "pelvis"
            if item == "pref", ku = kappaTab(O, 'none');
            else, kl = kappaTab(O, 'linear'); kq = kappaTab(O, 'quadratic'); end
        end
        R = [R; table(region, item, n, exact, within1, kl, kq, ku, ...
            'VariableNames', {'region', 'measure', 'n', 'exact', 'within_one', ...
            'kappa_linear', 'kappa_quadratic', 'kappa_unweighted'})]; %#ok<AGROW>
        txt(end+1) = sprintf('%s %s (rows first read, columns repeat; categories %s):', region, item, strjoin(string(cats), '/')); %#ok<AGROW>
        txt(end+1) = mat2str(O); %#ok<AGROW>
    end
end
disp(R);

% ---------- reproduction check against the in-memory calculation of 25 Nov 2026 ----------
P = R(R.region == "pelvis", :);
st(end+1) = chk('pelvis A2 exact',            P.exact(P.measure == "A2"),   14);
st(end+1) = chk('pelvis A3a exact',           P.exact(P.measure == "A3a"),  16);
st(end+1) = chk('pelvis pref exact',          P.exact(P.measure == "pref"), 14);
st(end+1) = chk('pelvis A2 linear kappa 3dp', round(P.kappa_linear(P.measure == "A2"), 3),    0.679);
st(end+1) = chk('pelvis A2 quadratic 3dp',    round(P.kappa_quadratic(P.measure == "A2"), 3), 0.718);
st(end+1) = chk('pelvis A3a linear kappa',    round(P.kappa_linear(P.measure == "A3a"), 3),   0);
st(end+1) = chk('pelvis pref kappa 3dp',      round(P.kappa_unweighted(P.measure == "pref"), 3), 0.617);
H = R(R.region == "head", :);
st(end+1) = chk('head all exact (3 of 3)',    all(H.exact == 3), true);
fprintf('Checks passed: %d of %d\n', sum([st.pass]), numel(st));
if ~all([st.pass]), warning('Some checks failed; see manifest.'); end

% ---------- outputs ----------
outDir = fullfile(root, '06_output', ['A3_kappa_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
assert(~exist(outDir, 'dir')); mkdir(outDir);
f1 = fullfile(outDir, 'A3_intrareader_agreement.csv');   writetable(R, f1);
f2 = fullfile(outDir, 'A3_intrareader_contingency.txt'); writelines(txt, f2);
M.script = [mfilename '.m'];  M.script_md5 = md5file([mfilename('fullpath') '.m']);
M.matlab_version = version;   M.run_time = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
M.input = struct('file', inCSV, 'md5', gotMD5);
M.weights_note = 'Plan did not specify kappa weights; linear = main value (decision 25 Nov 2026); quadratic for completeness.';
M.checks = st;
M.outputs = struct('file', {f1, f2}, 'md5', {md5file(f1), md5file(f2)});
fid = fopen(fullfile(outDir, 'A3_kappa_manifest.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(M, 'PrettyPrint', true)); fclose(fid);
fprintf('Results written to %s\n', outDir);
end

% =====================================================================================
function O = tabulateTab(a, b, cats)
k = numel(cats); O = zeros(k);
for i = 1:numel(a)
    O(find(cats == a(i), 1), find(cats == b(i), 1)) = O(find(cats == a(i), 1), find(cats == b(i), 1)) + 1;
end
assert(sum(O(:)) == numel(a), 'Unrecognised category in input.');
end

function k = kappaTab(O, w)
% Cohen's kappa from a square contingency table with disagreement weights:
% none (0/1), linear |i-j|/(K-1), quadratic ((i-j)/(K-1))^2.
K = size(O, 1); n = sum(O(:)); [I, J] = ndgrid(1:K, 1:K);
switch w
    case 'none',      W = double(I ~= J);
    case 'linear',    W = abs(I - J) / (K - 1);
    case 'quadratic', W = ((I - J) / (K - 1)).^2;
end
E = sum(O, 2) * sum(O, 1) / n;
den = sum(W(:) .* E(:));
if den == 0, k = NaN; else, k = 1 - sum(W(:) .* O(:)) / den; end
end

function s = chk(name, v, e)
s = struct('name', name, 'value', v, 'expected', e, 'pass', isequal(v, e) || (isnumeric(v) && abs(v - e) < 1e-12));
end

function h = md5file(p)
fid = fopen(p, 'r'); assert(fid > 0, 'Cannot open %s', p);
b = fread(fid, inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(b(:), 'int8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
