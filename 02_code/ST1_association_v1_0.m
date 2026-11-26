function ST1_association_v1_0(mode)
% ST1_association_v1_0 - PAPER 3: ST1A (patient level) and ST1B (lesion level) associations between organ displacement
% and correspondence status (MATLAB R2025a). Exploratory, associational; no hypothesis tests, no p-values.
%
% Implements: soft-tissue specification LOCKED (MD5 055199c4fd5d46039dcc65b6ab4f2416) SP4, SP5, SP7; dated implementation
% clarification 00_protocol\StepSAP_ST1_implementation_clarification_20261112.md (MD5 1021a7266020138fc7d0e13b2b4af6a8;
% author decisions A1-A7 and two safeguards).
% Organ displacement is an anatomical measurement; it does not establish that a specific component was mismatched.
%
% ST1A: Spearman rho (corr, 'Type','Spearman'; average ranks) between O1 displacement (centroid, ASSD, HD95; eligible
%   records) and primary DMTVF_ref, per structure (6) x metric (3) and for the patient-level median over eligible
%   structures (3); patient bootstrap over evaluable patients, 2,000, percentile 95 % CI, rng(20261111,'twister') reset
%   per correlation; invalid replicates (undefined rho) counted; CI only if >= 90 % valid.
% ST1B: eligible components V_M / V_U (primary); assigned to a structure if > 50 % of native voxel centres lie in its
%   FDG CT mask (ST2a counts) and the structure is displacement-eligible; logistic regression (fitglm, binomial, logit)
%   of V_U on ASSD (per mm) + structure (fixed order; reference = first present); failure criteria (i)-(iv) (SP5) and
%   category loss; patient-cluster bootstrap 2,000, rng(20261112,'twister'); fallback if the full fit fails or > 10 %
%   of replicates fail. The descriptive table is always reported.
%
% Usage (MATLAB, current folder = 02_code):
%   ST1_association_v1_0('selftest')
%   ST1_association_v1_0('run')
% The function calls restoredefaultpath (MATLAB default path for this session; the saved path is not changed).

if nargin < 1, mode = 'selftest'; end
restoredefaultpath;
assertMathWorks({'corr','fitglm','prctile','readtable','writetable','jsonencode','tiedrank','randi','categorical'});
switch lower(mode)
    case 'selftest'
        if ~selftest(), error('ST1:selftest', 'self-test failed'); end
    case 'run'
        runAll();
    otherwise
        error('ST1:mode', 'mode must be selftest or run');
end
end

function S = STRUCT()
S = ["spleen","kidney_right","kidney_left","liver","lung_left","lung_right"];
end

% =============================== run ===============================
function runAll()
t0 = tic; runId = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out = 'D:\DEEP PSMA\PAPER3_Discordance\06_output'; proto = 'D:\DEEP PSMA\PAPER3_Discordance\00_protocol';
P3d = fullfile(out, 'P3_conn26_20261112_001915'); S2d = fullfile(out, 'ST2a_20261112_072150');
inputs = {fullfile(out, 'O1_organ_displacement_20261111_232613.csv'), 'daaae3134932b6267807a128cb8f85dd';
          fullfile(P3d, 'P3_components_20261112_001915.csv'), 'ac0048cfe059c714b3c293fbe7e88279';
          fullfile(P3d, 'P3_patients_20261112_001915.csv'),   '2f88309e1be6eff7c499c6c1d36b50ed';
          fullfile(S2d, 'ST2a_manifest_20261112_072150.json'), 'bca65313ea0d550286fbf9c38bcfc33b';
          fullfile(S2d, 'ST2a_components_20261112_072150.csv'), 'a2fbe377c4081a64de98eec2b6b5f533';
          fullfile(proto, 'StepSAP_soft_tissue_ST_spec_LOCKED_20261111.md'), '055199c4fd5d46039dcc65b6ab4f2416';
          fullfile(proto, 'StepSAP_ST1_implementation_clarification_20261112.md'), '1021a7266020138fc7d0e13b2b4af6a8'};
md = cell(size(inputs,1), 2);
for r = 1:size(inputs,1)
    h = md5file(inputs{r,1}); md(r,:) = {inputs{r,1}, h};
    if ~strcmp(h, inputs{r,2}), error('ST1:C1', 'STOP C1: MD5 mismatch %s (%s)', inputs{r,1}, h); end
end
ms = jsondecode(fileread(inputs{4,1}));
if ~strcmp(ms.status, 'COMPLETE') || ~strcmp(ms.mode, 'FULL') || ...
        ~strcmp(ms.outputs.(matlab.lang.makeValidName('ST2a_components_20261112_072150.csv')), inputs{5,2})
    error('ST1:C1', 'STOP: ST2a component table not verified by its manifest');
end
fprintf('C1 MD5: %d/%d match\n', size(inputs,1), size(inputs,1));
O1 = readtable(inputs{1,1}, 'TextType', 'string'); O1.elig = tobool_v(O1.eligible);
PC = readtable(inputs{2,1}, 'TextType', 'string'); PP = readtable(inputs{3,1}, 'TextType', 'string');
S2 = readtable(inputs{5,1}, 'TextType', 'string');
od = fullfile(out, ['ST1_' runId]); mkdir(od);
metrics = ["centroid_mm","ASSD_mm","HD95_mm"];
SS = STRUCT();

% ---------------- ST1A ----------------
A = cell(0, 11);
for s = SS
    for m = metrics
        E = O1(O1.structure == s & O1.elig, :);
        [tf, loc] = ismember(E.case_id, PP.case_id);
        if ~all(tf) || any(~isfinite(E.(m))), error('ST1:ST1A', 'STOP: O1 record without DMTVF or non-finite metric (%s %s)', s, m); end
        r = rhoBoot(E.(m), PP.DMTVF_ref(loc), 20261111, 2000);
        A(end+1,:) = {s, m, numel(loc), r.rho, r.lo, r.hi, r.nValid, r.nInvalid, r.ciText, NaN, NaN}; %#ok<AGROW>
    end
end
G = O1(O1.elig, :); [u, ~, g] = unique(G.case_id);
PS = table(u, accumarray(g, 1), 'VariableNames', {'case_id','n_structures'});
for m = metrics
    PS.(m + "_median") = accumarray(g, G.(m), [], @median);
end
[tf, loc] = ismember(PS.case_id, PP.case_id); if ~all(tf), error('ST1:ST1A', 'STOP: patient without DMTVF'); end
PS.DMTVF_ref = PP.DMTVF_ref(loc);
for m = metrics
    r = rhoBoot(PS.(m + "_median"), PS.DMTVF_ref, 20261111, 2000);
    A(end+1,:) = {"patient median", m, height(PS), r.rho, r.lo, r.hi, r.nValid, r.nInvalid, r.ciText, ...
        median(PS.n_structures), min(PS.n_structures)}; %#ok<AGROW>
end
TA = cell2table(A, 'VariableNames', {'structure','metric','n_patients','rho','ci_low','ci_high','n_boot_valid', ...
    'n_boot_invalid','ci_text','n_structures_median','n_structures_min'});
fprintf('ST1A (Spearman rho vs primary DMTVF_ref; 95%% patient-bootstrap CI):\n');
for k = 1:height(TA)
    fprintf('   %-15s %-12s n %3d  rho %7.3f  CI %s  (valid %d, invalid %d)\n', TA.structure(k), TA.metric(k), ...
        TA.n_patients(k), TA.rho(k), TA.ci_text(k), TA.n_boot_valid(k), TA.n_boot_invalid(k));
end
fprintf('   patient-level summary: %d patients, contributing structures median %g (min %g)\n', height(PS), ...
    median(PS.n_structures), min(PS.n_structures));

% ---------------- ST1B ----------------
K = PC(PC.status == "V_M" | PC.status == "V_U", :);
key = @(c, q) c + "#" + string(q);
[tf, loc] = ismember(key(K.case_id, K.component_id), key(S2.case_id, S2.component_id));
if ~all(tf), error('ST1:ST1B', 'STOP: component missing in the ST2a table'); end
S2k = S2(loc, :);
if any(S2k.voxels_native <= 0), error('ST1:ST1B', 'STOP: invalid native voxel count'); end
N = zeros(height(K), numel(SS));
for j = 1:numel(SS), N(:,j) = S2k.("n_" + SS(j)); end
eligM = false(size(N)); assd = NaN(size(N));
for j = 1:numel(SS)
    E = O1(O1.structure == SS(j) & O1.elig, :);
    [tfe, le] = ismember(K.case_id, E.case_id);
    eligM(:,j) = tfe; assd(tfe, j) = E.ASSD_mm(le(tfe));
end
[sIdx, why] = assignComp(N, S2k.voxels_native, eligM);
CB = table(K.case_id, K.component_id, K.status, S2k.voxels_native, why, 'VariableNames', ...
    {'case_id','component_id','status','voxels_native','assignment'});
CB.structure = repmat("", height(K), 1); CB.ASSD_mm = NaN(height(K), 1); CB.max_fraction = max(N ./ S2k.voxels_native, [], 2);
a = sIdx > 0;
CB.structure(a) = SS(sIdx(a))';
CB.ASSD_mm(a) = assd(sub2ind(size(assd), find(a), sIdx(a)));
T = table(CB.case_id(a), CB.component_id(a), CB.structure(a), CB.ASSD_mm(a), double(CB.status(a) == "V_U"), ...
    'VariableNames', {'case_id','component_id','structure','ASSD','y'});
fprintf('ST1B: components V_M/V_U %d | assigned and analysed %d (%d patients) | not assigned (<= 50 %%) %d | assigned structure not displacement-eligible %d\n', ...
    height(K), height(T), numel(unique(T.case_id)), sum(why == "not assigned"), sum(why == "structure not displacement-eligible"));
% descriptive table (always)
D = cell(0, 8);
for s = [SS, "all"]
    for yv = [0 1]
        w = T.y == yv; if s ~= "all", w = w & T.structure == s; end
        q = T.ASSD(w);
        D(end+1,:) = {s, tern(yv == 1, "unmatched (V_U)", "matched (V_M)"), sum(w), numel(unique(T.case_id(w))), ...
            median(q), prctile0(q, 25), prctile0(q, 75), NaN}; %#ok<AGROW>
    end
end
TD = cell2table(D, 'VariableNames', {'structure','outcome','n_components','n_patients','ASSD_median_mm', ...
    'ASSD_q25_mm','ASSD_q75_mm','unused'}); TD.unused = [];
fprintf('ST1B descriptive (always reported):\n');
for k = 1:height(TD)
    fprintf('   %-13s %-16s n %4d (patients %3d)  ASSD median %6.2f (IQR %6.2f-%6.2f) mm\n', TD.structure(k), ...
        TD.outcome(k), TD.n_components(k), TD.n_patients(k), TD.ASSD_median_mm(k), TD.ASSD_q25_mm(k), TD.ASSD_q75_mm(k));
end
% model
cats = SS(ismember(SS, unique(T.structure)));
full = fitST1B(T, cats);
refc = ""; if ~isempty(cats), refc = cats(1); end
model = struct('locked_categories', cats, 'reference_category', refc, ...
    'full_fit_ok', full.ok, 'full_fit_reason', full.reason, 'full_fit_warnings', full.warnings, ...
    'OR_per_mm', full.OR, 'logOR', full.logOR, 'n_components', height(T), 'n_patients', numel(unique(T.case_id)), ...
    'bootstrap_run', false, 'n_boot', 2000, 'n_boot_failed', NaN, 'n_boot_category_loss', NaN, 'boot_fail_reasons', "", ...
    'ci_low', NaN, 'ci_high', NaN, 'model_reported', false, 'fallback', true, 'statement', "");
TB = table();
if full.ok
    pats = unique(T.case_id); np = numel(pats);
    idx = cell(np, 1); for i = 1:np, idx{i} = find(T.case_id == pats(i)); end
    rng(20261112, 'twister'); B = 2000;
    okb = false(B,1); lb = NaN(B,1); rs = strings(B,1);
    for b = 1:B
        d = randi(np, np, 1);
        rb = fitST1B(T(vertcat(idx{d}), :), cats);
        okb(b) = rb.ok; lb(b) = rb.logOR; rs(b) = rb.reason;
    end
    TB = table((1:B)', okb, lb, rs, 'VariableNames', {'replicate','ok','logOR','reason'});
    nf = sum(~okb); [ur, ~, gr] = unique(rs(~okb));
    model.bootstrap_run = true; model.n_boot_failed = nf; model.n_boot_category_loss = sum(rs == "category loss");
    if nf > 0, model.boot_fail_reasons = strjoin(ur + ": " + string(accumarray(gr(:), 1)), "; "); else, model.boot_fail_reasons = "none"; end
    if nf <= 0.10 * B
        q = prctile(exp(lb(okb)), [2.5 97.5]); model.ci_low = q(1); model.ci_high = q(2);
        model.model_reported = true; model.fallback = false;
        model.statement = sprintf("model reported: OR per mm %.3f (95%% CI %.3f-%.3f); %d of %d replicates failed", ...
            full.OR, q(1), q(2), nf, B);
    else
        model.statement = sprintf("model NOT reported (SP5 fallback): %d of %d bootstrap fits failed (> 10 %%)", nf, B);
    end
else
    model.statement = "model NOT reported (SP5 fallback): full-data fit failed (" + full.reason + "); bootstrap not run";
end
fprintf('ST1B model: locked categories [%s]; %s\n', strjoin(cats, ', '), model.statement);
if model.bootstrap_run
    fprintf('   bootstrap failures by reason: %s (category loss %d)\n', model.boot_fail_reasons, model.n_boot_category_loss);
end
if strlength(full.warnings) > 0, fprintf('   other warnings recorded (full fit): %s\n', full.warnings); end

% ---------------- outputs ----------------
f = {fullfile(od, ['ST1A_correlations_' runId '.csv']), fullfile(od, ['ST1A_patient_summary_' runId '.csv']), ...
     fullfile(od, ['ST1B_components_' runId '.csv']), fullfile(od, ['ST1B_descriptive_' runId '.csv']), ...
     fullfile(od, ['ST1B_model_' runId '.json'])};
writetable(TA, f{1}); writetable(PS, f{2}); writetable(CB, f{3}); writetable(TD, f{4});
fid = fopen(f{5}, 'w'); fwrite(fid, jsonencode(model, 'PrettyPrint', true)); fclose(fid);
if ~isempty(TB), f{end+1} = fullfile(od, ['ST1B_bootstrap_' runId '.csv']); writetable(TB, f{end}); end
outs = cell(numel(f), 2); for k = 1:numel(f), outs(k,:) = {f{k}, md5file(f{k})}; end
fns = {'corr','fitglm','prctile','tiedrank','randi','categorical'};
wp = struct(); for k = 1:numel(fns), wp.(fns{k}) = which(fns{k}); end
man = struct('stage', 'ST1A / ST1B association of organ displacement with correspondence status (exploratory)', ...
    'run', runId, 'script', 'ST1_association_v1_0.m', 'script_md5', md5file([mfilename('fullpath') '.m']), ...
    'inputs', {md}, 'outputs', {outs}, 'seeds', struct('ST1A', 20261111, 'ST1B', 20261112), 'matlab', version, ...
    'function_paths', wp, 'runtime_s', round(toc(t0)), 'status', 'COMPLETE');
mf = fullfile(od, ['ST1_manifest_' runId '.json']);
fid = fopen(mf, 'w'); fwrite(fid, jsonencode(man, 'PrettyPrint', true)); fclose(fid);
fprintf('done: %s; manifest MD5 %s\n', od, md5file(mf));
end

% =============================== core (tested) ===============================
function r = rhoBoot(x, y, seed, B)
% Spearman rho with patient bootstrap; undefined replicates invalid; CI only if >= 90 % valid
x = x(:); y = y(:); n = numel(x);
r = struct('rho', NaN, 'lo', NaN, 'hi', NaN, 'nValid', 0, 'nInvalid', B, 'ciText', "");
if n >= 2, r.rho = corr(x, y, 'Type', 'Spearman'); end
if isnan(r.rho)
    r.ciText = "not estimable: point estimate undefined"; return;
end
rng(seed, 'twister'); kb = NaN(B, 1);
for b = 1:B
    i = randi(n, n, 1); kb(b) = corr(x(i), y(i), 'Type', 'Spearman');
end
r.nValid = sum(~isnan(kb)); r.nInvalid = B - r.nValid;
if r.nValid >= 0.9 * B
    q = prctile(kb(~isnan(kb)), [2.5 97.5]); r.lo = q(1); r.hi = q(2); r.ciText = sprintf("%.3f to %.3f", q(1), q(2));
else
    r.ciText = sprintf("not estimable: only %d of %d replicates valid (< 90 %%)", r.nValid, B);
end
end

function [sIdx, why] = assignComp(N, nvox, eligM)
% structure index (0 = unassigned) if > 50 % of native voxel centres lie in one structure's mask and it is eligible
F = N ./ nvox(:); [fm, im] = max(F, [], 2);
sIdx = zeros(size(fm)); why = repmat("not assigned", size(fm));
over = fm > 0.5;
el = false(size(fm)); el(over) = eligM(sub2ind(size(eligM), find(over), im(over)));
sIdx(over & el) = im(over & el);
why(over & el) = "assigned"; why(over & ~el) = "structure not displacement-eligible";
end

function res = fitST1B(T, cats)
% locked model y ~ ASSD + structure (fixed category order); SP5 failure criteria (i)-(iv) and category loss
res = struct('ok', false, 'reason', "", 'OR', NaN, 'logOR', NaN, 'warnings', "");
present = STRUCT(); present = present(ismember(present, unique(T.structure)));
if ~isequal(present, cats), res.reason = "category loss"; return; end
if isempty(cats) || all(T.y == T.y(1)), res.reason = "(iii) no outcome variation"; return; end
for c = cats
    yy = T.y(T.structure == c);
    if all(yy == yy(1)), res.reason = "(iii) category with a single outcome"; return; end
end
tbl = table(T.y, T.ASSD, categorical(T.structure, cats), 'VariableNames', {'y','ASSD','structure'});
if numel(cats) > 1, frm = 'y ~ ASSD + structure'; else, frm = 'y ~ ASSD'; end
lastwarn('', ''); mdl = [];
txt = evalc('mdl = fitglm(tbl, frm, ''Distribution'', ''binomial'', ''Link'', ''logit'');');
[wm, wi] = lastwarn;
wt = lower(string(txt) + " " + string(wm) + " " + string(wi));
if contains(wt, "iteration") || contains(wt, "separat")
    res.reason = "(i) iteration-limit or separation warning"; return;
end
if strlength(strtrim(string(wm))) > 0, res.warnings = string(wi) + ": " + string(wm); end
co = mdl.Coefficients;
if any(~isfinite(co.Estimate)) || any(~isfinite(co.SE)), res.reason = "(ii) non-finite coefficient or SE"; return; end
p = mdl.Fitted.Probability;
if any(p < 1e-8 | p > 1 - 1e-8), res.reason = "(iv) fitted probability beyond 1e-8"; return; end
res.logOR = co.Estimate(strcmp(co.Properties.RowNames, 'ASSD')); res.OR = exp(res.logOR); res.ok = true;
end

% =============================== self-test ===============================
function ok = selftest()
res = cell(0, 2);
x = [1;2;2;3;4;4;5]; y = [0;0;1;1;2;0;3];
res(end+1,:) = {'Spearman with ties = Pearson of average ranks', abs(corr(x, y, 'Type', 'Spearman') - corr(tiedrank(x), tiedrank(y))) < 1e-12};
r1 = rhoBoot(ones(6,1), (1:6)', 1, 50);
res(end+1,:) = {'constant variable -> rho undefined, reported as not estimable', isnan(r1.rho) && contains(r1.ciText, "undefined")};
r2 = rhoBoot([0;0;0;0;0;0;0;1], (1:8)', 1, 200);
res(end+1,:) = {sprintf('undefined bootstrap replicates counted as invalid (%d of 200); CI withheld below 90 %% valid', r2.nInvalid), ...
    r2.nInvalid > 20 && r2.nValid + r2.nInvalid == 200 && isnan(r2.lo) && contains(r2.ciText, "not estimable")};
r3 = rhoBoot((1:30)', (1:30)' + 0.1 * sin(1:30)', 1, 200);
res(end+1,:) = {'well-defined rho -> CI reported, all replicates valid', ~isnan(r3.lo) && r3.nValid == 200};
[sI, wy] = assignComp([5 0; 6 0; 3 3; 0 7], [10; 10; 10; 10], logical([1 1; 1 1; 1 1; 1 0]));
res(end+1,:) = {'assignment: exactly 50 % not assigned; 60 % assigned; split 30/30 not; ineligible structure excluded', ...
    isequal(sI, [0; 1; 0; 0]) && wy(4) == "structure not displacement-eligible" && wy(1) == "not assigned"};
S = STRUCT();
T1 = table(repmat("x", 8, 1), (1:8)', [repmat(S(1), 4, 1); repmat(S(4), 4, 1)], [1;2;3;4;1;2;3;4], [0;1;0;1;0;0;0;0], ...
    'VariableNames', {'case_id','component_id','structure','ASSD','y'});
f1 = fitST1B(T1, [S(1) S(4)]);
res(end+1,:) = {'criterion (iii): a category with a single outcome fails', ~f1.ok && startsWith(f1.reason, "(iii)")};
f2 = fitST1B(T1(5:8, :), [S(1) S(4)]);
res(end+1,:) = {'category loss vs the locked set fails (bootstrap safeguard)', ~f2.ok && f2.reason == "category loss"};
a = (1:20)'; T3 = table(repmat("x", 40, 1), (1:40)', [repmat(S(1), 20, 1); repmat(S(4), 20, 1)], [a; a], ...
    double([a > 10; a > 10]), 'VariableNames', {'case_id','component_id','structure','ASSD','y'});
f3 = fitST1B(T3, [S(1) S(4)]);
res(end+1,:) = {sprintf('perfect separation by ASSD fails by (i), (ii) or (iv) [%s]', f3.reason), ~f3.ok && ~startsWith(f3.reason, "(iii)") ...
    && f3.reason ~= "category loss"};
rng(7, 'twister'); n = 6000; as = 10 * rand(n, 1); c = randi(2, n, 1);
yy = double(rand(n, 1) < 1 ./ (1 + exp(-(-2 + 0.2 * as + 0.5 * (c == 2)))));
T4 = table(repmat("x", n, 1), (1:n)', S(1 + 3 * (c == 2))', as, yy, 'VariableNames', {'case_id','component_id','structure','ASSD','y'});
f4 = fitST1B(T4, [S(1) S(4)]);
res(end+1,:) = {sprintf('known log-odds slope 0.20 recovered (%.3f)', f4.logOR), f4.ok && abs(f4.logOR - 0.2) < 0.05};
for k = 1:size(res,1), fprintf('  [%s] %s\n', tern(res{k,2}, 'PASS', 'FAIL'), res{k,1}); end
ok = all([res{:,2}]);
fprintf('selftest: %d/%d passed\n', sum([res{:,2}]), size(res,1));
end

% =============================== utilities ===============================
function v = prctile0(q, p)
if isempty(q), v = NaN; else, v = prctile(q, p); end
end

function assertMathWorks(fns)
root = lower(matlabroot);
for k = 1:numel(fns)
    w = regexprep(which(fns{k}), '^built-in \((.*)\)$', '$1');
    if isempty(w) || ~startsWith(lower(w), root)
        error('ST1:path', 'STOP: %s resolves to "%s" (not the MathWorks installation)', fns{k}, w);
    end
end
end

function h = md5file(f)
md = java.security.MessageDigest.getInstance('MD5');
fid = fopen(f, 'r'); if fid < 0, error('ST1:io', 'cannot open %s', f); end
data = fread(fid, Inf, '*uint8'); fclose(fid);
md.update(data);
h = lower(reshape(dec2hex(typecast(md.digest, 'uint8'), 2)', 1, []));
end

function b = tobool_v(x)
if islogical(x), b = x; elseif isnumeric(x), b = ~isnan(x) & x ~= 0; else, b = strcmpi(strtrim(string(x)), "true"); end
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
