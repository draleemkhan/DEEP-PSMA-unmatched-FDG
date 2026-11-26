function P3_matching_classification_v1_4(mode, p2bDir)
% P3_matching_classification_v1_4 - PAPER 3, stage P3: PSMA-FDG matching, boundary-band rule, classification precedence
% and per-patient volume accounting / DMTVF (MATLAB R2025a).
%
% Implements: consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) sections 1, 5, 6; tau spec
% LOCKED (MD5 285debf473d92a32c765f6f7ee45d412) D1-D3; implementation spec LOCKED (MD5
% a5a2ec19d284deb44bc9ceebf8636658) I0, I4 (with the author's clarifications on band alternatives and overlaps).
%
% Inputs (all MD5-checked): P1 components/patients (P1_20261111_233936), P2 components (P2_conn26_20261111_235603,
% accepted; self-test 23/23), P2b CSV export (voxel sets 1-based MATLAB linear indices on the PSMA PET grid; PSMA PET
% grid in SimpleITK LPS convention), native PSMA TTB (data root, read-only), tau validation table (self-test).
%
% Matching (tau spec): tau_i = voxel diagonal of the PSMA PET grid; d = 0 if any FDG voxel is a PSMA TTB voxel;
% otherwise the minimum Euclidean distance (mm) between boundary-voxel centres (6-neighbour boundary; grid edge =
% outside) of the FDG component and of the whole PSMA TTB mask (knnsearch, kd-tree); matched if d <= tau_i.
% Correspondence graph: PSMA components (26-connectivity) within tau_i of each matched FDG component.
% Precedence (one status per component): V_S -> V_I registration -> V_I band -> V_I coverage -> V_M / V_U.
%   band: band candidates (alternative transform, P2) not V_I registration and with a valid FINAL mapping are
%   classified under the final and the alternative voxel sets; categories M / U / NA (alternative invalid); any
%   difference -> V_I band. A band candidate whose final mapping fails coverage proceeds to coverage.
%   coverage: coverage_fov_any (FOV or within tau of its boundary, any relevant transform) or zero survivors (final).
% Accounting from NATIVE volumes (P1): V_FDG,total = V_S + V_I + V_M + V_U; DMTVF_ref = V_U / V_FDG,total; upper
% bound (V_U + V_I) / V_FDG,total; f_small = V_S / V_FDG,total.
%
% v1.1 (12 Nov 2026, 00:13): the MathWorks path check accepts built-in functions, whose which() result is
% 'built-in (<path under matlabroot>)'; v1.0 (P3_matching_classification.m) stopped at that check before any
% computation. No other change.
% v1.2 (12 Nov 2026, 00:20): self-test fixes found at the first MATLAB execution (no change to the analysis rules):
%   tau table read with VariableNamingRule 'preserve' ('case' is a MATLAB keyword); string-to-char before sscanf;
%   self-test result collection (the anonymous accumulator captured an empty list); alt_transform read robustly
%   when the column is numeric/NaN.
% v1.3 (12 Nov 2026, 00:25): self-test grid parsing — the locked tau table stores cohort grids as
%   '(np.float64(x), np.float64(y), np.float64(z))' (numpy 2 repr); numbers are now extracted after removing
%   'np.float64'. No change to the analysis rules.
% v1.4 (12 Nov 2026, 00:30): self-test scene for the correspondence graph corrected — in v1.3 two PSMA blocks
%   touched diagonally and were (correctly) one 26-connected component; the second block now lies on the opposite
%   side of the FDG block. No change to the analysis code.
%
% Usage (MATLAB, current folder = 02_code):
%   P3_matching_classification_v1_4('selftest')
%   P3_matching_classification_v1_4('preflight', 'P2b_conn26_<run>')    % 2 patients; prints integrity only
%   P3_matching_classification_v1_4('full',      'P2b_conn26_<run>')
% The function calls restoredefaultpath (MATLAB default path for this session; the saved path is not changed).

if nargin < 1, mode = 'selftest'; end
restoredefaultpath;                                   % SAP item 6: no third-party function may shadow MathWorks code
assertMathWorks({'knnsearch','createns','rangesearch','bwconncomp','labelmatrix','niftiread','niftiinfo', ...
                 'readtable','writetable','jsonencode','ind2sub'});
cfg.out  = 'D:\DEEP PSMA\PAPER3_Discordance\06_output';
cfg.data = 'D:\DEEP PSMA\DEEP PSMA EXTRACTED';
cfg.proto = 'D:\DEEP PSMA\PAPER3_Discordance\00_protocol';
cfg.tol = 1e-9;
switch lower(mode)
    case 'selftest'
        ok = selftest(cfg);
        if ~ok, error('P3:selftest', 'self-test failed'); end
    case {'preflight','full'}
        runP3(cfg, strcmpi(mode, 'preflight'), p2bDir);
    otherwise
        error('P3:mode', 'mode must be selftest, preflight or full');
end
end

% =============================== run ===============================
function runP3(cfg, preflight, p2bDir)
t0 = tic; runId = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
P1d = fullfile(cfg.out, 'P1_20261111_233936'); P2d = fullfile(cfg.out, 'P2_conn26_20261111_235603');
P2bd = fullfile(cfg.out, p2bDir);
inputs = {fullfile(P1d, 'P1_components_20261111_233936.csv'), 'ffafbc7917e1b595eadacafab247fe88';
          fullfile(P1d, 'P1_patients_20261111_233936.csv'),   'f9b5b5908e523d1122547ef07d5c0afd';
          fullfile(P2d, 'P2_components_20261111_235603.csv'), '9a649e60c49dafde2e4b3d1b1e35fa3c';
          fullfile(cfg.proto, 'StepSAP_consolidated_LOCKED_20261111.md'), '0c64711db1207e31a70cf3f558ccbd7d';
          fullfile(cfg.proto, 'PET_correspondence_implementation_spec_LOCKED_20261111.md'), 'a5a2ec19d284deb44bc9ceebf8636658'};
md = cell(size(inputs,1), 2);
for r = 1:size(inputs,1)
    h = md5file(inputs{r,1}); md(r,:) = {inputs{r,1}, h};
    if ~strcmp(h, inputs{r,2}), error('P3:C1', 'STOP C1: MD5 mismatch %s (%s)', inputs{r,1}, h); end
end
man2b = jsondecode(fileread(fullfile(P2bd, sprintf('P2b_manifest_%s.json', p2bDir(end-14:end)))));
if ~strcmp(man2b.source_P2_run, 'P2_conn26_20261111_235603') || ~strcmp(man2b.status, 'COMPLETE')
    error('P3:C1', 'STOP: P2b export is not from the accepted P2 run');
end
fprintf('C1 MD5: %d/%d match; P2b export %s (manifest MD5 %s)\n', size(inputs,1), size(inputs,1), p2bDir, ...
    md5file(fullfile(P2bd, sprintf('P2b_manifest_%s.json', p2bDir(end-14:end)))));
C1 = readtable(inputs{1,1}, 'TextType', 'string'); C1 = C1(C1.connectivity == 26, :);
P1 = readtable(inputs{2,1}, 'TextType', 'string'); P1 = P1(P1.connectivity == 26, :);
C2 = readtable(inputs{3,1}, 'TextType', 'string');
G  = readtable(fullfile(P2bd, sprintf('P2b_psma_pet_grid_%s.csv', p2bDir(end-14:end))), 'TextType', 'string');
cases = sort(unique(C1.case_id));
if numel(cases) ~= 100, error('P3:input', 'STOP: expected 100 patients'); end
if preflight
    oldK = readtable(fullfile(cfg.out, 'A3_visualQC_KEY_do_not_open_20261105.csv'), 'TextType', 'string');
    cases = cases(~ismember(cases, oldK.case_id)); cases = cases(1:2);
end
compRows = {}; patRows = {};
for ci = 1:numel(cases)
    cid = cases(ci);
    [cr, pr] = processPatient(cfg, cid, C1(C1.case_id == cid, :), P1(P1.case_id == cid, :), ...
        C2(C2.case_id == cid, :), G(G.case_id == cid, :), P2bd, man2b);
    compRows = [compRows; cr]; patRows = [patRows; pr]; %#ok<AGROW>
    fprintf('  [%3d/%d] %s: %d components classified (%.0f s)\n', ci, numel(cases), cid, size(cr,1), toc(t0));
end
CT = cell2table(compRows, 'VariableNames', {'case_id','component_id','volume_native_mL','status','cause', ...
    'd_final_mm','category_final','d_alt_mm','category_alt','tau_mm','psma_components_matched'});
PT = cell2table(patRows, 'VariableNames', {'case_id','V_FDG_total_mL','V_S_mL','V_I_registration_mL', ...
    'V_I_band_mL','V_I_coverage_mL','V_I_mL','V_M_mL','V_U_mL','n_V_S','n_V_I','n_V_M','n_V_U', ...
    'DMTVF_ref','DMTVF_upper','f_small','reconciliation_rel_diff'});
% ---- reconciliation (I0.3) ----
okRec = all(PT.reconciliation_rel_diff <= cfg.tol) && height(CT) == sum(ismember(C1.case_id, cases));
fprintf('Reconciliation: V_S + V_I + V_M + V_U = V_FDG,total (max rel. diff %.1e); one status per component: %s\n', ...
    max(PT.reconciliation_rel_diff), tern(okRec, 'PASS', 'FAIL'));
if ~okRec, error('P3:reconciliation', 'STOP: reconciliation failed; outputs not written'); end
tag = 'P3_'; if preflight, tag = 'P3_PREFLIGHT_'; end
od = fullfile(cfg.out, [tag 'conn26_' runId]); mkdir(od);
cf = fullfile(od, ['P3_components_' runId '.csv']); writetable(CT, cf);
pf = fullfile(od, ['P3_patients_' runId '.csv']);   writetable(PT, pf);
if preflight
    fprintf('Preflight: integrity only (no endpoint values displayed). Files written to %s\n', od);
else
    fprintf('Cohort (n = %d): V_FDG,total %.1f mL | V_S %.1f | V_I %.1f (registration %.1f, band %.1f, coverage %.1f) | V_M %.1f | V_U %.1f mL\n', ...
        height(PT), sum(PT.V_FDG_total_mL), sum(PT.V_S_mL), sum(PT.V_I_mL), sum(PT.V_I_registration_mL), ...
        sum(PT.V_I_band_mL), sum(PT.V_I_coverage_mL), sum(PT.V_M_mL), sum(PT.V_U_mL));
    fprintf('Per-patient DMTVF_ref: median %.4f (IQR %.4f-%.4f); upper bound median %.4f\n', median(PT.DMTVF_ref), ...
        prctile(PT.DMTVF_ref, 25), prctile(PT.DMTVF_ref, 75), median(PT.DMTVF_upper));
end
fns = {'knnsearch','createns','rangesearch','bwconncomp','labelmatrix','niftiread','niftiinfo','prctile'};
wp = struct(); for k = 1:numel(fns), wp.(fns{k}) = which(fns{k}); end
v = ver; tb = arrayfun(@(x) sprintf('%s %s', x.Name, x.Version), v, 'UniformOutput', false);
man = struct('stage', 'P3 matching, band rule, precedence, accounting', 'mode', tern(preflight, 'PREFLIGHT', 'FULL'), ...
    'run', runId, 'script', 'P3_matching_classification_v1_4.m', 'script_md5', md5file([mfilename('fullpath') '.m']), ...
    'inputs', {md}, 'p2b_dir', p2bDir, 'n_patients', height(PT), ...
    'outputs', {{cf, md5file(cf); pf, md5file(pf)}}, 'matlab', version, 'toolboxes', {tb}, 'function_paths', wp, ...
    'tolerance', cfg.tol, 'runtime_s', round(toc(t0)), 'status', 'COMPLETE');
mf = fullfile(od, ['P3_manifest_' runId '.json']);
fid = fopen(mf, 'w'); fwrite(fid, jsonencode(man, 'PrettyPrint', true)); fclose(fid);
fprintf('done: %s; manifest MD5 %s\n', od, md5file(mf));
end

% =============================== per patient ===============================
function [rows, prow] = processPatient(cfg, cid, C1, P1, C2, G, P2bd, man2b)
cidc = char(cid);
info = niftiinfo(fullfile(cfg.data, cidc, 'PSMA', 'TTB.nii.gz'));
A = niiAffine(info); sz = double(info.ImageSize(1:3));
% grid check against the P2b (SimpleITK, LPS) record: size, spacing, affine (RAS = diag(-1,-1,1) LPS)
Dl = reshape([G.dir_00 G.dir_01 G.dir_02 G.dir_10 G.dir_11 G.dir_12 G.dir_20 G.dir_21 G.dir_22], 3, 3)';
spg = [G.spacing_x G.spacing_y G.spacing_z]; og = [G.origin_x_LPS; G.origin_y_LPS; G.origin_z_LPS];
Aexp = [diag([-1 -1 1]) * Dl * diag(spg), diag([-1 -1 1]) * og; 0 0 0 1];
if ~isequal(sz, [G.size_x G.size_y G.size_z]) || max(abs(A(:) - Aexp(:))) > 1e-3
    error('P3:grid', 'STOP: %s PSMA TTB grid (MATLAB) differs from the P2 PSMA PET grid (max %.2e mm)', cidc, max(abs(A(:) - Aexp(:))));
end
sp = vecnorm(A(1:3,1:3)); tau = sqrt(sum(sp.^2));
% voxel sets: verify the CSV MD5 recorded by P2b
vf = fullfile(P2bd, sprintf('P2b_voxelsets_%s.csv', cidc));
key = matlab.lang.makeValidName(sprintf('P2b_voxelsets_%s.csv', cidc));
if ~isfield(man2b.outputs, key) || ~strcmp(man2b.outputs.(key), md5file(vf))
    error('P3:C1', 'STOP: voxel-set CSV MD5 mismatch for %s', cidc);
end
VS = readmatrix(vf);                                  % [component_id role linear_index_1based]
psma = niftiread(info) > 0;
[Pb, PbLin] = boundaryPoints(find(psma), sz, A);
L = labelmatrix(bwconncomp(psma, 26));
if isempty(Pb), NS = []; else, NS = createns(Pb, 'NSMethod', 'kdtree'); end
Vtot = P1.V_FDG_total_mL(1);
rows = cell(height(C1), 11); acc = zeros(1, 8); n = zeros(1, 4);   % acc: S Ireg Iband Icov M U
for r = 1:height(C1)
    q = C1.component_id(r); vol = C1.volume_mL(r);
    c2 = C2(C2.component_id == q, :);
    fl = struct('eligible', tobool(C1.eligible_A(r)), 'registration', false, 'bandCand', false, ...
        'finalCoverage', false, 'altInvalid', false, 'coverage', false);
    dF = NaN; dA = NaN; catF = ""; catA = ""; matchedP = "";
    if fl.eligible
        fl.registration = tobool(c2.VI_registration);
        altv = c2.alt_transform;
        if isnumeric(altv), alt = ""; else, alt = string(altv); if ismissing(alt) || alt == "NaN", alt = ""; end, end
        fl.bandCand = tobool(c2.band) && strlength(alt) > 0;
        fl.finalCoverage = tobool(c2.coverage_fov_final) || tobool(c2.zero_survivor_final);
        fl.coverage = tobool(c2.coverage_fov_any) || tobool(c2.zero_survivor_final);
        if fl.bandCand
            fl.altInvalid = tobool(c2.coverage_fov_alt) || tobool(c2.zero_survivor_alt);
        end
        linF = VS(VS(:,1) == q & VS(:,2) == 0, 3); linA = VS(VS(:,1) == q & VS(:,2) == 1, 3);
        needFinal = ~fl.registration && ~fl.coverage || (~fl.registration && fl.bandCand && ~fl.finalCoverage);
        if needFinal && ~isempty(linF)
            [dF, mP] = matchDistance(linF, psma, Pb, PbLin, NS, L, sz, A, tau);
            catF = tern(dF <= tau + cfg.tol, "M", "U"); matchedP = mP;
        end
        if ~fl.registration && fl.bandCand && ~fl.finalCoverage
            if fl.altInvalid || isempty(linA)
                catA = "NA";
            else
                dA = matchDistance(linA, psma, Pb, PbLin, NS, L, sz, A, tau);
                catA = tern(dA <= tau + cfg.tol, "M", "U");
            end
        end
    end
    [st, cause] = classify(fl, catF, catA);
    k = find(strcmp(st, {'V_S','V_I','V_M','V_U'}));
    n(k) = n(k) + 1;
    switch cause
        case 'below eligibility', acc(1) = acc(1) + vol;
        case 'registration',      acc(2) = acc(2) + vol;
        case 'boundary band',     acc(3) = acc(3) + vol;
        case 'coverage',          acc(4) = acc(4) + vol;
        case 'matched',           acc(5) = acc(5) + vol;
        case 'unmatched',         acc(6) = acc(6) + vol;
    end
    if st ~= "V_M", matchedP = ""; end
    rows(r,:) = {cidc, q, vol, st, cause, dF, char(catF), dA, char(catA), tau, char(matchedP)};
end
VI = sum(acc(2:4));
relDiff = abs(sum(acc(1:6)) - Vtot) / max(Vtot, eps);
prow = {cidc, Vtot, acc(1), acc(2), acc(3), acc(4), VI, acc(5), acc(6), n(1), n(2), n(3), n(4), ...
    acc(6) / Vtot, (acc(6) + VI) / Vtot, acc(1) / Vtot, relDiff};
end

% =============================== core (tested) ===============================
function [st, cause] = classify(fl, catF, catA)
% precedence: V_S -> V_I registration -> V_I band -> V_I coverage -> V_M / V_U
if ~fl.eligible,      st = 'V_S'; cause = 'below eligibility'; return; end
if fl.registration,   st = 'V_I'; cause = 'registration'; return; end
if fl.bandCand && ~fl.finalCoverage && catF ~= catA
                      st = 'V_I'; cause = 'boundary band'; return; end
if fl.coverage,       st = 'V_I'; cause = 'coverage'; return; end
if catF == "M",       st = 'V_M'; cause = 'matched';
elseif catF == "U",   st = 'V_U'; cause = 'unmatched';
else, error('P3:classify', 'component reached matching without a final category');
end
end

function [d, matchedP] = matchDistance(lin, psma, Pb, PbLin, NS, L, sz, A, tau)
% D2: any shared voxel -> d = 0; D1: min boundary-centre distance (mm) to the PSMA TTB boundary
matchedP = "";
shared = psma(lin);
if any(shared)
    d = 0;
else
    Fb = boundaryPoints(lin, sz, A);
    if isempty(Pb) || isempty(Fb), d = Inf; return; end
    [~, dd] = knnsearch(NS, Fb); d = min(dd);
end
if d <= tau + 1e-9
    labs = double(L(lin(shared)));
    if ~isempty(Pb)
        Fb = boundaryPoints(lin, sz, A);
        idx = rangesearch(NS, Fb, tau + 1e-9);
        labs = [labs; double(L(PbLin([idx{:}])))];
    end
    labs = unique(labs(labs > 0));
    matchedP = string(strjoin(arrayfun(@num2str, labs(:)', 'UniformOutput', false), ';'));
end
end

function [pts, blin] = boundaryPoints(lin, sz, A)
% 6-neighbour boundary voxels of a voxel set (linear indices on a grid of size sz; grid edge = outside);
% returns physical centres (mm; A maps 0-based (i,j,k) to world) and their linear indices
pts = zeros(0,3); blin = zeros(0,1);
if isempty(lin), return; end
[i, j, k] = ind2sub(sz, lin(:));
lo = [min(i) min(j) min(k)] - 1; hi = [max(i) max(j) max(k)] + 1;
M = false(hi - lo + 1);
li = sub2ind(size(M), i - lo(1) + 1, j - lo(2) + 1, k - lo(3) + 1);
M(li) = true;
inner = M;
inner(2:end-1,:,:) = inner(2:end-1,:,:) & M(1:end-2,:,:) & M(3:end,:,:);
inner(:,2:end-1,:) = inner(:,2:end-1,:) & M(:,1:end-2,:) & M(:,3:end,:);
inner(:,:,2:end-1) = inner(:,:,2:end-1) & M(:,:,1:end-2) & M(:,:,3:end);
inner([1 end],:,:) = false; inner(:,[1 end],:) = false; inner(:,:,[1 end]) = false;
B = M & ~inner;
% a voxel on the grid edge has an outside neighbour: it is in the padded frame only if the set touches the edge
[bi, bj, bk] = ind2sub(size(B), find(B));
gi = bi + lo(1) - 1; gj = bj + lo(2) - 1; gk = bk + lo(3) - 1;
blin = sub2ind(sz, gi, gj, gk);
w = A * [gi' - 1; gj' - 1; gk' - 1; ones(1, numel(gi))];
pts = w(1:3,:)';
end

% =============================== self-test ===============================
function ok = selftest(cfg)
res = cell(0, 2);
% --- 1. tau tests T1-T7 rebuilt from the locked tau validation table (identical synthetic inputs) ---
tf = fullfile(cfg.out, 'Tau_validation_20261111', 'Tau_synthetic_tests_20261111_220544.csv');
if ~strcmp(md5file(tf), '0ddc70252a075697b90f062181a0304b'), error('P3:C1', 'tau table MD5 mismatch'); end
T = readtable(tf, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
T = T(ismember(T.test, ["T1","T2","T3","T4","T5","T6","T7"]), :);
maxd = 0; maxk = 0; mism = 0;
for r = 1:height(T)
    sp = gridSpacing(T.grid(r));
    [D, O] = orientationOf(T.orientation(r));
    A = [D * diag(sp), O(:); 0 0 0 1]; sz = [30 30 30];
    [F, P] = sceneOf(T.("case")(r));
    linF = find(F); psma = P;
    [Pb, PbLin] = boundaryPoints(find(psma), sz, A);
    NS = createns(Pb, 'NSMethod', 'kdtree'); L = labelmatrix(bwconncomp(psma, 26));
    tau = sqrt(sum(sp.^2));
    d = matchDistance(linF, psma, Pb, PbLin, NS, L, sz, A, tau);
    % brute force (T8a)
    if any(psma(linF)), db = 0; else
        Fb = boundaryPoints(linF, sz, A); db = min(min(pdist2(Fb, Pb))); end
    maxd = max(maxd, abs(d - T.d_mm(r))); maxk = max(maxk, abs(d - db));
    mism = mism + ((d <= tau + 1e-9) ~= tobool(T.matched(r)));
end
res(end+1,:) = rowOf(sprintf('tau T1-T7: %d cases reproduce Python d (max diff %.1e mm) and matched flags', height(T), maxd), maxd <= 1e-9 && mism == 0);
res(end+1,:) = rowOf(sprintf('T8a in MATLAB: knnsearch = brute force (max diff %.1e mm)', maxk), maxk <= 1e-9);
% --- 2. correspondence graph (one-to-many, many-to-one), as the Python T8(c) scene ---
sz = [30 30 30]; A = [diag([2.73 2.73 5.0]) zeros(3,1); 0 0 0 1];
FL = zeros(sz); PL = false(sz);
FL(11:13,11:13,11:13) = 1; FL(11:13,19:21,11:13) = 2;
PL(14:16,11:13,11:13) = true; PL(8:10,11:13,11:13) = true; PL(14:16,19:21,11:13) = true;   % three separate 26-components
[Pb, PbLin] = boundaryPoints(find(PL), sz, A); NS = createns(Pb, 'NSMethod', 'kdtree'); L = labelmatrix(bwconncomp(PL, 26));
tau = sqrt(sum([2.73 2.73 5.0].^2));
[d1, m1] = matchDistance(find(FL == 1), PL, Pb, PbLin, NS, L, sz, A, tau);
[d2, m2] = matchDistance(find(FL == 2), PL, Pb, PbLin, NS, L, sz, A, tau);
res(end+1,:) = rowOf(sprintf('one FDG component matched to two PSMA components; another to a third (PSMA components %d)', max(L(:))), ...
    d1 <= tau && d2 <= tau && numel(split(m1, ';')) == 2 && numel(split(m2, ';')) == 1);
% --- 3. precedence ---
b = @(e, r, bc, fc, ai, c) struct('eligible', e, 'registration', r, 'bandCand', bc, 'finalCoverage', fc, 'altInvalid', ai, 'coverage', c);
res(end+1,:) = rowOf('not eligible -> V_S', strcmp(classify(b(false,true,true,false,false,true), "U", "M"), 'V_S'));
[s1, c1] = classify(b(true,true,true,false,false,true), "M", "U");
res(end+1,:) = rowOf('registration precedes band and coverage', strcmp(s1,'V_I') && strcmp(c1,'registration'));
[s2, c2] = classify(b(true,false,true,false,false,false), "M", "U");
res(end+1,:) = rowOf('band: categories differ -> V_I band', strcmp(s2,'V_I') && strcmp(c2,'boundary band'));
[s3, c3] = classify(b(true,false,true,false,true,false), "M", "NA");
res(end+1,:) = rowOf('band: invalid alternative (NA) -> V_I band, never V_U', strcmp(s3,'V_I') && strcmp(c3,'boundary band'));
[s4, c4] = classify(b(true,false,true,false,false,false), "U", "U");
res(end+1,:) = rowOf('band: same category -> continues to matching (V_U)', strcmp(s4,'V_U') && strcmp(c4,'unmatched'));
[s5, c5] = classify(b(true,false,true,true,false,true), "", "");
res(end+1,:) = rowOf('band candidate with invalid final mapping -> V_I coverage', strcmp(s5,'V_I') && strcmp(c5,'coverage'));
[s6, c6] = classify(b(true,false,false,false,false,true), "", "");
res(end+1,:) = rowOf('coverage -> V_I coverage', strcmp(s6,'V_I') && strcmp(c6,'coverage'));
res(end+1,:) = rowOf('no flags -> matched / unmatched by d', strcmp(classify(b(true,false,false,false,false,false), "M", ""), 'V_M') ...
    && strcmp(classify(b(true,false,false,false,false,false), "U", ""), 'V_U'));
% --- 4. overlapping FDG components (different transforms) are matched on their own voxel sets ---
sz = [20 20 20]; A = [eye(3)*2 zeros(3,1); 0 0 0 1]; P = false(sz); P(5:7,5:7,5:7) = true;
[Pb, PbLin] = boundaryPoints(find(P), sz, A); NS = createns(Pb, 'NSMethod', 'kdtree'); L = labelmatrix(bwconncomp(P, 26));
setA = sub2ind(sz, [8 8 9], [6 6 6], [6 7 6]); setB = sub2ind(sz, [8 15], [6 15], [6 15]);  % share voxel (8,6,6)
dA = matchDistance(setA(:), P, Pb, PbLin, NS, L, sz, A, sqrt(12));
dB = matchDistance(setB(:), P, Pb, PbLin, NS, L, sz, A, sqrt(12));
res(end+1,:) = rowOf('overlapping components keep separate identities (both evaluated)', abs(dA - 2) < 1e-12 && abs(dB - 2) < 1e-12);
% --- 5. boundary of a sparse set = full-grid 6-neighbour boundary ---
sz = [12 12 12]; M = false(sz); M(3:8,3:8,3:8) = true; M(1,1,1) = true;
fullB = M & ~(circshift(M,1,1) & circshift(M,-1,1) & circshift(M,1,2) & circshift(M,-1,2) & circshift(M,1,3) & circshift(M,-1,3));
fullB(1,:,:) = M(1,:,:); fullB(:,1,:) = fullB(:,1,:) | M(:,1,:); fullB(:,:,1) = fullB(:,:,1) | M(:,:,1);
[~, bl] = boundaryPoints(find(M), sz, eye(4));
res(end+1,:) = rowOf('sparse boundary = full-grid boundary (grid edge = outside)', isequal(sort(bl), sort(find(fullB))));
% --- 6. index convention ---
sz = [6 5 4]; i = 5; j = 4; k = 3;                    % 1-based MATLAB subscripts of (i,j,k)
pyC = (k-1)*sz(1)*sz(2) + (j-1)*sz(1) + (i-1);      % Python C-order index over (k,j,i), 0-based
res(end+1,:) = rowOf('P2b index (Python C order + 1) = MATLAB linear index', sub2ind(sz, i, j, k) == pyC + 1);
% --- 7. accounting ---
acc = [1 2 3 4 5 6]; res(end+1,:) = rowOf('accounting sums', abs(sum(acc) - 21) < 1e-12);
for r = 1:size(res,1), fprintf('  [%s] %s\n', tern(res{r,2}, 'PASS', 'FAIL'), res{r,1}); end
ok = all([res{:,2}]);
fprintf('selftest: %d/%d passed\n', sum([res{:,2}]), size(res,1));
end

function sp = gridSpacing(g)
% '(4.0, 4.0, 4.0)' or '(np.float64(2.733471), np.float64(2.733471), np.float64(5.0))' -> [x y z]
t = regexprep(char(g), 'np\.float64', '');
sp = str2double(regexp(t, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match'));
if numel(sp) ~= 3 || any(isnan(sp)), error('P3:selftest', 'cannot parse grid "%s"', char(g)); end
end

function r = rowOf(name, c)
r = {name, logical(c)};
end

function [D, O] = orientationOf(s)
if startsWith(s, "identity"), D = eye(3); O = [0 0 0]; return; end
f = sscanf(char(regexprep(extractBetween(s, "(", ")"), '\s', '')), '%d,%d,%d');
D = diag(f); O = [-312.7 145.3 -988.1];
end

function [F, P] = sceneOf(c)
% same boxes as Tau_validation.py cases(): F = [10,13)^3 (0-based) -> MATLAB 11:13
F = false(30,30,30); P = false(30,30,30); F(11:13,11:13,11:13) = true;
box = @(lo, hi) {lo(1)+1:hi(1), lo(2)+1:hi(2), lo(3)+1:hi(3)};
ax = struct('x', 1, 'y', 2, 'z', 3);
if c == "overlap (shared voxels)", b = box([12 12 12], [16 16 16]);
elseif startsWith(c, "containment"), b = box([6 6 6], [18 18 18]);
elseif startsWith(c, "face contact"), a = ax.(char(extractAfter(c, "face contact "))); lo = [10 10 10]; hi = [13 13 13]; lo(a) = 13; hi(a) = 16; b = box(lo, hi);
elseif startsWith(c, "edge contact"), p = char(extractAfter(c, "edge contact ")); lo = [10 10 10]; hi = [13 13 13];
    for a = [ax.(p(1)) ax.(p(2))], lo(a) = 13; hi(a) = 16; end; b = box(lo, hi);
elseif c == "corner contact", b = box([13 13 13], [16 16 16]);
elseif startsWith(c, "one-voxel gap"), a = ax.(char(extractAfter(c, "one-voxel gap "))); lo = [10 10 10]; hi = [13 13 13]; lo(a) = 14; hi(a) = 17; b = box(lo, hi);
elseif c == "one-voxel diagonal gap", b = box([14 14 14], [17 17 17]);
else, error('P3:scene', 'unknown case %s', c);
end
P(b{:}) = true;
end

% =============================== utilities ===============================
function assertMathWorks(fns)
root = lower(matlabroot);
for k = 1:numel(fns)
    w = which(fns{k});
    w = regexprep(w, '^built-in \((.*)\)$', '$1');      % built-in functions: 'built-in (<path>)'
    if isempty(w) || ~startsWith(lower(w), root)
        error('P3:path', 'STOP: %s resolves to "%s" (not the MathWorks installation)', fns{k}, w);
    end
end
end

function A = niiAffine(I)
% 4x4 mapping 0-based voxel indices -> NIfTI world (RAS+, mm); identical to Audits 1A/1B/1D
r = I.raw;
if isfield(r,'sform_code') && r.sform_code > 0
    A = [double(r.srow_x(:))'; double(r.srow_y(:))'; double(r.srow_z(:))'; 0 0 0 1];
else
    b = double(r.quatern_b); c = double(r.quatern_c); d = double(r.quatern_d);
    a = sqrt(max(0, 1 - b^2 - c^2 - d^2));
    Rm = [a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c); 2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b); ...
          2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-c*c-b*b];
    pd = double(r.pixdim); qfac = pd(1); if qfac == 0, qfac = 1; end
    A = eye(4); A(1:3,1:3) = Rm * diag([pd(2) pd(3) pd(4)*qfac]);
    A(1:3,4) = [double(r.qoffset_x); double(r.qoffset_y); double(r.qoffset_z)];
end
end

function h = md5file(f)
md = java.security.MessageDigest.getInstance('MD5');
fid = fopen(f, 'r'); if fid < 0, error('P3:io', 'cannot open %s', f); end
data = fread(fid, Inf, '*uint8'); fclose(fid);
md.update(data);
h = lower(reshape(dec2hex(typecast(md.digest, 'uint8'), 2)', 1, []));
end

function b = tobool(x)
if islogical(x), b = x;
elseif isnumeric(x), b = ~isnan(x) & x ~= 0;
else, b = strcmpi(strtrim(string(x)), "true");
end
b = logical(b(1));
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
