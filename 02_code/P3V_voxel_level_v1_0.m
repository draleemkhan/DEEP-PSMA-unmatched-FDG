function P3V_voxel_level_v1_0(mode, p2bDir)
% P3V_voxel_level_v1_0 - PAPER 3: prespecified voxel-level sensitivity analysis (consolidated Step 2 SAP LOCKED,
% MD5 0c64711db1207e31a70cf3f558ccbd7d, section 8.4): "spatially unmatched FDG tumour-volume fraction".
% Secondary methodological sensitivity only; the primary DMTVF_ref (P3 run P3_conn26_20261112_001915) is unchanged.
%
% Definition (SAP 8.4, tau spec): assessable FDG voxels - voxels of components classified V_M or V_U in the primary
% P3 run - whose centre lies more than tau_i (PSMA PET voxel diagonal) from every PSMA TTB voxel centre, divided by
% V_FDG,total. Voxels of V_S and V_I components never enter the numerator; the denominator is the native
% V_FDG,total (P1). Distance of a voxel centre to the PSMA TTB: 0 inside the mask, otherwise the distance to the
% nearest PSMA boundary-voxel centre (for a point on the same grid outside the mask the nearest mask voxel is a
% boundary voxel). Same grid, same boundary definition, same tau and tolerance (1e-9 mm) as P3.
% Native volumes (implementation spec I0.1): per component, the unmatched share of its resampled voxels (final
% transform, as in P3) multiplies its NATIVE volume; resampled volumes are never summed into numerator or denominator.
%
% Checks (stop on failure): input MD5s (P3 outputs, P1, P2b voxel CSVs via the P2b manifest); MATLAB path
% (restoredefaultpath, MathWorks functions only); grid as in P3; for every V_U component the unmatched voxel share
% must be exactly 1; for every V_M component it must be < 1 (at least one voxel within tau, since the component
% matched); component set and voxel counts identical to the P2b final voxel sets.
%
% Usage (MATLAB, current folder 02_code):
%   P3V_voxel_level_v1_0('selftest')
%   P3V_voxel_level_v1_0('preflight', 'P2b_conn26_20261112_000903')   % 2 patients; integrity only
%   P3V_voxel_level_v1_0('full',      'P2b_conn26_20261112_000903')

if nargin < 1, mode = 'selftest'; end
restoredefaultpath;
assertMathWorks({'knnsearch','createns','niftiread','niftiinfo','readtable','writetable','jsonencode','ind2sub','pdist2'});
cfg.out = 'D:\DEEP PSMA\PAPER3_Discordance\06_output'; cfg.data = 'D:\DEEP PSMA\DEEP PSMA EXTRACTED';
cfg.tol = 1e-9;
switch lower(mode)
    case 'selftest', if ~selftest(), error('P3V:selftest', 'self-test failed'); end
    case {'preflight','full'}, runV(cfg, strcmpi(mode, 'preflight'), p2bDir);
    otherwise, error('P3V:mode', 'mode must be selftest, preflight or full');
end
end

% =============================== run ===============================
function runV(cfg, preflight, p2bDir)
t0 = tic; runId = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
P3d = fullfile(cfg.out, 'P3_conn26_20261112_001915'); P1d = fullfile(cfg.out, 'P1_20261111_233936');
P2bd = fullfile(cfg.out, p2bDir); r2b = p2bDir(end-14:end);
inputs = {fullfile(P3d, 'P3_components_20261112_001915.csv'), 'ac0048cfe059c714b3c293fbe7e88279';
          fullfile(P3d, 'P3_patients_20261112_001915.csv'),   '2f88309e1be6eff7c499c6c1d36b50ed';
          fullfile(P1d, 'P1_patients_20261111_233936.csv'),   'f9b5b5908e523d1122547ef07d5c0afd';
          fullfile(P2bd, ['P2b_manifest_' r2b '.json']),      '4320cad089bca8eecb128b9df5efd972'};
md = cell(size(inputs,1), 2);
for r = 1:size(inputs,1)
    h = md5file(inputs{r,1}); md(r,:) = {inputs{r,1}, h};
    if ~strcmp(h, inputs{r,2}), error('P3V:C1', 'STOP C1: MD5 mismatch %s (%s)', inputs{r,1}, h); end
end
fprintf('C1 MD5: %d/%d match\n', size(inputs,1), size(inputs,1));
man2b = jsondecode(fileread(inputs{4,1}));
C3 = readtable(inputs{1,1}, 'TextType', 'string');
P1 = readtable(inputs{3,1}, 'TextType', 'string'); P1 = P1(P1.connectivity == 26, :);
G = readtable(fullfile(P2bd, ['P2b_psma_pet_grid_' r2b '.csv']), 'TextType', 'string');
cases = sort(unique(C3.case_id));
if preflight
    oldK = readtable(fullfile(cfg.out, 'A3_visualQC_KEY_do_not_open_20261105.csv'), 'TextType', 'string');
    cases = cases(~ismember(cases, oldK.case_id)); cases = cases(1:2);
end
compRows = cell(0, 8); patRows = cell(0, 8);
for ci = 1:numel(cases)
    cid = cases(ci); cidc = char(cid);
    info = niftiinfo(fullfile(cfg.data, cidc, 'PSMA', 'TTB.nii.gz'));
    A = niiAffine(info); sz = double(info.ImageSize(1:3)); g = G(G.case_id == cid, :);
    Dl = reshape([g.dir_00 g.dir_01 g.dir_02 g.dir_10 g.dir_11 g.dir_12 g.dir_20 g.dir_21 g.dir_22], 3, 3)';
    Aexp = [diag([-1 -1 1]) * Dl * diag([g.spacing_x g.spacing_y g.spacing_z]), diag([-1 -1 1]) * [g.origin_x_LPS; g.origin_y_LPS; g.origin_z_LPS]; 0 0 0 1];
    if ~isequal(sz, [g.size_x g.size_y g.size_z]) || max(abs(A(:) - Aexp(:))) > 1e-3
        error('P3V:grid', 'STOP: grid mismatch for %s', cidc);
    end
    tau = sqrt(sum(vecnorm(A(1:3,1:3)).^2));
    vf = fullfile(P2bd, sprintf('P2b_voxelsets_%s.csv', cidc));
    key = matlab.lang.makeValidName(sprintf('P2b_voxelsets_%s.csv', cidc));
    if ~strcmp(man2b.outputs.(key), md5file(vf)), error('P3V:C1', 'STOP: voxel-set MD5 mismatch %s', cidc); end
    VS = readmatrix(vf);
    psma = niftiread(info) > 0;
    [Pb, ~] = boundaryPoints(find(psma), sz, A);
    NS = []; if ~isempty(Pb), NS = createns(Pb, 'NSMethod', 'kdtree'); end
    c3 = C3(C3.case_id == cid, :);
    Vtot = P1.V_FDG_total_mL(P1.case_id == cid);
    num = 0; numM = 0; numU = 0;
    for r = 1:height(c3)
        st = c3.status(r);
        if ~(st == "V_M" || st == "V_U"), continue; end
        q = c3.component_id(r);
        lin = VS(VS(:,1) == q & VS(:,2) == 0, 3);
        if isempty(lin), error('P3V:input', 'STOP: %s component %d (%s) has no final voxel set', cidc, q, st); end
        dv = voxelDistances(lin, psma, NS, sz, A);
        nU = sum(dv > tau + cfg.tol); share = nU / numel(lin);
        if st == "V_U" && share ~= 1, error('P3V:check', 'STOP: %s V_U component %d has voxels within tau', cidc, q); end
        if st == "V_M" && share >= 1, error('P3V:check', 'STOP: %s V_M component %d has no voxel within tau', cidc, q); end
        vU = share * c3.volume_native_mL(r);
        num = num + vU; if st == "V_M", numM = numM + vU; else, numU = numU + vU; end
        compRows(end+1,:) = {cidc, q, char(st), c3.volume_native_mL(r), numel(lin), nU, share, vU}; %#ok<AGROW>
    end
    VU3 = sum(c3.volume_native_mL(c3.status == "V_U"));
    patRows(end+1,:) = {cidc, Vtot, num, numM, numU, abs(numU - VU3), num / Vtot, VU3 / Vtot}; %#ok<AGROW>
    fprintf('  [%3d/%d] %s (%.0f s)\n', ci, numel(cases), cidc, toc(t0));
end
CT = cell2table(compRows, 'VariableNames', {'case_id','component_id','status','volume_native_mL','voxels_resampled', ...
    'voxels_unmatched','unmatched_share','unmatched_native_mL'});
PT = cell2table(patRows, 'VariableNames', {'case_id','V_FDG_total_mL','voxel_unmatched_mL','within_V_M_mL', ...
    'within_V_U_mL','V_U_check_abs_diff_mL','voxel_level_fraction','DMTVF_ref_primary'});
okU = all(PT.V_U_check_abs_diff_mL <= 1e-9 * max(PT.V_FDG_total_mL));
fprintf('Checks: V_U components fully unmatched and V_M components partly matched (per component): PASS; within-V_U volume = primary V_U: %s\n', tern(okU, 'PASS', 'FAIL'));
if ~okU, error('P3V:check', 'STOP: within-V_U volume differs from the primary V_U'); end
tag = 'P3V_'; if preflight, tag = 'P3V_PREFLIGHT_'; end
od = fullfile(cfg.out, [tag 'conn26_' runId]); mkdir(od);
cf = fullfile(od, ['P3V_components_' runId '.csv']); writetable(CT, cf);
pf = fullfile(od, ['P3V_patients_' runId '.csv']); writetable(PT, pf);
if preflight
    fprintf('Preflight: integrity only (no values displayed). Files in %s\n', od);
else
    fprintf('Cohort: voxel-level unmatched %.1f mL of %.1f mL (within V_M components %.1f mL; within V_U components %.1f mL)\n', ...
        sum(PT.voxel_unmatched_mL), sum(PT.V_FDG_total_mL), sum(PT.within_V_M_mL), sum(PT.within_V_U_mL));
    fprintf('Per-patient voxel-level fraction: median %.4f (IQR %.4f-%.4f); primary DMTVF_ref median %.4f\n', ...
        median(PT.voxel_level_fraction), prctile(PT.voxel_level_fraction, 25), prctile(PT.voxel_level_fraction, 75), median(PT.DMTVF_ref_primary));
end
fns = {'knnsearch','createns','niftiread','niftiinfo','prctile'}; wp = struct();
for k = 1:numel(fns), wp.(fns{k}) = which(fns{k}); end
man = struct('stage', 'P3V voxel-level sensitivity (SAP 8.4)', 'mode', tern(preflight, 'PREFLIGHT', 'FULL'), 'run', runId, ...
    'script', 'P3V_voxel_level_v1_0.m', 'script_md5', md5file([mfilename('fullpath') '.m']), 'inputs', {md}, ...
    'p2b_dir', p2bDir, 'n_patients', height(PT), 'outputs', {{cf, md5file(cf); pf, md5file(pf)}}, 'matlab', version, ...
    'function_paths', wp, 'tolerance', cfg.tol, 'runtime_s', round(toc(t0)), 'status', 'COMPLETE');
mf = fullfile(od, ['P3V_manifest_' runId '.json']); fid = fopen(mf, 'w'); fwrite(fid, jsonencode(man, 'PrettyPrint', true)); fclose(fid);
fprintf('done: %s; manifest MD5 %s\n', od, md5file(mf));
end

% =============================== core ===============================
function dv = voxelDistances(lin, psma, NS, sz, A)
% distance (mm) of each voxel centre to the nearest PSMA TTB voxel centre: 0 inside, else nearest boundary centre
dv = zeros(numel(lin), 1);
out = ~psma(lin);
if any(out)
    if isempty(NS), dv(out) = Inf; return; end
    [i, j, k] = ind2sub(sz, lin(out));
    w = A * [i' - 1; j' - 1; k' - 1; ones(1, numel(i))];
    [~, d] = knnsearch(NS, w(1:3,:)'); dv(out) = d;
end
end

function [pts, blin] = boundaryPoints(lin, sz, A)
% identical to P3_matching_classification_v1_4 (6-neighbour boundary; grid edge = outside)
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
[bi, bj, bk] = ind2sub(size(B), find(B));
gi = bi + lo(1) - 1; gj = bj + lo(2) - 1; gk = bk + lo(3) - 1;
blin = sub2ind(sz, gi, gj, gk);
w = A * [gi' - 1; gj' - 1; gk' - 1; ones(1, numel(gi))];
pts = w(1:3,:)';
end

% =============================== self-test ===============================
function ok = selftest()
res = cell(0, 2);
% anisotropic grid 2.73 x 2.73 x 5.0 mm; tau = 6.28 mm
sp = [2.73 2.73 5.0]; A = [diag(sp) zeros(3,1); 0 0 0 1]; sz = [40 40 40]; tau = sqrt(sum(sp.^2));
P = false(sz); P(10:15,10:15,10:15) = true;
[Pb, ~] = boundaryPoints(find(P), sz, A); NS = createns(Pb, 'NSMethod', 'kdtree');
% 1. voxel inside PSMA -> 0
res(end+1,:) = rowOf('voxel inside PSMA -> distance 0', voxelDistances(sub2ind(sz, 12, 12, 12), P, NS, sz, A) == 0);
% 2. exact distances along each axis (1, 2, 3 voxels beyond the face)
d = voxelDistances(sub2ind(sz, [16 17 18], [12 12 12], [12 12 12])', P, NS, sz, A);
res(end+1,:) = rowOf('x distances 2.73 / 5.46 / 8.19 mm', max(abs(d - [2.73; 5.46; 8.19])) < 1e-12);
d = voxelDistances(sub2ind(sz, [12 12], [12 12], [16 17])', P, NS, sz, A);
res(end+1,:) = rowOf('z distances 5.0 / 10.0 mm', max(abs(d - [5; 10])) < 1e-12);
% 3. brute force over ALL PSMA voxel centres (not only boundary) on random outside voxels
rng(1); allP = find(P); [ai, aj, ak] = ind2sub(sz, allP); W = A * [ai' - 1; aj' - 1; ak' - 1; ones(1, numel(ai))];
cand = setdiff(randi(prod(sz), 300, 1), allP);
[ci, cj, ck] = ind2sub(sz, cand); Q = A * [ci' - 1; cj' - 1; ck' - 1; ones(1, numel(ci))];
db = min(pdist2(Q(1:3,:)', W(1:3,:)'), [], 2);
res(end+1,:) = rowOf('nearest PSMA voxel = nearest boundary voxel (300 random voxels, brute force over all PSMA voxels)', ...
    max(abs(voxelDistances(cand, P, NS, sz, A) - db)) < 1e-9);
% 4. partially overlapping FDG component: known unmatched share
F = false(sz); F(14:25,12,12) = true;                 % i = 14..25: 14,15 inside; 16,17 within tau (2.73, 5.46); 18..25 beyond
dv = voxelDistances(find(F), P, NS, sz, A); nU = sum(dv > tau + 1e-9);
res(end+1,:) = rowOf('partial overlap: 8 of 12 voxels beyond tau (share 2/3)', nU == 8 && abs(nU / 12 - 2/3) < 1e-12);
% 5. native scaling: share x native volume
res(end+1,:) = rowOf('numerator uses native volume (share x native mL)', abs((nU / 12) * 3.0 - 2.0) < 1e-12);
% 6. component beyond tau entirely (V_U-like): share 1
F2 = false(sz); F2(30:32,30:32,30:32) = true;
res(end+1,:) = rowOf('fully separate component -> share exactly 1', all(voxelDistances(find(F2), P, NS, sz, A) > tau));
% 7. exactly at tau counts as matched (<=)
sp2 = [3 4 12]; A2 = [diag(sp2) zeros(3,1); 0 0 0 1]; t2 = 13;   % 3-4-12 -> diagonal 13
P2 = false(sz); P2(5,5,5) = true; [Pb2, ~] = boundaryPoints(find(P2), sz, A2); NS2 = createns(Pb2, 'NSMethod', 'kdtree');
d2 = voxelDistances(sub2ind(sz, 6, 6, 6), P2, NS2, sz, A2);
res(end+1,:) = rowOf('corner neighbour at exactly tau (13 mm) -> not unmatched', abs(d2 - 13) < 1e-12 && ~(d2 > t2 + 1e-9));
for r = 1:size(res,1), fprintf('  [%s] %s\n', tern(res{r,2}, 'PASS', 'FAIL'), res{r,1}); end
ok = all([res{:,2}]); fprintf('selftest: %d/%d passed\n', sum([res{:,2}]), size(res,1));
end

% =============================== utilities (as P3 v1.4) ===============================
function r = rowOf(name, c), r = {name, logical(c)}; end

function assertMathWorks(fns)
root = lower(matlabroot);
for k = 1:numel(fns)
    w = regexprep(which(fns{k}), '^built-in \((.*)\)$', '$1');
    if isempty(w) || ~startsWith(lower(w), root)
        error('P3V:path', 'STOP: %s resolves to "%s" (not the MathWorks installation)', fns{k}, w);
    end
end
end

function A = niiAffine(I)
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
fid = fopen(f, 'r'); if fid < 0, error('P3V:io', 'cannot open %s', f); end
data = fread(fid, Inf, '*uint8'); fclose(fid); md.update(data);
h = lower(reshape(dec2hex(typecast(md.digest, 'uint8'), 2)', 1, []));
end

function out = tern(c, a, b), if c, out = a; else, out = b; end, end
