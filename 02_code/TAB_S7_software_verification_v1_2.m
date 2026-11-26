function TAB_S7_software_verification_v1_2
% TAB_S7_software_verification_v1_2  Online Resource 1, Table S7: software verification for every
% analysis run reported in the manuscript, MATLAB R2025a. Base MATLAB only. Read-only on all inputs.
%
% Version 1.2, 25 November 2026. Change from v1.1: the anatomical figure moved to the main text as Fig. 4
% (FIG4_anatomy_v1_0.m, folder FIG4_*); the Fig. S1 run is kept in the table as superseded.
%
% Version 1.1, 25 November 2026. Change from v1.0: the recorded script MD5 is also read from the
% field 'script_MD5' (upper case), used in the two registration manifests (Audit1C_e, Amendment 3);
% v1.0 read only 'script_md5' and reported those two runs as 'not recorded'. Nothing else changed.
%
% Version 1.0, 25 November 2026. Approved by the first author (M.A. Khan): Table S7 compiled in MATLAB.
%
% For each reported run (curated list below) the script
%   - reads the run manifest (JSON) and takes the recorded script name, script MD5 and, where present,
%     the reproduction gate and software versions;
%   - recomputes the MD5 of that script as it exists now in 02_code and compares it with the record;
%   - follows upstream P2/P2b runs referenced anywhere in the manifest and verifies them the same way;
%   - records the MD5 of the manifest and of the dated completion note(s) for the analysis.
% Output: new folder 06_output\TAB_S7_<yyyyMMdd_HHmmss>\ with TableS7_software_verification.csv,
% TableS7_notes.csv and TAB_S7_manifest.json. Nothing existing is modified.

root = fileparts(fileparts(mfilename('fullpath')));
out6 = fullfile(root, '06_output'); prot = fullfile(root, '00_protocol'); code = fullfile(root, '02_code');

% ---------- curated list: analysis label, run folders or manifest files (in 06_output), notes ----------
L = {
 'Registration: trunk translation correction (T_A)', {'Audit1C_e_manifest_20261031_141128.json'}, ...
     {'00_protocol/PAPER3_PSMA_FDG_discordance_protocol_v2.10.md'}
 'Registration: regional head and pelvic refinement', {'Amendment3_manifest_20261102_073040.json'}, ...
     {'00_protocol/Amendment3_regional_rigid_refinement_spec_LOCKED_v2_1_20261101.md'}
 'Registration QC: head/pelvis intra-reader agreement', {'A3_kappa_20261125_125041'}, ...
     {'06_output/A3_intrareader_kappa_note_20261125.md', '06_output/A3_visualQC_second_reader_note_20261109.md'}
 'Native FDG components (P1)', {'P1_20261111_233936'}, {'06_output/P1_completion_note_20261111.md'}
 'Primary analysis: mapping, matching, classification', {'P3_conn26_20261112_001915'}, ...
     {'06_output/P2_completion_note_20261112.md', '06_output/P2_completion_note_addendum_20261112.md', ...
      '06_output/P3_completion_note_20261112.md'}
 'Voxel-level analysis', {'P3V_conn26_20261112_002845'}, ...
     {'06_output/P3V_completion_note_20261112.md', '06_output/P3V_completion_note_addendum_20261112.md'}
 'Matching rules (reference, strict, permissive)', {'P3T_ref_conn26_20261112_003908', ...
     'P3T_strict_conn26_20261112_003925', 'P3T_perm_conn26_20261112_003954'}, ...
     {'06_output/P3T_tau_sensitivity_completion_note_20261112.md'}
 'Matching rule: overlap only', {'P3T_ovl_conn26_20261112_073943'}, {'06_output/P3T_OVL_completion_note_20261112.md'}
 'Connectivity (6, 18)', {'P2_conn6_20261112_004502', 'P2_conn18_20261112_004353', 'P3C_conn6_20261112_004741', ...
     'P3C_conn18_20261112_004711'}, {'06_output/P3C_connectivity_sensitivity_completion_note_20261112.md'}
 'Eligibility rule B; registration-QC handling S1, S2, S4', {'P3S_EB_20261112_005937', 'P3S_S1_20261112_010017', ...
     'P3S_S2_20261112_010111', 'P3S_S4_20261112_010156'}, {'06_output/P3S_EB_S1_S2_S4_completion_note_20261112.md'}
 'Registration-QC handling S3b, S3c', {'P3S_S3b_20261112_060859', 'P3S_S3c_20261112_060928'}, ...
     {'06_output/P3S_S3b_S3c_completion_note_20261112.md'}
 'Visceral organs (ST2)', {'ST2a_20261112_072150', 'P3S_ST2_20261112_072438'}, ...
     {'06_output/P3S_ST2_completion_note_20261112.md', '06_output/P3S_ST2_completion_note_addendum_20261112.md', ...
      '06_output/ST2_overlap_audit_note_20261112.md'}
 'Organ displacement (O1, ST1)', {'O1_manifest_20261111_232613.json', 'ST1_20261112_080558'}, ...
     {'06_output/O1_completion_note_20261111.md', '06_output/ST1_completion_note_20261112.md', ...
      '06_output/ST1_completion_note_addendum_20261112.md'}
 'Geometric stress test (perturbations)', {'PERT_perturbations_20261112_082222', 'P3P_ZERO_20261112_090859', ...
     'P3P_L1_20261112_101030', 'P3P_L2_20261112_103235'}, ...
     {'06_output/P3P_perturbation_completion_note_20261112.md', ...
      '06_output/P3P_perturbation_completion_note_addendum_20261112.md', ...
      '00_protocol/StepSAP_perturbation_imagecheck_note_20261112.md'}
 'PSMA uptake at unmatched locations', {'PU_20261112_115324'}, {'06_output/PU_psma_uptake_completion_note_20261112.md'}
 'Stability x uptake cross-tabulation (post-hoc)', {'PX_20261112_121119'}, {'06_output/PX_crosstab_completion_note_20261112.md'}
 'Registration-strategy comparator (post-hoc)', {'P3S_SUPPLIED_20261112_122924'}, ...
     {'00_protocol/StepSAP_supplied_registration_comparator_note_20261112.md', ...
      '06_output/P3S_supplied_registration_completion_note_20261112.md'}
 'Anatomical distribution and correspondence flags', {'AD_20261112_125432'}, ...
     {'06_output/AD_anatomical_distribution_completion_note_20261112.md'}
 'Figures 1-3', {'FIG_20261125_123846'}, {}
 'Figure 4 (anatomical distribution)', {latestDir(out6, 'FIG4_*')}, {}
 'Superseded: Fig. S1 (same figure, before the move to the main text)', {latestDir(out6, 'FIG_S1_*')}, {}
 };

% ---------- verify runs ----------
R = table('Size', [0 10], 'VariableTypes', repmat("string", 1, 10), 'VariableNames', ...
    {'analysis', 'run', 'role', 'manifest_md5', 'script', 'script_md5_recorded', 'script_md5_now', ...
     'script_match', 'gate', 'software'});
N = table('Size', [0 3], 'VariableTypes', repmat("string", 1, 3), 'VariableNames', {'analysis', 'note', 'note_md5'});
seen = strings(0, 1);
for a = 1:size(L, 1)
    queue = string(L{a, 2}); roles = repmat("reported run", size(queue));
    while ~isempty(queue)
        r = queue(1); role = roles(1); queue(1) = []; roles(1) = [];
        if r == "" || any(seen == L{a, 1} + "|" + r), continue; end
        seen(end+1) = L{a, 1} + "|" + r; %#ok<AGROW>
        [mf, M] = readManifest(out6, r);
        if isempty(mf)
            R(end+1, :) = {string(L{a, 1}), r, role, "manifest not found", "", "", "", "NO", "", ""}; %#ok<AGROW>
            continue
        end
        [scr, recMD5] = scriptInfo(M);
        nowMD5 = "script file not found";
        if scr ~= "" && isfile(fullfile(code, scr)), nowMD5 = string(md5file(fullfile(code, scr))); end
        match = "NO"; if recMD5 ~= "not recorded" && recMD5 == nowMD5, match = "yes"; end
        if recMD5 == "not recorded" && nowMD5 ~= "script file not found", match = "not recorded"; end
        R(end+1, :) = {string(L{a, 1}), r, role, string(md5file(mf)), scr, recMD5, nowMD5, match, gateInfo(M), softInfo(M)}; %#ok<AGROW>
        % upstream P2 / P2b runs referenced anywhere in the manifest
        allUp = regexp(fileread(mf), 'P2b?_[A-Za-z0-9_]*?\d{8}_\d{6}', 'match');
        up = unique(string(allUp)); up(up == r) = [];
        up = up(arrayfun(@(u) isfolder(fullfile(out6, u)), up));
        queue = [queue, up]; roles = [roles, repmat("upstream of " + r, size(up))]; %#ok<AGROW>
    end
    for k = 1:numel(L{a, 3})
        p = fullfile(root, L{a, 3}{k});
        if isfile(p), m5 = string(md5file(p)); else, m5 = "not found"; end
        N(end+1, :) = {string(L{a, 1}), string(L{a, 3}{k}), m5}; %#ok<AGROW>
    end
end

nRuns = height(R); nMatch = sum(R.script_match == "yes"); nNR = sum(R.script_match == "not recorded");
nGate = sum(R.gate ~= "");
fprintf('Runs verified: %d; script MD5 match: %d; not recorded in manifest: %d; mismatch or missing: %d\n', ...
    nRuns, nMatch, nNR, nRuns - nMatch - nNR);
fprintf('Runs with a recorded reproduction gate: %d\n', nGate);
fprintf('Completion notes found: %d of %d\n', sum(N.note_md5 ~= "not found"), height(N));
disp(R(R.script_match ~= "yes", {'analysis', 'run', 'script', 'script_match'}));

% ---------- outputs ----------
outDir = fullfile(out6, ['TAB_S7_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
assert(~exist(outDir, 'dir')); mkdir(outDir);
f1 = fullfile(outDir, 'TableS7_software_verification.csv'); writetable(R, f1);
f2 = fullfile(outDir, 'TableS7_notes.csv');                 writetable(N, f2);
S.script = [mfilename '.m']; S.script_md5 = md5file([mfilename('fullpath') '.m']);
S.matlab_version = version; S.run_time = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
S.summary = struct('runs', nRuns, 'script_md5_match', nMatch, 'not_recorded', nNR, ...
    'mismatch_or_missing', nRuns - nMatch - nNR, 'runs_with_gate', nGate, ...
    'notes_found', sum(N.note_md5 ~= "not found"), 'notes_listed', height(N));
S.outputs = struct('file', {f1, f2}, 'md5', {md5file(f1), md5file(f2)});
fid = fopen(fullfile(outDir, 'TAB_S7_manifest.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(S, 'PrettyPrint', true)); fclose(fid);
fprintf('Table S7 written to %s\n', outDir);
end

% =====================================================================================
function [mf, M] = readManifest(out6, r)
mf = ''; M = [];
p = fullfile(out6, r);
if endsWith(r, '.json')
    if isfile(p), mf = p; end
elseif isfolder(p)
    d = dir(fullfile(p, '*manifest*.json'));
    if ~isempty(d), mf = fullfile(p, d(1).name); end
end
if ~isempty(mf), M = jsondecode(fileread(mf)); end
end

function [scr, md] = scriptInfo(M)
scr = ""; md = "not recorded";
if isfield(M, 'script')
    [~, nm, ex] = fileparts(char(string(M.script))); scr = string([nm ex]);
end
if isfield(M, 'script_md5') && ~isempty(M.script_md5), md = string(M.script_md5);
elseif isfield(M, 'script_MD5') && ~isempty(M.script_MD5), md = string(M.script_MD5); end
end

function g = gateInfo(M)
g = "";
for f = ["gate", "fallback_gate"]
    if isfield(M, f)
        v = M.(f); if ~ischar(v) && ~isstring(v), v = jsonencode(v); end
        g = strtrim(g + " " + f + ": " + string(v));
    end
end
end

function s = softInfo(M)
s = "";
for f = ["software", "matlab", "matlab_version"]
    if isfield(M, f)
        v = M.(f); if ~ischar(v) && ~isstring(v), v = jsonencode(v); end
        s = string(v); if strlength(s) > 160, s = extractBefore(s, 161); end; return
    end
end
end

function d = latestDir(base, pat)
x = dir(fullfile(base, pat)); x = x([x.isdir]);
if isempty(x), d = ''; else, [~, i] = max([x.datenum]); d = x(i).name; end
end

function h = md5file(p)
fid = fopen(p, 'r'); assert(fid > 0, 'Cannot open %s', p);
b = fread(fid, inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(b(:), 'int8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
