%% Audit1B_totalsegmentator_labels.m
% =========================================================================
% PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - AUDIT 1B
% TotalSegmentator labels and anatomical coverage audit (all 100 cases,
% PSMA CT and FDG CT = 200 CT examinations).
% Protocol: PAPER3_PSMA_FDG_discordance_protocol_v2.5 (Section 8 Audit B, 8.1)
%
% READ-ONLY on:   D:\DEEP PSMA\DEEP PSMA EXTRACTED  (totseg_24.nii.gz, CT header)
% Writes ONLY to: D:\DEEP PSMA\PAPER3_Discordance\06_output
% No spatial transformation, no resampling, no PSMA<->FDG mapping. World
% coordinates (mm) come from each image's own NIfTI header only.
%
% Label dictionary: TotalSegmentator v2.4, task = total, 117 classes
%   Sources: DEEP-PSMA dataset page (https://deep-psma.grand-challenge.org/dataset/)
%            states "TotalSegmentator (v2.4, task=total)";
%            label IDs/names: TotalSegmentator README "Class details" /
%            totalsegmentator/map_to_binary.py (github.com/wasserth/TotalSegmentator)
%
% Per CT and per label (117):
%   presence, voxel count, volume (mL), centroid (world mm), craniocaudal
%   extent, contact with each CT image boundary (S/I and in-plane edges,
%   anatomical direction from the affine), connected components (26) and
%   fraction in the largest component.
% Structure categories (boundary contact SUGGESTS possible truncation, it
% does not prove it; absence is NOT interpreted clinically):
%   PRESENT_COMPLETE                  present, no boundary contact
%   PRESENT_BOUNDARY_Z                touches superior/inferior CT boundary
%   PRESENT_BOUNDARY_INPLANE          touches an in-plane CT edge
%   PRESENT_BOUNDARY_Z_AND_INPLANE    both
%   ABSENT_EXPECTED_OUTSIDE_CT        expected position entirely outside CT z-range
%   ABSENT_EXPECTED_AT_CT_BOUNDARY    expected extent overlaps a CT z-boundary
%   ABSENT_UNEXPLAINED                expected inside CT z-range (segmentation
%                                     failure, anatomical absence, in-plane
%                                     coverage or other - not interpreted)
%   ABSENT_POSITION_UNKNOWN           expected position cannot be estimated
% Expected position of an absent label: cohort median normalised
% craniocaudal position between vertebrae L5 and T1 (from complete
% structures), applied to the scan's own L5/T1 centroids.
%
% Independent verification of the dictionary (descriptive):
%   laterality (right-labelled structure lies to patient right: larger
%   NIfTI x, RAS+), vertebral ordering (sacrum..C1 ascending), craniocaudal
%   order of selected organs, robust volume outliers, fragmentation of
%   single-piece structures, labels outside 1..117.
% Paired coverage: completeness recorded separately for PSMA CT and FDG CT;
% "paired candidate" = PRESENT_COMPLETE on both. Final suitability criteria
% for paired registration QC are locked at Step 2 (not here).
%
% Requires: MATLAB R2017b+ with Image Processing Toolbox.
% =========================================================================

clear; clc;

%% ---------------------------- configuration ----------------------------
cfg.dataRoot   = 'D:\DEEP PSMA\DEEP PSMA EXTRACTED';
cfg.outDir     = 'D:\DEEP PSMA\PAPER3_Discordance\06_output';
cfg.tracers    = {'PSMA','FDG'};
cfg.nClasses   = 117;
cfg.geomTolMM  = 1e-3;          % totseg vs CT grid identity tolerance
cfg.refLow     = 27;            % vertebrae_L5 (normalisation reference)
cfg.refHigh    = 43;            % vertebrae_T1 (normalisation reference)
cfg.minScansForExpected = 10;   % min complete scans to estimate expected position
cfg.robustZ    = 3.5;           % robust z (log volume) for volume outliers
cfg.minScansForVolume = 10;
cfg.fragLabels = [1 2 3 5 21 22 25:50 51 75 76 77 78 90];  % expected single-piece
cfg.fragMinLargestFrac = 0.95;
cfg.dictionary = struct('software','TotalSegmentator','version','2.4','task','total', ...
    'n_classes',117, ...
    'source_version','DEEP-PSMA dataset page: https://deep-psma.grand-challenge.org/dataset/ ("TotalSegmentator (v2.4, task=total)")', ...
    'source_labels','TotalSegmentator README Class details; totalsegmentator/map_to_binary.py (https://github.com/wasserth/TotalSegmentator)');

names = {'spleen','kidney_right','kidney_left','gallbladder','liver','stomach','pancreas', ...
 'adrenal_gland_right','adrenal_gland_left','lung_upper_lobe_left','lung_lower_lobe_left', ...
 'lung_upper_lobe_right','lung_middle_lobe_right','lung_lower_lobe_right','esophagus','trachea', ...
 'thyroid_gland','small_bowel','duodenum','colon','urinary_bladder','prostate','kidney_cyst_left', ...
 'kidney_cyst_right','sacrum','vertebrae_S1','vertebrae_L5','vertebrae_L4','vertebrae_L3', ...
 'vertebrae_L2','vertebrae_L1','vertebrae_T12','vertebrae_T11','vertebrae_T10','vertebrae_T9', ...
 'vertebrae_T8','vertebrae_T7','vertebrae_T6','vertebrae_T5','vertebrae_T4','vertebrae_T3', ...
 'vertebrae_T2','vertebrae_T1','vertebrae_C7','vertebrae_C6','vertebrae_C5','vertebrae_C4', ...
 'vertebrae_C3','vertebrae_C2','vertebrae_C1','heart','aorta','pulmonary_vein', ...
 'brachiocephalic_trunk','subclavian_artery_right','subclavian_artery_left', ...
 'common_carotid_artery_right','common_carotid_artery_left','brachiocephalic_vein_left', ...
 'brachiocephalic_vein_right','atrial_appendage_left','superior_vena_cava','inferior_vena_cava', ...
 'portal_vein_and_splenic_vein','iliac_artery_left','iliac_artery_right','iliac_vena_left', ...
 'iliac_vena_right','humerus_left','humerus_right','scapula_left','scapula_right', ...
 'clavicula_left','clavicula_right','femur_left','femur_right','hip_left','hip_right', ...
 'spinal_cord','gluteus_maximus_left','gluteus_maximus_right','gluteus_medius_left', ...
 'gluteus_medius_right','gluteus_minimus_left','gluteus_minimus_right','autochthon_left', ...
 'autochthon_right','iliopsoas_left','iliopsoas_right','brain','skull', ...
 'rib_left_1','rib_left_2','rib_left_3','rib_left_4','rib_left_5','rib_left_6','rib_left_7', ...
 'rib_left_8','rib_left_9','rib_left_10','rib_left_11','rib_left_12', ...
 'rib_right_1','rib_right_2','rib_right_3','rib_right_4','rib_right_5','rib_right_6', ...
 'rib_right_7','rib_right_8','rib_right_9','rib_right_10','rib_right_11','rib_right_12', ...
 'sternum','costal_cartilages'};
assert(numel(names) == cfg.nClasses, 'Dictionary must contain 117 names.');

grp = repmat({''}, 1, 117);
grp([1 4 5 6 7 8 9 2 3 23 24]) = {'abdominal_organ'};
grp(10:14) = {'lung'}; grp([15 16 17]) = {'neck_thorax_organ'};
grp([18 19 20]) = {'bowel'}; grp([21 22]) = {'pelvic_organ'};
grp(25:50) = {'spine'}; grp(51:68) = {'heart_vessel'};
grp([69:78 116 117]) = {'bone_other'}; grp(92:115) = {'rib'};
grp(79) = {'spinal_cord'}; grp(80:89) = {'muscle'}; grp([90 91]) = {'head'};

% laterality checks: [label expected more to patient RIGHT, label expected more LEFT]
latPairs = [2 3; 8 9; 12 10; 14 11; 55 56; 57 58; 60 59; 66 65; 68 67; 70 69; 72 71; ...
            74 73; 76 75; 78 77; 81 80; 83 82; 85 84; 87 86; 89 88; (104:115)' (92:103)'; ...
            5 1; 5 6; 63 52];
% craniocaudal checks: [label expected SUPERIOR, label expected INFERIOR]
ccPairs = [90 51; 17 51; 51 5; 5 21; 2 21; 3 21; 21 22; 77 75; 78 76];
spineSeq = 25:50;   % sacrum, S1, L5 ... C1 (centroid z must ascend)

%% ------------------------------ set-up ---------------------------------
assert(exist(cfg.dataRoot,'dir')==7, 'Data root not found: %s', cfg.dataRoot);
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
ts = datestr(now,'yyyymmdd_HHMMSS');
fn = @(name,ext) fullfile(cfg.outDir, sprintf('Audit1B_%s_%s.%s', name, ts, ext));
logFile = fn('run_log','txt'); diary(logFile); diary on; tStart = tic;
fprintf('Audit 1B started %s\nData root: %s\nOutput: %s\n\n', datestr(now), cfg.dataRoot, cfg.outDir);

d = dir(fullfile(cfg.dataRoot,'train_*')); d = d([d.isdir]); caseIDs = sort({d.name});
nC = numel(caseIDs); nT = numel(cfg.tracers); nS = nC*nT; K = cfg.nClasses;
fprintf('Case folders: %d | CT examinations: %d\n\n', nC, nS);

CNT = zeros(nS,K); VOL = nan(nS,K); CX = nan(nS,K); CY = nan(nS,K); CZ = nan(nS,K);
ZLO = nan(nS,K); ZHI = nan(nS,K); CONT = zeros(nS,K,6); NCOMP = nan(nS,K); FRACL = nan(nS,K);
sc = struct('case_id',{},'tracer',{},'CT_size',{},'CT_spacing_mm',{},'CT_z_min_mm',{}, ...
    'CT_z_max_mm',{},'totseg_on_CT_grid',{},'unexpected_labels',{},'edge_dirs',{},'read_error',{});

%% ------------------------------ pass 1 ---------------------------------
s = 0;
for c = 1:nC
    for t = 1:nT
        s = s + 1; cid = caseIDs{c}; tr = cfg.tracers{t};
        dTr = fullfile(cfg.dataRoot, cid, tr);
        fprintf('[%3d/%3d] %s %s\n', s, nS, cid, tr);
        sc(s).case_id = cid; sc(s).tracer = tr; sc(s).read_error = '';
        sc(s).unexpected_labels = ''; sc(s).edge_dirs = '';
        sc(s).CT_size = ''; sc(s).CT_spacing_mm = ''; sc(s).CT_z_min_mm = NaN; sc(s).CT_z_max_mm = NaN;
        sc(s).totseg_on_CT_grid = false;
        try
            tsI = niftiinfo(fullfile(dTr,'totseg_24.nii.gz'));
            ctI = niftiinfo(fullfile(dTr,'CT.nii.gz'));       % header only
            sc(s).totseg_on_CT_grid = sameGrid(tsI, ctI, cfg.geomTolMM);
            A  = niiAffine(tsI); sz = double(tsI.ImageSize(1:3));
            vv = abs(det(A(1:3,1:3))) / 1000;                 % mL per voxel
            sc(s).CT_size = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false), 'x');
            sc(s).CT_spacing_mm = strjoin(arrayfun(@(x) sprintf('%.4f',x), vecnorm(A(1:3,1:3)), 'UniformOutput', false), 'x');
            [bMn, bMx] = bbox(A, sz); sc(s).CT_z_min_mm = bMn(3); sc(s).CT_z_max_mm = bMx(3);
            eNames = {edgeLabel(A,1,-1), edgeLabel(A,1,+1), edgeLabel(A,2,-1), edgeLabel(A,2,+1), ...
                      edgeLabel(A,3,-1), edgeLabel(A,3,+1)};
            sc(s).edge_dirs = strjoin(eNames, ',');
            if ~any(strcmp(eNames{5},{'S','I'})) || ~any(strcmp(eNames{6},{'S','I'}))
                sc(s).read_error = 'slice axis is not craniocaudal';
            end

            TS = niftiread(tsI);
            cnt = zeros(K,1); sI = cnt; sJ = cnt; sK = cnt;
            mnI = inf(K,1); mxI = -inf(K,1); mnJ = inf(K,1); mxJ = -inf(K,1);
            mnK = inf(K,1); mxK = -inf(K,1); cont = zeros(K,6); unexp = [];
            M = sz(1) + 1;
            for k = 1:sz(3)
                Lk = TS(:,:,k); nz = find(Lk);
                if isempty(nz), continue; end
                lab = double(Lk(nz));
                bad = lab < 1 | lab > K | lab ~= round(lab);
                if any(bad)
                    unexp = union(unexp, unique(lab(bad))); nz = nz(~bad); lab = lab(~bad);
                    if isempty(nz), continue; end
                end
                [ii, jj] = ind2sub(sz(1:2), nz);
                [labS, ord] = sort(lab);                       % stable sort
                st = find([true; diff(labS) > 0]); en = [st(2:end)-1; numel(labS)]; ul = labS(st);
                ck = zeros(K,1); ck(ul) = en - st + 1;
                jjS = jj(ord);                                 % jj non-decreasing within label
                mnJ(ul) = min(mnJ(ul), jjS(st)); mxJ(ul) = max(mxJ(ul), jjS(en));
                keyS = sort(lab*M + ii); iiS = keyS - labS*M;
                mnI(ul) = min(mnI(ul), iiS(st)); mxI(ul) = max(mxI(ul), iiS(en));
                mnK(ul) = min(mnK(ul), k);      mxK(ul) = max(mxK(ul), k);
                cnt = cnt + ck;
                sI = sI + accumarray(lab, ii, [K 1]); sJ = sJ + accumarray(lab, jj, [K 1]); sK = sK + ck*k;
                cont(:,1) = cont(:,1) + accumarray(lab, double(ii == 1),     [K 1]);
                cont(:,2) = cont(:,2) + accumarray(lab, double(ii == sz(1)), [K 1]);
                cont(:,3) = cont(:,3) + accumarray(lab, double(jj == 1),     [K 1]);
                cont(:,4) = cont(:,4) + accumarray(lab, double(jj == sz(2)), [K 1]);
                if k == 1,     cont(:,5) = cont(:,5) + ck; end
                if k == sz(3), cont(:,6) = cont(:,6) + ck; end
            end
            if ~isempty(unexp), sc(s).unexpected_labels = strjoin(arrayfun(@num2str, unexp(:)', 'UniformOutput', false), ';'); end

            P = find(cnt > 0)';
            CNT(s,:) = cnt'; CONT(s,:,:) = reshape(cont, [1 K 6]);
            for l = P
                ci = [sI(l)/cnt(l) - 1; sJ(l)/cnt(l) - 1; sK(l)/cnt(l) - 1; 1];
                w = A * ci; CX(s,l) = w(1); CY(s,l) = w(2); CZ(s,l) = w(3);
                z1 = A(3,:) * [ci(1); ci(2); mnK(l)-1; 1]; z2 = A(3,:) * [ci(1); ci(2); mxK(l)-1; 1];
                ZLO(s,l) = min(z1,z2); ZHI(s,l) = max(z1,z2);
                VOL(s,l) = cnt(l) * vv;
                sub = TS(mnI(l):mxI(l), mnJ(l):mxJ(l), mnK(l):mxK(l)) == l;
                CCl = bwconncomp(sub, 26); nn = cellfun(@numel, CCl.PixelIdxList);
                NCOMP(s,l) = CCl.NumObjects; FRACL(s,l) = max(nn) / sum(nn);
            end
            clear TS
        catch ME
            sc(s).read_error = ME.message;
            fprintf('   ERROR: %s\n', ME.message);
        end
    end
end

%% ------------------------------ pass 2 ---------------------------------
present = CNT > 0;
contZ   = (CONT(:,:,5) + CONT(:,:,6)) > 0;
contIP  = sum(CONT(:,:,1:4), 3) > 0;
complete = present & ~contZ & ~contIP;
caseOf = reshape(repmat(1:nC, nT, 1), [], 1);  tracerOf = repmat((1:nT)', nC, 1);

% expected normalised craniocaudal position (complete structures, refs complete)
zL = CZ(:,cfg.refLow); zH = CZ(:,cfg.refHigh);
refOK = complete(:,cfg.refLow) & complete(:,cfg.refHigh) & (zH - zL) > 0;
U = (CZ - zL) ./ (zH - zL);
uMed = nan(1,K); hMed = nan(1,K);
for l = 1:K
    m = complete(:,l) & refOK;
    if nnz(m) >= cfg.minScansForExpected
        uMed(l) = median(U(m,l)); hMed(l) = median((ZHI(m,l) - ZLO(m,l)) / 2);
    end
end
refAny = present(:,cfg.refLow) & present(:,cfg.refHigh) & (zH - zL) > 0;

catg = repmat({''}, nS, K);
for s = 1:nS
    for l = 1:K
        if present(s,l)
            if contZ(s,l) && contIP(s,l), catg{s,l} = 'PRESENT_BOUNDARY_Z_AND_INPLANE';
            elseif contZ(s,l),            catg{s,l} = 'PRESENT_BOUNDARY_Z';
            elseif contIP(s,l),           catg{s,l} = 'PRESENT_BOUNDARY_INPLANE';
            else,                         catg{s,l} = 'PRESENT_COMPLETE';
            end
        elseif ~isempty(sc(s).read_error) && isnan(sc(s).CT_z_min_mm)
            catg{s,l} = 'NOT_ASSESSED';
        elseif refAny(s) && ~isnan(uMed(l))
            zE = zL(s) + uMed(l) * (zH(s) - zL(s)); lo = zE - hMed(l); hi = zE + hMed(l);
            if hi < sc(s).CT_z_min_mm || lo > sc(s).CT_z_max_mm, catg{s,l} = 'ABSENT_EXPECTED_OUTSIDE_CT';
            elseif lo < sc(s).CT_z_min_mm || hi > sc(s).CT_z_max_mm, catg{s,l} = 'ABSENT_EXPECTED_AT_CT_BOUNDARY';
            else, catg{s,l} = 'ABSENT_UNEXPLAINED';
            end
        else
            catg{s,l} = 'ABSENT_POSITION_UNKNOWN';
        end
    end
end

% volume outliers (complete structures only; robust z of log volume)
volOut = false(nS,K); volRZ = nan(nS,K);
for l = 1:K
    m = complete(:,l);
    if nnz(m) >= cfg.minScansForVolume
        lv = log(VOL(m,l)); med = median(lv); madv = median(abs(lv - med));
        if madv > 0
            volRZ(:,l) = 0.6745 * (log(VOL(:,l)) - med) / madv;
            volOut(:,l) = m & abs(volRZ(:,l)) > cfg.robustZ;
        end
    end
end
% fragmentation of expected single-piece structures
fragFlag = false(nS,K);
fragFlag(:,cfg.fragLabels) = present(:,cfg.fragLabels) & FRACL(:,cfg.fragLabels) < cfg.fragMinLargestFrac;

%% ------------------------- verification checks -------------------------
abn = {};  latFailLabel = zeros(1,K);
scanRows = struct([]);
for s = 1:nS
    R = rmfield(sc(s), 'edge_dirs'); R.edge_dirs_iMin_iMax_jMin_jMax_kMin_kMax = sc(s).edge_dirs;
    issues = {};
    if ~isempty(sc(s).read_error), issues{end+1} = ['read/geometry: ' sc(s).read_error]; end %#ok<AGROW>
    if ~sc(s).totseg_on_CT_grid, issues{end+1} = 'totseg not on CT grid'; end %#ok<AGROW>
    if ~isempty(sc(s).unexpected_labels), issues{end+1} = ['labels outside 1..117: ' sc(s).unexpected_labels]; end %#ok<AGROW>
    % laterality
    nLat = 0; latF = {};
    for p = 1:size(latPairs,1)
        a = latPairs(p,1); b = latPairs(p,2);
        if present(s,a) && present(s,b)
            nLat = nLat + 1;
            if ~(CX(s,a) > CX(s,b))
                latF{end+1} = sprintf('%s not right of %s', names{a}, names{b}); %#ok<AGROW>
                latFailLabel([a b]) = latFailLabel([a b]) + 1;
            end
        end
    end
    % vertebral order
    pv = spineSeq(present(s,spineSeq)); spF = {};
    for q = 1:numel(pv)-1
        if ~(CZ(s,pv(q+1)) > CZ(s,pv(q))), spF{end+1} = sprintf('%s not above %s', names{pv(q+1)}, names{pv(q)}); end %#ok<AGROW>
    end
    % craniocaudal organ order
    nCC = 0; ccF = {};
    for p = 1:size(ccPairs,1)
        a = ccPairs(p,1); b = ccPairs(p,2);
        if present(s,a) && present(s,b)
            nCC = nCC + 1;
            if ~(CZ(s,a) > CZ(s,b)), ccF{end+1} = sprintf('%s not above %s', names{a}, names{b}); end %#ok<AGROW>
        end
    end
    R.n_labels_present   = nnz(present(s,:));
    R.n_present_complete = nnz(strcmp(catg(s,:),'PRESENT_COMPLETE'));
    R.n_boundary_z       = nnz(present(s,:) & contZ(s,:));
    R.n_boundary_inplane = nnz(present(s,:) & contIP(s,:));
    R.n_absent_outside_CT     = nnz(strcmp(catg(s,:),'ABSENT_EXPECTED_OUTSIDE_CT'));
    R.n_absent_at_CT_boundary = nnz(strcmp(catg(s,:),'ABSENT_EXPECTED_AT_CT_BOUNDARY'));
    R.n_absent_unexplained    = nnz(strcmp(catg(s,:),'ABSENT_UNEXPLAINED'));
    R.absent_unexplained      = strjoin(names(strcmp(catg(s,:),'ABSENT_UNEXPLAINED')), ';');
    R.n_absent_position_unknown = nnz(strcmp(catg(s,:),'ABSENT_POSITION_UNKNOWN'));
    R.boundary_z_structures       = strjoin(names(present(s,:) & contZ(s,:)), ';');
    R.boundary_inplane_structures = strjoin(names(present(s,:) & contIP(s,:)), ';');
    R.laterality_checks = nLat; R.laterality_failures = strjoin(latF, ' | ');
    R.vertebral_order_failures = strjoin(spF, ' | ');
    R.craniocaudal_checks = nCC; R.craniocaudal_failures = strjoin(ccF, ' | ');
    R.volume_outliers = strjoin(names(volOut(s,:)), ';');
    R.fragmented_single_piece = strjoin(names(fragFlag(s,:)), ';');
    R.reference_L5_T1_complete = refOK(s);
    if ~isempty(latF), issues{end+1} = ['laterality: ' R.laterality_failures]; end %#ok<AGROW>
    if ~isempty(spF),  issues{end+1} = ['vertebral order: ' R.vertebral_order_failures]; end %#ok<AGROW>
    if ~isempty(ccF),  issues{end+1} = ['craniocaudal order: ' R.craniocaudal_failures]; end %#ok<AGROW>
    if any(volOut(s,:)),   issues{end+1} = ['volume outlier: ' R.volume_outliers]; end %#ok<AGROW>
    if any(fragFlag(s,:)), issues{end+1} = ['fragmented: ' R.fragmented_single_piece]; end %#ok<AGROW>
    if R.n_absent_unexplained > 0, issues{end+1} = ['absent unexplained: ' R.absent_unexplained]; end %#ok<AGROW>
    if isempty(issues), R.review = 'NO_ISSUE'; else, R.review = 'REVIEW'; end
    R.issues = strjoin(issues, ' || ');
    for q = 1:numel(issues), abn(end+1,:) = {sc(s).case_id, sc(s).tracer, issues{q}}; end %#ok<AGROW>
    if isempty(scanRows), scanRows = R; else, scanRows(end+1) = R; end %#ok<AGROW>
end

%% ------------------------------ tables ---------------------------------
% label level (one row per CT x label)
[lS, lL] = ndgrid(1:nS, 1:K); lS = lS(:); lL = lL(:);
edgeStr = cell(numel(lS),1);
for r = 1:numel(lS)
    e = squeeze(CONT(lS(r), lL(r), :))' > 0;
    if any(e), en_ = strsplit(sc(lS(r)).edge_dirs, ','); edgeStr{r} = strjoin(unique(en_(e)), ';'); else, edgeStr{r} = ''; end
end
LT = table({sc(lS).case_id}', {sc(lS).tracer}', lL, names(lL)', grp(lL)', present(:), catg(:), ...
    CNT(:), VOL(:), CX(:), CY(:), CZ(:), ZLO(:), ZHI(:), ...
    reshape(CONT(:,:,5),[],1), reshape(CONT(:,:,6),[],1), sum(reshape(CONT(:,:,1:4),[],4),2), edgeStr, ...
    NCOMP(:), FRACL(:), volRZ(:), volOut(:), fragFlag(:), ...
    'VariableNames', {'case_id','tracer','label_id','structure','group','present','category', ...
    'voxels','volume_mL','centroid_x_mm','centroid_y_mm','centroid_z_mm','z_low_mm','z_high_mm', ...
    'contact_first_slice_vox','contact_last_slice_vox','contact_inplane_edge_vox','boundary_contact_dirs', ...
    'n_components26','frac_largest_component','volume_robust_z','volume_outlier','fragmented'});

% paired coverage (one row per case x label)
iP = find(tracerOf == 1); iF = find(tracerOf == 2);
[pc, pl] = ndgrid(1:nC, 1:K); pc = pc(:); pl = pl(:);
catP = catg(iP,:); catF = catg(iF,:);
pairCand = complete(iP,:) & complete(iF,:);
PT = table(caseIDs(pc)', pl, names(pl)', grp(pl)', catP(:), catF(:), pairCand(:), ...
    'VariableNames', {'case_id','label_id','structure','group','PSMA_CT_category','FDG_CT_category','paired_candidate'});

% cohort summary per label
cats = {'PRESENT_COMPLETE','PRESENT_BOUNDARY_Z','PRESENT_BOUNDARY_INPLANE','PRESENT_BOUNDARY_Z_AND_INPLANE', ...
        'ABSENT_EXPECTED_OUTSIDE_CT','ABSENT_EXPECTED_AT_CT_BOUNDARY','ABSENT_UNEXPLAINED','ABSENT_POSITION_UNKNOWN'};
CS = table((1:K)', names', grp', 'VariableNames', {'label_id','structure','group'});
CS.n_present_of_200 = sum(present,1)';
CS.n_present_PSMA   = sum(present(iP,:),1)'; CS.n_present_FDG = sum(present(iF,:),1)';
for q = 1:numel(cats), CS.(['n_' cats{q}]) = sum(strcmp(catg, cats{q}),1)'; end
mv = nan(K,1); q1 = nan(K,1); q3 = nan(K,1);
for l = 1:K
    v = sort(VOL(complete(:,l), l));
    if ~isempty(v), mv(l) = median(v); q1(l) = v(max(1,round(0.25*numel(v)))); q3(l) = v(max(1,round(0.75*numel(v)))); end
end
CS.median_volume_mL_complete = mv; CS.Q1_volume_mL = q1; CS.Q3_volume_mL = q3;
CS.n_volume_outliers = sum(volOut,1)'; CS.n_fragmented = sum(fragFlag,1)';
CS.n_laterality_failures = latFailLabel';
CS.expected_norm_position_L5_T1 = uMed'; CS.median_half_extent_mm = hMed';
CS.n_patients_paired_candidate_of_100 = sum(pairCand,1)';

ST = struct2table(scanRows);
if isempty(abn), abn = cell(0,3); end
AT = cell2table(abn, 'VariableNames', {'case_id','tracer','issue'});

%% ------------------------------ outputs --------------------------------
fCS = fn('cohort_label_summary','csv');  writetable(CS, fCS);
fST = fn('scan_level','csv');            writetable(ST, fST);
fLT = fn('label_level','csv');           writetable(LT, fLT);
fPT = fn('paired_coverage','csv');       writetable(PT, fPT);
fAT = fn('abnormalities','csv');         writetable(AT, fAT);
fX  = fn('tables','xlsx');
writetable(CS, fX, 'Sheet', 'cohort_label_summary');
writetable(ST, fX, 'Sheet', 'scan_level');
writetable(AT, fX, 'Sheet', 'abnormalities');
writetable(PT, fX, 'Sheet', 'paired_coverage');

fprintf('\n===== SUMMARY =====\n');
fprintf('CT examinations: %d | read/geometry errors: %d | totseg not on CT grid: %d\n', nS, ...
    sum(~cellfun(@isempty, {sc.read_error})), sum(~[sc.totseg_on_CT_grid]));
fprintf('Scans with labels outside 1..117: %d\n', sum(~cellfun(@isempty, {sc.unexpected_labels})));
fprintf('Labels never present in any CT: %s\n', strjoin(names(CS.n_present_of_200 == 0), ', '));
fprintf('Scans with laterality failures: %d | vertebral order failures: %d | craniocaudal failures: %d\n', ...
    sum(~cellfun(@isempty, ST.laterality_failures)), sum(~cellfun(@isempty, ST.vertebral_order_failures)), ...
    sum(~cellfun(@isempty, ST.craniocaudal_failures)));
fprintf('Scans flagged REVIEW: %d of %d\n', sum(strcmp(ST.review,'REVIEW')), nS);
fprintf('Elapsed: %.1f min\n', toc(tStart)/60);

Mf.audit = 'Audit 1B - TotalSegmentator labels and anatomical coverage';
Mf.protocol = 'PAPER3_PSMA_FDG_discordance_protocol_v2.5';
Mf.script = [mfilename('fullpath') '.m']; Mf.timestamp = ts;
Mf.matlab_version = version; Mf.computer = computer;
v = ver('images'); if ~isempty(v), Mf.image_toolbox_version = v(1).Version; end
Mf.config = cfg; Mf.label_dictionary = names;
Mf.n_case_folders = nC; Mf.n_CT = nS;
Mf.outputs = {fCS, fST, fLT, fPT, fAT, fX, logFile};
Mf.note = ['Read-only. No spatial transformation or resampling. Boundary contact suggests possible ', ...
           'truncation, not proof. Absent structures are not interpreted clinically. Paired-candidate ', ...
           'status is descriptive; registration-QC suitability criteria are locked at Step 2.'];
try, js = jsonencode(Mf, 'PrettyPrint', true); catch, js = jsonencode(Mf); end
fid = fopen(fn('manifest','json'),'w'); fwrite(fid, js, 'char'); fclose(fid);
fprintf('Outputs written to %s (timestamp %s)\n', cfg.outDir, ts);
diary off;

%% ============================ local functions ===========================
function tf = sameGrid(A, B, tol)
    if isempty(A) || isempty(B), tf = false; return; end
    tf = isequal(A.ImageSize(1:3), B.ImageSize(1:3)) && ...
         max(abs(niiAffine(A) - niiAffine(B)), [], 'all') <= tol;
end

function A = niiAffine(I)
    % 4x4 mapping 0-based voxel indices -> NIfTI world (mm); identical to Audit 1A
    r = I.raw;
    if isfield(r,'sform_code') && r.sform_code > 0
        A = [double(r.srow_x(:))'; double(r.srow_y(:))'; double(r.srow_z(:))'; 0 0 0 1];
    else
        b = double(r.quatern_b); c = double(r.quatern_c); d = double(r.quatern_d);
        a = sqrt(max(0, 1 - b^2 - c^2 - d^2));
        Rm = [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c);
              2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b);
              2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-c*c-b*b];
        pd = double(r.pixdim); qfac = pd(1); if qfac == 0, qfac = 1; end
        A = eye(4); A(1:3,1:3) = Rm * diag([pd(2) pd(3) pd(4)*qfac]);
        A(1:3,4) = [double(r.qoffset_x); double(r.qoffset_y); double(r.qoffset_z)];
    end
end

function [mn, mx] = bbox(A, sz)
    [a,b,c] = ndgrid([0 sz(1)-1], [0 sz(2)-1], [0 sz(3)-1]);
    w = A * [a(:)'; b(:)'; c(:)'; ones(1,8)];
    mn = min(w(1:3,:), [], 2); mx = max(w(1:3,:), [], 2);
end

function lab = edgeLabel(A, ax, sgn)
    % anatomical direction of the outward normal of an image edge (NIfTI RAS+:
    % +x = Right, +y = Anterior, +z = Superior)
    v = sgn * A(1:3,ax); [~, d] = max(abs(v));
    L = {'L','R'; 'P','A'; 'I','S'};
    lab = L{d, (v(d) > 0) + 1};
end
