%% Audit1D_lesion_mask_audit.m
% =========================================================================
% PAPER 3 (PSMA-FDG discordance, DEEP-PSMA) - AUDIT 1D
% Tumour-mask (TTB) characteristics audit: all 100 patients, both tracers.
% Protocol: PAPER3_PSMA_FDG_discordance_protocol_v2.6 (Section 8, Audit D)
%
% READ-ONLY on:   D:\DEEP PSMA\DEEP PSMA EXTRACTED (PET, TTB, totseg_24, CT header)
% Writes ONLY to: D:\DEEP PSMA\PAPER3_Discordance\06_output (single folder)
%
% Operating rules
%  - supplied TTB masks never modified; connectivity runs on in-memory copies
%  - all volumes in NATIVE PET space (voxel count x voxel volume from header)
%  - each tracer analysed separately: NO intertracer registration, NO
%    correspondence, NO overlap, NO discordance
%  - thresholds are exploratory feasibility values, NOT eligibility criteria;
%    no final threshold is selected here
%  - a processing failure is never silent: every scan receives a status
%    (OK / PARTIAL / FAILED); the run is declared COMPLETE only if all
%    100 patients and 200 tracer examinations have status OK AND every
%    mandatory output (9 CSV, 1 XLSX, 3 PNG, manifest) was written and is
%    non-empty; otherwise INCOMPLETE
%  - failed scans keep their geometry row (useful for diagnosis), labelled
%    with the scan's final status (column scan_status)
%  - execution test: set cfg.testN = 2 to process the first 2 cases only;
%    such a run is labelled TEST_ONLY and is never an Audit 1D result.
%    Definitive run: cfg.testN = [] (default)
%
% Components
%  1 connectivity 6 / 18 / 26: component counts; merges (6->18, 6->26,
%    18->26) and % of total tumour volume involved in merges
%  2 lesion volumes: counts, volume distributions (mL and voxels), very small
%    components; per-component SUVmax / SUVmean (native PET, same tracer)
%  3 minimum-volume feasibility: 0.1/0.25/0.5/1/2 mL and 2/5/10 voxels:
%    components removed, tumour volume removed (PSMA and FDG separately),
%    patients losing all components
%  4 confluent disease: largest component (volume, fraction of TTB, physical
%    extent); anatomical lookup for large components using PHYSICAL-COORDINATE
%    mapping PET -> same-scan CT/TotalSegmentator grid with nearest-neighbour
%    sampling at PET voxel centres (same examination; no PSMA-FDG
%    registration). A component spanning several labelled structures is an
%    AUTOMATED OBSERVATION; it is not evidence of a single confluent metastasis.
%    Unlabelled tissue is not automatically nodal.
%  5 grid geometry: PET spacing, dimensions, voxel volume, orientation,
%    obliquity, handedness, full PET affine; PET-TTB and CT-totseg grid
%    verification; PSMA/FDG voxel-size compatibility per patient
%
% Geometry safeguards
%  - PET-TTB: identical dimensions and affine (tolerance cfg.gridTol_mm);
%    on failure the scan is NOT processed (status FAILED)
%  - CT-totseg_24: identical dimensions and affine; on failure the anatomical
%    lookup is NOT performed (status PARTIAL)
%  - coordinates are NIfTI world coordinates (RAS+, mm), not SimpleITK LPS
%  - physical extent = world-axis bounding box of the voxel parallelepipeds
%    (correct for oblique grids)
%
% Requires: MATLAB R2018b+ with Image Processing Toolbox.
% =========================================================================

clear; clc;

%% ---------------------------- configuration ----------------------------
cfg.protocol  = 'PAPER3_PSMA_FDG_discordance_protocol_v2.6';
cfg.dataRoot  = 'D:\DEEP PSMA\DEEP PSMA EXTRACTED';
cfg.outDir    = 'D:\DEEP PSMA\PAPER3_Discordance\06_output';
cfg.tracers   = {'PSMA','FDG'};
cfg.files     = struct('PET','PET.nii.gz','TTB','TTB.nii.gz','CT','CT.nii.gz','TS','totseg_24.nii.gz');
cfg.expectedCases = 100;
cfg.testN     = [];                           % [] = DEFINITIVE run (all cases). A number (e.g. 2) = execution test on
                                              % the first N cases only; a test run is never declared COMPLETE.
cfg.gridTol_mm = 1e-3;                        % max abs affine difference for "same grid"
cfg.conns     = [6 18 26];
cfg.thrML     = [0.1 0.25 0.5 1.0 2.0];      % exploratory feasibility only
cfg.thrVox    = [2 5 10];                     % exploratory feasibility only
cfg.smallVox  = [1 4 9];                      % reporting bins: 1, 2-4, 5-9 voxels
cfg.largeML   = 10;                           % anatomical lookup: components >= 10 mL ...
cfg.alsoLargest = true;                       % ... plus the largest component of each scan
cfg.lookupConn = 26;                          % connectivity used for anatomical lookup
cfg.spanMinFrac = 0.05;                       % group counted as spanned if >= 5% of component voxels

% TotalSegmentator v2.4 (task=total, 117 classes) label groups for anatomical lookup
G = struct('name', {}, 'labels', {});
G(end+1) = struct('name','cervical_spine',         'labels', 44:50);
G(end+1) = struct('name','thoracic_spine',         'labels', 32:43);
G(end+1) = struct('name','lumbar_spine',           'labels', 27:31);
G(end+1) = struct('name','sacrum_S1',              'labels', [25 26]);
G(end+1) = struct('name','pelvis_hips',            'labels', [77 78]);
G(end+1) = struct('name','ribs_sternum_cartilage', 'labels', 92:117);
G(end+1) = struct('name','shoulder_girdle_humeri', 'labels', 69:74);
G(end+1) = struct('name','femora',                 'labels', [75 76]);
G(end+1) = struct('name','skull',                  'labels', 91);
G(end+1) = struct('name','prostate',               'labels', 22);
G(end+1) = struct('name','liver',                  'labels', 5);
G(end+1) = struct('name','lungs',                  'labels', 10:14);
G(end+1) = struct('name','heart_vessels',          'labels', 51:68);
G(end+1) = struct('name','muscles_spinal_cord',    'labels', 79:89);
G(end+1) = struct('name','brain',                  'labels', 90);
G(end+1) = struct('name','other_organs',           'labels', [1:4 6:9 15:21 23:24]);
lab2grp = zeros(1,255);
for g = 1:numel(G), lab2grp(G(g).labels) = g; end
assert(all(lab2grp(1:117) > 0) && nnz(lab2grp) == 117, 'Label grouping does not cover labels 1-117 exactly once.');
grpNames = [{G.name}, {'unlabelled_tissue','outside_CT'}];   % unlabelled is NOT automatically nodal
nG = numel(G);

%% ------------------------------ set-up ---------------------------------
assert(exist(cfg.dataRoot,'dir')==7, 'Data root not found.');
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
ts = datestr(now,'yyyymmdd_HHMMSS');
fn = @(name,ext) fullfile(cfg.outDir, sprintf('Audit1D_%s_%s.%s', name, ts, ext));
logFile = fn('run_log','txt'); diary(logFile); diary on; tStart = tic;

scriptPath = [mfilename('fullpath') '.m'];
if exist(scriptPath,'file') ~= 2, scriptPath = ''; end
scriptMD5 = fileMD5(scriptPath);

fprintf('Audit 1D started %s\nProtocol: %s\nScript: %s (MD5 %s)\nData root: %s\nOutput: %s\n\n', ...
    datestr(now), cfg.protocol, scriptPath, scriptMD5, cfg.dataRoot, cfg.outDir);

d = dir(fullfile(cfg.dataRoot,'train_*')); d = d([d.isdir]); caseIDs = sort({d.name});
fprintf('Case folders found: %d (expected %d)\n', numel(caseIDs), cfg.expectedCases);
nCasesFound = numel(caseIDs);
isTest = ~isempty(cfg.testN);
if isTest
    caseIDs = caseIDs(1:min(cfg.testN, numel(caseIDs)));
    fprintf('*** EXECUTION TEST: processing the first %d case(s) only; this run cannot be declared COMPLETE. ***\n', numel(caseIDs));
end
fprintf('\n');

comp = {};  scan = {};  merge = {};  thr = {};  anat = {};  geo = {};  inv = {};  status = {};
compHdr = {'case_id','tracer','connectivity','component_id','voxels','volume_mL','SUVmax','SUVmean', ...
           'centroid_x_RAS_mm','centroid_y_RAS_mm','centroid_z_RAS_mm', ...
           'extent_x_mm','extent_y_mm','extent_z_mm','max_extent_mm'};
pctOf = @(a,b) ternary(b > 0, 100*a/max(b,eps), NaN);

%% ------------------------------ main loop ------------------------------
for c = 1:numel(caseIDs)
    cid = caseIDs{c};
    fprintf('[%3d/%3d] %s\n', c, numel(caseIDs), cid);
    gp = struct(); stTr = struct();
    for t = 1:numel(cfg.tracers)
        tr = cfg.tracers{t}; dTr = fullfile(cfg.dataRoot, cid, tr);
        st = ''; msg = ''; gi = 0;
        lComp = {}; lScan = {}; lThr = {}; lMerge = {}; lAnat = {};
        try
            % ---------- input inventory ----------
            fl    = {fullfile(dTr,cfg.files.PET), fullfile(dTr,cfg.files.TTB), fullfile(dTr,cfg.files.CT), fullfile(dTr,cfg.files.TS)};
            roles = {'PET','TTB','CT_header_only','totseg_24'};
            for k = 1:numel(fl)
                dd = dir(fl{k});
                if isempty(dd)
                    inv(end+1,:) = {cid, tr, roles{k}, fl{k}, NaN, 'MISSING'}; %#ok<SAGROW>
                    error('Missing input file: %s', fl{k});
                end
                inv(end+1,:) = {cid, tr, roles{k}, fl{k}, dd.bytes, datestr(dd.datenum,'yyyy-mm-dd HH:MM:SS')}; %#ok<SAGROW>
            end

            % ---------- headers and grid verification ----------
            petI = niftiinfo(fl{1}); ttbI = niftiinfo(fl{2}); ctI = niftiinfo(fl{3}); tsI = niftiinfo(fl{4});
            Ap = niiAffine(petI); Att = niiAffine(ttbI); Act = niiAffine(ctI); At = niiAffine(tsI);
            sz = double(petI.ImageSize(1:3));
            M  = Ap(1:3,1:3); sp = vecnorm(M); Dn = M ./ sp;
            vv = abs(det(M)) / 1000;                          % mL per PET voxel
            obl = max(acosd(min(1, max(abs(Dn), [], 1))));    % largest deviation of a voxel axis from a world axis
            halfW = 0.5 * sum(abs(M), 2);                     % world-axis half-extent of one voxel (oblique-safe)
            [okTTB, dTTB] = sameGrid(petI, Ap, ttbI, Att, cfg.gridTol_mm);
            [okTS,  dTS ] = sameGrid(ctI,  Act, tsI, At,  cfg.gridTol_mm);
            gp.(tr) = struct('sp', sp, 'vv', vv);
            geo(end+1,:) = {cid, tr, sizeStr(sz), sp(1), sp(2), sp(3), vv*1000, obl, sign(det(M)), axisCode(M), ...
                okTTB, dTTB, okTS, dTS, affStr(Ap), ''}; %#ok<SAGROW>
            gi = size(geo,1);                                 % scan status filled in once processing ends
            if ~okTTB
                error('PET-TTB grid mismatch (PET %s, TTB %s, max affine difference %.3g mm): scan not processed', ...
                    sizeStr(sz), sizeStr(ttbI.ImageSize(1:3)), dTTB);
            end

            % ---------- read (copies in memory; originals untouched) ----------
            PET = double(niftiread(petI));
            s = petI.MultiplicativeScaling; o = petI.AdditiveOffset;
            if s ~= 0 && (s ~= 1 || o ~= 0), PET = PET * s + o; end
            BW = niftiread(ttbI) > 0;
            if ~isequal([size(BW,1) size(BW,2) size(BW,3)], sz) || ~isequal([size(PET,1) size(PET,2) size(PET,3)], sz)
                error('Array size differs from PET header size: scan not processed');
            end
            totVox = nnz(BW); totML = totVox * vv;

            % ---------- 1 + 2 + 3. components per connectivity ----------
            CCs = cell(1, numel(cfg.conns)); Ls = cell(1, numel(cfg.conns));
            for k = 1:numel(cfg.conns)
                cn = cfg.conns(k);
                CC = bwconncomp(BW, cn); CCs{k} = CC; Ls{k} = labelmatrix(CC);
                nObj = CC.NumObjects;
                nv = cellfun(@numel, CC.PixelIdxList); nv = reshape(nv, 1, []);
                ml = nv * vv;
                for q = 1:nObj
                    idx = CC.PixelIdxList{q};
                    [ii,jj,kk] = ind2sub(sz, idx);
                    w = Ap * [ii'-1; jj'-1; kk'-1; ones(1,numel(ii))];
                    ext = max(w(1:3,:),[],2) - min(w(1:3,:),[],2) + 2*halfW;   % world-axis bounding box (mm)
                    pv = PET(idx);
                    lComp(end+1,:) = {cid, tr, cn, q, nv(q), ml(q), max(pv), mean(pv), ...
                        mean(w(1,:)), mean(w(2,:)), mean(w(3,:)), ext(1), ext(2), ext(3), max(ext)}; %#ok<SAGROW>
                end
                % scan-level summary (empty mask -> zero components, NaN statistics; no dummy component)
                if nObj == 0
                    lScan(end+1,:) = {cid, tr, cn, totVox, totML, 0, NaN, NaN, NaN, 0, 0, 0, NaN, NaN}; %#ok<SAGROW>
                else
                    [mx, imx] = max(ml);
                    lScan(end+1,:) = {cid, tr, cn, totVox, totML, nObj, median(ml), mx, mx/totML, ...
                        nnz(nv <= cfg.smallVox(1)), nnz(nv > cfg.smallVox(1) & nv <= cfg.smallVox(2)), ...
                        nnz(nv > cfg.smallVox(2) & nv <= cfg.smallVox(3)), sum(ml(nv <= cfg.smallVox(3)))/totML, imx}; %#ok<SAGROW>
                end
                % threshold feasibility (exploratory)
                for x = cfg.thrML
                    rm = ml < x;
                    lThr(end+1,:) = {cid, tr, cn, 'mL', x, nObj, nnz(rm), pctOf(nnz(rm),nObj), ...
                        sum(ml(rm)), pctOf(sum(ml(rm)),totML), nObj > 0 && all(rm), x/vv}; %#ok<SAGROW>
                end
                for x = cfg.thrVox
                    rm = nv < x;
                    lThr(end+1,:) = {cid, tr, cn, 'voxels', x, nObj, nnz(rm), pctOf(nnz(rm),nObj), ...
                        sum(ml(rm)), pctOf(sum(ml(rm)),totML), nObj > 0 && all(rm), x}; %#ok<SAGROW>
                end
            end

            % ---------- 1. merge mapping ----------
            pairs = [1 2; 1 3; 2 3];   % 6->18, 6->26, 18->26
            for p = 1:size(pairs,1)
                a = pairs(p,1); b = pairs(p,2);
                CCb = CCs{b}; La = Ls{a};
                nMergedComps = 0; mergedVox = 0; maxParts = 0;
                for q = 1:CCb.NumObjects
                    u = unique(La(CCb.PixelIdxList{q}));
                    if numel(u) > 1
                        nMergedComps = nMergedComps + 1; mergedVox = mergedVox + numel(CCb.PixelIdxList{q});
                        maxParts = max(maxParts, numel(u));
                    end
                end
                lMerge(end+1,:) = {cid, tr, cfg.conns(a), cfg.conns(b), CCs{a}.NumObjects, CCb.NumObjects, ...
                    CCs{a}.NumObjects - CCb.NumObjects, nMergedComps, maxParts, mergedVox*vv, ...
                    pctOf(mergedVox, totVox)}; %#ok<SAGROW>
            end

            % ---------- 4. anatomical extent of large components ----------
            kL = find(cfg.conns == cfg.lookupConn, 1);
            CC = CCs{kL}; nv = reshape(cellfun(@numel, CC.PixelIdxList), 1, []); ml = nv * vv;
            if CC.NumObjects > 0
                if ~okTS
                    st = 'PARTIAL';
                    msg = sprintf('totseg_24-CT grid mismatch (max affine difference %.3g mm): anatomical lookup not performed', dTS);
                else
                    sel = find(ml >= cfg.largeML);
                    if cfg.alsoLargest, [~, imx] = max(ml); sel = unique([sel, imx]); end
                    TS = niftiread(tsI); szT = [size(TS,1) size(TS,2) size(TS,3)];
                    T = At \ Ap;                              % PET index (0-based) -> CT index (0-based), physical mapping
                    for q = sel
                        idx = CC.PixelIdxList{q};
                        [ii,jj,kk] = ind2sub(sz, idx);
                        ci = round(T * [ii'-1; jj'-1; kk'-1; ones(1,numel(ii))]) + 1;   % nearest CT voxel at PET voxel centre
                        inside = all(ci(1:3,:) >= 1 & ci(1:3,:) <= szT(:), 1);
                        gidx = zeros(1, numel(idx)); gidx(~inside) = nG + 2;
                        li = sub2ind(szT, ci(1,inside), ci(2,inside), ci(3,inside));
                        lab = double(TS(li));
                        gg = zeros(size(lab)); gg(lab > 0) = lab2grp(lab(lab > 0)); gg(gg == 0) = nG + 1;
                        gidx(inside) = gg;
                        fr = accumarray(gidx(:), 1, [nG+2 1])' / numel(idx);
                        spanned = find(fr(1:nG) >= cfg.spanMinFrac);
                        obs = '';
                        if numel(spanned) >= 2, obs = 'spans_multiple_labelled_structures (automated observation)'; end
                        lAnat(end+1,:) = [{cid, tr, cfg.lookupConn, q, ml(q), ml(q)/totML, numel(spanned), ...
                            strjoin(grpNames(spanned), ';'), fr(nG+1), fr(nG+2), obs}, num2cell(fr(1:nG))]; %#ok<SAGROW>
                    end
                    clear TS
                end
            end
            if isempty(st), st = 'OK'; end
            fprintf('   %s [%s]: TTB %.1f mL | components 6/18/26: %d/%d/%d | largest(26) %.1f mL\n', tr, st, totML, ...
                CCs{1}.NumObjects, CCs{2}.NumObjects, CCs{3}.NumObjects, max([ml 0]));
            if ~isempty(msg), fprintf('   %s WARNING: %s\n', tr, msg); end
            % commit this scan's rows
            comp = [comp; lComp]; scan = [scan; lScan]; thr = [thr; lThr]; merge = [merge; lMerge]; anat = [anat; lAnat]; %#ok<AGROW>
            clear PET BW CCs Ls
        catch ME
            st = 'FAILED'; msg = ME.message;          % no rows from a failed scan are committed
            fprintf('   %s FAILED: %s\n', tr, msg);
        end
        status(end+1,:) = {cid, tr, st, msg}; %#ok<SAGROW>
        if gi > 0, geo{gi, end} = st; end              % geometry row carries the scan's final status
        stTr.(tr) = st;
    end
    % grid compatibility per patient (no resampling)
    if isfield(gp,'PSMA') && isfield(gp,'FDG')
        geo(end+1,:) = {cid, 'PAIR_FDG_over_PSMA', '', gp.FDG.sp(1)/gp.PSMA.sp(1), gp.FDG.sp(2)/gp.PSMA.sp(2), ...
            gp.FDG.sp(3)/gp.PSMA.sp(3), gp.FDG.vv/gp.PSMA.vv, NaN, NaN, ...
            ternary(gp.FDG.vv < gp.PSMA.vv, 'FDG_voxels_smaller_than_PSMA (higher small-lesion resampling risk)', ''), ...
            NaN, NaN, NaN, NaN, '', sprintf('PSMA:%s;FDG:%s', stTr.PSMA, stTr.FDG)}; %#ok<SAGROW>
    end
end

%% ------------------------- data completeness ---------------------------
% (final COMPLETE/INCOMPLETE status is decided only after all mandatory outputs are verified)
SS = mkT(status, {'case_id','tracer','status','message'});
nCases   = numel(caseIDs);
nExpScan = cfg.expectedCases * numel(cfg.tracers);
nOK      = nnz(strcmp(SS.status,'OK'));
nPartial = nnz(strcmp(SS.status,'PARTIAL'));
nFailed  = nnz(strcmp(SS.status,'FAILED'));
dataComplete = ~isTest && (nCasesFound == cfg.expectedCases) && (nCases == cfg.expectedCases) && ...
               (height(SS) == nExpScan) && (nOK == nExpScan);

%% ------------------------------ tables ---------------------------------
CT  = mkT(comp, compHdr);
ST  = mkT(scan, {'case_id','tracer','connectivity','TTB_voxels','TTB_mL','n_components', ...
        'median_component_mL','largest_component_mL','largest_fraction_of_TTB','n_1voxel','n_2to4voxels', ...
        'n_5to9voxels','fraction_TTB_in_components_le9voxels','largest_component_id'});
MT  = mkT(merge, {'case_id','tracer','conn_from','conn_to','n_from','n_to','n_reduction', ...
        'n_merged_components','max_parts_in_one_component','merged_volume_mL','pct_TTB_volume_in_merges'});
TT  = mkT(thr, {'case_id','tracer','connectivity','threshold_unit','threshold','n_components', ...
        'n_removed','pct_components_removed','volume_removed_mL','pct_TTB_volume_removed','all_components_removed', ...
        'threshold_in_voxels'});
anHdr = [{'case_id','tracer','connectivity','component_id','volume_mL','fraction_of_TTB','n_labelled_groups_spanned', ...
          'groups_spanned','frac_unlabelled_tissue','frac_outside_CT','observation'}, strcat('frac_', {G.name})];
AT  = mkT(anat, anHdr);
GT  = mkT(geo, {'case_id','tracer','PET_size','spacing_x_or_ratio','spacing_y_or_ratio','spacing_z_or_ratio', ...
        'voxel_volume_uL_or_ratio','obliquity_deg','handedness','axis_code_or_flag','TTB_on_PET_grid', ...
        'PET_TTB_max_affine_diff_mm','totseg_on_CT_grid','CT_totseg_max_affine_diff_mm','PET_affine_RAS_row_major', ...
        'scan_status'});
IT  = mkT(inv, {'case_id','tracer','role','path','bytes','modified'});

% cohort threshold summary (PSMA and FDG separately; scans with status OK or PARTIAL)
sumRows = {};
for t = 1:numel(cfg.tracers)
    for cn = cfg.conns
        for u = {'mL','voxels'}
            base = strcmp(TT.tracer,cfg.tracers{t}) & TT.connectivity==cn & strcmp(TT.threshold_unit,u{1});
            vals = unique(TT.threshold(base))';
            stv = ST(strcmp(ST.tracer,cfg.tracers{t}) & ST.connectivity==cn, :);
            for x = vals
                m = base & TT.threshold==x;
                sumRows(end+1,:) = {cfg.tracers{t}, cn, u{1}, x, nnz(m), median(TT.pct_components_removed(m),'omitnan'), ...
                    median(TT.pct_TTB_volume_removed(m),'omitnan'), prctile90(TT.pct_TTB_volume_removed(m)), ...
                    pctOf(sum(TT.volume_removed_mL(m)), sum(stv.TTB_mL)), nnz(TT.all_components_removed(m)), ...
                    min(TT.threshold_in_voxels(m)), max(TT.threshold_in_voxels(m))}; %#ok<SAGROW>
            end
        end
    end
end
TS2 = mkT(sumRows, {'tracer','connectivity','threshold_unit','threshold','n_scans', ...
        'median_pct_components_removed','median_pct_TTB_volume_removed','P90_pct_TTB_volume_removed', ...
        'pooled_pct_TTB_volume_removed','n_patients_losing_all_components','min_threshold_in_voxels','max_threshold_in_voxels'});

files = {}; otherErrors = {};
outs = {'processing_status',SS; 'component_level',CT; 'scan_level',ST; 'connectivity_merges',MT; ...
        'threshold_feasibility',TT; 'threshold_cohort_summary',TS2; 'large_component_anatomy',AT; ...
        'grid_geometry',GT; 'input_inventory',IT};
for k = 1:size(outs,1)
    f = fn(outs{k,1},'csv');
    try
        writetable(outs{k,2}, f); files{end+1} = f; %#ok<SAGROW>
    catch ME
        otherErrors{end+1} = sprintf('CSV %s: %s', outs{k,1}, ME.message); %#ok<SAGROW>
        fprintf('CSV export error (%s): %s\n', outs{k,1}, ME.message);
    end
end
fX = fn('tables','xlsx');
try
    for k = [1 3 4 6 7 8]
        writetable(outs{k,2}, fX, 'Sheet', outs{k,1});
    end
    files{end+1} = fX;
catch ME
    otherErrors{end+1} = ['XLSX export: ' ME.message];
    fprintf('XLSX export error: %s\n', ME.message);
end

%% ------------------------------ figures --------------------------------
try
    f1 = figure('Visible','off','Position',[50 50 1300 500],'Color','w');
    for t = 1:2
        subplot(1,2,t);
        v = CT.volume_mL(strcmp(CT.tracer,cfg.tracers{t}) & CT.connectivity==26);
        histogram(log10(v), 50); xlabel('log_{10} component volume (mL)'); ylabel('components');
        title(sprintf('%s components (26-connectivity), n = %d', cfg.tracers{t}, numel(v)));
        hold on; yl = ylim; for x = cfg.thrML, plot(log10([x x]), yl, 'r--'); end; hold off;
    end
    p1 = fn('fig_component_volumes','png'); print(f1, p1, '-dpng', '-r120'); close(f1);
    f2 = figure('Visible','off','Position',[50 50 1300 500],'Color','w');
    for t = 1:2
        subplot(1,2,t);
        s6 = ST(strcmp(ST.tracer,cfg.tracers{t}) & ST.connectivity==6, :);
        s26 = ST(strcmp(ST.tracer,cfg.tracers{t}) & ST.connectivity==26, :);
        [~, ia, ib] = intersect(s6.case_id, s26.case_id);
        scatter(s6.n_components(ia), s26.n_components(ib), 15, 'filled'); hold on;
        mx = max([s6.n_components; 1]); plot([0 mx],[0 mx],'k--'); hold off;
        xlabel('components, 6-connectivity'); ylabel('components, 26-connectivity');
        title(sprintf('%s: component count by connectivity', cfg.tracers{t}));
    end
    p2 = fn('fig_connectivity_counts','png'); print(f2, p2, '-dpng', '-r120'); close(f2);
    f3 = figure('Visible','off','Position',[50 50 1300 500],'Color','w');
    for t = 1:2
        subplot(1,2,t); hold on;
        for cn = cfg.conns
            m = strcmp(TS2.tracer,cfg.tracers{t}) & TS2.connectivity==cn & strcmp(TS2.threshold_unit,'mL');
            plot(TS2.threshold(m), TS2.pooled_pct_TTB_volume_removed(m), '-o', 'DisplayName', sprintf('%d-conn', cn));
        end
        hold off; set(gca,'XScale','log'); grid on; legend('Location','northwest');
        xlabel('minimum volume (mL)'); ylabel('% of total TTB volume removed (pooled)');
        title(sprintf('%s: exploratory minimum-volume feasibility', cfg.tracers{t}));
    end
    p3 = fn('fig_threshold_feasibility','png'); print(f3, p3, '-dpng', '-r120'); close(f3);
    files = [files, {p1, p2, p3}];
catch ME
    otherErrors{end+1} = ['Figures: ' ME.message];
    fprintf('Figure error: %s\n', ME.message);
end

% ---------- verify mandatory outputs (9 CSV + 1 XLSX + 3 PNG) ----------
nMandatory = size(outs,1) + 1 + 3;
for k = 1:numel(files)
    dd = dir(files{k});
    if isempty(dd) || dd.bytes == 0
        otherErrors{end+1} = ['Output missing or empty: ' files{k}]; %#ok<SAGROW>
    end
end
if numel(files) ~= nMandatory
    otherErrors{end+1} = sprintf('Mandatory outputs produced %d of %d', numel(files), nMandatory);
end

% ---------- final run status ----------
if isTest
    completion = 'TEST_ONLY';
elseif dataComplete && isempty(otherErrors)
    completion = 'COMPLETE';
else
    completion = 'INCOMPLETE';
end

%% ------------------------------ summary --------------------------------
fprintf('\n===== SUMMARY =====\n');
for t = 1:2
    for cn = cfg.conns
        s = ST(strcmp(ST.tracer,cfg.tracers{t}) & ST.connectivity==cn, :);
        if isempty(s), continue; end
        fprintf('%-4s %2d-conn: components/patient median %g (range %g-%g) | total components %d | median TTB in components <=9 voxels %.2f%%\n', ...
            cfg.tracers{t}, cn, median(s.n_components), min(s.n_components), max(s.n_components), sum(s.n_components), ...
            100*median(s.fraction_TTB_in_components_le9voxels,'omitnan'));
    end
    for p = {[6 26],[6 18],[18 26]}
        m = strcmp(MT.tracer,cfg.tracers{t}) & MT.conn_from==p{1}(1) & MT.conn_to==p{1}(2);
        if ~any(m), continue; end
        fprintf('%-4s merges %d->%d: median %% TTB volume in merges %.2f (max %.2f); median count reduction %g\n', ...
            cfg.tracers{t}, p{1}(1), p{1}(2), median(MT.pct_TTB_volume_in_merges(m),'omitnan'), ...
            max(MT.pct_TTB_volume_in_merges(m)), median(MT.n_reduction(m)));
    end
end
fprintf('Large components examined anatomically: %d (spanning >= 2 labelled groups: %d; automated observation)\n', ...
    height(AT), nnz(AT.n_labelled_groups_spanned >= 2));
fprintf('Elapsed: %.1f min\n', toc(tStart)/60);

fprintf('\n=================================================================\n');
fprintf('AUDIT 1D RUN STATUS: %s\n', completion);
fprintf('Patients found %d / expected %d | tracer examinations OK %d, PARTIAL %d, FAILED %d / expected %d\n', ...
    nCases, cfg.expectedCases, nOK, nPartial, nFailed, nExpScan);
bad = SS(~strcmp(SS.status,'OK'), :);
for k = 1:height(bad)
    fprintf('   %s %s %s: %s\n', bad.case_id{k}, bad.tracer{k}, bad.status{k}, bad.message{k});
end
for k = 1:numel(otherErrors), fprintf('   OUTPUT ERROR: %s\n', otherErrors{k}); end
if strcmp(completion,'INCOMPLETE')
    fprintf('*** INCOMPLETE RUN: results must NOT be reported as a complete Audit 1D. ***\n');
elseif strcmp(completion,'TEST_ONLY')
    fprintf('*** EXECUTION TEST ONLY (%d case(s)): not an Audit 1D result. Set cfg.testN = [] for the definitive run. ***\n', nCases);
end
fprintf('=================================================================\n');

%% ------------------------------ manifest -------------------------------
Mf.audit = 'Audit 1D - tumour-mask characteristics (technical; per tracer)';
Mf.protocol = cfg.protocol;
Mf.completion_status = completion;
Mf.execution_test = isTest; Mf.n_case_folders_found = nCasesFound;
Mf.n_patients_processed = nCases; Mf.n_patients_expected = cfg.expectedCases;
Mf.n_scans_expected = nExpScan; Mf.n_scans_OK = nOK; Mf.n_scans_PARTIAL = nPartial; Mf.n_scans_FAILED = nFailed;
Mf.scan_errors = table2struct(bad);
Mf.output_errors = otherErrors;
Mf.script = scriptPath; Mf.script_MD5 = scriptMD5; Mf.timestamp = ts;
Mf.elapsed_min = toc(tStart)/60;
Mf.matlab_version = version; Mf.matlab_release = version('-release'); Mf.computer = computer;
v = ver('images'); if ~isempty(v), Mf.image_toolbox_version = v(1).Version; end
Mf.config = cfg; Mf.anatomical_groups = G;
Mf.input_inventory_n_files = height(IT);
Mf.input_inventory = table2struct(IT);
Mf.outputs = files;
Mf.note = ['Read-only; TTB masks unchanged (in-memory copies). Native PET-space volumes. No intertracer registration, ', ...
           'correspondence, overlap or discordance. Thresholds exploratory. Coordinates in NIfTI world space (RAS+, mm). ', ...
           'Anatomical lookup via physical coordinates within the same examination, nearest-neighbour sampling at PET voxel ', ...
           'centres; spanning several labelled structures is an automated observation, not evidence of a single confluent ', ...
           'metastasis. Unlabelled tissue is not automatically nodal.'];
% Manifest safeguard: written first to a temporary file, read back and decoded,
% its completion_status checked; only a verified manifest is renamed to its
% final name. A failed or partial write therefore never leaves a file named
% "manifest" that could be read as a complete run.
fM = fn('manifest','json'); fTmp = [fM '.tmp']; manifestOK = false; manifestErr = '';
try
    try, js = jsonencode(Mf,'PrettyPrint',true); catch, js = jsonencode(Mf); end
    fid = fopen(fTmp,'w'); if fid < 0, error('cannot open %s', fTmp); end
    nW = fwrite(fid, js, 'char'); cs = fclose(fid);
    if nW ~= numel(js) || cs ~= 0, error('incomplete write (%d of %d bytes)', nW, numel(js)); end
    chk = jsondecode(fileread(fTmp));
    if ~isfield(chk,'completion_status') || ~strcmp(chk.completion_status, completion)
        error('manifest read-back does not match the run status');
    end
    [okMv, mvMsg] = movefile(fTmp, fM); if ~okMv, error('rename failed: %s', mvMsg); end
    manifestOK = true;
catch ME
    manifestErr = ME.message;
    fprintf('MANIFEST ERROR: %s\n', manifestErr);
    if exist(fTmp,'file') == 2, try, delete(fTmp); catch, end, end   % remove only the script's own partial temp file
    if strcmp(completion,'COMPLETE'), completion = 'INCOMPLETE'; end % a run without a verified manifest is not complete
    fprintf('*** INCOMPLETE RUN (no verified manifest). ***\n');
end

fid = fopen(fn(['RUN_' completion],'txt'),'w');
fprintf(fid, 'Audit 1D run %s: %s\nPatients %d/%d; scans OK %d, PARTIAL %d, FAILED %d of %d expected.\n', ...
    ts, completion, nCases, cfg.expectedCases, nOK, nPartial, nFailed, nExpScan);
if manifestOK
    fprintf(fid, 'Manifest: written and verified (%s).\n', fM);
else
    fprintf(fid, 'Manifest: NOT WRITTEN (%s). This run must not be reported as a complete Audit 1D.\n', manifestErr);
end
fclose(fid);
fprintf('Outputs written to %s (timestamp %s)\n', cfg.outDir, ts);
diary off;

%% ============================ local functions ===========================
function A = niiAffine(I)
    % 4x4 mapping 0-based voxel indices -> NIfTI world (RAS+, mm); identical to Audits 1A/1B
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

function [ok, dmax] = sameGrid(I1, A1, I2, A2, tol)
    % identical first three dimensions AND identical spatial affine (within tol, mm)
    s1 = double(I1.ImageSize(1:3)); s2 = double(I2.ImageSize(1:3));
    dmax = max(abs(A1(:) - A2(:)));
    ok = isequal(s1, s2) && dmax < tol;
end

function s = axisCode(M)
    % anatomical direction (NIfTI RAS+) of increasing voxel index along i, j, k
    L = {'L','R'; 'P','A'; 'I','S'}; s = '';
    for ax = 1:3
        v = M(:,ax); [~, d] = max(abs(v)); s = [s L{d, (v(d) > 0) + 1}]; %#ok<AGROW>
    end
end

function s = sizeStr(v)
    s = strjoin(arrayfun(@(x) sprintf('%d', x), double(v(:))', 'UniformOutput', false), 'x');
end

function s = affStr(A)
    % full 3x4 spatial affine, row-major
    s = strtrim(sprintf('%.6f ', A(1:3,:)'));
end

function T = mkT(cells, names)
    if isempty(cells), cells = cell(0, numel(names)); end
    T = cell2table(cells, 'VariableNames', names);
end

function h = fileMD5(f)
    h = '';
    if isempty(f), return; end
    fid = fopen(f, 'r'); if fid < 0, return; end
    b = fread(fid, Inf, '*uint8'); fclose(fid);
    md = java.security.MessageDigest.getInstance('MD5');
    if ~isempty(b), md.update(typecast(b, 'int8')); end
    dg = typecast(md.digest(), 'uint8');
    h = lower(reshape(dec2hex(dg, 2)', 1, []));
end

function v = prctile90(x)
    x = sort(x(~isnan(x))); if isempty(x), v = NaN; else, v = x(max(1, round(0.9*numel(x)))); end
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
