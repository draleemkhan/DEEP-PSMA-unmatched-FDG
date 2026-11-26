function PERT_generate_v1_0(mode)
% PERT_generate_v1_0 - PAPER 3: generate the registration perturbations of SAP section 8.1 ONCE (MATLAB R2025a).
% Implements: consolidated Step 2 SAP LOCKED (MD5 0c64711db1207e31a70cf3f558ccbd7d) section 8.1; dated implementation
% clarification 00_protocol\StepSAP_perturbation_implementation_clarification_20261112.md
% (MD5 5f280df40970eec504bb4f08d421d69d), Q1-Q2.
%   rng(20261113, 'twister'); D = randn(50,3) (translation directions), then X = randn(50,3) (rotation axes); rows
%   normalised to unit length (uniform on the sphere). Paired levels: L1 5.16 mm / 2.29 deg; L2 8.50 mm / 4.47 deg.
%   R = axis-angle (Rodrigues). Output: 06_output\PERT_perturbations_<run>\PERT_perturbations_<run>.csv (%.17g) + manifest.
% No patient data are read.
% Usage (current folder = 02_code): PERT_generate_v1_0('selftest'); PERT_generate_v1_0('run')

if nargin < 1, mode = 'selftest'; end
restoredefaultpath;
switch lower(mode)
    case 'selftest'
        if ~selftest(), error('PERT:selftest', 'self-test failed'); end
    case 'run'
        runGen();
    otherwise
        error('PERT:mode', 'mode must be selftest or run');
end
end

function L = levels()
L = struct('name', {'L1', 'L2'}, 'tmag', {5.16, 8.50}, 'angle_deg', {2.29, 4.47});
end

function [D, X] = draws()
rng(20261113, 'twister');
D = randn(50, 3); D = D ./ vecnorm(D, 2, 2);
X = randn(50, 3); X = X ./ vecnorm(X, 2, 2);
end

function R = rodrigues(ax, deg)
a = ax(:) / norm(ax); th = deg * pi / 180;
K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
end

function runGen()
runId = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out = 'D:\DEEP PSMA\PAPER3_Discordance\06_output'; proto = 'D:\DEEP PSMA\PAPER3_Discordance\00_protocol';
nf = fullfile(proto, 'StepSAP_perturbation_implementation_clarification_20261112.md');
if ~strcmp(md5file(nf), '5f280df40970eec504bb4f08d421d69d'), error('PERT:C1', 'STOP C1: clarification note MD5 mismatch'); end
[D, X] = draws(); L = levels();
od = fullfile(out, ['PERT_perturbations_' runId]); mkdir(od);
cf = fullfile(od, ['PERT_perturbations_' runId '.csv']);
fid = fopen(cf, 'w');
fprintf(fid, 'level,k,tmag_mm,angle_deg,dir_x,dir_y,dir_z,axis_x,axis_y,axis_z,t_x,t_y,t_z,R11,R12,R13,R21,R22,R23,R31,R32,R33\n');
for l = 1:numel(L)
    for k = 1:50
        R = rodrigues(X(k,:), L(l).angle_deg); t = L(l).tmag * D(k,:);
        fprintf(fid, '%s,%d,%.17g,%.17g', L(l).name, k, L(l).tmag, L(l).angle_deg);
        fprintf(fid, ',%.17g', [D(k,:), X(k,:), t, reshape(R', 1, 9)]);
        fprintf(fid, '\n');
    end
end
fclose(fid);
man = struct('stage', 'PERT generation of registration perturbations (SAP 8.1)', 'run', runId, ...
    'script', 'PERT_generate_v1_0.m', 'script_md5', md5file([mfilename('fullpath') '.m']), ...
    'clarification_md5', '5f280df40970eec504bb4f08d421d69d', 'seed', 20261113, 'generator', 'twister', ...
    'n_per_level', 50, 'levels', L, 'csv', ['PERT_perturbations_' runId '.csv'], 'csv_md5', md5file(cf), ...
    'R_storage', 'row-major R11..R33; perturbed point = R (q - c) + c + t', 'matlab', version, 'status', 'COMPLETE');
mf = fullfile(od, ['PERT_manifest_' runId '.json']);
fid = fopen(mf, 'w'); fwrite(fid, jsonencode(man, 'PrettyPrint', true)); fclose(fid);
fprintf('done: %s; CSV MD5 %s; manifest MD5 %s\n', od, md5file(cf), md5file(mf));
end

function ok = selftest()
res = cell(0, 2);
[D, X] = draws(); [D2, X2] = draws();
res(end+1,:) = {'draws reproducible with the locked seed', isequal(D, D2) && isequal(X, X2)};
res(end+1,:) = {'directions and axes unit length', max(abs(vecnorm(D, 2, 2) - 1)) < 1e-12 && max(abs(vecnorm(X, 2, 2) - 1)) < 1e-12};
L = levels(); worst = 0; angErr = 0; axErr = 0; tErr = 0;
for l = 1:2
    for k = 1:50
        R = rodrigues(X(k,:), L(l).angle_deg);
        worst = max([worst, norm(R' * R - eye(3)), abs(det(R) - 1)]);
        angErr = max(angErr, abs(acosd(min(1, (trace(R) - 1) / 2)) - L(l).angle_deg));
        axErr = max(axErr, norm(R * X(k,:)' - X(k,:)'));
        tErr = max(tErr, abs(norm(L(l).tmag * D(k,:)) - L(l).tmag));
    end
end
res(end+1,:) = {sprintf('rotation matrices orthonormal, det 1 (max error %.1e)', worst), worst < 1e-12};
res(end+1,:) = {sprintf('rotation angle recovered (max error %.1e deg)', angErr), angErr < 1e-6};
res(end+1,:) = {'rotation axis invariant', axErr < 1e-12};
res(end+1,:) = {'translation magnitude exact', tErr < 1e-12};
Rz = rodrigues([0 0 1], 90);
res(end+1,:) = {'Rodrigues: 90 deg about z maps x to y', norm(Rz * [1;0;0] - [0;1;0]) < 1e-12};
c = [10; -5; 3]; q = c + [4; 0; 0]; t = [1; 2; 3];
p = Rz * (q - c) + c + t; pinv = Rz' * (p - c - t) + c;
res(end+1,:) = {'perturbation about a centre and its inverse', norm(pinv - q) < 1e-12 && norm(Rz * (c - c) + c - c) == 0};
rng(1, 'twister'); U = randn(1e5, 3); U = U ./ vecnorm(U, 2, 2);
res(end+1,:) = {'normalised randn is uniform on the sphere (moments)', max(abs(mean(U))) < 0.01 && max(abs(var(U) - 1/3)) < 0.01};
for k = 1:size(res,1), fprintf('  [%s] %s\n', tern(res{k,2}, 'PASS', 'FAIL'), res{k,1}); end
ok = all([res{:,2}]);
fprintf('selftest: %d/%d passed\n', sum([res{:,2}]), size(res,1));
end

function h = md5file(f)
md = java.security.MessageDigest.getInstance('MD5');
fid = fopen(f, 'r'); if fid < 0, error('PERT:io', 'cannot open %s', f); end
data = fread(fid, Inf, '*uint8'); fclose(fid);
md.update(data);
h = lower(reshape(dec2hex(typecast(md.digest, 'uint8'), 2)', 1, []));
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
