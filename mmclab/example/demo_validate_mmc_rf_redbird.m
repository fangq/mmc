%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-Domain comparison: Redbird (FEM diffusion) vs MMC (RF forward MC)
%
% Redbird solves the frequency-domain diffusion equation:
%   -div( D * grad(phi) ) + (mua + i*omega/(c/n)) * phi = q
% MMC RF forward uses complex photon weights per Eq.9:
%   w <- w * exp[ -(mua + i*omega*n/c0) * ds ]
%
% Both produce complex-valued fluence. We compare amplitude and phase.
%
% MMC uses method='grid' with a bounding-box tetrahedral mesh.
% Raw output dim=[61,61,31] for meshgrid6([0 60],[0 60],[0 30]);
% trim 1 boundary voxel per axis to get [60,60,30] matching the domain.
%
% See also: demo_validate_mmc_adjoint_jacobian.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg mcfg;

addpath('../../matlab/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Shared optical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq  = 200e6;              % modulation frequency [Hz]
omega = 2 * pi * freq;      % angular frequency [rad/s]

mua  = 0.005;               % absorption coeff [1/mm]
musp = 1.0;                 % reduced scattering [1/mm]
g    = 0;                   % anisotropy (isotropic, for diffusion comparison)
nmed = 1.37;                % refractive index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Redbird FEM setup  (domain [0,60]x[0,60]x[0,30] mm)
%%   Uses a fine uniform tetrahedral mesh (1 mm resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cfg.node, cfg.elem] = meshgrid6(0:1:60, 0:1:60, 0:1:30);
cfg.elem(:, 1:4) = meshreorient(cfg.node(:, 1:3), cfg.elem(:, 1:4));
cfg.face = volface(cfg.elem);

cfg.seg = ones(size(cfg.elem, 1), 1);

cfg.srcpos = [29.5 29.5 0];   % source at top surface [mm]
cfg.srcdir = [0    0    1];   % pointing into medium (+z)
cfg.detpos = [40   30   0];   % detector position for comparison
cfg.detdir = [0    0    1];

cfg.prop  = [0 0 1 1; mua musp g nmed];
cfg.omega = omega;             % non-zero omega: frequency-domain FEM

cfg = rbmeshprep(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Redbird FD solve
%%   phi_rb is [nnodes x 1]: forward fluence from source
%%   rbfemgetdet interpolates detector value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Building Redbird FD system (f = %.0f MHz)...\n', freq / 1e6);
[Amat, ~] = rbfemlhs(cfg);
[rhs, ~, ~] = rbfemrhs(cfg);

tic;
fprintf('Solving Redbird FD...\n');
phi_rb = rbfemsolve(Amat, rhs, 'qmr', 1e-7, 1000);
toc;

fprintf('Redbird output is complex: %d\n', ~isreal(phi_rb));

% Redbird CW for reference
cfg0       = cfg;
cfg0.omega = 0;
cfg0       = rbmeshprep(cfg0);
[Amat0]    = rbfemlhs(cfg0);
[rhs0]     = rbfemrhs(cfg0);
phi_rb0    = rbfemsolve(Amat0, rhs0, 'qmr', 1e-8, 1000);
phi_rb0(phi_rb0 < 0) = 0;

detval_rb = rbfemgetdet(phi_rb, cfg, rhs);
fprintf('Redbird detector: amplitude = %.6e, phase = %.2f deg\n', ...
        abs(detval_rb), angle(detval_rb) * 180 / pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MMC RF Forward setup  (same physical domain)
%%
%%  omega > 0 activates complex weight accumulation per segment.
%%  Output flux.data is complex double [dim.x, dim.y, dim.z, maxgate].
%%  For meshgrid6([0 60],[0 60],[0 30]): dim=[61,61,31].
%%  Multiply by tstep to convert flux rate to fluence.
%%  Trim boundary voxel to [60,60,30] for domain-interior comparison.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('mmclab', 'file')
    error('mmclab not found. Add mmclab to your MATLAB path.');
end

mcfg.nphoton = 1e7;

% Bounding-box tetrahedral mesh for grid-mode output
[mcfg.node, mcfg.elem] = meshgrid6([0 60], [0 60], [0 30]);
mcfg.elem(:, 1:4) = meshreorient(mcfg.node(:, 1:3), mcfg.elem(:, 1:4));
mcfg.elemprop = ones(size(mcfg.elem, 1), 1);
mcfg.method   = 'grid';       % BLBadouelGrid: uniform voxel output at 1 mm step

mcfg.prop     = [0 0 1 1; mua musp g nmed];
mcfg.tstart   = 0;
mcfg.tend     = 5e-9;
mcfg.tstep    = 5e-9;
mcfg.isreflect = 1;

mcfg.srctype  = 'pencil';
mcfg.srcpos   = [29.5 29.5 0];  % match Redbird source position [mm]
mcfg.srcdir   = [0    0    1];  % +z into medium

mcfg.omega    = omega;           % RF modulation angular frequency (omega>0 → complex weights)

fprintf('Running MMC RF forward (f = %.0f MHz, %.0e photons)...\n', freq / 1e6, mcfg.nphoton);
tic;
flux_rf = mmclab(mcfg);
toc;

% RF forward output: complex double [dim.x, dim.y, dim.z, maxgate] = [61, 61, 31, 1]
% Multiply by tstep to convert fluence rate to fluence (same as CW workflow)
phi_mmc = flux_rf.data * mcfg.tstep;  % complex double [61, 61, 31, 1]
phi_mmc = squeeze(phi_mmc);           % [61, 61, 31]
% Trim boundary voxel: mesh nodes 0..60 → dim=61; interior voxels cover [0,60)
phi_mmc = phi_mmc(1:end - 1, 1:end - 1, 1:end - 1);  % [60, 60, 30]
fprintf('MMC RF output is complex: %d, size: %s\n', ~isreal(phi_mmc), mat2str(size(phi_mmc)));

% MMC CW for reference
mcfg_cw          = mcfg;
mcfg_cw          = rmfield(mcfg_cw, 'omega');
mcfg_cw.outputtype = 'flux';

fprintf('Running MMC CW (%.0e photons)...\n', mcfg_cw.nphoton);
tic;
flux_cw = mmclab(mcfg_cw);
toc;

phi_cw = squeeze(flux_cw.data) * mcfg_cw.tstep;  % [61, 61, 31]
phi_cw = phi_cw(1:end - 1, 1:end - 1, 1:end - 1);       % [60, 60, 30]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract cross-sections at y = 29.5 mm
%%   MMC trimmed array [60,60,30]: voxel iy=29 (MATLAB index 30)
%%   covers y in [29,30), center at y=29.5.
%%   Redbird: qmeshcut slice at y=29.5.
%%   Both grids: x_mmc = 0.5:1:59.5, z_mmc = 0.5:1:29.5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ym_idx = 30;               % MATLAB 1-based, center at y=29.5 mm
x_mmc  = (0.5:1:59.5);    % [1 x 60]
z_mmc  = (0.5:1:29.5);    % [1 x 30]

% MMC: extract xz slice [Nx=60, Nz=30], transpose to [Nz x Nx] for imagesc
mmc_amp_yz   = squeeze(abs(phi_mmc(:, ym_idx, :))).';    % [30 x 60]
mmc_phase_yz = squeeze(angle(phi_mmc(:, ym_idx, :))).' * 180 / pi;
cw_yz        = squeeze(real(phi_cw(:,  ym_idx, :))).';   % [30 x 60]

% Redbird: cut at y=29.5, grid onto regular grid
[cutpos_rb, cutval_amp]   = qmeshcut(cfg.elem, cfg.node, abs(phi_rb(:, 1)),                  'y=29.5');
[~,         cutval_phase] = qmeshcut(cfg.elem, cfg.node, angle(phi_rb(:, 1)) * 180 / pi,     'y=29.5');
[~,         cutval_cw]    = qmeshcut(cfg.elem, cfg.node, real(phi_rb0(:, 1)),                 'y=29.5');

[xi, zi] = meshgrid(x_mmc, z_mmc);  % [30 x 60]
rb_amp_grid   = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_amp,   xi, zi);
rb_phase_grid = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_phase, xi, zi);
rb_cw_grid    = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_cw,    xi, zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 1: Side-by-side 2D maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amp_clim = [-6 0];
phs_clim = [-60 5];

figure('Name', 'FD Comparison: Maps', 'Position', [50 100 1400 800]);

subplot(2, 3, 1);
imagesc(x_mmc, z_mmc, log10(cw_yz + eps));
axis equal tight;
colorbar;
caxis(amp_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title('MMC CW log_{10}(fluence)');

subplot(2, 3, 2);
imagesc(x_mmc, z_mmc, log10(mmc_amp_yz + eps));
axis equal tight;
colorbar;
caxis(amp_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC RF Amp log_{10} (f=%g MHz)', freq / 1e6));

subplot(2, 3, 3);
imagesc(x_mmc, z_mmc, log10(rb_amp_grid + eps));
axis equal tight;
colorbar;
caxis(amp_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FD Amp log_{10} (f=%g MHz)', freq / 1e6));

subplot(2, 3, 4);
imagesc(x_mmc, z_mmc, zeros(length(z_mmc), length(x_mmc)));
axis equal tight;
colorbar;
caxis(phs_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title('MMC CW phase (0 deg)');

subplot(2, 3, 5);
imagesc(x_mmc, z_mmc, mmc_phase_yz);
axis equal tight;
colorbar;
caxis(phs_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC RF Phase (deg, f=%g MHz)', freq / 1e6));

subplot(2, 3, 6);
imagesc(x_mmc, z_mmc, rb_phase_grid);
axis equal tight;
colorbar;
caxis(phs_clim);
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FD Phase (deg, f=%g MHz)', freq / 1e6));

sgtitle_compat(gcf, sprintf('Frequency-Domain: Redbird (FEM) vs MMC RF Forward  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 2: Contour overlay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Contour Comparison', 'Position', [50 50 1200 500]);

clines_amp   = -7:0.5:0;
clines_phase = -60:5:0;

subplot(1, 2, 1);
contour(x_mmc, z_mmc, log10(rb_amp_grid + eps),  clines_amp, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, log10(mmc_amp_yz  + eps),  clines_amp, 'b--', 'LineWidth', 2);
legend('Redbird FD', 'MMC RF Forward', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}(Amplitude)  f=%g MHz', freq / 1e6));
axis equal tight;

subplot(1, 2, 2);
contour(x_mmc, z_mmc, rb_phase_grid,  clines_phase, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, mmc_phase_yz,   clines_phase, 'b--', 'LineWidth', 2);
legend('Redbird FD', 'MMC RF Forward', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Phase (deg)  f=%g MHz', freq / 1e6));
axis equal tight;

sgtitle_compat(gcf, 'Contour overlay: Redbird (red solid) vs MMC RF (blue dashed)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 3: Depth profiles at (x=29.5, y=29.5)
%%   Trimmed array index (30, 30, :) → x=29.5, y=29.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Depth Profiles', 'Position', [100 100 1000 700]);

z_line = z_mmc;   % voxel centers 0.5, 1.5, ..., 29.5

mmc_amp_line   = squeeze(abs(phi_mmc(30, ym_idx, :)));
mmc_phase_line = squeeze(angle(phi_mmc(30, ym_idx, :))) * 180 / pi;
cw_mmc_line    = squeeze(real(phi_cw(30, ym_idx, :)));

rb_amp_line   = rb_amp_grid(:, 30)';
rb_phase_line = rb_phase_grid(:, 30)';
rb_cw_line    = rb_cw_grid(:, 30)';

subplot(2, 2, 1);
semilogy(z_line, cw_mmc_line, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(z_line, rb_cw_line,  'r--', 'LineWidth', 1.5);
legend('MMC CW', 'Redbird CW', 'Location', 'northeast');
ylabel('Fluence');
xlabel('Depth z (mm)');
title('CW Fluence');
grid on;

subplot(2, 2, 2);
semilogy(z_line, mmc_amp_line, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(z_line, rb_amp_line,  'r--', 'LineWidth', 1.5);
legend('MMC RF', 'Redbird FD', 'Location', 'northeast');
ylabel('Amplitude');
xlabel('Depth z (mm)');
title(sprintf('FD Amplitude (f=%g MHz)', freq / 1e6));
grid on;

subplot(2, 2, 3);
plot(z_line, mmc_phase_line, 'b-', 'LineWidth', 1.5);
hold on;
plot(z_line, rb_phase_line,  'r--', 'LineWidth', 1.5);
legend('MMC RF', 'Redbird FD', 'Location', 'southwest');
ylabel('Phase (deg)');
xlabel('Depth z (mm)');
title(sprintf('FD Phase (f=%g MHz)', freq / 1e6));
grid on;

subplot(2, 2, 4);
mmc_ratio = mmc_amp_line(:) ./ max(cw_mmc_line(:), eps);
rb_ratio  = rb_amp_line(:)  ./ max(rb_cw_line(:),  eps);
plot(z_line, mmc_ratio, 'b-', 'LineWidth', 1.5);
hold on;
plot(z_line, rb_ratio,  'r--', 'LineWidth', 1.5);
legend('MMC', 'Redbird', 'Location', 'southwest');
ylabel('AC/DC ratio');
xlabel('Depth z (mm)');
title('Amplitude demodulation ratio');
ylim([0 1.1]);
grid on;

sgtitle_compat(gcf, sprintf('Depth profiles at (x=29.5, y=29.5)  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm, n=%.2f', ...
                            freq / 1e6, mua, musp, nmed));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Print summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n========== Summary ==========\n');
fprintf('Frequency:  %.0f MHz\n', freq / 1e6);
fprintf('mua=%.4f/mm,  musp=%.1f/mm,  n=%.2f\n', mua, musp, nmed);
fprintf('\nRedbird detector (at [40,30,0]): amp=%.4e, phase=%.2f deg\n', ...
        abs(detval_rb), angle(detval_rb) * 180 / pi);
% MMC: voxel (40,30,1) in trimmed array → center at (39.5, 29.5, 0.5)
mmc_det = phi_mmc(40, ym_idx, 1);
fprintf('MMC voxel (40,30,1):             amp=%.4e, phase=%.2f deg\n', ...
        abs(mmc_det), angle(mmc_det) * 180 / pi);
fprintf('Amp ratio (MMC/Redbird):   %.4f\n', abs(mmc_det) / abs(detval_rb));
fprintf('Phase diff (MMC-Redbird):  %.2f deg\n', ...
        angle(mmc_det) * 180 / pi - angle(detval_rb) * 180 / pi);
fprintf('\nNote: MMC (transport) and Redbird (diffusion) agree well in\n');
fprintf('the diffusive regime (far from source/boundary).\n');
fprintf('================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sgtitle_compat(figh, str)
    if exist('sgtitle', 'builtin') || exist('sgtitle', 'file')
        sgtitle(str);
    else
        oldax = gca;
        axes('Parent', figh, 'Position', [0 0.96 1 0.04], 'Visible', 'off');
        text(0.5, 0.5, str, 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'tex');
        axes(oldax);
    end
end
