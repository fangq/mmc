%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-Domain D-coefficient Jacobian comparison:
% Redbird (FEM diffusion) vs MMC (RF adjoint_dcoeff MC)
%
% Computes J_D(r) = -grad(phi_src) . grad(phi_det) * dV per voxel.
%
% Redbird: J_D_rb(r) from rbjac, complex [1 x nnodes]
% MMC:     J_D_mmc(r) from outputtype='adjoint_dcoeff' + omega>0, method='grid'
%          raw output is complex single [dim.x,dim.y,dim.z,maxgate,Ns*Nd]
%          where dim=[61,61,31] for meshgrid6([0 60],[0 60],[0 30]);
%          trim 1 boundary voxel per axis to get [60,60,30] matching domain.
%
% Agreement is expected in the diffusive regime (away from source/boundary).
%
% Geometry: slab domain [0,60]x[0,60]x[0,30] mm, 1 source, 1 detector.
% See also: demo_validate_mmc_adjoint_jacobian.m
%           demo_validate_mmc_adjoint_mua_d.m
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

cfg.srcpos = [20  29.5  0];   % source at top surface [mm]
cfg.srcdir = [0     0     1]; % pointing into medium (+z)
cfg.detpos = [40    30    0]; % detector at top surface [mm]
cfg.detdir = [0     0     1]; % pointing into medium (+z)

cfg.prop  = [0 0 1 1; mua musp g nmed];
cfg.omega = omega;            % non-zero omega: frequency-domain FEM

cfg = rbmeshprep(cfg);
sd  = rbsdmap(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Redbird: solve FD forward + compute D-coeff Jacobian
%%
%%  rbrunforward solves from BOTH source and detector positions.
%%  phi_rb is [nnodes x 2]: col 1 = phi_src, col 2 = phi_det.
%%  rbjac 3rd output = J_D, complex [1 x nnodes].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Building Redbird FD system (f = %.0f MHz)...\n', freq / 1e6);
tic;
[detphi_rb, phi_rb] = rbrunforward(cfg, 'solverflag', {'qmr', 1e-7, 1000});
toc;

[~, ~, J_D_rb] = rbjac(sd, phi_rb, cfg.deldotdel, cfg.elem, cfg.evol);  % [1 x nnodes] complex

fprintf('Redbird: detphi = %.4e + %.4ei  (amp=%.4e, phase=%.2f deg)\n', ...
        real(detphi_rb(1)), imag(detphi_rb(1)), ...
        abs(detphi_rb(1)), angle(detphi_rb(1)) * 180 / pi);
fprintf('J_D_rb size: %s,  complex: %d\n', mat2str(size(J_D_rb)), ~isreal(J_D_rb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MMC RF adjoint_dcoeff setup  (same physical domain)
%%
%%  outputtype='adjoint_dcoeff' + omega>0:
%%    - computes J_D(r) = -grad(phi_src) . grad(phi_det) * dV per voxel
%%    - output is complex single [Nx,Ny,Nz,maxgate,Ns*Nd] = [60,60,30,1,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('mmclab', 'file')
    error('mmclab not found. Add mmclab to your MATLAB path.');
end

mcfg.nphoton    = 1e8;
mcfg.gpuid = 2;
% Bounding-box tetrahedral mesh for grid-mode output
[mcfg.node, mcfg.elem] = meshgrid6([0 60], [0 60], [0 30]);
mcfg.elem(:, 1:4) = meshreorient(mcfg.node(:, 1:3), mcfg.elem(:, 1:4));
mcfg.elemprop   = ones(size(mcfg.elem, 1), 1);
mcfg.method     = 'grid';

mcfg.prop       = [0 0 1 1; mua musp g nmed];
mcfg.tstart     = 0;
mcfg.tend       = 5e-9;
mcfg.tstep      = 5e-9;
mcfg.isreflect  = 1;

mcfg.srctype    = 'pencil';
mcfg.srcpos     = [20  29.5  0];
mcfg.srcdir     = [0     0     1];

mcfg.detpos     = [40    30    0    1];  % position + radius 1 mm
mcfg.detdir     = [0     0     1    0];  % +z, focal length=0

mcfg.omega      = omega;
mcfg.outputtype = 'adjoint_dcoeff';    % compute J_D = -grad(phi_src).grad(phi_det)*dV

fprintf('Running MMC RF adjoint_dcoeff (f = %.0f MHz, %.0e photons)...\n', freq / 1e6, mcfg.nphoton);
tic;
flux = mmclab(mcfg);
toc;

% flux.data = forward fluence (complex double [dim.x,dim.y,dim.z,maxgate,srcnum+detnum])
%             last dim has srcnum forward-source slices then detnum detector-as-source slices
% flux.jd   = D-coefficient Jacobian (complex single [dim.x,dim.y,dim.z,maxgate,nsrcpairs])
% meshgrid6([0 60],[0 60],[0 30]) → dim=[61,61,31]; trim boundary to [60,60,30]
J_D_mmc = single(squeeze(flux.jd));           % [61, 61, 31] complex single
J_D_mmc = J_D_mmc(1:end - 1, 1:end - 1, 1:end - 1); % trim boundary voxel → [60, 60, 30]
fprintf('J_D_mmc size: %s,  complex: %d\n', mat2str(size(J_D_mmc)), ~isreal(J_D_mmc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract xz cross-section at y = 29.5 mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ym_idx = 30;
x_mmc  = (0.5:1:59.5);
z_mmc  = (0.5:1:29.5);

J_D_mmc_grid = squeeze(J_D_mmc(:, ym_idx, :)).';  % [Nz=30, Nx=60] complex

J_D_rb_nodes = J_D_rb(1, :).';   % [nnodes x 1] complex

[cutpos_rb, cutval_Dre] = qmeshcut(cfg.elem, cfg.node, real(J_D_rb_nodes), 'y=29.5');
[~,         cutval_Dim] = qmeshcut(cfg.elem, cfg.node, imag(J_D_rb_nodes), 'y=29.5');

[xi, zi] = meshgrid(x_mmc, z_mmc);
J_D_rb_grid = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_Dre, xi, zi) + ...
              1i * griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_Dim, xi, zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 1: Jacobian amplitude maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'RF D-coeff Jacobian Amplitude (MMC vs Redbird)', 'Position', [50 100 1400 430]);

subplot(1, 3, 1);
imagesc(x_mmc, z_mmc, log10(abs(J_D_mmc_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC adj |J_D| log_{10} (f=%g MHz)', freq / 1e6));

subplot(1, 3, 2);
imagesc(x_mmc, z_mmc, log10(abs(J_D_rb_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FEM |J_D| log_{10} (f=%g MHz)', freq / 1e6));

subplot(1, 3, 3);
imagesc(x_mmc, z_mmc, log10(abs(J_D_mmc_grid) + eps) - log10(abs(J_D_rb_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title('log_{10}|J_D| difference (MMC - Redbird)');

sgtitle_compat(gcf, sprintf('RF adjoint_dcoeff |J_D| comparison  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 2: Jacobian phase maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'RF D-coeff Jacobian Phase (MMC vs Redbird)', 'Position', [50 100 1400 430]);

subplot(1, 3, 1);
imagesc(x_mmc, z_mmc, angle(J_D_mmc_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC adj angle(J_D) (deg, f=%g MHz)', freq / 1e6));

subplot(1, 3, 2);
imagesc(x_mmc, z_mmc, angle(J_D_rb_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FEM angle(J_D) (deg, f=%g MHz)', freq / 1e6));

subplot(1, 3, 3);
imagesc(x_mmc, z_mmc, angle(J_D_mmc_grid) * 180 / pi - angle(J_D_rb_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title('Phase difference MMC - Redbird (deg)');

sgtitle_compat(gcf, sprintf('RF adjoint_dcoeff phase(J_D) comparison  |  f=%g MHz', freq / 1e6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 3: Contour overlay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'D-coeff Contour Overlay', 'Position', [50 50 1200 500]);

clines_amp = [-12:0.5:-6 -5.5:0.25:-4];

subplot(1, 2, 1);
contour(x_mmc, z_mmc, log10(abs(J_D_rb_grid) + eps), clines_amp, 'r-', 'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, log10(abs(J_D_mmc_grid) + eps), clines_amp, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_dcoeff', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}|J_D| contour overlay (f=%g MHz)', freq / 1e6));
axis equal tight;

subplot(1, 2, 2);
clines_phs = -180:10:180;
contour(x_mmc, z_mmc, angle(J_D_rb_grid) * 180 / pi, clines_phs, 'r-', 'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, angle(J_D_mmc_grid) * 180 / pi, clines_phs, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_dcoeff', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title('angle(J_D) contour overlay (deg)');
axis equal tight;

sgtitle_compat(gcf, 'Contour overlay: Redbird (red solid) vs MMC adjoint\_dcoeff (blue dashed)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 4: Depth profiles at x = 29.5 mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'D-coeff Depth Profiles', 'Position', [100 100 900 500]);

J_D_mmc_line = J_D_mmc_grid(:, 30);
J_D_rb_line  = J_D_rb_grid(:,  30);

subplot(1, 2, 1);
semilogy(z_mmc, abs(J_D_mmc_line) / max(abs(J_D_mmc_line)), 'b-',  'LineWidth', 1.5);
hold on;
semilogy(z_mmc, abs(J_D_rb_line)  / max(abs(J_D_rb_line)),  'r--', 'LineWidth', 1.5);
legend('MMC RF adjoint\_dcoeff', 'Redbird FEM', 'Location', 'northeast');
xlabel('Depth z (mm)');
ylabel('|J_D| (normalised)');
title('J_D amplitude profile at x=29.5 mm');
grid on;

subplot(1, 2, 2);
plot(z_mmc, angle(J_D_mmc_line) * 180 / pi, 'b-',  'LineWidth', 1.5);
hold on;
plot(z_mmc, angle(J_D_rb_line)  * 180 / pi, 'r--', 'LineWidth', 1.5);
legend('MMC RF adjoint\_dcoeff', 'Redbird FEM', 'Location', 'southwest');
xlabel('Depth z (mm)');
ylabel('Phase (deg)');
title('J_D phase profile at x=29.5 mm');
grid on;

sgtitle_compat(gcf, sprintf('Depth profiles at x=29.5 mm  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Print summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n========== Summary ==========\n');
fprintf('Frequency:  %.0f MHz\n', freq / 1e6);
fprintf('mua=%.4f/mm,  musp=%.1f/mm,  n=%.2f\n', mua, musp, nmed);
fprintf('\nRedbird detector value:  %.4e + %.4ei  (amp=%.4e, phase=%.2f deg)\n', ...
        real(detphi_rb(1)), imag(detphi_rb(1)), abs(detphi_rb(1)), angle(detphi_rb(1)) * 180 / pi);
fprintf('\nJ_D_rb  sum(real) = %.4e\n', sum(real(J_D_rb(1, :))));
fprintf('J_D_mmc sum(real) = %.4e\n', sum(real(J_D_mmc(:))));
fprintf('================================\n');
fprintf('Note: MMC (RTE) and Redbird (diffusion) agree well in the diffusive regime.\n');
fprintf('Deviations near the source and boundaries are expected.\n');

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
