%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dual-Jacobian comparison: Redbird (FEM diffusion) vs MMC (RF adjoint_mua_d)
%
% Demonstrates outputtype='adjoint_mua_d', which computes BOTH J_mua and J_D
% in a single MMC session (one set of forward/adjoint photon propagations).
%
%   J_mua(r) = -phi_src(r) * phi_det(r) * dV      [mua sensitivity]
%   J_D(r)   = -grad(phi_src) . grad(phi_det) * dV [diffusion-coeff sensitivity]
%
% Output: forward fluence in flux.data; Jacobians in separate fields:
%   flux.jmua -> J_mua (raw dim=[61,61,31], trim to [60,60,30])
%   flux.jd   -> J_D
%
% Both are validated against Redbird FEM diffusion. Agreement is expected in
% the diffusive regime (away from source/boundary).
%
% Geometry: slab domain [0,60]x[0,60]x[0,30] mm, 1 source, 1 detector.
% See also: demo_validate_mmc_adjoint_jacobian.m
%           demo_validate_mmc_adjoint_jacobian_dcoeff.m
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
cfg.omega = omega;

cfg = rbmeshprep(cfg);
sd  = rbsdmap(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Redbird: one FD solve -> extract both J_mua and J_D
%%
%%  phi_rb is [nnodes x 2]: col 1 = phi_src, col 2 = phi_det.
%%  J_mua_rb = rbfemmatrix(...)         -> [1 x nnodes] complex
%%  J_D_rb   = rbjac(...) 3rd output   -> [1 x nnodes] complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Building Redbird FD system (f = %.0f MHz)...\n', freq / 1e6);
tic;
[detphi_rb, phi_rb] = rbrunforward(cfg, 'solverflag', {'qmr', 1e-7, 1000});
toc;

J_mua_rb = rbfemmatrix(cfg, sd, phi_rb);                                    % [1 x nnodes] complex
[~, ~, J_D_rb] = rbjac(sd, phi_rb, cfg.deldotdel, cfg.elem, cfg.evol);     % [1 x nnodes] complex

fprintf('Redbird: detphi = %.4e + %.4ei  (amp=%.4e, phase=%.2f deg)\n', ...
        real(detphi_rb(1)), imag(detphi_rb(1)), ...
        abs(detphi_rb(1)), angle(detphi_rb(1)) * 180 / pi);
fprintf('J_mua_rb size: %s,  complex: %d\n', mat2str(size(J_mua_rb)), ~isreal(J_mua_rb));
fprintf('J_D_rb   size: %s,  complex: %d\n', mat2str(size(J_D_rb)),   ~isreal(J_D_rb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MMC: single session with outputtype='adjoint_mua_d'
%%
%%  One photon run computes both J_mua and J_D simultaneously.
%%  One photon run computes both J_mua and J_D simultaneously.
%%  Forward fluence is returned in flux.data (complex double [dim.x,dim.y,dim.z,maxgate]).
%%  Jacobians are in separate fields:
%%    flux.jmua  ->  J_mua   (point product of fluences, complex single [dim.x,dim.y,dim.z,maxgate,Ns*Nd])
%%    flux.jd    ->  J_D     (dot product of fluence gradients, complex single [...])
%%  meshgrid6([0 60],[0 60],[0 30]) → dim=[61,61,31]; trim 1 boundary voxel
%%  per axis to get [60,60,30] matching the physical domain [0,60]x[0,60]x[0,30]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('mmclab', 'file')
    error('mmclab not found. Add mmclab to your MATLAB path.');
end

mcfg.nphoton    = 1e8;

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
mcfg.outputtype = 'adjoint_mua_d';   % single session: output J_mua + J_D together

fprintf('Running MMC RF adjoint_mua_d (f = %.0f MHz, %.0e photons)...\n', freq / 1e6, mcfg.nphoton);
tic;
flux = mmclab(mcfg);
toc;

% flux.data = forward fluence (complex double [dim.x,dim.y,dim.z,maxgate,srcnum+detnum])
%             last dim has srcnum forward-source slices then detnum detector-as-source slices
% flux.jmua = mua Jacobian (complex single [61,61,31,1,1])
% flux.jd   = D Jacobian   (complex single [61,61,31,1,1])
fprintf('flux.jmua size: %s,  complex: %d\n', mat2str(size(flux.jmua)), ~isreal(flux.jmua));

% Trim boundary voxel (mesh nodes 0..60 → dim=61; interior voxels 1..60 = [0,60))
J_mua_mmc = single(squeeze(flux.jmua));                          % [61, 61, 31] complex
J_mua_mmc = J_mua_mmc(1:end - 1, 1:end - 1, 1:end - 1);         % trim to [60, 60, 30]
J_D_mmc   = single(squeeze(flux.jd));                            % [61, 61, 31] complex
J_D_mmc   = J_D_mmc(1:end - 1, 1:end - 1, 1:end - 1);           % trim to [60, 60, 30]

fprintf('J_mua_mmc sum(real) = %.4e\n', sum(real(J_mua_mmc(:))));
fprintf('J_D_mmc   sum(real) = %.4e\n', sum(real(J_D_mmc(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Interpolate Redbird onto MMC regular grid
%%   xz cross-section at y = 29.5 mm (MMC voxel index 30)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ym_idx = 30;
x_mmc  = (0.5:1:59.5);
z_mmc  = (0.5:1:29.5);
[xi, zi] = meshgrid(x_mmc, z_mmc);

% MMC xz slices [Nz x Nx]
J_mua_mmc_grid = squeeze(J_mua_mmc(:, ym_idx, :)).';  % [30 x 60] complex
J_D_mmc_grid   = squeeze(J_D_mmc(:,   ym_idx, :)).';  % [30 x 60] complex

% Redbird: cut at y=29.5 and interpolate Re/Im parts separately
J_mua_rb_nodes = J_mua_rb(1, :).';
J_D_rb_nodes   = J_D_rb(1,   :).';

[cutpos_rb, cutval_mua_re] = qmeshcut(cfg.elem, cfg.node, real(J_mua_rb_nodes), 'y=29.5');
[~,         cutval_mua_im] = qmeshcut(cfg.elem, cfg.node, imag(J_mua_rb_nodes), 'y=29.5');
[~,         cutval_D_re]   = qmeshcut(cfg.elem, cfg.node, real(J_D_rb_nodes),   'y=29.5');
[~,         cutval_D_im]   = qmeshcut(cfg.elem, cfg.node, imag(J_D_rb_nodes),   'y=29.5');

J_mua_rb_grid = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_mua_re, xi, zi) + ...
                1i * griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_mua_im, xi, zi);
J_D_rb_grid   = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_D_re, xi, zi) + ...
                1i * griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_D_im, xi, zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 1: Amplitude maps -- J_mua (top) and J_D (bottom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Dual-Jacobian amplitude: J_mua and J_D', 'Position', [50 80 1400 750]);

titles_mua = {sprintf('MMC |J_{mua}| log_{10}  (f=%g MHz)', freq / 1e6), ...
              sprintf('Redbird |J_{mua}| log_{10}  (f=%g MHz)', freq / 1e6), ...
              'log_{10}|J_{mua}| diff (MMC - Redbird)'};
titles_D   = {sprintf('MMC |J_D| log_{10}  (f=%g MHz)', freq / 1e6), ...
              sprintf('Redbird |J_D| log_{10}  (f=%g MHz)', freq / 1e6), ...
              'log_{10}|J_D| diff (MMC - Redbird)'};

data_mua = {log10(abs(J_mua_mmc_grid) + eps), log10(abs(J_mua_rb_grid) + eps), ...
            log10(abs(J_mua_mmc_grid) + eps) - log10(abs(J_mua_rb_grid) + eps)};
data_D   = {log10(abs(J_D_mmc_grid) + eps), log10(abs(J_D_rb_grid) + eps), ...
            log10(abs(J_D_mmc_grid) + eps) - log10(abs(J_D_rb_grid) + eps)};

for k = 1:3
    subplot(2, 3, k);
    imagesc(x_mmc, z_mmc, data_mua{k});
    axis equal tight;
    colorbar;
    xlabel('x (mm)');
    ylabel('z (mm)');
    title(titles_mua{k});

    subplot(2, 3, k + 3);
    imagesc(x_mmc, z_mmc, data_D{k});
    axis equal tight;
    colorbar;
    xlabel('x (mm)');
    ylabel('z (mm)');
    title(titles_D{k});
end

sgtitle_compat(gcf, sprintf('Dual-Jacobian amplitude  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 2: Phase maps -- J_mua (top) and J_D (bottom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Dual-Jacobian phase: J_mua and J_D', 'Position', [50 80 1400 750]);

titles_mua_ph = {sprintf('MMC angle(J_{mua}) (deg, f=%g MHz)', freq / 1e6), ...
                 sprintf('Redbird angle(J_{mua}) (deg, f=%g MHz)', freq / 1e6), ...
                 'Phase diff J_{mua}: MMC - Redbird (deg)'};
titles_D_ph   = {sprintf('MMC angle(J_D) (deg, f=%g MHz)', freq / 1e6), ...
                 sprintf('Redbird angle(J_D) (deg, f=%g MHz)', freq / 1e6), ...
                 'Phase diff J_D: MMC - Redbird (deg)'};

data_mua_ph = {angle(J_mua_mmc_grid) * 180 / pi, angle(J_mua_rb_grid) * 180 / pi, ...
               angle(J_mua_mmc_grid) * 180 / pi - angle(J_mua_rb_grid) * 180 / pi};
data_D_ph   = {angle(J_D_mmc_grid) * 180 / pi, angle(J_D_rb_grid) * 180 / pi, ...
               angle(J_D_mmc_grid) * 180 / pi - angle(J_D_rb_grid) * 180 / pi};

for k = 1:3
    subplot(2, 3, k);
    imagesc(x_mmc, z_mmc, data_mua_ph{k});
    axis equal tight;
    colorbar;
    xlabel('x (mm)');
    ylabel('z (mm)');
    title(titles_mua_ph{k});

    subplot(2, 3, k + 3);
    imagesc(x_mmc, z_mmc, data_D_ph{k});
    axis equal tight;
    colorbar;
    xlabel('x (mm)');
    ylabel('z (mm)');
    title(titles_D_ph{k});
end

sgtitle_compat(gcf, sprintf('Dual-Jacobian phase  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 3: Contour overlay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Dual-Jacobian contour overlay', 'Position', [50 50 1200 700]);

clines_amp = [-12:0.5:-6 -5.5:0.25:-4];
clines_phs = -180:10:180;

subplot(2, 2, 1);
contour(x_mmc, z_mmc, log10(abs(J_mua_rb_grid) + eps), clines_amp, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, log10(abs(J_mua_mmc_grid) + eps), clines_amp, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_mua\_d', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}|J_{mua}| contour (f=%g MHz)', freq / 1e6));
axis equal tight;

subplot(2, 2, 2);
contour(x_mmc, z_mmc, log10(abs(J_D_rb_grid) + eps), clines_amp, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, log10(abs(J_D_mmc_grid) + eps), clines_amp, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_mua\_d', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}|J_D| contour (f=%g MHz)', freq / 1e6));
axis equal tight;

subplot(2, 2, 3);
contour(x_mmc, z_mmc, angle(J_mua_rb_grid) * 180 / pi, clines_phs, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, angle(J_mua_mmc_grid) * 180 / pi, clines_phs, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_mua\_d', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title('angle(J_{mua}) contour (deg)');
axis equal tight;

subplot(2, 2, 4);
contour(x_mmc, z_mmc, angle(J_D_rb_grid) * 180 / pi, clines_phs, 'r-',  'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, angle(J_D_mmc_grid) * 180 / pi, clines_phs, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint\_mua\_d', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title('angle(J_D) contour (deg)');
axis equal tight;

sgtitle_compat(gcf, 'Contour overlay: Redbird (red solid) vs MMC adjoint\_mua\_d (blue dashed)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 4: Depth profiles at x = 29.5 mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Dual-Jacobian depth profiles', 'Position', [100 100 1100 700]);

J_mua_mmc_line = J_mua_mmc_grid(:, 30);
J_mua_rb_line  = J_mua_rb_grid(:,  30);
J_D_mmc_line   = J_D_mmc_grid(:, 30);
J_D_rb_line    = J_D_rb_grid(:,  30);

subplot(2, 2, 1);
semilogy(z_mmc, abs(J_mua_mmc_line) / max(abs(J_mua_mmc_line)), 'b-',  'LineWidth', 1.5);
hold on;
semilogy(z_mmc, abs(J_mua_rb_line)  / max(abs(J_mua_rb_line)),  'r--', 'LineWidth', 1.5);
legend('MMC adjoint\_mua\_d', 'Redbird FEM', 'Location', 'northeast');
xlabel('Depth z (mm)');
ylabel('|J_{mua}| (normalised)');
title('J_{mua} amplitude profile at x=29.5 mm');
grid on;

subplot(2, 2, 2);
plot(z_mmc, angle(J_mua_mmc_line) * 180 / pi, 'b-',  'LineWidth', 1.5);
hold on;
plot(z_mmc, angle(J_mua_rb_line)  * 180 / pi, 'r--', 'LineWidth', 1.5);
legend('MMC adjoint\_mua\_d', 'Redbird FEM', 'Location', 'southwest');
xlabel('Depth z (mm)');
ylabel('Phase (deg)');
title('J_{mua} phase profile at x=29.5 mm');
grid on;

subplot(2, 2, 3);
semilogy(z_mmc, abs(J_D_mmc_line) / max(abs(J_D_mmc_line)), 'b-',  'LineWidth', 1.5);
hold on;
semilogy(z_mmc, abs(J_D_rb_line)  / max(abs(J_D_rb_line)),  'r--', 'LineWidth', 1.5);
legend('MMC adjoint\_mua\_d', 'Redbird FEM', 'Location', 'northeast');
xlabel('Depth z (mm)');
ylabel('|J_D| (normalised)');
title('J_D amplitude profile at x=29.5 mm');
grid on;

subplot(2, 2, 4);
plot(z_mmc, angle(J_D_mmc_line) * 180 / pi, 'b-',  'LineWidth', 1.5);
hold on;
plot(z_mmc, angle(J_D_rb_line)  * 180 / pi, 'r--', 'LineWidth', 1.5);
legend('MMC adjoint\_mua\_d', 'Redbird FEM', 'Location', 'southwest');
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
fprintf('\nJ_mua_rb  sum(real) = %.4e\n', sum(real(J_mua_rb(1, :))));
fprintf('J_mua_mmc sum(real) = %.4e\n', sum(real(J_mua_mmc(:))));
fprintf('\nJ_D_rb    sum(real) = %.4e\n', sum(real(J_D_rb(1, :))));
fprintf('J_D_mmc   sum(real) = %.4e\n', sum(real(J_D_mmc(:))));
fprintf('================================\n');
fprintf('Note: MMC (RTE) and Redbird (diffusion) agree well in the diffusive regime.\n');
fprintf('Both Jacobians were computed from a SINGLE MMC session using\n');
fprintf('outputtype=''adjoint_mua_d'' (forward in flux.data, Jacobians in flux.jmua and flux.jd).\n');

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
