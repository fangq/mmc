%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-Domain Jacobian comparison: Redbird (FEM diffusion) vs MMC (RF adjoint MC)
%
% Both compute the mua Jacobian J(r) = phi_src(r) * phi_det(r) per voxel,
% where phi_src is the forward fluence from the source and phi_det is the
% (adjoint) forward fluence launched from the detector position.
%
% Redbird: J_rb(r) from rbfemmatrix, complex [1 x nnodes]
% MMC:     J_mmc(r) from outputtype='adjoint' + omega>0, method='grid'
%          raw output is complex single [dim.x,dim.y,dim.z,maxgate,Ns*Nd]
%          where dim = [61,61,31] for meshgrid6([0 60],[0 60],[0 30]);
%          trim 1 boundary voxel per axis to get [60,60,30] matching domain.
%
% Agreement is expected in the diffusive regime (away from source/boundary).
%
% Geometry: slab domain [0,60]x[0,60]x[0,30] mm, 1 source, 1 detector.
% See also: demo_validate_mmc_adjoint_jacobian_dcoeff.m
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
%%   Redbird: solve FD forward + compute mua Jacobian
%%
%%  rbrunforward solves from BOTH source and detector positions.
%%  phi_rb is [nnodes x 2]: col 1 = phi_src, col 2 = phi_det.
%%  J_rb = rbfemmatrix(cfg, sd, phi_rb) returns [1 x nnodes], complex.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Building Redbird FD system (f = %.0f MHz)...\n', freq / 1e6);
tic;
[detphi_rb, phi_rb] = rbrunforward(cfg, 'solverflag', {'qmr', 1e-7, 1000});
toc;

J_rb = rbfemmatrix(cfg, sd, phi_rb);   % [1 x nnodes] complex

fprintf('Redbird: detphi = %.4e + %.4ei  (amp=%.4e, phase=%.2f deg)\n', ...
        real(detphi_rb(1)), imag(detphi_rb(1)), ...
        abs(detphi_rb(1)), angle(detphi_rb(1)) * 180 / pi);
fprintf('J_rb size: %s,  complex: %d\n', mat2str(size(J_rb)), ~isreal(J_rb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MMC RF adjoint setup  (same physical domain [0,60]x[0,60]x[0,30] mm)
%%
%%  Uses method='grid' for voxel-based Jacobian output.
%%  The tetrahedral mesh is just the bounding box; the grid output
%%  is at 1 mm resolution (default cfg.steps = [1 1 1] mm).
%%
%%  outputtype='adjoint' + omega>0:
%%    - sources (forward + detector-as-reversed-source) are simulated together
%%    - adjoint kernel computes J(r) = phi_src(r) * phi_det(r) per voxel
%%    - raw output is complex single [dim.x,dim.y,dim.z,1,1]
%%    - meshgrid6([0 60],[0 60],[0 30]) gives nodes 0..60 so dim=[61,61,31];
%%      trim boundary voxel (1:end-1 per axis) to get [60,60,30]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('mmclab', 'file')
    error('mmclab not found. Add mmclab to your MATLAB path.');
end

mcfg.nphoton    = 1e7;

% Bounding-box tetrahedral mesh: meshgrid6 with coarse spacing creates just
% the 8 corner nodes of the box; grid output fills the interior at 1 mm steps.
[mcfg.node, mcfg.elem] = meshgrid6([0 60], [0 60], [0 30]);
mcfg.elem(:, 1:4) = meshreorient(mcfg.node(:, 1:3), mcfg.elem(:, 1:4));
mcfg.elemprop   = ones(size(mcfg.elem, 1), 1);
mcfg.method     = 'grid';    % BLBadouelGrid: uniform voxel output

mcfg.prop       = [0 0 1 1; mua musp g nmed];
mcfg.tstart     = 0;
mcfg.tend       = 5e-9;
mcfg.tstep      = 5e-9;
mcfg.isreflect  = 1;

mcfg.srctype    = 'pencil';
mcfg.srcpos     = [20  29.5  0];    % matches Redbird srcpos [mm]
mcfg.srcdir     = [0     0     1];  % +z direction

mcfg.detpos     = [40    30    0    1]; % detector position + radius 1 mm
mcfg.detdir     = [0     0     1    0]; % +z (reversed photons enter medium), focal length=0

mcfg.omega      = omega;              % RF modulation frequency [rad/s]
mcfg.outputtype = 'adjoint';          % compute J_mua = phi_src * phi_det per voxel

fprintf('Running MMC RF adjoint (f = %.0f MHz, %.0e photons)...\n', freq / 1e6, mcfg.nphoton);
tic;
flux = mmclab(mcfg);
toc;

% flux.data  = forward fluence (complex double [dim.x,dim.y,dim.z,maxgate,srcnum+detnum])
%              last dim has srcnum forward-source slices then detnum detector-as-source slices
% flux.jmua  = mua Jacobian (complex single [dim.x,dim.y,dim.z,maxgate,nsrcpairs])
%   real part = Re(J) = Re(phi_src)*Re(phi_det) - Im(phi_src)*Im(phi_det)
%   imag part = Im(J) = Re(phi_src)*Im(phi_det) + Im(phi_src)*Re(phi_det)
% meshgrid6([0 60],[0 60],[0 30]) → dim=[61,61,31]; trim boundary to [60,60,30]
J_mmc = single(squeeze(flux.jmua));        % [61, 61, 31] complex single
J_mmc = J_mmc(1:end - 1, 1:end - 1, 1:end - 1); % trim boundary voxel → [60, 60, 30]
fprintf('J_mmc size: %s,  complex: %d\n', mat2str(size(J_mmc)), ~isreal(J_mmc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Extract xz cross-section at y = 29.5 mm
%%   MMC grid voxel index 30 (1-based) has center at y = 29.5 mm
%%   (voxel k center = k - 0.5 mm since nmin = [0,0,0])
%%   Redbird: qmeshcut slices at y = 29.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ym_idx = 30;               % MMC 1-based voxel index with center at y=29.5 mm

% Voxel center coordinates
x_mmc = (0.5:1:59.5);     % [1 x 60]
z_mmc = (0.5:1:29.5);     % [1 x 30]

% MMC: extract xz slice, transpose to [Nz x Nx] for imagesc
J_mmc_xz   = squeeze(J_mmc(:, ym_idx, :));   % [60, 30] complex
J_mmc_grid = J_mmc_xz.';                     % [30, 60] complex (non-conj transpose)

% Redbird: cut mesh at y=29.5, interpolate Re/Im separately onto regular grid
J_rb_nodes = J_rb(1, :).';   % [nnodes x 1] complex

[cutpos_rb, cutval_Jre] = qmeshcut(cfg.elem, cfg.node, real(J_rb_nodes), 'y=29.5');
[~,         cutval_Jim] = qmeshcut(cfg.elem, cfg.node, imag(J_rb_nodes), 'y=29.5');

[xi, zi] = meshgrid(x_mmc, z_mmc);
J_rb_grid = griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_Jre, xi, zi) + ...
            1i * griddata(cutpos_rb(:, 1), cutpos_rb(:, 3), cutval_Jim, xi, zi);
% J_rb_grid: [Nz=30, Nx=60] complex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 1: Jacobian amplitude maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'RF Jacobian Amplitude (MMC vs Redbird)', 'Position', [50 100 1400 430]);

subplot(1, 3, 1);
imagesc(x_mmc, z_mmc, log10(abs(J_mmc_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC adj |J_{mua}| log_{10} (f=%g MHz)', freq / 1e6));

subplot(1, 3, 2);
imagesc(x_mmc, z_mmc, log10(abs(J_rb_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FEM |J_{mua}| log_{10} (f=%g MHz)', freq / 1e6));

subplot(1, 3, 3);
imagesc(x_mmc, z_mmc, log10(abs(J_mmc_grid) + eps) - log10(abs(J_rb_grid) + eps));
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title('log_{10}|J_{mua}| difference (MMC - Redbird)');

sgtitle_compat(gcf, sprintf('RF adjoint |J_{mua}| comparison  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 2: Jacobian phase maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'RF Jacobian Phase (MMC vs Redbird)', 'Position', [50 100 1400 430]);

subplot(1, 3, 1);
imagesc(x_mmc, z_mmc, angle(J_mmc_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('MMC adj angle(J_{mua}) (deg, f=%g MHz)', freq / 1e6));

subplot(1, 3, 2);
imagesc(x_mmc, z_mmc, angle(J_rb_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Redbird FEM angle(J_{mua}) (deg, f=%g MHz)', freq / 1e6));

subplot(1, 3, 3);
imagesc(x_mmc, z_mmc, angle(J_mmc_grid) * 180 / pi - angle(J_rb_grid) * 180 / pi);
axis equal tight;
colorbar;
xlabel('x (mm)');
ylabel('z (mm)');
title('Phase difference MMC - Redbird (deg)');

sgtitle_compat(gcf, sprintf('RF adjoint phase(J_{mua}) comparison  |  f=%g MHz', freq / 1e6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 3: Contour overlay -- banana shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Contour Overlay', 'Position', [50 50 1200 500]);

clines_amp = [-12:0.5:-6 -5.5:0.25:-4];

subplot(1, 2, 1);
contour(x_mmc, z_mmc, log10(abs(J_rb_grid) + eps), clines_amp, 'r-', 'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, log10(abs(J_mmc_grid) + eps), clines_amp, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('log_{10}|J_{mua}| contour overlay (f=%g MHz)', freq / 1e6));
axis equal tight;

subplot(1, 2, 2);
clines_phs = 0:10:270;
contour(x_mmc, z_mmc, angle(J_rb_grid) * 180 / pi, clines_phs, 'r-', 'LineWidth', 2);
hold on;
contour(x_mmc, z_mmc, angle(J_mmc_grid) * 180 / pi, clines_phs, 'b--', 'LineWidth', 2);
legend('Redbird FEM', 'MMC adjoint', 'Location', 'northeast');
xlabel('x (mm)');
ylabel('z (mm)');
title('angle(J_{mua}) contour overlay (deg)');
axis equal tight;

sgtitle_compat(gcf, 'Contour overlay: Redbird (red solid) vs MMC adjoint (blue dashed)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Figure 4: Depth profiles at x = 29.5 mm
%%   Grid column index 30 = x_mmc(30) = 29.5 mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Depth Profiles', 'Position', [100 100 900 500]);

J_mmc_line = J_mmc_grid(:, 30);   % [Nz x 1] complex
J_rb_line  = J_rb_grid(:,  30);   % [Nz x 1] complex

subplot(1, 2, 1);
semilogy(z_mmc, abs(J_mmc_line) / max(abs(J_mmc_line)), 'b-',  'LineWidth', 1.5);
hold on;
semilogy(z_mmc, abs(J_rb_line)  / max(abs(J_rb_line)),  'r--', 'LineWidth', 1.5);
legend('MMC RF adjoint', 'Redbird FEM', 'Location', 'northeast');
xlabel('Depth z (mm)');
ylabel('|J_{mua}| (normalised)');
title('Amplitude depth profile at x=29.5 mm');
grid on;

subplot(1, 2, 2);
plot(z_mmc, angle(J_mmc_line) * 180 / pi, 'b-',  'LineWidth', 1.5);
hold on;
plot(z_mmc, angle(J_rb_line)  * 180 / pi, 'r--', 'LineWidth', 1.5);
legend('MMC RF adjoint', 'Redbird FEM', 'Location', 'southwest');
xlabel('Depth z (mm)');
ylabel('Phase (deg)');
title('Phase depth profile at x=29.5 mm');
grid on;

sgtitle_compat(gcf, sprintf('Depth profiles at x=29.5 mm  |  f=%g MHz, \\mu_a=%.3f/mm, \\mu_s''=%.1f/mm', ...
                            freq / 1e6, mua, musp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Print summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n========== Summary ==========\n');
fprintf('Frequency:  %.0f MHz\n', freq / 1e6);
fprintf('mua=%.4f/mm,  musp=%.1f/mm,  n=%.2f\n', mua, musp, nmed);
fprintf('\nRedbird detector value (src->det):  %.4e + %.4ei\n', ...
        real(detphi_rb(1)), imag(detphi_rb(1)));
fprintf('  |amp| = %.4e,  phase = %.2f deg\n', ...
        abs(detphi_rb(1)), angle(detphi_rb(1)) * 180 / pi);
fprintf('\nJ_rb  sum(real) = %.4e\n', sum(real(J_rb(1, :))));
fprintf('J_mmc sum(real) = %.4e\n', sum(real(J_mmc(:))));
fprintf('================================\n');
fprintf('Note: MMC (transport, RTE) and Redbird (diffusion equation)\n');
fprintf('agree well in the diffusive regime (far from source/boundary).\n');
fprintf('Deviations near the source are expected.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sgtitle_compat(figh, str)
    %% Octave/MATLAB-compatible sgtitle substitute
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
