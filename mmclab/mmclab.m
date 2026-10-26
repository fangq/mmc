function varargout = mmclab(varargin)
%
% #############################################################################%
%                         Mesh-based Monte Carlo (MMC)                        %
%          Copyright (c) 2010-2026 Qianqian Fang <q.fang at neu.edu>          %
%              https://mcx.space/#mmc  &  https://neurojson.io                %
%                                                                             %
% Computational Optics & Translational Imaging (COTI) Lab [https://fanglab.org]%
%   Department of Bioengineering, Northeastern University, Boston, MA, USA    %
% #############################################################################%
%    The MCX Project is funded by the NIH/NIGMS under grant R01-GM114365      %
% #############################################################################%
%  Open-source codes and reusable scientific data are essential for research, %
% MCX proudly developed human-readable JSON-based data formats for easy reuse.%
%                                                                             %
% Please visit our free scientific data sharing portal at https://neurojson.io%
% and consider sharing your public datasets in standardized JSON/JData format %
% #############################################################################%
%
% Format:
%    fluence = mmclab(cfg);
%    [fluence, detphoton, ncfg, seeds, traj] = mmclab(cfg);
%    [fluence, detphoton, ncfg, seeds, traj] = mmclab(cfg, option);
%    ncfg = mmclab(cfg, 'prep');
%    info = mmclab('gpuinfo');
%
%    Request only the outputs you need: detphoton, seeds and traj are only
%    produced when they are requested, and [~, detphoton] = mmclab(cfg) skips
%    the fluence in Octave (set cfg.issave2pt=0 to do the same in MATLAB).
%
% Input:
%    cfg: a struct, or a struct array; each element is one simulation.
%         cfg='gpuinfo' returns a struct array describing the OpenCL devices.
%    option: (optional) a string
%         'prep':    only preprocess cfg (fill in the fields marked "-" below)
%                    and return it; reuse it to skip the preprocessing later
%         'preview': plot the mesh, sources and detectors with mcxpreview(cfg)
%         'sse', 'opencl', 'cuda' or 'optix': run on that backend, the same as
%                    setting cfg.compute
%
%    A variable USE_MCXCL in the base workspace set to an integer N>1 replaces
%    cfg.gpuid by N in all following mmclab calls, to pick one device on a
%    multi-GPU computer without editing scripts (as it does for mcxlabcl).
%
%    Fields marked with * are required, fields marked with - are computed when
%    absent (giving them saves time), and values in [] are the defaults.
%    Lengths are in the unit of the mesh nodes, which cfg.unitinmm gives in mm.
%
% == Mesh and optical properties ==
%     *cfg.node:       mesh nodes, an Nnx3 array of [x y z]
%     *cfg.elem:       tetrahedral elements, an Nex4 array of 1-based node indices;
%                      a 5th column, if present, is used as cfg.elemprop
%     *cfg.elemprop:   the medium label (from 1) of each element, Nex1
%      cfg.face:       OptiX only: a surface mesh given instead of cfg.elem, an
%                      Nfx5 array of [n1 n2 n3 inner_label outer_label], whose
%                      triangle normals (right-hand rule) point from the inner to
%                      the outer medium; it is ignored when cfg.elem is given
%     *cfg.prop:       an Nx4 array, one row per medium, [mua(1/mm), mus(1/mm), g, n];
%                      the first row is label 0 (outside of the mesh), the second
%                      row is label 1, and so on
%     -cfg.facenb:     element face neighbors, from faceneighbors(cfg.elem)
%     -cfg.evol:       element volumes, from elemvolume(cfg.node, cfg.elem)
%     -cfg.e0:         the element enclosing the source, from tsearchn(); set it to
%                      '>', '<' or '-' to search along srcdir, against srcdir, or
%                      both, from srcpos to the mesh surface and move the source
%                      there
%      cfg.unitinmm:   the length unit of the mesh in mm [1.0]; e.g. 10 for a mesh
%                      in cm
%      cfg.nout:       refractive index outside of the mesh [1.0]
%      cfg.nodemua:    per-node mua, Nnx1, replacing the mua of cfg.prop (used by
%                      the redbird reconstruction); each element uses the mean of
%                      its nodes
%      cfg.nodemusp:   per-node mus', Nnx1, used together with cfg.nodemua; g and n
%                      still come from cfg.prop
%      implicit MMC (iMMC), CPU only (cfg.gpuid=-1), see <demo_immc_basic.m> and
%      <demo_immc_vessel.m>:
%      cfg.noderoi:    a vector of sphere radii, one per node
%      cfg.edgeroi:    an Nex6 array of cylinder radii, one per element edge, in
%                      the nchoosek(1:4,2) order
%      cfg.faceroi:    an Nex4 array of slab thicknesses, one per element face
%
% == Time gates ==
%     *cfg.tend:       end of the simulated time window (s)
%     *cfg.tstart:     start of the time window (s)
%     *cfg.tstep:      width of each time gate (s)
%
% == Source ==
%     *cfg.srcpos:     [x y z] position of the source
%     *cfg.srcdir:     [vx vy vz] unit launch direction; an optional 4th value is
%                      the focal length: 0 collimated, >0 converging to a focus at
%                      that distance along srcdir, <0 diverging from a virtual
%                      point behind the source (no effect on pencil, isotropic
%                      and cone sources)
%                      several sources: give srcpos/srcdir (and srcparam1/2) one
%                      row per source
%      cfg.srctype:    source type ['pencil']; cfg.srcparam1/srcparam2 hold the
%                      parameters. Wide-field sources must start outside of the
%                      mesh: mmclab then adds the source domain to the mesh.
%                      'pencil'      a pencil beam
%                      'isotropic'   isotropic point source
%                      'arcsine'     like isotropic, with a uniform zenith angle
%                      'cone'        uniform cone; srcparam1(1): half-angle (rad)
%                      'zgaussian'   Gaussian zenith angle; srcparam1(1): variance
%                      'gaussian'    Gaussian beam; srcparam1(1): waist radius; with
%                                    a focal length, the waist is at srcpos+focus*srcdir
%                                    and srcparam1(2), the wavelength in mesh units,
%                                    makes the radius grow away from it
%                      'disk'        uniform disk normal to srcdir; srcparam1(1):
%                                    radius, srcparam1(2): inner radius of an annulus
%                      'planar'      uniform quadrilateral with corners srcpos,
%                                    srcpos+srcparam1(1:3), srcpos+srcparam2(1:3)
%                      'pattern'     the 'planar' quadrilateral with intensity
%                                    cfg.srcpattern (0-1) of srcparam1(4) x srcparam2(4)
%                                    pixels; cfg.srcnum patterns at once are given as
%                                    a [srcnum, srcparam1(4), srcparam2(4)] array.
%                                    Example: <demo_photon_sharing.m>
%                      'pattern3d'   a 3D intensity volume cfg.srcpattern of size
%                                    srcparam1(1:3), starting at srcpos
%                      'fourier'     spatial-frequency source over the 'planar'
%                                    quadrilateral, S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)]:
%                                    integer parts of srcparam1(4)/srcparam2(4) are
%                                    fx/fy, fraction of srcparam1(4) is phi0/(2*pi),
%                                    1 minus the fraction of srcparam2(4) is M
%                      'fourierx'    general Fourier source; srcparam1: [v1x v1y v1z
%                                    |v2|], srcparam2: [kx ky phi0 M], with
%                                    v2=srcdir x v1; the phase is phi0*2*pi
%                      'fourierx2d'  general 2D Fourier basis; srcparam1 as for
%                                    'fourierx', srcparam2: [kx ky phix phiy]
%                      'line'        line source from srcpos to srcpos+srcparam1(1:3),
%                                    emitting normal to the line
%                      'slit'        collimated slit beam along srcdir from the same
%                                    segment
%                      CPU (sse) and OptiX only:
%                      'pencilarray' an array of pencil beams over the 'planar'
%                                    quadrilateral; srcparam1(4)/srcparam2(4) are the
%                                    beam counts
%                      'hyperboloid' one-sheet hyperboloid Gaussian beam; srcparam1:
%                                    [waist radius, focal distance, Rayleigh range]
%                      'ring'        ring or ring sector; srcparam1: [outer radius,
%                                    inner radius, start angle, end angle (rad)]
%      cfg.srcparam1, cfg.srcparam2: 1x4 vectors, see cfg.srctype
%      cfg.srcpattern: the pattern of 'pattern'/'pattern3d' sources
%      cfg.srcnum:     number of patterns simulated together by a 'pattern'
%                      source [1]
%      cfg.srcid:      with several sources: [0] all together in one output, -1 each
%                      in its own output, -2 also run each detector as a source
%                      (forward fluence only), N>0 only the N-th source
%      cfg.voidtime:   for a source outside of the mesh, [1]: start the clock at
%                      launch, 0: when the photon enters the mesh
%      cfg.isspecular: 1: apply the specular reflection at entry for a source
%                      outside of the mesh; [0]: do not
%      cfg.omega:      modulation angular frequency 2*pi*f (rad/s) [0]; >0 makes a
%                      frequency-domain (RF) run, see cfg.outputtype
%
% == Detectors ==
%      cfg.detpos:     an Nx4 array, one row per disk detector: [x y z radius]
%      cfg.detdir:     an Nx4 array, one row per detector: [vx vy vz focal length],
%                      the inward direction of each detector, required by the
%                      adjoint outputs (the detectors are launched as sources)
%      cfg.sdmap:      an Nch x 2 (or wider) array selecting the channels of the
%                      adjoint Jacobians: each row is [srcid detid], both 1-based
%                      (srcid counts the forward sources, detid the rows of
%                      cfg.detpos); extra columns are ignored. The Jacobians then
%                      have Nch channels, in the row order of cfg.sdmap, instead
%                      of all Ns*Nd pairs. Default: all pairs, with the detector
%                      index varying the fastest, i.e. [srcid detid] =
%                      [1 1; 1 2; ... 1 Nd; 2 1; ...]
%      cfg.detparam1, cfg.detparam2: [ux uy uz Nx] and [vx vy vz Ny], the edges and
%                      pixel counts of a wide-field camera detector spanning
%                      detpos, detpos+u, detpos+u+v and detpos+v; used with
%                      cfg.issaveexit=2. Example: <demo_wide_det.m>
%      cfg.maxdetphoton: the most detected photons kept [1000000]
%      cfg.issaveexit: 1: also save the exit position and direction of each
%                      detected photon [0]; 2: instead, add up the exiting weight
%                      on the wide-field detector of cfg.detparam1/detparam2, an
%                      Nx x Ny image (detector radii are then ignored)
%      cfg.ismomentum: 1: also save the momentum transfer [0]
%
% == Monte Carlo settings ==
%      cfg.seed:       RNG seed, an integer [1648335518]; a uint8 array (the seeds
%                      output of an earlier run) replays those detected photons,
%                      see "Replay" below
%      cfg.isreflect:  boundary condition at the mesh surface: [1] Fresnel
%                      reflection, 0 index-matched (photons leave), 2 total
%                      absorption, 3 mirror (total reflection)
%      cfg.isnormalized: [1]: normalize the output to a unit source; 0: raw output
%      cfg.minenergy:  weight below which a photon enters Russian roulette [1e-6]
%      cfg.roulettesize: the Russian roulette survival factor [10]
%      cfg.mcmethod:   [0]: MCX-style (Beer-Lambert absorption along the path);
%                      1: MCML-style (albedo weight)
%      cfg.method:     ray-tracing algorithm:
%                      ['grid']: dual-grid MMC (DMMC), a voxelated output
%                      'elem':  branch-less Badouel (SSE), a mesh output
%                      'plucker', 'havel', 'badouel': other mesh tracers, CPU only
%                      (GPU runs use 'elem' in their place)
%      cfg.steps:      [dx dy dz], the voxel size of the 'grid' output, all equal
%                      [1 1 1]
%      cfg.basisorder: mesh output only: [1] values at the nodes (linear basis),
%                      0 values in the elements (constant basis)
%
% == Output ==
%      cfg.outputtype: what the 1st output holds ['flux']
%                      'flux'      fluence rate (1/mm^2/s)
%                      'fluence'   fluence (1/mm^2), integrated over each time gate
%                      'energy'    energy deposited
%                      'jacobian' or 'wl'  replay: mua Jacobian (weighted path length)
%                      'wp'        replay: weighted scattering counts (mus Jacobian)
%                      'rf'        with cfg.omega>0 and no replay: the complex
%                                  fluence rate
%                      adjoint Jacobians (need cfg.detpos and cfg.detdir; every
%                      detector is also run as a source, and the Jacobian of
%                      each source-detector pair is returned in its own field of
%                      the 1st output, [Nx,Ny,Nz,Ngate,Nch] in grid mode or
%                      [Nn,Nch] in mesh mode, where Nch is Ns*Nd or the rows of
%                      cfg.sdmap; complex if cfg.omega>0):
%                      'adjoint'          mua Jacobian, fluence.jmua
%                      'adjoint_dcoeff'   diffusion coefficient D=1/(3mus')
%                                         Jacobian, fluence.jd
%                      'adjoint_mus'      mus Jacobian, fluence.jmus
%                      'adjoint_musp'     mus' Jacobian, fluence.jmusp
%                      'adjoint_mua_d'    both fluence.jmua and fluence.jd
%                      'adjoint_mua_musp' both fluence.jmua and fluence.jmusp
%                      Examples: <demo_mmclab_mesh_adjoint.m>,
%                      <demo_mmclab_mesh_rf_adjoint.m>
%      cfg.adjointmode: mesh-output adjoint mua Jacobian: [0] the full FEM form,
%                      1 the nodal approximation -Vn*phi_src*phi_det
%      cfg.issave2pt:  [1]: compute the 1st output's data; 0: skip it, which is
%                      faster when only detphoton is needed
%      cfg.issaveref:  1: save the diffuse reflectance on the mesh surface in
%                      fluence.dref [0], [#surface triangles, #gates]; the
%                      triangles are faceneighbors(cfg.elem,'rowmajor') (iso2mesh).
%                      Example: <demo_mmclab_basic.m>
%      cfg.session:    a name for the output files
%
% == Replay ==
%      Replaying the detected photons of a baseline run, with cfg.seed set to its
%      seeds output, gives Jacobians with cfg.outputtype. Example:
%      <demo_example_replay.m>
%      cfg.replaydet:  [0]: replay all detectors into one output; N>0: only
%                      detector N
%      cfg.detphotons: the detphoton.data output of the baseline run (single)
%      cfg.replayweight, cfg.replaytime: the exit weight and time of flight of
%                      each replayed photon (single), as an alternative to
%                      cfg.detphotons
%
% == Compute backend ==
%      cfg.compute:    'sse' (CPU, OpenMP), ['opencl'], 'cuda' or 'optix'; 'cuda'
%                      and 'optix' need an mmc mex file built with them
%                      (make cudamex/cudaoct or optixmex/optixoct), and give an
%                      error otherwise. Example: <demo_optix_mmc.m>
%      cfg.gpuid:      which device to use, from 1 (mmclab('gpuinfo') lists the
%                      OpenCL devices) [1]; -1 runs on the CPU (sse); a string of
%                      0s and 1s selects several, e.g. '1101'
%      cfg.workload:   the relative share of photons for each selected device,
%                      e.g. [50 20 30]
%      cfg.nthread:    total GPU thread count [8192]
%      cfg.nblocksize: GPU threads per block [64]
%      cfg.optlevel:   OpenCL kernel optimization level, 0 (none) to 3 [3]
%      cfg.compileropt: extra OpenCL JIT compiler options
%      cfg.kernelfile: an OpenCL kernel file to run instead of the built-in one
%      cfg.isatomic:   [1]: accumulate the output with atomic operations; 0: not
%
% == Debugging ==
%      cfg.debuglevel: a string of debug flags (case insensitive):
%                        'P' progress bar           'T' timing
%                        'M' record trajectories (the traj output)
%                        'S' photon moves           'C' ray-triangle tests
%                        'B' barycentric coordinates 'W' weight changes
%                        'D' distances              'I'/'O' entering/leaving a triangle
%                        'X' hitting an edge        'A' accumulation into the mesh
%                        'R' reflections            'E' exiting photons
%      cfg.debugphoton: print the debug information of this photon only [-1: all]
%      cfg.maxjumpdebug: the most trajectory positions recorded [10000000]
%
%      accepted but without effect: cfg.shapes, cfg.isref3, cfg.outputformat,
%      cfg.isreoriented (set by mmclab itself)
%
% Output:
%      fluence: a struct array, one element per cfg element, with the subfields
%            data:  the quantity selected by cfg.outputtype, normalized as set by
%                   cfg.isnormalized: with method 'grid', a volume of voxels of
%                   size cfg.steps covering the mesh, [Nx,Ny,Nz,#gates]; with
%                   the mesh tracers, [#nodes,#gates] (basisorder 1) or
%                   [#elements,#gates] (basisorder 0); with several patterns or
%                   source slots, an extra dimension; complex for RF runs
%            dref:  the surface diffuse reflectance, when cfg.issaveref=1
%            jmua, jd, jmus, jmusp: the adjoint Jacobians, see cfg.outputtype
%      detphoton: a struct array of the detected photons, one row per photon:
%            detid: the detector index (from 1)
%            nscat: scattering counts in each medium
%            ppath: partial path lengths in each medium, in mesh units
%            mom:   momentum transfer in each medium (cfg.ismomentum=1)
%            p, v:  exit position and direction (cfg.issaveexit=1)
%            w0:    initial weight
%            prop, unitinmm: copies of cfg.prop and cfg.unitinmm
%            data:  all of the above stored, transposed, in one array
%            with cfg.issaveexit=2, the wide-field detector image instead
%      ncfg:  the preprocessed cfg, with the fields marked "-" filled in; pass it
%            to later runs to skip the preprocessing
%      seeds: the RNG seed of each detected photon, one uint8 column per photon,
%            for replay
%      traj:  photon trajectories, sorted by photon, with subfields id (photon
%            index), pos (x/y/z), eid (enclosing element) and data, all of them
%            as rows [id x y z weight eid] of one array
%
% Example:
%      cfg.nphoton = 1e5;
%      [cfg.node, face, cfg.elem] = meshabox([0 0 0], [60 60 30], 6);
%      cfg.elemprop = ones(size(cfg.elem, 1), 1);
%      cfg.srcpos = [30 30 0];
%      cfg.srcdir = [0 0 1];
%      cfg.prop = [0 0 1 1; 0.005 1 0 1.37];
%      cfg.tstart = 0;
%      cfg.tend = 5e-9;
%      cfg.tstep = 5e-10;
%      cfg.debuglevel = 'TP';
%      % preprocess once, to reuse below
%      ncfg = mmclab(cfg, 'prep');
%
%      cfgs(1) = ncfg;   % all elements of a struct array have the same fields
%      cfgs(2) = ncfg;
%      cfgs(1).isreflect = 0;
%      cfgs(2).isreflect = 1;
%      cfgs(2).detpos = [30 20 0 1; 30 40 0 1; 20 30 1 1; 40 30 0 1];
%      % run both, returning the fluence and the detected photons
%      [fluxs, detps] = mmclab(cfgs);
%
% This function is part of Mesh-based Monte Carlo (MMC) URL: https://mcx.space/#mmc
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

try
    defaultocl = evalin('base', 'USE_MCXCL');
catch
    defaultocl = 1;
end

type = '';
if (nargin == 2 && ischar(varargin{2}))
    type = varargin{2};
    if (~any(strcmp(type, {'prep', 'preview', 'sse', 'omp', 'opencl', 'cuda', 'optix'})))
        error('option ''%s'' is not recognized', type);
    end
end

if (isstruct(varargin{1}))
    for i = 1:length(varargin{1})
        % USE_MCXCL = N with |N| > 1 in the base workspace selects the
        % default device for every cfg in this call (parallel to
        % mcxcl/mcxlabcl/mcxlabcl.m). USE_MCXCL = 0 / 1 leaves cfg.gpuid
        % alone (gpuid stays at whatever the cfg sets, or the mex default).
        if (defaultocl < 0)
            varargin{1}(i).gpuid = -defaultocl;
        elseif (defaultocl > 1)
            varargin{1}(i).gpuid = defaultocl;
        end

        % option='sse'/'opencl'/'cuda'/'optix' picks the backend, as cfg.compute
        % does ('omp' is the old name of the CPU backend)
        if (any(strcmp(type, {'sse', 'omp', 'opencl', 'cuda', 'optix'})))
            varargin{1}(i).compute = strrep(type, 'omp', 'sse');
        end

        castlist = {'srcpattern', 'srcpos', 'detpos', 'prop', 'workload', 'srcdir', 'detdir', 'sdmap'};
        for j = 1:length(castlist)
            if (isfield(varargin{1}(i), castlist{j}))
                varargin{1}(i).(castlist{j}) = double(varargin{1}(i).(castlist{j}));
            end
        end
    end
end

if (nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'gpuinfo'))
    varargout{1} = mmc('gpuinfo');
    return
end

if (nargin == 0)
    return
end

cfg = varargin{1};

if (~isstruct(cfg))
    error('cfg must be a struct or struct array');
end

len = length(cfg);
for i = 1:len
    % a surface-only mesh for OptiX: cfg.node plus a 5-column cfg.face, [n1 n2 n3 inner outer],
    % without cfg.elem; the element preprocessing below does not apply
    if (isfield(cfg(i), 'face') && size(cfg(i).face, 2) == 5 && isfield(cfg(i), 'node') && ...
        (~isfield(cfg(i), 'elem') || isempty(cfg(i).elem)))
        if (~(isfield(cfg(i), 'compute') && strcmp(cfg(i).compute, 'optix')))
            error('a surface-only mesh (cfg.face with 5 columns, no cfg.elem) requires cfg.compute=''optix''');
        end
        if (~isfield(cfg(i), 'nphoton'))
            error('cfg.nphoton field is missing');
        end
        if (~isfield(cfg(i), 'prop') || size(cfg(i).prop, 1) < max(max(cfg(i).face(:, 4:5))) + 1)
            error('cfg.prop field is missing or insufficient');
        end
        if (~isfield(cfg(i), 'method') || isempty(cfg(i).method))
            cfg(i).method = 'grid';
        end
        continue
    end
    if (~isfield(cfg(i), 'node') || ~isfield(cfg(i), 'elem'))
        error('cfg.node or cfg.elem is missing');
    end
    if (size(cfg(i).elem, 2) > 4)
        cfg(i).elemprop = cfg(i).elem(:, 5);
    end
    if (~isfield(cfg(i), 'isreoriented') || isempty(cfg(i).isreoriented) || cfg(i).isreoriented == 0)
        if (isfield(cfg(i), 'edgeroi') || isfield(cfg(i), 'faceroi'))
            [cfg(i).elem, evol, idx] = meshreorient(cfg(i).node, cfg(i).elem(:, 1:4));
        else
            [cfg(i).elem, evol] = meshreorient(cfg(i).node, cfg(i).elem(:, 1:4));
        end
        if (isfield(cfg(i), 'edgeroi'))
            cfg(i).edgeroi(idx, :) = cfg(i).edgeroi(idx, [1 3 2 5 4 6]);
        end
        if (isfield(cfg(i), 'faceroi'))
            cfg(i).faceroi(idx, :) = cfg(i).faceroi(idx, [2 1 3 4]);
        end
        cfg(i).isreoriented = 1;
    end
    if (~isfield(cfg(i), 'facenb') || isempty(cfg(i).facenb))
        cfg(i).facenb = faceneighbors(cfg(i).elem);
    end
    if (~isfield(cfg(i), 'evol') || isempty(cfg(i).evol))
        cfg(i).evol = elemvolume(cfg(i).node, cfg(i).elem);
    end
    if (~(isfield(cfg(i), 'compute') && strcmp(cfg(i).compute, 'optix')))
        if (find(cfg(i).evol == 0))
            fprintf(1, ['degenerated elements are detected: [' sprintf('%d ', find(cfg(i).evol == 0)) ']\n']);
            error(['input mesh can not contain degenerated elements, ' ...
                   'please double check your input mesh; if you use a ' ...
                   'widefield source, please rerun mmcsrcdomain and setting ' ...
                   '''Expansion'' option to a larger value (default is 1)']);
        end
    end
    if (~isfield(cfg(i), 'srcpos'))
        error('cfg.srcpos field is missing');
    end
    if (~isfield(cfg(i), 'srcdir'))
        error('cfg.srcdir field is missing');
    end
    if (~isfield(cfg(i), 'e0') || isempty(cfg(i).e0))
        % For multi-source (Mx{3,4} srcpos / srcdir), only the first row
        % populates cfg.srcpos in the mex container; the remaining rows
        % go into cfg.srcdata[] and are looked up per-photon at launch
        % time inside the GPU kernel. Pass only the first 3 coordinates
        % to tsearchn so the e0 hint is for the main source slot.
        cfg(i).e0 = tsearchn(cfg(i).node, cfg(i).elem, cfg(i).srcpos(1, 1:3));
    end
    if ((isnan(cfg(i).e0) && (isfield(cfg(i), 'srctype') && strcmp(cfg(i).srctype, 'pencil'))) || ischar(cfg(i).e0))
        disp('searching initial element ...');
        [cfg(i).srcpos(1, 1:3), cfg(i).e0] = mmcraytrace(cfg(i).node, cfg(i).elem, cfg(i).srcpos(1, 1:3), cfg(i).srcdir(1, 1:3), cfg(i).e0);
    end
    if ((isfield(cfg(i), 'srctype') && strcmp(cfg(i).srctype, 'pattern')) && (ndims(cfg(i).srcpattern) == 2))
        cfg(i).srcpattern = reshape(cfg(i).srcpattern, ...
                                    size(cfg(i).srcpattern, 1), size(cfg(i).srcpattern, 2), 1);
    end
    % Wide-field sources pick a new launch position on the source plane for
    % every photon, so the launch element must be looked up per-photon. The GPU
    % kernels (mmc_core.cl/mmc_core.cu, "for(is=0;is<srcelemlen;is++)") scan
    % only the srcelem[] list that mesh_srcdetelem() builds from elements
    % labeled elemprop<0, and they have no fallback: when srcelemlen==0 the
    % loop body never runs and every photon silently keeps r->eid=cfg.e0,
    % depositing the whole beam into that one tetrahedron instead of the
    % requested pattern. (The CPU tracer starts the same loop at is=-1 and at
    % least raises "initial element does not enclose the source!".)
    %
    % Keying the retessellation below off isnan(cfg.e0) alone is not enough: a
    % cfg.srcpos placed exactly on (or inside) the mesh surface makes tsearchn
    % succeed, so the whole mmcsrcdomain/mmcaddsrc block used to be skipped and
    % the simulation ran on with srcelemlen==0. Trigger it off the source type
    % instead.
    %
    % fourierx/fourierx2d are wide-field too: mmcsrcdomain builds their v2 as
    % srcparam1(4)*(srcdir x v1)/|v1|, matching the launch code in
    % mmc_raytrace.c and mmc_core.cl.
    %
    % optix is excluded because it does its own source handling, never gets the
    % added source domain (see the mmcaddsrc guard below), and reads cfg.e0 to
    % pick launchParams.mediumid0 in mmc_optix_utils.cpp -- entering this block
    % would reset a perfectly good e0 to -1 for it.
    iswidesrc = isfield(cfg(i), 'srctype') && ...
        any(strcmp(cfg(i).srctype, {'planar', 'pattern', 'fourier', 'fourierx', 'fourierx2d'})) && ...
        ~(isfield(cfg(i), 'compute') && strcmp(cfg(i).compute, 'optix'));
    if (isnan(cfg(i).e0) || iswidesrc)  % widefield source
        if (~isfield(cfg(i), 'srcparam1') || ~isfield(cfg(i), 'srcparam2'))
            error('for wide-field sources, you must provide srcparam1 and srcparam2');
        end
        if (~isfield(cfg(i), 'srctype'))
            cfg(i).srctype = 'pencil';
        end
        srcdef = struct('srctype', cfg.srctype, 'srcpos', cfg.srcpos, 'srcdir', cfg.srcdir, ...
                        'srcparam1', cfg.srcparam1, 'srcparam2', cfg.srcparam2);
        sdom = mmcsrcdomain(srcdef, [min(cfg.node); max(cfg.node)]);
        isinside = ismember(round(sdom * 1e10) * 1e-10, round(cfg(i).node * 1e10) * 1e-10, 'rows');
        if (all(~isinside) && ~(isfield(cfg(i), 'compute') && strcmp(cfg(i).compute, 'optix')))
            if (size(cfg(i).elem, 2) == 4)
                cfg(i).elem(:, 5) = cfg(i).elemprop;
            end
            [cfg(i).node, cfg(i).elem] = mmcaddsrc(cfg(i).node, cfg(i).elem, sdom);
            cfg(i).elemprop = cfg(i).elem(:, 5);
            [cfg(i).elem, evol, idx] = meshreorient(cfg(i).node, cfg(i).elem(:, 1:4));
            if (isfield(cfg(i), 'edgeroi'))
                cfg(i).edgeroi(idx, :) = cfg(i).edgeroi(idx, [1 3 2 5 4 6]);
            end
            if (isfield(cfg(i), 'faceroi'))
                cfg(i).faceroi(idx, :) = cfg(i).faceroi(idx, [2 1 3 4]);
            end
            cfg(i).facenb = faceneighbors(cfg(i).elem);
            cfg(i).evol = elemvolume(cfg(i).node, cfg(i).elem);
            cfg(i).isreoriented = 1;
        end
        if (strcmp(cfg(i).srctype, 'pencil') || strcmp(cfg(i).srctype, 'isotropic') || ...
            strcmp(cfg(i).srctype, 'cone')   || strcmp(cfg(i).srctype, 'zgaussian'))
            cfg(i).e0 = tsearchn(cfg(i).node, cfg(i).elem, cfg(i).srcpos);
            if (isnan(cfg(i).e0))
                cfg(i).e0 = -1;
            end
        elseif (length(find(cfg(i).elemprop < 0)) == 1)
            cfg(i).e0 = find(cfg(i).elemprop < 0);
        else
            cfg(i).e0 = -1;
        end
        % refuse to hand the GPU a wide-field source with an empty srcelem[]
        % list -- see the comment above (iswidesrc is already false for optix)
        if (iswidesrc && ~any(cfg(i).elemprop < 0))
            error(['no source-domain element was tessellated into the mesh for ' ...
                   'the ''' cfg(i).srctype ''' wide-field source, so the ' ...
                   'per-photon launch-element search has nothing to search and ' ...
                   'the simulation would deposit the entire beam into element ' ...
                   'cfg.e0. this normally means cfg.srcpos sits on or inside ' ...
                   'the mesh surface: move the source outside the mesh, ' ...
                   'against cfg.srcdir, so that mmcaddsrc can attach the ' ...
                   'source domain (see mmclab/example/demo_compare_mmc_mcx.m)']);
        end
    end
    if (~isfield(cfg(i), 'elemprop'))
        error('cfg.elemprop field is missing');
    end
    if (~isfield(cfg(i), 'nphoton'))
        error('cfg.nphoton field is missing');
    end
    if (~isfield(cfg(i), 'prop') || size(cfg(i).prop, 1) < max(cfg(i).elemprop) + 1 || min(cfg(i).elemprop <= 0))
        error('cfg.prop field is missing or insufficient');
    end
end

if (strcmp(type, 'prep'))
    varargout{1} = cfg;
    return
elseif (strcmp(type, 'preview'))
    [varargout{1:nargout}] = mcxpreview(cfg);
    return
end

% must do a fflush, otherwise octave buffers the output until complete
if (exist('OCTAVE_VERSION'))
    fflush(stdout);
end

if (exist('maxNumCompThreads', 'file'))
    warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
    maxNumCompThreads('automatic');
    maxNumCompThreads(maxNumCompThreads * 2);
end

mmcout = nargout;
if (nargout >= 3)
    mmcout = nargout - 1;
    varargout{nargout} = cfg;
end

for i = 1:length(cfg)
    if (isfield(cfg(i), 'gpuid') && ~ischar(cfg(i).gpuid) && cfg(i).gpuid < -1)
        cfg(i).gpuid = 1;
    end
end

% detect [~,...]=mmclab() in Octave. unfortunately matlab does not have isargout()
if (nargout >= 1 && exist('isargout', 'builtin') && isargout(1) == 0)
    for i = 1:length(cfg)
        cfg(i).issave2pt = 0;
    end
end
[varargout{1:mmcout}] = mmc(cfg);

if (mmcout >= 2)
    for i = 1:length(varargout{2})
        if (~isfield(cfg(i), 'issaveexit') || cfg(i).issaveexit ~= 2)
            medianum = size(cfg(i).prop, 1) - 1;
            detp = varargout{2}(i).data;
            if (isempty(detp))
                continue
            end
            newdetp.detid = int32(detp(1, :))';
            newdetp.nscat = int32(detp(2:medianum + 1, :))';    % 1st medianum block is num of scattering
            newdetp.ppath = detp(medianum + 2:2 * medianum + 1, :)'; % 2nd medianum block is partial path
            if (isfield(cfg(i), 'ismomentum') && cfg(i).ismomentum)
                newdetp.mom = detp(2 * medianum + 2:3 * medianum + 1, :)'; % 3rd medianum block is the momentum transfer
            end
            if (isfield(cfg(i), 'issaveexit') && cfg(i).issaveexit)
                newdetp.p = detp(end - 6:end - 4, :)';             % columns 7-5 from the right store the exit positions*/
                newdetp.v = detp(end - 3:end - 1, :)';      % columns 4-2 from the right store the exit dirs*/
            end
            newdetp.w0 = detp(end, :)';  % last column is the initial packet weight
            if (isfield(cfg(i), 'srcpattern') && ndims(cfg(i).srcpattern) == 3)
                newdetp.w0 = typecast(newdetp.w0, 'uint32');
            end
            newdetp.prop = cfg(i).prop;
            if (isfield(cfg(i), 'unitinmm'))
                newdetp.unitinmm = cfg(i).unitinmm;
            end
            newdetp.data = detp;      % enable this line for compatibility
            newdetpstruct(i) = newdetp;
        else
            newdetpstruct(i) = varargout{2}(i);
        end
    end
    if (exist('newdetpstruct', 'var'))
        varargout{2} = newdetpstruct;
    end
end

if (mmcout > 3)  % the trajectory output only exists when it is requested, even with debuglevel 'M'
    for i = 1:length(varargout{4})
        data = varargout{4}.data;
        if (isempty(data))
            continue
        end
        traj.pos = data(2:4, :).';
        traj.id = typecast(data(1, :), 'uint32').';
        [traj.id, idx] = sort(traj.id);
        traj.pos = traj.pos(idx, :);
        traj.eid = typecast(data(6, :), 'uint32').';
        traj.data = [single(traj.id)'; data(2:end, idx)];
        newtraj(i) = traj;
    end
    if (exist('newtraj', 'var'))
        varargout{4} = newtraj;
    end
end

if (nargout >= 4)
    [varargout{3:end}] = deal(varargout{[end 3:end - 1]});
end
