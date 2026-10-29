function test_mmclab(varargin)
%
% test_mmclab
% test_mmclab('compute', 'opencl', 'gpuid', 1)
%
% Unit tests of the mmclab interface, for MATLAB and GNU Octave: the mesh
% preparation, the errors, the outputs and their sizes, and short absorbed-
% fraction end-to-end checks. The simulations run on the CPU (cfg.compute=
% 'sse') with a few thousand photons, so no GPU is needed; mmclab needs
% iso2mesh on the path. An error is raised if any test fails
%
% input:
%     optional name/value pairs added to the cfg of every test
%
% This file is part of Mesh-based Monte Carlo (MMC) URL: https://mcx.space/#mmc
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
if (~exist('meshabox', 'file') || ~exist('faceneighbors', 'file'))
    error('test_mmclab needs iso2mesh (https://iso2mesh.sf.net) on the path');
end

extra = struct(varargin{:});
tests = {@test_info, @test_errors, @test_prep, @test_outputs, ...
         @test_detected_photons, @test_struct_array, @test_rf, @test_replay, ...
         @test_option_sets_compute, @test_widefield_source, @test_flog_silent, ...
         @test_end_to_end};
failed = {};
for i = 1:numel(tests)
    name = func2str(tests{i});
    try
        tests{i}(extra);
        fprintf(1, '  %-26s ok\n', name);
    catch err
        fprintf(1, '  %-26s FAILED: %s\n', name, err.message);
        failed{end + 1} = [name ': ' err.message];
    end
end

if (~isempty(failed))
    error('%d of %d mmclab tests failed:\n  %s', numel(failed), numel(tests), strjoin(failed, '\n  '));
end
fprintf(1, 'all %d mmclab tests passed\n', numel(tests));

%% -------------------------------------------------------------------------
function cfg = mkcfg(extra, varargin)
% a 20 mm cube of 6 tetrahedra with a pencil beam; name/value pairs override
node = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1] * 20;
elem = [1 2 4 8; 1 3 4 8; 1 2 6 8; 1 5 6 8; 1 3 7 8; 1 5 7 8];
cfg = struct('nphoton', 2000, 'node', node, 'elem', elem, 'elemprop', ones(6, 1), ...
             'srcpos', [10 10 0], 'srcdir', [0 0 1], 'tstart', 0, 'tend', 5e-9, ...
             'tstep', 5e-9, 'prop', [0 0 1 1; 0.05 2 0.8 1], 'compute', 'sse', 'seed', 12345);
if (ispc)
    cfg.flog = 'NUL';  % silence the log
else
    cfg.flog = '/dev/null';
end
opt = [fieldnames(extra), struct2cell(extra); reshape(varargin, 2, [])'];
for i = 1:size(opt, 1)
    cfg.(opt{i, 1}) = opt{i, 2};
end

function cfg = boxcfg(extra, varargin)
% the same cube as a finer tetrahedral mesh, for the wide-field sources
[node, face, elem] = meshabox([0 0 0], [20 20 20], 10);
cfg = mkcfg(extra, 'node', node, 'elem', elem(:, 1:4), 'elemprop', ones(size(elem, 1), 1), varargin{:});

function e = energy(flux)
e = sum(flux.data(:));

function assert_error(fun, msg)
try
    fun();
catch
    return
end
error('no error raised: %s', msg);

%% -------------------------------------------------------------------------
function test_info(extra)
try
    info = mmclab('gpuinfo');
catch
    return  % no OpenCL runtime
end
assert(isstruct(info) || isempty(info), 'gpuinfo must return a struct array');

function test_errors(extra)
cfg = mkcfg(extra);
assert_error(@() mmclab(rmfield(cfg, 'node')), 'missing node');
assert_error(@() mmclab(rmfield(cfg, 'elemprop')), 'missing elemprop');
assert_error(@() mmclab(rmfield(cfg, 'prop')), 'missing prop');
assert_error(@() mmclab(mkcfg(extra, 'method', 'xyz')), 'unknown method');
assert_error(@() mmclab(cfg, 'xyz'), 'unknown option');

function test_prep(extra)
% 'prep' fills in the element data and the initial element, without running
cfg = mkcfg(extra);
ncfg = mmclab(cfg, 'prep');
assert(isequal(size(ncfg.facenb), [6 4]) && isequal(numel(ncfg.evol), 6), 'facenb/evol');
assert(all(ncfg.evol > 0), 'positive element volumes');
assert(ncfg.e0 >= 1 && ncfg.e0 <= 6, 'the source element');
% a prepared cfg runs as it is
flux = mmclab(ncfg);
assert(isequal(size(flux.data), [21 21 21]), 'running a prepared cfg');

function test_outputs(extra)
% 'grid' gives voxels covering the mesh; the mesh tracers give nodal or elemental values
flux = mmclab(mkcfg(extra));
assert(isequal(size(flux.data), [21 21 21]), sprintf('grid output size %s', mat2str(size(flux.data))));
flux = mmclab(mkcfg(extra, 'tstep', 1e-9));
assert(isequal(size(flux.data), [21 21 21 5]), 'time gates');
flux = mmclab(mkcfg(extra, 'method', 'elem', 'issaveref', 1));
assert(isequal(size(flux.data), [8 1]), 'nodal output');
assert(isequal(size(flux.dref), [12 1]) && sum(flux.dref(:)) > 0, 'reflectance on the 12 surface triangles');
flux = mmclab(mkcfg(extra, 'method', 'elem', 'basisorder', 0));
assert(isequal(size(flux.data), [6 1]), 'elemental output');

function test_detected_photons(extra)
cfg = mkcfg(extra, 'detpos', [10 8 0 2], 'ismomentum', 1, 'issaveexit', 1, 'nphoton', 1e4);
[flux, detp, ncfg, seeds, traj] = mmclab(cfg);
ndet = numel(detp.detid);
assert(ndet > 0, 'no photon detected');
assert(all(detp.detid == 1), 'detector id');
assert(isequal(size(detp.ppath), [ndet 1]) && all(detp.ppath > 0), 'ppath: one column per medium');
assert(isequal(size(detp.p), [ndet 3]) && isequal(size(detp.v), [ndet 3]), 'exit position and direction');
assert(isequal(size(detp.mom), [ndet 1]), 'momentum transfer');
assert(isa(seeds.data, 'uint8') && size(seeds.data, 2) == ndet, 'one seed column per detected photon');
assert(all(isfield(traj, {'id', 'pos', 'eid', 'data'})) && size(traj.pos, 2) == 3, 'trajectory fields');
assert(issorted(double(traj.id)), 'trajectories sorted by photon id');
assert(isfield(ncfg, 'facenb') && isfield(ncfg, 'e0'), 'the 3rd output is the prepared cfg');

function test_struct_array(extra)
% cfgs(1).isreflect is empty, which keeps the default (1); it was read as 0
cfgs = mkcfg(extra, 'prop', [0 0 1 1; 0.05 2 0.8 1.37], 'outputtype', 'energy', 'nphoton', 2e4);
cfgs(2) = cfgs(1);
cfgs(2).isreflect = 0;
flux = mmclab(cfgs);
assert(numel(flux) == 2, 'one output per cfg element');
% Fresnel reflection keeps more light inside the mesh
assert(energy(flux(1)) > energy(flux(2)) + 0.05, sprintf('isreflect of each element: %.4f %.4f', energy(flux(1)), energy(flux(2))));

function test_rf(extra)
flux = mmclab(mkcfg(extra, 'method', 'elem', 'omega', 2 * pi * 1e8));
assert(~isreal(flux.data) && isequal(size(flux.data), [8 1]), 'complex RF output');

function test_replay(extra)
cfg = mkcfg(extra, 'detpos', [10 8 0 2], 'nphoton', 2e4);
[flux, detp, ncfg, seeds] = mmclab(cfg);
cfg.seed = seeds.data;
cfg.detphotons = detp.data;
cfg.outputtype = 'jacobian';
jac = mmclab(cfg);
assert(all(isfinite(jac.data(:))) && sum(abs(jac.data(:))) > 0, 'Jacobian values');

function test_option_sets_compute(extra)
% mmclab(cfg, 'sse') runs on the CPU, overriding cfg.compute
cfg = mkcfg(extra, 'outputtype', 'energy', 'compute', 'xyz');
assert_error(@() mmclab(cfg), 'unknown cfg.compute');
flux = mmclab(cfg, 'sse');
assert(energy(flux) > 0.5, 'option sse');

function test_widefield_source(extra)
% a planar source outside the mesh gets a source domain (label -1) added
src = {'srctype', 'planar', 'srcpos', [5 5 -2], 'srcparam1', [10 0 0 0], 'srcparam2', [0 10 0 0]};
cfg = boxcfg(extra, src{:});
[flux, detp, ncfg] = mmclab(cfg);
assert(any(ncfg.elemprop < 0), 'no source domain added');
assert(size(ncfg.elem, 1) > size(cfg.elem, 1), 'the source domain adds elements');
assert(sum(flux.data(:)) > 0, 'no light from the wide-field source');
% a planar source on the mesh surface leaves no room for a source domain: an error
src{4} = [5 5 0];
assert_error(@() mmclab(boxcfg(extra, src{:})), 'a planar source on the surface');

function test_flog_silent(extra)
% flog set to the null device prints nothing, not even the cfg echo
if (~exist('evalc'))
    return
end
cfg = mkcfg(extra, 'nphoton', 200);
out = evalc('mmclab(cfg);');
assert(isempty(strtrim(out)), ['unexpected output: ' out]);

function test_end_to_end(extra)
% the absorbed fraction of the cube: 0.574 (+/-0.005) from mcx with the same medium;
% the same energy whether accumulated on the grid or in the elements
grid = energy(mmclab(mkcfg(extra, 'outputtype', 'energy', 'nphoton', 2e4)));
elem = energy(mmclab(mkcfg(extra, 'outputtype', 'energy', 'nphoton', 2e4, 'method', 'elem', 'basisorder', 0)));
assert(grid > 0.54 && grid < 0.60, sprintf('absorbed fraction %.4f', grid));
assert(abs(grid - elem) < 0.005, sprintf('grid %.4f vs elements %.4f', grid, elem));
