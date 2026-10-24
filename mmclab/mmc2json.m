function varargout = mmc2json(cfg, filestub, varargin)
%
% Format:
%    mmc2json(cfg,filestub)
%
% Save MCXLAB simulation configuration to a JSON file for MCX binary
%
% Author: Qianqian Fang <q.fang at neu.edu>
%
% Input:
%    cfg: a struct defining the parameters associated with a simulation.
%         Please run 'help mcxlab' or 'help mmclab' to see the details.
%         mcxpreview supports the cfg input for both mcxlab and mmclab.
%    filestub: the filestub is the name stub for all output files,including
%         filestub.json: the JSON input file
%
%         if filestub ends with '.json', then mmc2json saves the mesh data
%         in the single-file mode, and the mesh information is stored
%         inside the "MeshNode" and "MeshElem" JSON keys inside the "Mesh"
%         root-level object.
%
% Dependency:
%    this function depends on the savejson/saveubjson functions from the
%    Iso2Mesh toolbox (http://iso2mesh.sf.net) or JSONlab toolbox
%    (http://iso2mesh.sf.net/jsonlab)
%
% This function is part of Monte Carlo eXtreme (MCX) URL: https://mcx.space
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

singlefile = (nargin == 1 || ~isempty(regexp(filestub, '\.json$', 'once')));

if (nargin > 1)
    [fpath, fname, fext] = fileparts(filestub);
    filestub = fullfile(fpath, fname);
end

%% define the optodes: sources and detectors

Optode.Source = struct();
Optode.Source = copycfg(cfg, 'srcpos', Optode.Source, 'Pos');
Optode.Source = copycfg(cfg, 'srcdir', Optode.Source, 'Dir');
Optode.Source = copycfg(cfg, 'srcparam1', Optode.Source, 'Param1');
Optode.Source = copycfg(cfg, 'srcparam2', Optode.Source, 'Param2');
Optode.Source = copycfg(cfg, 'srctype', Optode.Source, 'Type');
Optode.Source = copycfg(cfg, 'srcnum', Optode.Source, 'SrcNum');
Optode.Source = copycfg(cfg, 'srcid', Optode.Source, 'ID');

if (isfield(cfg, 'detpos') && ~isempty(cfg.detpos))
    Optode.Detector = struct();
    Optode.Detector = cell2struct(mat2cell(cfg.detpos, ones(1, size(cfg.detpos, 1)), [3 1]), {'Pos', 'R'}, 2);
    Optode.Detector = Optode.Detector(:)';
    if (length(Optode.Detector) == 1)
        Optode.Detector = {Optode.Detector};
    end
end
if (isfield(cfg, 'srcpattern') && ~isempty(cfg.srcpattern))
    % a row-major JData array, [srcnum, Nx, Ny] for multiple patterns; mmc only reads
    % the compressed form, and takes the pattern count from its first dimension
    Optode.Source.Pattern = jdataencode(single(cfg.srcpattern), 'AnnotateArray', 1, 'Compression', 'zlib', 'CompressArraySize', 0, 'Base64', 1);
end

%% define the domain and optical properties

Mesh = struct();
Mesh = copycfg(cfg, 'unitinmm', Mesh, 'LengthUnit');

if (isfield(cfg, 'node') && ~isempty(cfg.node) && isfield(cfg, 'elem') && ~isempty(cfg.elem))
    node = cfg.node;
    elem = cfg.elem;
    if (isfield(cfg, 'elemprop') && size(elem, 2) == 4)
        elem = [elem, cfg.elemprop];
    end
    if (~singlefile)
        Mesh.MeshID = fname;
        savemmcmesh(fname, node, elem);
    else
        Mesh.MeshNode = single(node);
        Mesh.MeshElem = int32(elem);  % signed: uint32 clamps the -1/-2 labels of wide-field source/detector domains to 0
        if (isfield(cfg, 'edgeroi'))
            Mesh.MeshROI = single(cfg.edgeroi);
        elseif (isfield(cfg, 'noderoi'))
            Mesh.MeshROI = single(cfg.noderoi);
        elseif (isfield(cfg, 'faceroi'))
            Mesh.MeshROI = single(cfg.faceroi);
        end
    end

    if (~isfield(cfg, 'e0'))
        cfg.e0 = tsearchn(node, elem(:, 1:4), cfg.srcpos(1, 1:3));  % mmc locates any other sources itself
    end
    Mesh = copycfg(cfg, 'e0', Mesh, 'InitElem');

else
    Mesh.MeshID = fname;
    warning('mesh is missing!');
end

Domain = struct();
Domain = copycfg(cfg, 'steps', Domain, 'Step');
if (isfield(cfg, 'prop'))
    prop = [(1:size(cfg.prop, 1) - 1)' cfg.prop(2:end, :)];
    if (~singlefile)
        fid = fopen(['prop_', fname, '.dat'], 'wt');
        fprintf(fid, '1 %d\n', size(prop, 1));
        fprintf(fid, '%d %e %e %e %e\n', prop');
        fclose(fid);
    else
        Domain.Media = cell2struct(num2cell(cfg.prop), {'mua', 'mus', 'g', 'n'}, 2)';
    end
else
    warning('prop is missing');
end

%% define the simulation session flags

Session = struct();
if (exist('fname', 'var'))
    Session.ID = fname;
else
    Session.ID = inputname(1);
end
Session = copycfg(cfg, 'isreflect', Session, 'DoMismatch');
Session = copycfg(cfg, 'issave2pt', Session, 'DoSaveVolume');
Session = copycfg(cfg, 'issavedet', Session, 'DoPartialPath');
Session = copycfg(cfg, 'issaveexit', Session, 'DoSaveExit');
Session = copycfg(cfg, 'issaveseed', Session, 'DoSaveSeed');
Session = copycfg(cfg, 'isnormalized', Session, 'DoNormalize');
Session = copycfg(cfg, 'ismomentum', Session, 'DoDCS');
Session = copycfg(cfg, 'isspecular', Session, 'DoSpecular');
Session = copycfg(cfg, 'outputformat', Session, 'OutputFormat');
Session = copycfg(cfg, 'debuglevel', Session, 'DebugFlag');
Session = copycfg(cfg, 'autopilot', Session, 'DoAutoThread');
Session = copycfg(cfg, 'basisorder', Session, 'BasisOrder');

if (isfield(cfg, 'method'))
    methodmap = struct('plucker', 'p', 'havel', 'h', 'badouel', 'b', 'elem', 's', 'grid', 'g');
    if (isfield(methodmap, cfg.method))
        Session.RayTracer = methodmap.(cfg.method);
    else
        error('cfg.method is invalid');
    end
end

if (isfield(cfg, 'outputtype'))
    outputtypemap = struct('flux', 'x', 'fluence', 'f', 'energy', 'e', 'jacobian', 'j', 'wl', 'l', 'wp', 'p', ...
                           'rf', 'r', 'rfmus', 's', 'adjoint', 'a', 'adjoint_dcoeff', 'd', 'adjoint_mus', 'u', ...
                           'adjoint_musp', 'v', 'adjoint_mua_d', 'w', 'adjoint_mua_musp', 'q');
    if (isfield(outputtypemap, cfg.outputtype))
        Session.OutputType = outputtypemap.(cfg.outputtype);
    else
        error('cfg.outputtype is invalid');
    end
end

if (isfield(cfg, 'seed') && numel(cfg.seed) == 1)
    Session.RNGSeed = cfg.seed;
end
Session = copycfg(cfg, 'nphoton', Session, 'Photons');
% Session=copycfg(cfg,'rootpath',Session,'RootPath');

%% define the forward simulation settings

Forward.T0 = cfg.tstart;
Forward.T1 = cfg.tend;
Forward.Dt = cfg.tstep;
if (isfield(cfg, 'omega') && ~isempty(cfg.omega) && cfg.omega > 0)
    Forward.Omega = cfg.omega;  % modulation angular frequency (rad/s) of RF simulations
end
Forward = copycfg(cfg, 'nout', Forward, 'N0');

%% assemble the complete input, save to a JSON or UBJSON input file

mmcsession = struct('Session', Session, 'Domain', Domain, 'Mesh', Mesh, 'Forward', Forward, 'Optode', Optode);

if (nargin == 1)
    [varargout{1:nargout}] = savejson('', mmcsession, varargin{:});
elseif (strcmp(fext, 'ubj'))
    [varargout{1:nargout}] = saveubjson('', mmcsession, 'filename', [filestub, '.ubj'], varargin{:});
else
    [varargout{1:nargout}] = savejson('', mmcsession, 'filename', [filestub, '.json'], varargin{:});
end

function outdata = copycfg(cfg, name, outroot, outfield, defaultval)
if (nargin >= 5)
    outroot.(outfield) = defaultval;
end
if (isfield(cfg, name))
    outroot.(outfield) = cfg.(name);
end
outdata = outroot;
