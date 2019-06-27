function mmc2json(cfg,filestub)
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
%         filestub_vol.bin: the volume file if cfg.vol is defined
%         filestub_shapes.json: the domain shape file if cfg.shapes is defined
%         filestub_pattern.bin: the domain shape file if cfg.pattern is defined
%
% Dependency:
%    this function depends on the savejson/saveubjson functions from the 
%    Iso2Mesh toolbox (http://iso2mesh.sf.net) or JSONlab toolbox 
%    (http://iso2mesh.sf.net/jsonlab)
%
% This function is part of Monte Carlo eXtreme (MCX) URL: http://mcx.space
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

[fpath, fname, fext]=fileparts(filestub);
filestub=fullfile(fpath,fname);

%% define the optodes: sources and detectors

Optode.Source=struct();
Optode.Source=copycfg(cfg,'srcpos',Optode.Source,'Pos');
Optode.Source=copycfg(cfg,'srcdir',Optode.Source,'Dir');
Optode.Source=copycfg(cfg,'srcparam1',Optode.Source,'Param1');
Optode.Source=copycfg(cfg,'srcparam2',Optode.Source,'Param2');
Optode.Source=copycfg(cfg,'srctype',Optode.Source,'Type');
Optode.Source=copycfg(cfg,'srcnum',Optode.Source,'SrcNum');

if(isfield(cfg,'detpos') && ~isempty(cfg.detpos))
    Optode.Detector=struct();
    Optode.Detector=cell2struct(mat2cell(cfg.detpos, ones(1,size(cfg.detpos,1)),[3 1]), {'Pos','R'} ,2);
    if(length(Optode.Detector)==1)
        Optode.Detector={Optode.Detector};
    end
end
if(isfield(cfg,'srcpattern') && ~isempty(cfg.srcpattern))
    Optode.Source.Pattern.Nx=size(cfg.srcpattern,1);
    Optode.Source.Pattern.Ny=size(cfg.srcpattern,2);
    Optode.Source.Pattern.Nz=size(cfg.srcpattern,3);
    Optode.Source.Pattern.Data=[filestub '_pattern.bin'];
    fid=fopen(Optode.Source.Pattern.Data,'wb');
    fwrite(fid,cfg.srcpattern,'float32');
    fclose(fid);
end

%% define the domain and optical properties

Mesh=struct();
Mesh=copycfg(cfg,'unitinmm',Mesh,'LengthUnit');
Mesh.MeshID=fname;

if(isfield(cfg,'node') && ~isempty(cfg.node) && isfield(cfg,'elem') && ~isempty(cfg.elem))
    node=cfg.node;
    elem=cfg.elem;
    if(isfield(cfg,'elemprop') && size(elem,2)==4)
        elem=[elem, cfg.elemprop];
    end
    savemmcmesh(fname,node,elem);

    if(~isfield(cfg,'e0'))
        cfg.e0=tsearchn(node,elem(:,1:4),cfg.srcpos);
    end
    Mesh=copycfg(cfg,'e0',Mesh,'InitElem');

else
    warning('mesh is missing!')
end

if(isfield(cfg,'prop'))
    prop=[(1:size(cfg.prop,1)-1)' cfg.prop(2:end,:)];
    fid=fopen(['prop_',fname,'.dat'],'wt');
    fprintf(fid,'1 %d\n',size(prop,1));
    fprintf(fid,'%d %e %e %e %e\n',prop');
    fclose(fid);
else
    warning('prop is missing');
end

%% define the simulation session flags

Session=struct();
Session.ID=fname;
Session=copycfg(cfg,'isreflect',Session,'DoMismatch');
Session=copycfg(cfg,'issave2pt',Session,'DoSaveVolume');
Session=copycfg(cfg,'issavedet',Session,'DoPartialPath');
Session=copycfg(cfg,'issaveexit',Session,'DoSaveExit');
Session=copycfg(cfg,'issaveseed',Session,'DoSaveSeed');
Session=copycfg(cfg,'isnormalize',Session,'DoNormalize');
Session=copycfg(cfg,'ismomentum',Session,'DoDCS');
Session=copycfg(cfg,'DoSpecular',Session,'DoSpecular');
Session=copycfg(cfg,'outputformat',Session,'OutputFormat');
Session=copycfg(cfg,'outputtype',Session,'OutputType');
Session=copycfg(cfg,'debuglevel',Session,'DebugFlag');
Session=copycfg(cfg,'autopilot',Session,'DoAutoThread');
Session=copycfg(cfg,'basisorder',Session,'BasisOrder');
Session=copycfg(cfg,'method',Session,'RayTracer');

if(isfield(cfg,'seed') && numel(cfg.seed)==1)
    Session.RNGSeed=cfg.seed;
end
Session=copycfg(cfg,'nphoton',Session,'Photons');
%Session=copycfg(cfg,'rootpath',Session,'RootPath');

%% define the forward simulation settings

Forward.T0=cfg.tstart;
Forward.T1=cfg.tend;
Forward.Dt=cfg.tstep;
Forward=copycfg(cfg,'nout',Forward,'N0');

%% assemble the complete input, save to a JSON or UBJSON input file

mmcsession=struct('Session', Session, 'Mesh', Mesh, 'Forward', Forward, 'Optode',Optode);

if(strcmp(fext,'ubj'))
    saveubjson('',mmcsession,[filestub,'.ubj']);
else
    savejson('',mmcsession,[filestub,'.json']);
end


function outdata=copycfg(cfg,name,outroot,outfield,defaultval)
if(nargin>=5)
    outroot.(outfield)=defaultval;
end
if(isfield(cfg,name))
    outroot.(outfield)=cfg.(name);
end
outdata=outroot;
