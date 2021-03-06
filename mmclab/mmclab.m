function varargout=mmclab(varargin)
%
%#############################################################################%
%         MMCLAB - Mesh-based Monte Carlo (MMC) for MATLAB/GNU Octave         %
%          Copyright (c) 2010-2020 Qianqian Fang <q.fang at neu.edu>          %
%                            http://mcx.space/#mmc                            %
%                                                                             %
% Computational Optics & Translational Imaging (COTI) Lab- http://fanglab.org %
%            Department of Bioengineering, Northeastern University            %
%                                                                             %
%               Research funded by NIH/NIGMS grant R01-GM114365               %
%#############################################################################%
%$Rev::      $v2020 $Date::                       $ by $Author::             $%
%#############################################################################%
%
% Format:
%    [fluence,detphoton,ncfg,seeds]=mmclab(cfg);
%          or
%    fluence=mmclab(cfg);
%    newcfg=mmclab(cfg,'prep');
%    [fluence,detphoton,ncfg,seeds]=mmclab(cfg, options);
%
% Input:
%    cfg: a struct, or struct array. Each element in cfg defines 
%         a set of parameters for a simulation. 
%
%    option: (optional), options is a string, specifying additional options
%         option='preview': this plots the domain configuration using mcxpreview(cfg)
%         option='opencl':  force using OpenCL (set cfg.gpuid=1 if not set)
%                           instead of SSE on CPUs/GPUs that support OpenCL
%
%
%    cfg may contain the following fields:
%
%== Required ==
%     *cfg.nphoton:     the total number of photons to be simulated (integer)
%     *cfg.prop:        an N by 4 array, each row specifies [mua, mus, g, n] in order.
%                       the first row corresponds to medium type 0 which is 
%                       typically [0 0 1 1]. The second row is type 1, and so on.
%     *cfg.node:        node array for the input tetrahedral mesh, 3 columns: (x,y,z)
%     *cfg.elem:        element array for the input tetrahedral mesh, 4 columns
%     *cfg.elemprop:    element property index for input tetrahedral mesh
%     *cfg.tstart:      starting time of the simulation (in seconds)
%     *cfg.tstep:       time-gate width of the simulation (in seconds)
%     *cfg.tend:        ending time of the simulation (in second)
%     *cfg.srcpos:      a 1 by 3 vector, the position of the source in mesh node length unit
%     *cfg.srcdir:      if defined as [vx, vy, vy], it specifies the incident vector
%                       if defined as [vx, vy, vy, focus], the first 3 elements define
%                       the incident vector; focus controls the convergence or 
%                       divergence of the beam:
%                       focus=0: collimated beam
%                       focus<0: diverging beam from an imaginary src at c0-|focus|*[vx vy vz]
%                       focus>0: converging beam, focusing to a point at c0+|focus|*[vx vy vz]
%                       where c0 is the centroid of the source domain. Setting focus does 
%                       not impact pencil/isotropic/cone sources.
%
%== MC simulation settings ==
%      cfg.seed:        seed for the random number generator (integer)
%                       if set to a uint8 array, the binary data in each column is used 
%                       to seed a photon (i.e. the "replay" mode), default value: 1648335518
%      cfg.isreflect:   [1]-consider refractive index mismatch, 0-matched index
%                       2 - total absorption on exterior surface
%                       3 - prefect reflection (mirror) on exterior surface
%      cfg.isnormalized:[1]-normalize the output fluence to unitary source, 0-no reflection
%      cfg.isspecular:  [1]-calculate specular reflection if source is outside
%      cfg.ismomentum:  [0]-save momentum transfer for each detected photon
%      cfg.method:      ray-tracing method, ["plucker"]:Plucker, "havel": Havel (SSE4),
%                       "badouel": partial Badouel, "elem": branchless Badouel (SSE), 
%                       "grid": dual-grid MMC
%      cfg.mcmethod:    0 use MCX-styled MC method, 1 use MCML style MC
%      cfg.nout:        [1.0] refractive index for medium type 0 (background)
%      cfg.minenergy:   terminate photon when weight less than this level (float) [0.0]
%      cfg.roulettesize:[10] size of Russian roulette
%      cfg.unitinmm:    defines the default length unit (to interpret mesh nodes, src/det positions 
%                       the default value is 1.0 (mm). For example, if the mesh node length unit is 
%                       in cm, one should set unitinmm to 10.
%      cfg.basisorder:  [1]-linear basis, 0-piece-wise constant basis
%
%== Source-detector parameters ==
%      cfg.detpos:      an N by 4 array, each row specifying a detector: [x,y,z,radius]
%      cfg.srctype:     source type, the parameters of the src are specified by cfg.srcparam{1,2}
%                      'pencil' - default, pencil beam, no param needed
%                      'isotropic' - isotropic source, no param needed
%                      'cone' - uniform cone beam, srcparam1(1) is the half-angle in radian
%                      'gaussian' - a gaussian beam, srcparam1(1) specifies the waist radius 
%                                (in default length unit); if one specifies a non-zero focal length
%                                using cfg.srcdir, the gaussian beam can be converging to or 
%                                diverging from the waist center, which is located at srcpos+focus*srcdir;
%                                optionally, one can specify the wavelength lambda (in cfg.unitinmm mm), 
%                                using srcparam1(2). This will rescale the Gaussian profile according 
%                                to w(z)=w0*sqrt(1-(z/z0)^2), where w0 is the waist radius, z is the 
%                                distance (in mm) to the waist center (focus), and z0 is the Rayleigh 
%                                range (in mm), and z0 is related to w0 by z0=w0^2*pi/lambda
%                      'planar' - a 3D quadrilateral uniform planar source, with three corners specified 
%                                by srcpos, srcpos+srcparam1(1:3) and srcpos+srcparam2(1:3)
%                      'pattern' - a 3D quadrilateral pattern illumination, same as above, except
%                                srcparam1(4) and srcparam2(4) specify the pattern array x/y dimensions,
%                                and srcpattern is a floating-point pattern array, with values between [0-1]. 
%                                if cfg.srcnum>1, srcpattern must be a floating-point array with 
%                                a dimension of [srcnum srcparam1(4) srcparam2(4)]
%                                Example: <demo_photon_sharing.m>
%                      'fourier' - spatial frequency domain source, similar to 'planar', except
%                                the integer parts of srcparam1(4) and srcparam2(4) represent
%                                the x/y frequencies; the fraction part of srcparam1(4) multiplies
%                                2*pi represents the phase shift (phi0); 1.0 minus the fraction part of
%                                srcparam2(4) is the modulation depth (M). Put in equations:
%                                    S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
%                      'arcsine' - similar to isotropic, except the zenith angle is uniform
%                                distribution, rather than a sine distribution.
%                      'disk' - a uniform disk source pointing along srcdir; the radius is 
%                               set by srcparam1(1) (in default length unit)
%                      'fourierx' - a general Fourier source, the parameters are 
%                               srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phi0,M]
%                               normalized vectors satisfy: srcdir cross v1=v2
%                               the phase shift is phi0*2*pi
%                      'fourierx2d' - a general 2D Fourier basis, parameters
%                               srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phix,phiy]
%                               the phase shift is phi{x,y}*2*pi
%                      'zgaussian' - an angular gaussian beam, srcparam1(1) specifies the variance in  
%                               the zenith angle
%      cfg.{srcparam1,srcparam2}: 1x4 vectors, see cfg.srctype for details
%      cfg.srcpattern: see cfg.srctype for details
%      cfg.srcnum:     the number of source patterns that are
%                      simultaneously simulated; only works for 'pattern'
%                      source, see cfg.srctype='pattern' for details
%                      Example <demo_photon_sharing.m>
%      cfg.replaydet:  only works when cfg.outputtype is 'jacobian', 'wl', 'nscat', or 'wp' and cfg.seed is an array
%                       0 replay all detectors and sum all Jacobians into one volume
%                       a positive number: the index of the detector to replay and obtain Jacobians
%      cfg.voidtime:   for wide-field sources, [1]-start timer at launch, 0-when entering
%                      the first non-zero voxel
%
%      by default, mmc assumes the mesh and source position settings are all in mm unit.
%      if the mesh coordinates/source positions are not in mm unit, one needs to define
%      cfg.unitinmm  (in mm) to specify the actual length unit.
%
%== Optional mesh data ==
%     -cfg.facenb:      element face neighbohood list (calculated by faceneighbors())
%     -cfg.evol:        element volume (calculated by elemvolume() with iso2mesh)
%     -cfg.e0:          the element ID enclosing the source, if not defined,
%                       it will be calculated by tsearchn(node,elem,srcpos);
%                       if cfg.e0 is set as one of the following characters,
%                       mmclab will do an initial ray-tracing and move
%                       srcpos to the first intersection to the surface:
%                       '>': search along the forward (srcdir) direction
%                       '<': search along the backward direction
%                       '-': search both directions
%
%== Output control ==
%      cfg.issaveexit: [0]-save the position (x,y,z) and (vx,vy,vz) for a detected photon
%      cfg.issaveref:  [0]-save diffuse reflectance/transmittance on the exterior surfaces.
%                      The output is stored as flux.dref in a 2D array of size [#Nf,  #time_gate]
%                      where #Nf is the number of triangles on the surface; #time_gate is the
%                      number of total time gates. To plot the surface diffuse reflectance, the output
%                      triangle surface mesh can be extracted by faces=faceneighbors(cfg.elem,'rowmajor');
%                      where 'faceneighbors' can be found in the iso2mesh toolbox.
%                      Example: see <demo_mmclab_basic.m>
%      cfg.issaveseed:  [0]-save the RNG seed for a detected photon so one can replay
%      cfg.isatomic:    [1]-use atomic operations for saving fluence, 0-no atomic operations
%      cfg.outputtype:  'flux' - output fluence-rate
%                       'fluence' - fluence, 
%                       'energy' - energy deposit, 
%                       'jacobian' - mua Jacobian (replay mode)
%                       'wl'- weighted path lengths to build mua Jacobian (replay mode)
%                       'wp'- weighted scattering counts to build mus Jacobian (replay mode)
%      cfg.debuglevel:  debug flag string, a subset of [MCBWDIOXATRPE], no space
%      cfg.debugphoton: print the photon movement debug info only for a specified photon ID
%
%      fields marked with * are required; options in [] are the default values
%      fields marked with - are calculated if not given (can be faster if precomputed)
%
%
%    type: omit or 'omp' for multi-threading version; 'sse' for the SSE4 MMC,
%          the SSE4 version is about 25% faster, but requires newer CPUs; 
%          if type='prep' with a single output, mmclab returns ncfg only.
%
% Output:
%      fluence: a struct array, with a length equals to that of cfg.
%            For each element of fluence, fluence(i).data is a 2D array with
%            dimensions [size(cfg.node,1), total-time-gates] if cfg.basisorder=1,
%            or [size(cfg.elem,1), total-time-gates] if cfg.basisorder=0. 
%            The content of the array is the normalized fluence-rate (or others 
%            depending on cfg.outputtype) at each mesh node and time-gate.
%
%            If cfg.issaveref is set to 1, fluence(i).dref is not empty, and stores
%            the surface diffuse reflectance (normalized by default). The surface mesh
%            that the dref output is attached can be obtained by faces=faceneighbors(cfg.elem,'rowmajor');
%      detphoton: (optional) a struct array, with a length equals to that of cfg.
%            Starting from v2016.5, the detphoton contains the below subfields:
%              detphoton.detid: the ID(>0) of the detector that captures the photon
%              detphoton.nscat: cummulative scattering event counts in each medium
%              detphoton.ppath: cummulative path lengths in each medium (partial pathlength)
%                   one need to multiply cfg.unitinmm with ppath to convert it to mm.
%              detphoton.mom: cummulative cos_theta for momentum transfer in each medium
%              detphoton.p or .v: exit position and direction, when cfg.issaveexit=1
%              detphoton.w0: photon initial weight at launch time
%              detphoton.prop: optical properties, a copy of cfg.prop
%              detphoton.data: a concatenated and transposed array in the order of
%                    [detid nscat ppath mom p v w0]'
%              "data" is the is the only subfield in all MMCLAB before 2016.5
%      ncfg: (optional), if given, mmclab returns the preprocessed cfg structure,
%            including the calculated subfields (marked by "-"). This can be
%            used in the subsequent simulations to avoid repetitive preprocessing.
%      seeds: (optional), if give, mmclab returns the seeds, in the form of
%            a byte array (uint8) for each detected photon. The column number
%            of seed equals that of detphoton.
%
% Example:
%      cfg.nphoton=1e5;
%      [cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);
%      cfg.elemprop=ones(size(cfg.elem,1),1);
%      cfg.srcpos=[30 30 0];
%      cfg.srcdir=[0 0 1];
%      cfg.prop=[0 0 1 1;0.005 1 0 1.37];
%      cfg.tstart=0;
%      cfg.tend=5e-9;
%      cfg.tstep=5e-10;
%      cfg.debuglevel='TP';
%      % preprocessing to populate the missing fields to save computation
%      ncfg=mmclab(cfg,'prep');
%
%      cfgs(1)=ncfg;   % when using struct array input, all fields must be defined
%      cfgs(2)=ncfg;
%      cfgs(1).isreflect=0;
%      cfgs(2).isreflect=1;
%      cfgs(2).detpos=[30 20 0 1;30 40 0 1;20 30 1 1;40 30 0 1];
%      % calculate the fluence and partial path lengths for the two configurations
%      [fluxs,detps]=mmclab(cfgs);
%
%
% This function is part of Mesh-based Monte Carlo (MMC) URL: http://mcx.space/#mmc
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

try
    defaultocl=evalin('base','USE_MCXCL');
catch
    defaultocl=1;
end

useopencl=defaultocl;

if(nargin==2 && ischar(varargin{2}))
    if(strcmp(varargin{2},'preview') || strcmp(varargin{2},'prep') || strcmp(varargin{2},'cuda'))
        useopencl=0;
    end
end

if(isstruct(varargin{1}))
    for i=1:length(varargin{1})
        castlist={'srcpattern','srcpos','detpos','prop','workload','srcdir'};
        for j=1:length(castlist)
            if(isfield(varargin{1}(i),castlist{j}))
                varargin{1}(i).(castlist{j})=double(varargin{1}(i).(castlist{j}));
            end
        end
    end
end

if(nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'gpuinfo'))
    varargout{1}=mmc('gpuinfo');
    return;
end

if(nargin==0)
    return;
end

cfg=varargin{1};
if(length(varargin)>=2)
    type=varargin{2};
end
    
if(~isstruct(cfg))
    error('cfg must be a struct or struct array');
end

len=length(cfg);
for i=1:len
    if(~isfield(cfg(i),'node') || ~isfield(cfg(i),'elem'))
        error('cfg.node or cfg.elem is missing');
    end
    if(size(cfg(i).elem,2)>4)
        cfg(i).elemprop=cfg(i).elem(:,5);
    end
    if(~isfield(cfg(i),'isreoriented') || isempty(cfg(i).isreoriented) || cfg(i).isreoriented==0)
        [cfg(i).elem, evol, idx]=meshreorient(cfg(i).node,cfg(i).elem(:,1:4));
	if(isfield(cfg(i),'edgeroi'))
	    cfg(i).edgeroi(idx,:)=cfg(i).edgeroi(idx,[1 3 2 5 4 6]);
	end
	if(isfield(cfg(i),'faceroi'))
	    cfg(i).faceroi(idx,:)=cfg(i).faceroi(idx,[2 1 3 4]);
	end
        cfg(i).isreoriented=1;
    end
    if(~isfield(cfg(i),'facenb') || isempty(cfg(i).facenb))
        cfg(i).facenb=faceneighbors(cfg(i).elem);
    end
    if(~isfield(cfg(i),'evol') || isempty(cfg(i).evol))
        cfg(i).evol=elemvolume(cfg(i).node,cfg(i).elem);
    end
    if(find(cfg(i).evol==0))
        fprintf(1,['degenerated elements are detected: [' sprintf('%d ',find(cfg(i).evol==0)) ']\n']);
        error(['input mesh can not contain degenerated elements, ' ...
            'please double check your input mesh; if you use a ' ...
            'widefield source, please rerun mmcsrcdomain and setting ' ...
            '''Expansion'' option to a larger value (default is 1)']);
    end
    if(~isfield(cfg(i),'srcpos'))
        error('cfg.srcpos field is missing');
    end
    if(~isfield(cfg(i),'srcdir'))
        error('cfg.srcdir field is missing');
    end
    if(~isfield(cfg(i),'e0') || isempty(cfg(i).e0))
        cfg(i).e0=tsearchn(cfg(i).node,cfg(i).elem,cfg(i).srcpos);
    end
    if((isnan(cfg(i).e0) && (isfield(cfg(i),'srctype') ...
            && strcmp(cfg(i).srctype,'pencil')) )|| ischar(cfg(i).e0))
        disp('searching initial element ...');
        [cfg(i).srcpos,cfg(i).e0]=mmcraytrace(cfg(i).node,cfg(i).elem,cfg(i).srcpos,cfg(i).srcdir,cfg(i).e0);
    end
    if((isfield(cfg(i),'srctype') && strcmp(cfg(i).srctype,'pattern')) && (ndims(cfg(i).srcpattern)==2))
        cfg(i).srcpattern=reshape(cfg(i).srcpattern,...
            size(cfg(i).srcpattern,1),size(cfg(i).srcpattern,2),1);
    end
    if(isnan(cfg(i).e0))  % widefield source
        if(~isfield(cfg(i),'srcparam1') || ~isfield(cfg(i),'srcparam2'))
            error('for wide-field sources, you must provide srcparam1 and srcparam2');
        end
        if(~isfield(cfg(i),'srctype'))
            cfg(i).srctype='pencil';
        end
        srcdef=struct('srctype',cfg.srctype,'srcpos',cfg.srcpos,'srcdir',cfg.srcdir,...
            'srcparam1',cfg.srcparam1,'srcparam2',cfg.srcparam2);
        sdom=mmcsrcdomain(srcdef,[min(cfg.node);max(cfg.node)]);
        isinside=ismember(round(sdom*1e10)*1e-10,round(cfg(i).node*1e10)*1e-10,'rows');
        if(all(~isinside))
            if(size(cfg(i).elem,2)==4)
                cfg(i).elem(:,5)=1;
            end
            [cfg(i).node,cfg(i).elem] = mmcaddsrc(cfg(i).node,cfg(i).elem,sdom);
            cfg(i).elemprop=cfg(i).elem(:,5);
            [cfg(i).elem, evol, idx]=meshreorient(cfg(i).node,cfg(i).elem(:,1:4));
            if(isfield(cfg(i),'edgeroi'))
	        cfg(i).edgeroi(idx,:)=cfg(i).edgeroi(idx,[1 3 2 5 4 6]);
	    end
	    if(isfield(cfg(i),'faceroi'))
	        cfg(i).faceroi(idx,:)=cfg(i).faceroi(idx,[2 1 3 4]);
	    end
            cfg(i).facenb=faceneighbors(cfg(i).elem);
            cfg(i).evol=elemvolume(cfg(i).node,cfg(i).elem);
            cfg(i).isreoriented=1;
        end
        if(strcmp(cfg(i).srctype,'pencil') || strcmp(cfg(i).srctype,'isotropic') ...
                || strcmp(cfg(i).srctype,'cone')   || strcmp(cfg(i).srctype,'zgaussian'))
            cfg(i).e0=tsearchn(cfg(i).node,cfg(i).elem,cfg(i).srcpos);
            if(isnan(cfg(i).e0))
                cfg(i).e0=-1;
            end
        elseif(length(find(cfg(i).elemprop<0))==1)
            cfg(i).e0=find(cfg(i).elemprop<0);
        else
            cfg(i).e0=-1;
        end
    end
    if(~isfield(cfg(i),'elemprop'))
        error('cfg.elemprop field is missing');
    end
    if(~isfield(cfg(i),'nphoton'))
        error('cfg.nphoton field is missing');
    end
    if(~isfield(cfg(i),'prop') || size(cfg(i).prop,1)<max(cfg(i).elemprop)+1 || min(cfg(i).elemprop<=0))
        error('cfg.prop field is missing or insufficient');
    end
end

% must do a fflush, otherwise octave buffers the output until complete
if(exist('OCTAVE_VERSION'))
    fflush(stdout);
end

if(exist('maxNumCompThreads','file'))
    warning('off','MATLAB:maxNumCompThreads:Deprecated');
    maxNumCompThreads('automatic');
    maxNumCompThreads(maxNumCompThreads*2);
end

mmcout=nargout;
if(nargout>=3)
    mmcout=nargout-1;
    varargout{nargout}=cfg;
end

if(useopencl==1)
    if(isfield(cfg,'gpuid') && ~ischar(cfg.gpuid) && cfg.gpuid<-1)
	    cfg.gpuid=1;
    end
    [varargout{1:mmcout}]=mmc(cfg);
elseif(length(varargin)<2)
    [varargout{1:mmcout}]=mmc(cfg);
elseif(strcmp(type,'omp'))
    [varargout{1:mmcout}]=mmc(cfg);
elseif(strcmp(type,'sse'))
    [varargout{1:mmcout}]=mmc_sse(cfg);
elseif(strcmp(type,'prep') && nargout==1)
    varargout{1}=cfg;
elseif(strcmp(type,'preview') && nargout==1)
    varargout{1}=mcxpreview(cfg);
else
    error('type is not recognized');
end

if(mmcout>=2)
    for i=1:length(varargout{2})
        if(~isfield(cfg(i),'issaveexit') || cfg(i).issaveexit~=2)
            medianum=size(cfg(i).prop,1)-1;
            detp=varargout{2}(i).data;
            if(isempty(detp))
                continue;
            end
            newdetp.detid=int32(detp(1,:))';
            newdetp.nscat=int32(detp(2:medianum+1,:))';    % 1st medianum block is num of scattering
            newdetp.ppath=detp(medianum+2:2*medianum+1,:)';% 2nd medianum block is partial path
            if(isfield(cfg(i),'ismomentum') && cfg(i).ismomentum)
                newdetp.mom=detp(2*medianum+2:3*medianum+1,:)'; % 3rd medianum block is the momentum transfer
            end
            if(isfield(cfg(i),'issaveexit') && cfg(i).issaveexit)
                newdetp.p=detp(end-6:end-4,:)';             %columns 7-5 from the right store the exit positions*/
                newdetp.v=detp(end-3:end-1,:)';	     %columns 4-2 from the right store the exit dirs*/
            end
            newdetp.w0=detp(end,:)';  % last column is the initial packet weight
            newdetp.prop=cfg(i).prop;
            if(isfield(cfg(i),'unitinmm'))
                newdetp.unitinmm=cfg(i).unitinmm;
            end
            newdetp.data=detp;      % enable this line for compatibility
            newdetpstruct(i)=newdetp;
        else
            newdetpstruct(i)=varargout{2}(i);
        end
    end
    if(exist('newdetpstruct','var'))
        varargout{2}=newdetpstruct;
    end
end

if(nargout>=4)
    [varargout{3:end}]=deal(varargout{[end 3:end-1]});
end

