%%-----------------------------------------------------------------
%% simple demonstration for debug flags
%%-----------------------------------------------------------------
%
% In this example, we run photon migration in a toy-problem to show
% the basic steps and behaviors of the simulation code.
% The mesh used in this case is a simple cube, splitted into
% 6 tetrahedra. The edge length of the cube is 10mm.
% A total of 20 photons are simulated starting
% from a source position [2,8,0] along an incident
% vector [0 0 1]. We run the simulation and print out
% the moving trajectories of the photon, and the exit and
% accumulation sites. Using a matlab script, one can visualize
% the photon moving from the output of the simulation.
%
%%-----------------------------------------------------------------

%%-----------------------------------------------------------------
%% defining mesh and input structure
%%-----------------------------------------------------------------

addpath('../../matlab/');

[cfg.node,cfg.elem]=genT6mesh(0:10:10,0:10:10,0:10:10);
cfg.elem=sortrows(cfg.elem);
cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.prop=[0 0 1 1;0.005 1.0101010101 0.01 1.0];
cfg.srcpos=[2,8,0];
cfg.srcdir=[0 0 1];
cfg.seed=1648335518;
cfg.nphoton=100;
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;

%%-----------------------------------------------------------------
%% run simulation and print out photon movement
%%-----------------------------------------------------------------

cfg.debuglevel='M';
flux=mmclab(cfg);

cfg.debuglevel='A';
flux=mmclab(cfg);


