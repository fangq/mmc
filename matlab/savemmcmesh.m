function savemmcmesh(key,node,elem,face,evol,facenb)
%
% savemmcmesh(key,node,elem,face,evol,facenb)
%
% Export a tetrahedral mesh to the MMC mesh format
%
% author: Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%     key: a string included in all exported mesh file names, the 
%          format of the files are {node,elem,face,facenb,evol}_key.dat
%      node: a node coordinate list, 3 columns for x/y/z
%      elem: a tetrahedral element list
%      face: a triangular surface face list, if missing, it will be calculated
%      evol: the volumes of all tetrahedra, if missing, it will be calculated
%      facenb: the 4 neighboring elements for each element, if missing, it will
%             be calculated
%
% example:
%     savemmcmesh('sph1',node,elem);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if(nargin<5)
  evol=elemvolume(node,elem(:,1:4));
end

if(~isempty(node))
  fid=fopen(['node_',key,'.dat'],'wt');
  fprintf(fid,'%d\t%d\n',1,size(node,1));
  fprintf(fid,'%d\t%16.8e\t%16.8e\t%16.8e\n',[1:length(node);node']);
  fclose(fid);
end

if(~isempty(elem))
  elem(:,1:4)=meshreorient(node,elem(:,1:4));

  fid=fopen(['elem_',key,'.dat'],'wt');
  fprintf(fid,'%d\t%d\n',1,size(elem,1));
  if(size(elem,2)==4)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t1\n', [1:length(elem);elem']);
  elseif(size(elem,2)==5)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n',[1:length(elem);elem']);
  else
    fclose(fid);
    error('wrong elem input: must be 4 or 5 columns');
  end
  fclose(fid);
end

if(nargin<6)
  facenb=faceneighbors(elem(:,1:4));
  if(nargin<4)
     face=faceneighbors(elem(:,1:4),'rowmajor');
  end
end

if(~isempty(face))
  fid=fopen(['face_',key,'.dat'],'wt');
  fprintf(fid,'%d\t%d\n',1,size(face,1));
  if(size(face,2)==3)
    fprintf(fid,'%d\t%d\t%d\t%d\t1\n',[1:length(face);face']);
  elseif(size(face,2)==4)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',[1:length(face);face']);
  else
    fclose(fid);
    error('wrong face input: must be 3 or 4 columns');
  end
  fclose(fid);
end

if(~isempty(evol))
  fid=fopen(['velem_' key '.dat'],'wt');
  fprintf(fid,'%d %d\n',1,size(elem,1));
  fprintf(fid,'%d %e\n',[(1:size(elem,1))',evol]');
  fclose(fid);
end

if(~isempty(facenb))
  fid=fopen(['facenb_' key '.dat'],'wt');
  fprintf(fid,'%d %d\n',1,size(elem,1));
  fprintf(fid,'%d\t%d\t%d\t%d\n',facenb');
  fclose(fid);
end
