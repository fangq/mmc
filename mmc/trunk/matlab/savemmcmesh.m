function savemmcmesh(key,node,elem,face)

if(nargin<4)
  face=volface(elem); % volface is part of iso2mesh toolbox,http://iso2mesh.sf.net
end

if(~isempty(node))
  fid=fopen(['node_',key,'.dat'],'wt');
  fprintf(fid,'%d\t%d\n',1,size(node,1));
  fprintf(fid,'%d\t%16.8e\t%16.8e\t%16.8e\n',[1:length(node);node']);
  fclose(fid);
end

if(~isempty(elem))
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

if(~isempty(face))
  fid=fopen(['face_',key,'.dat'],'wt');
  fprintf(fid,'1 %d\n',size(face,1));
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

