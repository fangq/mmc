function savemmcmesh(key, node, elem, varargin)
%
% savemmcmesh(key,node,elem,facenb,roi)
% savemmcmesh(key,node,elem,'facenb',facenb,'roi',roi)
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
%      facenb: the 4 neighboring elements for each element, if missing, it will
%             be calculated
%      roi: if 1 colume, roi defines noderoi; if 4 column, roi defines
%             faceroi, if 6 columns, roi defines edgeroi
%
% example:
%     savemmcmesh('sph1',node,elem);
%
% this file is part of Mesh-based Monte Carlo (MMC)
%
% License: GPLv3, see http://mcx.sf.net/mmc/ for details
%

if (~isempty(varargin) && ischar(varargin{1}))
    opt = varargin2struct(varargin{:});
else
    opt = struct;
end

if (~isempty(varargin) && ~ischar(varargin{1}))
    facenb = varargin{1};
    if (length(varargin) > 1)
        roi = varargin{2};
    end
end

if (~isempty(node))
    fid = fopen(['node_', key, '.dat'], 'wt');
    fprintf(fid, '%d\t%d\n', 1, size(node, 1));
    fprintf(fid, '%d\t%16.8e\t%16.8e\t%16.8e\n', [1:length(node); node']);
    fclose(fid);
end

if (~isempty(elem))
    elem(:, 1:4) = meshreorient(node, elem(:, 1:4));

    fid = fopen(['elem_', key, '.dat'], 'wt');
    fprintf(fid, '%d\t%d\n', 1, size(elem, 1));
    if (size(elem, 2) == 4)
        fprintf(fid, '%d\t%d\t%d\t%d\t%d\t1\n', [1:length(elem); elem']);
    elseif (size(elem, 2) == 5)
        fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\n', [1:length(elem); elem']);
    else
        fclose(fid);
        error('wrong elem input: must be 4 or 5 columns');
    end
    fclose(fid);
end

if (~isempty(jsonopt('facenb', [], opt)))
    facenb = jsonopt('facenb', [], opt);
end

if (nargin < 4 || isempty(jsonopt('facenb', [], opt)))
    facenb = faceneighbors(elem(:, 1:4));
end

if (~isempty(facenb))
    fid = fopen(['facenb_' key '.dat'], 'wt');
    fprintf(fid, '%d %d\n', 1, size(elem, 1));
    fprintf(fid, '%d\t%d\t%d\t%d\n', facenb');
    fclose(fid);
end

if (~exist('roi', 'var'))
    roi = jsonopt('roi', [], opt);
end

if (~isempty(roi))
    fid = fopen(['roi_' key '.dat'], 'wt');
    fprintf(fid, '%d %d\n', size(roi, 2), size(roi, 1));
    format = [repmat('%d\t', 1, size(roi, 2)) '\n'];
    fprintf(fid, format, roi');
    fclose(fid);
end
