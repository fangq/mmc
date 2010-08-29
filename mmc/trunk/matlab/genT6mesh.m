function [vertices,tess]=genT6mesh(varargin)
% tess_lat: simplicial tessellation of a rectangular lattice
% usage: [tessellation,vertices]=genT6mesh(v1,v2,v3,...)
%
% URL: http://www.mathkb.com/Uwe/Forum.aspx/matlab/50484/Constant-surfaces
%
% arguments: input
% v1,v2,v3,... - numeric vectors defining the lattice in
% each dimension.
% Each vector must be of length >= 1
%
% arguments: (output)
% vertices - factorial lattice created from (v1,v2,v3,...)
% Each row of this array is one node in the lattice
% tess - integer array defining simplexes as references to
% rows of "vertices".

% dimension of the lattice
n = length(varargin);

% create a single n-d hypercube
% list of vertices of the cube itself
vhc=('1'==dec2bin(0:(2^n-1)));
% permutations of the integers 1:n
p=perms(1:n);
nt=factorial(n);
thc=zeros(nt,n+1);
for i=1:nt
  thc(i,:)=find(all(diff(vhc(:,p(i,:)),[],2)>=0,2))';
end

% build the complete lattice
nodecount = cellfun('length',varargin);
if any(nodecount<2)
  error 'Each dimension must be of size 2 or more.'
end
vertices = lattice(varargin{:});

% unrolled index into each hyper-rectangle in the lattice
ind = cell(1,n);
for i=1:n
ind{i} = 0:(nodecount(i)-2);
end
ind = lattice(ind{:});
k = cumprod([1,nodecount(1:(end-1))]);
ind = 1+ind*k';
nind = length(ind);

offset=vhc*k';
tess=zeros(nt*nind,n+1);
L=(1:nind)';
for i=1:nt
  tess(L,:)=repmat(ind,1,n+1)+repmat(offset(thc(i,:))',nind,1);
L=L+nind;
end

% ======== subfunction ========
function g = lattice(varargin)
% generate a factorial lattice in n variables
n=nargin;
sizes = cellfun('length',varargin);
c=cell(1,n);
[c{1:n}]=ndgrid(varargin{:});
g=zeros(prod(sizes),n);
for i=1:n
g(:,i)=c{i}(:);
end
