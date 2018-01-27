function [tau,g1]=generate_g1(fhist,tau, disp_model, DV, lambda, format, varargin)
%
%   [tau,g1]=generate_g1(fhist,tau, disp_model, DV, lambda, format)
%
%   Compute simulated electric-field auto-correlation function using
%   simulated photon pathlengths and scattering momentum transfer
%
%   author: Stefan Carp (carp <at> nmr.mgh.harvard.edu)
%
%   input:
%       fhist:      the file name of the output .mch file 
%       tau:        correlation times at which to compute g1 
%                   (default: 1e-7 to 1e-1 seconds, log equidistant)
%       disp_model: displacement model ('brownian', 'random_flow', <custom>)
%                   (default: brownian, see further explanation below)
%       disp_var:   value of displacement variable using mm as unit of
%                   length and s as unit of time
%                   (default: 1e-7 mm^2/s, see further explanation below)
%       lambda:     wavelenght of light used in nm
%                   (default: 785)
%       format:     the format used to save the .mch file 
%                   (default: 'float')
%
%   output:
%
%       tau:        correlation times at which g1 was computed provided for
%                   convenience (copied from input if set, otherwise 
%                   outputs default)
%       g1:         field auto-correlation curves, one for each detector
%
%   The displacement model indicates the formula used to compute the root
%   mean square displacement of scattering particles during a given delay
%   
%   brownian:       RMS= 6 * DV * tau; 
%                   DV(displacement variable)=Db (brownian diffusion coeff)
%   random_flow:    RMS= DV^2 * tau^2; 
%                   DV = V (first moment of velocity distribution)
%   <custom>:       any string other than 'brownian' or 'random_flow' will
%                   be evaluate as is using Matlab evalf, make sure it uses
%                   'DV' as the flow related independent variable, tau is
%                   indexed as tau(J). Any additional parameters can be 
%                   sent via "varargin"
%
%   This file is part of Mesh-Based Monte Carlo
%   License: GPLv3, see http://mcx.sf.net for details
%

if nargin<6, format='float'; end
if nargin<5, lambda=785; end
if nargin<4, DV=1e-7; end
if nargin<3, disp_model='brownian'; end
if nargin<2, tau=logspace(-7,-1,200); end

[mch_data,mch_header]=loadmch(fhist,format);
temp=strfind(fhist,filesep);
if isempty(temp), fhist=[pwd filesep fhist]; end
temp=strfind(fhist,filesep);
lastslash=temp(end);
sim_label=fhist(lastslash+1:end-4);

[mua,mus,g,n]=load_mc_prop([fhist(1:lastslash) filesep 'prop_' sim_label '.dat']);

if (mch_header.recordnum-2)~=(2*mch_header.medianum),
    fprintf('History file does not contain momentum transfer information \n');
    g1=-1;
    return;
end

if strcmp(disp_model,'brownian'),
    disp_str='rmsdisp=6*DV.*tau(J);';
elseif strcmp(disp_model,'random_flow'),
    disp_str='rmsdisp=DV.^2.*tau(J).^2;';
else 
    disp_str=['rmsdisp=' disp_model ';'];
end

k0=2*pi*n/(lambda*1e-6);

g1=zeros(mch_header.detnum,length(tau));

for I=1:mch_header.detnum,
    idx= find(mch_data(:,1)==I);     
    fprintf('Processing detector %.0f: %.0f photons\n',I,length(idx));

    for J=1:length(tau),
        eval(disp_str);
        g1(I,J)=sum(exp(-(k0.^2.*rmsdisp/3)*mch_data(idx,(3+mch_header.medianum):end)'-mua*mch_data(idx,3:(3+mch_header.medianum-1))'));
    end
    g1_norm=sum(exp(-mua*mch_data(idx,3:(3+mch_header.medianum-1))'));
    g1(I,:)=g1(I,:)./g1_norm;
end
    


