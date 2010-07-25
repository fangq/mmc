= Matlab Toolbox for the Analytical Diffusion/Helmholtz Solutions of a Sphere =

* Version 1.0 (Release date: 07/25/2010)

The "Sphere-Diffusion" toolbox solves for the analytical solutions for the diffusion 
inside and outside a sphere. This toolbox is part of an open-source software, MMC, 
developed by Qianqian Fang. The license of the toolbox is GNU General Public License
version 3. Please see LICENSE.txt for details.

The solutions are defined in the 3D space on a user specified grid; the sphere can 
have different absorption/scattering/refractive index to the background media. This 
toolbox can be useful when evaluating new algorithms in the heterogeneous media.
The toolbox is compatible with GNU Octave.

    1. Download
    2. Core functions
    3. Contour plots of the sample solutions

            3.1.1. A 10mm sphere inside an infinite homogeneous medium
            3.1.2. A 10mm sphere inside a semi-infinite homogeneous medium
            3.1.3. A 10mm sphere inside an infinite homogeneous slab

    4. Other functions
    5. Examples
    6. Acknowledgement

1. Download

Please browse the online page of this toolbox to proceed with downloading.

 http://mcx.sourceforge.net/cgi-bin/index.cgi?SphereDiffusion

The author of the toolbox is appreciated if you can cite the reference
[Fang2010] listed at the end of this document if you choose to use
this toolbox in your publication.


2. Core functions

[res,xi,yi,zi] = sphdiffusioninfinite(xrange,yrange,zrange,cfg)
    diffusion solution for a sphere inside the infinite homogeneous medium 
[res,xi,yi,zi] = sphdiffusionsemi(Reff,xrange,yrange,zrange,cfg)
    diffusion solution for a sphere inside a semi-infinite homogeneous medium 
[res,xi,yi,zi] = sphdiffusionslab(Reff1,Reff2,h,xrange,yrange,zrange,cfg)
    diffusion solution for a sphere inside an infinite homogeneous slab 

3. Contour plots of the sample solutions

3.1.1. A 10mm sphere inside an infinite homogeneous medium

    upload:infinite_sphere.png 

3.1.2. A 10mm sphere inside a semi-infinite homogeneous medium

    upload:semiinfinite_sphere.png 

3.1.3. A 10mm sphere inside an infinite homogeneous slab

    upload:slab_sphere.png 

4. Other functions
function 	description
besselhprime.m	Hankel function first order derivative
besseljprime.m	Bessel function (Bessel first kind) first order derivative
besselyprime.m	Neumann function (Bessel second kind) first order derivative
spbesselh.m	Spherical Hankel function
spbesselhprime.m	Spherical Hankel function first order derivative
spbesselj.m	Spherical Bessel function
spbesseljprime.m	Spherical Bessel function first order derivative
spbessely.m	Spherical Neumann function
spbesselyprime.m	Spherical Neumann function first order derivative
spharmonic.m	Spherical harmonics
sphdiffAcoeff.m	Sphere exterior field coefficients
sphdiffBcoeff.m	Sphere exterior field coefficients
sphdiffCcoeff.m	Sphere interior field coefficient
sphdiffexterior.m	Sphere exterior total field
sphdiffincident.m	Incident field
sphdiffinterior.m	Sphere interior field
sphdiffscatter.m	Sphere exterior scattered field

5. Examples

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % setting up problem domain
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cfg.v=299792458000;  % speed-of-light in mm/s
 cfg.a=10;            % radius of the sphere, mm
 cfg.omua=0.002;      % outside (background) mua 1/mm
 cfg.omusp=0.990;     % outside (background) mus' 1/mm
 cfg.imua=0.050;      % inside (sphere) mua 1/mm
 cfg.imusp=0.500;     % inside (sphere) mus' 1/mm
 %cfg.imua=0.002;
 %cfg.imusp=0.990;
 cfg.src=[30,pi,0];   % source position in spherical coordinates
 cfg.maxl=20;         % maximum orders for the series expansion
 cfg.omega=0;         % modulation frequency

 cfg.Din=cfg.v/(3*cfg.imusp);  % Diffusion coefficient in the sphere
 cfg.Dout=cfg.v/(3*cfg.omusp); % Diffusion coefficient outside the sphere
 cfg.kin=sqrt((-cfg.v*cfg.imua+i*cfg.omega)/cfg.Din);   % complex-wavenumber in the sphere
 cfg.kout=sqrt((-cfg.v*cfg.omua+i*cfg.omega)/cfg.Dout); % complex-wavenumber outside the sphere

 % solution of the sphere in an infinite medium
 [phi_ana,xa,ya,za]=sphdiffusioninfinite(-30:0.8:30,0,-30:0.8:30,cfg);
 figure;contourf(xa,za,log10(abs(phi_ana)),40);axis equal;

 % solution of the sphere in an infinite slab with a height of 60mm
 [phi_ana,xa,ya,za]=sphdiffusionslab(0,0,60,-30:0.8:30,0,-30:0.8:30,cfg);
 figure;contourf(xa,za,log10(abs(phi_ana)),40);axis equal;

6. Acknowledgement
Qianqian Fang would like to thank David Boas for the helpful 
discussions on the semi-infinite solutions. 

7. Reference

* Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing \
in Plücker coordinates," Biomed. Opt. Express 1, 165-175 (2010) 
