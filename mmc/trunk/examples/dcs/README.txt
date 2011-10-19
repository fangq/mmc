*Author: Stefan Carp (carp <at> nmr.mgh.harvard.edu or carp <at> alum.mit.edu)
*Date:   October 18, 2011

= Validation of MMCM simulation of diffuse correlation spectroscopy auto-correlation measurements =

== Introduction ==

Monte Carlo simulations of light propagation are generally used to estimate fluence distributions, 
or the intensity of diffusely reflected light from a turbid medium. By storing the history of 
scattering angles on a per-photon basis (momentum transfer) in addition to the photon pathlenghts
it is possible to derive the expected electric field auto-correlation on the surface of the medium
for a given distribution of scatterer dynamics in the medium. 
    
Diffuse correlation spectroscopy (DCS) as Diffusing Wave Spectroscopy (DWS) is often called in the
biomedical arena, is a non-invasive method for measuring tissue blood flow based on the extensions
of the dynamic light scattering (DLS) to the multiple scattering regime. 

DCS measurements consist of measuring diffusely reflected light using photon counting photo detectors,
and computing the temporal intensity auto-correlation curve g2(tau) for typical correlation delays from
100 ns to 1 s. These curves can be fitted using analytical expressions derived from the correlation
diffusion equation to estimate scatterer motion parameters. These expression give the field auto-
correlation function g1, related to the measured g2 by the Siegert relationship g2=1+beta*g1^2,
where beta is a factor that depends on the detection geometry and light polarization (usually lost
in tissue due to multiple scattering). 

== This example ==

Here we use an MMCM simulation of photon propagation in a thick slab (considered semi-infinite) for 
a symmetrical pattern of 4 2mm diameter detectors 15 mm away from a source. We turn on momentum 
transfer recording (the "-m" parameter in the command line), and compute the field auto-correlation
curve g1 assuming Brownian diffusion like scatterer motion (shown to approximate biological tissue
measurements well). We also simulate g2 by assuming a beta of 0.4 (typical for room light 
experimental conditions). We then fit both simulated curves using the analytical model and attempt
to recover the Brownian diffusion coefficient. The 2nd example where a simulated intensity auto-
correlation is fit represents the way experimental DCS data would be processed to obtain a blood
flow index. 

== Steps ==

1. First, you need to create the mesh to run this example. If you have not done it already,
install iso2mesh from Sourceforge:

    http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
    
2. Start Matlab, run createmesh and make a note of the value of the element ID that encloses
the source (reported by the createmesh script). Edit "dcs.inp" to set line 7 to this number.
If you receive errors, make sure iso2mesh is in your Matlab/Octave path 

3. Run the simulation shell script "run_test.sh" from the command line. It will run 20 million photons.

4. Run "dcs_example.m" in Matlab to produce the test-figures and verify the reported fitting error is
less than 5%. 

== Reference ==

Boas D.A., "Diffuse Photon Probes of Structural and Dynamical Properties of Turbid Media: Theory and 
Biomedical Applications", Ph.D. Thesis, Univ. of Pennsylvania, 1996

