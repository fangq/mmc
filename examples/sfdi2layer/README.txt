= Validation of SFDI source =

== Introduction ==

In this example, we test a spatial frequency domain (SFDI) 
source -- on a two-layer cubic model.

Three corners of the 3D quadrilateral are specified by 
srcpos, srcpos+srcparam1(0:2), srcpos+srcparam2(0:2), 
k-numbers of the spatial frequency along two directions 
are specified by srcparam1(3) and srcparam2(3).

The original mesh is saved as "**_media.dat' and the 
retessellated mesh is saved as "**_sfdi.dat".

The optical properties mimic that of the skull and CSF 
of a human brain model, and can be found at 

http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh#Tissue_optical_properties.

Run "createmesh.m" first to generate mesh (before/after 
retessellation), then run "run_test.sh" to perform the 
simulation. Finally, run "plot_result.m" to see the 
continuous wave (CW) fluence distribution inside the model.
