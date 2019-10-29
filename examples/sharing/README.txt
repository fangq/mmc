= Simulate multiple wide-field illumination patterns through "photon sharing" =

== Introduction ==

In this example, we again use a two-layer cubic model.
Mesh retessellation is performed based on the origin model 
to accomodate for the wide-feild configuration

Five illumination patterns of 40 by 40 mm2 are simulated

The optical properties mimic that of the skull and CSF 
of a human brain model, and can be found at 

http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh#Tissue_optical_properties.

Run "createmesh.m" first to generate mesh (before/after retessellation), 
then run "run_test.sh" to perform the simulation. Finally, run "plot_result.m" 
to see the continuous wave (CW) fluence distribution of five source patterns.
