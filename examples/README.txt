= README for examples =

In this folder, you will find a number of examples 
to validate, test, or demonstrate the functionalities 
of MMC.

Please read the README.txt file under each sub-folder
to understand the purpose, procedures to run the simulation 
and the interpretations of the results.

The examples contained in each subfolders are explained below:

examples/
|-validation - homogeneous cubic domain 60x60x60, i.e. Fig. 2 in Fang2010
|-meshtest   - 10mm sphere embedded in the cubic domain, Fig. 3 in Fang2010
|-mcxsph     - scripts to run MCX for the above two cases
|-onecube    - debugging flag show case 1, simple mesh without reflection
|-reftest    - debugging flag show case 2, simple mesh with reflection
|-rngtest    - random number generator test
|-dcs        - validation of momentum transfer in DCS simulation
|-statnoise  - test statistical noise against launched photon numbers
|-misctest   - other misc. tests
|-planar     - validation of uniform planar source in a simple cube
|-sphere     - validation of uniform planar source use a spherical domain
|-sfdi2layer - validation SFDI illumination on a two-layer brain model
|-regression/exitangle - regression testing of photon exit angles

A complex human brain atlas mesh, i.e. Fig. 4 in [Fang2010], can 
be downloaded separately at the following URL:

 http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh


Reference: 

[Fang2010] Qianqian Fang, "Mesh-based Monte Carlo method using fast ray-tracing \
           in Plücker coordinates," Biomed. Opt. Express 1, 165-175 (2010) 
[Yao2016]  Ruoyang Yao, Xavier Intes, Qianqian Fang, "Generalized mesh-based \
           Monte Carlo for wide-field illumination and detection via mesh \
           retessellation," Biomed. Optics Express, 7(1), 171-184 (2016)
