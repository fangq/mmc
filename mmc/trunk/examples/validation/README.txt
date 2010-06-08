1. First, please download and install the svn version of iso2mesh from 

 svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/iso2mesh/trunk/iso2mesh iso2mesh

the password is anonymous_user

2. Start matlab, run createmesh to create all mesh files

3. run the simulation bash script run_test.sh

4. To plot results from MCX, you need to download it from 

 svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/mcextreme/mcextreme_cuda/trunk/ mcx

the password is anonymous_user

5. compile and install MCX, and run MCX simulation

5.1 go to examples/mcxsph run createmcxbin from matlab
5.2 open benchbox.sh, edit the thread number (-t) and thread block size (-T)
based on the compute capability of your card. For 8800GT/9800GT, the -T number
can not exceed 128; for 280/285/295, -T can not be more than 256; for 470, 
-T can not be more than 576
5.3 run benchbox.sh to generate the MCX output box.mc2

6. when the simulation is complete, you can start matlab again, and 
run plotcuberes to generate the plots.
