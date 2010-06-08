1. First, please download and install the svn version of iso2mesh from 

 svn checkout --username anonymous_user https://orbit.nmr.mgh.harvard.edu/svn/iso2mesh/trunk/iso2mesh iso2mesh

the password is anonymous_user

2. Start matlab, run createmesh to create all mesh files

3. Open initial_elem.txt, and type the integer in each line to replace
the second integer at the 7-th line of sph{1,2,3}.inp 
(this sets the initial element ID, because the mesh generator may generate
a different mesh on your system, you have to update this manually)

4. run the simulation bash script run_test.sh

5. when all 3 simulations are complete, you can start matlab again, and 
run plotmmcsph to generate the plots.
