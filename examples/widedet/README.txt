= Output time-resolved images when using wide-field detector =

== Notice ==

1. In input file, the number of detectors is 1 (one wide-field) and the
    radius of that detector is 0, so that the detparams can be read.

2. The last two lines of input file are detparam1 and detparam1,
    they have the same meaning as srcparam1&2 when srctype is pattern.

3. In the shell file, set -x to be 2 and to get an "img" file as output,
    the file is binary and the size equlas to nx*ny*nTG, where nx and ny
    are number of pixels of detector defined in detparam1 and detparam2
    and nTG is the total number of time gates.

4. The order of row/column is different between order of C and Matlab, 
    so we need to do a matrix transpose before plotting the image.