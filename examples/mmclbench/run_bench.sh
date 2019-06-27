#!/bin/bash

################################################################
# B1
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "1" ]]
then
    cd ../validation

    echo testing benchmark B1 mmc mode with mmc
    ./run_test.sh -s b1_mmc -F bin -n 1e6 -G -1 -M S -D T -l $@

    echo testing benchmark B1 mmc mode with mmcl
    ./run_test.sh -s b1_mmcl -F bin -n 1e6 -G 1 -M S -D T -l $@

    echo testing benchmark B1 dmmc mode with mmc
    ./run_test.sh -s b1_dmmc -F bin -n 1e6 -G -1 -M G -D T -l $@

    echo testing benchmark B1 dmmc mode with mmcl
    ./run_test.sh -s b1_dmmcl -F bin -n 1e6 -G 1 -M G -D T -l $@

    cp b1_*mmc*.log ../mmclbench

fi

################################################################
# B2
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "2" ]]
then
    cd ../sphshells

    echo testing benchmark B2 mmc mode with mmc
    ./run_mmc.sh -s b2_mmc -F bin -n 1e6 -G -1 -M S -D T -l $@

    echo testing benchmark B2 mmc mode with mmcl
    ./run_mmc.sh -s b2_mmcl -F bin -n 1e6 -G 1 -M S -D T -l $@

    echo testing benchmark B2 dmmc mode with mmc
    ./run_dmmc.sh -s b2_dmmc -F bin -n 1e6 -G -1 -M G -D T -l $@

    echo testing benchmark B2 dmmc mode with mmcl
    ./run_dmmc.sh -s b2_dmmcl -F bin -n 1e6 -G 1 -M G -D T -l $@

    cp b2_*mmc*.log ../mmclbench
fi


################################################################
# B3
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "3" ]]
then
    cd ../colin27

    echo testing benchmark B3 mmc mode with mmc
    ./run_mmc.sh -s b3_mmc -F bin -n 1e6 -G -1 -M S -D T -l $@

    echo testing benchmark B3 mmc mode with mmcl
    ./run_mmc.sh -s b3_mmcl -F bin -n 1e6 -G 1 -M S -D T -l $@

    echo testing benchmark B3 dmmc mode with mmc
    ./run_dmmc.sh -s b3_dmmc -F bin -n 1e6 -G -1 -M G -D T -l $@

    #echo testing benchmark B3 dmmc mode with mmcl
    #./run_dmmc.sh -s b3_dmmcl -F bin -n 1e6 -G 1 -M G -D T -l $@

    cp b3_*mmc*.log ../mmclbench
fi

################################################################
# B4
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "4" ]]
then
    cd ../skinvessel

    echo testing benchmark B4 mmc mode with mmc
    ./run_mmc.sh -s b4_mmc -F bin -n 1e6 -G -1 -M S -D T -l $@

    echo testing benchmark B4 mmc mode with mmcl
    ./run_mmc.sh -s b4_mmcl -F bin -n 1e6 -G 1 -M S -D T -l $@

    echo testing benchmark B4 dmmc mode with mmc
    ./run_dmmc.sh -s b4_dmmc -F bin -n 1e6 -G -1 -M G -D T -l $@

    echo testing benchmark B4 dmmc mode with mmcl
    ./run_dmmc.sh -s b4_dmmcl -F bin -n 1e6 -G 1 -M G -D T -l $@

    cp b4_*mmc*.log ../mmclbench
fi