#!/bin/bash

# format:
#    ./run_bench.sh  1234cpugpu <more mmcl parameters>

OPT="-F bin -n 1e6 -D T -l"
EXTRA=${@:2}
CPUEXTRA="-n 1e7"
#CPUEXTRA=-L

hostid=$(hostname -s)

echo OPT=$OPT
echo EXTRA=$EXTRA

OUTPUTDIR=../mmclbench/$hostid

mkdir -p $OUTPUTDIR

################################################################
# B1
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "1" ]]
then
    cd ../validation

    if [ $# -eq 0 ] || [[ $1 =~ "cpu" ]]
    then    
	echo testing benchmark B1 mmc mode with mmc
	./run_test.sh -s b1_mmc $OPT -G -1 -M S $EXTRA $CPUEXTRA

	echo testing benchmark B1 dmmc mode with mmc
	./run_tess.sh -s b1_dmmc $OPT -G -1 -M G -F nii $EXTRA $CPUEXTRA
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B1 mmc mode with mmcl
	./run_test.sh -s b1_mmcl $OPT -G 1 -M S $EXTRA

	echo testing benchmark B1 dmmc mode with mmcl
	./run_tess.sh -s b1_dmmcl $OPT -G 1 -M G -F nii $EXTRA
    fi

    cp b1_*mmc*.log $OUTPUTDIR
fi

################################################################
# B2
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "2" ]]
then
    cd ../sphshells

    if [ $# -eq 0 ] || [[ $1 =~ "cpu" ]]
    then   
	echo testing benchmark B2 mmc mode with mmc
	./run_mmc.sh -s b2_mmc $OPT -G -1 -M S $EXTRA $CPUEXTRA

	echo testing benchmark B2 dmmc mode with mmc
	./run_dmmc.sh -s b2_dmmc $OPT -G -1 -M G -F nii $EXTRA $CPUEXTRA
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B2 mmc mode with mmcl
	./run_mmc.sh -s b2_mmcl $OPT -G 1 -M S $EXTRA

	echo testing benchmark B2 dmmc mode with mmcl
	./run_dmmc.sh -s b2_dmmcl $OPT -G 1 -M G -F nii $EXTRA
    fi

    cp b2_*mmc*.log $OUTPUTDIR
fi


################################################################
# B3
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "3" ]]
then
    cd ../colin27

    if [ $# -eq 0 ] || [[ $1 =~ "cpu" ]]
    then   
	echo testing benchmark B3 mmc mode with mmc
	./run_test.sh -s b3_mmc $OPT -G -1 -M S $EXTRA $CPUEXTRA

	#echo testing benchmark B3 dmmc mode with mmc
	#./run_test.sh -s b3_dmmc $OPT -G -1 -M G -F nii $EXTRA $CPUEXTRA
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B3 mmc mode with mmcl
	./run_test.sh -s b3_mmcl $OPT -G 1 -M S $EXTRA

	#echo testing benchmark B3 dmmc mode with mmcl
	#./run_dmmc.sh -s b3_dmmcl $OPT -G 1 -M G -F nii $EXTRA
    fi
    cp b3_*mmc*.log $OUTPUTDIR
fi

################################################################
# B4
################################################################

if [ $# -eq 0 ] || [[ $1 =~ "4" ]]
then
    cd ../skinvessel

    if [ $# -eq 0 ] || [[ $1 =~ "cpu" ]]
    then   
	echo testing benchmark B4 mmc mode with mmc
	./run_mmc.sh -s b4_mmc $OPT -G -1 -M S $EXTRA $CPUEXTRA

	echo testing benchmark B4 dmmc mode with mmc
	./run_dmmc.sh -s b4_dmmc $OPT -G -1 -M G -F nii $EXTRA $CPUEXTRA
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B4 mmc mode with mmcl
	./run_mmc.sh -s b4_mmcl $OPT -G 1 -M S $EXTRA

	echo testing benchmark B4 dmmc mode with mmcl
	./run_dmmc.sh -s b4_dmmcl $OPT -G 1 -M G -F nii $EXTRA
    fi

    cp b4_*mmc*.log $OUTPUTDIR
fi

cd ../mmclbench
