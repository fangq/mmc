#!/bin/bash

# format:
#    ./run_bench.sh  1234cpugpu <more mmc parameters>

OPT="-F nii -n 1e8 -D T"
EXTRA=${@:2}
CPUEXTRA="-n 1e8 -G -1"
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
	./run_test.sh -s b1_mmc $OPT -M S $EXTRA $CPUEXTRA > $OUTPUTDIR/b1_mmc.log

	echo testing benchmark B1 dmmc mode with mmc
	./run_tess.sh -s b1_dmmc $OPT -M G -F nii $EXTRA $CPUEXTRA > $OUTPUTDIR/b1_dmmc.log
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B1 mmc mode with mmc
	./run_test.sh -s b1_mmc $OPT -G 1 -M S $EXTRA > $OUTPUTDIR/b1_mmc.log

	echo testing benchmark B1 dmmc mode with mmc
	./run_tess.sh -s b1_dmmc $OPT -G 1 -M G -F nii $EXTRA > $OUTPUTDIR/b1_dmmc.log
    fi
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
	./run_mmc.sh -s b2_mmc $OPT -M S $EXTRA $CPUEXTRA > $OUTPUTDIR/b2_mmc.log

	echo testing benchmark B2 dmmc mode with mmc
	./run_dmmc.sh -s b2_dmmc $OPT -M G -F nii $EXTRA $CPUEXTRA > $OUTPUTDIR/b2_dmmc.log
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B2 mmc mode with mmc
	./run_mmc.sh -s b2_mmc $OPT -G 1 -M S $EXTRA > $OUTPUTDIR/b2_mmc.log

	echo testing benchmark B2 dmmc mode with mmc
	./run_dmmc.sh -s b2_dmmc $OPT -G 1 -M G -F nii $EXTRA > $OUTPUTDIR/b2_dmmc.log
    fi
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
	./run_test.sh -s b3_mmc $OPT -M S $EXTRA $CPUEXTRA > $OUTPUTDIR/b3_mmc.log

	echo testing benchmark B3 dmmc mode with mmc
	./run_test.sh -s b3_dmmc $OPT -M G -F nii $EXTRA $CPUEXTRA > $OUTPUTDIR/b3_dmmc.log
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B3 mmc mode with mmc
	./run_test.sh -s b3_mmc $OPT -G 1 -M S $EXTRA > $OUTPUTDIR/b3_mmc.log

	echo testing benchmark B3 dmmc mode with mmc
	./run_test.sh -s b3_dmmc $OPT -G 1 -M G -F nii $EXTRA > $OUTPUTDIR/b3_dmmc.log
    fi
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
	./run_mmc.sh -s b4_mmc $OPT -M S $EXTRA $CPUEXTRA > $OUTPUTDIR/b4_mmc.log

	echo testing benchmark B4 dmmc mode with mmc
	./run_dmmc.sh -s b4_dmmc $OPT -M G -F nii $EXTRA $CPUEXTRA > $OUTPUTDIR/b4_dmmc.log
    fi

    if [ $# -eq 0 ] || [[ $1 =~ "gpu" ]]
    then    
	echo testing benchmark B4 mmc mode with mmc
	./run_mmc.sh -s b4_mmc $OPT -G 1 -M S $EXTRA > $OUTPUTDIR/b4_mmc.log

	#echo testing benchmark B4 dmmc mode with mmc
	#./run_dmmc.sh -s b4_dmmc $OPT -G 1 -M G -F nii $EXTRA > $OUTPUTDIR/b4_dmmc.log
    fi
fi

cd ../mmclbench
