#!/bin/bash
analysis_tag="PPCG_primaryOnly"
date_tag="2026-06-26"

## make sure the output directory is updated

overwrite=true
myfolder=/.mounts/labs/reimandlab/private/users/jreimand/PPCG_drivers_pipeline/data/${date_tag}/GISTIC2_run1
mkdir -p $myfolder

basedirPaM=${myfolder}/GISTIC2_results/
mkdir $basedirPaM

#define this file with CNA segments upstream
segfilePaM=${myfolder}/CNA_segments_for_gistic.tsv

# reference gene file from hg19
refgenefile=/PATH_TO_WORKING_DIR//data/GISTIC_INPUT_Diogo_2024-12-02/hg19.mat



###################################
#OICR's HPC specific part to deal with cluster.
MODULEPATH=:/PATH_TO_MODULES/Ubuntu18.04
# To load the libraries of the cluster and the module command
source /etc/profile.d/uge_settings.sh
# To populate your MODULEPATH
source /etc/profile.d/modulator_env.sh

module load matlab

# jr: another one to solve libcurses
ln -s /usr/lib/x86_64-linux-gnu/libncursesw.so.6 libncurses.so.5
export LD_LIBRARY_PATH=/PATH_TO_GISTIC:$LD_LIBRARY_PATH
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
export LD_LIBRARY_PATH=/PATH_TO_MATLAB:$LD_LIBRARY_PATH


echo 'loaded matlab'
#### end of cluster specific part
###################################

cd /PATH_TO_GISTIC

#this folder is problematic!
if [ -d ~/.mcrCache8.3 ] || [[ -d ~/.mcrCache9.7 ]]; then
	echo Super Cache warning! Maybe "rm -rf ~/.mcrCach*"
fi


if [ ! -e ${basedirPaM}*.scores.gistic ] || $overwrite; then
	## diogo orig call
	./gistic2 -b $basedirPaM -seg $segfilePaM -refgene $refgenefile  -rx 0 -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -fname $analysis_tag -saveseg 0 -savedata 0 
fi

 
echo finished
