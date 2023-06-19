#!/bin/bash -f
#$ -N run-mwi-stitching
#$ -l h_rss=10G,mem_free=10G,h_data=10G
#$ -cwd -b n
#$ -S /bin/bash
#$ -pe shmem-1 1
#$ -m a
#$ -j y

export LD_LIBRARY_PATH=/home/emilyjd/osisaf_mwi/osi_base/OSI_HL_AUX/libs/PROJ.4/lib
export NC_BLKSZ=80M
source /modules/rhel8/conda/install/etc/profile.d/conda.sh
conda activate production-12-2021

python /home/emilyjd/osisaf_mwi/mwi-decode-osisaf/mwi_check_and_stitch.py -i /lustre/storeB/project/fou/fd/project/osisaf/eps-sg-test-data/osisaf_data_emilyjd/mwi-EUM-decoded -o /lustre/storeB/project/fou/fd/project/osisaf/eps-sg-test-data/osisaf_data_emilyjd/mwi-EUM-stitched -O
