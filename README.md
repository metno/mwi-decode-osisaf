# mwi-decode-osisaf
Decoder from OSI SAF for the Metop-SG MWI L1 files distributed by EUMETSAT

## How to build the required conda environment to run the code
```
conda env create -f reqs_mwi-decode.yml
conda activate mwi-decode
```

## How to test the decoder on the included granuale test file
```
python process_mwi.py -i ./W_XX-EUMETSAT-Darmstadt,SAT,SGB1-MWI-1B-RAD_C_EUMT_20221213180200_G_D_20080223093905_20080223094135_T_N____.nc -o ./test_mwi_small.nc
```

## How to run the decoder on all the test data, in orbit form
Download test data from https://www.eumetsat.int/eps-sg-user-test-data and set an environment variable MWIDATA
```
python process_mwi.py -i ${MWIDATA}/MWI_L1B_Test_Scenario_011/W_XX-EUMETSAT-Darmstadt,SAT,SGB1-MWI-1B-RAD_C_EUMT_20220615170800_G_D_20070912084321_20070912102224_T_N____.nc+${MWIDATA}/MWI_L1B_Test_Scenario_012/W_XX-EUMETSAT-Darmstadt,SAT,SGB1-MWI-1B-RAD_C_EUMT_20220615151100_G_D_20070912102225_20070912120114_T_N____.nc -o ./test_mwi_out.nc
```
