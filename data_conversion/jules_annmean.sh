#!/bin/bash
#Average the JULES monthly files to create annual ones
#T. Pugh

for file in cleaf croot cwood evapotrans evspsblveg gpp_gb gpp lai landCoverFrac leaf_litC lit_c_dyn_leaf lit_c_dyn_root lit_c_dyn_wood lit_c mrro nbp npp_gb npp rh tran wood_litC wood_turnover;
do
  cdo yearmean JULESC2_CRUNCEP_${file}_Monthly_1901_2014.nc JULESC2_CRUNCEP_${file}_Annual_1901_2014.nc
done

