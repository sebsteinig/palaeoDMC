export dir="/Volumes/external_Samsung-SSD/documents/coding_github/paleoDMC"

# time_slice_list="eeco petm lp"
time_slice_list="lp"
cut_offs="2.5"
export REF_FRAME="mantle" # or "pmag"

convert_hollis=0
calculate_GP=1
plot_GP_results=0
thin=1 # thinning factor for GP analysis (1=use all data)

# create output directories
mkdir -p ${dir}/results/netcdf
mkdir -p ${dir}/results/plots

# loop over cut_offs list for sensitivity tests
for cut_off in ${cut_offs}; do
  export CUT_OFF=${cut_off}

  # create input observations from Hollis et al. (2019) CSV file
  if [ ${convert_hollis} -eq 1 ]; then
    cd ${dir}
    ncl convert_hollis2PaleoDMC.ncl
  fi

  if [ ${calculate_GP} -eq 1 ]; then
    for time_slice in ${time_slice_list}; do
      cd ${dir}
      rm -f ${dir}/Output/O-DeepMIP_${time_slice}_M-no_model_Th-${thin}/*
      python run_GP_Palaeo_DMC.py no_model DeepMIP_${time_slice} ${thin}
      cd ${dir}/Output/O-DeepMIP_${time_slice}_M-no_model_Th-${thin}
      cp *.nc ${dir}/results/netcdf/
    done
  fi

  if [ ${plot_GP_results} -eq 1 ]; then
    cd ${dir}
    ncl plot_GP_results.ncl
  fi

done
