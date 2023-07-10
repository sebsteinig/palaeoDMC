# palaeoDMC for Eocene surface temperatures
PalaeoDMC is a python package that performs PALAEOclimate Data-Model_Comparison using a Gaussian process (GP) to make a statistical estimate of continuous fields of the mean and standard deviation of a geological climate quantity derived from scattered point observations. The uncertainty of each data point is accounted for in the process and must be supplied with the input. The GP can then be compared with regular gridded general circulation model (GCM, “climate model”) output and several comparison metrics are produced and plotted.

The code in this repository is intended to document the reconstruction of Eocene surface temperature fields from the [DeepMIP proxy compilation](https://gmd.copernicus.org/articles/12/3149/2019/). The current code is identical to the method D<sub>surf-3</sub> described in (Inglis et al., 2020)[https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html] and should reproduce the associated surface temperature distributions and figures.

## Copyright
The PalaeoDMC package is still work in progress and has not been formally published. The original PalaeoDMC code was written by Fran Bragg (University of Bristol). Richard Wilkinson (University of Sheffield) wrote all code in `./core`. The PalaeoDMC code in this respository has been adapted and new NCL code has been added by Sebastian Steinig (University of Bristol). Please contact the respective person(s) if you want to use any of this code in your own work.

## Getting started
Reproducing the [Inglis et al. (2020) analysis](https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html) can be done by following these five steps:

### 1. Install dependencies
The analysis uses the Python and NCL languages. All necessary packages can be installed with  [conda](https://conda.io/projects/conda/en/latest/index.html). Create an environment `paleoDMC` with 

```
conda env create --name paleoDMC --file=environment.yml
``` 

using the `environment.yml` file from this repository to install all necessary packages. This environment can then be activated with

```
conda activate paleoDMC
```

### 2. Set working directory
In `GP_analysis_wrapper`, set the working directory `dir` in the very first line to the full path of your local repo directory.

### 3. Create input proxy data
The full DeepMIP proxy compilation `data/full_deepMIP_temperature_compilation_1.5_subsampled.csv` is used to create CSV subsets for the individual time periods and subsampling experiments described in [Inglis et al. (2020)](https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html). In `GP_analysis_wrapper`, set `convert_hollis=1` and run the script (after activating the new conda environment) with `bash GP_analysis_wrapper.sh`.

The newly created CSV files are saved to `Observation_Data/`. The main input for this step is a minimum uncertainty threshold (`cut_offs`) that gets applied to all proxy data. The reasoning is that some entries in the DeepMIP proxy collection have unrealistically low uncertainties, which would increase their relevant influence on the GP results. A values of `cut_offs="2.5"` has been used for [Inglis et al. (2020)](https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html), but you can enter a list of values to check the sensitivity of the results. Original and updated uncertainties of the full compilation (using a threshold of 2.5 degC) are shown in the Figure below (which gets automatically created by the script).

![alt text](https://github.com/sebsteinig/paleoDMC/blob/main/example_output/DeepMIP_proxy_SD.png?raw=true)


### 4. Run the Gaussian process regression
In `GP_analysis_wrapper`, set `calculate_GP=1` and run the script (again, after activating the new conda environment) with `bash GP_analysis_wrapper.sh`.

This will run the GP analysis on each subsampling experiment individually and saves the output as netCDF files in the `Output/` directory. Required computation time mainly depends on the size of the chosen output grid (`nc_m`). Results in [Inglis et al. (2020)](https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html) use `nc_m = GCM_dir+'mask_144x72.nc'` (which can take several hours to be computed for all time slices and experiments). For testing, `nc_m = GCM_dir+'mask_36x18.nc'` can be used instead.

### 5. Plot the Gaussian process regression results
In `GP_analysis_wrapper`, set `plot_GP_results=1` and run the script (again, after activating the new conda environment) with `bash GP_analysis_wrapper.sh`.

This will recreate the maps from [Inglis et al. (2020)](https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html) based on the results produced above. Flags for which plots to produce can be set in `plot_GP_results.ncl`. All plots are saved to `results/plots/`.