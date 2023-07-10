# palaeoDMC for Eocene surface temperatures
PalaeoDMC is a python package that performs PALAEOclimate Data-Model_Comparison using a Gaussian process (GP) to make a statistical estimate of continuous fields of the mean and standard deviation of a geological climate quantity derived from scattered point observations. The uncertainty of each data point is accounted for in the process and must be supplied with the input. The GP can then be compared with regular gridded general circulation model (GCM, “climate model”) output and several comparison metrics are produced and plotted.

The code in this repository is intended to document the reconstruction of Eocene surface temperature fields from the [DeepMIP proxy compilation](https://gmd.copernicus.org/articles/12/3149/2019/). The current code is identical to the method D<sub>surf-3</sub> described in (Inglis et al., 2020)[https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html] and should reproduce the associated surface temperature distributions and figures.

## Copyright
The PalaeoDMC package is still work in progress and has not been formally published. The original PalaeoDMC code was written by Fran Bragg (University of Bristol). Richard Wilkinson (University of Sheffield) wrote all code in `./core`. The PalaeoDMC code in this respository has been adapted and new NCL code has been added by Sebastian Steinig (University of Bristol). Please contact the respective person(s) if you want to use any of this code in your own work.

## Getting started
Reproducing the (Inglis et al., 2020 analysis)[https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html] can be done by following these four steps:

#### 1. Install dependencies
The analysis uses the Python and NCL languages. All necessary packages can be installed with  [conda](https://conda.io/projects/conda/en/latest/index.html). Create an environment `paleoDMC` with 

```
conda env create --name paleoDMC --file=environment.yml
``` 

using the `environment.yml` file from this repository to install all necessary packages. This environment can then be activated with

```
conda activate paleoDMC
```

#### 2. Create input proxy data
The full DeepMIP proxy compilation `data/full_deepMIP_temperature_compilation_1.5_subsampled.csv` is used to create CSV subsets for the individual time periods and subsampling experiments described in (Inglis et al., 2020)[https://cp.copernicus.org/articles/16/1953/2020/cp-16-1953-2020.html]. In `GP_analysis_wrapper`, set `convert_hollis=1` and run the script (after activating the new conda environment) with 

```
bash GP_analysis_wrapper.sh
```

The newly created CSV files will be saved to `Observation_Data/`.

![alt text](https://github.com/sebsteinig/paleoDMC/blob/main/example_output/DeepMIP_proxy_SD?raw=true)


#### 3. Run the Gaussian process regression