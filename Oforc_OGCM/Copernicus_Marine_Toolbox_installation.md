# Copernicus Marine Toolbox installation

The Copernicus Marine Toolbox is used to download and extract data from the Copernicus Marine Data Store (Global Ocean Physics Reanalysis, Global Ocean Physics Analysis and Forecast, ...) to create oceanic open boundary and initial conditions. [ https://data.marine.copernicus.eu/products ]


All the instructions can be found here :
   *  1) https://help.marine.copernicus.eu/en/collections/4060068-copernicus-marine-toolbox
   *  2) https://help.marine.copernicus.eu/en/articles/7970514-copernicus-marine-toolbox-installation


Prerequisites:
   * python version >= 3.9 & <3.12


---
## 2 ways to download it:

   1) typing :
      ```console
      python -m pip install copernicusmarine
      ```


   2) Use of conda/mamba package (replace conda by mamba command below)

      * Install a dedicated environment, by default named cmt_1.0 . First, you need to copy the file copernicusmarine_env.yml from the directory Forecast_tools/CopernicusMarineToolbox/

         ```console
         conda env create -f copernicusmarine_env.yml
         ```

        Note that you can use micromamba instead of mamba to install the python environment


         Firstly, you need to activate the environment cmt_1.0 :

         ```console
         conda activate cmt_1.0
         ```

         The location of the executable (here after pathCMC) will be found typing :

         ```console
         ls $CONDA_PREFIX/bin/copernicusmarine
         ```
        Note that the value returned by your terminal here, will be your pathCMC to fill in your
           crocotools_param.m and download_glorys_data.sh script


---
## When it's installed :

You will have access to the copernicusmarine executable with various sub-command, very useful to get, extract and download data from the Copernicus Marine Data Store :

```console
copernicusmarine -h
```

Or alternatively:

```console
$CONDA_PREFIX/bin/copernicusmarine -h
```


To be used in the Matlab croco_tools, in the crocotools_param.m, you need to :

- define the path to the copernicusmarine executable (`pathCMC`) in crocotools_param.m at *section 7, Option for make_OGCM_frcst or make_OGCM_mercator*
- define your copernicus marine `login` and `password` in crocotools_param.m at *section 7, Options for for make_OGCM_frcst or make_OGCM_mercator*





