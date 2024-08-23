# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco_tools/-/releases

## [2.0.0] - 2024-04-22

### Added

- crocotools_param.m: deal with different mercator datasets, see issue [#23](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/23)
    - Reanalysis Glorys/Analysis:  [ *see make_OGCM_mercator.m* ]
    - Forecast/Mediterranean high resolution Forecasts : [ *see make_OGCM_mercator_frcst.m* ]

- Forecast_tools, see issue [#23](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/23)
  - Add pause instance to download GFS datasets to avoid over rate limit violation
  - Add dependency on time step to write on template croco_forecast.in (filled within run_croco_forecast.bash)
  - Remove hindcast runs (and iterations for nudging OGCM data), keep only one forecast run on hindcast/forecast days.

- Aforc_ERA5, see issue [#23](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/23)
  - Add python request for extraction in parallel mode
  - Add boundary creation with WKB model

- Coupling_tools: updates for machines, and improvements [#25](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/25)
  -  Re-organize and improve configure.namelist.wps
  -  Add leftraru and wchpc in run_wps
  -  Move job* names according to MACHINES as in WRF_IN

## Fixed

- Nesting_tools : fix in nested_grid, see issue [#15](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/15)

- Rivers : wrong runoff positioning after first guess, fixed see issue [#18](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/18)

## Changed

- Preprocessing_tools :  Update creation of forcing files for dust and nitrogen deposition in line with the PISCES code. See issue [#19](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/19)

- Oforc_OGCM : see issue [#23](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/23)
  - Add new Copernicus Marine Toolbox to deal with mercator datasets (Please see Copernicus_Marine_Toolbox_installation.md)
  - Launch SODA preprocessing with make_OGCM_SODA.m
  - Lauch MERCATOR/CMEMS preprocessing with make_OGCM_mercator.m

## Removed
- No more use of ECCO datasets for oceanic reanalysis, see issue [#23](https://gitlab.inria.fr/croco-ocean/croco_tools/-/issues/23)
- Remove croco_pyvisu which is now hosted in croco_pytools repository

### Deprecated

### Other

