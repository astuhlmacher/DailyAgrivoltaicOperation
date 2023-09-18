# DailyAgrivoltaicOperation
Source code and data files for optimizing the daily operation of dual-axis solar photovoltaic (PV) panels within an agrivoltaic system, created by Anna Stuhlmacher. The code is demonstrated with a case study located at the University of Michigan's Campus Farm in Ann Arbor, Michigan (42.3N, 83.7W) on July 14, 2021. 




## Description 
This repository contains:

- A CSV file with solar irradiance data (GHI, DHI, and DNI) for a 24-hour time horizon with a time resolution of 30 minutes. The data was obtained from the [National Solar Radiation Database
(NSRDB)](https://nsrdb.nrel.gov/). The file is named ```AnnArbor_2021_simple_july14_EDT.csv```.
- A MAT file of the solar position (azimuth and altitude) and shading factor approximations with a time resolution of 10 minutes. The solar position was calcuated using the Python package [pvlib](https://pvlib-python.readthedocs.io/en/stable/). The file is named ```AdjustTilt_90range_13density_fullday.mat```.
- Julia code used to run optimization problem, i.e., ```script.jl``` and ```apv_daily_problem.jl```.
- Supplemental Matlab and Python files used to generate solar position and shading factor approximations

# Publication
If you would like to use this data, please cite our publication:
```
@article{stuhlmacher2024agrivoltaics,
  title={Optimizing Dual-Axis Solar Panel Operation in an Agrivoltaic System and Implications for Power Systems},
  author={Anna Stuhlmacher and Johanna L. Mathieu and Peter Seiler},
  journal={Proceedings of the 57th Hawaii International Conference on System Sciences (HICSS)},
  year={2024}
} 
```

## Abstract
The concept of agrivoltaics, or co-locating photovoltaic panels and crops, is viewed as a potential solution to competing land demands for food and energy production. In this paper, we propose an optimal dual-axis photovoltaic panel formulation that adjusts the panel position to maximize power generation subject to crop requirements. Through convex relaxations and shading factor approximations, we reformulate the problem as a convex second-order cone program and solve for the panel position adjustments away from the sun-tracking trajectory. We demonstrate our approach in a case study by comparing our approach with an approach that maximizes solar power capture and a scenario in which there are only crops. We found that we are able to successfully adjust the panel position while accounting for the trade-offs between the photovoltaic panels' energy production and the crop health. Additionally, optimizing the operation of an agrivoltaic system allows us to better understand agrivoltaic systems as a resource connected to the power grid.
