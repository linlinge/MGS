# MMOD: A Flexible Module-based Outliers Detection Method for 3D Point Clouds
This is the source code for the paper "MMOD: A Flexible Module-based Outliers Detection Method for 3D Point Clouds". The file structure description for this project is presented in the following.  
```
|---- main.cpp : It contains the entry of the main function, and you can choose different patterns in this file.  
|---- glib : This folder contains some basic functions write by myself.
|---- OLEModule.h, OLEModule.cpp :
|---- HybridMethods.h,.cpp : It realized the function of modules' definition and hybrid methods generation.  
```

## Prerequisities  
* PCL
* Eigen
* OpenMP

## Setup  
## Run
## Experiments
You can find the comparison methods PMR [1], LoOP [2] and PointCleanNet [3] in the reference section.
## Reference  
[1] https://pdal.io/stages/filters.reciprocity.html  
[2] MeshLab -> Filters -> Selection -> Select Outliers  
[3] https://github.com/mrakotosaon/pointcleannet  