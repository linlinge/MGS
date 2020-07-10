# SGODF: Style-guided outlier detection framework
This is the source code for the paper "MGS: A Flexible Methods Generation Scheme for Point Cloud Outliers Detection". The file structure description for this project is presented in the following.  
```
|---- main.cpp : It contains the entry of the main function, and you can choose different patterns in this file.  
|---- glib : This folder contains some basic functions write by myself.
|---- OLEModule.h, OLEModule.cpp :
|---- HybridMethods.h,.cpp : It realized the function of modules' definition and hybrid methods generation.  
```

# How to use
## Prerequisities   
* PCL
* Eigen
* OpenMP
```
// step one
// step 1-1: download opencv from https://opencv.org/releases/
// step 1-2: install opencv
# mkdir build
# cd build
# cmake ..
# make -j8
# make install

// step two : install pcl
# sudo add-apt-repository ppa:v-launchpad-jochen-sprickerhof-de/pcl
# sudo apt-get update
# sudo apt-get install libpcl-dev

// step three: install eigen
# sudo apt-get install libeigen3-dev
# sudo cp -r /usr/local/include/eigen3 /usr/include 
```

## Compile  
```
# mkdir build
# cmake ..
# make -j12
```
# Experiments
You can find the comparison methods PMR [1], LoOP [2] and PointCleanNet [3] in the reference section.
# Reference  
[1] https://pdal.io/stages/filters.reciprocity.html  
[2] MeshLab -> Filters -> Selection -> Select Outliers  
[3] https://github.com/mrakotosaon/pointcleannet  