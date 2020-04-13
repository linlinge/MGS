/*
	Author: linlinge
	Date: 2020.04.13
*/
#include <iostream>	
#include "PCLExtend.h"
#include "V3.hpp"
#include <Eigen/Dense>
#include <omp.h>
#include "Table.h"
#include <algorithm>
#include "VectorExtend.h"
#include "SignalProcessing.h"
#include "OLEModule.h"
#include "FEModule.h"
#include<time.h>
#include "stdlib.h"
#include "time.h"    
int main(int argc,char** argv)
{
	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);	
	if (pcl::io::loadPLYFile<PointType>(argv[1], *cloud) == -1){
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return (-1);
	}
	pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
	kdtree->setInputCloud(cloud);

	OLEModule ole(cloud,kdtree);
	ole.GetMinorEval();
	cout<<"og_: "<<ole.og_<<endl;

	FEModule fem(cloud,kdtree);
	// fem.ApplyHierarchical(ole,100,0.95,0.002);
	// fem.ApplySlope(ole,200,0.8);
	fem.ApplyMEval(ole);
	fem.DemonstrateResult("Result/rst_color.ply");
	return 0;
}
