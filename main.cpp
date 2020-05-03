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
#include "HybridModule.h" 
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
	// ole.GetMinorEval();
	FEModule(cloud,kdtree,&ole);

	// ApplySlope(100,15,default_layer);
	// ApplyRegionGrowth(1,0.001,0.3);
	// ApplyMEval(40,10.0,0);
	// D();
	// ApplySlope(300,6,0);
	// H();
	// ApplyRegionGrowth(1,0.001,0.3);
	// D();

	ApplyMEval(10,10.0,default_layer);
	ApplyMajorityVote(1,8000);
	D();
	ApplySlope(400,3,0);
	H();
	ApplyRegionGrowth(1,0.001,0.3);
	D();

	// H(ApplyMEval,40,10.0,default_layer,ApplySlope,200,10,0);
	// H(ApplyMEval,10,3.0,default_layer,ApplySlope,800,3.0,0);
	// H(ApplyMEval,10,3.0,default_layer,ApplySlope,800,3.0,0);
	DemonstrateResult("Result/rst_color.ply");
	return 0;
}