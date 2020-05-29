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
#include<time.h>
#include "stdlib.h"
#include "time.h"   
#include "HybridMethods.h" 
#include "Statistics.h"
#include "StringExtend.h"

int main(int argc,char** argv)
{
	cout<<argv[1]<<endl;
	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);	
	if (pcl::io::loadPLYFile<PointType>(argv[1], *cloud) == -1){
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return (-1);
	}
	

	pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
	kdtree->setInputCloud(cloud);
	OLEModule ole(cloud,kdtree);
	cout<<"dmean_:"<<ole.dmean_<<endl;

	/* Hybrid Methods Generation */
	string str="MGS3";
	if("MGS1_1"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(150,5);
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult("Result/M1_1");
	}
	else if("MGS1_2"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_DB2(500,5,"1,2");
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult("Result/M1_2",hrd01.rst_slope_);
	}
	else if("MGS1_3"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,"-1");
		hrd01.DemonstrateResult("Result/M1_3",hrd01.rst_meval_);
	}
	else if("MGS2"==str){
		double startTime = omp_get_wtime();
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(150,5);
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.FM_NID(100,0.95,"1,2");
		double stopTime = omp_get_wtime();
		double secsElapsed = stopTime - startTime;
		cout<<"Elapsed:"<<secsElapsed<<endl;
		hrd01.DemonstrateResult("Result/test");
	}
	else if("MGS3"==str){
		HybridMethods hrd01(cloud,kdtree,&ole); 
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,"-1");
		hrd01.FM_DB2(500,5,"1,2");
		hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		hrd01.DemonstrateResult("Result/M3");
	}
	return 0; 
}