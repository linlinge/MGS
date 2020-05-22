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

	/* Hybrid Methods Generation */
	string str="M1_2";
	if("M1_1"==str){
		/* Prox -> NID */
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(150,15);
		hrd01.FM_NID(100,0.8,REGULAR_DOMAIN);
		hrd01.D(IRREGULAR_DOMAIN);
		hrd01.DemonstrateResult("Result/M1_1");
	}
	else if("M1_2"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(100,5);
		hrd01.FM_MEval(80,2,IRREGULAR_DOMAIN);
		hrd01.D(REGULAR_DOMAIN);
		hrd01.DemonstrateResult("Result/M1_2");
	}
	else if("M1_3"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_NID(100,0.8);
		
		hrd01.DemonstrateResult("Result/M1_2");
	}
	else if("M2_1"==str){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(100,10);
		hrd01.FM_RegionGrowth(5.0,80.0,150.0,IRREGULAR_DOMAIN);
		hrd01.D(IRREGULAR_DOMAIN);
		hrd01.FM_NID(100,0.95,REGULAR_DOMAIN);
		hrd01.D(IRREGULAR_DOMAIN);
		hrd01.DemonstrateResult("Result/M2_1");
	}
	else if("M3_1"==str){
		HybridMethods hrd01(cloud,kdtree,&ole); 
		hrd01.FM_Prox(300,5);
		hrd01.FM_MEval(30,1.5,REGULAR_DOMAIN);
		hrd01.J(IRREGULAR_DOMAIN);
		hrd01.FM_MajorityVote(500,IRREGULAR_DOMAIN);
		hrd01.D(REGULAR_DOMAIN);
		hrd01.FM_Prox(300,5,REGULAR_DOMAIN);
		hrd01.J(IRREGULAR_DOMAIN);
		hrd01.DemonstrateResult("Result/M3_1");
	}
	return 0;
}