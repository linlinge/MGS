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
	string Mtype,input_path,output_path,json_path;
	string filename;

	/* parameters */
	for(int i=1;i<argc;i++){
		string param=argv[i];
		if("-Mtype"==param || "-M"==param){
			Mtype=argv[i+1];
		}
		else if("-i"==param){
			input_path=argv[i+1];
			vector<string> str_tmp;
			StrSplit(input_path,"/",str_tmp);
			filename=str_tmp[str_tmp.size()-1];
			cout<<"**********************************"<<endl;
			cout<<filename<<endl;
		}
		else if("-o"==param){
			output_path=argv[i+1];
		}
		else if("-json"==param || "-J"==param){
			json_path=argv[i+1];
		}
		else if("-show"==param){
			is_show_progress=atoi(argv[i+1]);
		}
		else if("--help"==param || "-H"== param){
			cout<<"MGS -i path/to/input -o path/to/output -Mtype Methods_Type -json path/of/json/file"<<endl;
			cout<<"Mtype: BSC1_1 BSC1_2 BSC1_3 MGS1 MGS2"<<endl;
			return 0;
		}
	}


	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);	
	if (pcl::io::loadPLYFile<PointType>(input_path, *cloud) == -1){
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return (-1);
	}	

	pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
	kdtree->setInputCloud(cloud);
	OLEModule ole(cloud,kdtree);
	// cout<<"dmean_:"<<ole.dmean_<<endl;

	/* Hybrid Methods Generation */
	if("BSC1_1"==Mtype){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Prox(150,5);
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult("Result/M1_1");
	}
	else if("BSC1_2"==Mtype){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_D4(500,5,"1,2");
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult("Result/M1_2",hrd01.rst_slope_);
	}
	else if("BSC1_3"==Mtype){
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,"-1");
		hrd01.DemonstrateResult("Result/M1_3",hrd01.rst_meval_);
	}
	else if("MGS1"==Mtype){
		/* Tanks and Temples */
		// HybridMethods hrd01(cloud,kdtree,&ole);
		// hrd01.store_path_=output_path;
		// hrd01.FM_Prox(150,5);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);

		/* DTU scan*/
		// HybridMethods hrd01(cloud,kdtree,&ole);
		// hrd01.store_path_=output_path;
		// hrd01.FM_Prox(800,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);


		/* DTU furu*/
		// HybridMethods hrd01(cloud,kdtree,&ole);
		// hrd01.store_path_=output_path;
		// hrd01.FM_Prox(150,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);

		/* THU */
		HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.store_path_=output_path;
		hrd01.FM_Prox(150,3);
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.FM_NID(100,0.999,"1,2");
		hrd01.DemonstrateResult(output_path);

		/* Human */
		// HybridMethods hrd01(cloud,kdtree,&ole,output_path);
		// hrd01.FM_Prox(400,6);
		// // hrd01.FM_RegionGrowth(2,80.0,4000.0,"-1");
		// hrd01.FM_NID(100,0.7,"1,2");
		// hrd01.DemonstrateResult(output_path);
	}
	else if("MGS2"==Mtype){
		/* THU */
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.FM_MEval(400,1.5);
		// hrd01.FM_MajorityVote(4000,"-1");
		// hrd01.FM_D4(500,30,"1,2");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		// hrd01.FM_NID(30,0.999,"1,2,3,4");
		// hrd01.DemonstrateResult("Result/MGS3");

		/* Tanks and Temples*/
		HybridMethods hrd01(cloud,kdtree,&ole); 
		hrd01.store_path_=output_path;
		hrd01.FM_MEval(100,9);
		hrd01.FM_MajorityVote(4000,"-1");
		hrd01.FM_D4(500,30,"1,2");
		hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		hrd01.FM_NID(30,0.7,"1,2,3,4");
		hrd01.DemonstrateResult(output_path);

		/* Human */
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.FM_MEval(200,3);
		// hrd01.FM_D4(500,30,"1");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-2");
		// hrd01.FM_NID(30,0.6,"1,2,3");
		// hrd01.DemonstrateResult(output_path);
	}
	cout<<endl<<endl;
	return 0; 
}