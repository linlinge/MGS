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
#include "StyleEvaluation.h"
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
	string fileNoType;

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
			str_tmp.clear();
			StrSplit(input_path,".",str_tmp);
			fileNoType=str_tmp[0];
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
			cout<<"Mtype: OGEM UEM SEM BSC1_1 BSC1_2 BSC1_3 MGS1 MGS2"<<endl;
			return 0;
		}
	}

	/* Load Model */
	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);	
	if (pcl::io::loadPLYFile<PointType>(input_path, *cloud) == -1){
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return (-1);
	}	
	pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
	kdtree->setInputCloud(cloud);


	/* Hybrid Methods Generation */	
	StyleEvaluation ole(cloud,kdtree);
	HybridMethods hrd01(cloud,kdtree,&ole);	
	if("OutlierGrade"==Mtype){		
		ole.OutlierGrade(output_path,50);
	}
	else if("meval"==Mtype){
		hrd01.FM_MEval(50,3.0);
		hrd01.DemonstrateResult_Color("Result/color.ply",hrd01.rst_meval_,"nonlinear");
	}
	else if("Homogeneity"==Mtype){		
		ole.Homogeneity(output_path);
	}
	else if("Singularity"==Mtype){
		ole.Sigularity(output_path,60,2000);
	}
	else if("OR-BG"==Mtype){
		/* Homework ply */
		hrd01.FM_Slope(150,5);
		// hrd01.rst_db2_.Normalize_Min_Max();
		// hrd01.rst_db2_.Normalize_Tanh();
		hrd01.DemonstrateResult_Color("Result/color.ply",hrd01.rst_slope_);
		// hrd01.DemonstrateResult_Color("Result/color.ply",hrd01.rst_slope_,"nonlinear");
	}
	else if("OR-D4"==Mtype){
		/* Homework ply */
		hrd01.FM_D4(3,1);
		// hrd01.FM_D4(30,1);
		// hrd01.FM_D4(30,10);
		// hrd01.FM_D4(70,3); /* S */		
		// hrd01.FM_D4(200,3);
		// hrd01.FM_D4(40,3);
		hrd01.DemonstrateResult_Color("Result/color.ply",hrd01.rst_db2_);
	}
	else if("BSC1_1"==Mtype){
		// HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_Slope(150,5);
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult("Result/M1_1");
	}
	else if("BSC1_2"==Mtype){
		// HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_D4(500,5,"1,2");
		hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		hrd01.DemonstrateResult_Color("Result/M1_2",hrd01.rst_slope_);
	}
	else if("BSC1_3"==Mtype){
		// HybridMethods hrd01(cloud,kdtree,&ole);
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,0.01,"-1");
		hrd01.DemonstrateResult_Color("Result/M1_3",hrd01.rst_meval_);
	}
	else if("MGS1"==Mtype){
		hrd01.FM_D4(100,4);
		hrd01.FM_LoOP();
		hrd01.DemonstrateResult(output_path);
	}
	else if("MGS2"==Mtype){
		// hrd01.FM_Slope(200,3);
		// hrd01.FM_NID(100,2);
		
		hrd01.FM_MEval(80,30);
		hrd01.FM_LoOP("1");
		// hrd01.FM_RegionGrowth(0.1,20.0,3.0,"-1");
		hrd01.DemonstrateResult(output_path);
	}
	// else if("MGS2"==Mtype){
	// 	// hrd01.FM_Slope(60,10);
	// 	hrd01.FM_D4(70,10);
	// 	hrd01.rst_db2_.Normalize_Min_Max();
	// 	hrd01.FM_Density(80,1);
	// 	hrd01.rst_density_.Normalize_Min_Max();
	// 	Table<Rrd1> rst_tmp;
	// 	rst_tmp.Resize(cloud->points.size());
	// 	for(int i=0;i<cloud->points.size();i++){
	// 		rst_tmp.records_[i].item1_=hrd01.rst_db2_.records_[i].item1_*hrd01.rst_density_.records_[i].item1_;
	// 		// rst_tmp.records_[i].item1_=hrd01.rst_density_.records_[i].item1_;
	// 		// return 0;
	// 	}
	// 	double IQR=rst_tmp.GetQuantile(0.75)-rst_tmp.GetQuantile(0.25);
    // 	double thresh=rst_tmp.GetQuantile(0.75)+IQR*1;
	// 	for(int i=0;i<rst_tmp.records_.size();i++){
	// 		if(rst_tmp.records_[i].item1_>thresh){
	// 			hrd01.status_[i]=-1;
	// 		}
	// 		else{
	// 			hrd01.status_[i]=1;
	// 		}
	// 	}
	// 	// hrd01.DemonstrateResult_Color("Result/Ignatius_color.ply",rst_tmp,"linear");
	// 	hrd01.DemonstrateResult("Result/Ignatius_MGS2.ply");
	// }	
	else if("MGS3"==Mtype){
		// cout<<"MGS3"<<endl;
		/* Tanks and Temples */
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(150,5);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");	
		// hrd01.DemonstrateResult(output_path);

		/* DTU scan*/
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(800,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);


		/* DTU furu*/		
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(150,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.8,"1,2");
		// hrd01.DemonstrateResult(output_path);

		/* THU */		
		// hrd01.store_path_=output_path;
		// hrd01.FM_Slope(150,3);
		// hrd01.FM_RegionGrowth(1.5,80.0,200.0,"-1");
		// hrd01.FM_NID(100,0.999,"1,2");
		// hrd01.DemonstrateResult(output_path);

		/* Human (OPENMVS)*/		
		// hrd01.FM_Slope(400,6);
		// hrd01.FM_RegionGrowth(2,80.0,4000.0,"-1");
		// hrd01.FM_NID(100,0.7,"1,2");
		// hrd01.DemonstrateResult(fileNoType+"_MGS1");

		/* Human (COLMAP)*/		
		// hrd01.FM_Slope(400,7);
		hrd01.FM_Slope(400,6);
		// hrd01.FM_RegionGrowth(2,80.0,4000.0,"-1");		
		// // hrd01.FM_NID(100,0.9,"1,2");
		hrd01.FM_LoOP("1");
		hrd01.DemonstrateResult(output_path);
	}
	else if("MGS4"==Mtype){
		/* THU */		
		hrd01.FM_MEval(400,1.5);
		hrd01.FM_MajorityVote(4000,0.01,"-1");
		hrd01.FM_D4(500,30,"1,2");
		hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		hrd01.FM_LoOP("1,2,3,4");
		hrd01.DemonstrateResult(output_path);

		/* Tanks and Temples*/
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.store_path_=output_path;
		// hrd01.FM_MEval(100,9);
		// hrd01.FM_MajorityVote(4000,"-1");
		// hrd01.FM_D4(500,30,"1,2");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		// hrd01.FM_NID(30,0.7,"1,2,3,4");
		// hrd01.DemonstrateResult(fileNoType+"_MGS2");

		/* Human (OpenMVS)*/
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.FM_MEval(200,3);
		// hrd01.FM_D4(500,30,"1");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-2");
		// hrd01.FM_NID(30,0.6,"1,2,3");
		// hrd01.DemonstrateResult(fileNoType+"_MGS2");

		/* Human (COLMAP)*/
		// HybridMethods hrd01(cloud,kdtree,&ole); 
		// hrd01.FM_MEval(150,20);
		// hrd01.FM_MajorityVote(400,0.002,"1");
		// hrd01.FM_D4(100,40,"1,2");
		// hrd01.FM_RegionGrowth(2.0,80.0,200.0,"-3");
		// hrd01.FM_LoOP("1,2,3,4");
		// hrd01.DemonstrateResult(output_path);
	}
	else if("nohair"==Mtype){
		hrd01.FM_MEval(250,6);
		hrd01.DemonstrateResult(output_path);
	}
	return 0;
}


