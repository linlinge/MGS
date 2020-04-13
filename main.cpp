/*
result= A-(A and B)
*/
#include <iostream>	
#include "PCLExtend.h"
#include "V3.hpp"
#include "DataMining.h"
#include "ImprovedLoop.h"
#include "EdgeDetection.h"
#include <Eigen/Dense>
#include <omp.h>
#include "AngleBasedOutlier.h"
#include "PointSetFeatures.h"
#include "Table.h"
#include <algorithm>
#include "VectorExtend.h"
#include "SignalProcessing.h"
#include "OLEModule.h"
#include "FEModule.h"
#include<time.h>

void Method01_pulse(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,100,"pulse");	
	// psf.rst_pulse_.Standardize_Zscore();
	psf.rst_pulse_.Normalize_Tanh();
	psf.rst_pulse_.Write("Result/rst_color.csv");

	psf.rst_pulse_.GetCorrespondingColor();
	for(int i=0;i<cloud->points.size();i++)
	{
		V3 ctmp=psf.rst_pulse_.color_[i];
		cloud->points[i].r=ctmp.r;
		cloud->points[i].g=ctmp.g;
		cloud->points[i].b=ctmp.b;
	}
	pcl::io::savePLYFileBinary("Result/rst_color.ply",*cloud);
}
void Method02_Slope(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,100,"slope");
	psf.rst_slope_.Write("Result/rst_color.csv");	
	psf.rst_slope_.Standardize_Zscore();
	psf.rst_slope_.Normalize_Tanh();
	

	psf.rst_slope_.GetCorrespondingColor();
	for(int i=0;i<cloud->points.size();i++)
	{
		V3 ctmp=psf.rst_slope_.color_[i];
		cloud->points[i].r=ctmp.r;
		cloud->points[i].g=ctmp.g;
		cloud->points[i].b=ctmp.b;
	}
	pcl::io::savePLYFileBinary("Result/rst_color.ply",*cloud);
}
void Method02_Slope_PDF(pcl::PointCloud<PointType>::Ptr cloud,string str, double p)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,200,"slope");
	// psf.rst_slope_.Write("Result/slope.csv");
	// psf.rst_slope_.WritePDF();
	double t=psf.rst_slope_.ReversePDF(p);

	for(int i=0;i<psf.rst_slope_.records_.size();i++){
		if(psf.rst_slope_.records_[i].item1_>t){
			int itmp=psf.rst_slope_.records_[i].id_;
			cloud->points[itmp].r=255;
			cloud->points[itmp].g=0;
			cloud->points[itmp].b=0;
		}
	}
	pcl::io::savePLYFileBinary(str,*cloud);
}
void Method02_Slope_db_PDF(pcl::PointCloud<PointType>::Ptr cloud,string str,double p)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,200,"slope_db2");
	double t=psf.rst_slope_.ReversePDF(p);

	for(int i=0;i<psf.rst_slope_.records_.size();i++){
		if(psf.rst_slope_.records_[i].item1_>t){
			int itmp=psf.rst_slope_.records_[i].id_;
			cloud->points[itmp].r=255;
			cloud->points[itmp].g=0;
			cloud->points[itmp].b=0;
		}
	}
	pcl::io::savePLYFileBinary(str,*cloud);
}

void Method02_Slope_CLT_PDF(pcl::PointCloud<PointType>::Ptr cloud,string str, double p)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,100,"slope_CLT");
	// psf.rst_slope_.Write("Result/slope.csv");
	// psf.rst_slope_.WritePDF();
	double t=psf.rst_slope_.ReversePDF(p);

	for(int i=0;i<psf.rst_slope_.records_.size();i++){
		if(psf.rst_slope_.records_[i].item1_>t){
			int itmp=psf.rst_slope_.records_[i].id_;
			cloud->points[itmp].r=255;
			cloud->points[itmp].g=0;
			cloud->points[itmp].b=0;
		}
	}
	pcl::io::savePLYFileBinary(str,*cloud);
}
void Method03_Jump_PDF(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointSetFeatures psf;
	psf.ApplykNN(cloud,200,"gap_IQR");
}
void Method03_MinorEigenvalue(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointSetFeatures psf;
	psf.ApplyMinorEigenvalue(cloud);
	// psf.rst_MinorEigenvalue_.LocalFilter("std",cloud,80);
	psf.rst_MinorEigenvalue_.Standardize_Zscore();
	psf.rst_MinorEigenvalue_.Normalize_Tanh();
	psf.rst_MinorEigenvalue_.GetCorrespondingColor();
	#pragma omp parallel for
	for(int i=0;i<cloud->points.size();i++){
		V3 ctmp=psf.rst_MinorEigenvalue_.color_[i];
		cloud->points[i].r=ctmp.r;
		cloud->points[i].g=ctmp.g;
		cloud->points[i].b=ctmp.b;
	}
	pcl::io::savePLYFileBinary("Result/rst_color.ply",*cloud);
}
void Method03_MinorEigenvalue_PDF(pcl::PointCloud<PointType>::Ptr cloud,string str,double p)
{
	PointSetFeatures psf;
	psf.ApplyMinorEigenvalue(cloud);
	double t=psf.rst_MinorEigenvalue_.ReversePDF(p);
	for(int i=0;i<psf.rst_MinorEigenvalue_.records_.size();i++){
		if(psf.rst_MinorEigenvalue_.records_[i].item1_>t){
			int itmp=psf.rst_MinorEigenvalue_.records_[i].id_;
			cloud->points[itmp].r=255;
			cloud->points[itmp].g=0;
			cloud->points[itmp].b=0;
		}
	}
	pcl::io::savePLYFileBinary(str,*cloud);
}
void Method04_ThinCounter(pcl::PointCloud<PointType>::Ptr cloud,string str)
{
	PointSetFeatures psf;
	psf.ApplyThinPlaneAngle(cloud);
	psf.rst_ThinPlaneAngle_.GetCorrespondingColor();
	for(int i=0;i<cloud->points.size();i++){
		V3 ctmp=psf.rst_ThinPlaneAngle_.color_[i];
		cloud->points[i].r=ctmp.r;
		cloud->points[i].g=ctmp.g;
		cloud->points[i].b=ctmp.b;
	}
	pcl::io::savePLYFileBinary(str,*cloud);
}
void Blending01_MinorEigenvalue_Slope(pcl::PointCloud<PointType>::Ptr cloud)
{
	PointSetFeatures psf;
	psf.ApplyMinorEigenvalue(cloud,50);
	psf.rst_MinorEigenvalue_.Standardize_Zscore();
	psf.rst_MinorEigenvalue_.Normalize_Tanh();
	psf.rst_MinorEigenvalue_.SetActiveIndex("1000","believable");

	psf.ApplykNN(cloud,80,"slope");	
	psf.rst_slope_.Standardize_Zscore();
	psf.rst_slope_.Normalize_Tanh();
	psf.rst_slope_.SetActiveIndex("1110","believable");

	/*
		Step 03: Blending outlier
	*/
	vector<int> oidx;
	std::set_difference(psf.rst_slope_.unbelievable_idx_.begin(),
						psf.rst_slope_.unbelievable_idx_.end(),
						psf.rst_MinorEigenvalue_.believable_idx_.begin(),
						psf.rst_MinorEigenvalue_.believable_idx_.end(),
						std::back_inserter(oidx));
	vector<int> nidx;
	vector<int> idx;

	// #pragma omp parallel for
	for(int i=0;i<cloud->points.size();i++){
		idx.push_back(i);
	}
	std::set_difference(idx.begin(),
						idx.end(),
						oidx.begin(),
						oidx.end(),
						std::back_inserter(nidx));

	/*
		Step 04: output result
	*/
	pcl::PointCloud<PointType>::Ptr rst_cloud(new pcl::PointCloud<PointType>);

	// #pragma omp parallel for
	for(int i=0;i<nidx.size();i++){
		int tmp_idx=nidx[i];
		rst_cloud->points.push_back(cloud->points[tmp_idx]);
	}
	pcl::io::savePLYFileBinary("Result/rst.ply",*rst_cloud);

	#pragma omp parallel for
	for(int i=0;i<oidx.size();i++){
		int tmp_idx=oidx[i];
		cloud->points[tmp_idx].r=255;
		cloud->points[tmp_idx].g=0;
		cloud->points[tmp_idx].b=0;
	}
	pcl::io::savePLYFileBinary("Result/rst_color.ply",*cloud);

}
// void Blending01_MinorEigenvalue_Slope_PDF(pcl::PointCloud<PointType>::Ptr cloud,string str)
// {
// 	PointSetFeatures psf;
// 	vector<int> bidx,uidx;
// 	psf.ApplyMinorEigenvalue(cloud,50);
// 	psf.rst_MinorEigenvalue_.GetSpecifiedIndex(bidx,"<",0.9);
	
// 	psf.ApplykNN(cloud,100,"slope");	
// 	psf.rst_slope_.GetSpecifiedIndex(uidx,">",0.9);

// 	/*
// 		Step 03: Blending outlier
// 	*/
// 	vector<int> oidx;
// 	std::set_difference(uidx.begin(),uidx.end(),bidx.begin(),bidx.end(),std::back_inserter(oidx));

// 	/*
// 		Step 04: output result
// 	*/
// 	// #pragma omp parallel for
// 	for(int i=0;i<oidx.size();i++){
// 		int tmp_idx=oidx[i];
// 		cloud->points[tmp_idx].r=255;
// 		cloud->points[tmp_idx].g=0;
// 		cloud->points[tmp_idx].b=0;
// 	}
// 	pcl::io::savePLYFileBinary(str,*cloud);

// }

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

	// Blending01_MinorEigenvalue_Slope(cloud);
	// Blending01_MinorEigenvalue_Slope_PDF(cloud,"Result/rst_color.ply");

	// Method02_Slope(cloud);
	// Method02_Slope_PDF(cloud,"Result/rst_color.ply",0.95);
	// Method02_Slope_CLT_PDF(cloud,"Result/slope_pdf.ply",0.98);
	// Method02_Slope_db_PDF(cloud,"Result/slope_pdf.ply",0.98);

	// Method03_MinorEigenvalue(cloud);
	// Method03_MinorEigenvalue_PDF(cloud,"Result/rst_color.ply",0.75);
	// Method03_Jump_PDF(cloud);
	// Method04_ThinCounter(cloud,"Result/thincounter.ply");

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
