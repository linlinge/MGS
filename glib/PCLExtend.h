#pragma once
#include <iostream>
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/octree/octree.h>
#include <pcl/point_cloud.h>
#include <pcl/common/centroid.h>
#include <pcl/features/boundary.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/segmentation/conditional_euclidean_clustering.h>
#include <omp.h>
#include <vector>
#include <Eigen/Dense>
#include "V3.hpp"
#include "VectorExtend.h"
#ifndef PointType
#define PointType pcl::PointXYZRGBNormal
#endif
using namespace std;
/*
        Compute
*/
double ComputeMeanDistance(const pcl::PointCloud<PointType>::ConstPtr cloud);
double ComputeMaxDistance(const pcl::PointCloud<PointType>::ConstPtr cloud);
double ComputeNearestDistance(pcl::search::KdTree<PointType>::Ptr kdtree,int index);
double ComputeDB2(pcl::search::KdTree<PointType>::Ptr kdtree,int K,int index);
PointType ComputeCentroid(const pcl::PointCloud<PointType>::ConstPtr cloud);
void TransformPointCloud(pcl::PointCloud<PointType>::Ptr cloud,
                         pcl::PointCloud<PointType>::Ptr cloud_tf,Eigen::Affine3f tf);

/*
        Get
*/
int GetIndex(pcl::search::KdTree<PointType>::Ptr kdtree,PointType ptmp);

/* 
        Eigenvector and Eigenvalue 
*/
double GetMEval(pcl::PointCloud<PointType>::Ptr cloud,int k,
		pcl::search::KdTree<PointType>::Ptr kdtree,int index);
double GetMEval(pcl::PointCloud<PointType>::Ptr cloud,double radius,
		pcl::search::KdTree<PointType>::Ptr kdtree,int index);
double GetMEval2(pcl::PointCloud<PointType>::Ptr cloud,int k,
		pcl::search::KdTree<PointType>::Ptr kdtree,int index);
void GetEvalAndEvec(pcl::PointCloud<PointType>::Ptr cloud,int k,
		pcl::search::KdTree<PointType>::Ptr kdtree,int index,
                vector<double>& eval,vector<V3>& evec);
void GetEvalAndEvec(pcl::PointCloud<PointType>::Ptr ctmp,vector<double>& eval,vector<V3>& evec);
void GetEval(pcl::PointCloud<PointType>::Ptr cloud,int k,
	     pcl::search::KdTree<PointType>::Ptr kdtree,int index);

/*
        Region Growth
*/
bool customRegionGrowing (const PointType& point_a, const PointType& point_b, float squared_distance);
class RegionGrowth
{
    public:
        // pcl::PointCloud<PointType>::Ptr cloud_;
        // pcl::search::KdTree<PointType>::Ptr kdtree_;
        pcl::IndicesClustersPtr clusters,small_clusters,large_clusters;
        vector<int> oidx_;
        RegionGrowth(pcl::PointCloud<PointType>::Ptr cloud, 
                     pcl::PointCloud<PointType>::Ptr rg_cloud, 
                     pcl::search::KdTree<PointType>::Ptr kdtree,
                     double eclidean_thresh,
                     double small_cluster_ratio);
};