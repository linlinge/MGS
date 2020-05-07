/*
    Description: Feature Extraction Module
    Author: linlinge
    Date: 2020.04.05
*/
#pragma once
#include "Table.h"
#include "OLEModule.h"
#include "SignalProcessing.h"
#define SLOPE_BIN   0x01
#define SLOPE_COL   0x02
#define MEVAL_BIN   0x04
#define MEVAL_COL   0x08
#define DENSITY_COL 0x10
#define STATUS      0xFF
#define default_layer -1
extern pcl::PointCloud<PointType>::Ptr cloud_;
extern pcl::search::KdTree<PointType>::Ptr kdtree_;
extern OLEModule* pOle_;
extern vector<int> status_;
void FEModule(pcl::PointCloud<PointType>::Ptr cloud,
              pcl::search::KdTree<PointType>::Ptr kdtree,
              OLEModule* pOle_);
void ApplySlope(int K, double p,int active_layer);
void ApplyMEval(int K, double p,int active_layer);
void ApplyDB2(int K,double P,int active_layer);
void ApplyDensity(int K, double alpha,int active_layer);
void ApplyRegionGrowth(int active_layer,double minimum_clusters_ratio=0.001,double maximum_clusters_ratio=0.1);
void ApplyMajorityVote(int active_layer,int K);

/*
    Status Operations
*/
void GetScopeIndices(int st,vector<int>& cIdx);
void DemonstrateResult(string str,int mode=STATUS);