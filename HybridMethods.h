/*
    Description: Hybrid Methods Generation
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
#define FULL_DOMAIN -1
#define REGULAR_DOMAIN 0
#define IRREGULAR_DOMAIN 1
extern pcl::PointCloud<PointType>::Ptr cloud_;
extern pcl::search::KdTree<PointType>::Ptr kdtree_;
extern OLEModule* pOle_;
extern vector<int> status_;
extern int accumulator;

class HybridMethods
{
    public:
        /* Variables */
        Table<Rrd1> rst_meval_;
        Table<Rrd1> rst_slope_;
        Table<Rrd1> rst_db2_;
        Table<Rrd1> rst_density_;
        Table<Rrd1> rst_nid_;
        pcl::PointCloud<PointType>::Ptr cloud_;
        pcl::search::KdTree<PointType>::Ptr kdtree_;
        OLEModule* pOle_;
        int N_;
        int accumulator_;
        int domain_xi_1_;
        int cid_;

        /* Construction */
        HybridMethods(pcl::PointCloud<PointType>::Ptr cloud,pcl::search::KdTree<PointType>::Ptr kdtree,OLEModule* pOle){
            cloud_=cloud;
            kdtree_=kdtree;
            pOle_=pOle;
            N_=cloud_->points.size();
            status_.resize(N_);
            for(int i=0;i<N_;i++){
                status_[i]=0;
            }
            accumulator_=1;
            domain_xi_1_=FULL_DOMAIN;
            cid_=0;
        }

        /* Function Modules*/
        void FM_Prox(int K, double kIQR,int domain_xi_1=FULL_DOMAIN);          // Outliers detection based on proximity
        void FM_MEval(int K, double kIQR,int domain_xi_1=FULL_DOMAIN);      // IRPC: Irregular and Regular Points Classfication based on Minor Eigenvalue
        void FM_NID(int K, double kIQR, int domain_xi_1=FULL_DOMAIN);
        void FM_DB2(int K,double P,int domain_xi_1=FULL_DOMAIN);
        void FM_Density(int K, double alpha,int domain_xi_1=FULL_DOMAIN);
        void FM_RegionGrowth(double thresh_eclidean, double thresh_tolerance,double thresh_kIQR, int domain_xi_1);
        void FM_MajorityVote(int K,int domain_xi_1=FULL_DOMAIN);

        /* Combination Operator */
        vector<int> status_;
        void U(int domain_xi_2);
        void D(int domain_xi_2);
        void J(int domain_xi_2);

        /* Status Operations */
        void GetScopeIndices(int st,vector<int>& cIdx);
        void DemonstrateResult(string str,int mode=STATUS);
};
