/*
    Description: Feature Extraction Module
    Author: linlinge
    Date: 2020.04.05
*/
#include "Table.h"
#include "OLEModule.h"
class FEModule
{
    public:
        Table<Rrd1> rst_slope_;
        vector<int> idx_rg_in_,oidx_;
        pcl::PointCloud<PointType>::Ptr cloud_;
        pcl::search::KdTree<PointType>::Ptr kdtree_;
        FEModule(pcl::PointCloud<PointType>::Ptr cloud,pcl::search::KdTree<PointType>::Ptr kdtree){
            cloud_=cloud;
            kdtree_=kdtree;
        };
        void ApplyHierarchical(OLEModule& ole,int K=100, double p=0.95,double thresh=0.001);
        void ApplySlope(OLEModule& ole,int K, double p);
        void ApplyMEval(OLEModule& ole);
        void DemonstrateResult(string str);
};

void FEModule::ApplyHierarchical(OLEModule& ole,int K, double p, double small_clusters_quantity_thresh)
{
    rst_slope_.Resize(cloud_->points.size());
    Eigen::MatrixXd xtmp(K,1);
    for(int i=0;i<K;i++) xtmp(i,0)=i;

    // Init ytmp
    #pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);
        Eigen::MatrixXd ytmp(K,1);
        kdtree_->nearestKSearch(cloud_->points[i], K, idx, dist);
        for(int j=0;j<dist.size();j++) ytmp(j,0)=dist[j];
        Eigen::MatrixXd atmp=(xtmp.transpose()*xtmp).inverse()*xtmp.transpose()*ytmp;
        rst_slope_.records_[i].id_=i;
        rst_slope_.records_[i].item1_=atmp(0,0);
    }

    // Get all outliers indices with slope
    double t=rst_slope_.ReversePDF(p);
    Table<Rrd1> rst_lv1;
    for(int i=0;i<rst_slope_.records_.size();i++){
		if(rst_slope_.records_[i].item1_>t){
            double meval_tmp=GetMEval(cloud_,30,kdtree_,i);
            rst_lv1.push_back(Rrd1(i,meval_tmp));
            cloud_->points[i].r=255;
            cloud_->points[i].g=0;
            cloud_->points[i].b=0;
        }            		
	}

    // double thresh=ole.meval_Q3_+ole.meval_IQR_*3.0;
    // for(int i=0;i<rst_lv1.records_.size();i++){
    //     if(rst_lv1.records_[i].item1_<thresh)
    //         idx_rg_in_.push_back(rst_lv1.records_[i].id_);
    //     else
    //         oidx_.push_back(rst_lv1.records_[i].id_);
    // }



    // // lv2 store
    // pcl::PointCloud<PointType>::Ptr lv2_cloud(new pcl::PointCloud<PointType>);
    // for(int i=0;i<idx_rg_in_.size();i++)
    //     lv2_cloud->points.push_back(cloud_->points[idx_rg_in_[i]]);
    
    // RegionGrowth rg(cloud_,lv2_cloud,kdtree_,ole.dmean_,small_clusters_quantity_thresh);
    // // for(int i=0;i<rg.oidx_.size();i++)
    // //     oidx_.push_back(rg.oidx_[i]);
    // for(int i=0;i<oidx_.size();i++){
    //     cloud_->points[oidx_[i]].r=255;
    //     cloud_->points[oidx_[i]].g=0;
    //     cloud_->points[oidx_[i]].b=0;
    // }

    // for(int i=0;i<rg.oidx_.size();i++){
    //     cloud_->points[rg.oidx_[i]].r=0;
    //     cloud_->points[rg.oidx_[i]].g=255;
    //     cloud_->points[rg.oidx_[i]].b=0;
    // }

}


void FEModule::ApplySlope(OLEModule& ole,int K, double p)
{
    rst_slope_.Resize(cloud_->points.size());
    Eigen::MatrixXd xtmp(K,1);
    for(int i=0;i<K;i++) xtmp(i,0)=i;

    // Init ytmp
    #pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i++){
        vector<int> idx(K);
        vector<float> dist(K);
        Eigen::MatrixXd ytmp(K,1);
        kdtree_->nearestKSearch(cloud_->points[i], K, idx, dist);
        for(int j=0;j<dist.size();j++) ytmp(j,0)=dist[j];
        Eigen::MatrixXd atmp=(xtmp.transpose()*xtmp).inverse()*xtmp.transpose()*ytmp;
        rst_slope_.records_[i].id_=i;
        rst_slope_.records_[i].item1_=atmp(0,0);
    }

    // Get all outliers indices with slope
    // double t=rst_slope_.ReversePDF(p);
    // Table<Rrd1> rst_lv1;
    // for(int i=0;i<rst_slope_.records_.size();i++){
	// 	if(rst_slope_.records_[i].item1_>t){
    //         cloud_->points[i].r=255;
    //         cloud_->points[i].g=0;
    //         cloud_->points[i].b=0;
    //     }            		
	// }

    rst_slope_.Standardize_Zscore();
    rst_slope_.Normalize_Tanh();
    rst_slope_.GetCorrespondingColor();
    for(int i=0;i<rst_slope_.color_.size();i++){
        V3 ctmp=rst_slope_.color_[i];
        cloud_->points[i].r=ctmp.r;
        cloud_->points[i].g=ctmp.g;
        cloud_->points[i].b=ctmp.b;
    }
}

void FEModule::ApplyMEval(OLEModule& ole)
{
    Table<Rrd1> rst_meval_;
    rst_meval_.Resize(cloud_->points.size());
	#pragma omp parallel for
	for(int i=0;i<cloud_->points.size();i++){
        double eval_tmp=GetMEval(cloud_,ole.dmean_*10,kdtree_,i);
		rst_meval_.records_[i].id_=i;
		rst_meval_.records_[i].item1_=eval_tmp;
	}

    // rst_meval_.Standardize_Zscore();
    // rst_meval_.Multiply(1000);
    // rst_meval_.Normalize_Tanh();
    rst_meval_.GetCorrespondingColor();
    for(int i=0;i<rst_meval_.color_.size();i++){
        V3 ctmp=rst_meval_.color_[i];
        cloud_->points[i].r=ctmp.r;
        cloud_->points[i].g=ctmp.g;
        cloud_->points[i].b=ctmp.b;
    }
}

void FEModule::DemonstrateResult(string str)
{
	pcl::io::savePLYFileBinary(str,*cloud_);
}