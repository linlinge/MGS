#include "FEModule.h"
pcl::PointCloud<PointType>::Ptr cloud_;
pcl::search::KdTree<PointType>::Ptr kdtree_;
OLEModule* pOle_;
vector<int> idx_rg_in_;
vector<int> status_;
int N_=0;

/*
    Private Functions and Variables
*/
Table<Rrd1> rst_meval_;
Table<Rrd1> rst_slope_;
Table<Rrd1> rst_db2_;
Table<Rrd1> rst_density_;
double cthresh=INT_MAX;
bool customRegionGrowing (const PointType& point_a, const PointType& point_b, float squared_distance)
{
  if (squared_distance < cthresh)
    return true;
  else
    return false;
}


/*
    Public Function
*/
void FEModule(pcl::PointCloud<PointType>::Ptr cloud,pcl::search::KdTree<PointType>::Ptr kdtree,OLEModule* pOle){
    cloud_=cloud;
    kdtree_=kdtree;
    pOle_=pOle;
    N_=cloud_->points.size();
    status_.resize(N_);
    for(int i=0;i<N_;i++){
        status_[i]=0;

    }
};
/*
    Indices Operation
*/
void GetScopeIndices(int st,vector<int>& cIdx)
{
    cIdx.clear();
    for(int i=0;i<status_.size();i++){
        if(status_[i]==st)
            cIdx.push_back(i);
    }     
}
void ApplySlope(int K, double p,int active_layer)
{
    // Default Scope
    if(active_layer==default_layer){
        cout<<"Slope 1 Start"<<endl;
        rst_slope_.Clear();
        rst_slope_.Resize(cloud_->points.size());
        Eigen::MatrixXd xtmp(K,1);
        for(int i=0;i<K;i++) xtmp(i,0)=i;
        // calculate y
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
    
        // double t=rst_slope_.ReversePDF(p);
        // for(int i=0;i<N_;i++){
        //     if(rst_slope_.records_[i].item1_>t){
        //         status_[i]+=1;
        //     }
        // }
        double IQR=rst_slope_.GetQuantile(0.75)-rst_slope_.GetQuantile(0.25);
        double thresh=rst_slope_.GetQuantile(0.75)+IQR*p;
        for(int i=0;i<rst_slope_.GetSize();i++){
            if(rst_slope_.records_[i].item1_>thresh){
                status_[i]=1;
            }
        }
        cout<<"Slope 1 end!"<<endl;
    }
    else{
        cout<<"2"<<endl;
        vector<int> scope;
        GetScopeIndices(active_layer,scope);

        rst_slope_.Clear();
        rst_slope_.Resize(scope.size());
        pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
        pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++)
            cloud_tmp->points.push_back(cloud_->points[scope[i]]);
        kdtree_tmp->setInputCloud(cloud_tmp);

        Eigen::MatrixXd xtmp(K,1);
        for(int i=0;i<K;i++) xtmp(i,0)=i;
        // calculate y
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            vector<int> idx(K);
            vector<float> dist(K);
            Eigen::MatrixXd ytmp(K,1);
            int itmp=scope[i];
            kdtree_tmp->nearestKSearch(cloud_->points[itmp], K, idx, dist);
            for(int j=0;j<dist.size();j++) ytmp(j,0)=dist[j];
            Eigen::MatrixXd atmp=(xtmp.transpose()*xtmp).inverse()*xtmp.transpose()*ytmp;
            rst_slope_.records_[i].id_=itmp;
            rst_slope_.records_[i].item1_=atmp(0,0);
        }
        // double t=rst_slope_.ReversePDF(p);
        // for(int i=0;i<scope.size();i++){
        //     if(rst_slope_.records_[i].item1_>t){
        //         status_[scope[i]]+=1;
        //     }
        // }
        double IQR=rst_slope_.GetQuantile(0.75)-rst_slope_.GetQuantile(0.25);
        double thresh=rst_slope_.GetQuantile(0.75)+IQR*p;
        for(int i=0;i<scope.size();i++){
            if(rst_slope_.records_[i].item1_>thresh){
                status_[scope[i]]+=2;
            }
        }
    }
}
void ApplyMEval(int K, double alpha,int active_layer)
{
    if(active_layer==default_layer){
        rst_meval_.Clear();
        rst_meval_.Resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<cloud_->points.size();i++){
            double eval_tmp=GetMEvalRadius(cloud_,K*pOle_->dmean_,kdtree_,i);
            rst_meval_.records_[i].id_=i;
            rst_meval_.records_[i].item1_=eval_tmp;
        }
        // rst_meval_.Standardize_Zscore();
        // rst_meval_.Normalize_Tanh();
        double IQR=rst_meval_.GetQuantile(0.75)-rst_meval_.GetQuantile(0.25);
        double thresh=rst_meval_.GetQuantile(0.75)+IQR*alpha;
        #pragma omp parallel for
        for(int i=0;i<N_;i++){
            if(rst_meval_.records_[i].item1_>thresh){
                status_[i]=1;
            }
        }
    }
    else{
        cout<<"meval 2"<<endl;
        vector<int> scope;
        GetScopeIndices(active_layer,scope);
        rst_meval_.Clear();
        rst_meval_.Resize(scope.size());
        pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
        pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++)
            cloud_tmp->points.push_back(cloud_->points[scope[i]]);
        kdtree_tmp->setInputCloud(cloud_tmp);

        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            int itmp=scope[i];
            double eval_tmp=GetMEvalRadius(cloud_tmp,alpha*pOle_->dmean_,kdtree_tmp,i);
            rst_meval_.records_[i].id_=itmp;
            rst_meval_.records_[i].item1_=eval_tmp;
        }
        // rst_meval_.Standardize_Zscore();
        // rst_meval_.Normalize_Tanh();
        double IQR=rst_meval_.GetQuantile(0.75)-rst_meval_.GetQuantile(0.25);
        double thresh=rst_meval_.GetQuantile(0.75)+IQR*alpha;
        #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            if(rst_meval_.records_[i].item1_>thresh){
                status_[scope[i]]+=2;
            }
        }
    }
}

void ApplyDB2(int K, double P,int active_layer)
{
    rst_db2_.Clear();
    if(active_layer!=default_layer){
        rst_db2_.Resize(cloud_->points.size());
        // Init ytmp
        // #pragma omp parallel for
        for(int i=0;i<cloud_->points.size();i++){
            vector<int> idx(K);
            vector<float> dist(K);
            kdtree_->nearestKSearch(i, K, idx, dist); 
            vector<double> db;
            DaubechiesWavelet(dist,db);
            db.pop_back();
            rst_db2_.records_[i].id_=i;
            rst_db2_.records_[i].item1_=VectorMaximum(db);
        }

        double IQR=rst_db2_.GetQuantile(0.75)-rst_db2_.GetQuantile(0.25);
        double thresh=rst_db2_.GetQuantile(0.75)+IQR*P;
        for(int i=0;i<rst_db2_.records_.size();i++){
            if(rst_db2_.records_[i].item1_>thresh){
                status_[i]=1;
            }
        }        
    }
    else{
        vector<int> scope;
        GetScopeIndices(active_layer,scope);
        rst_db2_.Resize(scope.size());
        pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
        pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
        for(int i=0;i<scope.size();i++)
            cloud_tmp->points.push_back(cloud_->points[scope[i]]);
        kdtree_tmp->setInputCloud(cloud_tmp);

        // Init ytmp
        // #pragma omp parallel for
        for(int i=0;i<scope.size();i++){
            vector<int> idx(K);
            vector<float> dist(K);
            kdtree_tmp->nearestKSearch(cloud_->points[scope[i]], K, idx, dist); 
            vector<double> db;
            DaubechiesWavelet(dist,db);
            db.pop_back();
            rst_db2_.records_[i].id_=i;
            rst_db2_.records_[i].item1_=VectorMaximum(db);
        }
        double IQR=rst_db2_.GetQuantile(0.75)-rst_db2_.GetQuantile(0.25);
        double thresh=rst_db2_.GetQuantile(0.75)+IQR*P;
        for(int i=0;i<scope.size();i++){
            if(rst_db2_.records_[i].item1_>thresh){
                status_[scope[i]]+=2;
            }
        }
    }
}
void ApplyDensity(int K, double alpha,int active_layer)
{
    cout<<"1"<<endl;
    rst_density_.Resize(cloud_->points.size());
    #pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i++){
        double eval_tmp=GetCounterAmongRadius(cloud_,K*pOle_->dmean_,kdtree_,i);
        rst_density_.records_[i].id_=i;
        rst_density_.records_[i].item1_=eval_tmp;
    }

    double IQR=rst_density_.GetQuantile(0.75)-rst_density_.GetQuantile(0.25);
    double thresh=rst_density_.GetQuantile(0.75)+IQR*alpha;
    #pragma omp parallel for
    for(int i=0;i<N_;i++){
        if(rst_density_.records_[i].item1_>thresh){
            status_[i]+=2;
        }
    }
}
void ApplyRegionGrowth(int active_layer,double minimum_clusters_ratio,double maximum_clusters_ratio)
{
    vector<int> scope;
    pcl::IndicesClustersPtr clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr small_clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr large_clusters (new pcl::IndicesClusters);
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    GetScopeIndices(active_layer,scope);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);
    cthresh=pOle_->dmean_*2.0;

    // Set up a Conditional Euclidean Clustering class
    pcl::ConditionalEuclideanClustering<PointType> cec (true);
    cec.setInputCloud (cloud_tmp);
    cec.setConditionFunction (&customRegionGrowing);
    cec.setClusterTolerance (0.1);
    cec.setMinClusterSize (cloud_tmp->points.size () *minimum_clusters_ratio);
    cec.setMaxClusterSize (cloud_tmp->points.size () *maximum_clusters_ratio);
    cec.segment (*clusters);
    cec.getRemovedClusters (small_clusters, large_clusters);

    for (int i = 0; i < large_clusters->size (); ++i){
        for(int j = 0; j < (*large_clusters)[i].indices.size (); ++j){
            int itmp=(*large_clusters)[i].indices[j];
            int global_index=GetIndex(kdtree_,cloud_tmp->points[itmp]);
            status_[global_index]+=2;
        }
    }
}

void ApplyMajorityVote(int active_layer,int K)
{
    vector<int> scope;
    GetScopeIndices(active_layer,scope);
    Table<Rrd1> rst_mv_;
    rst_mv_.Resize(scope.size());
    cout<<"wokaka"<<endl;

    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        int itmp=scope[i];
        vector<int> idx;
        vector<float> dist;
        // kdtree_->nearestKSearch(itmp,K,idx,dist);
        kdtree_->radiusSearch(itmp,pOle_->dmean_*200,idx,dist);

        double conf=0.0;
        for(int j=0;j<idx.size();j++){
            if(status_[idx[j]]==0)
                conf+=1.0;
        }
        conf=conf/K;
        rst_mv_.records_[i].id_=scope[i];
        rst_mv_.records_[i].item1_=conf;
    }

    double IQR=rst_mv_.GetQuantile(0.75)-rst_mv_.GetQuantile(0.25);
    // double thresh=rst_mv_.GetQuantile(0.75)+IQR*3.0;
    double thresh=(rst_mv_.GetMaximum()-rst_mv_.GetMinimum())*0.25+rst_mv_.GetMinimum();
    for(int i=0;i<scope.size();i++){
        if(rst_mv_.records_[i].item1_>thresh){
            status_[scope[i]]+=2;
        }            
    }

    // rst_mv_.GetCorrespondingColor();
    // for(int i=0;i<scope.size();i++){
    //     V3 ctmp=rst_mv_.color_[i];
    //     cloud_->points[scope[i]].r=ctmp.r;
    //     cloud_->points[scope[i]].g=ctmp.g;
    //     cloud_->points[scope[i]].b=ctmp.b;
    // }
    // pcl::io::savePLYFileBinary("Result/rst_cloud.ply",*cloud_);
}

void DemonstrateResult(string path,int mode)
{
    pcl::PointCloud<PointType>::Ptr rst_cloud(new pcl::PointCloud<PointType>);
    if(mode==STATUS){
        double max_tmp=*max_element(status_.begin(),status_.end());;
        if(max_tmp!=0){
             for(int i=0;i<status_.size();i++){
                if(status_[i]==max_tmp){
                    cloud_->points[i].r=255;
                    cloud_->points[i].g=0;
                    cloud_->points[i].b=0;
                }
                else{
                    rst_cloud->points.push_back(cloud_->points[i]);
                }
            }
            pcl::io::savePLYFileBinary("Result/rst_cloud.ply",*rst_cloud);
        }
        else
            cout<<"Status Error!"<<endl;
    }
    else if(mode==MEVAL_COL){
        rst_meval_.GetCorrespondingColor();
        for(int i=0;i<rst_meval_.color_.size();i++){
            V3 ctmp=rst_meval_.color_[i];
            cloud_->points[i].r=ctmp.r;
            cloud_->points[i].g=ctmp.g;
            cloud_->points[i].b=ctmp.b;
        }
    }
    else if(mode==SLOPE_COL){
        rst_slope_.GetCorrespondingColor();
        for(int i=0;i<rst_slope_.color_.size();i++){
            V3 ctmp=rst_slope_.color_[i];
            cloud_->points[i].r=ctmp.r;
            cloud_->points[i].g=ctmp.g;
            cloud_->points[i].b=ctmp.b;
        }
    }
    else if(mode==DENSITY_COL){
        rst_density_.GetCorrespondingColor();
        for(int i=0;i<rst_density_.color_.size();i++){
            V3 ctmp=rst_density_.color_[i];
            cloud_->points[i].r=ctmp.r;
            cloud_->points[i].g=ctmp.g;
            cloud_->points[i].b=ctmp.b;
        }
    }
    
	pcl::io::savePLYFileBinary(path,*cloud_);
}