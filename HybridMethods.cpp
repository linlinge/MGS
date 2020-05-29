#include "HybridMethods.h"
/* Private Functions and Variables */
double cthresh=INT_MAX;
bool customRegionGrowing(const PointType& point_a, const PointType& point_b, float squared_distance)
{
  if (squared_distance < cthresh)
    return true;
  else
    return false;
}


/* Public Function */
void HybridMethods::GetScopeIndices(string str,vector<int>& cIdx)
{
    if("full"==str){
        cIdx.resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<cIdx.size();i++)
            cIdx[i]=i;
    }
    else{
        vector<int> emnt; // elements
        Str2Vec(str,",",emnt);
        VecFindPos(status_,emnt,cIdx);
    }
}


void HybridMethods::FM_Prox(int K, double kIQR,string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    cout<<"Prox start!"<<endl;
    vector<int> scope;
    GetScopeIndices(domain,scope);

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
    double thresh=rst_slope_.GetQuantile(0.75)+IQR*kIQR;
    for(int i=0;i<scope.size();i++){
        if(rst_slope_.records_[i].item1_>thresh){
            status_[scope[i]]=-accumulator_; 
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            status_[scope[i]]=accumulator_; 
            rst0->points.push_back(cloud_->points[scope[i]]);
        }
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"Prox end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_MEval(int K, double alpha,string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    cout<<"MEval start!"<<endl;
    vector<int> scope;
    GetScopeIndices(domain,scope);
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
        // double eval_tmp=GetMEvalRadius(cloud_tmp,K*pOle_->dmean_,kdtree_tmp,i);
        double eval_tmp=GetMEvalKNeighours(cloud_tmp,K,kdtree_tmp,i);
        rst_meval_.records_[i].id_=itmp;
        rst_meval_.records_[i].item1_=eval_tmp;
    }
    // rst_meval_.Standardize_Zscore();
    // rst_meval_.Normalize_Tanh();
    double IQR=rst_meval_.GetQuantile(0.75)-rst_meval_.GetQuantile(0.25);
    double thresh=rst_meval_.GetQuantile(0.75)+IQR*alpha;
    // #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        if(rst_meval_.records_[i].item1_>thresh){
            status_[scope[i]]=-accumulator_;
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            status_[scope[i]]=accumulator_;
            rst0->points.push_back(cloud_->points[scope[i]]);
        }
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"MEval end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_NID(int K, double p, string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    cout<<"NID start!"<<endl;
    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_nid_.Clear();
    rst_nid_.Resize(scope.size());
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
        rst_nid_.records_[i].id_=itmp;
        rst_nid_.records_[i].item1_=atmp(0,0);
    }

    rst_nid_.LocalFilter("average",cloud_tmp,80);
    rst_nid_.Standardize_Zscore();

    for(int i=0;i<rst_nid_.GetSize();i++){
        double etmp=GaussErrorFunction(rst_nid_.records_[i].item1_);
        if(etmp>p){
            status_[scope[i]]=-accumulator_;
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            status_[scope[i]]=accumulator_;
            rst0->points.push_back(cloud_->points[scope[i]]);
        }            
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"NID end!"<<endl;
    accumulator_++;     
}


void HybridMethods::FM_DB2(int K, double P,string domain)
{
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    rst_db2_.Clear();
    cout<<"DB start!"<<endl;
    vector<int> scope;
    GetScopeIndices(domain,scope);
    rst_db2_.Resize(scope.size());
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);

    // Init ytmp
    #pragma omp parallel for
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
            status_[scope[i]]=-accumulator_;
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            status_[scope[i]]=accumulator_;
            rst0->points.push_back(cloud_->points[scope[i]]);
        }        
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"DB end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_Density(int K, double alpha,string domain)
{
    cout<<"Desisty 1 start!"<<endl;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
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
            rst1->points.push_back(cloud_->points[i]);
        }
        else{
            rst0->points.push_back(cloud_->points[i]);
        }
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"Desisty 1 end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_RegionGrowth(double thresh_eclidean, double thresh_tolerance,double thresh_kIQR, string domain)
{
    cout<<"RegionGrowth start!"<<endl;
    vector<int> scope;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);    
    pcl::IndicesClustersPtr clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr small_clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr large_clusters (new pcl::IndicesClusters);
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    GetScopeIndices(domain,scope);

    for(int i=0;i<scope.size();i++)
        cloud_tmp->points.push_back(cloud_->points[scope[i]]);
    kdtree_tmp->setInputCloud(cloud_tmp);
    cthresh=pOle_->dmean_*thresh_eclidean;

    // Set up a Conditional Euclidean Clustering class
    pcl::ConditionalEuclideanClustering<PointType> cec (true);
    cec.setInputCloud (cloud_tmp);
    cec.setConditionFunction (&customRegionGrowing);
    cec.setClusterTolerance (pOle_->dmean_*thresh_tolerance);
    cec.segment (*clusters);
    vector<int> cluster_size_set;
    for(int i=0;i<clusters->size();i++){
       cluster_size_set.push_back((*clusters)[i].indices.size());
    }
    cout<<"cluster size:"<<cluster_size_set.size()<<endl;

    double IQR=Quantile(cluster_size_set,0.75)-Quantile(cluster_size_set,0.25);
    double thresh=Quantile(cluster_size_set,0.75)+IQR*thresh_kIQR;
    for(int i=0;i<clusters->size();i++){
        int current_cluster_size=(*clusters)[i].indices.size();
        if(current_cluster_size<=thresh){            
            for(int j=0;j<(*clusters)[i].indices.size();j++){
                int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                status_[itmp]=-accumulator_;
                rst1->points.push_back(cloud_->points[itmp]);
            }
        }
        else{
             for(int j=0;j<(*clusters)[i].indices.size();j++){
                int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                status_[itmp]=accumulator_;
                rst0->points.push_back(cloud_->points[itmp]);
            }
        }     
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"RegionGrowth end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_MajorityVote(int K,string domain)
{
    cout<<"MajorityVote start!"<<endl;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    vector<int> scope;
    GetScopeIndices(domain,scope);
    Table<Rrd1> rst_mv_;
    rst_mv_.Resize(scope.size());

    #pragma omp parallel for
    for(int i=0;i<scope.size();i++){
        int itmp=scope[i];
        vector<int> idx;
        vector<float> dist;
        kdtree_->nearestKSearch(itmp,K,idx,dist);
        // kdtree_->radiusSearch(itmp,pOle_->dmean_*200,idx,dist);

        double conf=0.0;
        double T1=5*pOle_->dmean_;
        for(int j=0;j<idx.size();j++){
            if(status_[idx[j]]>0 && dist[j]<T1)
                conf+=1.0;
        }
        conf=conf/K;
        rst_mv_.records_[i].id_=scope[i];
        rst_mv_.records_[i].item1_=conf;
    }

    // double IQR=rst_mv_.GetQuantile(0.75)-rst_mv_.GetQuantile(0.25);
    // double thresh=rst_mv_.GetQuantile(0.75)+IQR*1.5;
    // double thresh=(rst_mv_.GetMaximum()-rst_mv_.GetMinimum())*0.25+rst_mv_.GetMinimum();
    double thresh=0.02;
    for(int i=0;i<scope.size();i++){
        if(rst_mv_.records_[i].item1_<thresh){
            status_[scope[i]]=-accumulator_;
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            status_[scope[i]]=accumulator_;
            rst0->points.push_back(cloud_->points[scope[i]]);
        }      
    }
    #if STORE
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+").ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst("+to_string(accumulator_)+"-).ply",*rst1);
    #endif
    cout<<"MajorityVote end!"<<endl;
    accumulator_++;
}

void HybridMethods::DemonstrateResult(string path)
{
    pcl::PointCloud<PointType>::Ptr rst_regular_cloud(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst_irregular_cloud(new pcl::PointCloud<PointType>);
    double max_tmp=*max_element(status_.begin(),status_.end());
    VecUnique(status_);

    if(max_tmp>0){
            for(int i=0;i<status_.size();i++){
            if(status_[i]>0){                    
                rst_regular_cloud->points.push_back(cloud_->points[i]);
            }
            else{
                rst_irregular_cloud->points.push_back(cloud_->points[i]);
            }
        }
        pcl::io::savePLYFileBinary(path+"_regular_cloud.ply",*rst_regular_cloud);
        pcl::io::savePLYFileBinary(path+"_irregular_cloud.ply",*rst_irregular_cloud);
    }
    else
        cout<<"Status Error!"<<endl;
}
void HybridMethods::DemonstrateResult(string path, Table<Rrd1>& tb)
{
    tb.GetCorrespondingColor();
    for(int i=0;i<tb.GetSize();i++){
        V3 ctmp=tb.color_[i];
        cloud_->points[i].r=ctmp.r;
        cloud_->points[i].g=ctmp.g;
        cloud_->points[i].b=ctmp.b;        
    }
    pcl::io::savePLYFileBinary(path+"_color.ply",*cloud_);
}