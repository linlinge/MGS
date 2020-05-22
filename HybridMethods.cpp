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
void HybridMethods::GetScopeIndices(int st,vector<int>& cIdx)
{
    cIdx.clear();
    for(int i=0;i<status_.size();i++){
        if(status_[i]==st)
            cIdx.push_back(i);
    }     
}

void HybridMethods::FM_Prox(int K, double kIQR,int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    // Default Scope
    if(domain_xi_1==FULL_DOMAIN){
        cout<<"Prox 1 Start"<<endl;
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

        /* Result Classification */
        double IQR=rst_slope_.GetQuantile(0.75)-rst_slope_.GetQuantile(0.25);
        double thresh=rst_slope_.GetQuantile(0.75)+IQR*kIQR;
        for(int i=0;i<rst_slope_.GetSize();i++){
            if(rst_slope_.records_[i].item1_>thresh){
                status_[i]=1;
                rst1->points.push_back(cloud_->points[i]);
            }
            else{
                rst0->points.push_back(cloud_->points[i]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"Prox 1 end!"<<endl;
    }
    else{
        cout<<"Prox 2 start!"<<endl;
        vector<int> scope;
        GetScopeIndices(domain_xi_1,scope);

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
                status_[scope[i]]+=accumulator_; 
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                rst0->points.push_back(cloud_->points[scope[i]]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"Prox 2 end!"<<endl;
    }
    accumulator_++;
}

void HybridMethods::FM_MEval(int K, double alpha,int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    if(domain_xi_1==FULL_DOMAIN){
        cout<<"MEval 1 start!"<<endl;
        rst_meval_.Clear();
        rst_meval_.Resize(cloud_->points.size());
        #pragma omp parallel for
        for(int i=0;i<cloud_->points.size();i++){
            double eval_tmp=GetMEvalRadius(cloud_,K*pOle_->dmean_,kdtree_,i);
            // double eval_tmp=GetMEvalKNeighours(cloud_,K,kdtree_,i);
            rst_meval_.records_[i].id_=i;
            rst_meval_.records_[i].item1_=eval_tmp;
        }
        // rst_meval_.Standardize_Zscore();
        // rst_meval_.Normalize_Tanh();
        double IQR=rst_meval_.GetQuantile(0.75)-rst_meval_.GetQuantile(0.25);
        double thresh=rst_meval_.GetQuantile(0.75)+IQR*alpha;
        // #pragma omp parallel for
        for(int i=0;i<N_;i++){
            if(rst_meval_.records_[i].item1_>thresh){
                status_[i]+=accumulator_;
                rst1->points.push_back(cloud_->points[i]);
            }
            else{
                rst0->points.push_back(cloud_->points[i]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"MEval 1 end!"<<endl;
    }
    else{
        cout<<"MEval 2 start!"<<endl;
        vector<int> scope;
        GetScopeIndices(domain_xi_1,scope);
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
                status_[scope[i]]+=accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                rst0->points.push_back(cloud_->points[scope[i]]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"MEval 2 end!"<<endl;
    }
    accumulator_++;
}

void HybridMethods::FM_NID(int K, double p, int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
   if(domain_xi_1==FULL_DOMAIN){
        cout<<"NID 1 start!"<<endl;
        rst_nid_.Clear();
        rst_nid_.Resize(cloud_->points.size());
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
            rst_nid_.records_[i].id_=i;
            rst_nid_.records_[i].item1_=atmp(0,0);
        }

        rst_nid_.LocalFilter("average",cloud_,30);
        rst_nid_.Standardize_Zscore();
        
        for(int i=0;i<rst_nid_.GetSize();i++){
            double etmp=GaussErrorFunction(rst_nid_.records_[i].item1_);
            if(etmp>p){
                status_[i]=1;                
                rst1->points.push_back(cloud_->points[i]);
            }
            else{
                rst0->points.push_back(cloud_->points[i]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"NID 1 end!"<<endl;
    }
    else{
        cout<<"NID 2 start!"<<endl;
        vector<int> scope;
        GetScopeIndices(domain_xi_1,scope);
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
                status_[scope[i]]+=accumulator_;
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                rst0->points.push_back(cloud_->points[scope[i]]);
            }            
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"NID 2 end!"<<endl;
    }
    accumulator_++;     
}


void HybridMethods::FM_DB2(int K, double P,int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    rst_db2_.Clear();
    if(domain_xi_1!=FULL_DOMAIN){
        cout<<"DB2 1 start!"<<endl;
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
                rst1->points.push_back(cloud_->points[i]);
            }
            else{
                rst0->points.push_back(cloud_->points[i]);
            }
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"DB2 1 end!"<<endl;
    }
    else{
        cout<<"DB2 2 start!"<<endl;
        vector<int> scope;
        GetScopeIndices(domain_xi_1,scope);
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
                rst1->points.push_back(cloud_->points[scope[i]]);
            }
            else{
                rst0->points.push_back(cloud_->points[scope[i]]);
            }        
        }
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
        pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
        cid_++;
        cout<<"DB2 2 end!"<<endl;
    }
}

void HybridMethods::FM_Density(int K, double alpha,int domain_xi_1)
{
    cout<<"Desisty 1 start!"<<endl;
    domain_xi_1_=domain_xi_1;
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
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
    cid_++;
    cout<<"Desisty 1 end!"<<endl;
}

void HybridMethods::FM_RegionGrowth(double thresh_eclidean, double thresh_tolerance,double thresh_kIQR, int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    cout<<"RegionGrowth start!"<<endl;
    vector<int> scope;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);    
    pcl::IndicesClustersPtr clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr small_clusters (new pcl::IndicesClusters);
    pcl::IndicesClustersPtr large_clusters (new pcl::IndicesClusters);
    pcl::search::KdTree<PointType>::Ptr kdtree_tmp(new pcl::search::KdTree<PointType>);
    pcl::PointCloud<PointType>::Ptr cloud_tmp(new pcl::PointCloud<PointType>);
    GetScopeIndices(domain_xi_1,scope);

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
    double IQR=Quantile(cluster_size_set,0.75)-Quantile(cluster_size_set,0.25);
    double thresh=Quantile(cluster_size_set,0.75)+IQR*thresh_kIQR;
    for(int i=0;i<clusters->size();i++){
        int current_cluster_size=(*clusters)[i].indices.size();
        if(current_cluster_size>thresh){            
            for(int j=0;j<(*clusters)[i].indices.size();j++){
                int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);
                status_[itmp]+=accumulator_;
                rst1->points.push_back(cloud_->points[itmp]);
            }
        }
        else{
             for(int j=0;j<(*clusters)[i].indices.size();j++){
                int itmp=GetIndex(kdtree_,cloud_tmp->points[(*clusters)[i].indices[j]]);                
                rst0->points.push_back(cloud_->points[itmp]);
            }
        }     
    }
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
    cid_++;
    cout<<"RegionGrowth end!"<<endl;
    accumulator_++;
}

void HybridMethods::FM_MajorityVote(int K,int domain_xi_1)
{
    domain_xi_1_=domain_xi_1;
    cout<<"MajorityVote start!"<<endl;
    pcl::PointCloud<PointType>::Ptr rst0(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst1(new pcl::PointCloud<PointType>);
    vector<int> scope;
    GetScopeIndices(domain_xi_1,scope);
    Table<Rrd1> rst_mv_;
    rst_mv_.Resize(scope.size());

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
        if(rst_mv_.records_[i].item1_<thresh){
            status_[scope[i]]+=2;
            rst1->points.push_back(cloud_->points[scope[i]]);
        }
        else{
            rst0->points.push_back(cloud_->points[scope[i]]);
        }      
    }
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_0.ply",*rst0);
    pcl::io::savePLYFileBinary("Result/rst_"+to_string(cid_)+"_1.ply",*rst1);
    cid_++;
    cout<<"MajorityVote end!"<<endl;

    // rst_mv_.GetCorrespondingColor();
    // for(int i=0;i<scope.size();i++){
    //     V3 ctmp=rst_mv_.color_[i];
    //     cloud_->points[scope[i]].r=ctmp.r;
    //     cloud_->points[scope[i]].g=ctmp.g;
    //     cloud_->points[scope[i]].b=ctmp.b;
    // }
    // pcl::io::savePLYFileBinary("Result/rst_cloud.ply",*cloud_);
}

void HybridMethods::DemonstrateResult(string path,int mode)
{
    pcl::PointCloud<PointType>::Ptr rst_regular_cloud(new pcl::PointCloud<PointType>);
    pcl::PointCloud<PointType>::Ptr rst_irregular_cloud(new pcl::PointCloud<PointType>);
    if(mode==STATUS){
        double max_tmp=*max_element(status_.begin(),status_.end());;
        if(max_tmp==1){
             for(int i=0;i<status_.size();i++){
                if(status_[i]==max_tmp){                    
                    rst_irregular_cloud->points.push_back(cloud_->points[i]);
                }
                else{
                    rst_regular_cloud->points.push_back(cloud_->points[i]);
                }
            }
            pcl::io::savePLYFileBinary(path+"_regular_cloud.ply",*rst_regular_cloud);
            pcl::io::savePLYFileBinary(path+"_irregular_cloud.ply",*rst_irregular_cloud);
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
}


void HybridMethods::J(int domain_xi_2)
{
  if(domain_xi_1_== FULL_DOMAIN){
      if(domain_xi_2==REGULAR_DOMAIN){
            #pragma omp parallel for
            for(int i=0;i<status_.size();i++){
                if(status_[i]!=2)
                    status_[i]=1;
                else
                    status_[i]=0;
            }
      }
      else{
           #pragma omp parallel for
           for(int i=0;i<status_.size();i++){
                if(status_[i]>0)
                    status_[i]=1;
                else
                    status_[i]=0;
           }
      }        
  }
  else if(domain_xi_1_==REGULAR_DOMAIN){
      if(domain_xi_2==REGULAR_DOMAIN){
         #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==0 || status_[i]==1)
                status_[i]=1;
            else
                status_[i]=0;
         }
      }
      else{
         #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==1 || status_[i]==2)
                status_[i]=1;
            else
                status_[i]=0;
         }
      }
  }
  else if(domain_xi_1_==IRREGULAR_DOMAIN && domain_xi_2==REGULAR_DOMAIN){
    #pragma omp parallel for
    for(int i=0;i<status_.size();i++){
        if(status_[i]==0 || status_[i]==1)
            status_[i]=1;
        else
            status_[i]=0;
    }
  }
  else
  {
      cout<<"J() domain error!"<<endl;
  }  
  accumulator_=1;
}

void HybridMethods::D(int domain_xi_2)
{
    if(domain_xi_1_==FULL_DOMAIN && domain_xi_2==REGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==3)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else if(domain_xi_1_==FULL_DOMAIN && domain_xi_2==IRREGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==1)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else if(domain_xi_1_==REGULAR_DOMAIN && domain_xi_2==REGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==2)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else if(domain_xi_1_==REGULAR_DOMAIN && domain_xi_2==IRREGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==1 || status_[i]==2)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else if(domain_xi_1_==IRREGULAR_DOMAIN && domain_xi_2==REGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==3)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else if(domain_xi_1_==IRREGULAR_DOMAIN && domain_xi_2==IRREGULAR_DOMAIN){
        #pragma omp parallel for
         for(int i=0;i<status_.size();i++){
            if(status_[i]==1)
                status_[i]=1;
            else
                status_[i]=0;
         }
    }
    else
    {
        cout<<"D() domain error!"<<endl;
    }    
    accumulator_=1;
}

void HybridMethods::U(int domain_xi_2)
{
    if(domain_xi_1_==FULL_DOMAIN && domain_xi_2==REGULAR_DOMAIN){
        #pragma omp parallel for
        for(int i=0;i<status_.size();i++){
            if(status_[i]==1)
            status_[i]=1;
            else
            status_[i]=0;
        } 
    }
    else if(domain_xi_1_==FULL_DOMAIN && domain_xi_2==IRREGULAR_DOMAIN){
        #pragma omp parallel for
        for(int i=0;i<status_.size();i++){
            if(status_[i]==3)
            status_[i]=1;
            else
            status_[i]=0;
        } 
    }
    else
    {
        cout<<"U() domain error!"<<endl;
    }
    accumulator_=1;
}

