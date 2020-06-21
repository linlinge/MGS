#include "StyleEvaluation.h"
StyleEvaluation::StyleEvaluation(pcl::PointCloud<PointType>::Ptr cloud,
					 pcl::search::KdTree<PointType>::Ptr kdtree)
{
	cloud_=cloud;
	N_=cloud_->points.size();
	if(kdtree==NULL){
		kdtree_=pcl::search::KdTree<PointType>::Ptr(new pcl::search::KdTree<PointType>);
		kdtree_->setInputCloud(cloud);
	}
	else
		kdtree_=kdtree;
	// dmean_=ComputeMeanDistance(cloud_);

	// fast mode parameter
	double z=2.58;	// confidence 0.99
	double MOE=0.001;
	double X=z*z*0.5*0.5/(MOE*MOE);
	n_=round((N_*X)/(X+N_-1));
}

void StyleEvaluation::OutlierGradeMetric(string output_path, int K)
{
	double ratio=0;
	#pragma omp parallel for
    for(int i=0;i<cloud_->points.size();i+=100){
		vector<V3> arrows;
		V3 mevec;
		double count=0;

		/* Get local point cloud */
		vector<int> idx(K);
		vector<float> dist(K);
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		ctmp->points.resize(K);
		kdtree_->nearestKSearch(i, K, idx, dist);        
		#pragma omp parallel for
		for(int i=0;i<idx.size();i++){		
			ctmp->points[i]=cloud_->points[idx[i]];
		}
        
		/* ratio */
		mevec=GetMEvec(ctmp);		
		arrows.resize(ctmp->points.size());
		for(int j=1;j<ctmp->points.size();j++){
			arrows.push_back(V3(ctmp->points[j].x-ctmp->points[0].x,ctmp->points[j].y-ctmp->points[0].y,ctmp->points[j].z-ctmp->points[0].z));
		}
		double thresh=cos(M_PI/24.0);
		for(int j=0;j<arrows.size();j++){
			double plen=Dot(arrows[j],mevec)/mevec.GetLength()/arrows[j].GetLength();
			if(abs(plen)<thresh)
				count++;
		}
		ratio+=count*1.0/arrows.size();
    }
	ratio=1-ratio/cloud_->points.size();
	ofstream fout(output_path);
    fout<<ratio<<endl;
	fout<<ratio<<endl;
	fout.close();
}

void StyleEvaluation::UniformityMetric(string output_path, int K)
{
	vector<double> density;
	density.resize(cloud_->points.size());
	#pragma omp parallel for
	for(int i=0;i<cloud_->points.size();i++){
		/* Get local point cloud */
		vector<int> idx(K);
		vector<float> dist(K);
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		ctmp->points.resize(K);
		kdtree_->nearestKSearch(i, K, idx, dist);        
		density[i]=K/dist[K-1];
	}
		
	VectorNormalize(density);
	VectorWrite(output_path,density,"cover");
}

void StyleEvaluation::SigularityMetric(string output_path, int K)
{	
	vector<double> sg;
	sg.resize(cloud_->points.size());
	int N=cloud_->points.size();

	#pragma omp parallel for
	for(int i=0;i<cloud_->points.size();i++){
		/* Get local point cloud */
		vector<int> idx(K);
		vector<float> dist(K);
		pcl::PointCloud<PointType>::Ptr ctmp(new pcl::PointCloud<PointType>);
		ctmp->points.resize(K);
		kdtree_->nearestKSearch(i, K, idx, dist);        		
		for(int i=0;i<idx.size();i++){		
			ctmp->points[i]=cloud_->points[idx[i]];
		}
		double meval=GetMEval(ctmp);
		sg[i]=N*pow(meval/dist[K-1],2);
	}
	VectorNormalize(sg);
	int count=0;
	double IQR=VectorQuantile(sg,0.75)-VectorQuantile(sg,0.25);
	// cout<<"IQR:"<<IQR<<endl;
	double thresh=IQR*1.5+VectorQuantile(sg,0.75);
	// cout<<"thresh:"<<thresh<<endl;
	// #pragma omp parallel for
	for(int i=0;i<N;i++){
		if(sg[i]>thresh)
			count++;
	}

	double ratio=count*1.0/N;
	ofstream fout(output_path);
    fout<<ratio<<endl;
	fout<<ratio<<endl;
	fout.close();
}

double StyleEvaluation::GetDmean(int K,string str)
{
	if("Common"==str){
		rst_dmean.Resize(cloud_->points.size());
		for(int i=0;i<cloud_->points.size();i++){
			double dmean_tmp=ComputeNearestDistance(kdtree_,i);
			rst_dmean.records_[i].id_=i;
			rst_dmean.records_[i].item1_=dmean_tmp;
		}
		dmean_=rst_dmean.GetMean();
		return dmean_;
	}
	else if("Fast"==str){
		srand((int)time(0));
		rst_dmean.Resize(n_);
		for(int i=0;i<n_;i++){
			int itmp=rand()%N_;
			double dmean_tmp=ComputeNearestDistance(kdtree_,itmp);
			rst_dmean.records_[i].id_=itmp;
			rst_dmean.records_[i].item1_=dmean_tmp;
		}
		dmean_=rst_dmean.GetMean();
		return dmean_;
	}
	else{
		cout<<"Dmean Mode Error!"<<endl;
		return -1;
	}
}

double StyleEvaluation::GetMinorEval(int K,string str)
{
	
	if("Common"==str){
		if(rst_meval.GetSize()==0)
			rst_meval.Resize(cloud_->points.size());
	
		#pragma omp parallel for
		for(int i=0;i<cloud_->points.size();i++){
			double meval_tmp=GetMEvalRadius(cloud_,dmean_*10.0,kdtree_,i);
			rst_meval.records_[i].id_=i;
			rst_meval.records_[i].item1_=meval_tmp;
		}
		rst_meval.Write("Result/meval.csv");
		meval_=rst_meval.GetMean();
		double meval_Q1=rst_meval.GetQuantile(0.25);
		meval_Q3_=rst_meval.GetQuantile(0.75);
		meval_IQR_=meval_Q3_-meval_Q1;
		og_=meval_/(dmean_*dmean_);
		return meval_;
	}
	else if("Fast"==str){
		if(rst_meval.GetSize()!=0){
			rst_meval.Clear();
		}
		rst_meval.Resize(n_);
		srand((int)time(0)); 

		#pragma omp parallel for
		for(int i=0;i<n_;i++){
			int itmp=rand()%N_;
			double meval_tmp=GetMEvalRadius(cloud_,dmean_*10.0,kdtree_,itmp);
			double dmean_tmp=ComputeNearestDistance(kdtree_,itmp);
			rst_meval.records_[i].id_=itmp;
			rst_meval.records_[i].item1_=meval_tmp;
		}
		meval_=rst_meval.GetMean();
		return meval_;
	}
	else
		cout<<"Minor Eigenvalue Mode Error!"<<endl;
}

double StyleEvaluation::GetDB2(int K, string str)
{
	
	if("Common"==str){
		if(rst_db2.GetSize()==0)
			rst_db2.Resize(cloud_->points.size());
		
		#pragma omp parallel for
		for(int i=0;i<cloud_->points.size();i++){
			double db2_tmp=ComputeDB2(kdtree_,100,i);
			rst_db2.records_[i].id_=i;
			rst_db2.records_[i].item1_=db2_tmp;
		}
		meval_=rst_db2.GetMean();


		double db2_Q1=rst_db2.GetQuantile(0.25);
		db2_Q3_=rst_db2.GetQuantile(0.75);
		db2_IQR_=db2_Q3_-db2_Q1;

		return meval_;
	}
	else if("Fast"==str){
		if(rst_db2.GetSize()!=0)
			rst_db2.Clear();
		rst_db2.Resize(n_);

		#pragma omp parallel for
		for(int i=0;i<n_;i++){
			int itmp=rand()%N_;
			double db2_tmp=ComputeDB2(kdtree_,100,itmp);
			rst_db2.records_[i].id_=i;
			rst_db2.records_[i].item1_=db2_tmp;
		}
		meval_=rst_db2.GetMean();
		return meval_;	
	}
	else{
		cout<<"DB2 Mode Error!"<<endl;
	}	
}