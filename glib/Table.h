#pragma once
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string>
#include "Statistics.h"
#include "PCLExtend.h"
#include "VectorExtend.h"
#include "Color.h"
using namespace std;
class Rrd1
{
	public:
		int id_;
		double item1_;
		Rrd1(){}
		Rrd1(int id,double item1){
			id_=id;
			item1_=item1;
		}
};
class Rrd2
{
    public:
        int id_;
        double item1_;
        double item2_;
        Rrd2(){id_=0;item1_=0;item2_=0;};
        Rrd2(int id,double item1,double item2){
            id_=id;
            item1_=item1_;
            item2_=item2_;
        }
};
class Cluster
{
    public:
        vector<int> c_;
};

template <class Rrd>
class Table
{
	public:
		vector<Rrd> records_;
		vector<Rrd> records_sorted_;
		vector<int> records_hist_;
		double mean_,sigma_;		
		double min_=INT_MAX;
		double max_=-INT_MAX;

		// use on hierarchical method
		vector<int> believable_idx_;	// believable area index for upper level
		vector<int> unbelievable_idx_;	// unbelivable area index for upper level
		vector<int> lower_unbelievable_idx_;	// index for outlier below lower bound
		vector<int> upper_unbelievable_idx_;	// index for outlier greater than upper bound
		vector<V3> color_;
		

		double Q1_,Q3_;
		double median_;	
		double upper_inner_limit_,lower_inner_limit_;
		int flag_is_active_=0;

		// Quantile
		double Quantile(double p);
		double Median();
		double Mean();
		double Std();
		double Minimum();
		double Maximum();
        int Size();
        void Clear();

		// 
		void EnableActive();
		void DisableActive();
		int Rows();
		void Resize(int n);
		void push_back(Rrd e);
		void Ascending(int item_index=1);
		void Descending(int item_index=1);

		void GetBoxplot(double lamda=20.0,int item_index=1);
		void SortBackup(int item_index=1);		
		void Print();
		double PDF(double t,double n,double h);
		double WritePDF();
		double ReversePDF(double P);

        // Find
		void Find(vector<int>& out,string str=">",double thresh1=0, double thresh2=0);

		// Standarize with Z-score method
		void Standardize_Zscore(int item_index=1);
		// Normalize with tanh
		void Normalize_Tanh(double scale=1.0, int item_index=1);		
		void Normalize_Min_Max(int item_index=1);
		// Write File
		void Write(string path);
		// Multiply
		void MultiplyQ6(int m1,int m2,double scale);
        void Multiply(double);
		//
		void SetQ6(int m1,int m2);
		// Overturn the table
		void Overturn();
		// error
		void GetNormalDistributionError();
		//
		void LocalFilter(string str,pcl::PointCloud<PointType>::Ptr cloud,int K);
		void GetCorrespondingColor(int item_index=1);
		//
		void GetHistogram(int k=30);
		void GaussianOutlier(int k=1);

		void nPLOF(pcl::PointCloud<PointType>::Ptr cloud,int K);
		void SetActiveIndex(string erea_color, string status);

		// Get 
		vector<int> GetIndex(string str,double v);
        void GetMeanAndVariance(int item_index=1);
		void GetMinimumAndMaximum(int item_index=1);
        double GetMininum(int item_index=1);
        double GetMean(int item_index=1);

		// Fast Remove 
		void FastRemove(int index);
		void FastRemove(vector<int> indices);
		
};
void TableHierarchy();



template <class Rrd>
void Table<Rrd>::Resize(int n)
{
    records_.resize(n);
}
template <class Rrd>
void Table<Rrd>::Clear()
{
    records_.clear();
    records_sorted_.clear();
    records_hist_.clear();
}

template <class Rrd>
int Table<Rrd>::Rows()
{
    return records_.size();
}

template <class Rrd>
int Table<Rrd>::Size()
{
    return records_.size();
}

template <class Rrd>
void Table<Rrd>::EnableActive()
{
    flag_is_active_=1;
}

template <class Rrd>
void Table<Rrd>::DisableActive()
{
    flag_is_active_=0;
}

template <class Rrd>
void Table<Rrd>::GetMeanAndVariance(int item_index)
{    
    if(1==item_index){
        if(flag_is_active_==0){
            double sum=0;
            // calculate average
            for(int i=0;i<records_.size();i++){
                sum+=records_[i].item1_;
            }
            mean_=sum/records_.size();

            // calculate variance
            sum=0;
            for(int i=0;i<records_.size();i++){
                sum+=pow(records_[i].item1_-mean_,2);
            }
            sigma_=sqrt(sum/(records_.size()-1));
        }
        else if(flag_is_active_==1){
            double sum=0;
            // calculate average
            for(int i=0;i<unbelievable_idx_.size();i++){
                sum+=records_[unbelievable_idx_[i]].item1_;
            }
            mean_=sum/unbelievable_idx_.size();

            // calculate variance
            sum=0;
            for(int i=0;i<unbelievable_idx_.size();i++){
                sum+=pow(records_[unbelievable_idx_[i]].item1_-mean_,2);
            }
            sigma_=sqrt(sum/(unbelievable_idx_.size()-1));
        }
        else 
            cout<<"Error 01: active mode error!"<<endl;
    }
}

template <class Rrd>
void Table<Rrd>::GetMinimumAndMaximum(int item_index)
{
    if(1==item_index){
        if(flag_is_active_==0){
            min_=INT_MAX;
            max_=-INT_MAX;
            // #pragma omp parallel for
            for(int i=0;i<records_.size();i++){
                min_=min_<records_[i].item1_ ? min_:records_[i].item1_;
                max_=max_>records_[i].item1_ ? max_:records_[i].item1_;
            }
        }
        else if(flag_is_active_==1){
            min_=INT_MAX;
            max_=-INT_MAX;
            for(int i=0;i<unbelievable_idx_.size();i++){
                min_=min_<records_[unbelievable_idx_[i]].item1_ ? min_:records_[unbelievable_idx_[i]].item1_;
                max_=max_>records_[unbelievable_idx_[i]].item1_ ? max_:records_[unbelievable_idx_[i]].item1_;
            }
        }
        else
            cout<<"Error 02: GetMinimumAndMaximum error!"<<endl;
    }    
}

template<class Rrd>
double Table<Rrd>::GetMininum(int item_index)
{
    min_=INT_MAX;
    if(1==item_index){        
        // #pragma omp parallel for
        for(int i=0;i<records_.size();i++){
            min_=min_<records_[i].item1_ ? min_:records_[i].item1_;
        }
    }
    else if(2==item_index){
        for(int i=0;i<records_.size();i++){
            min_=min_<records_[i].item2_ ? min_:records_[i].item2_;
        }
    }

    return min_;
}

template<class Rrd>
double Table<Rrd>::GetMean(int item_index)
{
    double sum=0;
    // calculate average
    for(int i=0;i<records_.size();i++){
        sum+=records_[i].item1_;
    }
    mean_=sum/records_.size();
    return mean_;
}

template <class Rrd>
void Table<Rrd>::Standardize_Zscore(int item_index)
{
    if(1==item_index){
        if(flag_is_active_==0){
            GetMeanAndVariance();
            #pragma omp parallel for
            for(int i=0;i<records_.size();i++){
                records_[i].item1_=(records_[i].item1_-mean_)/sigma_;
            }
        }
        else if(flag_is_active_==1){
            GetMeanAndVariance();
            #pragma omp parallel for
            for(int i=0;i<unbelievable_idx_.size();i++){
                records_[unbelievable_idx_[i]].item1_=(records_[unbelievable_idx_[i]].item1_-mean_)/sigma_;
            }
        }
        else
            cout<<"Error 03: Standardize Zscore error!"<<endl;        
    }
}
template <class Rrd>
void Table<Rrd>::Normalize_Min_Max(int item_index)
{
    if(1==item_index){
        if(flag_is_active_==0){
            GetMinimumAndMaximum();
            #pragma omp parallel for
            for(int i=0;i<records_.size();i++)
                records_[i].item1_=(records_[i].item1_-min_)/(max_-min_);        
        }
        else if(flag_is_active_==1){
            GetMinimumAndMaximum();
            #pragma omp parallel for
            for(int i=0;i<unbelievable_idx_.size();i++)
                records_[unbelievable_idx_[i]].item1_=(records_[unbelievable_idx_[i]].item1_-min_)/(max_-min_); 
        }
        else
            cout<<"Error 04: Normalize min max error!"<<endl;
    }
}
template <class Rrd>
void Table<Rrd>::Normalize_Tanh(double scale,int item_index)
{
    if(1==item_index){        
        if(flag_is_active_==0){
            GetMeanAndVariance();

            #pragma omp parallel for
            for(int i=0;i<records_.size();i++){
                records_[i].item1_=tanh(scale*records_[i].item1_);
            }    
        }
        else if(flag_is_active_==1){
            GetMeanAndVariance();

            #pragma omp parallel for
            for(int i=0;i<unbelievable_idx_.size();i++){
                records_[unbelievable_idx_[i]].item1_=tanh(scale*records_[unbelievable_idx_[i]].item1_);
            }
        }
        else
           cout<<"Error 05: Normalize Tanh error!"<<endl; 
    }
}
template <class Rrd>
void Table<Rrd>::Ascending(int item_index)
{
    if(1==item_index)
        sort(records_.begin(),records_.end(),[](Rrd& e1,Rrd& e2){ return e1.item1_<e2.item1_;});
}
template <class Rrd>
void Table<Rrd>::Descending(int item_index)
{
    if(1==item_index)
        sort(records_.begin(),records_.end(),[](Rrd& e1, Rrd& e2){ return e1.item1_> e2.item1_;});
}
template <class Rrd>
void Table<Rrd>::Write(string path)
{
    ofstream fout(path);
    if(flag_is_active_==0){
        for(int i=0;i<records_.size();i++){
            fout<<records_[i].item1_<<endl;
        }
    }
    else if(flag_is_active_==1){
        GetBoxplot();
        for(int i=0;i<unbelievable_idx_.size();i++){
            fout<<records_[unbelievable_idx_[i]].item1_<<endl;
        }
    }
    fout.close();
}
template <class Rrd>
void Table<Rrd>::MultiplyQ6(int m1,int m2,double scale)
{
    GetMinimumAndMaximum();
    double step=(max_-min_)/6.0;
    double lower_bound=(m1-1)*step+min_;
    double upper_bound=m2*step+min_;
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        if(records_[i].item1_>= lower_bound && records_[i].item1_<= upper_bound){
            records_[i].item1_=records_[i].item1_*scale;
        }
    }
}
template <class Rrd>
void Table<Rrd>::Multiply(double val)
{
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++)
        records_[i].item1_*=val;
}
template <class Rrd>
void Table<Rrd>::Overturn()
{
    GetMinimumAndMaximum();
    // #pragma omp parallel for
    for(int i=0;i<records_.size();i++){        
        records_[i].item1_=max_-records_[i].item1_;        
    }
}
template <class Rrd>
void Table<Rrd>::SetQ6(int m1,int m2)
{
    GetMinimumAndMaximum();
    double step=(max_-min_)/6.0;
    double lower_bound=(m1-1)*step+min_;
    double upper_bound=m2*step+min_;
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        if(records_[i].item1_>= lower_bound && records_[i].item1_<= upper_bound){
            records_[i].item1_=1.0f;
        }
    }
}
template <class Rrd>
void Table<Rrd>::GetNormalDistributionError()
{
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++){
        records_[i].item1_=Erf(records_[i].item1_);
    }
}
template <class Rrd>
void Table<Rrd>::LocalFilter(string str, pcl::PointCloud<PointType>::Ptr cloud,int K)
{
    vector<double> lf;
    lf.resize(cloud->points.size());
    pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
    kdtree->setInputCloud(cloud);

    if("average"==str){
         // #pragma omp parallel for
        for(int i=0;i<cloud->points.size();i++){
            vector<int> idx(K+1);
            vector<float> dist(K+1);        
            kdtree->nearestKSearch(cloud->points[i], K+1, idx, dist);
            double sum=0;
            for(int j=1;j<K+1;j++){
                sum+=records_[idx[j]].item1_;
            }
            sum/=K;
            lf[i]=records_[i].item1_/sum-1.0f;
        }

        for(int i=0;i<cloud->points.size();i++){
            records_[i].item1_=lf[i];
        }
    }
    else if("std"==str){
        for(int i=0;i<cloud->points.size();i++){
            vector<int> idx(K+1);
            vector<float> dist(K+1);        
            kdtree->nearestKSearch(cloud->points[i], K+1, idx, dist);
            vector<double> tmp;
            for(int j=1;j<K+1;j++){
                tmp.push_back(abs(records_[i].item1_-records_[idx[j]].item1_));
            }
            lf[i]=VectorSum(tmp);
            // lf[i]=VectorStd(tmp);
        }

        for(int i=0;i<cloud->points.size();i++){
            records_[i].item1_=lf[i];
        }
    }
   
}
template <class Rrd>
void Table<Rrd>::GetHistogram(int k)
{
    records_hist_.resize(k);
    for(int i=0;i<records_.size();i++){
        records_hist_[records_[i].item1_]+=1;
    }
}
template <class Rrd>
void Table<Rrd>::GaussianOutlier(int k)
{
    GetMeanAndVariance();
    double lower_thresh=max((double)0.0,mean_-k*sigma_);
    double upper_thresh=mean_+2*sigma_;
    for(int i=0;i<records_.size();i++){
        if(records_[i].item1_<=lower_thresh){
            unbelievable_idx_.push_back(records_[i].id_);
        }

        if(records_[i].item1_>=upper_thresh){
            believable_idx_.push_back(records_[i].id_);
        }
    }
}
template <class Rrd>
double Table<Rrd>::Quantile(double p)
{   
    SortBackup(1);
    double Q_idx=1+(records_sorted_.size()-1)*p;
    int Q_idx_integer=(int)Q_idx;
    double Q_idx_decimal=Q_idx-Q_idx_integer;
    double Q=records_sorted_[Q_idx_integer-1].item1_+(records_sorted_[Q_idx_integer].item1_-records_sorted_[Q_idx_integer-1].item1_)*Q_idx_decimal;    
    return Q;
}
template <class Rrd>
double Table<Rrd>::Median()
{
    SortBackup(1);
    double median=0;
    if(records_sorted_.size()%2==1){// odd
        median= records_sorted_[(records_sorted_.size()+1.0)/2.0-1].item1_;
    }
    else{ // even        
        median=(records_sorted_[records_sorted_.size()/2-1].item1_+records_sorted_[records_sorted_.size()/2].item1_)/2.0f;
    }
    median_=median;
    return median;
}
template <class Rrd>
void Table<Rrd>::SortBackup(int item_index)
{
    // Step 1: Does it have backup ?
    if(records_sorted_.size()==0){
        records_sorted_.resize(records_.size());
        for(int i=0;i<records_.size();i++){
            records_sorted_[i].id_=records_[i].id_;
            records_sorted_[i].item1_=records_[i].item1_;
        }
    }

    // Step 2: Sort Backup Records
    if(item_index==1){
        sort(records_sorted_.begin(),records_sorted_.end(),[](Rrd& e1, Rrd& e2){ return e1.item1_<e2.item1_;});
    }
    else if(item_index==0){
        sort(records_sorted_.begin(),records_sorted_.end(),[](Rrd& e1, Rrd& e2){ return e1.id_<e2.id_;});
    }
}
template <class Rrd>
void Table<Rrd>::GetBoxplot(double lamda, int item_index)
{
    if(item_index==1){
        GetMinimumAndMaximum();
        GetMeanAndVariance();

        // Get Median
        median_=Median();
        // Q1
        Q1_=Quantile(0.25);
        // Q3
        Q3_=Quantile(0.75);
        // IQR
        double IQR=Q3_-Q1_;

        // Part 01: upper inner limit
        upper_inner_limit_=Q3_+lamda*IQR;
        int i=0;
        for(;i<records_sorted_.size();i++){
            if(records_sorted_[i].item1_>upper_inner_limit_){
                upper_inner_limit_=records_sorted_[i-1].item1_;
                break;
            }
        }
        if(i==records_sorted_.size())
            upper_inner_limit_=min(upper_inner_limit_,records_sorted_[i-1].item1_);

        // Part 02: lower inner limit
        lower_inner_limit_=Q1_-lamda*IQR;
        i=0;
        for(;i<records_sorted_.size();i++){
            if(records_sorted_[i].item1_>lower_inner_limit_){
                lower_inner_limit_=records_sorted_[i].item1_;
                break;
            }
        }

        if(i==0)
            lower_inner_limit_=max(lower_inner_limit_,records_sorted_[0].item1_);          
    }

    // Get all inactive index
    for(int i=0;i<records_.size();i++){
        if(records_[i].item1_< lower_inner_limit_){
            lower_unbelievable_idx_.push_back(i);
            unbelievable_idx_.push_back(i);
        }
        else if(records_[i].item1_ > upper_inner_limit_){
            upper_unbelievable_idx_.push_back(i);
            unbelievable_idx_.push_back(i);
        }
        else{
            believable_idx_.push_back(i);
        }
    }
}
template <class Rrd>
void Table<Rrd>::push_back(Rrd e)
{
    records_.push_back(e);
}
template <class Rrd>
void Table<Rrd>::Print()
{
    for(int i=0;i<records_.size();i++){
        cout<<records_[i].item1_<<" ";
    }
    cout<<endl;
}
template <class Rrd>
void Table<Rrd>::GetCorrespondingColor(int item_index)
{
    if(item_index==1){
        color_.resize(records_.size());

        if(flag_is_active_==1){
            // Get Minimum and Maximum
            GetMinimumAndMaximum();

            // Lower inactive value
            for(int i=0;i<lower_unbelievable_idx_.size();i++){
                int itmp=lower_unbelievable_idx_[i];
                color_[itmp].r=0;
                color_[itmp].g=0;
                color_[itmp].b=255;
            }

            // Upper inactive value
            for(int i=0;i<upper_unbelievable_idx_.size();i++){
                int itmp=upper_unbelievable_idx_[i];
                color_[itmp].r=255;
                color_[itmp].g=0;
                color_[itmp].b=0;
            }

            // 
            for(int i=0;i<unbelievable_idx_.size();i++){
                int itmp=unbelievable_idx_[i];
                V3 ctmp=get_color(min_,max_,records_[itmp].item1_);                  
                color_[itmp].r=ctmp.r;
                color_[itmp].g=ctmp.g;
                color_[itmp].b=ctmp.b;
            }
        }
        else if(flag_is_active_==0){
            GetMinimumAndMaximum();

            for(int i=0;i<records_.size();i++){
                V3 ctmp=get_color(min_,max_,records_[i].item1_);
                color_[i].r=ctmp.r;
                color_[i].g=ctmp.g;
                color_[i].b=ctmp.b;

                // if(records_[i].item1_<27){
                //     color_[i].r=255;
                //     color_[i].g=0;
                //     color_[i].b=0;
                // }
            }
        }
        else{
            cout<<"Error: flag status error!"<<endl;
        }
    }
}
template <class Rrd>
void Table<Rrd>::nPLOF(pcl::PointCloud<PointType>::Ptr cloud,int K)
{
    pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>());
    kdtree->setInputCloud(cloud);
    
    vector<double> plof;    
    plof.resize(cloud->points.size());
	// #pragma omp parallel for
	for (int i = 0; i < cloud->points.size(); i++){        
        vector<int> idx(K+1);
		vector<float> dist(K+1);
		kdtree->nearestKSearch (cloud->points[i], K+1, idx, dist);
        double sum = 0;
        for (int j = 1; j < K+1; j++)
          sum += records_[idx[j]].item1_;
        sum /= K;
        plof[i] = records_[i].item1_ / sum  - 1.0f;
    }

    for(int i=0;i<cloud->points.size();i++){
        records_[i].item1_=plof[i];
    }
}
// blue -> green -> yellow -> red
template <class Rrd>
void Table<Rrd>::SetActiveIndex(string erea_color, string status)
{
    GetMinimumAndMaximum();
    double n=4.0;
    double step=(max_-min_)/n;
    double Q1=min_+step;
    double Q2=min_+2*step;
    double Q3=min_+3*step;
    if(erea_color=="1000" && status == "believable"){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<=Q1){
                believable_idx_.push_back(i);
            }
            else{
                unbelievable_idx_.push_back(i);
            }
        }
    }
    else if(erea_color=="1100" && status == "believable"){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<=Q2){
                believable_idx_.push_back(i);
            }
            else{
                unbelievable_idx_.push_back(i);
            }
        }
    }
    else if(erea_color=="1110" && status == "believable"){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<=Q3){
                believable_idx_.push_back(i);
            }
            else{
                unbelievable_idx_.push_back(i);
            }
        }
    }
     
}
template <class Rrd>
double Table<Rrd>::PDF(double t,double n,double h)
{
    double rst=0;
    for(int i=0;i<records_.size();i++){
        double tmp=sqrt(2)/2.0*(t-records_[i].item1_)/h;
        rst+=Erf(tmp);
    }
    rst=rst/(2.0*n)+1/2.0;
    return rst;
}
template <class Rrd>
double Table<Rrd>::Mean()
{
    double rst=0;
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++) 
        rst+=records_[i].item1_;
    rst=rst/records_.size();
}
template <class Rrd>
double Table<Rrd>::Std()
{
    double rst=0;
    double vmean=Mean();
    #pragma omp parallel for
    for(int i=0;i<records_.size();i++) 
        rst+=pow(records_[i].item1_-vmean,2);
    return sqrt(rst/(records_.size()-1));
}
template <class Rrd>
double Table<Rrd>::ReversePDF(double P)
{
    // double IQR=Quantile(0.75)-Quantile(0.25);
    // vector<int> new_indices;
    // new_indices=GetIndex(">=",Quantile(0.75)+3.0*IQR);
    // FastRemove(new_indices);

    // update parameters
    double n=records_.size();
    double IQR=Quantile(0.75)-Quantile(0.25);
    double h=0.9*pow(n,-0.2)*min(Std(),IQR/1.34);

    double t0=Minimum();
    double t1=Mean();
    double t2=Maximum();
    double ptmp=PDF(t1,n,h);
    while(abs(ptmp-P)>0.001){
        if(ptmp>P)
            t2=t1;
        else
            t0=t1;
        t1=(t0+t2)/2.0;
        ptmp=PDF(t1,n,h);
    }
    return t1;
}
template <class Rrd>
void Table<Rrd>::Find(vector<int>& out, string str,double thresh1, double thresh2)
{
    if(">"==str){
       double t1=ReversePDF(thresh1);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>t1)
                out.push_back(i);
        }
    }
    else if("<"==str){
        double t1=ReversePDF(thresh1);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<t1)
                out.push_back(i);
        }
    }
    else if("t1<x<t2"==str){
        double t1=ReversePDF(thresh1);
        double t2=ReversePDF(thresh2);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>t1 && records_[i].item1_<t2){
                out.push_back(i);
            }
        }
    }
    else if("x<t1 || x>t2"==str){
        double t1=ReversePDF(thresh1);
        double t2=ReversePDF(thresh2);
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<t1 || records_[i].item1_>t2){
                out.push_back(i);
            }
        }
    }
}
template <class Rrd>
double Table<Rrd>::Minimum()
{
    double rst=INT_MAX;
    for(int i=0;i<records_.size();i++)
        rst=rst<records_[i].item1_ ? rst:records_[i].item1_;
    return rst;
}
template <class Rrd>
double Table<Rrd>::Maximum()
{
    double rst=-INT_MAX;
    for(int i=0;i<records_.size();i++)
        rst=rst>records_[i].item1_ ? rst:records_[i].item1_;
    return rst;
}
template <class Rrd>
double Table<Rrd>::WritePDF()
{
    double IQR=Quantile(0.75)-Quantile(0.25);
    vector<int> new_indices;
    double thresh=Quantile(0.75)+3.0*IQR;
    new_indices=GetIndex(">=",thresh);
    FastRemove(new_indices);
   

    // update parameters
    double n=records_.size();
    IQR=Quantile(0.75)-Quantile(0.25);
    double h=0.9*pow(n,-0.2)*min(Std(),IQR/1.34);
    double vmin=Minimum();
    double vmax=Maximum();
    double step=(vmax-vmin)/55.0;
    ofstream fout("Result/PDF.csv");
    double ptmp=0;
    for(double t=vmin;t<=vmax;t+=step){
        ptmp=PDF(t,n,h);
        fout<<ptmp<<endl;
    }
    fout.close();
    cout<<ptmp<<endl;
}
template <class Rrd>
vector<int> Table<Rrd>::GetIndex(string str, double thresh)
{
    vector<int> rst;
    if(str==">"){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>thresh)
                rst.push_back(i);
        }
    }
    else if(str==">="){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_>=thresh)
                rst.push_back(i);
        }
    }
    else if(str=="<"){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<thresh)
                rst.push_back(i);
        }
    }
    else if(str=="<="){
        for(int i=0;i<records_.size();i++){
            if(records_[i].item1_<=thresh)
                rst.push_back(i);
        }
    }
    return rst;
}
template <class Rrd>
void Table<Rrd>::FastRemove(int index)
{
    int tail_idx=records_.size()-1;
    records_[index].id_=records_[tail_idx].id_;
    records_[index].item1_=records_[tail_idx].item1_;
    records_.pop_back();
}
template <class Rrd>
void Table<Rrd>::FastRemove(vector<int> indices)
{
    for(int i=0;i<indices.size();i++){
        FastRemove(indices[i]);
    }
}