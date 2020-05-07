#include "VectorExtend.h"
double VectorMean(vector<double>& dat)
{
    double rst=0;
    for(int i=0;i<dat.size();i++)
        rst+=dat[i];
    
    return rst/dat.size();
}

double VectorMinimum(vector<double>& dat)
{
    auto vmin=std::min_element(dat.begin(),dat.end());
    return *vmin;
}
double VectorMaximum(vector<double>& dat)
{
    auto vmax=std::max_element(dat.begin(),dat.end());
    return *vmax;
}
double VectorStd(vector<double>& dat)
{
    double rst=0;
    double dat_mean=VectorMean(dat);
    for(int i=0;i<dat.size();i++){
        rst+=(dat[i]-dat_mean)*(dat[i]-dat_mean);
    }
    return rst/(dat.size()-1);
}

double VectorQuantile(vector<double>& dat,double p)
{
    sort(dat.begin(),dat.end());
    double Q_idx=1+(dat.size()-1)*p;
    int Q_idx_integer=(int)Q_idx;
    double Q_idx_decimal=Q_idx-Q_idx_integer;
    double Q=dat[Q_idx_integer-1]+(dat[Q_idx_integer]-dat[Q_idx_integer-1])*Q_idx_decimal;    
    return Q;
}

int VectorIQR(vector<double>& dat)
{
    double IQR=VectorQuantile(dat,0.75)-VectorQuantile(dat,0.25);
    double thresh=VectorQuantile(dat,0.75)+10*IQR;
    for(int i=0;i<dat.size();i++)
    {
        if(dat[i]>thresh){
            return 1;
        }
    }
    return 0;
}

double VectorSum(vector<double>& dat)
{
    double sum=0;
    for(int i=0;i<dat.size();i++)
    {
        sum+=dat[i];
    }
    return sum;
}
void VectorDelete(vector<int>&rdat,int index)
{
    rdat[index]=rdat[rdat.size()-1];
    rdat.pop_back();
}
void VectorDelete(vector<int>&raw_dat, vector<int>& delete_dat)
{
    if(delete_dat.size()!=0){
        vector<int> target;
        sort(raw_dat.begin(),raw_dat.end());
        sort(delete_dat.begin(),delete_dat.end());
        set_difference(raw_dat.begin(),raw_dat.end(),delete_dat.begin(),delete_dat.end(),std::back_inserter(target));
        raw_dat.clear();
        for(int i=0;i<target.size();i++) raw_dat.push_back(target[i]);
    }
}
void VectorDifference(vector<int>& dat1,vector<int>& dat2,vector<int>& out)
{
    sort(dat1.begin(),dat1.end());
    sort(dat2.begin(),dat2.end());
    set_difference(dat1.begin(),dat1.end(),dat2.begin(),dat2.end(),std::back_inserter(out));
}
void VectorUnion(vector<int>& dat1,vector<int>& dat2,vector<int>& out)
{
    sort(dat1.begin(),dat1.end());
    sort(dat2.begin(),dat2.end());
    set_union(dat1.begin(),dat1.end(),dat2.begin(),dat2.end(),std::back_inserter(out));
}
void VectorIntersection(vector<int>& dat1,vector<int>& dat2,vector<int>& out)
{
    sort(dat1.begin(),dat1.end());
    sort(dat2.begin(),dat2.end());
    set_intersection(dat1.begin(),dat1.end(),dat2.begin(),dat2.end(),std::back_inserter(out));
}