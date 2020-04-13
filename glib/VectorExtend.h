#pragma once
#include <iostream>
#include <vector>
#include <limits.h>
#include <algorithm>
using namespace std;
double VectorMean(vector<double>& dat);
double VectorMaximum(vector<double>& dat);
double VectorMinimum(vector<double>& dat);
double VectorStd(vector<double>& dat);
double VectorQuantile(vector<double>& dat,double p);
double VectorSum(vector<double>& dat);
void VectorDelete(vector<int>&rdat,int index);
void VectorDelete(vector<int>&raw_dat, vector<int>& delete_dat);
int VectorIQR(vector<double>& dat);

/*
    Logical Operation
*/
void VectorDifference(vector<int>& dat1,vector<int>& dat2,vector<int>& out);
void VectorUnion(vector<int>& dat1,vector<int>& dat2,vector<int>& out);
void VectorIntersection(vector<int>& dat1,vector<int>& dat2,vector<int>& out);