#pragma once
#include<fstream>
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
void VectorNormalize(vector<double>& dat);

/*
    Logical Operation
*/
void VectorDifference(vector<int>& dat1,vector<int>& dat2,vector<int>& out);
void VectorUnion(vector<int>& dat1,vector<int>& dat2,vector<int>& out);
void VectorIntersection(vector<int>& dat1,vector<int>& dat2,vector<int>& out);

/* Find */
int VecFindPos(vector<int>& buf,int dat);
/* Get the indices of all elements of buf2 in buf1*/
void VecFindPos(vector<int>& buf1,vector<int>& buf2, vector<int>& out);
void VecPrint(vector<int>& dat);
void VecUnique(vector<int>& dat);
/* write */
void VectorWrite(string filename,vector<double>& dat,string mode="append");
void VectorWrite(string filename,vector<float>& dat,string mode="append");