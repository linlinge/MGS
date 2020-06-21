/*
    Name:   Outliers Level Evaluation Module
    Author: linlinge
    Date:   2020.04.02
*/
#pragma once
#include "Table.h"
#include "PCLExtend.h"
#include <stdlib.h>
#include<time.h>
class StyleEvaluation
{
    public:
		pcl::PointCloud<PointType>::Ptr cloud_;
		pcl::search::KdTree<PointType>::Ptr kdtree_;
		Table<Rrd1> rst_dmean;
		Table<Rrd1> rst_meval;
		Table<Rrd1> rst_db2;
		double dnst_,dmean_,meval_,db2_;
		double db2_Q3_,db2_IQR_;
		double meval_Q3_,meval_IQR_;
		double og_;
		int n_,N_;
		/* Style Evaluation Metrics */
		void OutlierGradeMetric(string path, int K=30);
		void UniformityMetric(string path, int K=30);
		void SigularityMetric(string path, int K=30);

        double ApplyEigenvalue(int K=32);
		double GetMinorEval(int K=32, string str="Common");
		double GetDmean(int K=32,string str="Common");
		double GetDB2(int K=100, string str="Common");
		StyleEvaluation(pcl::PointCloud<PointType>::Ptr cloud,pcl::search::KdTree<PointType>::Ptr kdtree=NULL);
};