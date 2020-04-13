#include "SignalProcessing.h"
void DaubechiesWavelet(vector<double>& dat,vector<double>& output)
{
    double h0=(1+sqrt(3))/(4*sqrt(2));
    double h1=(3+sqrt(3))/(4*sqrt(2));
    double h2=(1+sqrt(3))/(4*sqrt(2));
    double h3=(1+sqrt(3))/(4*sqrt(2));
    double g0=h3;
    double g1=-h2;
    double g2=h1;
    double g3=-h0;
    output.resize(dat.size()/2);
    for(int i=0;i<output.size();i++){
        output[i]=g0*dat[2*i]+g1*dat[2*i+1]+g2*dat[2*i+2]+g3*dat[2*i+3];
    }
}

void DaubechiesWavelet(vector<float>& dat,vector<double>& output)
{
    double h0=(1+sqrt(3))/(4*sqrt(2));
    double h1=(3+sqrt(3))/(4*sqrt(2));
    double h2=(1+sqrt(3))/(4*sqrt(2));
    double h3=(1+sqrt(3))/(4*sqrt(2));
    double g0=h3;
    double g1=-h2;
    double g2=h1;
    double g3=-h0;
    output.resize(dat.size()/2);
    for(int i=0;i<output.size();i++){
        output[i]=g0*dat[2*i]+g1*dat[2*i+1]+g2*dat[2*i+2]+g3*dat[2*i+3];
    }
}