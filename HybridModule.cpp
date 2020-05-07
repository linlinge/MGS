#include "HybridModule.h"
void H()
{
  int max_status=*std::max_element(status_.begin(),status_.end());
  if(max_status==2){
    #pragma omp parallel for
    for(int i=0;i<status_.size();i++){
      if(status_[i]>0)
        status_[i]=1;
      else
        status_[i]=0;
    } 
    cout<<"wokaka"<<endl;
  }
  else if(max_status==3){
    #pragma omp parallel for
    for(int i=0;i<status_.size();i++){
      if(status_[i]==3)
        status_[i]=1;
      else
        status_[i]=0;
    } 
    cout<<"wokaka2"<<endl;
  }
  else 
    cout<<"H status error!"<<endl;
}

void U()
{
  #pragma omp parallel for
  for(int i=0;i<status_.size();i++){
    if(status_[i]==3)
      status_[i]=1;
    else
      status_[i]=0;
  } 
}
void D()
{
  #pragma omp parallel for
  for(int i=0;i<status_.size();i++){
    if(status_[i]==1)
      status_[i]=1;
    else
      status_[i]=0;
  } 
}
void J()
{
  #pragma omp parallel for
  for(int i=0;i<status_.size();i++){
    if(status_[i]>0)
      status_[i]=1;
    else
      status_[i]=0;
  }
}
// void H(void (*f1)(int K,double P,int active_layer), int K1, double P1,int active_layer1,
//        void (*f2)(int K,double P,int active_layer), int K2, double P2,int active_layer2)
// {
//     f1(K1,P1,active_layer1);
//     StatusToBinary();
//     f2(K2,P2,active_layer2);
//     StatusToBinary();
// }
// void U(void (*f1)(int K,double P,int active_layer), int K1, double P1,int active_layer1,
//        void (*f2)(int K,double P,int active_layer), int K2, double P2,int active_layer2)
// {
//     f1(K1,P1,active_layer1);
//     StatusToBinary();
//     f2(K2,P2,active_layer2);
//     StatusToBinary();
// }
// void D(void (*f1)(int K,double P,int active_layer), int K1, double P1,int active_layer1,
//        void (*f2)(int K,double P,int active_layer), int K2, double P2,int active_layer2)
// {
// //     f1(K1,P1,active_layer1);
// //     StatusToBinary();
// //     f2(K2,P2,active_layer2);
// //     StatusToBinary();
// }
// void J(void (*f1)(int K,double P,int active_layer), int K1, double P1,int active_layer1,
//        void (*f2)(int K,double P,int active_layer), int K2, double P2,int active_layer2)
// {
//     f1(K1,P1,active_layer1);
//     StatusToBinary();
//     f2(K2,P2,active_layer2);
//     StatusToBinary();
// }