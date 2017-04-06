/*************************************************************************
	> File Name: kernel_calEva.h
	> Author:
	> Mail:
	> Created Time: 2017年03月21日 星期二 18时20分53秒
 ************************************************************************/

#ifndef _KERNEL_CALEVA_H
#define _KERNEL_CALEVA_H
#include <vector>
//double *d_float_pp=NULL;
//#include "cu_PWA_PARAS.h"
    //int initialize_data(std::vector<PWA_PARAS>&, DataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //int data_distribution(DataPointers&, CudaDataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //double calEva(const PWA_PARAS &pp, int idp);
    //double kernel_calEva(const PWA_PARAS &pp,int idp);
    //void func(DataPointers& cpu_data_pointers);
    void cu_malloc_h_pp(double *,double *&,int,int);
    int host_store_fx(std::vector<double *>d_float_pp,int *h_parameter,double *h_paraList,int para_size, double *h_fx,double * h_mlk,int numElements,int begin);
#endif
