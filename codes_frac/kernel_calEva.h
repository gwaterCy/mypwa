/*************************************************************************
	> File Name: kernel_calEva.h
	> Author: 
	> Mail: 
	> Created Time: 2017年03月21日 星期二 18时20分53秒
 ************************************************************************/

#ifndef _KERNEL_CALEVA_H
#define _KERNEL_CALEVA_H
//float *d_float_pp=NULL;
//#include "cu_PWA_PARAS.h"
    //int initialize_data(std::vector<PWA_PARAS>&, DataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //int data_distribution(DataPointers&, CudaDataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //float calEva(const PWA_PARAS &pp, int idp);
    //float kernel_calEva(const PWA_PARAS &pp,int idp);
    //void func(DataPointers& cpu_data_pointers);
    void cu_malloc_h_pp(float *,float *&,int);
    int host_store_fx(float *d_float_pp,int *h_parameter,float *h_paraList,int para_size, float *h_fx,float * h_mlk,int numElements,int begin);
#endif
