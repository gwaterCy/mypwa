#ifndef CU_PWA_PARAS_H
#define CU_PWA_PARAS_H
//定义数据结构PWA_PARAS，其中包括所有需要用于计算分波振幅，但在拟合过程中不需要改变的参数

typedef struct {
    double wu[4],w0p22[4],ak23w[4],w2p2[4],w2p1[4],w2p3[2],w2p4[2],w2p5[2];
    double b2qf2xx,b4qjvf2,b2qjv2,b2qjv3,b2qbv2,b2qbv3,b1qjv2,b1qjv3,b1qbv2,b1qbv3,b1q2r23;
    double sv,s23,sv2,sv3;
    double wpf22[2];
    double b2qjvf2;
    double w1p12_1[4],w1p13_1[4],w1p12_2[4],w1p13_2[4];
    double w1p12_3[2],w1p13_3[2],w1p12_4[2],w1p13_4[2];
    double w1m12[2],w1m13[2];
} cu_PWA_PARAS;

#endif

