#ifndef PWA_PARAS_H
#define PWA_PARAS_H

// 定义数据结构PWA_PARAS，其中包括所有需要用于计算分波振幅，但在拟合过程中不需要改变的参数
typedef struct {
    Double_t wu[4],w0p22[4],ak23w[4],w2p2[4],w2p1[4],w2p3[2],w2p4[2],w2p5[2];
    Double_t b2qf2xx,b4qjvf2,b2qjv2,b2qjv3,b2qbv2,b2qbv3,b1qjv2,b1qjv3,b1qbv2,b1qbv3,b1q2r23;
    Double_t sv,s23,sv2,sv3;
    Double_t wpf22[2];
    Double_t b2qjvf2;
    Double_t w1p12_1[4],w1p13_1[4],w1p12_2[4],w1p13_2[4];
    Double_t w1p12_3[2],w1p13_3[2],w1p12_4[2],w1p13_4[2];
    Double_t w1m12[2],w1m13[2];
} PWA_PARAS;


#endif

