#include <cuda_runtime.h>
#include "cuComplex.h"
#include <iostream>
#include "cu_PWA_PARAS.h"
#include <vector>
#include <fstream>
#include <math.h>
#include "cu_DPFPropogator.h"
#include "kernel_calEva.h"
#include <assert.h>

using namespace std;

#define CUDA_CALL(x) {const cudaError_t a=(x); if(a != cudaSuccess) {printf("\nerror in line:%d CUDAError:%s(err_num=%d)\n",__LINE__,cudaGetErrorString(a),a); cudaDeviceReset(); assert(0); }}
 


 __device__ double calEva(const cu_PWA_PARAS *pp, const int * parameter , const double2 * complex_para ,const double * d_paraList,double *d_mlk,int idp) 
    ////return square of complex amplitude
{
    //	static int A=0;
    //	A++;
    
    int _N_spinList     =parameter[0];
    int _N_massList     =parameter[1];
    int _N_mass2List    =parameter[2];
    int _N_widthList    =parameter[3];
    int _N_g1List       =parameter[4];
    int _N_g2List       =parameter[5];
    int _N_b1List       =parameter[6];
    int _N_b2List       =parameter[7];
    int _N_b3List       =parameter[8];
    int _N_b4List       =parameter[9];
    int _N_b5List       =parameter[10];
    int _N_rhoList      =parameter[11];
    int _N_fracList     =parameter[12];
    int _N_phiList      =parameter[13];
    int _N_propList     =parameter[14];
    const int const_nAmps=parameter[15];
    double value = 0.;
    //double2 fCF[const_nAmps][4];
    double2 *fCF=complex_para; 
    //double2 (*fCF)[4]=(double2 (*)[4])malloc(sizeof(double2)*const_nAmps*4);
    //double2 fCP[const_nAmps];
    //double2 * fCP=(double2 *)malloc(sizeof(double2)*const_nAmps);
    double2 * fCP=&complex_para[4*const_nAmps];
    double2 * crp1=&complex_para[5*const_nAmps];
    double2 * crp11=&complex_para[6*const_nAmps];


    //double2 pa[const_nAmps][const_nAmps];
    double2 * pa=&complex_para[7*const_nAmps];
    double2 * fu=&complex_para[(7+const_nAmps)*const_nAmps];


    /*double2 **pa,**fu;
    pa=(double2 **)malloc(sizeof(double2 *)*const_nAmps);
    fu=(double2 **)malloc(sizeof(double2 *)*const_nAmps);
    for(int i=0;i<const_nAmps;i++)
    {
        pa[i]=(double2 *)malloc(sizeof(double2)*const_nAmps);
        fu[i]=(double2 *)malloc(sizeof(double2)*const_nAmps);
    }
    //double2 fu[const_nAmps][const_nAmps];
    //double2 crp1[const_nAmps];
    double2 * crp1=(double2 *)malloc(sizeof(double2)*const_nAmps);
    //double2 crp11[const_nAmps];
    double2 * crp11=(double2 *)malloc(sizeof(double2)*const_nAmps);
    */
    double2 cr0p11;
    //double2 ca2p1;
    double2 cw2p11;
    double2 cw2p12;
    double2 cw2p15;
    double2 cw;
    double2 c1p12_12,c1p13_12,c1p12_13,c1p13_13,c1p12_14,c1p13_14;
    double2 cr1m12_1,cr1m13_1;
    double2 crpf1,crpf2;

    for(int index=0; index<const_nAmps; index++) {
        double rho0 = d_paraList[_N_rhoList++];
        double frac0 = d_paraList[_N_fracList++];
        double phi0 = d_paraList[_N_phiList++];
        int spin_now = d_paraList[_N_spinList++];
        int propType_now = d_paraList[_N_propList++];
    //cout<<"haha: "<< __LINE__ << endl;

        rho0 *= std::exp(frac0);
        fCP[index]=make_cuDoubleComplex(rho0*std::cos(phi0),rho0*std::sin(phi0));
        //        //cout<<"fCP[index]="<<fCP[index]<<endl;
        //std::cout << __FILE__ << __LINE__ << " : " << propType_now << std::endl;
        switch(propType_now)
        {
            //  //cout<<"haha: "<< __LINE__ << endl;
            //                     ordinary  Propagator  Contribution
            case 1:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    double mass0 = d_paraList[_N_massList++];
                    double width0 = d_paraList[_N_widthList++];
                    //					//cout<<"mass0="<<mass0<<endl;
                    //					//cout<<"width0="<<width0<<endl;
                    crp1[index]=propogator(mass0,width0,pp->s23);
                }
                break;
            //	Flatte   Propagator Contribution
            case 2:
                {
                    //RooRealVar *g1 = (RooRealVar*)_g1IterV[omp_id]->Next();
                    //RooRealVar *g2 = (RooRealVar*)_g2IterV[omp_id]->Next();
                    double mass980 = d_paraList[_N_massList++];
                    double g10 = d_paraList[_N_g1List++];
                    double g20 = d_paraList[_N_g2List++];
                    //double g10=g1->getVal();
                    //double g20=g2->getVal();
     //               			//cout<<"mass980="<<mass980<<endl;
     //               			//cout<<"g10="<<g10<<endl;
     //               			//cout<<"g20="<<g20<<endl;
     //                           //cout<<"pp.s23="<<pp.s23<< endl;
                    crp1[index]=propogator980(mass980,g10,g20,pp->s23);
     //               			//cout<<"crp1[index]="<<crp1[index]<<endl;
                }
                break;
                // sigma  Propagator Contribution
            case 3:
                {
                    //RooRealVar *b1 = (RooRealVar*)_b1IterV[omp_id]->Next();
                    //RooRealVar *b2 = (RooRealVar*)_b2IterV[omp_id]->Next();
                    //RooRealVar *b3 = (RooRealVar*)_b3IterV[omp_id]->Next();
                    //RooRealVar *b4 = (RooRealVar*)_b4IterV[omp_id]->Next();
                    //RooRealVar *b5 = (RooRealVar*)_b5IterV[omp_id]->Next();
                    //double mass600=mass->getVal();
                    //double b10=b1->getVal();
                    //double b20=b2->getVal();
                    //double b30=b3->getVal();
                    //double b40=b4->getVal();
                    //double b50=b5->getVal();
                    double mass600 = d_paraList[_N_massList++];
                    double b10 = d_paraList[_N_b1List++];
                    double b20 = d_paraList[_N_b2List++];
                    double b30 = d_paraList[_N_b3List++];
                    double b40 = d_paraList[_N_b4List++];
                    double b50 = d_paraList[_N_b5List++];
                    crp1[index]=propogator600(mass600,b10,b20,b30,b40,b50,pp->s23);
                    //			//cout<<"crp1[index]3="<<crp1[index]<<endl;
                }
                break;
                // 1- or 1+  Contribution
            case 4:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //double mass0=mass->getVal();
                    //double width0=width->getVal();
                    double mass0 = d_paraList[_N_massList++];
                    double width0 = d_paraList[_N_widthList++];
                    crp1[index]=propogator(mass0,width0,pp->sv2);
                    crp11[index]=propogator(mass0,width0,pp->sv3);
                }
                break;
                //  phi(1650) f0(980) include flatte and ordinary Propagator joint Contribution
            case 5:
                {
                    //RooRealVar *mass2  = (RooRealVar*)_mass2IterV[omp_id]->Next();
                    //RooRealVar *g1 = (RooRealVar*)_g1IterV[omp_id]->Next();
                    //RooRealVar *g2 = (RooRealVar*)_g2IterV[omp_id]->Next();
                    //double mass980=mass2->getVal();
                    //double g10=g1->getVal();
                    //double g20=g2->getVal();
                    double mass980 = d_paraList[_N_mass2List++];
                    double g10 = d_paraList[_N_g1List++];
                    double g20 = d_paraList[_N_g2List++];
                    //					//cout<<"mass980="<<mass980<<endl;
                    //					//cout<<"g10="<<g10<<endl;
                    //					//cout<<"g20="<<g20<<endl;
                    crp1[index]=propogator980(mass980,g10,g20,pp->sv);
                    //					//cout<<"crp1[index]="<<crp1[index]<<endl;
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //double mass1680=mass->getVal();
                    //double width1680=width->getVal();
                    double mass1680 = d_paraList[_N_massList++];
                    double width1680 = d_paraList[_N_widthList++];
                    //					//cout<<"mass1680="<<mass1680<<endl;
                    //					//cout<<"width1680="<<width1680<<endl;
                    crp11[index]=propogator(mass1680,width1680,pp->s23);
                    //					//cout<<"crp11[index]="<<crp11[index]<<endl;
                }
                break;
            case 6:
                {
                    //RooRealVar *width = (RooRealVar*)_widthIterV[omp_id]->Next();
                    //double mass0=mass->getVal();
                    //double width0=width->getVal();
                    double mass0 = d_paraList[_N_massList++];
                    double width0 = d_paraList[_N_widthList++];
                    //					//cout<<"mass0="<<mass0<<endl;
                    //					//cout<<"width0="<<width0<<endl;
                    crp1[index]=propogator1270(mass0,width0,pp->s23);
                    //			//cout<<"crp1[index]6="<<crp1[index]<<endl;
                }
            default :
                ;
        }
    //cout << "LINE: " << __LINE__ << endl;
        for(int i=0;i<2;i++){
            ////cout<<"haha: "<< __LINE__ << endl;
            //		//cout<<"spin_now="<<spin_now<<endl;
            switch(spin_now)
            {
                case 11:
                    //1+_1 contribution
                    //fCF[index][i]=pp.w1p12_1[i]*crp1[index]+pp.w1p13_1[i]*crp11[i];
                    fCF[index*4+i]=cuCadd( cuCmuldc(pp->w1p12_1[i],crp1[index]),cuCmuldc(pp->w1p13_1[i],crp11[i]) );

                    break;
                case 12:
                    //1+_2 contribution
                    //c1p12_12=crp1[index]/pp.b2qbv2;
                    c1p12_12=cuCdivcd(crp1[index],pp->b2qbv2);
                    //c1p13_12=crp11[index]/pp.b2qbv3;
                    c1p13_12=cuCdivcd(crp11[index],pp->b2qbv3);
                    //fCF[index][i]=pp.w1p12_2[i]*c1p12_12+pp.w1p13_2[i]*c1p13_12;
                    fCF[index*4+i]=cuCadd( cuCmuldc(pp->w1p12_2[i],c1p12_12) , cuCmuldc(pp->w1p13_2[i],c1p13_12) );
                
                    break;
                case 13:
                    //1+_3 contribution
                    //c1p12_13=crp1[index]/pp.b2qjv2;
                    c1p12_13=cuCdivcd(crp1[index],pp->b2qjv2);
                    //c1p13_13=crp11[index]/pp.b2qjv3;
                    c1p13_13=cuCdivcd(crp11[index],pp->b2qjv3);
                    //fCF[index][i]=pp.w1p12_3[i]*c1p12_13+pp.w1p13_3[i]*c1p13_13;
                    fCF[index*4+i]=cuCadd( cuCmuldc(pp->w1p12_3[i],c1p12_13) , cuCmuldc(pp->w1p13_3[i],c1p13_13) );

                    break;
                case 14:
                    //1+_4 contribution
                    //c1p12_12=crp1[index]/pp.b2qbv2;
                    c1p12_12=cuCdivcd(crp1[index],pp->b2qbv2);
                    
                    c1p13_12=cuCdivcd(crp11[index],pp->b2qbv3);
                    c1p12_14=cuCdivcd(c1p12_12,pp->b2qjv2);
                    c1p13_14=cuCdivcd(c1p13_12,pp->b2qjv3);
                    fCF[index*4+i]=cuCadd( cuCmuldc(pp->w1p12_4[i],c1p12_14), cuCmuldc(pp->w1p13_4[i],c1p13_14));

                    break;
                case 111:
                    //1-__1 contribution
                    cr1m12_1=cuCdivcd( cuCdivcd(crp1[index],pp->b1qjv2) , pp->b1qbv2);
                    cr1m13_1=cuCdivcd( cuCdivcd(crp11[index],pp->b1qjv3) , pp->b1qbv3);
                    fCF[index*4+i]=cuCadd( cuCmuldc(pp->w1m12[i],cr1m12_1), cuCmuldc(pp->w1m13[i],cr1m13_1));

                    break;
                case 191:
                    //phi(1650)f0(980)_1 contribution
                    //		//cout<<"b1q2r23="<<b1q2r23<<endl;
                    crpf1=cuCdivcd( cuCmul(crp1[index],crp11[index]),pp->b1q2r23 );
                    //		//cout<<"crpf1="<<crpf1<<endl;
                    fCF[index*4+i]=cuCmuldc(pp->ak23w[i],crpf1);
                    //	//cout<<"fCF[index][i]="<<fCF[index][i]<<endl;

                    break;
                case 192:
                    //phi(1650)f0(980)_2 contribution
                    crpf1=cuCdivcd( cuCmul(crp1[index],crp11[index]) , pp->b1q2r23);
                    crpf2=cuCdivcd(crpf1,pp->b2qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->wpf22[i],crpf2);

                    break;
                case 1:
                    //  //cout<<"haha: "<< __LINE__ << endl;
                    //01 contribution
                    //	//cout<<"wu[i]="<<wu[i]<<endl;
                    //	//cout<<"crp1[index]="<<crp1[index]<<endl;
                    //	//cout<<"index="<<index<<endl;
                    fCF[index*4+i]=cuCmuldc(pp->wu[i],crp1[index]);
                    //	//cout<<"fCF[index][i]="<<fCF[index][i]<<endl;
                    //	//cout<<"i="<<i<<endl;

                    break;
                case 2:
                    //02 contribution
                    cr0p11=cuCdivcd(crp1[index],pp->b2qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->w0p22[i],cr0p11);
                    //	//cout<<"fCF[index][i]02="<<fCF[index][i]<<endl;

                    break;
                case 21:
                    //21 contribution
                    //	//cout<<"b2qf2xx="<<b2qf2xx<<endl;
                    cw2p11=cuCdivcd(crp1[index],pp->b2qf2xx);
                    //	//cout<<"cw2p11="<<cw2p11<<endl;
                    //	//cout<<"w2p1[0]="<<w2p1[0]<<endl;
                    //	//cout<<"w2p1[1]="<<w2p1[1]<<endl;
                    fCF[index*4+i]=cuCmuldc(pp->w2p1[i],cw2p11);
                    //	//cout<<"fCF[index][i]21="<<fCF[index][i]<<endl;

                    break;
                case 22:
                    //22 contribution
                    cw2p11=cuCdivcd(crp1[index],pp->b2qf2xx);
                    cw2p12=cuCdivcd(cw2p11,pp->b2qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->w2p2[i],cw2p12);

                    break;
                case 23:
                    //23 contribution
                    cw2p11=cuCdivcd(crp1[index],pp->b2qf2xx);
                    cw2p12=cuCdivcd(cw2p11,pp->b2qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->w2p3[i],cw2p12);

                    break;
                case 24:
                    //24 contribution
                    cw2p11=cuCdivcd(crp1[index],pp->b2qf2xx);
                    cw2p12=cuCdivcd(cw2p11,pp->b2qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->w2p4[i],cw2p12);

                    break;
                case 25:
                    //25 contribution
                    cw2p11=cuCdivcd(crp1[index],pp->b2qf2xx);
                    cw2p15=cuCdivcd(cw2p11,pp->b4qjvf2);
                    fCF[index*4+i]=cuCmuldc(pp->w2p5[i],cw2p15);

                default:		;
            }
        }

    }
    double carry(0);
    //#pragmaint  omp parallel for reduction(+:value)
    for(int i=0;i<const_nAmps;i++){
        //  //cout<<"haha: "<< __LINE__ << endl;    int mlk_cro_size=sizeof(double)*numElements
        for(int j=0;j<const_nAmps;j++){
            cw=cuCmul(fCP[i],cuConj(fCP[j]));
            //    //cout<<"cw="<<cw<<endl;
            if(i==j) pa[i*const_nAmps+j]=make_cuDoubleComplex(cuCreal(cw),0.0);
            else if(i<j) pa[i*const_nAmps+j]=make_cuDoubleComplex(2*cuCreal(cw),0.0);
            else pa[i*const_nAmps+j]=make_cuDoubleComplex(0.0,2*cuCimag(cw));
            cw=make_cuDoubleComplex(0.0,0.0);
            for(int k=0;k<2;k++){
                cw=cuCadd(cw,cuCdivcd( cuCmul( fCF[i][k],cuConj(fCF[j][k]) ),(double)2.0) );
                //   //cout<<"cwfu="<<cw<<endl;

            }
            if(i<=j) fu[i*const_nAmps+j]=make_cuDoubleComplex(cuCreal(cw),0.0);
            if(i>j) fu[i*const_nAmps+j]=make_cuDoubleComplex(0.0,-cuCimag(cw));
            //      //cout<<"pa[i][j]="<<pa[i][j]<<endl;
            //      //cout<<"fu[i][j]="<<fu[i][j]<<endl;
            double temp = cuCreal( cuCmul(pa[i*const_nAmps+j],fu[i*const_nAmps+j]) );//i have a big change here 
            double y = temp - carry;
            double t = value + y;
            carry = (t - value) - y;

            value = t; // Kahan Summation
        }
    }

    for(int i=0;i<const_nAmps;i++){
        double2 cw=cuCmul(fCP[i],cuConj(fCP[i]));
        double pa=cuCreal(cw);

        cw=make_cuDoubleComplex(0.0,0.0);
        for(int k=0;k<2;k++){
            //cw+=fCF[i][k]*cuConj(fCF[i][k])/(double)2.0;
            cw=cuCadd(cw,cuCdivcd( cuCmul( fCF[i][k],cuConj(fCF[i][k]) ),(double)2.0) );
        }
        double fu=cuCreal(cw);
        d_mlk[idp*const_nAmps+i] = pa * fu;
    }
    /*
    free(fCP);
    for(int i=0;i<const_nAmps;i++)
    {
        free(pa[i]);
        free(fu[i]);
        //free(fCF[i]);
    } 
    free(fCF);
    free(pa);
    free(fu);
    free(crp1);
    free(crp11);
*/
    return (value <= 0) ? 1e-20 : value;
}

__global__ void kernel_store_fx(const double * float_pp,const int *parameter,const double *d_paraList,double * d_fx,double *d_mlk,int numElements,int begin)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i<numElements && i>= begin)
    {
        int pwa_paras_size = sizeof(cu_PWA_PARAS) / sizeof(double);
        cu_PWA_PARAS *pp = (cu_PWA_PARAS*)&float_pp[i*pwa_paras_size];
        d_fx[i]=calEva(pp,parameter,d_paraList,d_mlk,i);
        //printf("%dgpu :: %.7f\n",i,pp->wu[0]);
        //printf("\nfx[%d]:%f\n",i,d_fx[i]);
        //fx[i]=calEva(pp,parameter,d_paraList,i);
    }
    
    //if(i==1)
    //{
        //printf("pp[0]:%f pp[end]:%f parameter[0]:%d parameter[16]:%d paraList[0]:%f \n",float_pp[0],float_pp[numElements*sizeof(cu_PWA_PARAS)/sizeof(double)-1],parameter[0],parameter[16],d_paraList[0]);
    //}
}

int host_store_fx(double *d_float_pp,int *h_parameter,double *h_paraList,int para_size, double *h_fx,double * h_mlk,int numElements,int begin)
{
    double *d_fx;
    CUDA_CALL(cudaMalloc((void **)&(d_fx),numElements * sizeof(double)));
    //std::cout << __LINE__ << endl;
    int *d_parameter;
    CUDA_CALL(cudaMalloc((void **)&(d_parameter),18 * sizeof(int)));
    CUDA_CALL(cudaMemcpy(d_parameter , h_parameter, 18*sizeof(int), cudaMemcpyHostToDevice));
    //cout << "\nd_parameter[16]" <<h_parameter[16] << endl;
    //std::cout << __LINE__ << endl;
    //std::cout << "d_paralist[0]: "<< h_paraList[0] << std::endl;
    //std::cout << "paralist[0]: "<< paraList[0] << std::endl;
    double *d_paraList;
    CUDA_CALL(cudaMalloc((void **)&(d_paraList),para_size * sizeof(double)));
    CUDA_CALL(cudaMemcpy(d_paraList , h_paraList, para_size * sizeof(double), cudaMemcpyHostToDevice));
    //cout << "\nd_paraList : " <<h_paraList[0] << endl;
    //std::cout << __LINE__ << endl;

    //init mlk
    double *d_mlk=NULL;
    CUDA_CALL(cudaMalloc( (void **)&(d_mlk),(h_parameter[16]+h_parameter[17])*h_parameter[15]*sizeof(double) ));
    //ut << "nAmps="<< h_parameter[15] << "iEnd=" << (h_parameter[16]+h_parameter[17]) << endl;
    int threadsPerBlock = 256;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    printf("%d\n",sizeof(double2)*h_parameter[15]*(7+h_parameter[15])*numElements );
    kernel_store_fx<<<blocksPerGrid, threadsPerBlock>>>(d_float_pp, d_parameter,d_paraList,d_fx,d_mlk, numElements,begin);
     //std::cout << __LINE__ << endl;
    CUDA_CALL(cudaGetLastError());
    //CUDA_CALL(cudaMemcpy(h_fx , d_fx, numElements * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_mlk , d_mlk, (h_parameter[16]+h_parameter[17])*h_parameter[15]*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h_fx , d_fx, numElements * sizeof(double), cudaMemcpyDeviceToHost));

    //free memory
    //CUDA_CALL(cudaFree(d_float_pp));
    CUDA_CALL(cudaFree(d_fx));
    CUDA_CALL(cudaFree(d_parameter));
    CUDA_CALL(cudaFree(d_paraList));
    CUDA_CALL(cudaFree(d_mlk));

    //ofstream cout("data_fx_cal");
    //std::cout << __LINE__ << endl;
    //for(int i=begin;i<numElements;i++)
    //{
        //cout << h_fx[i] << endl;
    //}
    //cout.close();
    return 0;
}

void cu_malloc_h_pp(double *h_float_pp,double *&d_float_pp,int length)
{
    int array_size = sizeof(cu_PWA_PARAS) / sizeof(double) * length;
    int mem_size = array_size * sizeof(double);
    CUDA_CALL(cudaMalloc((void **)&d_float_pp, mem_size));
    CUDA_CALL(cudaMemcpy(d_float_pp , h_float_pp, mem_size, cudaMemcpyHostToDevice));
}

