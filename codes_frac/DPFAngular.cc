#include "DPFAngular.h"
#include "TComplex.h"

//ClassImp(DPFAngular);


Double_t DPFAngular::scalar(Double_t *a1,Double_t *a2)const
{
    //  cout<<"haha: "<< __LINE__ << endl;
    Double_t scal=0;
    //	Double_t fDel[4][4];
    //	for(Int_t i=0;i<4;i++){
    //		for(Int_t j=0;j<4;j++){
    //			if(i==j){
    //				if(i<3) fDel[i][j]=-1;
    //				else fDel[i][j]=1;
    //			}
    //			else fDel[i][j]=0;
    //		}
    //	}

    for(Int_t i=0;i<4;i++){
        scal+=a1[i]*a2[i]*_dp->fDel[i][i];
    }
    return scal;
}

Double_t DPFAngular::calculate0p(
        Double_t _p11, Double_t _p12, Double_t _p13, Double_t _p14,
        Double_t _p21, Double_t _p22, Double_t _p23, Double_t _p24,
        Double_t _p31, Double_t _p32, Double_t _p33, Double_t _p34,
        Double_t _p41, Double_t _p42, Double_t _p43, Double_t _p44,
        Double_t _p51, Double_t _p52, Double_t _p53, Double_t _p54,
        PWA_PARAS &pp) const {
    Double_t ap23[4],apv2[4],apv3[4];
    Double_t ak1[4],ak2[4],ak3[4],ak4[4],ak5[4];
    ak1[0]=_p11; ak1[1]=_p12; ak1[2]=_p13; ak1[3]=_p14;
    ak2[0]=_p21; ak2[1]=_p22; ak2[2]=_p23; ak2[3]=_p24;
    ak3[0]=_p31; ak3[1]=_p32; ak3[2]=_p33; ak3[3]=_p34;
    ak4[0]=_p41; ak4[1]=_p42; ak4[2]=_p43; ak4[3]=_p44;
    ak5[0]=_p51; ak5[1]=_p52; ak5[2]=_p53; ak5[3]=_p54;

    Double_t ak23[4],ak45[4],ar[4],ak23u[4],ak45u[4],aru[4];
    Double_t ak12[4],ak13[4],akv2m3[4],akv3m2[4],ak12u[4],ak13u[4],akv2m3u[4],akv3m2u[4];
    Double_t delv2w[4][4],delv3w[4][4],ak12w[4],ak13w[4],akv2m3w[4],akv3m2w[4];
    Double_t del45w[4][4],t2wvf[4][4],wt[4];
    Double_t arw[4];
    Double_t ak23wu[4];
    Double_t fPCMS[4];
    Double_t arwdarw,tmp0;
    fPCMS[0]=fPCMS[1]=fPCMS[2]=0; fPCMS[3]=psi_mass;
    Double_t fFUD=0.22/3.0;
    Double_t fMK2=kaon_mass*kaon_mass,fMP2=pion_mass*pion_mass, fMPsi2=psi_mass*psi_mass;
    Double_t q2r45,b1q2r45,qjvf2;
    Double_t qf2xx,tmpf2,qjv2,qjv3,qbv2,qbv3;
    Double_t q2r23;
    Double_t tmp12w,tmp13w,tmpv2m3,tmpv3m2;
    Double_t temp,tmp1,tmp2,tmp;
    Double_t (*fDel)[4],(*fGel)[4],(*E)[4][4][4],(*G1)[4][4][4],(*G3)[4][4][4][4][4],t2wf[4][4],del23w[4][4],t2wfu[4][4],w2p1u[4],ttfw[4][4][4],t2p3[4][4],t4wvf[4][4][4][4];
    Double_t t2v2[4][4],t2v3[4][4],t2b3[4][4],t2b2[4][4];
    Double_t w1p12_1u[4],w1p13_1u[4],w1p12_2u[4],w1p13_2u[4];

    fDel = _dp->fDel;
    fGel = _dp->fGel;
    E = _dp->E;
    G1 = _dp->G1;
    G3 = _dp->G3;

    //#include dpfangular_cc_commet.cc 影响程序的格式，因此这部分被注释的内容放到了dpfangular_cc_commet中


    for(Int_t i=0;i<4;i++){
        ap23[i]=ak2[i]+ak3[i];
        apv2[i]=ak1[i]+ak2[i];
        apv3[i]=ak1[i]+ak3[i];
        ak23[i]=ak2[i]-ak3[i];
        ak45[i]=ak4[i]-ak5[i];
        //  cout<<"haha: ak45[i] = "<< ak45[i] << endl;
        ar   [i]=ak1[i]-ap23[i];
        ak23u[i]=ak23[i]*fDel[i][i];
        ak45u[i]=ak45[i]*fDel[i][i];
        //  cout<<"haha: ak45u[i] = "<< ak45u[i] << endl;
        aru[i]=ar[i]*fDel[i][i];
        //    arw[i] =fDel[i][j]*aru[j];
        //    akx[i]=ak4[i]-ak5[i];
        //    q[i]=fPCMS[i]-ak1[i]-ak2[i];
        //    akd[i]=ak1[i]-ak2[i];
        //    akxq[i]=2.0*(ak1[i]+ak2[i])-fPCMS[i];
    }
    //  fSX=scalar(akx,akx);
    pp.sv=scalar(ak1,ak1);
    //	cout<<"pp.sv="<<pp.sv<<endl;
    pp.s23=scalar(ap23,ap23);
    pp.sv2=scalar(apv2,apv2);
    pp.sv3=scalar(apv3,apv3);
    // 0+ contribution
    q2r45=0.25*pp.sv-_dp->_m2[0];
    //	cout<<"q2r45="<<q2r45<<endl;
    b1q2r45=sqrt(q2r45+fFUD);
    //	cout<<"b1q2r45="<<b1q2r45<<endl;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            del45w[i][j]=fDel[i][j]-ak1[i]*ak1[j]/pp.sv;
            //  cout<<"haha: fDel = "<< fDel[i][j] << endl;
            //  cout<<"haha: ak1 = " << ak1[i] << endl;
            //  cout<<"haha: pp.sv = "  << pp.sv << endl;
            //  cout<<"haha: del45w[i][j] = "<< del45w[i][j] << endl;
        }
    }
    for(Int_t i=0;i<4;i++){
        pp.wu[i]=0.0;
        for(Int_t j=0;j<4;j++){
            pp.wu[i]=pp.wu[i]+del45w[i][j]*ak45u[j];
            //  cout<<"haha: wu[i] = "<< wu[i] << endl;
        }
        pp.wu[i]=pp.wu[i]/b1q2r45;
        wt[i]=pp.wu[i]*fDel[i][i];
    }
    for(Int_t i=0;i<4;i++){
        arw[i]=0.0;
        for(Int_t j=0;j<4;j++){
            arw[i]=arw[i]+fGel[i][j]*aru[j];
        }
    }
    arwdarw=scalar(arw,arw);
    tmp0=(-arwdarw)/3;
    //	cout<<"tmp0="<<tmp0<<endl;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            t2wvf[i][j]=arw[i]*arw[j]+tmp0*fGel[i][j];
            //	cout<<"t2wvf[i][j]="<<t2wvf[i][j]<<endl;
        }
    }
    qjvf2=0.25*(pow((fMPsi2+pp.sv-pp.s23),2))/fMPsi2-pp.sv;
    //	cout<<"qjvf2="<<qjvf2<<endl;
    pp.b2qjvf2=sqrt(pow(qjvf2,2)+3.*qjvf2*fFUD+9.*pow(fFUD,2));
    //	cout<<"b2qjvf2="<<b2qjvf2<<endl;
    //0+ contribution
    for(Int_t i=0;i<2;i++){
        pp.w0p22[i]=0.0;
        for(Int_t j=0;j<4;j++){
            pp.w0p22[i]=pp.w0p22[i]+t2wvf[i][j]*wt[j];
            //	cout<<"w0p22[i]="<<w0p22[i]<<endl;
        }
    }
    //  cr0p1=propogator(mass0,width0,fSX);
    //  pp.cr0p11=cr0p1/b2qjvf2;
    //   2+ contribution
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            del23w[i][j]=fDel[i][j]-ap23[i]*ap23[j]/pp.s23;
            //	cout<<"del23w[i][j]="<<del23w[i][j]<<endl;
        }
    }

    for(Int_t i=0;i<4;i++){
        pp.ak23w[i]=0.0;
        for(Int_t j=0;j<4;j++){
            pp.ak23w[i]=pp.ak23w[i]+del23w[i][j]*ak23u[j];
            //	cout<<"pp.ak23w[i]="<<pp.ak23w[i]<<endl;
        }
    }
    qf2xx=0.25*pp.s23-_dp->_m2[2];
    pp.b2qf2xx=sqrt(pow(qf2xx,2)+3.*qf2xx*fFUD+9.*pow(fFUD,2));
    pp.b4qjvf2=sqrt(pow(qjvf2,4)+10*pow(qjvf2,3)*fFUD+135*pow(qjvf2,2)*pow(fFUD,2)+1575*qjvf2*pow(fFUD,3)+11025*pow(fFUD,4));
    tmpf2=-scalar(pp.ak23w,pp.ak23w)/3;
    //	cout<<"tmpf2="<<tmpf2<<endl;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            t2wf[i][j]=pp.ak23w[i]*pp.ak23w[j]+tmpf2*del23w[i][j];
            //	cout<<"t2wf[i][j]="<<t2wf[i][j]<<endl;
            t2wfu[i][j]=t2wf[i][j]*fDel[i][i]*fDel[j][j];
        }
    }
    for(Int_t i=0;i<4;i++){
        pp.w2p1[i]=0.0;
        for(Int_t j=0;j<4;j++){
            pp.w2p1[i]=pp.w2p1[i]+t2wf[i][j]*wt[j];
            //	cout<<"w2p1[i]="<<w2p1[i]<<endl;
        }
        w2p1u[i]=pp.w2p1[i]*fDel[i][i];
    }
    for(Int_t i=0;i<2;i++){
        pp.w2p2[i]=0.0;
        for(Int_t j=0;j<4;j++){
            pp.w2p2[i]=pp.w2p2[i]+t2wvf[i][j]*w2p1u[j];
        }
    }
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            for(Int_t k=0;k<4;k++){
                ttfw[i][j][k]=t2wfu[i][j]*wt[k];
            }
        }
    }
    //  cout<<"haha: "<< __LINE__ << endl;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            t2p3[i][j]=0.0;
            for(Int_t k1=0;k1<3;k1++){
                for(Int_t k2=0;k2<3;k2++){
                    t2p3[i][j]=t2p3[i][j]+(E[i][k1][3][k2]*ttfw[k1][j][k2]+E[j][k1][3][k2]*ttfw[k1][i][k2]);
                }
            }
        }
    }
    for(Int_t i=0;i<2;i++){
        pp.w2p3[i]=0.0;
        for(Int_t j1=0;j1<3;j1++){
            for(Int_t j2=0;j2<3;j2++){
                for(Int_t j3=0;j3<4;j3++){
                    pp.w2p3[i]=pp.w2p3[i]+E[i][j1][j2][3]*t2wvf[j1][j3]*t2p3[j2][j3];
                }
            }
        }
    }
    for(Int_t i=0;i<2;i++){
        pp.w2p4[i]=0.0;
        for(Int_t j1=0;j1<4;j1++){
            for(Int_t j2=0;j2<4;j2++){
                for(Int_t j3=0;j3<4;j3++){
                    for(Int_t j4=0;j4<4;j4++){
                        for(Int_t j5=0;j5<4;j5++){
                            pp.w2p4[i]=pp.w2p4[i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*ttfw[j3][j4][j5];
                        }
                    }
                }
            }
        }
    }
    tmp1=(-arwdarw)/7;
    tmp2=(pow(arwdarw,2))/35;
    for(Int_t i=0;i<2;i++){
        for(Int_t j=0;j<4;j++){
            for(Int_t k=0;k<4;k++){
                for(Int_t l=0;l<4;l++){
                    t4wvf[i][j][k][l]=arw[i]*arw[j]*arw[k]*arw[l]+tmp1*(fGel[i][j]*arw[k]*arw[l]+fGel[j][k]*arw[i]*arw[l]+fGel[k][i]*arw[j]*arw[l]+fGel[i][l]*arw[j]*arw[k]+fGel[j][l]*arw[k]*arw[i]+fGel[k][l]*arw[i]*arw[j])+tmp2*G1[i][j][k][l];
                }
            }
        }
    }
    for(Int_t i=0;i<2;i++){
        pp.w2p5[i]=0.0;
        for(Int_t j1=0;j1<4;j1++){
            for(Int_t j2=0;j2<4;j2++){
                for(Int_t j3=0;j3<4;j3++){
                    pp.w2p5[i]=pp.w2p5[i]+t4wvf[i][j1][j2][j3]*ttfw[j2][j3][j1];
                }
            }
        }
    }
    // 1+ and 1- contribution
    qjv2=0.25*pow((fMPsi2+pp.sv2-_dp->_m2[2]),2)/fMPsi2-pp.sv2;
    qjv3=0.25*pow((fMPsi2+pp.sv3-_dp->_m2[2]),2)/fMPsi2-pp.sv3;
    qbv2=0.25*pow((pp.sv2+pp.sv-_dp->_m2[2]),2)/pp.sv2-pp.sv;
    qbv3=0.25*pow((pp.sv3+pp.sv-_dp->_m2[2]),2)/pp.sv3-pp.sv;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            delv2w[i][j]=fDel[i][j]-apv2[i]*apv2[j]/pp.sv2;
            delv3w[i][j]=fDel[i][j]-apv3[i]*apv3[j]/pp.sv3;
        }
    }
    for(int i=0;i<4;i++){
        ak12[i]=ak1[i]-ak2[i];
        ak13[i]=ak1[i]-ak3[i];
        akv2m3[i]=apv2[i]-ak3[i];
        akv3m2[i]=apv3[i]-ak2[i];
        ak12u[i]=ak12[i]*fDel[i][i];
        ak13u[i]=ak13[i]*fDel[i][i];
        akv2m3u[i]=akv2m3[i]*fDel[i][i];
        akv3m2u[i]=akv3m2[i]*fDel[i][i];
    }
    for(int i=0;i<4;i++){
        ak12w[i]=0.0;
        ak13w[i]=0.0;
        akv2m3w[i]=0.0;
        akv3m2w[i]=0.0;
        for(int j=0;j<4;j++){
            ak12w[i]=ak12w[i]+delv2w[i][j]*ak12u[j];
            ak13w[i]=ak13w[i]+delv3w[i][j]*ak13u[j];
            akv2m3w[i]=akv2m3w[i]+fGel[i][j]*akv2m3u[j];
            akv3m2w[i]=akv3m2w[i]+fGel[i][j]*akv3m2u[j];
        }
    }
    // 1+ contribution
    pp.b2qjv2=sqrt(pow(qjv2,2)+3.*qjv2*fFUD+9.*pow(fFUD,2));
    pp.b2qjv3=sqrt(pow(qjv3,2)+3.*qjv3*fFUD+9.*pow(fFUD,2));
    pp.b2qbv2=sqrt(pow(qbv2,2)+3.*qbv2*fFUD+9.*pow(fFUD,2));
    pp.b2qbv3=sqrt(pow(qbv3,2)+3.*qbv3*fFUD+9.*pow(fFUD,2));
    tmp12w=-scalar(ak12w,ak12w)/3.0;
    tmp13w=-scalar(ak13w,ak13w)/3.0;
    tmpv2m3=-scalar(akv2m3w,akv2m3w)/3.0;
    tmpv3m2=-scalar(akv3m2w,akv3m2w)/3.0;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            t2v2[i][j]=ak12w[i]*ak12w[j]+tmp12w*delv2w[i][j];
            t2v3[i][j]=ak13w[i]*ak13w[j]+tmp13w*delv3w[i][j];
            t2b3[i][j]=akv2m3w[i]*akv2m3w[j]+tmpv2m3*fGel[i][j];
            t2b2[i][j]=akv3m2w[i]*akv3m2w[j]+tmpv3m2*fGel[i][j];
        }
    }
    for(int i=0;i<4;i++){
        pp.w1p12_1[i]=0.0;
        pp.w1p13_1[i]=0.0;
        pp.w1p12_2[i]=0.0;
        pp.w1p13_2[i]=0.0;
        for(int j=0;j<4;j++){
            pp.w1p12_1[i]=pp.w1p12_1[i]+delv2w[i][j]*wt[j];
            pp.w1p13_1[i]=pp.w1p13_1[i]+delv3w[i][j]*wt[j];
            pp.w1p12_2[i]=pp.w1p12_2[i]+t2v2[i][j]*wt[j];
            pp.w1p13_2[i]=pp.w1p13_2[i]+t2v3[i][j]*wt[j];
        }
        w1p12_1u[i]=pp.w1p12_1[i]*fDel[i][i];
        w1p13_1u[i]=pp.w1p13_1[i]*fDel[i][i];
        w1p12_2u[i]=pp.w1p12_2[i]*fDel[i][i];
        w1p13_2u[i]=pp.w1p13_2[i]*fDel[i][i];
    }
    for(int i=0;i<2;i++){
        pp.w1p12_3[i]=0.0;
        pp.w1p13_2[i]=0.0;
        pp.w1p12_4[i]=0.0;
        pp.w1p13_4[i]=0.0;
        for(int j=0;j<4;j++){
            pp.w1p12_3[i]=pp.w1p12_3[i]+t2b3[i][j]*w1p12_1u[j];
            pp.w1p13_3[i]=pp.w1p13_3[i]+t2b2[i][j]*w1p13_1u[j];
            pp.w1p12_4[i]=pp.w1p12_4[i]+t2b3[i][j]*w1p12_2u[j];
            pp.w1p13_4[i]=pp.w1p13_4[i]+t2b2[i][j]*w1p13_2u[j];
        }
    }
    //  1- contribution
    pp.b1qjv2=sqrt(qjv2+fFUD);
    pp.b1qjv3=sqrt(qjv3+fFUD);
    pp.b1qbv2=sqrt(qbv2+fFUD);
    pp.b1qbv3=sqrt(qbv3+fFUD);
    for(int i=0;i<2;i++){
        pp.w1m12[i]=0.0;
        pp.w1m13[i]=0.0;
        for(int j1=0;j1<3;j1++){
            for(int j2=0;j2<3;j2++){
                for(int j3=0;j3<3;j3++){
                    for(int j4=0;j4<3;j4++){
                        tmp=E[i][j1][j2][4]*E[j2][j3][j4][4];
                        pp.w1m12[i]=pp.w1m12[i]+tmp*akv2m3w[j1]*ak12w[j3]*wt[j4];
                        pp.w1m13[i]=pp.w1m13[i]+tmp*akv3m2w[j1]*ak13w[j3]*wt[j4];
                    }
                }
            }
        }
    }
    // phi(1650)f0(980)  contribution
    q2r23=0.25*pp.s23-_dp->_m2[2];
    pp.b1q2r23=sqrt(q2r23+fFUD);
    for(int i=0;i<4;i++){
        ak23wu[i]=pp.ak23w[i]*fDel[i][i];
    }
    for(int i=0;i<2;i++){
        pp.wpf22[i]=0.0;
        for(int j=0;j<4;j++){
            pp.wpf22[i]=pp.wpf22[i]-t2wvf[i][j]*ak23wu[j];
        }
    }
}
