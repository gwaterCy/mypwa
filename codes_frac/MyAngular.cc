#include "MyAngular.h"
#include "DPFPWAPoint.h"
#include "DPFPWAPoint.cc"
#include "TComplex.h"
const Double_t mpsip=3.686,mka=0.493677,mpi=0.13957;
Double_t fDel[4][4], fGel[4][4],E[4][4][4][4],G1[4][4][4][4],G3[4][4][4][4][4][4];
Double_t   _M   = mpsip;
Double_t _m[4] = {0.493677, 0.493677, 0.13957, 0.13957};

//ClassImp(MyAngular);

Double_t MyAngular::scalar(Double_t *a1,Double_t *a2)const
 {
       for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            if(i==j){
                if(i<3) 
                    fDel[i][j]=-1;
                else 
                    fDel[i][j]=1;
            } else { 
                fDel[i][j]=0;
            }
        }
    }

	Double_t scal=0;
	for(Int_t i=0;i<4;i++){
		scal+=a1[i]*a2[i]*(fDel[i][i]);
	}
	return scal;
} 

void MyAngular::calculate0p(Double_t ak1[], Double_t ak2[], Double_t ak3[], Double_t ak4[], Double_t ak5[], 
			     Double_t wu[], Double_t w0p22[],
			     Double_t w2p1[], Double_t w2p2[], Double_t w2p3[], Double_t w2p4[], Double_t w2p5[],
			     Double_t ak23w[],Double_t wpf22[], 
			     Double_t& b2qjvf2, Double_t& b2qf2xx, Double_t& b4qjvf2, Double_t& b1q2r23,
			     Double_t& b2qbv2, Double_t& b2qbv3, Double_t& b2qjv2, Double_t& b2qjv3,
			     Double_t w1p12_2[], Double_t w1p13_2[], Double_t w1p12_3[], Double_t w1p13_3[],
			     Double_t w1p12_4[], Double_t w1p13_4[],
			     Double_t& b1qjv2, Double_t& b1qbv2, Double_t& b1qjv3, Double_t& b1qbv3,
			     Double_t w1m12[], Double_t w1m13[]) const{     

  Double_t ap23[4],apv2[4],apv3[4],ak23[4],ak45[4],ar[4],ak23u[4],ak45u[4],aru[4];
  Double_t ak12[4],ak13[4],akv2m3[4],akv3m2[4],ak12u[4],ak13u[4],akv2m3u[4],akv3m2u[4];
  Double_t delv2w[4][4],delv3w[4][4],ak12w[4],ak13w[4],akv2m3w[4],akv3m2w[4];
  Double_t del45w[4][4],t2wvf[4][4],t2wvfu[4][4],t2wvfuu[4][4],wt[4];
  Double_t arw[4];
  Double_t ak23wu[4];
  //,tf4[4][4][4][4];
  Double_t fPCMS[4],sv,s23,sv2,sv3,arwdarw,tmp0;
  fPCMS[0]=fPCMS[1]=fPCMS[2]=0; fPCMS[3]=psi_mass;
  Double_t fMK2=kaon_mass*kaon_mass,fMP2=pion_mass*pion_mass, fMPsi2=psi_mass*psi_mass;
  Double_t q2r45,b1q2r45,qjvf2;
  Double_t qf2xx,tmpf2,qjv2,qjv3,qbv2,qbv3,q2r23;
  Double_t tmp12w,tmp13w,tmpv2m3,tmpv3m2;
  Double_t temp,tmp1,tmp2,tmp;
  Double_t fDel[4][4], fGel[4][4],E[4][4][4][4],G1[4][4][4][4],G3[4][4][4][4][4][4];
  Double_t t2wf[4][4],del23w[4][4],t2wfu[4][4],t2wfuu[4][4],w2p1u[4],ttfw[4][4][4],t2p3[4][4],t4wvf[4][4][4][4];
  Double_t t2v2[4][4],t2v3[4][4],t2b3[4][4],t2b2[4][4];
  Double_t w1p12_1[4],w1p13_1[4],w1p12_1u[4],w1p13_1u[4],w1p12_2u[4],w1p13_2u[4];
  //Double_t w1p12_3[2],w1p13_3[2],w1p12_4[2],w1p13_4[2];
  //Double_t w1m12[2],w1m13[2];
  TComplex cr0p11;
  TComplex ca2p1;
  TComplex cw2p11;
  TComplex cw2p12;
  TComplex cw2p15;
  TComplex cd0p1;
  TComplex cd0p2;
  TComplex cd2p1;
  TComplex cw;
  TComplex c1p12_12,c1p13_12,c1p12_13,c1p13_13,c1p12_14,c1p13_14;
  TComplex cr1m12_1,cr1m13_1;
  TComplex crpf1,crpf2;

      for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            if(i==j){
                if(i<3) 
                    fDel[i][j]=-1;
                else 
                    fDel[i][j]=1;
            } else { 
                fDel[i][j]=0;
            }
        }
    }
    
    
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            if(i==j){
                if(i<3) 
                    fGel[i][j]=-1;
                else 
                    fGel[i][j]=0;
            } else {
                fGel[i][j]=0;
            }
        }
    }
    //  cout<<"haha: "<< __LINE__ << endl;
    E[0][1][2][3]= 1.0;
    E[0][1][3][2]=-1.0;
    E[0][2][1][3]=-1.0;
    E[0][2][3][1]= 1.0;
    E[0][3][1][2]= 1.0;
    E[0][3][2][1]=-1.0;
    E[1][0][2][3]=-1.0;
    E[1][0][3][2]= 1.0;
    E[1][2][0][3]= 1.0;
    E[1][2][3][0]=-1.0;
    E[1][3][0][2]=-1.0;
    E[1][3][2][0]= 1.0;
    E[2][0][1][3]= 1.0;
    E[2][0][3][1]=-1.0;
    E[2][1][0][3]=-1.0;
    E[2][1][3][0]= 1.0;
    E[2][3][0][1]= 1.0;
    E[2][3][1][0]=-1.0;
    E[3][0][1][2]=-1.0;
    E[3][0][2][1]= 1.0;
    E[3][1][0][2]= 1.0;
    E[3][1][2][0]=-1.0;
    E[3][2][0][1]=-1.0;
    E[3][2][1][0]= 1.0;
    for(Int_t I =0;I<4;I++){
        for(Int_t J=0;J<4;J++){
            for(Int_t K=0;K<4;K++){
                for(Int_t L=0;L<4;L++){
		 // G1[I][J][K][L] = (fGel[J][K]*fGel[I][L] + fGel[K][I]*fGel[J][L])/2 -(fGel[I][J]*fGel[K][L])/3 ;
           G1[I][J][K][L] = fGel[I][J]*fGel[K][L] + fGel[J][K]*fGel[I][L] +fGel[K][I]*fGel[J][L]; // last (...), i.e. combination of \tilde{g} in Eq.(37) of PRD48,1225
               }
            }
        }
    }// Eq.(20)
    for(Int_t I1=0;I1<4;I1++){
        for(Int_t I2=0;I2<4;I2++){
            for(Int_t I3=0;I3<4;I3++){
                for(Int_t I4=0;I4<4;I4++){
                    for(Int_t I5=0;I5<4;I5++){
                        for(Int_t I6=0;I6<4;I6++){
                            G3[I1][I2][I3][I4][I5][I6] =
                                ( fGel[I1][I2]*fGel[I4][I5]*fGel[I3][I6] +
                                  fGel[I1][I2]*fGel[I5][I6]*fGel[I3][I4] +
                                  fGel[I1][I2]*fGel[I4][I6]*fGel[I3][I5] +
                                  fGel[I1][I3]*fGel[I4][I6]*fGel[I2][I5] +
                                  fGel[I1][I3]*fGel[I4][I5]*fGel[I2][I6] +
                                  fGel[I1][I3]*fGel[I5][I6]*fGel[I2][I4] +
                                  fGel[I2][I3]*fGel[I5][I6]*fGel[I1][I4] +
                                  fGel[I2][I3]*fGel[I4][I5]*fGel[I1][I6] +
                                  fGel[I2][I3]*fGel[I4][I6]*fGel[I1][I5] )/15.0 -
                                ( fGel[I1][I4]*fGel[I2][I5]*fGel[I3][I6] +
                                  fGel[I1][I4]*fGel[I2][I6]*fGel[I3][I5] +
                                  fGel[I1][I5]*fGel[I2][I4]*fGel[I3][I6] +
                                  fGel[I1][I5]*fGel[I2][I6]*fGel[I3][I4] +
                                  fGel[I1][I6]*fGel[I2][I5]*fGel[I3][I4] +
                                  fGel[I1][I6]*fGel[I2][I4]*fGel[I3][I5] )/6.0;
                        }
                    }
                }
            }
        }
    }// Eq. (21)

  // 1 phi 2 pi+ 3 pi- 4 K+ 5 K-
  // ak1 ~ ak5, four momentum from input, ak1 is the sum of 3 and 4, etc. the phi's momentum
  for(Int_t i=0;i<4;i++){
    ap23[i]=ak2[i]+ak3[i]; // pi+ + pi-
    apv2[i]=ak1[i]+ak2[i]; // phi + pi+
    apv3[i]=ak1[i]+ak3[i]; // phi + pi-
    ak23[i]=ak2[i]-ak3[i]; // r_34 (relative difference between pi+ and pi-)
    ak45[i]=ak4[i]-ak5[i]; // r_12 (relative difference between K+ and K-)
    ar[i]=ak1[i]-ap23[i];  // r^\mu (relative difference between phi and two pions)
    // 
    ak23u[i]=ak23[i]*fDel[i][i];
    ak45u[i]=ak45[i]*fDel[i][i];
    aru[i]=ar[i]*fDel[i][i];// r_\mu(phi f) 
  }
  sv=scalar(ak1,ak1);      // p^2(phi)
  s23=scalar(ap23,ap23);   // p^2(f)
  sv2=scalar(apv2,apv2);  //  p^2(123)
  sv3=scalar(apv3,apv3); //   p^2(124)

  // 0+ contribution
  // prepare tensors for Eq.(36)
  q2r45=sv*0.25 - (_m[0])*(_m[0]); // q2r45 is Q^2_abc, Eq.(13), phi(a)-> K+(b) K-(c)
  b1q2r45=sqrt(2/(q2r45+fFUD)); // b1q2r45 is B1(Q_abc), Eq.(14), it can be ignored for phi is narrow, i.e. it should be always close to a constant 
  // when final states' momenta change 
  //cout << "b1q2r45 " << b1q2r45 << endl;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      del45w[i][j]=fDel[i][j]-ak1[i]*ak1[j]/sv;
      //cout << "del45w " << del45w[i][j] << endl;
    }
  }// del45w is \tilde{g}^{\mu\nu}(\phi)
  
  for(Int_t i=0;i<4;i++){
    wt[i]=0.0;
    for(Int_t j=0;j<4;j++){
      wt[i]=wt[i]+del45w[i][j]*ak45u[j];
    }
    wt[i]=wt[i]*b1q2r45; // Eq. (10), \tilde{r}^\mu = t^{(1)\mu}(12)
    wu[i]=wt[i]*fDel[i][i]; // \tilde{r}^\mu = t^{(1)}_\mu(12)
  }
  // prepare tensors for Eq.(37)
  for(Int_t i=0;i<4;i++){
    arw[i]=0.0;
    for(Int_t j=0;j<4;j++){
      arw[i]=arw[i]+fGel[i][j]*aru[j];  // \tilde(r)^\mu(phi f0) 
    }
  } 
  arwdarw=scalar(arw,arw); 
  tmp0=(-arwdarw)/3;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2wvf[i][j]=arw[i]*arw[j]+tmp0*fGel[i][j]; // \tilde{T}^{(2)\mu\nu}(phi f0), from Eq.(11)
      //cout << "t2wvf " << t2wvfuu[i][j] << endl;
    }
  }
  qjvf2=0.25*(pow((mpsip*mpsip+sv-s23),2))/(mpsip*mpsip)-sv; // Eq.(13), a is psi, b is phi, c is resonance
  b2qjvf2=sqrt(13./(pow(qjvf2,2)+3.*qjvf2*fFUD+9.*pow(fFUD,2))); // B2(Qabc) in Eq.(15), can be neglected since psi is very narrow
  //cout << "b2qjvf2 " << b2qjvf2 << endl;

  //0+ contribution
  for(Int_t i=0;i<2;i++){
    w0p22[i]=0.0;
    for(Int_t j=0;j<4;j++){
      w0p22[i]=w0p22[i]+t2wvf[i][j]*wu[j]; // the \tilde{T}^{(2)\mu\nv}(phi f0) \tilde{t}^{(1)}_{\nu}(12)  in Eq.(37), only prepare \mu = 1 and 2
    }
  }
  // for(int i=0; i<4; i++) cout << "w0pp2 " << w0p22[i] << endl;
  // Till now, the tensors in Eq.(36) [wt] and (37) [w0p22] have been prepared. 
  // The missing parts are propagators, etc. BW.

  //   2+ contribution
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      del23w[i][j]=fDel[i][j]-ap23[i]*ap23[j]/s23; // \tilde{g}^{\mu\nu}(f2) in Eq.(11)
    }
  }
  for(Int_t i=0;i<4;i++){
    ak23w[i]=0.0;
    for(Int_t j=0;j<4;j++){
      ak23w[i]=ak23w[i]+del23w[i][j]*ak23u[j]; // \tilde{r}^{\mu}(f2) in Eq.(11), r(a)-> pi+(b) pi-(c)
    }
  }
  qf2xx=0.25*s23-(_m[2])*(_m[2]);
  b2qf2xx=sqrt(13./(pow(qf2xx,2)+3.*qf2xx*fFUD+9.*pow(fFUD,2))); // B2 in Eq.(15) of pi+ pi-
  b4qjvf2=sqrt(12746./(pow(qjvf2,4)+10*pow(qjvf2,3)*fFUD+135*pow(qjvf2,2)*pow(fFUD,2)+1575*qjvf2*pow(fFUD,3)+11025*pow(fFUD,4))); // B4 in Eq.(17) of psi
  tmpf2=-scalar(ak23w,ak23w)/3;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2wf[i][j]=(ak23w[i]*ak23w[j]+tmpf2*del23w[i][j])*b2qf2xx; // \tilde{t}^{(2)\mu\nu}(f2) in Eq.(11), \tilde{t}^{(2)\mu\nu}(34) in Eq.(41)
      t2wfu[i][j]=t2wf[i][j]*fDel[j][j];  // \tilde{t}^{(2)\mu}_{\nu}(34) 
      //cout << "t2wfu " << t2wfu[i][j] << endl;
    }
  }
  for(Int_t i=0;i<4;i++){
    w2p1[i]=0.0;
    for(Int_t j=0;j<4;j++){
      w2p1[i]=w2p1[i]+t2wf[i][j]*wu[j]; // \tilde{t}^{(2)\mu\nu}(34) \tilde{t}^{(1)}_{\nu}(12) the second and third terms in Eq.(42)
    }
    w2p1u[i]=w2p1[i]*fDel[i][i]; // \tilde{t}^{(2)}_{\mu\nu}(34) \tilde{t}^{(1)\nu}(12)
  }
  for(Int_t i=0;i<2;i++){
    w2p2[i]=0.0;
    for(Int_t j=0;j<4;j++){
      w2p2[i]=w2p2[i]+t2wvf[i][j]*w2p1u[j]; // \tilde{T}^{(2)\mu\alpha}(phi f2)\tilde{t}^{(2)}_{\alpha\nu}(34) \tilde{t}^{(1)\nu}(12)£¬ the first three terms in Eq.(42)
    }
  }
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      for(Int_t k=0;k<4;k++){
    	ttfw[i][j][k]=t2wfu[i][j]*wt[k]; // \tilde{t}^{(2)\lambda}_\delta(34) \tilde{t}^{(1)\nu}(12)
      }
    }
  }
  //  cout<<"haha: "<< __LINE__ << endl;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2p3[i][j]=0.0;
      for(Int_t k1=0;k1<3;k1++){
	for(Int_t k2=0;k2<3;k2++){
	  t2p3[i][j]=t2p3[i][j]+(E[i][k1][3][k2]*ttfw[k1][j][k2]+E[j][k1][3][k2]*ttfw[k1][i][k2]); // [...]p^\sigma_psi \tilde{t}^{(1)\nu}(12)
	    }
      }
    }
  }
  for(Int_t i=0;i<2;i++){
    w2p3[i]=0.0;
    for(Int_t j1=0;j1<3;j1++){
      for(Int_t j2=0;j2<3;j2++){
	for(Int_t j3=0;j3<4;j3++){
	  w2p3[i]=w2p3[i]+E[i][j1][j2][3]*t2wvf[j1][j3]*t2p3[j2][j3]; // epsilon^{\mu\alpha\beta\gamma}p_{\psi\alpha}\tilde^{(2)\lambda}_{\beta}(\phi f2)[...] in Eq. (43)
      	}
      }
    }
  }
  double test[4][4][4];
  for(int i=0;i<4;i++){
    for(int j=0; j<4; j++){
        for(int k=0; k<4; k++){
            test[i][j][k]=0;
            for(int m=0; m<4; m++){
                for(int n=0; n<4; n++){
                    test[i][j][k] = test[i][j][k] + ttfw[i][m][n]*fDel[m][j]*fDel[n][k]; }}}}}
  for(Int_t i=0;i<2;i++){
    w2p4[i]=0.0;
    for(Int_t j1=0;j1<4;j1++){
      for(Int_t j2=0;j2<4;j2++){
	for(Int_t j3=0;j3<4;j3++){
	  for(Int_t j4=0;j4<4;j4++){
	    for(Int_t j5=0;j5<4;j5++){
	      w2p4[i]=w2p4[i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*ttfw[j3][j4][j5]; //  P^{(3)\tilde{T}^{(2)}(\phi f2)\tilde{t}^{(2)}(34)\tilde{t}^{(1)}(12) Eq.(44)
	     // w2p4[i]=w2p4[i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*test[j3][j4][j5];
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
	  t4wvf[i][j][k][l]=arw[i]*arw[j]*arw[k]*arw[l]+tmp1*(fGel[i][j]*arw[k]*arw[l]+fGel[j][k]*arw[i]*arw[l]+fGel[k][i]*arw[j]*arw[l]+fGel[i][l]*arw[j]*arw[k]+fGel[j][l]*arw[k]*arw[i]+fGel[k][l]*arw[i]*arw[j])+tmp2*G1[i][j][k][l]; // defined in Eq(37) of PRD48,1225
	}
      }
    }
  }
  for(Int_t i=0;i<2;i++){
    w2p5[i]=0.0;
    for(Int_t j1=0;j1<4;j1++){
      for(Int_t j2=0;j2<4;j2++){
	for(Int_t j3=0;j3<4;j3++){
	  w2p5[i]=w2p5[i]+t4wvf[i][j1][j2][j3]*ttfw[j2][j3][j1]; // \tilde{T}^{(4)(\phi f2) \tilde{t}^{(1)(12) \tilde{t}^{(2)}(34), in Eq.(45)
	}
      }
    }
  }

  // 1+ and 1- contribution
  qjv2=0.25*pow((fMPsi2+sv2-(_m[2])*(_m[2])),2)/fMPsi2-sv2; // Q^2(abc) = Q^2(\psi' 123 \pi4) 
  qjv3=0.25*pow((fMPsi2+sv3-(_m[2])*(_m[2])),2)/fMPsi2-sv3; // Q^2(\psi 124 \pi3)
  qbv2=0.25*pow((sv2+sv-(_m[2])*(_m[2])),2)/sv2-sv;         // Q^2(123 12 \pi3)
  qbv3=0.25*pow((sv3+sv-(_m[2])*(_m[2])),2)/sv3-sv;         // Q^2(124 12 \pi4)
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      delv2w[i][j]=fDel[i][j]-apv2[i]*apv2[j]/sv2;          // \tilde{g}^{\mu\nu}(123)
      delv3w[i][j]=fDel[i][j]-apv3[i]*apv3[j]/sv3;          // \tilde{g}^{\mu\nu}(124)
    }
  }
  for(int i=0;i<4;i++){
    ak12[i]=ak1[i]-ak2[i];        // r^\mu(phi pi3) i.e. p_{12}-p_3             
    ak13[i]=ak1[i]-ak3[i];        // r^\mu(phi pi4) i.e. p_{12}-p_4
    akv2m3[i]=apv2[i]-ak3[i];     // r^\mu(rho' pi4) i.e. p_{123}-p_4
    akv3m2[i]=apv3[i]-ak2[i];     // r^\mu(rho' pi3) i.e. p_{124}-p_3
    ak12u[i]=ak12[i]*fDel[i][i];  // r_\mu(phi pi3)
    ak13u[i]=ak13[i]*fDel[i][i];  // r_\mu(phi pi4)
    akv2m3u[i]=akv2m3[i]*fDel[i][i]; // r_\mu(rho' pi4)
    akv3m2u[i]=akv3m2[i]*fDel[i][i]; // r_\mu(rho' pi3)
  }
  for(int i=0;i<4;i++){
    ak12w[i]=0.0;
    ak13w[i]=0.0;
    akv2m3w[i]=0.0;
    akv3m2w[i]=0.0;
    for(int j=0;j<4;j++){
      ak12w[i]=ak12w[i]+delv2w[i][j]*ak12u[j];      // \tilde{r}^\mu (phi pi3) 
      ak13w[i]=ak13w[i]+delv3w[i][j]*ak13u[j];      // \tilde{r}^\mu (phi pi4)
      akv2m3w[i]=akv2m3w[i]+fGel[i][j]*akv2m3u[j];  // \tilde{r}^\mu (b pi4)
      akv3m2w[i]=akv3m2w[i]+fGel[i][j]*akv3m2u[j];  // \tilde{r}^\mu (b pi3)
    }
  }
  // 1+ contribution
  b2qjv2=sqrt(13./(pow(qjv2,2)+3.*qjv2*fFUD+9.*pow(fFUD,2)));  // B2(Q psi b pi4)
  b2qjv3=sqrt(13./(pow(qjv3,2)+3.*qjv3*fFUD+9.*pow(fFUD,2)));  // B2(Q psi b pi3)
  b2qbv2=sqrt(13./(pow(qbv2,2)+3.*qbv2*fFUD+9.*pow(fFUD,2)));  // B2(Q b phi pi3)
  b2qbv3=sqrt(13./(pow(qbv3,2)+3.*qbv3*fFUD+9.*pow(fFUD,2)));  // B2(Q b phi pi4)
  tmp12w=-scalar(ak12w,ak12w)/3.0;
  tmp13w=-scalar(ak13w,ak13w)/3.0;
  tmpv2m3=-scalar(akv2m3w,akv2m3w)/3.0;
  tmpv3m2=-scalar(akv3m2w,akv3m2w)/3.0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      t2v2[i][j]=(ak12w[i]*ak12w[j]+tmp12w*delv2w[i][j])*b2qbv2;    // \tilde{t}^{(2)\mu\nu}(phi pi3)
      t2v3[i][j]=(ak13w[i]*ak13w[j]+tmp13w*delv3w[i][j])*b2qbv3;    // \tilde{t}^{(2)\mu\nu}(phi pi4)
      t2b3[i][j]=(akv2m3w[i]*akv2m3w[j]+tmpv2m3*fGel[i][j])*b2qjv2; // \tilde{T}^{(2)\mu\nu}(b pi4)
      t2b2[i][j]=(akv3m2w[i]*akv3m2w[j]+tmpv3m2*fGel[i][j])*b2qjv3; // \tilde{T}^{(2)\mu\nu}(b pi3)
    }
  }
  for(int i=0;i<4;i++){
    w1p12_1[i]=0.0;
    w1p13_1[i]=0.0;
    w1p12_2[i]=0.0;
    w1p13_2[i]=0.0;
    for(int j=0;j<4;j++){
      w1p12_1[i]=w1p12_1[i]+delv2w[i][j]*wu[j]; // \tilde{g}^{\mu\nu}(123) \tilde{t}^{(1)}_\nu(12) Eq.(47)
      w1p13_1[i]=w1p13_1[i]+delv3w[i][j]*wu[j]; // \tilde{g}^{\mu\nu}(124) \tilde{t}^{(1)}_\nu(12) Eq.(47)
      w1p12_2[i]=w1p12_2[i]+t2v2[i][j]*wu[j];   // \tilde{t}_{(2)\mu\nu}(phi 3) \tilde{t}^{(1\nu}(12) Eq.(48)
      w1p13_2[i]=w1p13_2[i]+t2v3[i][j]*wu[j];   // \tilde{t}^{(2)\mu\nu}(phi 4) \tilde{t}^{(1)\nu}(12) Eq.(48)
    }
   // cout << "w1p12_1 " << w1p12_1[i] << "  w1p13_1 " << w1p13_1[i] << "  w1p12_2 " << w1p12_2[i] << "  w1p13_2 " << w1p13_2[i] << endl;
    w1p12_1u[i]=w1p12_1[i]*fDel[i][i];
    w1p13_1u[i]=w1p13_1[i]*fDel[i][i];
    w1p12_2u[i]=w1p12_2[i]*fDel[i][i];
    w1p13_2u[i]=w1p13_2[i]*fDel[i][i];
  }
  for(int i=0;i<2;i++){
    w1p12_3[i]=0.0;
    w1p13_3[i]=0.0;
    w1p12_4[i]=0.0;
    w1p13_4[i]=0.0;
    for(int j=0;j<4;j++){
      w1p12_3[i]=w1p12_3[i]+t2b3[i][j]*w1p12_1u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 4} \tilde{g}_{(123)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
      w1p13_3[i]=w1p13_3[i]+t2b2[i][j]*w1p13_1u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 3} \tilde{g}_{(124)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
      w1p12_4[i]=w1p12_4[i]+t2b3[i][j]*w1p12_2u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 4} \tilde{t}_{(phi 3)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
      w1p13_4[i]=w1p13_4[i]+t2b2[i][j]*w1p13_2u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 3} \tilde{g}_{(phi 4)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
    }
   // cout << "w1p12_3 " << w1p12_3[i] << "  w1p13_3 " << w1p13_3[i] << "  w1p12_4 " << w1p12_4[i] << "  w1p13_4 " << w1p13_4[i] << endl;
  }
  //  1- contribution
  b1qjv2=sqrt(2./(qjv2+fFUD)); // B1(Q psi rho' pi4), rho'(123) 
  b1qjv3=sqrt(2./(qjv3+fFUD)); // B1(Q psi rho' pi3), rho'(124)
  b1qbv2=sqrt(2./(qbv2+fFUD)); // B1(Q rho' phi pi3), rho'(123)
  b1qbv3=sqrt(2./(qbv3+fFUD)); // B1(Q rho' phi pi4), rho'(124)
  for(int i=0;i<2;i++){
    w1m12[i]=0.0;
    w1m13[i]=0.0;
    for(int j1=0;j1<3;j1++){
      for(int j2=0;j2<3;j2++){
	for(int j3=0;j3<3;j3++){
	  for(int j4=0;j4<3;j4++){
	    tmp=E[i][j1][j2][4]*E[j2][j3][j4][4];
	    w1m12[i]=w1m12[i]+tmp*akv2m3w[j1]*ak12w[j3]*wt[j4]*b1qbv2; // The second term in Eq.(46)
	    w1m13[i]=w1m13[i]+tmp*akv3m2w[j1]*ak13w[j3]*wt[j4]*b1qbv3; // The first term in Eq.(46)
	  }
	}
      }
    }
  }

  // phi(1650)f0(980)  contribution // mass 1680 +- 20 MeV, width 150 +- 50 MeV, dominant decay KK*(892) 
  q2r23=0.25*s23-(_m[2])*(_m[2]); // Q^2_{abc}(f->pi+ pi-)
  b1q2r23=sqrt(2./(q2r23+fFUD));  // B1(Q_abc) = B1(Q f pi+ pi-)
  for(int i=0;i<4;i++){
    ak23wu[i]=ak23w[i]*fDel[i][i]*b1q2r23; // \tilde{t}_mu(34)
  }
  for(int i=0;i<2;i++){
    wpf22[i]=0.0;
    for(int j=0;j<4;j++){
      wpf22[i]=wpf22[i]-t2wvf[i][j]*ak23wu[j]; // -\tilde{T}^{(2)\mu\nu}(phi f0) \tilde{t}_nu(34)
    }
  }
}
