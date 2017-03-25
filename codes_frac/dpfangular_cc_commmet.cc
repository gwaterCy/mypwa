for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
        if(i==j){
            if(i<3) fDel[i][j]=-1;
            else fDel[i][j]=1;
        }
        else fDel[i][j]=0;
    }
}
for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
        if(i==j){
            if(i<3) fGel[i][j]=-1;
            else fGel[i][j]=0;
        }
        else fGel[i][j]=0;
    }
}
cout<<"haha: "<< __LINE__ << endl;
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
                G1[I][J][K][L] =fGel[I][J]*fGel[K][L] + fGel[J][K]*fGel[I][L] + fGel[K][I]*fGel[J][L];
            }
        }
    }
}
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
}
