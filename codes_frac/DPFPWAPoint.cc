#include "DPFPWAPoint.h"
#include "TMath.h"

//ClassImp (DPFPWAPoint);
//using namespace DPFCoord;

DPFPWAPoint::DPFPWAPoint(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t M, TString phspfile, TString datafile)
{
    _M    = M;
    _m[0] = m1;
    _m[1] = m2;
    _m[2] = m3;
    _m[3] = m4;
    _M2    = M*M;
    _m2[0] = m1*m1;
    _m2[1] = m2*m2;
    _m2[2] = m3*m3;
    _m2[3] = m4*m4;
    _phspfile = phspfile;
    _datafile = datafile;
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
    _area  = 0;
    _sum2  = _M2 + _m2[0] + _m2[1] + _m2[2];
    _dalitzX  = DPFCoord::AB;   // default Dalitz plot
    _dalitzY  = DPFCoord::CA;
    _dalitzZ  = DPFCoord::BC;
    fm2Mmpi2  = _M2 - _m2[2];
    fm2Amds2  = _M2 + _m2[0];
    fm2Ampi2  = _M2 + _m2[2];
    mpi2Mmds2 = _m2[2] - _m2[0];
    D25Dfm2   = 0.25/_M2;
    xLower = lowLimit (_dalitzX);
    xUpper = highLimit(_dalitzX);
    xDiff  = xUpper - xLower;
    yLower = lowLimit (_dalitzY);
    yUpper = highLimit(_dalitzY);
    yDiff  = yUpper - yLower;
    cout<<"[DPFPWAPoint] ==> Creating DPFPWAPoint !" << endl;
}

DPFPWAPoint::DPFPWAPoint (const DPFPWAPoint &other)
{
    cout<<"[DPFPWAPoint] ==> Copying DPFPWAPoint !" << endl;
    _M     = other._M;
    _m[0]  = other._m[0];
    _m[1]  = other._m[1];
    _m[2]  = other._m[2];
    _M2    = other._M2;
    _m2[0] = other._m2[0];
    _m2[1] = other._m2[1];
    _m2[2] = other._m2[2];
    _area       = other._area;
    _sum2       = other._sum2;
    _dalitzX    = other._dalitzX;
    _dalitzY    = other._dalitzY;
    _dalitzZ    = other._dalitzZ;
    fm2Mmpi2   = other.fm2Mmpi2;
    fm2Amds2   = other.fm2Amds2;
    fm2Ampi2   = other.fm2Ampi2;
    mpi2Mmds2  = other.mpi2Mmds2;
    D25Dfm2    = other.D25Dfm2;
    xLower = other.xLower;
    xUpper = other.xUpper;
    xDiff  = other.xDiff;
    yLower = other.yLower;
    yUpper = other.yUpper;
    yDiff  = other.yDiff;
}

Double_t DPFPWAPoint::firstOfPair (DPFCoord::Pair pair) const
{return _m[firstMember(pair)];}

Double_t DPFPWAPoint::secondOfPair (DPFCoord::Pair pair) const
{return _m[secondMember(pair)];}

Double_t DPFPWAPoint::notOfPair (DPFCoord::Pair pair) const
{return _m[notAMember(pair)];}

Double_t DPFPWAPoint::lowLimit(DPFCoord::Pair pair) const
{
    Double_t low = firstOfPair(pair) + secondOfPair(pair);
    return low*low;
}
Double_t DPFPWAPoint::highLimit(DPFCoord::Pair pair) const
{
    Double_t smaller = firstOfPair(pair);
    smaller = smaller < secondOfPair(pair) ? smaller : secondOfPair(pair);
    Double_t high = _M - smaller;
    return high*high;
}

void DPFPWAPoint::calcArea()
{
    // returns the area of the Dalitz plot
    // using a very dumb numerical integral
    Double_t nSeg = 1e6;
    Double_t step = (xUpper-xLower)/nSeg;
    Double_t sum  = 0;
    Double_t ymax, ymin;
    for (Double_t dx=xLower+step/2.; dx<xUpper; dx+=step) {
        qMaxMin(dx,ymax, ymin);
        sum += (ymax-ymin);
    }
    sum *= step;
    _area = sum;
}

DPFCoord::Pair DPFPWAPoint::otherPair (const DPFCoord::Pair pair) const
{
    if      (pair==_dalitzX) return _dalitzY;
    else if (pair==_dalitzY) return _dalitzX;
    else {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

DPFCoord::Index DPFPWAPoint::commonMember (const DPFCoord::Pair pair1, const DPFCoord::Pair pair2) const
{
    if (pair1==pair2) {cout<<"[DPFPWAPoint] ==> pair1=pair2 ?"<<endl;return DPFCoord::A;}
    if      ((pair1==0&&pair2==1) || (pair1==1&&pair2==0)) {return DPFCoord::A;}
    else if ((pair1==0&&pair2==2) || (pair1==2&&pair2==0)) {return DPFCoord::B;}
    else if ((pair1==1&&pair2==2) || (pair1==2&&pair1==1)) {return DPFCoord::C;}
    else    {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

DPFCoord::Index DPFPWAPoint::notAMember (const DPFCoord::Pair pair) const
{
    if      (0==pair) {return DPFCoord::C;}
    else if (1==pair) {return DPFCoord::B;}
    else if (2==pair) {return DPFCoord::A;}
    else    {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

DPFCoord::Index DPFPWAPoint::notAMember (const DPFCoord::Index i1, const DPFCoord::Index i2) const
{
    if (DPFCoord::A!=i1 && DPFCoord::A!=i2) return DPFCoord::A;
    else if (DPFCoord::B!=i1 && DPFCoord::B!=i2) return DPFCoord::B;
    else if (DPFCoord::C!=i1 && DPFCoord::C!=i2) return DPFCoord::C;
    else    {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

DPFCoord::Index DPFPWAPoint::firstMember (const DPFCoord::Pair pair) const
{
    if      (0==pair) {return DPFCoord::A;}
    else if (1==pair) {return DPFCoord::C;}
    else if (2==pair) {return DPFCoord::B;}
    else    {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

DPFCoord::Index DPFPWAPoint::secondMember (const DPFCoord::Pair pair) const
{
    if      (0==pair) {return DPFCoord::B;}
    else if (1==pair) {return DPFCoord::A;}
    else if (2==pair) {return DPFCoord::C;}
    else    {cout<<"wrong combination...bailing out..."<<endl;assert(0);}
}

void DPFPWAPoint::qMaxMin(const Double_t q, Double_t &qmax, Double_t &qmin) const
{
    Double_t d2m  = 0.5/sqrt(q);
    Double_t e1CM = d2m*(q+mpi2Mmds2);
    Double_t e2CM = d2m*(fm2Mmpi2-q);

    Double_t e12CM2 = pow(e1CM+e2CM,2.);
    Double_t e1CM2  = e1CM*e1CM;
    Double_t e2CM2  = e2CM*e2CM;
    Double_t d1 = sqrt(e1CM2-_m2[2]);
    Double_t d2 = sqrt(e2CM2-_m2[2]);

    Double_t min = e12CM2-pow(d1+d2,2.);
    Double_t max = e12CM2-pow(d1-d2,2.);
    Double_t tmp = _sum2-q;
    qmax = tmp-min;
    qmin = tmp-max;
}

Bool_t DPFPWAPoint::valid(const Double_t dalX, const Double_t dalY) const
{
    Double_t dalYMin, dalYMax;
    qMaxMin(dalX, dalYMax, dalYMin);
    Bool_t val = ((dalY>dalYMin)&&(dalY<dalYMax)) ? kTRUE:kFALSE;
    return val;
}

void DPFPWAPoint::print()
{
    cout<<"[DPFPWAPoint] == > DPFPWAPoint defined with the following params:"<<endl;
    cout<<" M     : "<<_M<<endl;
    cout<<" m[0]  : "<<_m[0]<<endl;
    cout<<" m[1]  : "<<_m[1]<<endl;
    cout<<" m[2]  : "<<_m[2]<<endl;
    cout<<" M2    : "<<_M2<<endl;
    cout<<" m2[0] : "<<_m2[0]<<endl;
    cout<<" m2[1] : "<<_m2[1]<<endl;
    cout<<" m2[2] : "<<_m2[2]<<endl;
    cout<<" area  : "<<_area<<endl;
    cout<<" sum2  : "<<_sum2<<endl;
    cout<<" DalX Limit : " << xLower << " to " << xUpper << " = " << xDiff << endl;
    cout<<" DalY Limit : " << yLower << " to " << yUpper << " = " << yDiff << endl;
}

void DPFPWAPoint::calBack (Double_t dalX, Double_t dalY,
        Double_t dalQ[], Double_t dalSQ[], Double_t pPi2[], Double_t pPi1[],
        Double_t a2[], Double_t a1[]) const
{
    dalQ[0]  = dalX;
    dalQ[1]  = dalY;
    dalSQ[0] = sqrt(dalX);
    dalSQ[1] = sqrt(dalY);

    for (Int_t i=0;i<2;i++){
        pPi2[i] = p2Cal(dalQ[i]);
        pPi1[i] = p1Cal(dalQ[i]);
    }
    Double_t p0x   = p0Cal(dalX);
    Double_t p0y   = p0Cal(dalY);

    Double_t p1p2  = 0.5*(fm2Amds2-dalX-dalY);
    Double_t p1q2  = 0.5*(fm2Mmpi2-dalX);
    Double_t p2q2  = 0.5*(dalX+mpi2Mmds2);
    Double_t p1q2a = 0.5*(fm2Mmpi2-dalY);
    Double_t p2q2a = 0.5*(dalY+mpi2Mmds2);
    Double_t dalXTmpi2 = dalX*_m2[2];
    Double_t dalYTmpi2 = dalY*_m2[2];

    Double_t prodX = -2.0*p0x*pPi2[0];
    Double_t prodY = -2.0*p0y*pPi2[1];


    Double_t cthX = (p1q2*p2q2-p1p2*dalX)  /sqrt((p1q2*p1q2  -dalXTmpi2)*(p2q2*p2q2  -dalXTmpi2));
    Double_t cthY = (p1q2a*p2q2a-p1p2*dalY)/sqrt((p1q2a*p1q2a-dalYTmpi2)*(p2q2a*p2q2a-dalYTmpi2));

    a1[0] = prodX*cthX;
    a1[1] = prodY*cthY;
    a2[0] = prodX*prodX*(cthX*cthX-1.0/3.0);
    a2[1] = prodY*prodY*(cthY*cthY-1.0/3.0);

    return;
}

Double_t DPFPWAPoint::p2Cal(Double_t sx) const
{
    Double_t x = sx+mpi2Mmds2;
    Double_t t = 0.25*x*x/sx-_m2[2];
    Double_t p;
    if (t>0.0) {
        p = sqrt(t);
    } else {
        p = 0.04;
    }
    return p;
}

Double_t DPFPWAPoint::p1Cal(Double_t sx) const
{
    Double_t x = fm2Ampi2-sx;
    Double_t t = D25Dfm2*x*x-_m2[2];
    Double_t p;
    if (t>0.0) {
        p = sqrt(t);                 
    } else {
        p = 0.04;
    }
    return p;
}

Double_t DPFPWAPoint::p0Cal(Double_t sx) const
{
    Double_t x = fm2Mmpi2-sx;
    Double_t t = 0.25*x*x/sx-_m2[2];
    Double_t p;
    if (t>0.0) {
        p = sqrt(t);                 
    } else {
        cout << " hello, less than 0.0" << endl;
        p = 0.04;
    }
    return p;
}

Double_t DPFPWAPoint::calCos(Double_t x, Double_t y) const
{
    Double_t Ymax,Ymin;
    qMaxMin(x,Ymax,Ymin);
    if (y<Ymin||y>Ymax) return -999.9;
    Double_t coshel = (2.0*y-Ymax-Ymin)/(Ymax-Ymin);
    return coshel;
}

Double_t DPFPWAPoint::calY(Double_t x, Double_t coshel) const
{
    Double_t dx = sqrt(x);
    Double_t eD = (27.8494170251+x)/10.5582;
    Double_t ptmp = eD*eD-x;
    if (ptmp<0) return 0.0;
    Double_t p1 = sqrt(ptmp);
    Double_t e2 = 0.5*(x-3.4748027051)/dx;
    ptmp = e2*e2-0.0194797849;
    if (ptmp<0) return 0.0;
    Double_t p2 = sqrt(ptmp);
    Double_t w  = 10.5582*p1*p2/dx;
    Double_t w1 = 10.5582*eD*e2/dx;
    Double_t y  = 27.8883765949-(w1-coshel*w);
    return y;
}
