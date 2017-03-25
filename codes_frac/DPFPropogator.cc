#include "DPFPropogator.h"
#include "TComplex.h"

//ClassImp(DPFPropogator);


const double rk=0.493677;
const double rp=0.139556995;
TComplex DPFPropogator::cro(Double_t sx, Double_t am1, Double_t am2)const
  {
	           TComplex ci(0,1);
//  cout<<"haha: "<< __LINE__ << endl;
	           Double_t t1=pow((am1+am2),2);
//		   cout<<"t1="<<t1<<endl;
	           Double_t t2=pow((am1-am2),2);
//		   cout<<"t2="<<t2<<endl;
//		   cout<<"sx="<<sx<<endl;
	           Double_t st=(sx-t1)*(sx-t2);
//		   cout<<"st="<<st<<endl;
	           Double_t cro=sqrt(fabs(st))/sx;
//		   cout<<"cro="<<cro<<endl;
		   TComplex result = cro;
	           if (st<0.) result=cro*ci;
	           return  result;
  }
TComplex DPFPropogator::propogator980(Double_t mass, Double_t g11, Double_t g22,Double_t sx)const
  {
//  cout<<"haha: "<< __LINE__ << endl;
	           TComplex ci(0,1);
	           Double_t rm=mass*mass;
//		   cout<<"rm="<<rm<<endl;
	           TComplex propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
//		   cout<<"propogator980="<<propogator980<<endl;
	           return propogator980;
  }
TComplex DPFPropogator::pip(Double_t sx)const
  {
//  cout<<"haha: "<< __LINE__ << endl;
	           TComplex ci(0,1);
	           Double_t xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
		   if(xk2<=0.)xk2=0.0;
	           Double_t r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
	           return  r4pip;
  }
TComplex DPFPropogator::propogator600(Double_t mass, Double_t b1, Double_t b2, Double_t b3, Double_t b4, Double_t b5, Double_t sx)const
  {
//  cout<<"haha: "<< __LINE__ << endl;
//std::cout << __FILE__ << __LINE__ << std::endl;
	           TComplex ci(0,1);
              	   Double_t am1=mass;
//		   cout<<"am1="<<am1<<endl;
	           Double_t as=am1*am1;
//		   cout<<"as="<<as<<endl;
	           Double_t cgam1=am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(sx-0.0097)/(as-0.0097)*exp(-(sx-as)/b3);
//		   cout<<"cgam1="<<cgam1<<endl;
	           Double_t cgam2=am1*b4*pip(sx)/pip(as);
//		   cout<<"cgam2="<<cgam2<<endl;
	           TComplex propogator600=1.0/(as-sx-ci*b5*(cgam1+cgam2));
  //      cout<<"propogator600="<<propogator600<<endl;
	           return propogator600;
	    }

TComplex DPFPropogator::propogator(Double_t mass, Double_t width, Double_t sx)const
{
//  cout<<"haha: "<< __LINE__ << endl;
//std::cout << __FILE__ << __LINE__ << std::endl;
	TComplex ci(0,1);
	Double_t am=mass;
//        cout<<"am="<<am<<endl;
	Double_t g1=mass*width;
	TComplex prop=g1/(sx-pow(am,2)+ci*g1);
//        cout<<"prop="<<prop<<endl;
	return prop;
}
TComplex DPFPropogator::propogator1270(Double_t mass, Double_t width, Double_t sx)const
{
//  cout<<"haha: "<< __LINE__ << endl;
//std::cout << __FILE__ << __LINE__ << std::endl;
	TComplex ci(0,1);
	Double_t rm=mass*mass;
	Double_t gr=mass*width;
	Double_t q2r=0.25*rm-0.0194792;
	Double_t b2r=q2r*(q2r+0.1825)+0.033306;
	Double_t g11270=gr*b2r/pow(q2r,2.5);
	Double_t q2=0.25*sx-0.0194792;
	Double_t b2=q2*(q2+0.1825)+0.033306;
	Double_t g1=g11270*pow(q2,2.5)/b2;
	TComplex prop=gr/(sx-rm+ci*g1);
//        cout<<"prop1270="<<prop<<endl;
	return prop;
}
