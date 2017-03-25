#include <math.h>
#include <cuda_runtime.h>
#include "cuComplex.h"
#include "DPFPropogator.h"
#include "conf.h"
const double rk=0.493677;
const double rp=0.139556995;

__device__ TComplex cro(my_float sx, my_float am1, my_float am2)
  {
	           TComplex ci = make_complex(0.0,1.0);
	           my_float t1=pow((am1+am2),2);
	           my_float t2=pow((am1-am2),2);
	           my_float st=(sx-t1)*(sx-t2);
	           my_float cro=sqrt(fabs(st))/sx;
		   TComplex result = make_complex(cro,0.0);
	           if (st<0.) result=cuCmuldc(cro,ci);
	           return  result;
  }
__device__ TComplex propogator980(my_float mass, my_float g11, my_float g22,my_float sx)
  {
	           TComplex ci = make_complex(0.0,1.0);
	           my_float rm=mass*mass;
               //TComplex propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
	           TComplex propogator980=cuCdivdc(1.0,(cuCsubdc((rm-sx),cuCmul(ci,cuCadd( cuCmuldc(g11,cro(sx,rp,rp)),cuCmuldc(g22,cro(sx,rk,rk)) )) )) );
	           return propogator980;
  }
__device__ TComplex pip(my_float sx)
  {
	           TComplex ci = make_complex(0.0,1.0);
	           my_float xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
		   if(xk2<=0.)xk2=0.0;
	           my_float r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
	           return  make_complex(r4pip,0.0);
  }
__device__ TComplex propogator600(my_float mass, my_float b1, my_float b2, my_float b3, my_float b4, my_float b5, my_float sx)
  {
	           TComplex ci = make_complex(0.0,1.0);
      my_float am1=mass;
	           my_float as=am1*am1;
	           //my_float cgam1=(am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(my_float)(sx-0.0097)/(my_float)(as-0.0097)*(my_float)exp(-(sx-as)/b3)).real();
              my_float cgam1=cuCreal(cuCmulcd(cuCdiv(cuCmuldc(am1*(b1+b2*sx),cro(sx,rp,rp)) , cro(as,rp,rp)) ,((my_float)(sx-0.0097)/(my_float)(as-0.0097)*(my_float)exp(-(sx-as)/b3))) );
	           my_float cgam2=cuCreal(cuCdiv( cuCmuldc((am1*b4),pip(sx)) , pip(as) ));
	           //TComplex propogator600=make_complex(1.0,0.0)/(as-sx-ci*b5*(cgam1+cgam2));
               TComplex propogator600=cuCdivdc(1.0,cuCsubdc((as-sx),cuCmulcd(ci,(b5*(cgam1+cgam2))) ) );
	           return propogator600;
	    }

__device__ TComplex propogator(my_float mass, my_float width, my_float sx)
{
	TComplex ci = make_complex(0.0,1.0);
	my_float am=mass;
	my_float g1=mass*width;
    //TComplex prop=g1/(sx-pow(am,2)+ci*g1);
	TComplex prop=cuCdivdc(g1,cuCadddc((sx-pow(am,2)),cuCmulcd(ci,g1)) );
	return prop;
}
__device__ TComplex propogator1270(my_float mass, my_float width, my_float sx)
{
	TComplex ci = make_complex(0.0,1.0);
	my_float rm=mass*mass;
	my_float gr=mass*width;
	my_float q2r=0.25*rm-0.0194792;
	my_float b2r=q2r*(q2r+0.1825)+0.033306;
	my_float g11270=gr*b2r/pow(q2r,2.5);
	my_float q2=0.25*sx-0.0194792;
	my_float b2=q2*(q2+0.1825)+0.033306;
	my_float g1=g11270*pow(q2,2.5)/b2;
	//TComplex prop=gr/(sx-rm+ci*g1);
    TComplex prop=cuCdivdc(gr,cuCadddc( (sx-rm),cuCmulcd(ci,g1)) );
	return prop;
}
