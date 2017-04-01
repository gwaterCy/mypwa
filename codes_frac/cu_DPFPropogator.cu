#include <math.h>
#include <cuda_runtime.h>
#include "cuComplex.h"
#include "cu_DPFPropogator.h"

const float rk=0.493677;
const float rp=0.139556995;

__device__ float2 cro(float sx, float am1, float am2)
  {
	           float2 ci = make_cuFloatComplex(0.0,1.0);
	           float t1=pow((am1+am2),2);
	           float t2=pow((am1-am2),2);
	           float st=(sx-t1)*(sx-t2);
	           float cro=sqrt(fabs(st))/sx;
		   float2 result = make_cuFloatComplex(cro,0.0);
	           if (st<0.) result=cuCmulfc(cro,ci);
	           return  result;
  }
__device__ float2 propogator980(float mass, float g11, float g22,float sx)
  {
	           float2 ci = make_cuFloatComplex(0.0,1.0);
	           float rm=mass*mass;
               //float2 propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
	           float2 propogator980=cuCdivfc(1.0,(cuCsubfc((rm-sx),cuCmulf(ci,cuCaddf( cuCmulfc(g11,cro(sx,rp,rp)),cuCmulfc(g22,cro(sx,rk,rk)) )) )) );
	           return propogator980;
  }
__device__ float2 pip(float sx)
  {
	           float2 ci = make_cuFloatComplex(0.0,1.0);
	           float xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
		   if(xk2<=0.)xk2=0.0;
	           float r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
	           return  make_cuFloatComplex(r4pip,0.0);
  }
__device__ float2 propogator600(float mass, float b1, float b2, float b3, float b4, float b5, float sx)
  {
	           float2 ci = make_cuFloatComplex(0.0,1.0);
      float am1=mass;
	           float as=am1*am1;
	           //float cgam1=(am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(float)(sx-0.0097)/(float)(as-0.0097)*(float)exp(-(sx-as)/b3)).real();
              float cgam1=cuCrealf(cuCmulcf(cuCdivf(cuCmulfc(am1*(b1+b2*sx),cro(sx,rp,rp)) , cro(as,rp,rp)) ,((float)(sx-0.0097)/(float)(as-0.0097)*(float)exp(-(sx-as)/b3))) );
	           float cgam2=cuCrealf(cuCdivf( cuCmulfc((am1*b4),pip(sx)) , pip(as) ));
	           //float2 propogator600=make_cuFloatComplex(1.0,0.0)/(as-sx-ci*b5*(cgam1+cgam2));
               float2 propogator600=cuCdivfc(1.0,cuCsubfc((as-sx),cuCmulcf(ci,(b5*(cgam1+cgam2))) ) );
	           return propogator600;
	    }

__device__ float2 propogator(float mass, float width, float sx)
{
	float2 ci = make_cuFloatComplex(0.0,1.0);
	float am=mass;
	float g1=mass*width;
    //float2 prop=g1/(sx-pow(am,2)+ci*g1);
	float2 prop=cuCdivfc(g1,cuCaddfc((sx-pow(am,2)),cuCmulcf(ci,g1)) );
	return prop;
}
__device__ float2 propogator1270(float mass, float width, float sx)
{
	float2 ci = make_cuFloatComplex(0.0,1.0);
	float rm=mass*mass;
	float gr=mass*width;
	float q2r=0.25*rm-0.0194792;
	float b2r=q2r*(q2r+0.1825)+0.033306;
	float g11270=gr*b2r/pow(q2r,2.5f);
	float q2=0.25*sx-0.0194792;
	float b2=q2*(q2+0.1825)+0.033306;
	float g1=g11270*pow(q2,2.5f)/b2;
	//float2 prop=gr/(sx-rm+ci*g1);
    float2 prop=cuCdivfc(gr,cuCaddfc( (sx-rm),cuCmulcf(ci,g1)) );
	return prop;
}
