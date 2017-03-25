#include <math.h>
#include <cuda_runtime.h>
#include "cuComplex.h"
#include "cu_DPFPropogator.h"

const double rk=0.493677;
const double rp=0.139556995;

__device__ double2 cro(double sx, double am1, double am2)
  {
	           double2 ci = make_cuDoubleComplex(0.0,1.0);
	           double t1=pow((am1+am2),2);
	           double t2=pow((am1-am2),2);
	           double st=(sx-t1)*(sx-t2);
	           double cro=sqrt(fabs(st))/sx;
		   double2 result = make_cuDoubleComplex(cro,0.0);
	           if (st<0.) result=cuCmuldc(cro,ci);
	           return  result;
  }
__device__ double2 propogator980(double mass, double g11, double g22,double sx)
  {
	           double2 ci = make_cuDoubleComplex(0.0,1.0);
	           double rm=mass*mass;
               //double2 propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
	           double2 propogator980=cuCdivdc(1.0,(cuCsubdc((rm-sx),cuCmul(ci,cuCadd( cuCmuldc(g11,cro(sx,rp,rp)),cuCmuldc(g22,cro(sx,rk,rk)) )) )) );
	           return propogator980;
  }
__device__ double2 pip(double sx)
  {
	           double2 ci = make_cuDoubleComplex(0.0,1.0);
	           double xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
		   if(xk2<=0.)xk2=0.0;
	           double r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
	           return  make_cuDoubleComplex(r4pip,0.0);
  }
__device__ double2 propogator600(double mass, double b1, double b2, double b3, double b4, double b5, double sx)
  {
	           double2 ci = make_cuDoubleComplex(0.0,1.0);
      double am1=mass;
	           double as=am1*am1;
	           //double cgam1=(am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(double)(sx-0.0097)/(double)(as-0.0097)*(double)exp(-(sx-as)/b3)).real();
              double cgam1=cuCreal(cuCmulcd(cuCdiv(cuCmuldc(am1*(b1+b2*sx),cro(sx,rp,rp)) , cro(as,rp,rp)) ,((double)(sx-0.0097)/(double)(as-0.0097)*(double)exp(-(sx-as)/b3))) );
	           double cgam2=cuCreal(cuCdiv( cuCmuldc((am1*b4),pip(sx)) , pip(as) ));
	           //double2 propogator600=make_cuDoubleComplex(1.0,0.0)/(as-sx-ci*b5*(cgam1+cgam2));
               double2 propogator600=cuCdivdc(1.0,cuCsubdc((as-sx),cuCmulcd(ci,(b5*(cgam1+cgam2))) ) );
	           return propogator600;
	    }

__device__ double2 propogator(double mass, double width, double sx)
{
	double2 ci = make_cuDoubleComplex(0.0,1.0);
	double am=mass;
	double g1=mass*width;
    //double2 prop=g1/(sx-pow(am,2)+ci*g1);
	double2 prop=cuCdivdc(g1,cuCadddc((sx-pow(am,2)),cuCmulcd(ci,g1)) );
	return prop;
}
__device__ double2 propogator1270(double mass, double width, double sx)
{
	double2 ci = make_cuDoubleComplex(0.0,1.0);
	double rm=mass*mass;
	double gr=mass*width;
	double q2r=0.25*rm-0.0194792;
	double b2r=q2r*(q2r+0.1825)+0.033306;
	double g11270=gr*b2r/pow(q2r,2.5);
	double q2=0.25*sx-0.0194792;
	double b2=q2*(q2+0.1825)+0.033306;
	double g1=g11270*pow(q2,2.5)/b2;
	//double2 prop=gr/(sx-rm+ci*g1);
    double2 prop=cuCdivdc(gr,cuCadddc( (sx-rm),cuCmulcd(ci,g1)) );
	return prop;
}
