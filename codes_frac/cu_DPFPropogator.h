#ifndef CU_DPF_PROPOGATOR_HH
#define CU_DPF_PROPOGATOR_HH

#include <cuda_runtime.h>

  __device__ double2 cro(double sx, double am1, double am2);
  __device__ double2 propogator980(double mass, double g11, double g22,double sx);
  __device__ double2 pip(double sx);
  __device__ double2 propogator600(double mass, double b1, double b2, double b3, double b4, double b5, double sx);
  __device__ double2 propogator(double mass,double width,double sx) ;
  __device__ double2 propogator1270(double mass,double width,double sx) ;

#endif
