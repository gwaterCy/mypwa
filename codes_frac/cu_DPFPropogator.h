#ifndef CU_DPF_PROPOGATOR_HH
#define CU_DPF_PROPOGATOR_HH

#include <cuda_runtime.h>

  __device__ float2 cro(float sx, float am1, float am2);
  __device__ float2 propogator980(float mass, float g11, float g22,float sx);
  __device__ float2 pip(float sx);
  __device__ float2 propogator600(float mass, float b1, float b2, float b3, float b4, float b5, float sx);
  __device__ float2 propogator(float mass,float width,float sx) ;
  __device__ float2 propogator1270(float mass,float width,float sx) ;

#endif
