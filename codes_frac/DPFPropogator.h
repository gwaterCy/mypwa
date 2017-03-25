#ifndef DPF_PROPOGATOR_HH
#define DPF_PROPOGATOR_HH

#include "conf.h"
#include <cuda_runtime.h>

  __device__ TComplex cro(my_float sx, my_float am1, my_float am2);
  __device__ TComplex propogator980(my_float mass, my_float g11, my_float g22,my_float sx);
  __device__ TComplex pip(my_float sx);
  __device__ TComplex propogator600(my_float mass, my_float b1, my_float b2, my_float b3, my_float b4, my_float b5, my_float sx);
  __device__ TComplex propogator(my_float mass,my_float width,my_float sx) ;
  __device__ TComplex propogator1270(my_float mass,my_float width,my_float sx) ;

#endif
