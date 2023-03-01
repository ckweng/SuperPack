#ifndef SHARE_CONV_H
#define SHARE_CONV_H

#include "tp.h"
#include "evpoly_internal.h"

namespace tp {

  template<typename T>
  T PackedToAdditive(T y, std::size_t x, std::size_t dest, std::size_t degree) {
    return y *= EvPolynomialInternal::mLagrangePolynomialPos[degree][-1*dest][x];
    /*T ret = y;
    T t_x(x);
    T t_dest(dest);
    for(std::size_t i = 1; i < degree+2; ++i) {
      if(i == x) continue;
      T xm(i);
      ret *= (t_dest - xm) / (t_x - xm);
    }
    return ret;*/
  }

}

#endif
