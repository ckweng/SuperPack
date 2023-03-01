/**
 * @file evpoly.h
 *
 * SCL --- Secure Computation Library
 * Copyright (C) 2022 Anders Dalskov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 * USA
 */
#ifndef _TP_EVPOLY_INTERNAL_H
#define _TP_EVPOLY_INTERNAL_H

#include "tp.h"

namespace tp {

  /**
   * @brief Polynomials over finite fields, in evaluation representation.
   * 
   * @details A polynomial of degree d can be represented as a vector of
   * length d+1 containing its coefficients, but alternatively, it can
   * also be represented by a vector of its evaluations at d+1
   * points. This is the representation we make use of here.
   */
  class EvPolynomialInternal {
  public:
    /**
     * @brief Construct a constant polynomial with constant term 0.
     */
    EvPolynomialInternal() : mX{0}, mY{FF(0)}  {};

    /**
     * @brief Construct a constant polynomial.
     * @param constant the constant term of the polynomial
     */
    EvPolynomialInternal(const FF& constant) : mX{0}, mY({constant}) {};

    /**
     * @brief Construct a polynomial from a list of x and y points.
     * @param x_points the set of evaluation points
     * @param y_points the set of evaluations
     */
    EvPolynomialInternal(const std::vector<std::size_t>& x_points, const Vec& y_points) {
      if ((x_points.size() != y_points.Size()))
        throw std::invalid_argument("number of evaluation points and evaluations do not match");
      if ((x_points.size() == 0))
        throw std::invalid_argument("empty set cannot be used for initialization");
      mX = x_points;
      mY = y_points;
    };

    /**
     * @brief Construct a polynomial from a list of y points and a
     * starting point x. The x points are taken consecutive after this one.
     * @param x_start the first evaluation point
     * @param y_points the set of evaluations
     */
    EvPolynomialInternal(const std::size_t& x_start, const Vec& y_points) {
      if ((y_points.Size() == 0))
        throw std::invalid_argument("empty set cannot be used for initialization");
      std::vector<std::size_t> x_points;
      x_points.reserve(y_points.Size());
      for (std::size_t i = 0; i < y_points.Size(); i++){
        x_points.emplace_back(x_start + i);
      };
      mX = x_points;
      mY = y_points;
    };

    /**
     * @brief Construct a polynomial from a list of y points. The
     * starting point is taken to be 0
     * @param y_points the set of evaluations
     */
    EvPolynomialInternal(const Vec& y_points) : EvPolynomialInternal(0, y_points){};

    // from [1, ..., n] to [-k+1, ..., 0]
    static std::vector<std::vector<std::vector<FF>>> mLagrangePolynomialPos;

    // from [-k+1, ..., degree-k+1] to [degree-k+2, ..., n]
    static std::vector<std::vector<std::vector<FF>>> mLagrangePolynomialNeg;

    static void PreprocessLagrangePolynomial(
        std::size_t k,
        std::vector<std::size_t> degree_set,
        std::size_t max_degree) {
      // the d-th matrix contains lagrange polynomials for degree d
      mLagrangePolynomialPos.resize(max_degree+1);
      for(auto degree : degree_set) {
        // the jth vector contains polys evaluated at j in [0, -1, ..., -j+1]
        mLagrangePolynomialPos[degree].resize(k);
        for(std::size_t j = 0; j < k; ++j) {
          // the lth value is the poly evaluated from l to j
          // l is actually in [1, ..., n+1]
          FF x(-1*j);
          mLagrangePolynomialPos[degree][j].resize(degree+1);
          for(std::size_t l = 0; l < degree+1; ++l) { // source_idx = l
            FF ell(1);
            FF xl(l);
            for(std::size_t m = 0; m < degree+1; ++m) {
              if (m == l) continue;
              auto xm = FF(m);
              ell *= (x - xm) / (xl - xm);
            }
            mLagrangePolynomialPos[degree][j][l] = ell;
          }
        }
      }

      // the d-th matrix contains lagrange polynomials for degree d
      mLagrangePolynomialNeg.resize(max_degree+1);
      for(auto degree : degree_set) {
        std::size_t start_idx = degree+2-k;
        std::size_t num_dest_points = max_degree-degree+k;
        // the jth vector contains polys evaluated at j in [degree+2-k, ..., n]
        mLagrangePolynomialNeg[degree].resize(num_dest_points); // assume max_degree = n-1
        for(std::size_t j = 0; j < num_dest_points; ++j) {
          // the lth value is the poly evaluated from l to j
          // l is actually in [1, ..., n+1]
          FF x(j+start_idx);
          mLagrangePolynomialNeg[degree][j].resize(degree+1);
          for(std::size_t l = -k+1; l < start_idx; ++l) { // source_idx = l
            FF ell(1);
            FF xl(l);
            for(std::size_t m = -k+1; m < start_idx; ++m) {
              if (m == l) continue;
              auto xm = FF(m);
              ell *= (x - xm) / (xl - xm);
            }
            mLagrangePolynomialNeg[degree][j][l+k-1] = ell;
          }
        }
      }
    }

    FF EvaluatePos(const std::size_t &x) {
      FF z;
      if(mLagrangePolynomialPos.size() != 0 &&
          mLagrangePolynomialPos[Degree()].size() != 0) {
        for(std::size_t j = 0; j < Degree()+1; ++j) {
          z += mY[j] * mLagrangePolynomialPos[Degree()][x][j];
        }
      } else {
        throw std::invalid_argument("evpoly: evaluate positive undefined preprocessing");
      }
      return z;
    }

    Vec EvaluatePos(const std::vector<std::size_t>& points) {
      Vec output;
      output.Reserve(points.size());
      for(const auto& point : points) {
        output.Emplace(EvaluatePos(point));
      }
      return output;
    }

    FF Evaluate(const std::size_t &x) {
      FF z;
      if(x > Degree()+1) {
        size_t xp = -x;
        if(mLagrangePolynomialPos.size() != 0 &&
            mLagrangePolynomialPos[Degree()].size() != 0) {
          for(std::size_t j = 0; j < Degree()+1; ++j) {
            z += mY[j] * mLagrangePolynomialPos[Degree()][xp][j];
          }
        } else {
          throw std::invalid_argument("evpoly: evaluate positive undefined preprocessing");
        }
      } else {
        if(mLagrangePolynomialNeg.size() != 0 &&
            mLagrangePolynomialNeg[Degree()].size() != 0) {
          for(std::size_t j = 0; j < Degree()+1; ++j) {
            z += mY[j] * mLagrangePolynomialPos[Degree()][x][j];
          }
        } else {
          throw std::invalid_argument("evpoly: evaluate negative undefined preprocessing");
        }
      }
      return z;
    }

    Vec Evaluate(const std::vector<std::size_t>& points) {
      Vec output;
      output.Reserve(points.size());
      for(const auto& point : points) {
        output.Emplace(Evaluate(point));
      }
      return output;
    }

    static FF EvaluateShift(const std::size_t &x, Vec &y_values, std::size_t degree) {
      FF z;
      if(x > degree+1) {
        size_t xp = -x;
        if(mLagrangePolynomialPos.size() != 0 &&
            mLagrangePolynomialPos[degree].size() != 0) {
          for(std::size_t j = 0; j < degree+1; ++j) {
            z += y_values[j] * mLagrangePolynomialPos[degree][xp][j];
          }
        } else {
          throw std::invalid_argument("evpoly: shift positive undefined preprocessing");
        }
      } else {
        if(mLagrangePolynomialNeg.size() != 0 &&
            mLagrangePolynomialNeg[degree].size() != 0) {
          for(std::size_t j = 0; j < degree+2; ++j) {
            z += y_values[j] * mLagrangePolynomialNeg[degree][x][j];
          }
        } else {
          throw std::invalid_argument("evpoly: shift negative undefined preprocessing");
        }
      }
      return z;
    }


    /**
     * @brief Evaluate this polynomial on a supplied point.
     * @param x the point to evaluate this polynomial on
     * @return f(x) where \p x is the supplied point and f this polynomial.
     */
    FF Evaluate(const FF& x) {
      FF z;
      for (std::size_t j = 0; j < Degree()+1; ++j) {
        FF ell(1);
        auto xj = FF(mX[j]);
        for (std::size_t m = 0; m < Degree()+1; ++m) {
          if (m == j) continue;
          auto xm = FF(mX[m]);
          ell *= (x - xm) / (xj - xm);
        }
        z += mY[j] * ell;
      }
      return z;
    };

    /**
     * @brief Evaluate this polynomial on a supplied list of points.
     * @param points points to evaluate on
     * @return f(x) for x in points
     */
    Vec Evaluate(const Vec& points) {
      Vec output;
      output.Reserve(points.Size());
      for(const auto& point : points) {
        output.Emplace(Evaluate(point));
      }
      return output;
    }

    Vec GetY() const { return mY; }
    std::vector<std::size_t> GetX() const { return mX; }

    /**
     * @brief Add two polynomials.
     */
    EvPolynomialInternal Add(const EvPolynomialInternal& q) {
      const auto c = mY.Add(q.GetY());
      return EvPolynomialInternal(mX, c);
    }

    /**
     * @brief Subtraction two polynomials.
     */
    EvPolynomialInternal Subtract(const EvPolynomialInternal& q) {
      const auto c = mY.Subtract(q.GetY());
      return EvPolynomialInternal(mX, c);
    }

    /**
     * @brief Returns true if this is the 0 polynomial.
     */
    bool IsZero() const { return Degree() == 0 && mY[0] == FF(0); };

    /**
     * @brief Degree of this polynomial.
     */
    std::size_t Degree() const { return mY.Size() - 1; };

  private:
    std::vector<std::size_t> mX;
    Vec mY;
  };

}  // namespace tp

#endif /* _TP_EVPOLY_INTERNAL_H */
