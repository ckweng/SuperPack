/**
 * @file packed.h
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
#ifndef _TP_PACKED_INTERNAL_H
#define _TP_PACKED_INTERNAL_H

#include "evpoly_internal.h"

namespace tp {

  /**
   * @brief Creates an evpoly from a vector of secrets
   * @param secrets the vector of secrets
   * @param x_points the evaluation points. The first secret.Size()
   * points are used for the secrets
   * @param prg pseudorandom function for randomness
   * @return A polynomial
   */
  inline EvPolynomialInternal EvPolyFromSecretsAndPointsInternal(const Vec& secrets, const std::vector<std::size_t>& x_points, scl::PRG& prg) {
    if (secrets.Size() > x_points.size())
      throw std::invalid_argument("number of secrets is larger than number of x points");
    std::size_t n_secrets = secrets.Size();
    std::size_t n_points = x_points.size();

    // Populate y_points with random entries for indexes above n_secrets
    Vec y_points_ = Vec::Random(n_points, prg);
    for (std::size_t i = 0; i < n_secrets; i++) {
      y_points_[i] = secrets[i];
    }

    Vec y_points(n_points);
    for(std::size_t i = 0; i < n_points - n_secrets; ++i) {
      y_points[i] = y_points_[i+n_secrets];
    }
    for(std::size_t i = n_points - n_secrets; i < n_secrets; ++i) {
      y_points[i] = EvPolynomialInternal::EvaluateShift(i, y_points, n_points-1);
    }

    auto poly = EvPolynomialInternal(x_points, y_points);
    return poly;
  }

  /**
   * @brief Creates an evpoly from a vector of secrets, putting them
   * in positions [0,-1,...]
   * @param secrets the vector of secrets
   * @param degree the desired degree
   * @param prg pseudorandom function for randomness
   * @return A polynomial
   */
  inline EvPolynomialInternal EvPolyFromSecretsAndDegreeInternal(const Vec& secrets, std::size_t degree, scl::PRG& prg) {
    std::vector<std::size_t> x_points;
    x_points.reserve(degree+1);
    for (std::size_t i = 0; i < degree + 1; i++) { x_points.emplace_back(i+1); }

    return EvPolyFromSecretsAndPointsInternal(secrets, x_points, prg);
  }

  /**
   * @brief Returns the shares f(1)...f(n) of a polynomial f
   * @param poly the ev polynomial
   * @param n_shares how many shares to return
   * @return A vector of shares
   */
  inline Vec SharesFromEvPolyInternal(EvPolynomialInternal& poly, std::size_t n_shares) {
    if(n_shares == poly.Degree()+1) {
      return poly.GetY();
    } else {
      Vec y = poly.GetY();
      y.Reserve(n_shares);
      for(std::size_t i = poly.Degree()+2; i < n_shares; ++i) {
        y.Emplace(poly.Evaluate(FF(i)));
      }
      return y;
    }
  }

  /**
   * @brief Returns the polynomial f of degree n-1 given by the
   * shares f(1)...f(n)
   * @param shares the shares. The ev points are 1...n
   * @return An EvPolynomialInternal
   */
  inline EvPolynomialInternal EvPolyFromSharesInternal(const Vec& shares) {
    std::vector<std::size_t> x_points;
    x_points.reserve(shares.Size());
    for (std::size_t i = 0; i < shares.Size(); i++) { x_points.emplace_back(i+1); }
    return EvPolynomialInternal(x_points, shares);
  }

  /**
   * @brief Reconstructs the secrets given shares
   * @param x_secrets x points for the secrets
   * @param shares the shares. The ev points are 1...n
   * @return A vector with the secrets
   */
  inline Vec SecretsFromPointsAndSharesInternal(const std::vector<std::size_t>& x_secrets, const Vec& shares) {
    auto poly = EvPolyFromSharesInternal(shares);
    auto ret = poly.EvaluatePos(x_secrets);
    return ret;
  }

  /**
   * @brief Reconstructs a vector of secrets given shares. The
   * secrets are assumed to be in positions [0,-1,...-length]
   * @param shares the shares. The ev points are 1...n
   * @param length length of the secret
   * @return The secrets
   */
  inline Vec SecretsFromSharesAndLengthInternal(const Vec& shares, std::size_t length) {
    std::vector<std::size_t> x_secrets;
    x_secrets.reserve(length);
    for (std::size_t i = 0; i < length; i++) { x_secrets.emplace_back(i); }
    return SecretsFromPointsAndSharesInternal(x_secrets, shares);
  }

  /**
   * @brief Reconstructs a single secret given shares
   * @param x_secret x point for the secret
   * @param shares the shares. The ev points are 1...n
   * @return The secret
   */
  inline FF SecretFromPointAndSharesInternal(const std::size_t& x_secret, const Vec& shares) {
    auto poly = EvPolyFromSharesInternal(shares);
    return poly.Evaluate(x_secret);
  }

  /**
   * @brief Reconstructs a secret given shares. The
   * secret is assumed to be in position 0
   * @param shares the shares. The ev points are 1...n
   * @return The secret
   */
  inline FF SecretFromSharesInternal(const Vec& shares) {
    return SecretFromPointAndSharesInternal(std::size_t(0), shares);
  }



}  // namespace tp

#endif  // _TP_PACKED_INTERNAL_H
