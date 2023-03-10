/**
 * @file vec.h
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
#ifndef _SCL_MATH_VEC_H
#define _SCL_MATH_VEC_H

#include <cstring>
#include <functional>
#include <sstream>
#include <type_traits>
#include <vector>

#include "scl/math/mat.h"
#include "scl/prg.h"

namespace scl {

/**
 * @brief Vector.
 *
 * This class is a thin wrapper around std::vector meant only to provide some
 * functionality that makes it behave like other classes present in SCUtil.
 */
template <typename T>
class Vec {
 public:
  /**
   * @brief The type of vector elements.
   */
  using ValueType = T;

  /**
   * @brief Iterator type.
   */
  using iterator = typename std::vector<T>::iterator;

  /**
   * @brief Const iterator type.
   */
  using const_iterator = typename std::vector<T>::const_iterator;

  /**
   * @brief Reverse iterator type.
   */
  using reverse_iterator = typename std::vector<T>::reverse_iterator;

  /**
   * @brief Reverse const iterator type.
   */
  using const_reverse_iterator =
      typename std::vector<T>::const_reverse_iterator;

  /**
   * @brief Read a vec from a stream of bytes.
   * @param n the number of elements to read
   * @param src the buffer to read from
   */
  static Vec<T> Read(std::size_t n, const unsigned char* src);

  /**
   * @brief Create a partially random Vec.
   *
   * This method creates a Vec object where entries satisfying a predicate are
   * random. Entries for which the predicate is false are default initialized.
   *
   * The \p predicate argument should be a callable that accepts an std::size_t
   * value and outputs something which can be treated as a bool.
   *
   * @param n the size of the Vec
   * @param predicate a predicate indicating what entries to randomize
   * @param prg a prg for creating random entries
   */
  template <typename Pred>
  static Vec<T> PartialRandom(std::size_t n, Pred predicate, PRG& prg);

  /**
   * @brief Create a Vec and populate it with random elements.
   * @param n the size of the vector
   * @param prg a PRG used to generate random elements
   * @return a Vec with random elements.
   */
  static Vec<T> Random(std::size_t n, PRG& prg);

  /**
   * @brief Default constructor that creates an empty Vec.
   */
  Vec(){};

  /**
   * @brief Construct a new Vec of explicit size.
   * @param n the size
   */
  explicit Vec(std::size_t n) : mValues(n){};

  /**
   * @brief Construct a vector from an initializer_list.
   * @brief values an initializer_list.
   */
  Vec(std::initializer_list<T> values) : mValues(values){};

  /**
   * @brief Construct a vector from an STL vector.
   * @param values an STL vector
   */
  explicit Vec(const std::vector<T>& values) : mValues(values){};

  /**
   * @brief Construct a Vec from a pair of iterators.
   * @param first iterator pointing to the first element
   * @param last iterator pointing to the one past last element
   * @tparam It iterator type
   */
  template <typename It>
  explicit Vec(It first, It last) : mValues(first, last) {}

  /**
   * @brief The size of the Vec.
   */
  std::size_t Size() const { return mValues.size(); };

  /**
   * @brief Reserve space for vector
   */
  void Reserve(std::size_t size) { mValues.reserve(size); };

  /**
   * @brief Appends element at the end of the vector
   */
  void Emplace(T value) { mValues.emplace_back(value); };  

  /**
   * @brief Mutable access to vector elements.
   */
  T& operator[](std::size_t idx) { return mValues[idx]; };

  /**
   * @brief Read only access to vector elements.
   */
  T operator[](std::size_t idx) const { return mValues[idx]; };

  /**
   * @brief Add two Vec objects entry-wise.
   * @param other the other vector
   * @return the sum of this and \p other.
   */
  Vec Add(const Vec& other) const;

  /**
   * @brief Add two Vec objects entry-wise in-place.
   * @param other the other vector
   * @return the sum of this and \p other, assigned to this.
   */
  Vec& AddInPlace(const Vec& other) {
    EnsureCompatible(other);
    for (std::size_t i = 0; i < Size(); i++) {
      mValues[i] += other.mValues[i];
    }
    return *this;
  };

  /**
   * @brief Subtract two Vec objects entry-wise.
   * @param other the other vector
   * @return the difference of this and \p other.
   */
  Vec Subtract(const Vec& other) const;

  /**
   * @brief Subtract two Vec objects entry-wise in-place.
   * @param other the other vector
   * @return the difference of this and \p other, assigned to this.
   */
  Vec& SubtractInPlace(const Vec& other) {
    EnsureCompatible(other);
    for (std::size_t i = 0; i < Size(); i++) {
      mValues[i] -= other.mValues[i];
    }
    return *this;
  };

  /**
   * @brief Multiply two Vec objects entry-wise.
   * @param other the other vector
   * @return the product of this and \p other.
   */
  Vec MultiplyEntryWise(const Vec& other) const;

  /**
   * @brief Multiply two Vec objects entry-wise in-place.
   * @param other the other vector
   * @return the product of this and \p other, assigned to this.
   */
  Vec& MultiplyEntryWiseInPlace(const Vec& other) {
    EnsureCompatible(other);
    for (std::size_t i = 0; i < Size(); i++) {
      mValues[i] *= other.mValues[i];
    }
    return *this;
  };

  /**
   * @brief Compute a dot product between this and another vector.
   * @param other the other vector
   * @return the dot (or inner) product of this and \p other.
   */
  T Dot(const Vec& other) const {
    EnsureCompatible(other);
    T result;
    for (std::size_t i = 0; i < Size(); i++)
      result += mValues[i] * other.mValues[i];
    return result;
  };

  /**
   * @brief Compute the sum over entries of this vector.
   * @return the sum of the entries of this vector.
   */
  T Sum() const {
    T sum;
    for (const auto& v : mValues) sum += v;
    return sum;
  };

  /**
   * @brief Scale this vector by a constant.
   * @param scalar the scalar
   * @return a scaled version of this vector.
   */
  Vec ScalarMultiply(const T& scalar) const {
    std::vector<T> r;
    r.reserve(Size());
    for (const auto& v : mValues) r.emplace_back(scalar * v);
    return Vec(r);
  };

  /**
   * @brief Scale this vector in-place by a constant.
   * @param scalar the scalar
   * @return a scaled version of this vector.
   */
  Vec& ScalarMultiplyInPlace(const T& scalar) {
    for (auto& v : mValues) v *= scalar;
    return *this;
  };

  /**
   * @brief Test if this vector is equal to another vector in constant time.
   * @param other the other vector
   * @return true if this vector is equal to \p other and false otherwise.
   */
  bool Equals(const Vec& other) const;

  /**
   * @brief Convert this vector into a 1-by-N row matrix.
   */
  Mat<T> ToRowMatrix() const { return Mat<T>{1, Size(), mValues}; };

  /**
   * @brief Convert this vector into a N-by-1 column matrix.
   */
  Mat<T> ToColumnMatrix() const { return Mat<T>{Size(), 1, mValues}; };

  /**
   * @brief Convert this Vec object to an std::vector.
   */
  std::vector<T>& ToStlVector() { return mValues; };

  /**
   * @brief Convert this Vec object to a const std::vector.
   */
  const std::vector<T>& ToStlVector() const { return mValues; };

  /**
   * @brief Return a string representation of this vector.
   */
  std::string ToString() const;

  /**
   * @brief Write a string representation of this vector to a stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const Vec<T>& v) {
    return os << v.ToString();
  };

  /**
   * @brief Write this Vec to a buffer.
   * @param dest the buffer to write this Vec to
   */
  void Write(unsigned char* dest) const;

  /**
   * @brief Returns the number of bytes that Write writes.
   */
  std::size_t ByteSize() const { return Size() * T::ByteSize(); };

  /**
   * @brief Return an iterator pointing to the start of this Vec.
   */
  iterator begin() { return mValues.begin(); };

  /**
   * @brief Provides a const iterator to the start of this Vec.
   */
  const_iterator begin() const { return mValues.begin(); };

  /**
   * @brief Provides a const iterator to the start of this Vec.
   */
  const_iterator cbegin() const { return mValues.cbegin(); };

  /**
   * @brief Provides an iterator pointing to the end of this Vec.
   */
  iterator end() { return mValues.end(); };

  /**
   * @brief Provides a const iterator pointing to the end of this Vec.
   */
  const_iterator end() const { return mValues.end(); };

  /**
   * @brief Provides a const iterator pointing to the end of this Vec.
   */
  const_iterator cend() const { return mValues.cend(); };

  /**
   * @brief Provides a reverse iterator pointing to the end of this Vec.
   */
  reverse_iterator rbegin() { return mValues.rbegin(); };

  /**
   * @brief Provides a reverse const iterator pointing to the end of this Vec.
   */
  const_reverse_iterator rbegin() const { return mValues.rbegin(); };

  /**
   * @brief Provides a reverse const iterator pointing to the end of this Vec.
   */
  const_reverse_iterator crbegin() const { return mValues.crbegin(); };

  /**
   * @brief Provides a reverse iterator pointing to the start of this Vec.
   */
  reverse_iterator rend() { return mValues.rend(); };

  /**
   * @brief Provides a reverse const iterator pointing to the start of this Vec.
   */
  const_reverse_iterator rend() const { return mValues.rend(); };

  /**
   * @brief Provides a reverse const iterator pointing to the start of this Vec.
   */
  const_reverse_iterator crend() const { return mValues.crend(); };

 private:
  void EnsureCompatible(const Vec& other) const {
    if (Size() != other.Size())
      throw std::invalid_argument("Vec sizes mismatch");
  };

  std::vector<T> mValues;
};

template <typename T>
Vec<T> Vec<T>::Read(std::size_t n, const unsigned char* src) {
  std::vector<T> r;
  r.reserve(n);
  std::size_t offset = 0;
  while (n-- > 0) {
    r.emplace_back(T::Read(src + offset));
    offset += T::ByteSize();
  }
  return Vec<T>(r);
}

template <typename T>
template <typename Pred>
Vec<T> Vec<T>::PartialRandom(std::size_t n, Pred predicate, PRG& prg) {
  std::vector<T> v;
  v.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    if (predicate(i)) {
      v.emplace_back(T::Random(prg));
    } else {
      v.emplace_back(T{});
    }
  }
  return Vec<T>(v);
}

template <typename T>
Vec<T> Vec<T>::Random(std::size_t n, PRG& prg) {
  return Vec<T>::PartialRandom(
      n,
      [](std::size_t i) {
        (void)i;
        return true;
      },
      prg);
}

template <typename T>
Vec<T> Vec<T>::Add(const Vec<T>& other) const {
  EnsureCompatible(other);
  std::vector<T> r;
  auto n = Size();
  r.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    r.emplace_back(mValues[i] + other.mValues[i]);
  }
  return Vec(r);
}

template <typename T>
Vec<T> Vec<T>::Subtract(const Vec<T>& other) const {
  EnsureCompatible(other);
  std::vector<T> r;
  auto n = Size();
  r.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    r.emplace_back(mValues[i] - other.mValues[i]);
  }
  return Vec(r);
}

template <typename T>
Vec<T> Vec<T>::MultiplyEntryWise(const Vec<T>& other) const {
  EnsureCompatible(other);
  std::vector<T> r;
  auto n = Size();
  r.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    r.emplace_back(mValues[i] * other.mValues[i]);
  }
  return Vec(r);
}

template <typename T>
bool Vec<T>::Equals(const Vec<T>& other) const {
  if (Size() != other.Size()) {
    return false;
  }

  bool equal = true;
  for (std::size_t i = 0; i < Size(); i++) {
    equal &= mValues[i] == other.mValues[i];
  }

  return equal;
}

template <typename T>
std::string Vec<T>::ToString() const {
  if (Size()) {
    std::stringstream ss;
    ss << "[";
    std::size_t i = 0;
    for (; i < Size() - 1; i++) ss << mValues[i] << ", ";
    ss << mValues[i] << "]";
    return ss.str();
  } else {
    return "[ EMPTY_VECTOR ]";
  }
}

template <typename T>
void Vec<T>::Write(unsigned char* dest) const {
  unsigned char* p = dest;
  for (const auto& v : mValues) {
    v.Write(p);
    p += T::ByteSize();
  }
}

}  // namespace scl

#endif  // _SCL_MATH_VEC_H
