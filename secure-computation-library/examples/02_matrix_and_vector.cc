#include <iostream>

#include "scl/math.h"

int main() {
  /* This defines a finite field (as in the previous example) and then Matrices,
   * respectively vectors over said finite field.
   */
  using FF = scl::FF<32>;
  using Mat = scl::Mat<FF>;
  using Vec = scl::Vec<FF>;

  /* Like with field elements, it is possible to create random vectors and
   * matrices using a PRG. The first argument denotes the size for vectors, or
   * the row, column count for matrices.
   */
  scl::PRG prg;
  auto a = Vec::Random(5, prg);
  auto b = Vec::Random(5, prg);

  std::cout << "a: " << a << "\n";
  std::cout << "b: " << b << "\n";

  /* Addition, subtraction, etc can be performed on both vectors and matrices,
   * provided their sizes match.
   */
  auto c = a.Add(b);

  std::cout << "a + b = " << c << "\n";

  /* In-place operators are also supported.
   *
   * All operators (that return a vector/matrix of the same size) come in pairs
   * of "Op" and "OpInPlace".
   */
  auto d = a.AddInPlace(c);

  auto A = Mat::Random(3, 5, prg);

  std::cout << A << "\n";

  /* It is possible to coerce a vector into either a column or row matrix.
   */
  auto e = a.ToColumnMatrix();

  std::cout << "e: " << e << "\n";

  auto C = A.Multiply(e);
  std::cout << "A*e: " << C << "\n";
}
