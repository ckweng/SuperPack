#include <iostream>
#include <stdexcept>

#include "scl/math.h"
#include "scl/secret_sharing.h"

int main() {
  using FF = scl::FF<32>;
  using Vec = scl::Vec<FF>;
  scl::PRG prg;

  /* We can easily create an additive secret sharing of some secret value:
   */
  FF secret(12345);
  Vec shares = scl::CreateAdditiveShares(secret, 5, prg);

  std::cout << "additive shares:\n" << shares << "\n";

  /* We can of course also reconstruct the secret again, given all 5 shares:
   */
  auto reconstructed = scl::ReconstructAdditive(shares);
  std::cout << "secret: " << reconstructed << "\n";

  /* SCL comes with three different ways of creating secret-shares for an honest
   * majority: Without any checks, with error detection, and with error
   * correction. Lets see error detection at work first
   */

  /* We create 4 shamir shares with a threshold of 1.
   */
  auto shamir_shares = scl::CreateShamirShares(secret, 4, 1, prg);
  std::cout << shamir_shares << "\n";

  /* Of course, these can be reconstructed. The second parameter is the
   * threshold. This performs reconstruction with error detection.
   */
  auto shamir_reconstructed = scl::ReconstructShamir(shamir_shares, 1);
  std::cout << shamir_reconstructed << "\n";

  /* If we introduce an error, then reconstruction fails
   */
  shamir_shares[2] = FF(123);
  try {
    std::cout << scl::ReconstructShamir(shamir_shares, 1) << "\n";
  } catch (std::logic_error& e) {
    std::cout << e.what() << "\n";
  }

  /* On the other hand, we can use the robust reconstruction since the threshold
   * is low enough. I.e., because 4 >= 3*1 + 1.
   */
  auto r = scl::ReconstructShamirRobust(shamir_shares, 1);
  std::cout << r << "\n";

  /* With a bit of extra work, we can even learn which share had the error.
   */

  /* first we need the alphas that were used when generating the shares. By
   * default these are just the field elements 1 through 4.
   */
  Vec alphas = {FF(1), FF(2), FF(3), FF(4)};
  auto pe = scl::ReconstructShamirRobust(shamir_shares, alphas, 1);

  /* pe is a pair of polynomials. The first is the original polynomial used for
   * generating the shares and the second is a polynomial whose roots tell which
   * share had errors.
   *
   * The secret is embedded in the constant term.
   */
  std::cout << pe[0].Evaluate(FF(0)) << "\n";

  /* This will be 0, indicating that the share corresponding to party 3 had an
   * error.
   */
  std::cout << pe[1].Evaluate(FF(3)) << "\n";

  /* Lastly, if there's too many errors, then correction is not possible
   */
  shamir_shares[1] = FF(22);
  try {
    scl::ReconstructShamirRobust(shamir_shares, 1);
  } catch (std::logic_error& e) {
    std::cout << e.what() << "\n";
  }
}
