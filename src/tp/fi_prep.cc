#include "tp/correlator.h"

namespace tp {

  void Correlator::OLEDummyPrep() {
    // shared PRG
    if(scl::PRG::BlockSize() != 16) {
      throw std::invalid_argument("The seeds are not 16 bytes\n");
    }
    unsigned char seed[16] = {0x00, 0x01, 0x02, 0x03,
      0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f};
    scl::PRG prg(seed);

    // dummy delta and its additive shares 
    FF delta = FF::Random(prg);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      FF x_point = FF(-i);
      auto poly_delta = scl::details::EvPolyFromSecretAndPointAndDegree(
          delta, x_point, mParties-2*mBatchSize+1, prg);
      Vec shares_delta = scl::details::SharesFromEvPoly(poly_delta, mParties);
      Gate::SetPackedShrDelta(shares_delta[mID]);
    }

    Vec parties_add_shr_delta = Vec::Random(mParties, prg); // all deltas
    auto sum = parties_add_shr_delta.Sum();
    parties_add_shr_delta[mParties-1] += (delta - sum);
    mAddShrDelta = parties_add_shr_delta[mID];

    // dummy packed shares of a, b
    Vec vec_a = Vec::Random(mParties, mPRG);
    Vec vec_b = Vec::Random(mParties, mPRG);
    Vec vec_c = vec_a.MultiplyEntryWise(vec_b);
    auto poly_a = scl::details::EvPolyFromSecretsAndDegree(vec_a, mParties-1, prg);
    auto poly_b = scl::details::EvPolyFromSecretsAndDegree(vec_b, mParties-1, prg);
    Vec packed_shr_a = scl::details::SharesFromEvPoly(poly_a, mParties);
    Vec packed_shr_b = scl::details::SharesFromEvPoly(poly_b, mParties);
    mDummyPackedShrA = packed_shr_a[mID];
    mDummyPackedShrB = packed_shr_b[mID];

    // dummy prep vole
    // pair-wise additive shares of a^i * delta^j and b^i * delta^j
    mDummyAddVoleShrA.Reserve(mParties);
    mDummyAddVoleShrB.Reserve(mParties);
    for(std::size_t i = 0; i < mID; ++i) {
      mDummyAddVoleShrA.Emplace(FF(0));
      mDummyAddVoleShrB.Emplace(FF(0));
    }
    // this can be used to compute additive shares of authenticated a, b
    for(std::size_t i = mID; i < mParties; ++i) {
      mDummyAddVoleShrA.Emplace(mDummyPackedShrA * parties_add_shr_delta[i]);
      mDummyAddVoleShrB.Emplace(mDummyPackedShrB * parties_add_shr_delta[i]);
    }

    // dummy prep ole
    // store cross terms of a * b
    // can be used to compute shares of <c_1>, ..., <c_k>
    // how? multiply with the lagrange polynomials
    for(std::size_t i = 0; i < mID; ++i) {
      mDummyAddOleShrAiBj.Emplace(FF(0));
      mDummyAddOleShrAjBi.Emplace(FF(0));
    }
    for(std::size_t i = mID; i < mParties; ++i) {
      mDummyAddOleShrAiBj.Emplace(mDummyPackedShrA * packed_shr_b[i]);
      mDummyAddOleShrAjBi.Emplace(mDummyPackedShrB * packed_shr_a[i]);
    }
  }

  // A sender party generate random Shamir shares
  // [r|1]_t, ..., [r|k]_t and distribute shares
  void Correlator::GenIndRandShr(std::size_t owner_id) {
    if(mID != owner_id) return;
    Vec rand_r = Vec::Random(mBatchSize, mPRG);
    std::vector<Vec> shares_mat_rand_r(mParties);
    std::size_t t = mParties-2*mBatchSize+1;
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      // x: [-i, 1, 2, ..., t]
      // y: [r_i, random values]
      Vec x_points(t+1);
      x_points[0] = FF(-1*i);
      for(std::size_t j = 1; j <= t; ++j) {
        x_points[j] = FF(j);
      }
      Vec y_points = Vec::Random(t+1, mPRG);
      y_points[0] = rand_r[i];
      auto poly_rand_r = scl::details::EvPolynomial(x_points, y_points);
      Vec shamir_shr_rand_r = scl::details::SharesFromEvPoly(
          poly_rand_r, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        shares_mat_rand_r[j].Emplace(shamir_shr_rand_r[j]);
      }
    }

    for(std::size_t i = 0; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shares_mat_rand_r[i]);
    }
  }

  // Generate shares of [r|1]_t, ..., [r|k]_t
  Vec Correlator::GenIndRandShr() {
    std::vector<FF> zeros(mBatchSize, FF(0));
    Vec shr_r(zeros);
    for(std::size_t i = 0; i < mParties; ++i) {
      GenIndRandShr(i);
      Vec shr_single;
      mNetwork->Party(i)->Recv(shr_single);
      shr_r = shr_r.Add(shr_single);
    }
    return shr_r;
  }

  // generate individual Shamir secret shares of MAC key Delta
  // [Delta|1]_t, ..., [Delta|k]_t
  void Correlator::MacKeyGen() {
    Vec shr_ind_r = GenIndRandShr();

    FF add_shr_r = PackedToAdditive(shr_ind_r[0], mID+1,
        0, mParties-2*mBatchSize+1);
    FF add_shr_delta_p_r = mAddShrDelta + add_shr_r;
    mNetwork->Party(0)->Send(add_shr_delta_p_r);

    if(mID == 0) {
      FF delta_p_r(0);
      FF add_shr_delta_p_r;
      for(std::size_t i = 0; i < mParties; ++i) {
        mNetwork->Party(i)->Recv(add_shr_delta_p_r);
        delta_p_r = delta_p_r + add_shr_delta_p_r;
      }
      std::vector<FF> all_delta_p_r_buf(mBatchSize, delta_p_r);
      Vec all_delta_p_r(all_delta_p_r_buf);
      auto poly_all_delta_p_r = EvPolyFromSecretsAndDegreeInternal(
          all_delta_p_r, 2*(mBatchSize-1), mPRG);
      Vec packed_shrs_all_delta_p_r = SharesFromEvPolyInternal(
          poly_all_delta_p_r, mParties);
      for(std::size_t i = 0; i < mParties; ++i) {
        mNetwork->Party(i)->Send(packed_shrs_all_delta_p_r[i]);
      }
    }

    FF packed_shr_all_delta_p_r;
    mNetwork->Party(0)->Recv(packed_shr_all_delta_p_r);

    for(std::size_t i = 0; i < mBatchSize; ++i) {
      mShamirShrDelta.Emplace(packed_shr_all_delta_p_r-shr_ind_r[i]);
    }
  }

  // generate triples from dummy VOLEs and OLEs
  // unauthenticated c
  // both a, b in high degree
  void Correlator::UnAuthTriplesHighDeg(std::vector<Vec> &add_shr_delta_a,
      std::vector<Vec> &add_shr_delta_b,
      std::vector<Vec> &add_shr_c, std::size_t num) {

    // precompute lagrange polynomials
    std::vector<Vec> lagrange_poly(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      lagrange_poly[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF res(1);
        FF t_x(i+1);
        FF t_dest(-1*j);
        for(std::size_t k = 1; k < mParties+1; ++k) {
          if(k == (i+1)) continue;
          FF xm(k);
          res *= (t_dest - xm) / (t_x - xm);
        }
        lagrange_poly[i].Emplace(res);
      }
    }

    add_shr_delta_a.resize(num);
    add_shr_delta_b.resize(num);
    add_shr_c.resize(num);

    for(std::size_t i = 0; i < num; ++i) {
      add_shr_delta_a[i].Reserve(mBatchSize);
      add_shr_delta_b[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF a(0), b(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          a += lagrange_poly[mID][j] * mDummyAddVoleShrA[k];
          b += lagrange_poly[mID][j] * mDummyAddVoleShrB[k];
        }
        add_shr_delta_a[i].Emplace(a);
        add_shr_delta_b[i].Emplace(b);
      }

      add_shr_c[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF c(0);
        for(std::size_t k = 0; k < mID; ++k) {
          FF lag = lagrange_poly[mID][j] * lagrange_poly[k][j];
          c += mDummyAddOleShrAiBj[k] * lag;
          c += mDummyAddOleShrAjBi[k] * lag;
        }
        c += mDummyPackedShrA * mDummyPackedShrB *
          lagrange_poly[mID][j] * lagrange_poly[mID][j];
        for(std::size_t k = mID+1; k < mParties; ++k) {
          FF lag = lagrange_poly[mID][j] * lagrange_poly[k][j];
          c += mDummyAddOleShrAiBj[k] * lag;
          c += mDummyAddOleShrAjBi[k] * lag;
        }
        add_shr_c[j].Emplace(c);
      }
    }
  }

  // A sender party generate random packed Shamir shares
  // ([r]_{n-k}, [r]_{n-1}) and distribute shares
  // TODO check if need to generate |C|/n pairs
  void Correlator::GenPairRandShr(std::size_t owner_id,
      std::size_t n_rand_pairs) {

    if(mID != owner_id) return;

    std::vector<Vec> shr_r_low_deg(mParties);
    std::vector<Vec> shr_r_high_deg(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      shr_r_low_deg[i].Reserve(n_rand_pairs);
      shr_r_high_deg[i].Reserve(n_rand_pairs);
    }

    for(std::size_t i = 0; i < n_rand_pairs; ++i) {
      Vec rand_r = Vec::Random(mBatchSize, mPRG);
      auto poly_r_low_deg = EvPolyFromSecretsAndDegreeInternal(
          rand_r, mParties-mBatchSize, mPRG);
      Vec shares_r_low_deg = SharesFromEvPolyInternal(
          poly_r_low_deg, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        shr_r_low_deg[j].Emplace(shares_r_low_deg[j]);
      }
      auto poly_r_high_deg = EvPolyFromSecretsAndDegreeInternal(
          rand_r, mParties-1, mPRG);
      Vec shares_r_high_deg = SharesFromEvPolyInternal(
          poly_r_high_deg, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        shr_r_high_deg[j].Emplace(shares_r_high_deg[j]);
      }
    }

    for(std::size_t i = 0; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shr_r_low_deg[i]);
      mNetwork->Party(i)->Send(shr_r_high_deg[i]);
    }
  }

  // Generate shares of ([r]_{n-k}, [r]_{n-1})
  void Correlator::GenPairRandShr(
      std::vector<std::pair<FF, FF>> &rand_pair,
      std::size_t n_rand_pairs) {

    std::vector<Vec> shr_low(mParties), shr_high(mParties);
    rand_pair.resize(n_rand_pairs*(2*mBatchSize-1));

    for(std::size_t i = 0; i < mParties; ++i) {
      GenPairRandShr(i, n_rand_pairs);
      mNetwork->Party(i)->Recv(shr_low[i]);
      mNetwork->Party(i)->Recv(shr_high[i]);
    }

    std::size_t cnt = 0;
    for(std::size_t i = 0; i < n_rand_pairs; ++i) {
      for(std::size_t j = 0; j < 2*mBatchSize-1; ++j) {
        FF a(0), b(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          a += mVandermonde[k][j] * shr_low[k][i];
          b += mVandermonde[k][j] * shr_low[k][i];
        }
        rand_pair[cnt++] = std::make_pair(a, b);
      }
    }
  }

  // reduce degree of a, b from n-1 to n-k
  // TODO check number of preprocessed rand pairs
  void Correlator::DegreeReduce(std::vector<FF> &packed_shr_a,
      std::vector<FF> &packed_shr_b, std::size_t num_triples) {
    // need num_triples pairs of ([r]_{n-k},[r]_{n-1})
    std::size_t n_m_t = 2 * mBatchSize - 1;
    std::size_t n_rand_pairs = (num_triples + n_m_t - 1) / n_m_t;
    std::vector<std::pair<FF, FF>> rand_pair;
    GenPairRandShr(rand_pair, n_rand_pairs);

    Vec share_masked_a(num_triples);
    Vec share_masked_b(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      share_masked_a[i] = mDummyPackedShrA + rand_pair[i].second;
      share_masked_b[i] = mDummyPackedShrA + rand_pair[i].second;
    }

    mNetwork->Party(0)->Send(share_masked_a);
    mNetwork->Party(0)->Send(share_masked_b);

    packed_shr_a.resize(num_triples);
    packed_shr_b.resize(num_triples);

    // P1 recover a+r, b+r from degree n-1 masked shares
    if(mID == 0) {
      std::vector<Vec> parties_masked_shares_a(mParties);
      std::vector<Vec> parties_masked_shares_b(mParties);
      for(std::size_t i = 0; i < mParties; ++i) {
        mNetwork->Party(i)->Recv(parties_masked_shares_a[i]);
        mNetwork->Party(i)->Recv(parties_masked_shares_b[i]);
      }

      std::vector<Vec> new_shares_mask_a_all(mParties);
      std::vector<Vec> new_shares_mask_b_all(mParties);
      for(std::size_t i = 0; i < mParties; ++i) {
        new_shares_mask_a_all[i].Reserve(num_triples);
        new_shares_mask_b_all[i].Reserve(num_triples);
      }

      for(std::size_t i = 0; i < num_triples; ++i) {
        Vec masked_shares_a(mParties), masked_shares_b(mParties);
        for(std::size_t j = 0; j < mParties; ++j) {
          masked_shares_a[j] = parties_masked_shares_a[j][i];
          masked_shares_b[j] = parties_masked_shares_b[j][i];
        }
        Vec vec_mask_a = SecretsFromSharesAndLengthInternal(
            masked_shares_a, mBatchSize);
        Vec vec_mask_b = SecretsFromSharesAndLengthInternal(
            masked_shares_b, mBatchSize);
        auto poly_mask_a = EvPolyFromSecretsAndDegreeInternal(
            vec_mask_a, 2*(mBatchSize-1), mPRG);
        auto poly_mask_b = EvPolyFromSecretsAndDegreeInternal(
            vec_mask_b, 2*(mBatchSize-1), mPRG);
        Vec new_shares_mask_a = SharesFromEvPolyInternal(
            poly_mask_a, mParties);
        Vec new_shares_mask_b = SharesFromEvPolyInternal(
            poly_mask_b, mParties);

        for(std::size_t j = 0; j < mParties; ++j) {
          new_shares_mask_a_all[j].Emplace(new_shares_mask_a[j]);
          new_shares_mask_b_all[j].Emplace(new_shares_mask_b[j]);
        }

        packed_shr_a[i] = new_shares_mask_a[0] - rand_pair[i].first;
        packed_shr_b[i] = new_shares_mask_b[0] - rand_pair[i].first;
      }

      for(std::size_t j = 1; j < mParties; ++j) {
        mNetwork->Party(j)->Send(new_shares_mask_a_all[j]);
        mNetwork->Party(j)->Send(new_shares_mask_b_all[j]);
      }
    } else {
      Vec a, b;
      mNetwork->Party(0)->Recv(a);
      mNetwork->Party(0)->Recv(b);
      for(std::size_t i = 0; i < num_triples; ++i) {
        packed_shr_a[i] = a[i] - rand_pair[i].first;
        packed_shr_b[i] = b[i] - rand_pair[i].first;
      }
    }
  }

  // A sender party generate random packed Shamir shares
  // ([r]_{n-k}) and distribute shares
  // TODO check if need to generate |C|/n pairs
  void Correlator::GenRandPackedShr(std::size_t owner_id,
      std::size_t degree, std::size_t num_shares) {
    if(mID != owner_id) return;

    std::vector<Vec> shares_r(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      shares_r[i].Reserve(num_shares);
    }

    for(std::size_t i = 0; i < num_shares; ++i) {
      Vec rand_r = Vec::Random(mBatchSize, mPRG);
      auto poly_r = EvPolyFromSecretsAndDegreeInternal(
          rand_r, degree, mPRG);
      Vec ind_shares_r = SharesFromEvPolyInternal(
          poly_r, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        shares_r[j].Emplace(ind_shares_r[j]);
      }
    }

    for(std::size_t i = 0; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shares_r[i]);
    }
  }

  // Generate shares of [r]_{n-k}
  void Correlator::GenRandPackedShr(Vec &shr_r,
      std::size_t degree, std::size_t num_shares) {
    shr_r.Reserve(num_shares*(2*mBatchSize-1));

    std::vector<Vec> ind_shr_r(mParties);

    for(std::size_t i = 0; i < mParties; ++i) {
      GenRandPackedShr(i, degree, num_shares);
      mNetwork->Party(i)->Recv(ind_shr_r[i]);
    }

    for(std::size_t i = 0; i < num_shares; ++i) {
      for(std::size_t j = 0; j < 2*mBatchSize-1; ++j) {
        FF a(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          a += mVandermonde[k][j] * ind_shr_r[k][i];
        }
        shr_r.Emplace(a);
      }
    }
  }
  
  void Correlator::AdditiveShareTransform(
      std::vector<FF> &packed_shr_delta_a,
      std::vector<FF> &packed_shr_delta_b,
      std::vector<Vec> &add_shr_delta_a,
      std::vector<Vec> &add_shr_delta_b,
      std::size_t num_triples) {

    // generate packed shares of random mask
    Vec shr_r;
    std::size_t n_m_t = 2 * mBatchSize - 1;
    std::size_t num_rand_shr = (num_triples + n_m_t - 1) / n_m_t;
    GenRandPackedShr(shr_r,
        mParties-mBatchSize, num_rand_shr);

    // convert packed shares to additive shares
    std::vector<Vec> add_shr_r(num_triples);
    for(std::size_t j = 0; j < num_triples; ++j) {
      add_shr_r[j].Reserve(mBatchSize);
      for(std::size_t i = 0; i < mBatchSize; ++i) {
        add_shr_r[j].Emplace(PackedToAdditive(shr_r[j], mID+1,
            -1*i, mParties-mBatchSize));
      }
    }

    // mask additive shares of delta_a and delta_b
    std::vector<Vec> add_shr_masked_a(num_triples);
    std::vector<Vec> add_shr_masked_b(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      add_shr_masked_a[i].Reserve(mBatchSize);
      add_shr_masked_b[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        add_shr_masked_a[i].Emplace(
            add_shr_delta_a[i][j] + add_shr_r[i][j]);
        add_shr_masked_b[i].Emplace(
            add_shr_delta_b[i][j] + add_shr_r[i][j]);
      }
    }

    // send masked shares to P1
    if(mID != 0) {
      for(std::size_t i = 0; i < mBatchSize; ++i) {
        Vec a(num_triples), b(num_triples);
        for(std::size_t j = 0; j < num_triples; ++j) {
          a[j] = add_shr_masked_a[j][i];
          b[j] = add_shr_masked_b[j][i];
        }
        mNetwork->Party(0)->Send(a);
        mNetwork->Party(0)->Send(b);
      }
    } else {
      for(std::size_t i = 1; i < mParties; ++i) {
        std::vector<Vec> recv_a(mBatchSize), recv_b(mBatchSize);
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          mNetwork->Party(i)->Recv(recv_a[j]);
          mNetwork->Party(i)->Recv(recv_b[j]);
        }
        for(std::size_t j = 0; j < num_triples; ++j) {
          Vec a(mBatchSize), b(mBatchSize);
          for(std::size_t k = 0; k < mBatchSize; ++k) {
            a[k] = recv_a[k][j];
            b[k] = recv_b[k][j];
          }
          add_shr_masked_a[j] = add_shr_masked_a[j].Add(a);
          add_shr_masked_b[j] = add_shr_masked_b[j].Add(b);
        }
      }
    }

    packed_shr_delta_a.resize(num_triples);
    packed_shr_delta_b.resize(num_triples);
    // P1 reconstruct shares
    if(mID == 0) {
      std::vector<Vec> shares_masked_a_all(mParties);
      std::vector<Vec> shares_masked_b_all(mParties);
      for(std::size_t i = 1; i < mParties; ++i) {
        shares_masked_a_all[i].Reserve(num_triples);
        shares_masked_b_all[i].Reserve(num_triples);
      }

      for(std::size_t i = 0; i < num_triples; ++i) {
        auto poly_masked_a = EvPolyFromSecretsAndDegreeInternal(
            add_shr_masked_a[i], 2*(mBatchSize-1), mPRG);
        Vec shares_masked_a = SharesFromEvPolyInternal(
            poly_masked_a, mParties);
        auto poly_masked_b = EvPolyFromSecretsAndDegreeInternal(
            add_shr_masked_b[i], 2*(mBatchSize-1), mPRG);
        Vec shares_masked_b = SharesFromEvPolyInternal(
            poly_masked_b, mParties);

        packed_shr_delta_a[i] = shares_masked_a[0] - shr_r[i];
        packed_shr_delta_b[i] = shares_masked_b[0] - shr_r[i];

        for(std::size_t j = 1; j < mParties; ++j) {
          shares_masked_a_all[j].Emplace(shares_masked_a[j]);
          shares_masked_b_all[j].Emplace(shares_masked_b[j]);
        }
      }

      for(std::size_t j = 1; j < mParties; ++j) {
        mNetwork->Party(j)->Send(shares_masked_a_all[j]);
        mNetwork->Party(j)->Send(shares_masked_b_all[j]);
      }
    } else {
      Vec a, b;
      mNetwork->Party(0)->Recv(a);
      mNetwork->Party(0)->Recv(b);
      for(std::size_t i = 0; i < num_triples; ++i) {
        packed_shr_delta_a[i] = a[i] - shr_r[i];
        packed_shr_delta_b[i] = b[i] - shr_r[i];
      }
    }
  }

  // generate random authenticated shares from dummy VOLEs
  // sacrifice to get authenticated c value in triples
  void Correlator::AuthTriplesC(
      std::vector<FF> &packed_shr_c,
      std::vector<Vec> &add_shr_delta_c,
      std::vector<Vec> &add_shr_c,
      std::size_t num_triples) {

    // precompute lagrange polynomials
    std::vector<Vec> lagrange_poly(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      lagrange_poly[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF res(1);
        FF t_x(i+1);
        FF t_dest(-1*j);
        for(std::size_t k = 1; k < mParties+1; ++k) {
          if(k == (i+1)) continue;
          FF xm(k);
          res *= (t_dest - xm) / (t_x - xm);
        }
        lagrange_poly[i].Emplace(res);
      }
    }

    // compute [r]_{n-1}, {<delta*r_i>}i\in[k]
    std::vector<FF> packed_shr_r(num_triples);
    std::vector<Vec> add_shr_delta_r(num_triples);

    for(std::size_t i = 0; i < num_triples; ++i) {
      add_shr_delta_r[i].Reserve(mBatchSize);
      packed_shr_r[i] = mDummyPackedShrA;
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF r(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          r += lagrange_poly[mID][j] * mDummyAddVoleShrA[k];
        }
        add_shr_delta_r[i][j] = r;
      }
    }

    // convert [r]_{n-1} to additive shares
    std::vector<Vec> add_shr_r(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      add_shr_r[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        add_shr_r[i].Emplace(PackedToAdditive(
              packed_shr_r[i], mID+1, -1*j, mParties-1));
      }
    }

    // 1. mask <c> and send masked values to P1
    //
    // 2. compute [Delta*(c_i+r_i)|i] = [Delta|i]_t * [c+r]_{2k-2}
    //                                = mShrmirShrDelta[i] * shares_masked_c
    std::vector<EvPolynomialInternal>
      polys_masked_c(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      Vec add_shr_masked_c(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        add_shr_masked_c[j] = add_shr_c[i][j] +
          add_shr_r[i][j];
      }
      if(mID == 0) {
        for(std::size_t j = 1; j < mParties; ++j) {
          Vec recv;
          mNetwork->Party(j)->Recv(recv);
          add_shr_masked_c = add_shr_masked_c.Add(recv);
        } 
        polys_masked_c[i] = EvPolyFromSecretsAndDegreeInternal(
            add_shr_masked_c, 2*(mBatchSize-1), mPRG);
      } else {
        mNetwork->Party(0)->Send(add_shr_masked_c);
      }
    }

    packed_shr_c.resize(num_triples);
    std::vector<Vec> packed_shr_delta_masked_c(num_triples);
    if(mID == 0) {
      std::vector<Vec> shares_masked_c_all(mParties);
      for(std::size_t i = 1; i < mParties; ++i) {
        shares_masked_c_all[i].Reserve(num_triples);
      }

      for(std::size_t i = 0; i < num_triples; ++i) {
        Vec shares_masked_c = SharesFromEvPolyInternal(
            polys_masked_c[i], mParties);
        for(std::size_t j = 1; j < mParties; ++j) {
          shares_masked_c_all[j].Emplace(shares_masked_c[j]);
          //mNetwork->Party(j)->Send(shares_masked_c[j]);
        }
        packed_shr_c[i] = shares_masked_c[0] - packed_shr_r[i];

        packed_shr_delta_masked_c[i].Reserve(mBatchSize);
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          packed_shr_delta_masked_c[i].Emplace(
              shares_masked_c[0] * mShamirShrDelta[j]);
        }
      }

      for(std::size_t j = 1; j < mParties; ++j) {
        mNetwork->Party(j)->Send(shares_masked_c_all[j]);
      }
    } else {
      Vec recv;
      mNetwork->Party(0)->Recv(recv);

      for(std::size_t i = 0; i < num_triples; ++i) {
        packed_shr_c[i] = recv[i] - packed_shr_r[i];

        packed_shr_delta_masked_c[i].Reserve(mBatchSize);
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          packed_shr_delta_masked_c[i].Emplace(
              recv[i] * mShamirShrDelta[j]);
        }
      }
    }

    // transform packed shamir share of authenticated masked c
    // and remove the mask
    add_shr_delta_c.resize(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      add_shr_delta_c[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF add_delta_masked_c = PackedToAdditive(
            packed_shr_delta_masked_c[i][j], mID+1, -1*j, mParties-1);
        add_shr_delta_c[i].Emplace(
            add_delta_masked_c - add_shr_delta_r[i][j]);
      }
    }
  }

  void Correlator::Sacrifice(
      FF mDummyPackedShrA, FF mDummyPackedShrB,
      std::vector<FF> packed_shr_a,
      std::vector<FF> packed_shr_b,
      std::vector<Vec> packed_shr_delta_a,
      std::vector<Vec> packed_shr_delta_b,
      std::vector<FF> packed_shr_c,
      std::vector<Vec> add_shr_delta_c,
      std::size_t num_triples) {

    // TODO random coin
    FF rho = FF(12345);

    // compute degree n-1 shares of
    // [\tilde{al}-rho*al] and [\tilde{bl}-bl]
    Vec packed_shr_tilde_al_m_rho_al(num_triples);
    Vec packed_shr_tilde_bl_m_bl(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      packed_shr_tilde_al_m_rho_al[i] =
        mDummyPackedShrA - rho * packed_shr_a[i];
      packed_shr_tilde_bl_m_bl[i] =
        mDummyPackedShrB - packed_shr_b[i];
    }

    // send above values to P1, who reconstructs and
    // distributes as degree k-1 shares
    if(mID == 0) {
      std::vector<Vec> packed_diff_a(mParties);
      std::vector<Vec> packed_diff_b(mParties);
      for(std::size_t i = 1; i < mParties; ++i) {
        mNetwork->Party(i)->Recv(packed_diff_a[i]);
        mNetwork->Party(i)->Recv(packed_diff_b[i]);
      }

      std::vector<Vec> new_shares_diff_a_all(mParties);
      std::vector<Vec> new_shares_diff_b_all(mParties);
      for(std::size_t i = 0; i < mParties; ++i) {
        new_shares_diff_a_all[i].Reserve(num_triples);
        new_shares_diff_b_all[i].Reserve(num_triples);
      }

      for(std::size_t i = 0; i < num_triples; ++i) {
        Vec shares_diff_a(mParties), shares_diff_b(mParties);
        shares_diff_a[0] = packed_shr_tilde_al_m_rho_al[i];
        shares_diff_b[0] = packed_shr_tilde_bl_m_bl[i];
        for(std::size_t j = 1; j < mParties; ++j) {
          shares_diff_a[j] = packed_diff_a[j][i];
          shares_diff_b[j] = packed_diff_b[j][i];
        }

        Vec vec_diff_a = SecretsFromSharesAndLengthInternal(
            shares_diff_a, mBatchSize);
        Vec vec_diff_b = SecretsFromSharesAndLengthInternal(
            shares_diff_b, mBatchSize);
        auto poly_diff_a = EvPolyFromSecretsAndDegreeInternal(
            vec_diff_a, mBatchSize-1, mPRG);
        auto poly_diff_b = EvPolyFromSecretsAndDegreeInternal(
            vec_diff_b, mBatchSize-1, mPRG);
        Vec new_shares_diff_a = SharesFromEvPolyInternal(
            poly_diff_a, mParties);
        Vec new_shares_diff_b = SharesFromEvPolyInternal(
            poly_diff_b, mParties);

        for(std::size_t j = 1; j < mParties; ++j) {
          new_shares_diff_a_all[j].Emplace(new_shares_diff_a[j]);
          new_shares_diff_b_all[j].Emplace(new_shares_diff_b[j]);
        }
        
        packed_shr_tilde_al_m_rho_al[0] = new_shares_diff_a[0];
        packed_shr_tilde_bl_m_bl[0] = new_shares_diff_b[0];
      }
      for(std::size_t i = 1; i < mParties; ++i) {
        mNetwork->Party(i)->Send(new_shares_diff_a_all[i]);
        mNetwork->Party(i)->Send(new_shares_diff_b_all[i]);
      }
    } else {
      mNetwork->Party(0)->Send(packed_shr_tilde_al_m_rho_al);
      mNetwork->Party(0)->Send(packed_shr_tilde_bl_m_bl);

      Vec a, b;
      mNetwork->Party(0)->Recv(a);
      mNetwork->Party(0)->Recv(b);
      for(std::size_t i = 0; i < num_triples; ++i) {
        packed_shr_tilde_al_m_rho_al[i] = a[i];
        packed_shr_tilde_bl_m_bl[i] = b[i];
      }
    }

    // main step of sacrifice: compute theta
    // it combines two triples in a way similar to
    // the use of Beaver triple
    Vec theta(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      FF a = packed_shr_tilde_al_m_rho_al[i];
      FF b = packed_shr_tilde_bl_m_bl[i];
      theta[i] = a * b + a * packed_shr_a[i] +
        rho * b * packed_shr_b[i] +
        rho * packed_shr_c[i] - packed_shr_c[i+num_triples];
    }

    // check if authentication part is correct
    // by sacrifice
    std::vector<Vec> delta_theta(num_triples);
    for(std::size_t i = 0; i < num_triples; ++i) {
      delta_theta[i].Reserve(mBatchSize);
      FF a = packed_shr_tilde_al_m_rho_al[i];
      FF b = packed_shr_tilde_bl_m_bl[i];
      FF t1 = a * b;
      FF t2 = rho * b ;
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF v1 = PackedToAdditive(mShamirShrDelta[j] * t1,
            mID+1, -1*j, mParties-1);
        FF v2 = PackedToAdditive(a * packed_shr_delta_b[i][j],
            mID+1, -1*j, mParties-1);
        FF v3 = PackedToAdditive(t2 * packed_shr_delta_a[i][j],
            mID+2, -1*j, mParties-1);
        delta_theta[i].Emplace(v1 + v2 + v3 +
            rho * add_shr_delta_c[i][j]
            - add_shr_delta_c[i+num_triples][j]);
      }
    }

    // TODO coin tossing
    unsigned char seed[16] = {0x00, 0x01, 0x02, 0x03,
      0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f};
    scl::PRG prg(seed);
    std::vector<FF> zeros(mBatchSize, FF(0));
    Vec chi = Vec::Random(num_triples, prg);
    FF packed_sum(0);
    Vec add_sum(zeros);
    for(std::size_t i = 0; i < num_triples; ++i) {
      packed_sum += chi[i] * theta[i];
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        add_sum[j] += chi[i] * delta_theta[i][j];
      }
    }

    // TODO commit and open, not sent to P1
    if(mID == 0) {
      for(std::size_t i = 1; i < mParties; ++i) {
        FF recv;
        Vec recv_vec;
        mNetwork->Party(i)->Recv(recv);
        packed_sum += recv;
        mNetwork->Party(i)->Recv(recv_vec);
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          add_sum[j] += recv_vec[j];
        }
      }

      // TODO check if all values are 0
    } else {
      mNetwork->Party(0)->Send(packed_sum);
      mNetwork->Party(0)->Send(add_sum);
    }


    // TODO degree check by VerifyDeg protocol
    // Not needed
  }

  void Correlator::Triple() {
    std::size_t num_triples = mNInOutBatches + mNMultBatches;

    // convert VOLE and OLE materials into triples
    // double number of triples, save another half for sacrifice check
    std::vector<Vec> add_shr_delta_a;
    std::vector<Vec> add_shr_delta_b;
    std::vector<Vec> add_shr_c; 
    UnAuthTriplesHighDeg(add_shr_delta_a, add_shr_delta_b, add_shr_c, 2*num_triples);

    // reduce degree of packed shares of a, b
    std::vector<FF> packed_shr_a;
    std::vector<FF> packed_shr_b;
    DegreeReduce(packed_shr_a, packed_shr_b, num_triples);

    // convert additive shares of authenticated a, b to packed shares
    std::vector<FF> packed_shr_delta_a;
    std::vector<FF> packed_shr_delta_b;
    AdditiveShareTransform(packed_shr_delta_a, packed_shr_delta_b,
        add_shr_delta_a, add_shr_delta_b, num_triples);

    // convert additive shares of c to packed shares and its authentication
    // double the number of output, save half for sacrifice check
    std::vector<FF> packed_shr_c;
    std::vector<Vec> add_shr_delta_c;
    AuthTriplesC(packed_shr_c, add_shr_delta_c,
      add_shr_c, 2*num_triples);

    Sacrifice(mDummyPackedShrA, mDummyPackedShrB,
        packed_shr_a, packed_shr_b,
        add_shr_delta_a, add_shr_delta_b,
        packed_shr_c, add_shr_delta_c,
        num_triples);

    for(std::size_t i = 0; i < mNInOutBatches; ++i) {
      IOBatchFIPrep prep;
      prep.mShrA = packed_shr_a[i];
      prep.mShrDeltaA = packed_shr_delta_a[i];
      prep.mShrB = packed_shr_b[i];
      prep.mShrDeltaB = packed_shr_delta_b[i];
      prep.mShrC = packed_shr_c[i];
      prep.mAddShrDeltaC = add_shr_delta_c[i];
      prep.mShrO1 = FF(0);
      mIOBatchFIPrep.emplace_back(prep);
    }

    for(std::size_t i = mNInOutBatches; i < num_triples; ++i) {
      MultBatchFIPrep prep;
      prep.mShrA = packed_shr_a[i];
      prep.mShrDeltaA = packed_shr_delta_a[i];
      prep.mShrB = packed_shr_b[i];
      prep.mShrDeltaB = packed_shr_delta_b[i];
      prep.mShrC = packed_shr_c[i];
      prep.mAddShrDeltaC = add_shr_delta_c[i];
      prep.mShrO1 = FF(0);
      prep.mShrO2 = FF(0);
      prep.mShrO3 = FF(0);
      mMultBatchFIPrep.emplace_back(prep);
    }
  }


  // Genrate [0]_{n-1} packed shares
  // A sender party generate random packed Shamir shares
  // [0]_{n-1} and distribute shares
  // TODO check how many time the function should be invoked 
  void Correlator::GenZeroShr(std::size_t owner_id,
      std::size_t num_shares) {
    if(mID != owner_id) return;

    std::vector<Vec> shares_zeros(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      shares_zeros[i].Reserve(num_shares);
    }
    std::vector<FF> zeros_in(mBatchSize, FF(0));
    Vec zeros(zeros_in);

    for(std::size_t i = 0; i < num_shares; ++i) {
      auto poly_zeros = EvPolyFromSecretsAndDegreeInternal(
          zeros, mParties-1, mPRG);
      Vec shares_zeros_ = SharesFromEvPolyInternal(
          poly_zeros, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        shares_zeros[j].Emplace(shares_zeros_[j]);
      }
    }

    for(std::size_t i = 0; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shares_zeros[i]);
    }
  }

  // Generate shares of [0]_{n-1}
  void Correlator::GenZeroShr(
      std::vector<FF> &zero_shares,
      std::size_t num_shares) {
    zero_shares.resize(num_shares*(2*mBatchSize-1));

    std::vector<Vec> all_shares(mParties);

    for(std::size_t i = 0; i < mParties; ++i) {
      GenZeroShr(i, num_shares);
      mNetwork->Party(i)->Recv(all_shares[i]);
    }

    std::size_t cnt = 0;
    for(std::size_t i = 0; i < num_shares; ++i) {
      for(std::size_t j = 0; j < 2*mBatchSize-1; ++j) {
        FF a(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          a += mVandermonde[k][j] * all_shares[k][i];
        }
        zero_shares[cnt++] = a;
      }
    }
  }

  // generate and fill in zero shares
  // TODO how many executions needed
  void Correlator::ZeroShr() {

    std::size_t n_m_t = 2 * mBatchSize - 1;
    std::size_t n_shares = (mNInOutBatches +
        3*mNMultBatches + n_m_t - 1) / n_m_t;
    std::vector<FF> zero_shares;
    GenZeroShr(zero_shares, n_shares);

    for(std::size_t i = 0; i < mNInOutBatches; ++i) {
      mIOBatchFIPrep[i].mShrO1 = zero_shares[i];
    }

    for(std::size_t i = 0, j = mNInOutBatches; i < mNMultBatches; ++i, ++j) {
      mMultBatchFIPrep[i].mShrO1 = zero_shares[j];
      mMultBatchFIPrep[i].mShrO2 = zero_shares[j];
      mMultBatchFIPrep[i].mShrO3 = zero_shares[j];
    }
  }

  // Genrate [rho*1]_{n-k} packed shares
  // A sender party generate random packed Shamir shares
  // [rho*1]_{n-1} and distribute shares
  // TODO check how many time the function should be invoked 
  void Correlator::GenUniShr(
      std::size_t owner_id,
      std::size_t num_shares) {
    if(mID != owner_id) return;

    std::vector<Vec> rho(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      rho[i].Reserve(num_shares);
    }

    for(std::size_t i = 0; i < num_shares; ++i) {
      FF rho_ = FF::Random(mPRG);
      std::vector<FF> rho_in(mBatchSize, rho_);
      Vec vec_rho(rho_in);
      auto poly_rho = EvPolyFromSecretsAndDegreeInternal(
          vec_rho, mParties-mBatchSize, mPRG);
      Vec shares_rho = SharesFromEvPolyInternal(
          poly_rho, mParties);
      for(std::size_t j = 0; j < mParties; ++j) {
        rho[j].Emplace(shares_rho[j]);
      }
    }

    for(std::size_t i = 0; i < mParties; ++i) {
      mNetwork->Party(i)->Send(rho[i]);
    }
  }

  // Generate shares of [rho]_{n-1}
  void Correlator::GenUniShr(
      std::vector<FF> &packed_shr_uni,
      std::size_t num_shares) {
    packed_shr_uni.resize(num_shares*(2*mBatchSize-1));

    std::vector<Vec> packed_shares_uni(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      GenUniShr(i, num_shares);
      mNetwork->Party(i)->Recv(packed_shares_uni[i]);
    }

    std::size_t cnt = 0;
    for(std::size_t i = 0; i < num_shares; ++i) {
      for(std::size_t j = 0; j < 2*mBatchSize-1; ++j) {
        FF a(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          a += packed_shares_uni[k][i] * mVandermonde[k][j];
        }
        packed_shr_uni[cnt++] = a;
      }
    }
  }


  // generate random authenticated values for individual wire
  // output [r*1]_{n-k}, <Delta*r>
  void Correlator::RandomAuthenticatedShare() {
    std::size_t num_wire = mNIndShrs;
    std::size_t num_batch = (num_wire + mBatchSize - 1) / mBatchSize;

    // precompute lagrange polynomials
    std::vector<Vec> lagrange_poly(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      lagrange_poly[i].Reserve(mBatchSize);
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF res(1);
        FF t_x(i+1);
        FF t_dest(-1*j);
        for(std::size_t k = 1; k < mParties+1; ++k) {
          if(k == (i+1)) continue;
          FF xm(k);
          res *= (t_dest - xm) / (t_x - xm);
        }
        lagrange_poly[i].Emplace(res);
      }
    }

    std::vector<FF> packed_shr_r(num_batch);

    // fill in [r]_{n-1} and <delta*r_i> from VOLE
    for(std::size_t i = 0; i < num_batch; ++i) {
      packed_shr_r[i] = mDummyPackedShrA;
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        FF r(0);
        for(std::size_t k = 0; k < mParties; ++k) {
          r += lagrange_poly[mID][j] * mDummyAddVoleShrA[k];
        }
        mAuthIndAddShrs.emplace_back(r);
      }
    }


    // generate [rho*1]_{n-k} shares via RandSh
    // total number of num_batch * k
    std::size_t n_rand_uni_shr = num_batch * mBatchSize;
    std::vector<FF> packed_shr_uni_rho;
    GenUniShr(packed_shr_uni_rho, (n_rand_uni_shr+2*mBatchSize-2)/(2*mBatchSize-1));

    // compute [r_l+rho]_{n-1} = [r_l]_{n-1} + \sum_{i=0}^k e_i * [rho*1]_{n-k}
    // reveal shares to P1
    std::size_t cnt = 0;
    Vec packed_shr_rl_p_rho(num_batch);
    for(std::size_t i = 0; i < num_batch; ++i) {
      packed_shr_rl_p_rho[i] = packed_shr_r[i];
      for(std::size_t j = 0; j < mBatchSize; ++j) {
        packed_shr_rl_p_rho[i] += packed_shr_uni_rho[cnt++] * mSharesOfEi[j];
      }
    }

    std::vector<Vec> packed_shares_rl_p_rho;
    if(mID == 0) {
      packed_shares_rl_p_rho.resize(num_batch);
      std::vector<Vec> recv(mParties);
      for(std::size_t i = 1; i < mParties; ++i) {
        mNetwork->Party(i)->Recv(recv[i]);
      }

      for(std::size_t i = 0; i < num_batch; ++i) {
        packed_shares_rl_p_rho[i].Reserve(mParties);
        packed_shares_rl_p_rho[i].Emplace(packed_shr_rl_p_rho[i]);
        for(std::size_t j = 1; j < mParties; ++j) {
          packed_shares_rl_p_rho[i].Emplace(recv[j][i]);
        }
      }
    } else {
      mNetwork->Party(0)->Send(packed_shr_rl_p_rho);
    }

    // P1 reconstruct r_l+rho, send cleartext to all parties
    // unmask
    std::vector<Vec> vec_rl_p_rho_all(mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      vec_rl_p_rho_all[i].Reserve(num_batch);
    }

    if(mID == 0) {
      for(std::size_t i = 0; i < num_batch; ++i) {
        Vec vec_rl_p_rho = SecretsFromSharesAndLengthInternal(
            packed_shares_rl_p_rho[i], mBatchSize);
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          vec_rl_p_rho_all[j].Emplace(vec_rl_p_rho[j]);
        }

        // unmask
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          std::vector<FF> cvec_uni_rli_p_rho(mBatchSize, vec_rl_p_rho[j]);
          Vec vec_uni_rli_p_rho(cvec_uni_rli_p_rho);
          auto poly_uni_rli_p_rho = EvPolyFromSecretsAndDegreeInternal(
              vec_uni_rli_p_rho, mParties-mBatchSize, mPRG);
          mIndShrs.emplace_back(poly_uni_rli_p_rho.Evaluate(FF(mID+1)));
        }
      }

      for(std::size_t j = 1; j < mParties; ++j) {
        for(std::size_t i = 0; i < mBatchSize; ++i) {
          mNetwork->Party(j)->Send(vec_rl_p_rho_all[i]);
        }
      }
    } else {

      for(std::size_t i = 0; i < mBatchSize; ++i) {
        mNetwork->Party(0)->Recv(vec_rl_p_rho_all[i]);
      }

      for(std::size_t i = 0; i < num_batch; ++i) {
        // unmask
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          std::vector<FF> cvec_uni_rli_p_rho(mBatchSize, vec_rl_p_rho_all[j][i]);
          Vec vec_uni_rli_p_rho(cvec_uni_rli_p_rho);
          auto poly_uni_rli_p_rho = EvPolyFromSecretsAndDegreeInternal(
              vec_uni_rli_p_rho, mParties-mBatchSize, mPRG);
          mIndShrs.emplace_back(poly_uni_rli_p_rho.Evaluate(FF(mID+1)));
        }
      }
    }

  } 

}
