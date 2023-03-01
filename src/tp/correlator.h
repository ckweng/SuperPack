#ifndef CORRELATOR_H
#define CORRELATOR_H

#include <catch2/catch.hpp>
#include <iostream>
#include <map>
#include <assert.h>

#include "tp.h"
#include "mult_gate.h"
#include "input_gate.h"
#include "output_gate.h"
#include "util.h"

namespace tp {

  // Preprocessing materials for input/output gates
  struct IOBatchFIPrep {
    // Triples
    FF mShrA;
    FF mShrDeltaA;
    FF mShrB;
    FF mShrDeltaB;
    FF mShrC;
    Vec mAddShrDeltaC;

    // Shares of 0
    FF mShrO1;
  };

  // Preprocessing materials for multiplication gates
  // Inherit IOBatchFIPrep
  struct MultBatchFIPrep : IOBatchFIPrep {
    // additional 0 shares
    FF mShrO2;
    FF mShrO3;
  };  

  class Correlator {
  public:
    Correlator() {};

    Correlator(std::size_t n_ind_shares, std::size_t n_mult_batches, std::size_t n_inout_batches, std::size_t batch_size) :
      mNIndShrs(n_ind_shares), mNMultBatches(n_mult_batches), mNInOutBatches(n_inout_batches), \
      mCTRIndShrs(0), mCTRMultBatches(0), mCTRInOutBatches(0), mBatchSize(batch_size) {

      mIndShrs.reserve(mNIndShrs);
      mAuthIndAddShrs.reserve(mNIndShrs);
      mMultBatchFIPrep.reserve(mNMultBatches);
      mIOBatchFIPrep.reserve(mNInOutBatches);
    }

    // Generates dummy F.I. preprocessing
    void GenIndShrsDummy(FF lambda, FF delta) {
      mIndShrs = std::vector<FF>(mNIndShrs, lambda);
      mAuthIndAddShrs = std::vector<FF>(mNIndShrs, lambda * delta); // TODO might not be naively multiplying delta
    }
    void GenMultBatchDummy(MultBatchFIPrep &mult_batch_fi_prep) {
      mMultBatchFIPrep = std::vector<MultBatchFIPrep>(mNMultBatches, mult_batch_fi_prep);      
    }
    void GenIOBatchDummy(IOBatchFIPrep &io_batch_fi_prep) {
      mIOBatchFIPrep = std::vector<IOBatchFIPrep>(mNInOutBatches, io_batch_fi_prep);
    }

    void GenPrepDummy(FF lambda, FF delta, scl::PRG prg) {
      GenIndShrsDummy(lambda, delta);

      // Fix the value of Beaver triples
      Vec triple_a(mBatchSize);
      Vec triple_a_delta(mBatchSize);
      Vec triple_b(mBatchSize);
      Vec triple_b_delta(mBatchSize);
      Vec triple_c(mBatchSize);
      Vec triple_c_delta(mBatchSize);
      Vec dummy_vec(mBatchSize);
      for(std::size_t i = 0; i < mBatchSize; ++i) {
        triple_a[i] = FF(i+1);
        triple_a_delta[i] = delta * FF(i+1);
        triple_b[i] = FF(i+1);
        triple_b_delta[i] = delta * FF(i+1);
        triple_c[i] = FF(i+1) * FF(i+1);
        triple_c_delta[i] = delta * FF(i+1) * FF(i+1);
        dummy_vec[i] = FF(0);
      }
      auto poly_triple_a = scl::details::EvPolyFromSecretsAndDegree(
          triple_a, mParties-mBatchSize, prg);
      Vec triple_a_shares = scl::details::SharesFromEvPoly(
          poly_triple_a, mParties);
      auto poly_triple_a_delta = scl::details::EvPolyFromSecretsAndDegree(
          triple_a_delta, mParties-mBatchSize, prg);
      Vec triple_a_delta_shares = scl::details::SharesFromEvPoly(
          poly_triple_a_delta, mParties);
      auto poly_triple_b = scl::details::EvPolyFromSecretsAndDegree(
          triple_b, mParties-mBatchSize, prg);
      Vec triple_b_shares = scl::details::SharesFromEvPoly(
          poly_triple_b, mParties);
      auto poly_triple_b_delta = scl::details::EvPolyFromSecretsAndDegree(
          triple_b_delta, mParties-mBatchSize, prg);
      Vec triple_b_delta_shares = scl::details::SharesFromEvPoly(
          poly_triple_b_delta, mParties);
      auto poly_triple_c = scl::details::EvPolyFromSecretsAndDegree(
          triple_c, mParties-1, prg);
      Vec triple_c_shares = scl::details::SharesFromEvPoly(
          poly_triple_c, mParties);
      auto poly_zero = scl::details::EvPolyFromSecretsAndDegree(
          dummy_vec, mParties-1, prg);
      Vec zero_shares = scl::details::SharesFromEvPoly(
          poly_zero, mParties);
      
      Vec add_shr_triple_c_delta_choice = dummy_vec;
      if(mID == 0) {
        add_shr_triple_c_delta_choice = triple_c_delta;
      }
      IOBatchFIPrep io_batch_fi_prep = {
        triple_a_shares[mID], triple_a_delta_shares[mID],
        triple_b_shares[mID], triple_b_delta_shares[mID],
        triple_c_shares[mID], add_shr_triple_c_delta_choice,
        zero_shares[mID]};
      MultBatchFIPrep mult_batch_fi_prep = {
        triple_a_shares[mID], triple_a_delta_shares[mID],
        triple_b_shares[mID], triple_b_delta_shares[mID],
        triple_c_shares[mID], add_shr_triple_c_delta_choice,
        zero_shares[mID], zero_shares[mID], zero_shares[mID]};

      GenMultBatchDummy(mult_batch_fi_prep);
      GenIOBatchDummy(io_batch_fi_prep);
    }

    /*void GenPrepDummy() {
      GenPrepDummy(FF(0));
    }*/

    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    void SetThreshold(std::size_t threshold) {
      if ( mParties != threshold + 2*(mBatchSize - 1) + 1 )
	throw std::invalid_argument("It must hold that n = t + 2(k-1) + 1");
      if ( mParties >= 2*threshold )
	throw std::invalid_argument("It must hold that t > n/2");
      mThreshold = threshold;
    }

    // GENERATE F.I. PREPROCESSING

    void OLEDummyPrep();

    void GenIndRandShr(std::size_t owner_id);
    Vec GenIndRandShr();
    void MacKeyGen();

    void UnAuthTriplesHighDeg(std::vector<Vec> &add_shr_delta_a,
      std::vector<Vec> &add_shr_delta_b,
      std::vector<Vec> &add_shr_c, std::size_t num);

    void GenPairRandShr(std::size_t owner_id,
        std::size_t n_rand_pairs);
    void GenPairRandShr(
      std::vector<std::pair<FF, FF>> &rand_pair,
      std::size_t n_rand_pairs);
    void DegreeReduce(std::vector<FF> &packed_shr_a,
      std::vector<FF> &packed_shr_b, std::size_t num_triple);

    void GenRandPackedShr(std::size_t owner_id,
      std::size_t degree, std::size_t num_shares);
    void GenRandPackedShr(Vec &shr_r,
        std::size_t degree, std::size_t num_shares);
    void AdditiveShareTransform(
      std::vector<FF> &packed_shr_delta_a,
      std::vector<FF> &packed_shr_delta_b,
      std::vector<Vec> &add_shr_delta_a,
      std::vector<Vec> &add_shr_delta_b,
      std::size_t num_triple);

    void AuthTriplesC(std::vector<FF> &packed_shr_c,
      std::vector<Vec> &add_shr_delta_c,
      std::vector<Vec> &add_shr_c,
      std::size_t num_triples);
    
    void Sacrifice(
      FF mDummyPackedShrA, FF mDummyPackedShrB,
      std::vector<FF> packed_shr_a,
      std::vector<FF> packed_shr_b,
      std::vector<Vec> packed_shr_delta_a,
      std::vector<Vec> packed_shr_delta_b,
      std::vector<FF> packed_shr_c,
      std::vector<Vec> add_shr_delta_c,
      std::size_t num_triples);

    void Triple();

    void GenZeroShr(std::size_t owner_id,
        std::size_t num_shares);
    void GenZeroShr(std::vector<FF> &zero_shares,
        std::size_t num_shares);
    void ZeroShr();

    void GenUniShr(std::size_t owner_id,
        std::size_t num_shares);
    void GenUniShr(std::vector<FF> &packed_shr_uni,
        std::size_t num_shares);
    void RandomAuthenticatedShare();

    void PrepFromDummyOle() {
      OLEDummyPrep();
      MacKeyGen();
      Triple();
      ZeroShr();
      RandomAuthenticatedShare();
    }

    // Mapping gates to preprocessed data
    void PopulateIndvShrs(std::shared_ptr<MultGate> gate) {
      if ( !gate->IsPadding() ) {
        mMapIndShrs[gate] = std::make_pair(mIndShrs[mCTRIndShrs], mAuthIndAddShrs[mCTRIndShrs]);
        mCTRIndShrs++;
      }
    }
    void PopulateIndvShrs(std::shared_ptr<InputGate> gate) {
      if ( !gate->IsPadding() ) {
        mMapIndShrs[gate] = std::make_pair(mIndShrs[mCTRIndShrs], mAuthIndAddShrs[mCTRIndShrs]);
        mCTRIndShrs++;
      }
    }
    void PopulateIndvShrs(std::shared_ptr<AddGate> gate) {
      // This is done in topological order, so previous (left&right) gates are already handled
      auto left = mMapIndShrs[gate->GetLeft()];
      auto right = mMapIndShrs[gate->GetRight()];
      mMapIndShrs[gate] = std::make_pair(left.first+right.first, left.second+right.second);
    }
    void PopulateIndvShrs(std::shared_ptr<OutputGate> gate) {
      mMapIndShrs[gate] = mMapIndShrs[gate->GetLeft()];
    }

    void PopulateMultBatches(std::shared_ptr<MultBatch> mult_batch) {
      mMapMultBatch[mult_batch] = mMultBatchFIPrep[mCTRMultBatches++];
    }

    void PopulateInputBatches(std::shared_ptr<InputBatch> input_batch) {
      mMapInputBatch[input_batch] = mIOBatchFIPrep[mCTRInOutBatches++];
    }

    void PopulateOutputBatches(std::shared_ptr<OutputBatch> output_batch) {
      mMapOutputBatch[output_batch] = mIOBatchFIPrep[mCTRInOutBatches++];
    }

    // Generate FD Prep from FI Prep
    
    // PREP INPUT BATCH
    void PrepInput(std::shared_ptr<InputBatch> input_batch);

    // PREP OUTPUT BATCH
    void PrepOutput(std::shared_ptr<OutputBatch> output_batch);

    // PREP MULT BATCH
    void PrepMultPartiesSendP1(std::shared_ptr<MultBatch> mult_batch);
    
    void PrepMultP1Receives(std::shared_ptr<MultBatch> mult_batch);

    // Populate shares of e_i
    void PrecomputeEi() {
      for (std::size_t i = 0; i < mBatchSize; i++) {
	FF shr(1);
	for (std::size_t j = 0; j < mBatchSize; ++j) {
	  if (j == i) continue;
	  shr *= (FF(mID+1) + FF(j)) / (FF(j) - FF(i));
	}
	mSharesOfEi.emplace_back(shr);
      }
    }

    // Populate vandermonde matrix
    void PrecomputeVandermonde() {
      mVandermonde.reserve(mParties);
      for (std::size_t i = 0; i < mParties; i++) {
	mVandermonde.emplace_back(std::vector<FF>());
	mVandermonde[i].reserve(mThreshold + 1);
	FF entry(1);
	for (std::size_t j = 0; j < mThreshold + 1; ++j) {
	  mVandermonde[i].emplace_back(entry);
	  entry *= FF(i);
	}
      }
    }

    // Maps
    std::map<std::shared_ptr<Gate>, std::pair<FF, FF>> mMapIndShrs;
    std::map<std::shared_ptr<MultBatch>, MultBatchFIPrep> mMapMultBatch;
    std::map<std::shared_ptr<InputBatch>, IOBatchFIPrep> mMapInputBatch;
    std::map<std::shared_ptr<OutputBatch>, IOBatchFIPrep> mMapOutputBatch;

  private:
    // Sizes
    std::size_t mNIndShrs;
    std::size_t mNMultBatches;
    std::size_t mNInOutBatches;

    // Counters
    std::size_t mCTRIndShrs;
    std::size_t mCTRMultBatches;
    std::size_t mCTRInOutBatches;

    std::size_t mBatchSize;


    // Shares of e_i for current party
    std::vector<FF> mSharesOfEi; // len = batchsize
    std::vector<std::vector<FF>> mVandermonde; // mParties x (mThreshold + 1)

    // dummy ole
    FF mAddShrDelta;
    Vec mShamirShrDelta;
    FF mDummyPackedShrA;
    FF mDummyPackedShrB;
    FF mDummyPackedShrC;
    Vec mDummyAddVoleShrA;
    Vec mDummyAddVoleShrB;
    Vec mDummyAddOleShrAiBj;
    Vec mDummyAddOleShrAjBi;

    FF mDummyAddShrDeltaA;
    FF mDummyAddShrDeltaB;

    // SHARINGS

    // Individual sharings
    std::vector<FF> mIndShrs;
    std::vector<FF> mAuthIndAddShrs;

    // Multiplication batches
    std::vector<MultBatchFIPrep> mMultBatchFIPrep;

    // Input & Output batches
    std::vector<IOBatchFIPrep> mIOBatchFIPrep;

    // HELPERS TO GET THE F.I. PREP
    std::vector<std::vector<FF>> mUnpackedShrsA; // idx = packed index
    std::vector<std::vector<FF>> mUnpackedShrsB; // idx = packed index
    std::vector<std::vector<FF>> mUnpackedShrsMask; // idx = packed index
    std::vector<std::vector<FF>> mZeroProdShrs; // idx = packed index

    

    // NETWORK-RELATED
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    std::size_t mThreshold;

    scl::PRG mPRG;
  };

} // namespace tp

#endif  // CORRELATOR_H
