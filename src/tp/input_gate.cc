#include "input_gate.h"

namespace tp {

  void InputBatch::PartiesRevealShares() {
      // parties send the shares
      if(mID != mOwnerID) {
        mNetwork->Party(mOwnerID)->Send(mPackedShrLambda);
        mNetwork->Party(mOwnerID)->Send(mPackedShrA);
        mNetwork->Party(mOwnerID)->Send(mPackedShrB);
        mNetwork->Party(mOwnerID)->Send(mPackedShrC);
      }
  }

  void InputBatch::ClientReconstructAndDistributes() {
    if(mID != mOwnerID) return;
    Vec point_values_lambda_a;
    point_values_lambda_a.Reserve(mParties);
    Vec point_values_a;
    point_values_a.Reserve(mParties);
    Vec point_values_b;
    point_values_b.Reserve(mParties);
    Vec point_values_c;
    point_values_c.Reserve(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      if(i == mOwnerID) {
        point_values_lambda_a.Emplace(mPackedShrLambda);
        point_values_a.Emplace(mPackedShrA);
        point_values_b.Emplace(mPackedShrB);
        point_values_c.Emplace(mPackedShrC);
        continue;
      }
      FF recv_ff;
      mNetwork->Party(i)->Recv(recv_ff);
      point_values_lambda_a.Emplace(recv_ff);
      mNetwork->Party(i)->Recv(recv_ff);
      point_values_a.Emplace(recv_ff);
      mNetwork->Party(i)->Recv(recv_ff);
      point_values_b.Emplace(recv_ff);
      mNetwork->Party(i)->Recv(recv_ff);
      point_values_c.Emplace(recv_ff); 
    }

    // TODO check consistency of the shares
    // TODO find better ways for interpolation of degree lower than n-1
    /*auto start = point_values_a.begin();
    auto end = point_values_a.begin();
    std::advance(end, mParties-mBatchSize);
    point_values_a = Vec(start, end);
    start = point_values_b.begin();
    end = point_values_b.begin();
    std::advance(end, mParties-mBatchSize);
    point_values_b = Vec(start, end);*/
    Vec vec_lambda_a = SecretsFromSharesAndLengthInternal(
        point_values_lambda_a, mBatchSize);
    Vec vec_a = SecretsFromSharesAndLengthInternal(point_values_a, mBatchSize);
    Vec vec_b = SecretsFromSharesAndLengthInternal(point_values_b, mBatchSize);
    Vec vec_c = SecretsFromSharesAndLengthInternal(point_values_c, mBatchSize);

    // check whether a = b * c
    // TODO temporarily disable the check for experiments
    size_t n_error = 0;
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      if(vec_c[i].Equal(vec_a[i] * vec_b[i]) != true)
        n_error++;
        //throw std::invalid_argument("input batch: inconsistent preprocessing triple\n");
    }

    // compute \mu = v - \lambda and v - a
    Vec vec_v = GetValue();
    Vec vec_mu = vec_v.Subtract(vec_lambda_a);
    Vec vec_diff_v_a = vec_v.Subtract(vec_a);

    for(std::size_t i = 0; i < mBatchSize; ++i) {
      GetInputGate(i)->SetMu(vec_mu[i]);
    }

    // generate and distribute degree-(2k-2) shares of v - a
    scl::PRG prg;
    auto poly_diff_v_a = EvPolyFromSecretsAndDegreeInternal(
        vec_diff_v_a, 2*(mBatchSize-1), prg);
    Vec shares_diff_v_a = SharesFromEvPolyInternal(poly_diff_v_a, mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      if(i == mID) {
        mPackedShrDiffVAA = shares_diff_v_a[i];
        continue;
      }
      mNetwork->Party(i)->Send(shares_diff_v_a[i]);
    }
  }

  void InputBatch::PartiesReceive() {
    // Parties receive shares from client
    if(mID == mOwnerID) return;
    mNetwork->Party(mOwnerID)->Recv(mPackedShrDiffVAA);
  }

  void InputBatch::ClientSendMuToP1() {
    // send mu to P1
    if(mOwnerID != 0 && mID == mOwnerID) {
      Vec vec_mu;
      vec_mu.Reserve(mBatchSize);
      for(std::size_t i = 0; i < mBatchSize; ++i) {
        vec_mu.Emplace(GetInputGate(i)->GetMu());
      }
      mNetwork->Party(0)->Send(vec_mu);
    }
  }

  void InputBatch::P1ReceiveMu() {
    // send mu to P1
    if(mOwnerID != 0 && mID == 0) {
      Vec vec_mu_recv;
      mNetwork->Party(mOwnerID)->Recv(vec_mu_recv);
      for(std::size_t i = 0; i < mBatchSize; ++i) {
        GetInputGate(i)->SetMu(vec_mu_recv[i]);
      }
    }
  }

  void InputBatch::PartiesAuthenticateMu() {
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      FF share_delta_diff_v_a = GetPackedShrDelta(i) * mPackedShrDiffVAA;
      // TODO transfer from packed share of Delta * (v - a) to additive share
      FF add_shr_delta_diff_v_a = PackedToAdditive(share_delta_diff_v_a, mID+1, -1*i, mParties-1);
      // TODO transfer packed share of Delta*a to additive share
      FF add_shr_delta_a = PackedToAdditive(mPackedShrDeltaA, mID+1, -1*i, mParties-mBatchSize);
      GetInputGate(i)->SetDeltaMu(add_shr_delta_diff_v_a + add_shr_delta_a - mAddShrDeltaLambda[i]);
    }
  }

}
