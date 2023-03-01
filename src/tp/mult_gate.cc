#include "tp/mult_gate.h"

namespace tp {

  void MultBatch::P1Distributes() {
    if (mID != 0) return;

    // P1 computes v_\alpha - a, v_\beta - b
    Vec vec_mu_a(mBatchSize);    
    Vec vec_mu_b(mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      vec_mu_a[i] = mMultGatesPtrs[i]->GetLeft()->GetMu();
      vec_mu_b[i] = mMultGatesPtrs[i]->GetRight()->GetMu();
    }
    Vec vec_diff_va_a = vec_mu_a.Add(mVecDiffLambdaAA);
    Vec vec_diff_vb_b = vec_mu_b.Add(mVecDiffLambdaBB);

    // generate packed shares of v_\alpha - a, v_\beta - b
    scl::PRG prg;
    auto poly_diff_va_a = EvPolyFromSecretsAndDegreeInternal(
        vec_diff_va_a, mBatchSize-1, prg);
    auto poly_diff_vb_b = EvPolyFromSecretsAndDegreeInternal(
        vec_diff_vb_b, mBatchSize-1, prg);
    Vec shares_diff_va_a = SharesFromEvPolyInternal(
        poly_diff_va_a, mParties);
    Vec shares_diff_vb_b = SharesFromEvPolyInternal(
        poly_diff_vb_b, mParties); 
    
    // distribute shares to parties
    mPackedShrDiffVAA = shares_diff_va_a[0];
    mPackedShrDiffVBB = shares_diff_vb_b[0];
    for(std::size_t i = 1; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shares_diff_va_a[i]);
      mNetwork->Party(i)->Send(shares_diff_vb_b[i]);
    }
  }

  void MultBatch::PartiesReceive() {
    if(mID == 0) return;
    mNetwork->Party(0)->Recv(mPackedShrDiffVAA);
    mNetwork->Party(0)->Recv(mPackedShrDiffVBB);
  }

  void MultBatch::PartiesMultiplyAndSend() {
    if(mID == 0) return;
    // Compute share of mu_\gamma
    FF shr_mu_C = mPackedShrDiffVAA * mPackedShrDiffVBB
      + mPackedShrDiffVAA * mPackedShrB
      + mPackedShrDiffVBB * mPackedShrA
      + mPackedShrC - mPackedShrDeltaC;

    // Send share to P1
    mNetwork->Party(0)->Send(shr_mu_C);
  }

  void MultBatch::P1ReceivesAndReconstruct() {
    if(mID != 0) return;
    // P1 computes share of \mu_\gamma
    Vec shares_mu_c(mParties);
    shares_mu_c[0] = mPackedShrDiffVAA * mPackedShrDiffVBB
      + mPackedShrDiffVAA * mPackedShrB
      + mPackedShrDiffVBB * mPackedShrA
      + mPackedShrC - mPackedShrDeltaC;

    // P1 receives shares of \mu_\gamma
    for(std::size_t i = 1; i < mParties; ++i) {
      mNetwork->Party(i)->Recv(shares_mu_c[i]);
    }

    // P1 reconstruct secrets of \mu_\gamma
    Vec vec_mu_c = SecretsFromSharesAndLengthInternal(
        shares_mu_c, mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      GetMultGate(i)->SetMu(vec_mu_c[i]);
    }
  }

  void MultBatch::PartiesCheckAndComputeAuthMu(std::vector<FF> &theta_list) {
    // compute theta
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      // Get \mu_\alpha and \mu_\beta
      FF add_shr_delta_mu_a = mMultGatesPtrs[i]->GetLeft()->GetDeltaMu();
      FF add_shr_delta_mu_b = mMultGatesPtrs[i]->GetRight()->GetDeltaMu();
      // Compute additive share of \Delta*(v_\alpha-a), \Delta*(v_\beta-b)
      FF add_shr_delta_diff_va_a = add_shr_delta_mu_a + mAddShrDeltaDiffLambdaAA[i];
      FF add_shr_delta_diff_vb_b = add_shr_delta_mu_b + mAddShrDeltaDiffLambdaBB[i];

      FF shr_bar_delta_diff_va_a =
        mMultGatesPtrs[0]->mPackedShrDelta[i] * mPackedShrDiffVAA;
      FF shr_bar_delta_diff_vb_b =
        mMultGatesPtrs[0]->mPackedShrDelta[i] * mPackedShrDiffVBB;
      // TODO transform packed share shr_bar_delta_diff_va_a to ith additive share
      FF add_shr_bar_delta_diff_va_a = PackedToAdditive(
          shr_bar_delta_diff_va_a, mID+1, -1*i, mParties-mBatchSize);
      // TODO transform packed share shr_bar_delta_diff_vb_b to ith additive share
      FF add_shr_bar_delta_diff_vb_b = PackedToAdditive(
          shr_bar_delta_diff_vb_b, mID+1, -1*i, mParties-mBatchSize);
      theta_list.push_back(add_shr_delta_diff_va_a - add_shr_bar_delta_diff_va_a);
      theta_list.push_back(add_shr_delta_diff_vb_b - add_shr_bar_delta_diff_vb_b);
    }


    // compute authenticated mu
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      // compute additive share of \Delta*(v_\alpha-a)*(v_\beta-b)
      FF shr_delta_diff_va_a_diff_vb_b =
        mMultGatesPtrs[0]->mPackedShrDelta[i]
        * mPackedShrDiffVAA * mPackedShrDiffVBB;
      // TODO transform from packed sharing to additive shares
      FF add_shr_delta_diff_va_a_diff_vb_b = PackedToAdditive(
          shr_delta_diff_va_a_diff_vb_b, mID+1, -1*i, mParties-1);

      FF shr_cross_term = mPackedShrDiffVAA * mPackedShrDeltaB
        + mPackedShrDiffVBB * mPackedShrDeltaA;
      // TODO transform from packed sharing to additive shares
      FF add_shr_cross_term = PackedToAdditive(
          shr_cross_term, mID+1, -1*i, mParties-1);

      GetMultGate(i)->SetDeltaMu(add_shr_delta_diff_va_a_diff_vb_b
        + add_shr_cross_term + mAddShrDeltaC[i] - mAddShrDeltaDeltaC[i]);
    }
  }
    
} // namespace tp
  
