#include "tp/correlator.h"

namespace tp {
  


  // PREP INPUT BATCH
  void Correlator::PrepInput(std::shared_ptr<InputBatch> input_batch) {
    // collect share of [lambda_gamma]_n-1
    FF shr_lambda_gamma(0);
    Vec add_shr_delta_lambda_gamma(mBatchSize);
    for (std::size_t i = 0; i < mBatchSize; i++) {
      auto input_gate = input_batch->GetInputGate(i);
      shr_lambda_gamma += mSharesOfEi[i] * mMapIndShrs[input_gate].first;
      add_shr_delta_lambda_gamma[i] = mMapIndShrs[input_gate].second;
    }

    // compute shares of [lambda_gamma + 0]_n-1
    auto input_prep_info = mMapInputBatch[input_batch];
    shr_lambda_gamma += input_prep_info.mShrO1;

    // set preprocessing
    input_batch->SetPreprocessing(
        shr_lambda_gamma, add_shr_delta_lambda_gamma,
        input_prep_info.mShrA, input_prep_info.mShrDeltaA,
        input_prep_info.mShrB, input_prep_info.mShrDeltaB,
        input_prep_info.mShrC, input_prep_info.mAddShrDeltaC);
  }

  // PREP OUTPUT BATCH
  void Correlator::PrepOutput(std::shared_ptr<OutputBatch> output_batch) {
    // collect share of [lambda_gamma]_n-1
    FF shr_lambda_gamma(0);
    Vec add_shr_delta_lambda_gamma(mBatchSize);
    for (std::size_t i = 0; i < mBatchSize; i++) {
      auto output_gate = output_batch->GetOutputGate(i);
      shr_lambda_gamma += mSharesOfEi[i] * mMapIndShrs[output_gate].first;
      add_shr_delta_lambda_gamma[i] = mMapIndShrs[output_gate].second;
    }

    // compute shares of [lambda_gamma + 0]_n-1
    auto output_prep_info = mMapOutputBatch[output_batch];
    shr_lambda_gamma += output_prep_info.mShrO1;

    // set preprocessing
    output_batch->SetPreprocessing(
        shr_lambda_gamma, add_shr_delta_lambda_gamma,
        output_prep_info.mShrA, output_prep_info.mShrDeltaA,
        output_prep_info.mShrB, output_prep_info.mShrDeltaB,
        output_prep_info.mShrC, output_prep_info.mAddShrDeltaC);
  }

  // PREP MULT BATCH
  // Reveal (lambda_alpha - a), (lambda_beta - b) to P1
  // Generate their authenticated additive shares
  void Correlator::PrepMultPartiesSendP1(std::shared_ptr<MultBatch> mult_batch) {
    // collect share of [lambda_alpha]_n-1, [lambda_beta]_n-1
    // and [lambda_gamma]_n-1
    FF shr_diff_lambda_a_a(0);
    FF shr_diff_lambda_b_b(0);
    FF shr_lambda_gamma(0);
    Vec add_shr_delta_lambda_alpha(mBatchSize);
    Vec add_shr_delta_lambda_beta(mBatchSize);
    Vec add_shr_delta_lambda_gamma(mBatchSize);
    for (std::size_t i = 0; i < mBatchSize; i++) {
      auto mult_gate = mult_batch->GetMultGate(i);
      shr_lambda_gamma += mSharesOfEi[i] * mMapIndShrs[mult_gate].first;
      add_shr_delta_lambda_gamma[i] = mMapIndShrs[mult_gate].second;
      shr_diff_lambda_a_a += mSharesOfEi[i] * mMapIndShrs[mult_gate->GetLeft()].first;
      add_shr_delta_lambda_alpha[i] = mMapIndShrs[mult_gate->GetLeft()].second;
      shr_diff_lambda_b_b += mSharesOfEi[i] * mMapIndShrs[mult_gate->GetRight()].first;
      add_shr_delta_lambda_beta[i] = mMapIndShrs[mult_gate->GetRight()].second;
    }

    // compute shares of [lambda_alpha - a + 0]_n-1, [lambda_beta - b + 0]_n-1
    // and [lambda_gamma + 0]_n-1
    auto mult_prep_info = mMapMultBatch[mult_batch];
    shr_diff_lambda_a_a += (mult_prep_info.mShrO1 - mult_prep_info.mShrA);
    shr_diff_lambda_b_b += (mult_prep_info.mShrO2 - mult_prep_info.mShrB);
    shr_lambda_gamma += mult_prep_info.mShrO3;

    Vec add_shr_diff_lambda_alpha_a(mBatchSize);
    Vec add_shr_diff_lambda_beta_b(mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      // convert [Delta * a]_n-k to additive shares
      FF add_shr_delta_a = PackedToAdditive(mult_prep_info.mShrDeltaA, mID+1, -1*i, mParties-mBatchSize);
      // output <Delta * (lambda_alpha - a)>
      add_shr_diff_lambda_alpha_a[i] = add_shr_delta_lambda_alpha[i] - add_shr_delta_a;
      // convert [Delta * b]_n-k to additive shares
      FF add_shr_delta_b = PackedToAdditive(mult_prep_info.mShrDeltaB, mID+1, -1*i, mParties-mBatchSize);
      // output <Delta * (lambda_beta - b)>
      add_shr_diff_lambda_beta_b[i] = add_shr_delta_lambda_beta[i] - add_shr_delta_b;
    }

    // set preprocessing
    mult_batch->SetPreprocessing(add_shr_diff_lambda_alpha_a,
        add_shr_diff_lambda_beta_b, shr_lambda_gamma, add_shr_delta_lambda_gamma);
    mult_batch->SetPreprocessingTriple(
        mult_prep_info.mShrA, mult_prep_info.mShrDeltaA,
        mult_prep_info.mShrB, mult_prep_info.mShrDeltaB,
        mult_prep_info.mShrC, mult_prep_info.mAddShrDeltaC);

    // send shares to P1
    mNetwork->Party(0)->Send(shr_diff_lambda_a_a);
    mNetwork->Party(0)->Send(shr_diff_lambda_b_b);    
  }
    
  void Correlator::PrepMultP1Receives(std::shared_ptr<MultBatch> mult_batch) {
    if(mID != 0) return;

    Vec shares_diff_lambda_a_a;
    Vec shares_diff_lambda_b_b;
    shares_diff_lambda_a_a.Reserve(mParties);
    shares_diff_lambda_b_b.Reserve(mParties);

    // P1 receives
    for (std::size_t i = 0; i < mParties; i++) {
      FF buffer;
      mNetwork->Party(i)->Recv(buffer);
      shares_diff_lambda_a_a.Emplace(buffer);
      mNetwork->Party(i)->Recv(buffer);
      shares_diff_lambda_b_b.Emplace(buffer);
    }
    Vec vec_diff_lambda_a_a = SecretsFromSharesAndLengthInternal(
        shares_diff_lambda_a_a, mBatchSize);
    Vec vec_diff_lambda_b_b = SecretsFromSharesAndLengthInternal(
        shares_diff_lambda_b_b, mBatchSize);

    // set values to mult_batch
    mult_batch->SetP1Preprocessing(vec_diff_lambda_a_a,
        vec_diff_lambda_b_b);
  }
}

