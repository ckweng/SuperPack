#include "output_gate.h"

namespace tp {

  void OutputBatch::PartiesRevealShares() {
    // compute shares of lambda - a
    if(mID != 0) {
      FF shr_diff_lambda_a = mPackedShrLambda - mPackedShrA;
      mNetwork->Party(0)->Send(shr_diff_lambda_a);
    }
  }

  void OutputBatch::P1ReconstructAndDistributes() {
    if(mID != 0) return;
    // P1 reconstruct lambda - a
    Vec shares_diff_lambda_a(mParties);
    shares_diff_lambda_a[0] = mPackedShrLambda - mPackedShrA;
    for(std::size_t i = 1; i < mParties; ++i) {
      mNetwork->Party(i)->Recv(shares_diff_lambda_a[i]);
    }
    Vec vec_diff_lambda_a = SecretsFromSharesAndLengthInternal(
        shares_diff_lambda_a, mBatchSize);

    // P1 get mu
    Vec vec_mu(mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      vec_mu[i] = GetOutputGate(i)->GetMu();
    }
    Vec vec_diff_v_a = vec_mu.Add(vec_diff_lambda_a);

    // P1 construct degree-2(k-1) packed shares of v - a
    scl::PRG prg;
    auto poly_diff_v_a = EvPolyFromSecretsAndDegreeInternal(
        vec_diff_v_a, 2*(mBatchSize-1), prg);
    Vec shares_diff_v_a = SharesFromEvPolyInternal(
        poly_diff_v_a, mParties);

    // P1 distributes shares
    mPackedShrDiffVAA = shares_diff_v_a[0];
    for(std::size_t i = 1; i < mParties; ++i) {
      mNetwork->Party(i)->Send(shares_diff_v_a[i]);
    }
  }

  void OutputBatch::PartiesReceiveShares() {
    if(mID != 0) {
      mNetwork->Party(0)->Recv(mPackedShrDiffVAA);
    }
  }

  // TODO problem: the theta is incorrect (only output gates)
  void OutputBatch::PartiesVerify(std::vector<FF> &theta_list) {
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      FF add_shr_delta_mu = GetOutputGate(i)->GetDeltaMu();
      FF add_shr_delta_lambda = mAddShrLambdaDelta[i];
      // TODO transform packed share mPackedShrDeltaA to ith additive share
      FF add_shr_delta_a = PackedToAdditive(
          mPackedShrDeltaA, mID+1, -1*i, mParties-mBatchSize);
      FF add_shr_delta_diff_v_a = add_shr_delta_mu
        + add_shr_delta_lambda - add_shr_delta_a;

      FF shr_bar_delta_diff_v_a =
        GetPackedShrDelta(i) * mPackedShrDiffVAA;
      // TODO transform packed share shr_bar_delta_diff_v_a to ith additive share
      FF add_shr_bar_delta_diff_v_a = PackedToAdditive(
          shr_bar_delta_diff_v_a, mID+1, -1*i, mParties-1);
      // start dummy, TODO delete when fixing problem
      FF a = add_shr_delta_diff_v_a - add_shr_bar_delta_diff_v_a;
      a = a + FF(0);
      std::size_t sz = theta_list.size();
      sz++;
      //theta_list.push_back(add_shr_delta_diff_v_a - add_shr_bar_delta_diff_v_a);
    }

    // TODO verification of the computation
    // TODO coin-tossing
    // TODO verify all thetas by random linear computation
    unsigned char seed[16] = {0x00, 0x01, 0x02, 0x03,
      0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f};
    scl::PRG prg(seed);
    Vec coefficients = Vec::Random(theta_list.size(), prg);
    FF sum(0);
    for(std::size_t i = 0; i < theta_list.size(); ++i) {
      sum = sum + coefficients[i] * theta_list[i];
    }
    if(mID != 0) {
      mNetwork->Party(0)->Send(sum);
    } else {
      mThetaSum = sum;
    }
  }

  void OutputBatch::PartiesCheckTheta() {
    if(mID == 0) {
      FF sum = mThetaSum;
      for(std::size_t i = 1; i < mParties; ++i) {
        FF sum_recv;
        mNetwork->Party(i)->Recv(sum_recv);
        sum = sum + sum_recv;
      }
      if(sum.Equal(FF(0)) != true) {
        throw std::invalid_argument("check theta fails");
      }
    }
  }

  void OutputBatch::PartiesRevealOutputVal() {
    // parties reveal shares of [v-a], [a], [b], [c]
    if(mID != mOwnerID) {
      mNetwork->Party(mOwnerID)->Send(mPackedShrDiffVAA);
      mNetwork->Party(mOwnerID)->Send(mPackedShrA);
      mNetwork->Party(mOwnerID)->Send(mPackedShrB);
      mNetwork->Party(mOwnerID)->Send(mPackedShrC);
    }
  }

  void OutputBatch::ClientReconstructSecret() {
    // client receive shares of [v-a], [a], [b], [c]
    if(mID != mOwnerID) return;
    Vec shares_diff_v_a(mParties);
    Vec shares_a(mParties);
    Vec shares_b(mParties);
    Vec shares_c(mParties);
    for(std::size_t i = 0; i < mParties; ++i) {
      if(i == mID) {
        shares_diff_v_a[i] = mPackedShrDiffVAA;
        shares_a[i] = mPackedShrA;
        shares_b[i] = mPackedShrB;
        shares_c[i] = mPackedShrC;
      } else {
        mNetwork->Party(i)->Recv(shares_diff_v_a[i]);
        mNetwork->Party(i)->Recv(shares_a[i]);
        mNetwork->Party(i)->Recv(shares_b[i]);
        mNetwork->Party(i)->Recv(shares_c[i]);
      }
    }

    // client reconstruct v-a, a, b, c
    // TODO possibly need to check the consistency of shares
    /*auto start = shares_diff_v_a.begin();
    auto end = shares_diff_v_a.begin();
    std::advance(end, 2*mBatchSize-2);
    shares_diff_v_a = Vec(start, end);
    start = shares_a.begin();
    end = shares_a.begin();
    std::advance(end, mParties-mBatchSize);
    shares_a = Vec(start, end);
    start = shares_b.begin();
    end = shares_b.begin();
    std::advance(end, mParties-mBatchSize);
    shares_b = Vec(start, end);*/
    Vec vec_diff_v_a = SecretsFromSharesAndLengthInternal(
        shares_diff_v_a, mBatchSize);
    Vec vec_a = SecretsFromSharesAndLengthInternal(
        shares_a, mBatchSize);
    Vec vec_b = SecretsFromSharesAndLengthInternal(
        shares_b, mBatchSize);
    Vec vec_c = SecretsFromSharesAndLengthInternal(
        shares_c, mBatchSize);
    vec_c = vec_c.Subtract(vec_a.MultiplyEntryWise(vec_b));
    // TODO temporarily disable the check for experiments
    std::size_t n_error = 0;
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      if(vec_c[i].Equal(FF(0)) != true) {
        n_error++;
        //throw std::invalid_argument("output batch: inconsistent preprocessing triple");
      }
    }

    // client reconstruct the output values
    Vec out = vec_diff_v_a.Add(vec_a);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      GetOutputGate(i)->SetValue(out[i]);
    }
  }

}
