#ifndef MULT_GATE_H
#define MULT_GATE_H

#include <vector>
#include <assert.h>

#include "gate.h"

namespace tp {
  class MultGate : public Gate {
  public:
    MultGate() {}; // for the padding gates below
    MultGate(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
      mLeft = left;
      mRight = right;
    };

    FF GetMu() {
      if ( !mLearned ) {
	throw std::invalid_argument("mu is not set");
      }
      return mMu;
    };

    void SetMu(FF mu) {
      // TODO temporarily disable the check
      /*if ( mLearned ) {
	throw std::invalid_argument("already set mu");
      }*/
      mMu = mu;
      mLearned = true;
    };

    FF GetDeltaMu() {
      /*if ( !mAuthLearned ) {
	throw std::invalid_argument("auth mu is not set");
      }*/
      return mDeltaMu;
    };

    void SetDeltaMu(FF delta_mu) {
      /*if ( mAuthLearned ) {
	throw std::invalid_argument("already set auth mu");
      }*/
      mDeltaMu = delta_mu;
      mAuthLearned = true;
    };

    void SetDummyLambda(FF lambda) {
      mLambda = lambda;
      mLambdaSet = true;
    };

    FF GetDummyLambda() {
      if ( !mLambdaSet ) {
	throw std::invalid_argument("Lambda is not set in this multiplication gate");
      }
      return mLambda;
    };

    void SetIndvShrLambda(FF indv_shr) {
      mIndvShrLambdaC = indv_shr;
      mIndvShrLambdaCSet = true;
    }

    FF GetIndvShrLambda() {
      if ( !mIndvShrLambdaCSet )
	throw std::invalid_argument("IndvShrLambda is not set in this multiplication gate");
      return mIndvShrLambdaC;
    }

    FF GetDn07Share() {
      if ( !mDn07Set )
	throw std::invalid_argument("Dn07 shares is not set in this multiplication gate");
      return mDn07Share;
    }

    FF GetClear() {
      if ( !mEvaluated ) {
	mClear = mLeft->GetClear() * mRight->GetClear();
	mEvaluated = true;
      }
      return mClear;
    }

  private:
  };

  // Used for padding batched multiplications
  class PadMultGate : public MultGate {
  public:
    PadMultGate() : MultGate() { mIsPadding = true; };

    FF GetMu() override { return FF(0); }
    FF GetDummyLambda() override { return FF(0); }

    // Just a technicality needed to make the parents of this gate be
    // itself when instantiated
    void UpdateParents(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
      mLeft = left;
      mRight = right;
    }

  private:
  };
    
  class MultBatch {
  public:
    MultBatch(std::size_t batch_size) : mBatchSize(batch_size) {
      mMultGatesPtrs.reserve(mBatchSize);
    };

    // Adds a new mult_gate to the batch. It cannot add more gates than
    // the batch_size
    void Append(std::shared_ptr<MultGate> mult_gate) {
      if ( mMultGatesPtrs.size() == mBatchSize )
	throw std::invalid_argument("Trying to batch more than batch_size gates");
      mMultGatesPtrs.emplace_back(mult_gate); };

    // For testing purposes: sets the required preprocessing for this
    // batch to be just constant shares
    /*void _DummyPrep(FF lambda_A, FF lambda_B, FF lambda_C) {
      if ( mMultGatesPtrs.size() != mBatchSize )
	throw std::invalid_argument("The number of mult gates does not match the batch size");

      mPackedShrLambdaA = lambda_A;
      mPackedShrLambdaB = lambda_B;
      mPackedShrDeltaC = lambda_A * lambda_B - lambda_C;
    };*/

    /*void _DummyPrep() {
      _DummyPrep(FF(0), FF(0), FF(0));
    };*/
    
    // Generates the preprocessing from the lambdas of the inputs
    /*void PrepFromDummyLambdas() {
      Vec lambda_A;
      Vec lambda_B;
      Vec delta_C;

      for (std::size_t i = 0; i < mBatchSize; i++) {
	auto l_A = mMultGatesPtrs[i]->GetLeft()->GetDummyLambda();
	auto l_B = mMultGatesPtrs[i]->GetRight()->GetDummyLambda();
	auto l_C = mMultGatesPtrs[i]->GetDummyLambda();
	auto d_C = l_A * l_B - l_C;

	lambda_A.Emplace(l_A);
	lambda_B.Emplace(l_B);
	delta_C.Emplace(d_C);
      }
      
      // Using deg = BatchSize-1 ensures there's no randomness involved
      auto poly_A = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mBatchSize-1, mPRG);
      mPackedShrLambdaA = poly_A.Evaluate(FF(mID));
	
      auto poly_B = scl::details::EvPolyFromSecretsAndDegree(lambda_B, mBatchSize-1, mPRG);
      mPackedShrLambdaB = poly_B.Evaluate(FF(mID));

      auto poly_C = scl::details::EvPolyFromSecretsAndDegree(delta_C, mBatchSize-1, mPRG);
      mPackedShrDeltaC = poly_C.Evaluate(FF(mID));
    }*/

    // For cleartext evaluation: calls GetClear on all its gates to
    // populate their mClear. This could return a vector with these
    // values but we're not needing them
    void GetClear() {
      for (auto gate : mMultGatesPtrs) { gate->GetClear(); }
    }

    // Determines whether the batch is full
    bool HasRoom() { return mMultGatesPtrs.size() < mBatchSize; }
    
    // For fetching mult gates
    std::shared_ptr<MultGate> GetMultGate(std::size_t idx) { return mMultGatesPtrs[idx]; }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the batch since sometimes we just want
    // to evaluate in the clear and this won't be needed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    // Online protocol
    // P1 distributes the shares of v_\alpha - a and v_\beta - b
    void P1Distributes();

    // The parties receive the packed shares
    void PartiesReceive();

    // The parties compute \mu_\gamma and send the shares back
    void PartiesMultiplyAndSend();

    // P1 receives the shares from the parties, reconstructs the \mu_\gamma
    void P1ReceivesAndReconstruct();

    // Parties check and compute \mu for previous batch and this batch
    void PartiesCheckAndComputeAuthMu(std::vector<FF> &theta_list);

    /*void RunProtocol() {
      P1Sends();
      PartiesReceive();
      PartiesSend();
      P1Receives();
    }*/

    void SetPreprocessing(Vec &shr_delta_diff_lambda_A_A,
        Vec &shr_delta_diff_lambda_B_B,
        FF shr_lambda_C,
        Vec &shr_delta_lambda_C) {
      mAddShrDeltaDiffLambdaAA = Vec(shr_delta_diff_lambda_A_A.begin(), shr_delta_diff_lambda_A_A.end());
      mAddShrDeltaDiffLambdaBB = Vec(shr_delta_diff_lambda_B_B.begin(), shr_delta_diff_lambda_B_B.end());
      mPackedShrDeltaC = shr_lambda_C;
      mAddShrDeltaDeltaC = Vec(shr_delta_lambda_C.begin(), shr_delta_lambda_C.end());
    }

    void SetPreprocessingTriple(FF shr_a, FF shr_a_delta,
        FF shr_b, FF shr_b_delta,
        FF shr_c, Vec &shr_c_delta) {
      mPackedShrA = shr_a;
      mPackedShrDeltaA = shr_a_delta;
      mPackedShrB = shr_b;
      mPackedShrDeltaB = shr_b_delta;
      mPackedShrC = shr_c;
      mAddShrDeltaC = Vec(shr_c_delta.begin(), shr_c_delta.end());
    }

    void SetP1Preprocessing(Vec diff_lambda_a_a, Vec diff_lambda_b_b) {
      mVecDiffLambdaAA = Vec(diff_lambda_a_a.begin(), diff_lambda_a_a.end());
      mVecDiffLambdaBB = Vec(diff_lambda_b_b.begin(), diff_lambda_b_b.end());
    }

  private:
    std::size_t mBatchSize;

    // The mult gates that are part of this batch
    vec<std::shared_ptr<MultGate>> mMultGatesPtrs;

    // The packed sharings associated to this batch
    FF mPackedShrDeltaC; // \lambda_\gamma
    Vec mAddShrDeltaDeltaC;
    Vec mAddShrDeltaDiffLambdaAA;
    Vec mAddShrDeltaDiffLambdaBB;
    // The semi-authenticated triple
    FF mPackedShrA;
    FF mPackedShrDeltaA;
    FF mPackedShrB;
    FF mPackedShrDeltaB;
    FF mPackedShrC;
    Vec mAddShrDeltaC;
    // known to P1
    Vec mVecDiffLambdaAA;
    Vec mVecDiffLambdaBB;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    // Intermediate-protocol
    scl::PRG mPRG;
    FF mPackedShrMuA;    // Shares of mu_alpha
    FF mPackedShrMuB;    // Shares of mu_beta
    FF mPackedShrDiffVAA;
    FF mPackedShrDiffVBB;
  };

  // Basically a collection of batches
  class MultLayer {
  public:
    MultLayer(std::size_t batch_size) : mBatchSize(batch_size) {
      auto first_batch = std::make_shared<MultBatch>(mBatchSize);
      // Append a first batch
      mBatches.emplace_back(first_batch);
    };

    // Adds a new mult_gate to the layer. It checks if the current
    // batch is full and if so creates a new one.
    void Append(std::shared_ptr<MultGate> mult_gate) {
      auto current_batch = mBatches.back(); // accessing last elt
      if ( current_batch->HasRoom() ) {
	current_batch->Append(mult_gate);
      } else {
	auto new_batch = std::make_shared<MultBatch>(mBatchSize);
	new_batch->Append(mult_gate);
	mBatches.emplace_back(new_batch);
      }
    }

    std::shared_ptr<MultBatch> GetMultBatch(std::size_t idx) { return mBatches[idx]; }

    // Pads the current batch if necessary
    void Close() {
      auto padding_gate = std::make_shared<PadMultGate>();
      padding_gate->UpdateParents(padding_gate, padding_gate);

      assert(padding_gate->GetMu() == FF(0));
      assert(padding_gate->GetLeft()->GetMu() == FF(0));
      assert(padding_gate->GetRight()->GetMu() == FF(0));      

      auto last_batch = mBatches.back(); // accessing last elt
      while ( last_batch->HasRoom() ) {
	last_batch->Append(padding_gate);
      } 
      // assert(last_batch->HasRoom() == false); // PASSES
    }

    // For testing purposes: sets the required preprocessing for each
    // batch to be just 0 shares
    /*void _DummyPrep(FF lambda_A, FF lambda_B, FF lambda_C) {
      for (auto batch : mBatches) batch->_DummyPrep(lambda_A, lambda_B, lambda_C);
    }*/
    /*void _DummyPrep() {
      for (auto batch : mBatches) batch->_DummyPrep();
    }*/
    /*void PrepFromDummyLambdas() {
      for (auto batch : mBatches) batch->PrepFromDummyLambdas();
    }*/

    void ClearEvaluation() {
      for (auto batch : mBatches) batch->GetClear();
    }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the layer since sometimes we just want
    // to evaluate in the clear and this won't be needed
    // To be used after the layer is closed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto batch : mBatches) batch->SetNetwork(network, id);
    }

    // For testing purposes, to avoid blocking
    /*void P1Sends() { for (auto batch : mBatches) batch->P1Sends(); }
    void PartiesReceive() { for (auto batch : mBatches) batch->PartiesReceive(); }
    void PartiesSend() { for (auto batch : mBatches) batch->PartiesSend(); }
    void P1Receives() { for (auto batch : mBatches) batch->P1Receives(); }*/
    void P1Distributes() {
      for(auto batch : mBatches) {
        batch->P1Distributes();
      }
    }
    void PartiesReceive() {
      for(auto batch : mBatches) {
        batch->PartiesReceive();
      }
    }
    void PartiesMultiplyAndSend() {
      for(auto batch : mBatches) {
        batch->PartiesMultiplyAndSend();
      }
    }
    void P1ReceivesAndReconstruct() {
      for(auto batch : mBatches) {
        batch->P1ReceivesAndReconstruct();
      }
    }
    void PartiesCheckAndComputeAuthMu(std::vector<FF> theta_list) {
      for(auto batch : mBatches) {
        batch->PartiesCheckAndComputeAuthMu(theta_list);
      }
    }
    
    // Metrics
    std::size_t GetSize() { return mBatches.size(); }

  private:
    vec<std::shared_ptr<MultBatch>> mBatches;
    std::size_t mBatchSize;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    friend class Circuit;
  };
} // namespace tp

#endif  // MULT_GATE_H
