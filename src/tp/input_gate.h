#ifndef INP_GATE_H
#define INP_GATE_H

#include <vector>
#include <assert.h>

#include "gate.h"

namespace tp {
  class InputGate : public Gate {
  public:
    InputGate(std::size_t owner_id) : mOwnerID(owner_id) {
      // TODO sample at random, interactively. Get from correlator
      mIndvShrLambdaC = FF(0);
    }

    // Set cleartext inputs, for the case of cleartext evaluation
    void ClearInput(FF input) {
      mClear = input;
      mEvaluated = true;
    }

    std::size_t GetOwner() {
      return mOwnerID;
    }

    FF GetMu() {
      if ( !mLearned )
	throw std::invalid_argument("P1 hasn't learned this value yet");
      return mMu;
    }

    void SetMu(FF mu) {
      if ( mLearned )
	throw std::invalid_argument("P1 already learned this value");
      mLearned = true;
      mMu = mu;
    }

    FF GetDeltaMu() {
      if ( !mAuthLearned )
	throw std::invalid_argument("P1 hasn't learned this value yet");
      return mDeltaMu;
    }

    void SetDeltaMu(FF delta_mu) {
      if ( mAuthLearned )
	throw std::invalid_argument("P1 already learned this value");
      mDeltaMu = delta_mu;
      mAuthLearned = true;
    }

    static FF GetPackedShrDelta(std::size_t i) {
      if(mPackedShrDelta.Size() < i+1)
        throw std::invalid_argument("Packed Shares of Delta not set yet");
      return mPackedShrDelta[i];
    }

    void SetLambda(FF lambda) {
      mLambda = lambda;
      mLambdaSet = true;
    }

    FF GetDummyLambda() {
      if ( !mLambdaSet )
	throw std::invalid_argument("Lambda is not set in this input gate");
      return mLambda;
    }

    void SetIndvShrLambda(FF indv_shr) {
      mIndvShrLambdaC = indv_shr;
      mIndvShrLambdaCSet = true;
    }

    FF GetIndvShrLambda() {
      if ( !mIndvShrLambdaCSet )
	throw std::invalid_argument("IndvShrLambda is not set in this input gate");
      return mIndvShrLambdaC;
    }

    FF GetDn07Share() {
      if ( !mDn07Set )
	throw std::invalid_argument("Dn07 shares is not set in this input gate");
      return mDn07Share;
    }
    
    FF GetClear() {
      if ( !mEvaluated )
	throw std::invalid_argument("This input has not been provided yet");
      return mClear;
    };

    // Protocol-related
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    /*void _DummyPrep(FF lambda) {
      if (mID == mOwnerID) mLambda = lambda;
    }
    void _DummyPrep() {
      _DummyPrep(FF(0));
    }*/

    void SetInput(FF input) {
      if ( mID == mOwnerID ) mValue = input;
    }

    /*void OwnerSendsP1() {
      if (mID == mOwnerID) mNetwork->Party(0)->Send(mValue - mLambda);
    }

    void P1Receives() {
      if (mID == 0) {
	mNetwork->Party(mOwnerID)->Recv(mMu);
	mLearned = true;
      }
    }*/

    FF GetValue() { return mValue; }

  private:
    // ID of the party who owns this gate
    std::size_t mOwnerID;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    // Protocol-specific
    FF mLambda; // Lambda, learned by owner
    FF mValue; // Actual input, known by owner
  };

  // Used for padding batched inputs
  class PadInputGate : public InputGate {
  public:
    PadInputGate(std::size_t owner_id) : InputGate(owner_id) { mIsPadding = true; }

    FF GetMu() override { return FF(0); }
    FF GetDummyLambda() override { return FF(0); }
    FF GetIndvShrLambda() override { return FF(0); }

  private:
  };

  class InputBatch {
  public:
    InputBatch(std::size_t owner_id, std::size_t batch_size) : mOwnerID(owner_id), mBatchSize(batch_size) {
      mInputGatesPtrs.reserve(mBatchSize);
    };

    // Adds a new input_gate to the batch. It cannot add more gates than
    // the batch_size
    void Append(std::shared_ptr<InputGate> input_gate) {
      if ( mInputGatesPtrs.size() == mBatchSize )
	throw std::invalid_argument("Trying to batch more than batch_size gates");
      if ( input_gate->GetOwner() != mOwnerID )
	throw std::invalid_argument("Owner IDs do not match");
      mInputGatesPtrs.emplace_back(input_gate); }

    // For testing purposes: sets the required preprocessing for this
    // batch to be constant shares
    /*void _DummyPrep(FF lambda) {
      if ( mInputGatesPtrs.size() != mBatchSize )
	throw std::invalid_argument("The number of input gates does not match the batch size");

      mPackedShrLambda = lambda;
      for (auto input_gate : mInputGatesPtrs) input_gate->_DummyPrep(lambda);
    }
    void _DummyPrep() {
      _DummyPrep(FF(0));
    }*/

    std::size_t GetOwner() {
      return mOwnerID;
    };

                                            
    // Generates the preprocessing from the lambdas of the inputs
    /*void PrepFromDummyLambdas() {
      Vec lambda;

      for (std::size_t i = 0; i < mBatchSize; i++) {
	lambda.Emplace(mInputGatesPtrs[i]->GetDummyLambda());
      }
      // Using deg = BatchSize-1 ensures there's no randomness involved
      auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda, mBatchSize-1, mPRG);
      Vec shares = scl::details::SharesFromEvPoly(poly, mParties);

      mPackedShrLambda = shares[mID];
    }*/


    // For cleartext evaluation: calls GetClear on all its gates to
    // populate their mClear. This could return a vector with these
    // values but we're not needing them
    Vec GetClear() {
      Vec clear_val;
      for (auto gate : mInputGatesPtrs) { clear_val.Emplace(gate->GetClear()); }
      return clear_val;
    }

    Vec GetValue() {
      Vec clear_val;
      for (auto gate : mInputGatesPtrs) { clear_val.Emplace(gate->GetValue()); }
      return clear_val;
    }


    // Determines whether the batch is full
    bool HasRoom() { return mInputGatesPtrs.size() < mBatchSize; }
    
    // For fetching input gates
    std::shared_ptr<InputGate> GetInputGate(std::size_t idx) { return mInputGatesPtrs[idx]; }

    std::size_t GetNumGate() {
      return mInputGatesPtrs.size();
    }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the batch since sometimes we just want
    // to evaluate in the clear and this won't be needed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto input_gate : mInputGatesPtrs) input_gate->SetNetwork(network, id);
    }
    
    // Online protocol
    // All parties send to client the shares of masks and semi-authenticated triple
    void PartiesRevealShares();
    // Owner check the shares and distribute mu
    void ClientReconstructAndDistributes();
    // Parties receive the shares
    void PartiesReceive();
    void ClientSendMuToP1();
    void P1ReceiveMu();
    // Parties compute authenticated mu
    void PartiesAuthenticateMu();

    void SetPreprocessing(FF packed_shr_lambda,
        Vec &add_shr_delta_lambda,
        FF packed_shr_a, FF packed_shr_delta_a,
        FF packed_shr_b, FF packed_shr_delta_b,
        FF packed_shr_c, Vec &add_shr_delta_c) {
      mPackedShrLambda = packed_shr_lambda;
      mAddShrDeltaLambda = Vec(add_shr_delta_lambda.begin(),
          add_shr_delta_lambda.end());
      mPackedShrA = packed_shr_a;
      mPackedShrDeltaA = packed_shr_delta_a;
      mPackedShrB = packed_shr_b;
      mPackedShrDeltaB = packed_shr_delta_b;
      mPackedShrC = packed_shr_c;
      mAddShrDeltaC = Vec(add_shr_delta_c.begin(),
          add_shr_delta_c.end());
    }

    FF GetPackedShrLambda() { return mPackedShrLambda; }

    FF GetPackedShrDelta(std::size_t i) {
      return InputGate::GetPackedShrDelta(i);
    }

    void SetInput(std::vector<FF> inputs) {
      for(std::size_t i = 0; i < mInputGatesPtrs.size(); ++i) {
        GetInputGate(i)->SetInput(inputs[i]);
      }
    }

  private:
    // ID of the party who owns this batch
    std::size_t mOwnerID;

    std::size_t mBatchSize;


    // The input gates that are part of this batch
    vec<std::shared_ptr<InputGate>> mInputGatesPtrs;
    
    // The sharings associated to this batch
    // The authenticated lambda
    FF mPackedShrLambda;
    Vec mAddShrDeltaLambda;
    // The semi-authenticated triple
    FF mPackedShrA;
    FF mPackedShrDeltaA;
    FF mPackedShrB;
    FF mPackedShrDeltaB;
    FF mPackedShrC;
    Vec mAddShrDeltaC;
    // Output of this Batch
    FF mPackedShrDiffVAA;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    scl::PRG mPRG;
  };

  // Basically a collection of batches
  class InputLayer {
  public:
    InputLayer(std::size_t owner_id, std::size_t batch_size) : mOwnerID(owner_id), mBatchSize(batch_size) {
      auto first_batch = std::make_shared<InputBatch>(mOwnerID, mBatchSize);
      // Append a first batch
      mBatches.emplace_back(first_batch);
      mBatches.pop_back(); // TODO initialize empty
    };

    // Adds a new input_gate to the layer. It checks if the current
    // batch is full and if so creates a new one.
    void Append(std::shared_ptr<InputBatch> new_batch) {
      mBatches.emplace_back(new_batch);
    }

    std::shared_ptr<InputBatch> GetInputBatch(std::size_t idx) { return mBatches[idx]; }

    // Pads the current batch if necessary
    void Close() {
      auto padding_gate = std::make_shared<PadInputGate>(mOwnerID);
      assert(padding_gate->GetMu() == FF(0));

      if(mBatches.size() == 0) { return; }
      auto last_batch = mBatches.back(); // accessing last elt
      while ( last_batch->HasRoom() ) {
	last_batch->Append(padding_gate);
      } 
      // assert(last_batch->HasRoom() == false); // PASSES
    }

    // For testing purposes: sets the required preprocessing for each
    // batch to be just 0 shares
    /*void _DummyPrep(FF lambda, FF delta) {
      for (auto batch : mBatches) batch->_DummyPrep(lambda);
    }
    void _DummyPrep() {
      for (auto batch : mBatches) batch->_DummyPrep();
    }*/
    /*void PrepFromDummyLambdas() {
      for (auto batch : mBatches) batch->PrepFromDummyLambdas();
    }*/

    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto batch : mBatches) batch->SetNetwork(network, id);
    }

    void ClearEvaluation() {
      for (auto batch : mBatches) batch->GetClear();
    }

    // Metrics
    std::size_t GetSize() { return mBatches.size(); }

  private:
    std::size_t mOwnerID;
    vec<std::shared_ptr<InputBatch>> mBatches;
    std::size_t mBatchSize;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    friend class Circuit;
  };

} // namespace tp

#endif  // INP_GATE_H
