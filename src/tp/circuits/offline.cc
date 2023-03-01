#include "tp/circuits.h"

namespace tp {
  void Circuit::_DummyPrep(FF lambda) {
    // fix the seed of PRG
    if(scl::PRG::BlockSize() != 16) {
      throw std::invalid_argument("The seeds are not 16 bytes\n");
    }
    unsigned char seed[16] = {0x00, 0x01, 0x02, 0x03,
      0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f};
    scl::PRG prg(seed);
    // TODO generate shares of Gate::mPackedShrDelta
    FF delta = FF::Random(prg);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      FF x_point = FF(-i);
      auto poly_delta = scl::details::EvPolyFromSecretAndPointAndDegree(
          delta, x_point, mParties-2*mBatchSize+1, prg);
      Vec shares_delta = scl::details::SharesFromEvPoly(poly_delta, mParties);
      Gate::SetPackedShrDelta(shares_delta[mID]);
    }

    for (auto input_gate : mInputGates) {
      input_gate->SetLambda(lambda);
    }
    for (std::size_t layer = 0; layer < GetDepth(); layer++) {
      for (auto mult_gate : mFlatMultLayers[layer]) {
	mult_gate->SetDummyLambda(lambda);
      }
    }
    for (auto output_gate : mOutputGates) {
      output_gate->GetDummyLambda(); //populate outputs and add wires
    }

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

    // For each batch of input/output gates
    // Set degree-(n-1) packed shares lambda
    // Set additive shares for Delta * (lambda_1,...,lambda_n)
    // Set preprocessed semi-authenticated Beaver triple
    for (auto input_layer : mInputLayers) {
      for (auto input_batch : input_layer.mBatches) {
	Vec lambda_A;
	lambda_A.Reserve(mBatchSize);
	for (std::size_t i = 0; i < mBatchSize; i++) {
          FF get_lambda = input_batch->GetInputGate(i)->GetDummyLambda();
	  lambda_A.Emplace(get_lambda);
	}
	auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mParties-1, prg);
	Vec new_shares = scl::details::SharesFromEvPoly(poly, mParties);

	if(mID == 0) {
          Vec delta_lambda = lambda_A.ScalarMultiply(delta);
	  input_batch->SetPreprocessing(new_shares[mID], delta_lambda,
              triple_a_shares[mID], triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], triple_c_delta);
        } else {
          input_batch->SetPreprocessing(new_shares[mID], dummy_vec,
              triple_a_shares[mID], triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], dummy_vec);
        }
      }
    }

    for (auto output_layer : mOutputLayers) {
      for (auto output_batch : output_layer.mBatches) {
	Vec lambda_A;
	lambda_A.Reserve(mBatchSize);
	for (std::size_t i = 0; i < mBatchSize; i++) {
	  lambda_A.Emplace(output_batch->GetOutputGate(i)->GetDummyLambda());
	}
	auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mParties-1, prg);
	Vec new_shares = scl::details::SharesFromEvPoly(poly, mParties);

	if(mID == 0) {
          Vec delta_lambda = lambda_A.ScalarMultiply(delta);
	  output_batch->SetPreprocessing(new_shares[mID], delta_lambda,
              triple_a_shares[mID], triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], triple_c_delta);
        } else {
          output_batch->SetPreprocessing(new_shares[mID], dummy_vec,
              triple_a_shares[mID], triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], dummy_vec);
        }
      }
    }

    // Mults
    for (auto mult_layer : mMultLayers) {
      for (auto mult_batch : mult_layer.mBatches) {
        // lambda_alpha, lambda_beta, lambda_gamma
	Vec lambda_A;
	Vec lambda_B;
	Vec lambda_C;
	lambda_A.Reserve(mBatchSize);
	lambda_B.Reserve(mBatchSize);
	lambda_C.Reserve(mBatchSize);
	for (std::size_t i = 0; i < mBatchSize; i++) {
	  lambda_A.Emplace(mult_batch->GetMultGate(i)->GetLeft()->GetDummyLambda());
	  lambda_B.Emplace(mult_batch->GetMultGate(i)->GetRight()->GetDummyLambda());
	  lambda_C.Emplace(mult_batch->GetMultGate(i)->GetDummyLambda());
	}
        // delta*(lambda_alpha - a), delta*(lambda_beta - b)
        Vec diff_lambda_AA = lambda_A.Subtract(triple_a);
        Vec delta_diff_lambda_AA = diff_lambda_AA.ScalarMultiply(delta);
        Vec diff_lambda_BB = lambda_B.Subtract(triple_b);
        Vec delta_diff_lambda_BB = diff_lambda_BB.ScalarMultiply(delta);
        Vec delta_lambda_C = lambda_C.ScalarMultiply(delta);
        // interpolate polynomials
	auto poly_C = scl::details::EvPolyFromSecretsAndDegree(lambda_C, mParties-1, prg);
	Vec shares_C = scl::details::SharesFromEvPoly(poly_C, mParties);

        // set values
        if(mID == 0) {
          mult_batch->SetPreprocessing(delta_diff_lambda_AA,
              delta_diff_lambda_BB, shares_C[mID], delta_lambda_C);
          mult_batch->SetPreprocessingTriple(triple_a_shares[mID],
              triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], triple_c_delta);
          mult_batch->SetP1Preprocessing(diff_lambda_AA, diff_lambda_BB);
        } else {
          mult_batch->SetPreprocessing(dummy_vec,
              dummy_vec, shares_C[mID], dummy_vec);
          mult_batch->SetPreprocessingTriple(triple_a_shares[mID],
              triple_a_delta_shares[mID],
              triple_b_shares[mID], triple_b_delta_shares[mID],
              triple_c_shares[mID], dummy_vec);
        }
      }
    }
      
      
  }

  // Dummy prep of function-independent preprocessing
  // order of invocation
  // 1. GenCorrelator();
  // 2. _DummyPrepFI
  // 3. MapCorrToCircuit()
  // 4. other function-dependent preprocessing functions
  void Circuit::_DummyPrepFI(FF lambda) {
    // fix the seed of PRG
    if(scl::PRG::BlockSize() != 16) {
      throw std::invalid_argument("The seeds are not 16 bytes\n");
    }
    unsigned char seed[16] = {0x00, 0x01, 0x02, 0x03,
      0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f};
    scl::PRG prg(seed);
    // TODO generate shares of Gate::mPackedShrDelta
    FF delta = FF::Random(prg);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      FF x_point = FF(-i);
      auto poly_delta = scl::details::EvPolyFromSecretAndPointAndDegree(
          delta, x_point, mParties-2*mBatchSize+1, prg);
      Vec shares_delta = scl::details::SharesFromEvPoly(poly_delta, mParties);
      Gate::SetPackedShrDelta(shares_delta[mID]);
    }
    mCorrelator.GenPrepDummy(lambda, delta, prg);
  }


    // Populates each batch with dummy preprocessing (all zeros)
  /*void Circuit::_DummyPrep() {
      for (auto input_layer : mInputLayers) input_layer._DummyPrep();
      for (auto mult_layer : mMultLayers) mult_layer._DummyPrep();
      for (auto output_layer : mOutputLayers) output_layer._DummyPrep();
    }*/

    // Set all settable lambdas to a constant
    void Circuit::SetDummyLambdas(FF lambda) {
      // Output wires of multiplications
      for (std::size_t layer = 0; layer < GetDepth(); layer++) {
	for (auto mult_gate : mFlatMultLayers[layer]) {
	  mult_gate->SetDummyLambda(lambda);
	}
      }
      // Output wires of input gates
      for (auto input_gate : mInputGates) {
	input_gate->SetLambda(lambda);
      }
    }

    // Populates the lambdas of addition and output gates
    void Circuit::PopulateDummyLambdas() {
      // Addition gates
      for (auto add_gate : mAddGates) {
	(void)add_gate->GetDummyLambda();
      }
      // Output gates
      for (auto output_gate : mOutputGates) {
	(void)output_gate->GetDummyLambda();
      }      
    }

    // Sets the lambdas for the output wires of addition and output
    // gates based on the lambdas for multiplications and input gates
    /*void Circuit::PrepFromDummyLambdas() {
      for (auto input_layer : mInputLayers) input_layer.PrepFromDummyLambdas();
      for (auto output_layer : mOutputLayers) output_layer.PrepFromDummyLambdas();
      for (auto mult_layer : mMultLayers) mult_layer.PrepFromDummyLambdas();
    }*/

    void Circuit::GenCorrelator() {
      if ( !mIsNetworkSet )
	throw std::invalid_argument("Cannot set correlator without setting a network first");
      
      std::size_t n_ind_shares = GetNInputs() + GetSize();
      std::size_t n_mult_batches = GetNMultBatches();
      std::size_t n_inout_batches = GetNInputBatches() + GetNOutputBatches();

      mCorrelator = Correlator(n_ind_shares, n_mult_batches, n_inout_batches, mBatchSize);
      mCorrelator.SetNetwork(mNetwork, mID);
      mCorrelator.PrecomputeEi();
    }

    // Populates the mappings in the correlator so that the
    // F.I. preprocessing is mapped to different gates/batches
    void Circuit::MapCorrToCircuit() {
      // Individual shares
      for (auto input_gate : mInputGates) {
	mCorrelator.PopulateIndvShrs(input_gate);
      }
      for (std::size_t layer = 0; layer < GetDepth(); layer++) {
	for (auto mult_gate : mFlatMultLayers[layer]) {
	  mCorrelator.PopulateIndvShrs(mult_gate);
	}
      }
      for (auto add_gate : mAddGates) {
	mCorrelator.PopulateIndvShrs(add_gate);
      }
      for (auto output_gate : mOutputGates) {
	mCorrelator.PopulateIndvShrs(output_gate);
      }

      // Input batches
      for (auto input_layer : mInputLayers) {
	for (auto input_batch : input_layer.mBatches) {
	  mCorrelator.PopulateInputBatches(input_batch);
	}
      }
      // Output batches
      for (auto output_layer : mOutputLayers) {
	for (auto output_batch : output_layer.mBatches) {
	  mCorrelator.PopulateOutputBatches(output_batch);
	}
      }
      // Mult batches
      for (auto mult_layer : mMultLayers) {
	for (auto mult_batch : mult_layer.mBatches) {
	  mCorrelator.PopulateMultBatches(mult_batch);
	}
      }
    }

    // Prep inputs & outputs
    void Circuit::PrepMultPartiesSendP1() {
      for (auto mult_layer : mMultLayers) {
	for (auto mult_batch : mult_layer.mBatches) {
	  mCorrelator.PrepMultPartiesSendP1(mult_batch);
	}
      }
    }
    void Circuit::PrepMultP1Receives() {
      for (auto mult_layer : mMultLayers) {
	for (auto mult_batch : mult_layer.mBatches) {
	  (void) mult_batch;
	  mCorrelator.PrepMultP1Receives(mult_batch);
	}
      }      
    }

    void Circuit::PrepIO() {
      for (auto input_layer : mInputLayers) {
	for (auto input_batch : input_layer.mBatches) {
	  mCorrelator.PrepInput(input_batch);
	}
      }
      for (auto output_layer : mOutputLayers) {
	for (auto output_batch : output_layer.mBatches) {
	  mCorrelator.PrepOutput(output_batch);
	}
      }
    }

} // namespace tp
