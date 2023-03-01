#include "tp/circuits.h"

namespace tp {

  std::vector<std::vector<std::vector<FF>>> EvPolynomialInternal::mLagrangePolynomialPos;
  std::vector<std::vector<std::vector<FF>>> EvPolynomialInternal::mLagrangePolynomialNeg;

  void Circuit::Init(std::size_t n_clients, std::size_t batch_size) {

    std::vector<std::size_t> degree_set = {batch_size-1, 2*(batch_size-1), n_clients-2*batch_size+1, n_clients-batch_size, n_clients-batch_size+1, n_clients-1};
    EvPolynomialInternal::PreprocessLagrangePolynomial(batch_size, degree_set, n_clients-1);

    // Initialize input layer
    mInputLayers.reserve(n_clients);
    mFlatInputBatches.resize(n_clients);
    for (std::size_t i = 0; i < n_clients; i++) {
      auto input_layer = InputLayer(i, batch_size);
      mInputLayers.emplace_back(input_layer);
    }
      
    // Initialize output layer
    mOutputLayers.reserve(n_clients);
    mFlatOutputBatches.resize(n_clients);
    for (std::size_t i = 0; i < n_clients; i++) {
      auto output_layer = OutputLayer(i, batch_size);
      mOutputLayers.emplace_back(output_layer);
    }

    // Initialize first mult layer
    auto first_layer = MultLayer(mBatchSize);
    mMultLayers.emplace_back(first_layer);
    mFlatMultLayers.emplace_back(VecMultGates());

    mIsClosed = false;
    mIsNetworkSet = false;
  }

  std::shared_ptr<InputBatch> Circuit::Input(std::size_t owner_id) {
    auto new_batch = std::make_shared<InputBatch>(owner_id, mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      auto input_gate = std::make_shared<tp::InputGate>(owner_id);
      mInputGates.emplace_back(input_gate);
      new_batch->Append(input_gate);
    }
    mInputLayers[owner_id].Append(new_batch);
    mFlatInputBatches[owner_id].emplace_back(new_batch);
    return new_batch;
  }

  std::shared_ptr<OutputBatch> Circuit::Output(std::size_t owner_id, std::vector<std::shared_ptr<Gate>> output) {
    auto new_batch = std::make_shared<OutputBatch>(owner_id, mBatchSize);
    for(std::size_t i = 0; i < mBatchSize; ++i) {
      auto output_gate = std::make_shared<tp::OutputGate>(owner_id, output[i]);
      mOutputGates.emplace_back(output_gate);
      new_batch->Append(output_gate);
    }
    mOutputLayers[owner_id].Append(new_batch);
    mFlatOutputBatches[owner_id].emplace_back(new_batch);
    return new_batch;
  }

  std::shared_ptr<MultGate> Circuit::Mult(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
    auto mult_gate = std::make_shared<tp::MultGate>(left, right);
    mMultLayers.back().Append(mult_gate);
    mFlatMultLayers.back().emplace_back(mult_gate);
    mSize++;
    return mult_gate;
  }

  void Circuit::NewLayer() {
    // Update width
    if (mWidth < mFlatMultLayers.back().size()) mWidth = mFlatMultLayers.back().size();
    mFlatMultLayers.emplace_back(VecMultGates());
      
    // Pad the batches if necessary
    mMultLayers.back().Close();
    // Open space for the next layer
    auto next_layer = MultLayer(mBatchSize);
    mMultLayers.emplace_back(next_layer);
  }

  void Circuit::LastLayer() {
    // Update width
    if (mWidth < mFlatMultLayers.back().size()) mWidth = mFlatMultLayers.back().size();
    mFlatMultLayers.emplace_back(VecMultGates());

    mMultLayers.back().Close();
  }

} // namespace tp
