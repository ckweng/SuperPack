#include "tp/circuits.h"

namespace tp {
  void Circuit::SetClearInputs(std::vector<std::vector<FF>> inputs_per_client) {
    if ( inputs_per_client.size() != mClients )
      throw std::invalid_argument("Number of clients do not match");

    for (std::size_t i = 0; i < mClients; i++) {
      if ( inputs_per_client[i].size() != mBatchSize*mFlatInputBatches[i].size() )
	throw std::invalid_argument("Number of inputs provided for a client does not match its number of input gates");

      for (std::size_t j = 0; j < inputs_per_client[i].size()/mBatchSize; j++) {
        for(std::size_t k = 0; k < mBatchSize; ++k) {
          mFlatInputBatches[i][j]->GetInputGate(k)->ClearInput(inputs_per_client[i][j*mBatchSize+k]);
        }
      }
    }
  }
  
  void Circuit::SetClearInputsFlat(std::vector<FF> inputs) {
    std::size_t idx(0);
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto input_batch : mFlatInputBatches[i]) {
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          input_batch->GetInputGate(j)->ClearInput(inputs[idx]);
	  idx++;
        }
      }
    }
  }
  
  std::vector<std::vector<FF>> Circuit::GetClearOutputs() {
    std::vector<std::vector<FF>> output;
    output.reserve(mClients);
    for (std::size_t i = 0; i < mClients; i++) {
      std::vector<FF> output_i;
      output_i.reserve(mBatchSize*mFlatOutputBatches[i].size());
      for (auto output_batch : mFlatOutputBatches[i]) {
        for(std::size_t j = 0; j < mBatchSize; ++j) {
          output_i.emplace_back(output_batch->GetOutputGate(j)->GetClear());
        }
      }
      output.emplace_back(output_i);
    }
    return output;
  }

  std::vector<FF> Circuit::GetClearOutputsFlat() {
    std::vector<FF> output;
    for (auto output_gate : mOutputGates) output.emplace_back(output_gate->GetClear());
    return output;
  }
} // namespace tp
