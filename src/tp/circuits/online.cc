#include "tp/circuits.h"

namespace tp {
  void Circuit::SetInputs(std::vector<std::vector<FF>> inputs) {
    if ( inputs.size() != mFlatInputBatches[mID].size() )
      throw std::invalid_argument("Number of inputs do not match");
    // Set input and send to P1
    for (std::size_t i = 0; i < inputs.size(); i++) {
      mFlatInputBatches[mID][i]->SetInput(inputs[i]);
    }
  }

  // Input protocol
  void Circuit::InputPartiesRevealShares() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto input_batch : mFlatInputBatches[i]) {
	input_batch->PartiesRevealShares();
      }
    }
  }
  void Circuit::InputClientReconstructAndDistributes() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto input_batch : mFlatInputBatches[i]) {
	input_batch->ClientReconstructAndDistributes();
      }
    }
  }
  void Circuit::InputPartiesReceive() {
    for(std::size_t i = 0; i < mClients; ++i) {
      for(auto input_batch : mFlatInputBatches[i]) {
        input_batch->PartiesReceive();
      }
    }
  }
  void Circuit::InputClientSendMuToP1() {
    for(std::size_t i = 0; i < mClients; ++i) {
      for(auto input_batch : mFlatInputBatches[i]) {
        input_batch->ClientSendMuToP1();
      }
    }
  }
  void Circuit::InputP1ReceiveMu() {
    for(std::size_t i = 0; i < mClients; ++i) {
      for(auto input_batch : mFlatInputBatches[i]) {
        input_batch->P1ReceiveMu();
      }
    }
  }
  void Circuit::InputPartiesAuthenticateMu() {
    for(std::size_t i = 0; i < mClients; ++i) {
      for(auto input_batch : mFlatInputBatches[i]) {
        input_batch->PartiesAuthenticateMu();
      }
    }
  }

  void Circuit::RunInput() {
    InputPartiesRevealShares();
    InputClientReconstructAndDistributes();
    InputPartiesReceive();
    InputClientSendMuToP1();
    InputP1ReceiveMu();
    InputPartiesAuthenticateMu();
  }

  // Multiplications in the i-th layer
  void Circuit::MultP1Distribute(std::size_t layer) {
    mMultLayers[layer].P1Distributes();
  }
  void Circuit::MultPartiesReceive(std::size_t layer) {
    mMultLayers[layer].PartiesReceive();
  }
  void Circuit::MultPartiesMultiplyAndSend(std::size_t layer) {
    mMultLayers[layer].PartiesMultiplyAndSend();
  }
  void Circuit::MultP1ReceivesAndReconstruct(std::size_t layer) {
    mMultLayers[layer].P1ReceivesAndReconstruct();
  }
  void Circuit::MultPartiesCheckAndComputeAuthMu(std::size_t layer) {
    mMultLayers[layer].PartiesCheckAndComputeAuthMu(mTheta);
  }

  void Circuit::RunMult(std::size_t layer) {
    MultP1Distribute(layer);
    MultPartiesReceive(layer);
    MultPartiesMultiplyAndSend(layer);
    MultP1ReceivesAndReconstruct(layer);
    MultPartiesCheckAndComputeAuthMu(layer);
  }

  // Output layers
  void Circuit::OutputPartiesRevealShares() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->PartiesRevealShares();
      }
    }
  }
  void Circuit::OutputP1ReconstructAndDistributes() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->P1ReconstructAndDistributes();
      }
    }
  }
  void Circuit::OutputPartiesReceiveShares() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->PartiesReceiveShares();
      }
    }
  }
  void Circuit::OutputPartiesVerify() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->PartiesVerify(mTheta);
        mTheta.clear();
      }
    }
  }
  void Circuit::OutputPartiesCheckTheta() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->PartiesCheckTheta();
      }
    }
  }
  void Circuit::OutputPartiesRevealOutputVal() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->PartiesRevealOutputVal();
      }
    }
  }
  void Circuit::OutputClientReconstructSecret() {
    for (std::size_t i = 0; i < mClients; i++) {
      for (auto output_batch : mFlatOutputBatches[i]) {
	output_batch->ClientReconstructSecret();
      }
    }
  }

  void Circuit::RunOutput() {
    OutputPartiesRevealShares();
    OutputP1ReconstructAndDistributes();
    OutputPartiesReceiveShares();
    OutputPartiesVerify();
    OutputPartiesRevealOutputVal();
    OutputClientReconstructSecret();
  }

  void Circuit::RunProtocol() {
    RunInput();
    for (std::size_t layer = 0; layer < mMultLayers.size(); layer++) {
      RunMult(layer);
    }
    RunOutput();
  }
    
  // Returns a vector with the outputs after computation
  std::vector<FF> Circuit::GetOutputs() {
    std::vector<FF> output;
    output.reserve(mBatchSize*mFlatOutputBatches[mID].size());
    for (auto output_batches : mFlatOutputBatches[mID]) {
      for(size_t i = 0; i < mBatchSize; ++i) {
        output.emplace_back(output_batches->GetOutputGate(i)->GetValue());
      }
    }
    return output;
  }
} // namespace tp
