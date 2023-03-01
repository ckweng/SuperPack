#include <iostream>
#include <chrono>
#include <thread>

#include "tp/circuits.h"
#include "misc.h"

#define DELIM std::cout << "========================================\n"

#define DEBUG true
#define THREAD true

#define PRINT(x) if (DEBUG) std::cout << x << "\n";

inline std::size_t ValidateN(const std::size_t n) {
  assert(n > 3
	 && ((n % 16) == 0)
	 );
  return n;
}

inline std::size_t ValidateId(const std::size_t id, const std::size_t n) {
      if ( id >= n )
	throw std::invalid_argument("ID cannot be larger than number of parties");
  return id;
}

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cout << "usage: " << argv[0] << " [N] [id] [size] [depth] [percentage of corruption]\n";
    return 0;
  }

  std::size_t n = ValidateN(std::stoul(argv[1]));
  std::size_t id = ValidateId(std::stoul(argv[2]), n);
  std::size_t size = std::stoul(argv[3]);
  std::size_t depth = std::stoul(argv[4]);
  std::size_t width = size/depth;
  float curruption_percentage = std::stof(argv[5]);
  std::size_t t = n * curruption_percentage;
  if((t % 2) == 0) {
    t--;
  }

  std::size_t batch_size = (n+1-t)/2;
  std::size_t n_parties = t + 2*(batch_size - 1) + 1;

  DELIM;
  std::cout << "Running benchmark with N " << n << ", size " <<
    size << ", width " << width << " and depth " << depth << "\n";
  std::cout << "With corruption " << t << " and batch size " <<
    batch_size << std::endl;
  DELIM;

  auto config = scl::NetworkConfig::Load(id, "./client_info.txt");
  // std::cout << "Config:"
  //           << "\n";
  // for (const auto& party : config.Parties()) {
  //   std::cout << " -- " << party.id << ", " << party.hostname << ", "
  //             << party.port << "\n";
  // }

  std::cout << "Connecting ..."
            << "\n";
  auto network = scl::Network::Create(config);

  std::cout << "Done!\n";

  // std::size_t n_clients = n_parties;

  tp::CircuitConfig circuit_config;
  circuit_config.n_parties = n_parties;
  circuit_config.inp_gates = std::vector<std::size_t>(n_parties, 0);
  circuit_config.inp_gates[0] = 2;
  circuit_config.out_gates = std::vector<std::size_t>(n_parties, 0);
  circuit_config.out_gates[0] = batch_size;
  circuit_config.width = width;
  circuit_config.depth = depth;
  circuit_config.batch_size = batch_size;
    
  auto circuit = tp::Circuit::FromConfig(circuit_config);

  circuit.SetNetwork(std::make_shared<scl::Network>(network), id);

  circuit.GenCorrelator();
  circuit.SetThreshold(t);

  DELIM;
  std::cout << "Running function-independent preprocessing\n";
  
  START_TIMER(fi_prep);

  // only with dummy circuit-independent preprocessing
  PRINT("fi_prep");
  //circuit.FIPrepFromDummyOle();
  circuit._DummyPrepFI(tp::FF(1));

  STOP_TIMER(fi_prep);

  circuit.MapCorrToCircuit(); 

  DELIM;
  std::cout << "Running function-dependent preprocessing\n";
  
  START_TIMER(fd_prep);
  PRINT("fd_prep");
  circuit.PrepMultPartiesSendP1(); 
  circuit.PrepMultP1Receives(); 
  circuit.PrepIO(); 

  STOP_TIMER(fd_prep);

  std::vector<tp::FF> result;
  if (id == 0) {
    std::vector<tp::FF> inputs(batch_size, tp::FF(0432432));
    circuit.SetClearInputsFlat(inputs);
    result = circuit.GetClearOutputsFlat();
    circuit.SetInputs(std::vector<std::vector<tp::FF>>{inputs, inputs});
  }

  DELIM;
  std::cout << "Running online phase\n";
  START_TIMER(online);
  // INPUT
  PRINT("Input");

  circuit.InputPartiesRevealShares();
  circuit.InputClientReconstructAndDistributes();
  circuit.InputPartiesReceive();
  circuit.InputClientSendMuToP1();
  circuit.InputP1ReceiveMu();
  circuit.InputPartiesAuthenticateMu();


  // MULT
  for (std::size_t layer = 0; layer < circuit_config.depth; layer++) {
    // PRINT("Mult layer " << layer);
    circuit.MultP1Distribute(layer);
    circuit.MultPartiesReceive(layer);
    circuit.MultPartiesMultiplyAndSend(layer);
    circuit.MultP1ReceivesAndReconstruct(layer);
    circuit.MultPartiesCheckAndComputeAuthMu(layer);
  }

  // OUTPUT
  PRINT("Output");
  circuit.OutputPartiesRevealShares();
  circuit.OutputP1ReconstructAndDistributes();
  circuit.OutputPartiesReceiveShares();
  circuit.OutputPartiesVerify();
  circuit.OutputPartiesCheckTheta();
  circuit.OutputPartiesRevealOutputVal();
  circuit.OutputClientReconstructSecret();
  STOP_TIMER(online);
  
  // Check output
  if (id == 0) {
    // std::cout << "\nOUTPUT = " << circuit.GetOutputs()[0] << "\nREAL = " << result[0] << "\n";
    assert( circuit.GetOutputs() == result );
  }

  std::cout << "\nclosing the network ...\n";
  // technically not necessary as channels are destroyed when their dtor is
  // called.
  network.Close();

}
