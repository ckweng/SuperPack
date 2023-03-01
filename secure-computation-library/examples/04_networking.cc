#include <iostream>

#include "scl.h"

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " [id] [n]\n";
    return 0;
  }

  auto id = (unsigned)std::stoul(argv[1]);
  auto n = std::stoi(argv[2]);

  std::cout << "About to start a TCP network where I am party " << id
            << " and with a total of " << n << " parties\n";

  auto config = scl::NetworkConfig::Localhost(id, n);
  std::cout << "Config:"
            << "\n";
  for (const auto& party : config.Parties()) {
    std::cout << " -- " << party.id << ", " << party.hostname << ", "
              << party.port << "\n";
  }

  std::cout << "Connecting ..."
            << "\n";
  auto network = scl::Network::Create(config);

  std::cout << "Done!\n";

  std::cout << "Sending my ID to other everyone!\n";

  for (std::size_t i = 0; i < network.Size(); ++i) {
    network.Party(i)->Send(id);
  }

  unsigned received_id;
  for (std::size_t i = 0; i < network.Size(); ++i) {
    network.Party(i)->Recv(received_id);
    std::cout << "received " << received_id << " from " << i << "\n";
  }

  std::cout << "closing the network ...\n";
  // technically not necessary as channels are destroyed when their dtor is
  // called.
  network.Close();
}
