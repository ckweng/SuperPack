#include <iostream>

#include "scl.h"

scl::NetworkConfig RunServer(int n) {
  scl::DiscoveryServer server(n);
  scl::Party party{0, "127.0.0.1", 5000};

  return server.Run(party);
}

scl::NetworkConfig RunClient(unsigned id) {
  scl::DiscoveryClient client("127.0.0.1");

  return client.Run(id, 5000 + id * 1000);
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " [id] [n]\n";
    return 0;
  }

  /* Read the players ID and the number of players to coordinate from the
   * commandline.
   */
  auto id = (unsigned)std::stoul(argv[1]);
  auto n = std::stoi(argv[2]);

  /* Player 0 is designated as the coordination server.
   */
  scl::NetworkConfig config;
  if (id == 0) {
    config = RunServer(n);
  } else {
    config = RunClient(id);
  }

  std::cout << "Discovery done!\n";
  std::cout << config.ToString() << "\n";

  /* With the received config, everyone can connect to each other and send stuff
   * around.
   */

  auto network = scl::Network::Create(config);

  for (std::size_t i = 0; i < 3; ++i) {
    // similar to the TCP channel example, send our ID to everyone:
    network.Party(i)->Send(config.Id());
    unsigned received_id;
    network.Party(i)->Recv(received_id);
    std::cout << "Received " << received_id << " from " << i << "\n";
  }
}
