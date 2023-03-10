#include "scl/net/discovery/server.h"

#include <memory>
#include <vector>

#include "scl/net/tcp_utils.h"

using Server = scl::DiscoveryServer;

scl::NetworkConfig Server::Run(const scl::Party& me) {
  // one of the parties is us, which we do not connect to.
  auto ssock = scl::details::CreateServerSocket(mPort, mNumberOfParties - 1);
  std::vector<std::shared_ptr<scl::Channel>> channels;
  std::vector<std::string> hostnames;

  for (std::size_t i = 0; i < mNumberOfParties; ++i) {
    if (i == me.id) {
      channels.emplace_back(nullptr);
      hostnames.emplace_back(me.hostname);
    } else {
      auto ac = scl::details::AcceptConnection(ssock);
      auto hostname = scl::details::GetAddress(ac);
      auto channel = std::make_shared<scl::TcpChannel>(ac.socket);
      channels.emplace_back(channel);
      hostnames.emplace_back(hostname);
    }
  }

  Server::CollectIdsAndPorts discovery(hostnames);
  Server::Ctx ctx{me, scl::Network{channels}};
  return scl::Evaluate(discovery, ctx);
}

Server::SendNetworkConfig Server::CollectIdsAndPorts::Run(Server::Ctx& ctx) {
  auto my_id = ctx.me.id;
  std::vector<scl::Party> parties(mHostnames.size());
  parties[my_id] = ctx.me;

  for (std::size_t i = 0; i < mHostnames.size(); ++i) {
    if (my_id != i) {
      unsigned id;
      ctx.network.Party(i)->Recv(id);
      if (id >= parties.size()) {
        throw std::logic_error("received invalid party ID");
      }
      int port;
      ctx.network.Party(i)->Recv(port);
      parties[id] = scl::Party{id, mHostnames[id], port};
    }
  }

  scl::NetworkConfig cfg(ctx.me.id, parties);
  return Server::SendNetworkConfig(cfg);
}

static inline void SendHostname(scl::Channel* channel, std::string hostname) {
  std::size_t len = hostname.size();
  const unsigned char* ptr =
      reinterpret_cast<const unsigned char*>(hostname.c_str());

  channel->Send(len);
  channel->Send(ptr, len);
}

static inline void SendConfig(scl::Channel* channel,
                              const scl::NetworkConfig& config) {
  channel->Send(config.NetworkSize());
  for (std::size_t i = 0; i < config.NetworkSize(); ++i) {
    auto party = config.Parties()[i];
    channel->Send(party.id);
    channel->Send(party.port);
    SendHostname(channel, party.hostname);
  }
}

scl::NetworkConfig Server::SendNetworkConfig::Finalize(Server::Ctx& ctx) {
  std::size_t network_size = mConfig.NetworkSize();
  for (std::size_t i = 0; i < network_size; ++i) {
    if (i == mConfig.Id()) continue;

    auto channel = ctx.network.Party(i);
    SendConfig(channel, mConfig);
  }

  return mConfig;
}
