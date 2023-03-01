#include <catch2/catch.hpp>
#include <thread>

#include "scl/net/tcp_channel.h"
#include "scl/net/tcp_utils.h"
#include "scl/prg.h"

static int test_port = 14421;

static inline int GetPort() { return test_port++; }

static inline bool Eq(const unsigned char* p, const unsigned char* q, int n) {
  while (n-- > 0 && *p++ == *q++)
    ;
  return n < 0;
}

TEST_CASE("TcpChannel", "[network]") {
  SECTION("Connect and Close") {
    auto port = GetPort();

    std::shared_ptr<scl::TcpChannel> client, server;
    std::thread clt([&]() {
      int socket = scl::details::ConnectAsClient("0.0.0.0", port);
      client = std::make_shared<scl::TcpChannel>(socket);
    });
    std::thread srv([&]() {
      int ssock = scl::details::CreateServerSocket(port, 1);
      auto ac = scl::details::AcceptConnection(ssock);
      server = std::make_shared<scl::TcpChannel>(ac.socket);
      scl::details::CloseSocket(ssock);
    });

    clt.join();
    srv.join();

    REQUIRE(client->Alive());
    REQUIRE(server->Alive());

    client->Close();
    server->Close();

    REQUIRE(!server->Alive());
    REQUIRE(!client->Alive());
  }

  SECTION("Send Receive") {
    auto port = GetPort();

    std::shared_ptr<scl::TcpChannel> client, server;
    std::thread clt([&]() {
      int socket = scl::details::ConnectAsClient("0.0.0.0", port);
      client = std::make_shared<scl::TcpChannel>(socket);
    });

    std::thread srv([&]() {
      int ssock = scl::details::CreateServerSocket(port, 1);
      auto ac = scl::details::AcceptConnection(ssock);
      server = std::make_shared<scl::TcpChannel>(ac.socket);
      scl::details::CloseSocket(ssock);
    });

    clt.join();
    srv.join();

    scl::PRG prg;
    unsigned char send[200] = {0};
    unsigned char recv[200] = {0};
    prg.Next(send, 200);

    client->Send(send, 100);
    client->Send(send + 100, 100);
    server->Recv(recv, 20);
    server->Recv(recv + 20, 180);

    REQUIRE(Eq(send, recv, 200));
  }
}
