/**
 * @file mem_channel.h
 *
 * SCL --- Secure Computation Library
 * Copyright (C) 2022 Anders Dalskov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 * USA
 */
#ifndef _SCL_NET_MEM_CHANNEL_H
#define _SCL_NET_MEM_CHANNEL_H

#include <array>
#include <memory>

#include "scl/net/channel.h"
#include "scl/net/shared_deque.h"

namespace scl {

/**
 * @brief Channel that communicates through in-memory buffers.
 */
class InMemoryChannel final : public Channel {
 private:
  using Buffer = details::SharedDeque<std::vector<unsigned char>>;

 public:
  /**
   * @brief Create a pair of paired channels.
   *
   * This method returns a pair of channels that shared their buffers such that
   * what is sent on one can be retrieved on the other.
   */
  static std::array<std::shared_ptr<InMemoryChannel>, 2> CreatePaired() {
    auto buf0 = std::make_shared<Buffer>();
    auto buf1 = std::make_shared<Buffer>();
    return {std::make_shared<InMemoryChannel>(buf0, buf1),
            std::make_shared<InMemoryChannel>(buf1, buf0)};
  };

  /**
   * @brief Create a channel that sends to itself.
   */
  static std::shared_ptr<InMemoryChannel> CreateSelfConnecting() {
    auto buf = std::make_shared<Buffer>();
    return std::make_shared<InMemoryChannel>(buf, buf);
  }

  /**
   * @brief Create a new channel that sends and receives on in-memory buffers.
   * @param in_buffer the buffer to read incoming messages from
   * @param out_buffer the buffer to put outgoing messages
   */
  InMemoryChannel(std::shared_ptr<Buffer> in_buffer,
                  std::shared_ptr<Buffer> out_buffer)
      : mIn(in_buffer), mOut(out_buffer){};

  /**
   * @brief Flush the incomming buffer;
   */
  void Flush() {
    while (mIn->Size()) mIn->PopFront();
    mOverflow.clear();
  };

  void Send(const unsigned char* src, std::size_t n) override;
  void Recv(unsigned char* dst, std::size_t n) override;
  void Close() override{};

 private:
  std::shared_ptr<Buffer> mIn;
  std::shared_ptr<Buffer> mOut;
  std::vector<unsigned char> mOverflow;
};

}  // namespace scl

#endif /* _SCL_NET_MEM_CHANNEL_H */
