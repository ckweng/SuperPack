#include "scl/prg.h"

#include <cstring>
#include <stdexcept>

/* https://github.com/sebastien-riou/aes-brute-force */

using byte_t = unsigned char;
using block_t = __m128i;

using std::size_t;
using std::vector;

#define DO_ENC_BLOCK(m, k)              \
  do {                                  \
    m = _mm_xor_si128(m, k[0]);         \
    m = _mm_aesenc_si128(m, k[1]);      \
    m = _mm_aesenc_si128(m, k[2]);      \
    m = _mm_aesenc_si128(m, k[3]);      \
    m = _mm_aesenc_si128(m, k[4]);      \
    m = _mm_aesenc_si128(m, k[5]);      \
    m = _mm_aesenc_si128(m, k[6]);      \
    m = _mm_aesenc_si128(m, k[7]);      \
    m = _mm_aesenc_si128(m, k[8]);      \
    m = _mm_aesenc_si128(m, k[9]);      \
    m = _mm_aesenclast_si128(m, k[10]); \
  } while (0)

#define AES_128_key_exp(k, rcon) \
  aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

inline static block_t aes_128_key_expansion(block_t key, block_t keygened) {
  keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3, 3, 3, 3));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, keygened);
}

inline static void aes128_load_key(byte_t* enc_key, block_t* key_schedule) {
  key_schedule[0] = _mm_loadu_si128((const block_t*)enc_key);
  key_schedule[1] = AES_128_key_exp(key_schedule[0], 0x01);
  key_schedule[2] = AES_128_key_exp(key_schedule[1], 0x02);
  key_schedule[3] = AES_128_key_exp(key_schedule[2], 0x04);
  key_schedule[4] = AES_128_key_exp(key_schedule[3], 0x08);
  key_schedule[5] = AES_128_key_exp(key_schedule[4], 0x10);
  key_schedule[6] = AES_128_key_exp(key_schedule[5], 0x20);
  key_schedule[7] = AES_128_key_exp(key_schedule[6], 0x40);
  key_schedule[8] = AES_128_key_exp(key_schedule[7], 0x80);
  key_schedule[9] = AES_128_key_exp(key_schedule[8], 0x1B);
  key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);
}

inline static void aes128_enc(block_t* key_schedule, byte_t* pt, byte_t* ct) {
  block_t m = _mm_loadu_si128((block_t*)pt);
  DO_ENC_BLOCK(m, key_schedule);
  _mm_storeu_si128((block_t*)ct, m);
}

scl::PRG::PRG() { Init(); }

scl::PRG::PRG(const unsigned char* seed) {
  memcpy(mSeed, seed, SeedSize());
  Init();
}

void scl::PRG::Update() { mCounter += 1; }

void scl::PRG::Init() { aes128_load_key(mSeed, mState); }

void scl::PRG::Reset() {
  Init();
  mCounter = PRG_INITIAL_COUNTER;
}

static inline auto create_mask(const long counter) {
  return _mm_set_epi64x(PRG_NONCE, counter);
}

void scl::PRG::Next(byte_t* dest, size_t nbytes) {
  if (!nbytes) return;

  size_t nblocks = nbytes / BlockSize();

  if (nbytes % BlockSize()) nblocks++;

  block_t mask = create_mask(mCounter);
  byte_t* out = (byte_t*)malloc(nblocks * BlockSize());
  byte_t* p = out;

  // LCOV_EXCL_START
  if (!out) throw std::runtime_error("Could not allocate memory for PRG.");
  // LCOV_EXCL_STOP

  for (size_t i = 0; i < nblocks; i++) {
    aes128_enc(mState, (byte_t*)(&mask), p);
    Update();
    mask = create_mask(mCounter);
    p += BlockSize();
  }

  memcpy(dest, out, nbytes);
  free(out);
}
