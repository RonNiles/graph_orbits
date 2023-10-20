#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <unordered_map>

/*
 * Each bit element represents an edge
 */
static struct BitElement {
  uint8_t from_node : 4;
  uint8_t to_node : 4;
} bit_elements[64];

static void MakeBitElements() {
  unsigned next_element = 0;
  for (uint8_t n = 1;; ++n) {
    for (uint8_t k = 0; k < n; ++k) {
      bit_elements[next_element++] = BitElement{k, n};
      if (next_element >= 64) return;
    }
  }
}

/*
 * Precompute a table. This value holds 16 nibbles, each representing how many edges
 * emanate from up to 16 nodes
 */
struct EdgeCount {
  uint64_t val;
  uint8_t get(int i) const { return (val >> (i * 4)) & 0xf; }
  void set(int i, uint8_t item) {
    val &= ~(uint64_t(0xf) << (i * 4));
    val |= uint64_t(item) << (i * 4);
  }
};

/*
 * Make edge count tables. 4 tables which handle 9 bits each for the 36-bit case.
 */
static EdgeCount edge_counts[4][512] = {0};
static void MakeEdgeCountTables() {
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 512; ++j) {
      EdgeCount &ec = edge_counts[i][j];
      unsigned k = 0;
      for (unsigned mask = 1; mask < 512; mask <<= 1) {
        if (j & mask) {
          BitElement &be = bit_elements[i * 9 + k];
          ec.set(be.from_node, ec.get(be.from_node) + 1);
          ec.set(be.to_node, ec.get(be.to_node) + 1);
        }
        ++k;
      }
    }
  }
}

static uint64_t nibble_table[256];

/*
 * This table will allow us to lookup byte-wise and accumulate a nibble-count nibble-wise,
 * i.e. lowest nibble will have number of zero-nibbles in the original, highest nibble
 * will have number of 0xf nibbles in the original.
 */
static void MakeNibbleCountTable() {
  unsigned next = 0;
  for (uint64_t i = 1; i; i <<= 4)
    for (uint64_t j = 1; j; j <<= 4) {
      nibble_table[next++] = i + j;
    }
}

static uint64_t SortNibbles(uint64_t word) {
  uint64_t nibble_counts = 0;

  for (int i = 0; i < 8; i++) nibble_counts += nibble_table[(word >> 8 * i) & 0xff];
  /* The table has allowed us to translate nibbles to counts. Now we have to reconstruct
   * the nibbles according to the counts */
  uint64_t output = 0;
  uint64_t nibbles = 0x1111111111111111;
  for (int i = 0; i < 16; i++) {
    unsigned current_count = (nibble_counts >> (4 * i)) & 0xf;
    nibbles >>= (current_count * 4);
    output += nibbles;
  }
  return output;
}

std::unordered_map<uint64_t, uint64_t> itemcounts;

std::unordered_map<uint64_t, uint64_t> trcounts;

static uint64_t CountEdges(uint64_t x) {
  uint64_t sav_x = x;
  uint64_t val = edge_counts[0][x % 512].val;
  x >>= 9;
  val += edge_counts[1][x % 512].val;
  x >>= 9;
  val += edge_counts[2][x % 512].val;
  x >>= 9;
  val += edge_counts[3][x % 512].val;
  EdgeCount ec{val};
  val = SortNibbles(val);
  ++itemcounts[val];
  x = ~sav_x;
  uint64_t tv = 0;
  if ((x & 0x7) == 0) tv +=  0x111;
  if ((x & 0x19) == 0) tv +=  0x1011;
  if ((x & 0x2a) == 0) tv +=  0x1101;
  if ((x & 0x34) == 0) tv +=  0x1110;
  if ((x & 0xc1) == 0) tv +=  0x10011;
  if ((x & 0x142) == 0) tv +=  0x10101;
  if ((x & 0x184) == 0) tv +=  0x10110;
  if ((x & 0x248) == 0) tv +=  0x11001;
  if ((x & 0x290) == 0) tv +=  0x11010;
  if ((x & 0x320) == 0) tv +=  0x11100;
  if ((x & 0xc01) == 0) tv +=  0x100011;
  if ((x & 0x1402) == 0) tv +=  0x100101;
  if ((x & 0x1804) == 0) tv +=  0x100110;
  if ((x & 0x2408) == 0) tv +=  0x101001;
  if ((x & 0x2810) == 0) tv +=  0x101010;
  if ((x & 0x3020) == 0) tv +=  0x101100;
  if ((x & 0x4440) == 0) tv +=  0x110001;
  if ((x & 0x4880) == 0) tv +=  0x110010;
  if ((x & 0x5100) == 0) tv +=  0x110100;
  if ((x & 0x6200) == 0) tv +=  0x111000;
  if ((x & 0x18001) == 0) tv +=  0x1000011;
  if ((x & 0x28002) == 0) tv +=  0x1000101;
  if ((x & 0x30004) == 0) tv +=  0x1000110;
  if ((x & 0x48008) == 0) tv +=  0x1001001;
  if ((x & 0x50010) == 0) tv +=  0x1001010;
  if ((x & 0x60020) == 0) tv +=  0x1001100;
  if ((x & 0x88040) == 0) tv +=  0x1010001;
  if ((x & 0x90080) == 0) tv +=  0x1010010;
  if ((x & 0xa0100) == 0) tv +=  0x1010100;
  if ((x & 0xc0200) == 0) tv +=  0x1011000;
  if ((x & 0x108400) == 0) tv +=  0x1100001;
  if ((x & 0x110800) == 0) tv +=  0x1100010;
  if ((x & 0x121000) == 0) tv +=  0x1100100;
  if ((x & 0x142000) == 0) tv +=  0x1101000;
  if ((x & 0x184000) == 0) tv +=  0x1110000;
  if ((x & 0x600001) == 0) tv +=  0x10000011;
  if ((x & 0xa00002) == 0) tv +=  0x10000101;
  if ((x & 0xc00004) == 0) tv +=  0x10000110;
  if ((x & 0x1200008) == 0) tv +=  0x10001001;
  if ((x & 0x1400010) == 0) tv +=  0x10001010;
  if ((x & 0x1800020) == 0) tv +=  0x10001100;
  if ((x & 0x2200040) == 0) tv +=  0x10010001;
  if ((x & 0x2400080) == 0) tv +=  0x10010010;
  if ((x & 0x2800100) == 0) tv +=  0x10010100;
  if ((x & 0x3000200) == 0) tv +=  0x10011000;
  if ((x & 0x4200400) == 0) tv +=  0x10100001;
  if ((x & 0x4400800) == 0) tv +=  0x10100010;
  if ((x & 0x4801000) == 0) tv +=  0x10100100;
  if ((x & 0x5002000) == 0) tv +=  0x10101000;
  if ((x & 0x6004000) == 0) tv +=  0x10110000;
  if ((x & 0x8208000) == 0) tv +=  0x11000001;
  if ((x & 0x8410000) == 0) tv +=  0x11000010;
  if ((x & 0x8820000) == 0) tv +=  0x11000100;
  if ((x & 0x9040000) == 0) tv +=  0x11001000;
  if ((x & 0xa080000) == 0) tv +=  0x11010000;
  if ((x & 0xc100000) == 0) tv +=  0x11100000;
  if ((x & 0x30000001) == 0) tv +=  0x100000011;
  if ((x & 0x50000002) == 0) tv +=  0x100000101;
  if ((x & 0x60000004) == 0) tv +=  0x100000110;
  if ((x & 0x90000008) == 0) tv +=  0x100001001;
  if ((x & 0xa0000010) == 0) tv +=  0x100001010;
  if ((x & 0xc0000020) == 0) tv +=  0x100001100;
  if ((x & 0x110000040) == 0) tv +=  0x100010001;
  if ((x & 0x120000080) == 0) tv +=  0x100010010;
  if ((x & 0x140000100) == 0) tv +=  0x100010100;
  if ((x & 0x180000200) == 0) tv +=  0x100011000;
  if ((x & 0x210000400) == 0) tv +=  0x100100001;
  if ((x & 0x220000800) == 0) tv +=  0x100100010;
  if ((x & 0x240001000) == 0) tv +=  0x100100100;
  if ((x & 0x280002000) == 0) tv +=  0x100101000;
  if ((x & 0x300004000) == 0) tv +=  0x100110000;
  if ((x & 0x410008000) == 0) tv +=  0x101000001;
  if ((x & 0x420010000) == 0) tv +=  0x101000010;
  if ((x & 0x440020000) == 0) tv +=  0x101000100;
  if ((x & 0x480040000) == 0) tv +=  0x101001000;
  if ((x & 0x500080000) == 0) tv +=  0x101010000;
  if ((x & 0x600100000) == 0) tv +=  0x101100000;
  if ((x & 0x810200000) == 0) tv +=  0x110000001;
  if ((x & 0x820400000) == 0) tv +=  0x110000010;
  if ((x & 0x840800000) == 0) tv +=  0x110000100;
  if ((x & 0x881000000) == 0) tv +=  0x110001000;
  if ((x & 0x902000000) == 0) tv +=  0x110010000;
  if ((x & 0xa04000000) == 0) tv +=  0x110100000;
  if ((x & 0xc08000000) == 0) tv +=  0x111000000;
  ++trcounts[tv];
  return val;
}

int main(int argc, char *argv[]) {
  int k = -1, n = -1;
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) k = atoi(argv[2]);
  if (k <= 0) {
    fprintf(stderr, "k (%d) must be greater than 0\n", k);
    return EXIT_FAILURE;
  }
  if (n <= k) {
    fprintf(stderr, "n (%d) must be greater than k(%d)\n", n, k);
    return EXIT_FAILURE;
  }
  if (n > 64) {
    fprintf(stderr, "n (%d) must be less than or equal to 64\n", n);
    return EXIT_FAILURE;
  }
  MakeBitElements();
  MakeEdgeCountTables();
  MakeNibbleCountTable();

  uint64_t bitperm = (uint64_t(1) << k) - 1;
  uint64_t endval = bitperm << (n - k);

  uint64_t count = 0;
  time_t start = time(nullptr);
  for (;;) {
    ++count;

    EdgeCount ec{CountEdges(bitperm)};
    if (bitperm == endval) break;

    /* These next two lines of code efficiently compute the lexicographically next bit
       permutation. This is from the famous "bit-twiddling hacks" by Sean Eron Anderson at
       Stanford University, which are in the public domain. Credit is given to acknowledge
       that this method is not original work.
    */
    uint64_t t = bitperm | (bitperm - 1);
    bitperm = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzl(bitperm) + 1));
  }
  time_t finish = time(nullptr);
  printf("count: %ld in %d seconds\n", count, (int)(finish - start));
  printf("itemcounts: %lu\n", itemcounts.size());
  std::unordered_map<unsigned, unsigned> validate;
  for (const auto &val : itemcounts) {
    EdgeCount ec{val.first};
    unsigned sum = 0;
    for (int i = 0; i < 9; ++i) sum += ec.get(i);
    ++validate[sum];
  }
  for (const auto &val : validate) {
    printf("%d: %d\n", val.first, val.second);
  }
  for (const auto &val : itemcounts) {
    EdgeCount ec{val.first};
    unsigned sum = 0;
    for (int i = 0; i < 9; ++i) {
      printf("%d ", ec.get(i));
      sum += ec.get(i);
    }
    printf(" %lX (%lu)\n", val.first, val.second);
    ++validate[sum];
  }
  printf("Triangles %ld\n", trcounts.size());
#if 0
  for (const auto &val : trcounts) {
    EdgeCount ec{val.first};
    unsigned sum = 0;
    for (int i = 0; i < 9; ++i) {
      printf("%d ", ec.get(i));
      sum += ec.get(i);
    }
    printf(" %lX (%lu)\n", val.first, val.second);
  }
#endif
  return 0;
}
