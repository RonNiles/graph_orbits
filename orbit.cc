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

std::unordered_map<uint64_t, uint64_t> itemcounts;

static uint64_t CountEdges(uint64_t x) {
  uint64_t val = edge_counts[0][x % 512].val;
  x >>= 9;
  val += edge_counts[1][x % 512].val;
  x >>= 9;
  val += edge_counts[2][x % 512].val;
  x >>= 9;
  val += edge_counts[3][x % 512].val;
  ++itemcounts[val];
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

  uint64_t bitperm = (uint64_t(1) << k) - 1;
  uint64_t endval = bitperm << (n - k);

  uint64_t count = 0;
  time_t start = time(nullptr);
  for (;;) {
    ++count;
    if (bitperm == endval) break;

    EdgeCount ec{CountEdges(bitperm)};
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
  return 0;
}
