#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ORDERED_CONTAINERS

#ifdef ORDERED_CONTAINERS
#include <map>
#include <set>

#define STD_MAP_TYPE std::map
#define STD_SET_TYPE std::set
#else
#include <unordered_map>
#include <unordered_set>
#define STD_MAP_TYPE std::unordered_map
#define STD_SET_TYPE std::unordered_set
#endif

template <class T>
void EnumerateBitPermutations(unsigned n, unsigned k, const T &t) {
  uint64_t bitperm = (uint64_t(1) << k) - 1;
  uint64_t endval = bitperm << (n - k);

  for (;;) {
    if (!t(bitperm)) return;
    if (bitperm == endval) break;

    /* These next two lines of code efficiently compute the lexicographically next bit
       permutation. This is from the famous "bit-twiddling hacks" by Sean Eron Anderson at
       Stanford University, which are in the public domain. Credit is given to acknowledge
       that this method is not original work.
    */
    uint64_t t = bitperm | (bitperm - 1);
    bitperm = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzl(bitperm) + 1));
  }
}

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

STD_MAP_TYPE<uint64_t, uint64_t> itemcounts;

static uint64_t CountEdges(uint64_t x) {
  uint64_t val = edge_counts[0][x % 512].val;
  x >>= 9;
  val += edge_counts[1][x % 512].val;
  x >>= 9;
  val += edge_counts[2][x % 512].val;
  x >>= 9;
  val += edge_counts[3][x % 512].val;
  uint64_t sorted_val = SortNibbles(val);
  ++itemcounts[sorted_val];
  return val;
}

struct TriangleCount {
  uint64_t val[3];
  /* 84 triangles times 2-bits each */
  void IncTriangle(unsigned tr) {
    div_t res = div(tr, 32);
    val[res.quot] += (uint64_t(1) << (2 * res.rem));
  }
  TriangleCount &operator+=(const TriangleCount &rhs) {
    val[0] += rhs.val[0];
    val[1] += rhs.val[1];
    val[2] += rhs.val[2];
    return *this;
  }
  TriangleCount &operator&=(const TriangleCount &rhs) {
    val[0] &= rhs.val[0];
    val[1] &= rhs.val[1];
    val[2] &= rhs.val[2];
    return *this;
  }

  void MaskTriangles() {
    for (uint64_t &v : val) {
      v &= ((v >> 1) & 0x5555555555555555);
    }
  }
  unsigned CountTriangles() const {
    return __builtin_popcountl(val[0]) + __builtin_popcountl(val[1]) +
           __builtin_popcountl(val[2]);
  }
};

static TriangleCount triangle_node_mask[9] = {0};
static TriangleCount triangle_counts[4][512] = {0};

static TriangleCount CountTriangles(uint64_t x) {
  TriangleCount res = {0};
  res += triangle_counts[0][x % 512];
  x >>= 9;
  res += triangle_counts[1][x % 512];
  x >>= 9;
  res += triangle_counts[2][x % 512];
  x >>= 9;
  res += triangle_counts[3][x % 512];

  res.MaskTriangles();
  return res;
}

static void MakeTriangleTable() {
  STD_MAP_TYPE<uint64_t, STD_SET_TYPE<unsigned>> triangles;
  EnumerateBitPermutations(36, 3, [&](uint64_t bitperm) {
    if (SortNibbles(CountEdges(bitperm)) == 0x222) {
      triangles.emplace(bitperm, STD_SET_TYPE<unsigned>());
    }
    return true;
  });
  int n = 0;
  for (auto &tr : triangles) {
    EdgeCount ec{CountEdges(tr.first)};
    for (unsigned i = 0; i < 9; ++i) {
      if (ec.get(i) != 0) tr.second.insert(i);
    }
  }

  unsigned tr_index = 0;
  for (const auto tr : triangles) {
    assert(tr.second.size() == 3);
    for (unsigned from = 0; from < 9; ++from) {
      for (unsigned to = from + 1; to < 9; ++to) {
        if (tr.second.find(from) != tr.second.end() &&
            tr.second.find(to) != tr.second.end()) {
          triangle_node_mask[from].IncTriangle(tr_index);
          triangle_node_mask[to].IncTriangle(tr_index);
        }
      }
    }
    ++tr_index;
  }
  for (auto &tr : triangle_node_mask) {
    for (unsigned index = 0; index < 84; ++index) tr.IncTriangle(index);
    tr.MaskTriangles();
  }

  for (const auto &tr : triangle_node_mask) assert(tr.CountTriangles() == 28);

  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 512; ++j) {
      TriangleCount &tc = triangle_counts[i][j];
      unsigned k = 0;
      for (unsigned mask = 1; mask < 512; mask <<= 1) {
        BitElement &be = bit_elements[i * 9 + k];
        if (j & mask) {
          unsigned tr_index = 0;
          unsigned found = 0;
          for (const auto tr : triangles) {
            if (tr.second.find(be.from_node) != tr.second.end() &&
                tr.second.find(be.to_node) != tr.second.end()) {
              tc.IncTriangle(tr_index);
              ++found;
            }
            ++tr_index;
          }
          assert(found == 7);
        }
        ++k;
      }
    }
  }
}

void TestTriangleTable() {
  /* Test individual triangles */
  EnumerateBitPermutations(36, 3, [&](uint64_t bitperm) {
    uint64_t val = SortNibbles(CountEdges(bitperm));
    TriangleCount tc = CountTriangles(bitperm);
    unsigned num = tc.CountTriangles();
    if (val == 0x222) {
      assert(num == 1); /* it is a triangle */
      EdgeCount ec{CountEdges(bitperm)};
      for (unsigned node = 0; node < 9; ++node) {
        TriangleCount t = tc;
        t &= triangle_node_mask[node];
        unsigned count = t.CountTriangles();
        assert(ec.get(node) / 2 == count);
      }
    } else {
      assert(num == 0);
    }
    return true;
  });

  /* Test tetrahedra */
  EnumerateBitPermutations(36, 6, [&](uint64_t bitperm) {
    uint64_t val = SortNibbles(CountEdges(bitperm));
    if (val == 0x3333) {
      TriangleCount tc = CountTriangles(bitperm);
      unsigned num = tc.CountTriangles();
      assert(num == 4); /* it is a tetrahedron */
      EdgeCount ec{CountEdges(bitperm)};
      for (unsigned node = 0; node < 9; ++node) {
        TriangleCount t = tc;
        t &= triangle_node_mask[node];
        unsigned count = t.CountTriangles();
        assert(ec.get(node) == count);
      }
    }
    return true;
  });
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
  MakeTriangleTable();
  TestTriangleTable();

  itemcounts.clear();
  uint64_t count = 0;
  time_t start = time(nullptr);
  EnumerateBitPermutations(n, k, [&](uint64_t bitperm) {
    ++count;

    EdgeCount ec{CountEdges(bitperm)};
    return true;
  });
  time_t finish = time(nullptr);
  printf("count: %ld in %d seconds\n", count, (int)(finish - start));
  printf("itemcounts: %lu\n", itemcounts.size());
  STD_MAP_TYPE<unsigned, unsigned> validate;
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
    printf("(%lu)\n", val.second);
    ++validate[sum];
  }
  return 0;
}
