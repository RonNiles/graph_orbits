#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>

#include "btree_map.h"
#include "btree_set.h"

/**
 * 1. machine representation
 * 2. permutations (swap adjacent)
 * 3. Filter 1: number of edges
 * 4. Filter 2: number of triangles
 * 5. Filtered code via bubble sort
 * 6. Orderly permutations of filtered code
 * 7. Maximal code
 */

/*
 * The segment table is organized so that swapping consecutive characters such as 'B' and
 * 'C' causes a number of bits with constant width to be swapped. Swapping 'A' and 'B'
 *  causes segments one-entry apart to be exchanged, swapping 'B' and 'C' causes segments
 *  two-entries apart to be exchanged, etc. This enables efficient graph node swapping.
 */
static char segment_table[][2] = {
    {'A', 'B'}, {'B', 'C'}, {'A', 'C'}, {0},        {'B', 'D'}, {'A', 'D'}, {'C', 'D'},
    {'D', 'E'}, {'B', 'E'}, {'A', 'E'}, {'C', 'E'}, {0},        {'D', 'F'}, {'B', 'F'},
    {'A', 'F'}, {'C', 'F'}, {'E', 'F'}, {'F', 'G'}, {'D', 'G'}, {'B', 'G'}, {'A', 'G'},
    {'C', 'G'}, {'E', 'G'}, {0},        {'F', 'H'}, {'D', 'H'}, {'B', 'H'}, {'A', 'H'},
    {'C', 'H'}, {'E', 'H'}, {'G', 'H'}, {'H', 'I'}, {'F', 'I'}, {'D', 'I'}, {'B', 'I'},
    {'A', 'I'}, {'C', 'I'}, {'E', 'I'}, {'G', 'I'}, {0},        {'H', 'J'}, {'F', 'J'},
    {'D', 'J'}, {'B', 'J'}, {'A', 'J'}, {'C', 'J'}, {'E', 'J'}, {'G', 'J'}, {'I', 'J'},
    {'J', 'K'}, {'H', 'K'}, {'F', 'K'}, {'D', 'K'}, {'B', 'K'}, {'A', 'K'}, {'C', 'K'},
    {'E', 'K'}, {'G', 'K'}, {'I', 'K'}};

template <unsigned N>
static constexpr unsigned EndIndex() {
  unsigned index = 0;
  constexpr char endchar = 'A' + N;
  for (const char *segment : segment_table) {
    if (segment[0] >= endchar || segment[1] >= endchar) break;
    ++index;
  }
  return index;
}

template <unsigned N>
class MachineRepresentation {
 public:
  MachineRepresentation() {
    int index = 0;
    constexpr char endchar = 'A' + N;
    for (const char *segment : segment_table) {
      if (segment[0] != '\0') {
        if (segment[0] >= endchar || segment[1] >= endchar) break;
        segment_bits_.emplace(std::string_view(segment, 2), index);
        bit_segments_.emplace(index, std::string_view(segment, 2));
      }
      ++index;
    }
    BuildMasks();
  }

  int SegmentBit(const char *seg) const {
    auto it = segment_bits_.find(seg);
    return it == segment_bits_.end() ? -1 : it->second;
  }
  std::string BitSegment(int bit) const { return std::string(segment_table[bit], 2); }

  std::string Segments(uint64_t val) {
    std::string s;
    const char *sep = "";
    while (val) {
      int bit = __builtin_ctzl(val);
      val &= ~(uint64_t(1) << bit);
      s += sep;
      s += BitSegment(bit);
      sep = " ";
    }
    return s;
  }
  /*
   * Because of the property of the table described above, we can apply a swap of nodes
   * on the graph with a sequence of bitwise operations. This uses the common trick of
   * XOR swapping (i.e.  a ^= b; b ^= a; a ^= b; swaps a and b) combined with the bitmask
   * that isolates the bits that need to be xored with bits (n + 1) bits away.
   */
  uint64_t ApplySwap(uint64_t val, unsigned n) {
    uint64_t bits = (val & masks_[n]) << (n + 1);
    val ^= bits;
    bits = (val >> (n + 1)) & masks_[n];
    val ^= bits;
    bits = (val & masks_[n]) << (n + 1);
    val ^= bits;
    return val;
  }

 private:
  void BuildMasks() {
    std::map<std::string, int> index;
    int idx = 0;
    for (const auto &t : segment_table) {
      ++idx;
      if (t[0] == '\0') continue;
      if (t[0] - 'A' >= N) continue;
      if (t[1] - 'A' >= N) continue;
      index.emplace(std::string(t, 2), idx);
      std::string s(t, 2);
      std::swap(s[0], s[1]);
      index.emplace(s, idx);
    }
    for (const auto &i : index) printf("%s: %d\n", i.first.c_str(), i.second);
    for (int i = 0; i < N - 1; ++i) {
      printf("Swap %c<->%c\n", 'A' + i, 'B' + i);
      uint64_t mask = 0;
      for (const auto &item : index) {
        std::string test(item.first);
        if (test[0] == 'A' + i)
          test[0] = 'B' + i;
        else if (test[0] == 'B' + i)
          test[0] = 'A' + i;
        if (test[1] == 'A' + i)
          test[1] = 'B' + i;
        else if (test[1] == 'B' + i)
          test[1] = 'A' + i;
        if (index[test] == 0) abort();
        if (item.second == index[test]) continue;
        int diff = item.second - index[test];
        if (diff < 0) {
          diff = -diff;
        } else {
          printf("%d -> %d\n", index[test], item.second);
          mask |= uint64_t(1) << (index[test] - 1);
        }
        if (diff != i + 1) abort();
      }
      printf("Mask for %d is %016lX\n", i + 1, mask);
      masks_[i] = mask;
    }
  }
  uint64_t masks_[N];
  std::map<std::string_view, int> segment_bits_;
  std::map<int, std::string_view> bit_segments_;
};

template <unsigned N>
class EdgeCount {
 public:
  EdgeCount() {
    for (int i = 0; i < N; ++i) {
      int bit = 0;
      uint64_t mask = 0;
      for (const char *segment : segment_table) {
        if (segment[0] == 'A' + i || segment[1] == 'A' + i) mask |= uint64_t(1) << bit;
        ++bit;
      }
      masks_[i] = mask;
    }
  }
  void Count(uint64_t val, uint8_t *totals) const {
    for (int i = 0; i < N; ++i) totals[i] = __builtin_popcountl(val & masks_[i]);
  }

 private:
  uint64_t masks_[N];
};

template <unsigned N>
class TriangleCount {
 public:
  TriangleCount() {
    for (char v1 = 'A'; v1 < 'A' + N; ++v1) {
      for (char v2 = v1 + 1; v2 < 'A' + N; ++v2) {
        for (char v3 = v2 + 1; v3 < 'A' + N; ++v3) {
          triangles_.emplace_back(std::move(Triangle{0}));
          uint8_t *counts = reinterpret_cast<uint8_t *>(triangles_.back().counts);
          counts[v1 - 'A'] = 1;
          counts[v2 - 'A'] = 1;
          counts[v3 - 'A'] = 1;
          uint64_t bit = 1;
          for (const char *segment : segment_table) {
            if ((segment[0] == v1 && segment[1] == v2) ||
                (segment[0] == v1 && segment[1] == v3) ||
                (segment[0] == v2 && segment[1] == v3)) {
              triangles_.back().mask |= bit;
            }
            bit <<= 1;
          }
        }
      }
    }
  }
  void Count(uint64_t val, uint8_t *totals) const {
    uint64_t totals64[2] = {0};
    for (const Triangle &t : triangles_) {
      if ((val & t.mask) == t.mask) {
        totals64[0] += t.counts[0];
        totals64[1] += t.counts[1];
      }
    }
    memcpy(totals, &totals64, N);
  }

  void PrintCount(uint64_t val) {
    uint8_t totals[N];
    Count(val, totals);
    for (int i = 0; i < N; ++i) {
      printf("%c: %02u ", 'A' + i, totals[i]);
    }
    printf("\n");
  }

 private:
  struct Triangle {
    uint64_t mask;
    uint64_t counts[2];
  };
  std::vector<Triangle> triangles_;
};

class AdjacentSwapPermutation {
 public:
  void Initialize(int order) {
    order_ = count_ = order;
    if (order > 2) (this + 1)->Initialize(order - 1);
  }
  int Next() {
    --count_;
    if (count_ == -order_) {
      count_ = order_;
      return order_ > 2 ? (this + 1)->Next() : -1;
    }
    if (count_ == 0) return order_ > 2 ? (this + 1)->Next() + 1 : -1;
    return count_ > 0 ? count_ - 1 : -count_ - 1;
  }
  int Order() const { return order_; }

 private:
  int order_;
  int count_;
};

static constexpr uint64_t Factorial(uint64_t x) {
  return x <= 1 ? 1 : x * Factorial(x - 1);
}

template <int N>
class LevelManager {
 public:
  LevelManager() { elements_[0] = 1; }
  void NextLevel() {
    ++level_;
    btree::btree_map<uint64_t, unsigned> old = std::move(elements_);
    elements_.clear();
    for (auto [graph, count] : old) {
      uint64_t mask = uint64_t(1);
      unsigned end_index = EndIndex<N>();
      for (unsigned i = 0; i < end_index; ++i) {
        const char *segment = segment_table[i];
        if (segment[0] != 0 && (mask & graph) == 0) {
          TestNewItem(mask | graph);
        }
        mask <<= 1;
      }
    }
  }

  void DumpElements() {
    uint64_t total = 0;
    for (auto [graph, count] : elements_) {
      total += count;
      printf("%016lX; %u %s\n", graph, count, mr_.Segments(graph).c_str());
    }
    uint64_t expected = 1;
    if (level_) {
      uint64_t num = N * (N - 1) / 2;
      uint64_t denom = 1;
      for (int i = 1; i <= level_; ++i) {
        expected *= num;
        expected /= denom;
        ++denom;
        --num;
      }
    }
    printf("Total: %lu %lu\n", total, expected);
    if (total != expected) abort();
  }

  size_t Size() const { return elements_.size(); }
  unsigned MaxLevel() const { return N * (N - 1) / 4; }

 private:
  void TestNewItem(uint64_t graph) {
    uint8_t edges[N];
    uint8_t triangles[N];
    uint8_t combined[N];
    ec_.Count(graph, edges);
    tc_.Count(graph, triangles);
    for (int i = 0; i < sizeof(edges); ++i) {
      int val = N * triangles[i] + edges[i];
      if (val > 255) abort();
      combined[i] = uint8_t(val);
    }

    // printf("Test %s\n", mr_.Segments(graph).c_str());
    /* bubble sort */
    int n = sizeof(combined);
    do {
      int new_n = 0;
      for (int i = 1; i < n; ++i) {
        if (combined[i - 1] < combined[i]) {
          std::swap(combined[i - 1], combined[i]);
          graph = mr_.ApplySwap(graph, i - 1);
          new_n = i;
        }
      }
      n = new_n;
    } while (n > 1);

    // for (uint8_t c : combined) printf("%u ", c);
    // printf("\n");

    AdjacentSwapPermutation asp[N];
    std::vector<int> bases;

    uint64_t min_code = graph;
    if (elements_.find(min_code) != elements_.end()) return;
    std::set<uint64_t> all_codes;
    for (int i = 0; i < N; ++i) {
      if (bases.empty() || combined[bases.back()] != combined[i]) bases.push_back(i);
    }
    if (combined[bases.back()] != 0) bases.push_back(sizeof(combined));
    for (int i = 0; i < bases.size() - 1; ++i)
      asp[bases[i]].Initialize(bases[i + 1] - bases[i]);
    bases.pop_back();
    int singletons = 0;
    for (int i = 0; i < bases.size();) {
      int order = asp[bases[i]].Order();
      if (order <= 1) {
        if (order == 1) ++singletons;
        bases.erase(bases.begin() + i);
      } else {
        ++i;
      }
    }

    if (bases.empty()) {
      uint64_t val = Factorial(N) / Factorial(N - singletons);
      elements_[graph] = val;
      return;
    }
    all_codes.insert(graph);
    for (;;) {
      int p = bases.size() - 1;
      int col = bases[p];
      int swap = asp[col].Next();
      while (swap < 0) {
        graph = mr_.ApplySwap(graph, col);
        if (p == 0) {
          uint64_t factor = Factorial(N);
          int remain = N;
          for (int base : bases) {
            int order = asp[base].Order();
            remain -= order;
            factor /= Factorial(order);
          }
          remain -= singletons;
          if (remain < 0) abort();
          factor /= Factorial(remain);
          elements_[min_code] = all_codes.size() * factor;

          //          printf("insert %lx\n", min_code);
          return;
        }
        --p;
        col = bases[p];
        swap = asp[col].Next();
      }
      graph = mr_.ApplySwap(graph, swap + col);
      all_codes.insert(graph);
      if (graph < min_code) {
        min_code = graph;
        if (elements_.find(min_code) != elements_.end()) return;
      }
    }
  }
  MachineRepresentation<N> mr_;
  EdgeCount<N> ec_;
  TriangleCount<N> tc_;
  unsigned level_ = 0;
  btree::btree_map<uint64_t, unsigned> elements_;
};

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s <num_nodes>\n", argv[0]);
    return 1;
  }
  const int num_nodes = atoi(argv[1]);
  if (num_nodes < 2 || num_nodes > 11) {
    fprintf(stderr, "Nodes (%d) must be between 2 and 11\n", num_nodes);
    return 1;
  }

  auto enumerate = [](auto *lm) {
    printf("Level 0 count %lu\n", lm->Size());
    lm->DumpElements();
    time_t prev = time(nullptr);
    for (int i = 1; i <= lm->MaxLevel(); ++i) {
      lm->NextLevel();
      time_t now = time(nullptr);
      printf("Level %d count %lu elapsed %u\n", i, lm->Size(), unsigned(now - prev));
      prev = now;
      lm->DumpElements();
    }
  };
  switch (num_nodes) {
    case 2: {
      LevelManager<2> lm;
      enumerate(&lm);
      break;
    }
    case 3: {
      LevelManager<3> lm;
      enumerate(&lm);
      break;
    }
    case 4: {
      LevelManager<4> lm;
      enumerate(&lm);
      break;
    }
    case 5: {
      LevelManager<5> lm;
      enumerate(&lm);
      break;
    }
    case 6: {
      LevelManager<6> lm;
      enumerate(&lm);
      break;
    }
    case 7: {
      LevelManager<7> lm;
      enumerate(&lm);
      break;
    }
    case 8: {
      LevelManager<8> lm;
      enumerate(&lm);
      break;
    }
    case 9: {
      LevelManager<9> lm;
      enumerate(&lm);
      break;
    }
    case 10: {
      LevelManager<10> lm;
      enumerate(&lm);
      break;
    }
    case 11: {
      LevelManager<11> lm;
      enumerate(&lm);
      break;
    }
  }
  return 0;
}
