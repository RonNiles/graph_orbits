#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>

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
    {'A', 'I'}, {'C', 'I'}, {'E', 'I'}, {'G', 'I'},
};

class MachineRepresentation {
 public:
  MachineRepresentation() {
    int index = 0;
    for (const char *segment : segment_table) {
      if (segment[0] != '\0') {
        segment_bits_.emplace(std::string_view(segment, 2), index);
        bit_segments_.emplace(index, std::string_view(segment, 2));
      }
      ++index;
    }
    BuildMasks(9);
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
  void BuildMasks(int nodes) {
    std::map<std::string, int> index;
    int idx = 0;
    for (const auto &t : segment_table) {
      ++idx;
      if (t[0] == '\0') continue;
      if (t[0] - 'A' >= nodes) continue;
      if (t[1] - 'A' >= nodes) continue;
      index.emplace(std::string(t, 2), idx);
      std::string s(t, 2);
      std::swap(s[0], s[1]);
      index.emplace(s, idx);
    }
    for (const auto &i : index) printf("%s: %d\n", i.first.c_str(), i.second);
    for (int i = 0; i < nodes - 1; ++i) {
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
  uint64_t masks_[9];
  std::map<std::string_view, int> segment_bits_;
  std::map<int, std::string_view> bit_segments_;
};

class EdgeCount {
 public:
  EdgeCount() {
    for (int i = 0; i < 9; ++i) {
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
    for (int i = 0; i < 9; ++i) totals[i] = __builtin_popcountl(val & masks_[i]);
  }

 private:
  uint64_t masks_[9];
};

class TriangleCount {
 public:
  TriangleCount() {
    int nodes = 9;
    for (char v1 = 'A'; v1 < 'A' + nodes; ++v1) {
      for (char v2 = v1 + 1; v2 < 'A' + nodes; ++v2) {
        for (char v3 = v2 + 1; v3 < 'A' + nodes; ++v3) {
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
    *reinterpret_cast<uint64_t *>(totals) = totals64[0];
    totals[8] = *reinterpret_cast<uint8_t *>(&totals64[1]);
  }

  void PrintCount(uint64_t val) {
    uint8_t totals[9];
    Count(val, totals);
    for (int i = 0; i < 9; ++i) {
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

 private:
  int order_;
  int count_;
};

class LevelManager {
public:
  void NextLevel() {
    std::set<uint64_t> old = std::move(elements_);
    elements_.clear();
    for (uint64_t graph : old) {
      uint64_t mask = uint64_t(1);
      for (const char *segment : segment_table) {
        if (segment[0] != 0 && (mask & graph) == 0) {
          TestNewItem(mask | graph);
        }
        mask <<= 1;
      }
    }
  }
private:
  void TestNewItem(uint64_t graph) {

  }
  MachineRepresentation mr_;
  EdgeCount ec;
  TriangleCount tc;
  std::set<uint64_t> elements_;
};

int main(int argc, char *argv[]) {
  MachineRepresentation mr;
  EdgeCount ec;
  TriangleCount tc;
  uint64_t rep = uint64_t(1) << mr.SegmentBit("AB");
  rep |= uint64_t(1) << mr.SegmentBit("AC");
  rep |= uint64_t(1) << mr.SegmentBit("BC");
  rep |= uint64_t(1) << mr.SegmentBit("GH");
  rep |= uint64_t(1) << mr.SegmentBit("GI");
  rep |= uint64_t(1) << mr.SegmentBit("HI");

  tc.PrintCount(rep);

  uint8_t edgecount[9];
  ec.Count(rep, edgecount);
  for (int i = 0; i < 9; ++i) {
    printf("%c: %d\n", 'A' + i, edgecount[i]);
  }

  printf("%s\n", mr.Segments(rep).c_str());
  for (int i = 0; i < 8; ++i) {
    rep = mr.ApplySwap(rep, i);
    printf("%s\n", mr.Segments(rep).c_str());
    tc.PrintCount(rep);
  }
  return 0;
}
