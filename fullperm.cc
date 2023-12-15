#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

void ReserveMemory() {
  size_t num_bytes = size_t(1) * 1024 * 1024 * 1024;
  void *res = mmap(nullptr, num_bytes, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (res == MAP_FAILED) {
    perror("mmap");
    exit(1);
  }
  constexpr size_t kPageSize = 4096;
  uint8_t *ptr = static_cast<uint8_t *>(res);
  struct timespec ts0, ts1;
  clock_gettime(CLOCK_MONOTONIC, &ts0);
  for (size_t i = 0; i < num_bytes; i += kPageSize) {
    ptr[i] = '\0';
  }
  clock_gettime(CLOCK_MONOTONIC, &ts1);
  uint64_t ns0 = ts0.tv_sec * 1000000000 + ts0.tv_nsec;
  uint64_t ns1 = ts1.tv_sec * 1000000000 + ts1.tv_nsec;
  uint64_t diff = ns1 - ns0;
  ldiv_t dt = ldiv(diff, 1000000000);
  printf("Populate %lu.%09lu\n", dt.quot, dt.rem);

  uint8_t data[4096];
  for (int i = 0; i < sizeof(data); ++i) data[i] = (i * 93) % 71;
  clock_gettime(CLOCK_MONOTONIC, &ts0);
  for (size_t i = 0; i < num_bytes; i += kPageSize) {
    memcpy(data, ptr + i, kPageSize);
  }
  clock_gettime(CLOCK_MONOTONIC, &ts1);
  ns0 = ts0.tv_sec * 1000000000 + ts0.tv_nsec;
  ns1 = ts1.tv_sec * 1000000000 + ts1.tv_nsec;
  diff = ns1 - ns0;
  dt = ldiv(diff, 1000000000);
  printf("Write %lu.%09lu\n", dt.quot, dt.rem);
  dt = ldiv(num_bytes * 1000, diff);
  printf("Write %lu Mbytes/sec\n", dt.quot);

  clock_gettime(CLOCK_MONOTONIC, &ts0);
  for (size_t i = 0; i < num_bytes; i += kPageSize) {
    memcpy(ptr + i, data, kPageSize);
  }
  clock_gettime(CLOCK_MONOTONIC, &ts1);
  ns0 = ts0.tv_sec * 1000000000 + ts0.tv_nsec;
  ns1 = ts1.tv_sec * 1000000000 + ts1.tv_nsec;
  diff = ns1 - ns0;
  dt = ldiv(diff, 1000000000);
  printf("Read %lu.%09lu\n", dt.quot, dt.rem);
  dt = ldiv(num_bytes * 1000, diff);
  printf("Read %lu Mbytes/sec\n", dt.quot);
}

/* the initial set 1...n, specify n as "nsetsize" */
unsigned nsetsize = 8;

static std::vector<unsigned char> swaps;

/**
 * This uses Steinhaus–Johnson–Trotter algorithm to generate each permutation exactly
 * once via a sequence of adjacent transpositions. We pre-compute this sequence and apply
 * it to generate orbits
 */
void GenerateSwaps(int n) {
  if (n == 2) {
    swaps.clear();
    swaps.push_back(0);
    return;
  }
  std::vector<unsigned char> old = std::move(swaps);
  swaps.clear();
  auto down = [&]() {
    for (int i = 0; i < n - 1; ++i) swaps.push_back(n - 2 - i);
  };
  auto up = [&]() {
    for (int i = 0; i < n - 1; ++i) swaps.push_back(i);
  };
  int next = 0;
  int ofs = 1;
  while (next < old.size()) {
    ofs ? down() : up();
    swaps.push_back(old[next++] + ofs);
    ofs = !ofs;
  }
  ofs ? down() : up();
}

/*
 * The segment table is organized so that swapping consecutive characters such as 'B' and
 * 'C' causes a number of bits with constant width to be swapped. Swapping 'A' and 'B'
 *  causes segments one-entry apart to be exchanged, swapping 'B' and 'C' causes segments
 *  two-entries apart to be exchanged, etc. This enables efficient graph node swapping.
 */
char segment_table[][2] = {
    {'A', 'B'}, {'B', 'C'}, {'A', 'C'}, {0},        {'B', 'D'}, {'A', 'D'}, {'C', 'D'},
    {'D', 'E'}, {'B', 'E'}, {'A', 'E'}, {'C', 'E'}, {0},        {'D', 'F'}, {'B', 'F'},
    {'A', 'F'}, {'C', 'F'}, {'E', 'F'}, {'F', 'G'}, {'D', 'G'}, {'B', 'G'}, {'A', 'G'},
    {'C', 'G'}, {'E', 'G'}, {0},        {'F', 'H'}, {'D', 'H'}, {'B', 'H'}, {'A', 'H'},
    {'C', 'H'}, {'E', 'H'}, {'G', 'H'}, {'H', 'I'}, {'F', 'I'}, {'D', 'I'}, {'B', 'I'},
    {'A', 'I'}, {'C', 'I'}, {'E', 'I'}, {'G', 'I'},
};

void PrintSegments(uint64_t bits) {
  while (bits) {
    unsigned bit = __builtin_ctzl(bits);
    printf("%s ", std::string(segment_table[bit], 2).c_str());
    bits ^= (uint64_t(1) << bit);
  }
};

static uint64_t masks[8] = {0};

/*
 * Because of the property of the table described above, we can apply a swap of nodes
 * on the graph with a sequence of bitwise operations. This uses the common trick of
 * XOR swapping (i.e.  a ^= b; b ^= a; a ^= b; swaps a and b) combined with the bitmask
 * that isolates the bits that need to be xored with bits (n + 1) bits away.
 */
uint64_t ApplySwap(uint64_t val, unsigned n) {
  uint64_t bits = (val & masks[n]) << (n + 1);
  val ^= bits;
  bits = (val >> (n + 1)) & masks[n];
  val ^= bits;
  bits = (val & masks[n]) << (n + 1);
  val ^= bits;
  return val;
}

/**
 * Validate that the segment table has the property that swapping nodes causes segments
 * a fixed distance from each other in the segment_table to be swapped
 */
static void AnalyzeSegmentTable() {
  std::map<std::string, int> index;
  int idx = 0;
  for (const auto &t : segment_table) {
    ++idx;
    if (t[0] == '\0') continue;
    if (t[0] - 'A' >= nsetsize) continue;
    if (t[1] - 'A' >= nsetsize) continue;
    index.emplace(std::string(t, 2), idx);
    std::string s(t, 2);
    std::swap(s[0], s[1]);
    index.emplace(s, idx);
  }
  for (const auto &i : index) printf("%s: %d\n", i.first.c_str(), i.second);
  for (int i = 0; i < nsetsize - 1; ++i) {
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
    masks[i] = mask;
  }

  /* Do a little test and printout */
  uint64_t start = 0;
  start |= (uint64_t(1) << (index["AB"] - 1));
  start |= (uint64_t(1) << (index["AC"] - 1));
  PrintSegments(start);
  uint64_t end = ApplySwap(start, 2);
  PrintSegments(end);
}

class Orbits {
 public:
  Orbits() {
    int n = 0;
    for (const auto &s : segment_table) {
      if (s[0] != 0 && s[0] - 'A' < nsetsize && s[1] - 'A' < nsetsize) {
        single_items_.push_back(uint64_t(1) << n);
        printf("single_item %016lX\n", single_items_.back());
      }
      ++n;
    }
  }
  /* basically what we're going to do here is to take the previous level of orbits, add in
     one of the single_items_ one by one, and permute it each time through all
     permutations to generate all of the orbit members , and then we remember the smallest
     representative and the size of the orbit for this member of the next level.
     */
  uint64_t NextImageOrbitLevel(void) {
    ++imagesize_;
    uint64_t num_analyzed = 0;
    std::set<uint64_t> oldorbits = std::move(prevorbits_);
    prevorbits_.clear();
    for (uint64_t old : oldorbits) {
      for (const auto &mask : single_items_) {
        if ((mask & old) != 0) continue;
        std::set<uint64_t> sorter;
        uint64_t next = mask | old;
        sorter.insert(next);
        ++num_analyzed;
        for (const auto &s : swaps) {
          next = ApplySwap(next, s);
          sorter.insert(next);
          ++num_analyzed;
        }
        if (prevorbits_.emplace(*sorter.begin()).second) {
          //          printf("found %016lX %lu\n",  *sorter.begin(), sorter.size());
          orbit_sizes_[*sorter.begin()] = sorter.size();
        }
      }
    }
    return num_analyzed;
  }

  void FindImageOrbits(void) {
    /* put an empty placeholder for index zero which is unused in this scheme */
    prevorbits_.emplace(0);

    /* compute orbits for each possible image size. The final divide by two is because
      with no non-repeated elements in image it can only use at most half of the full
      image */
    time_t start = time(nullptr);
    time_t level_start = start;
    for (unsigned i = 0; i < nsetsize * (nsetsize - 1) / 2 / 2; ++i) {
      printf("Level %d\n", i + 1);
      uint64_t num_analyzed = NextImageOrbitLevel();
      time_t now = time(nullptr);
      printf("%u seconds for %lu analyzed\n", (unsigned)(now - level_start),
             num_analyzed);
      level_start = now;
    }
    for (const auto &item : orbit_sizes_) {
      printf("%016lX; %u ", item.first, item.second);
      PrintSegments(item.first);
      printf("\n");
    }
    printf("%u seconds\n", (unsigned)(level_start - start));
  };

 private:
  std::set<uint64_t> prevorbits_;
  std::map<uint64_t, unsigned> orbit_sizes_;
  std::vector<uint64_t> single_items_;
  int imagesize_ = 0;
};

int main(int argc, char *argv[]) {
  ReserveMemory();
  nsetsize = argc > 1 ? atoi(argv[1]) : 0;
  if (nsetsize < 4) nsetsize = 4;
  for (int i = 2; i <= nsetsize; ++i) GenerateSwaps(i);
  AnalyzeSegmentTable();
  class Orbits o;
  o.FindImageOrbits();
  return 0;
}
