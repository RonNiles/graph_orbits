#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {
  int k = -1, n = -1;
  ;
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
  uint64_t bitperm = (uint64_t(1) << k) - 1;
  uint64_t endval = bitperm << (n - k);

  uint64_t count = 0;
  time_t start = time(nullptr);
  for (;;) {
    ++count;
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
  return 0;
}
