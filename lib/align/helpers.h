#ifndef HELPERS_H_INCLUDED
#define HELPERS_H_INCLUDED

#define SWAP(x, y) do { typeof(x) _v = x; x = y; y = _v; } while(0)

// ============================================================================
static inline void to_raw(const char *in_seq,
                          unsigned char *out_seq,
                          int seq_len,
                          const unsigned char *map) {
  int i;
  for (i = 0; i < seq_len; i++) out_seq[i] = (unsigned char)map[(unsigned)in_seq[i]];
}

#endif
