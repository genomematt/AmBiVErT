#include "align.h"

int align_frag_count(AlignFrag *f) {
  int i = 0;
  while (f) { i++; f=f->next; }
  return i;
}

void align_frag_free(AlignFrag *f) {
  while (f) {
    AlignFrag *n = f->next;
    free(f);
    f = n;
  }
}

void alignment_free(Alignment *alignment) {
  align_frag_free(alignment->align_frag);
  free(alignment);
}

Alignment *alignment_new(AlignFrag *align_frag, int score) {
  Alignment *alignment = malloc(sizeof(Alignment));
  if (alignment) {
    alignment->align_frag = align_frag;
    alignment->frag_count = align_frag_count(align_frag);
    alignment->score = score;
  }
  return alignment;
}

