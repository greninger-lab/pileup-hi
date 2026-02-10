/*  sam.c -- SAM and BAM file I/O and manipulation.

    Copyright (C) 2008-2010, 2012-2025 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "htslib/khash.h"
#include "htslib/sam.h"

//---------------------------------
//---  Tweak overlapping reads
//---------------------------------

/**
 *  cigar_iref2iseq_set()  - find the first CMATCH setting the ref and the read
 * index cigar_iref2iseq_next() - get the next CMATCH base
 *  @cigar:       pointer to current cigar block (rw)
 *  @cigar_max:   pointer just beyond the last cigar block
 *  @icig:        position within the current cigar block (rw)
 *  @iseq:        position in the sequence (rw)
 *  @iref:        position with respect to the beginning of the read (iref_pos -
 * b->core.pos) (rw)
 *
 *  Returns BAM_CMATCH, -1 when there is no more cigar to process or the
 * requested position is not covered, or -2 on error.
 */
inline int cigar_iref2iseq_set(const uint32_t **cigar,
                               const uint32_t *cigar_max, hts_pos_t *icig,
                               hts_pos_t *iseq, hts_pos_t *iref) {
  hts_pos_t pos = *iref;
  if (pos < 0)
    return -1;
  *icig = 0;
  *iseq = 0;
  *iref = 0;
  while (*cigar < cigar_max) {
    int cig = (**cigar) & BAM_CIGAR_MASK;
    int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

    if (cig == BAM_CSOFT_CLIP) {
      (*cigar)++;
      *iseq += ncig;
      *icig = 0;
      continue;
    }
    if (cig == BAM_CHARD_CLIP || cig == BAM_CPAD) {
      (*cigar)++;
      *icig = 0;
      continue;
    }
    if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
      pos -= ncig;
      if (pos < 0) {
        *icig = ncig + pos;
        *iseq += *icig;
        *iref += *icig;
        return BAM_CMATCH;
      }
      (*cigar)++;
      *iseq += ncig;
      *icig = 0;
      *iref += ncig;
      continue;
    }
    if (cig == BAM_CINS) {
      (*cigar)++;
      *iseq += ncig;
      *icig = 0;
      continue;
    }
    if (cig == BAM_CDEL || cig == BAM_CREF_SKIP) {
      pos -= ncig;
      if (pos < 0)
        pos = 0;
      (*cigar)++;
      *icig = 0;
      *iref += ncig;
      continue;
    }
    hts_log_error("Unexpected cigar %d", cig);
    return -2;
  }
  *iseq = -1;
  return -1;
}
static inline int cigar_iref2iseq_next(const uint32_t **cigar,
                                       const uint32_t *cigar_max,
                                       hts_pos_t *icig, hts_pos_t *iseq,
                                       hts_pos_t *iref) {
  while (*cigar < cigar_max) {
    int cig = (**cigar) & BAM_CIGAR_MASK;
    int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

    if (cig == BAM_CMATCH || cig == BAM_CEQUAL || cig == BAM_CDIFF) {
      if (*icig >= ncig - 1) {
        *icig = -1;
        (*cigar)++;
        continue;
      }
      (*iseq)++;
      (*icig)++;
      (*iref)++;
      return BAM_CMATCH;
    }
    if (cig == BAM_CDEL || cig == BAM_CREF_SKIP) {
      (*cigar)++;
      (*iref) += ncig;
      *icig = -1;
      continue;
    }
    if (cig == BAM_CINS) {
      (*cigar)++;
      *iseq += ncig;
      *icig = -1;
      continue;
    }
    if (cig == BAM_CSOFT_CLIP) {
      (*cigar)++;
      *iseq += ncig;
      *icig = -1;
      continue;
    }
    if (cig == BAM_CHARD_CLIP || cig == BAM_CPAD) {
      (*cigar)++;
      *icig = -1;
      continue;
    }
    hts_log_error("Unexpected cigar %d", cig);
    return -2;
  }
  *iseq = -1;
  *iref = -1;
  return -1;
}

// Given overlapping read 'a' (left) and 'b' (right) on the same
// template, adjust quality values to zero for either a or b.
// Note versions 1.12 and earlier always removed quality from 'b' for
// matching bases.  Now we select a or b semi-randomly based on name hash.
// Returns 0 on success,
//        -1 on failure
int tweak_overlap_quality(bam1_t *a, bam1_t *b) {
  const uint32_t *a_cigar = bam_get_cigar(a),
                 *a_cigar_max = a_cigar + a->core.n_cigar;
  const uint32_t *b_cigar = bam_get_cigar(b),
                 *b_cigar_max = b_cigar + b->core.n_cigar;
  hts_pos_t a_icig = 0, a_iseq = 0;
  hts_pos_t b_icig = 0, b_iseq = 0;
  uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
  uint8_t *a_seq = bam_get_seq(a), *b_seq = bam_get_seq(b);

  hts_pos_t iref = b->core.pos;
  hts_pos_t a_iref = iref - a->core.pos;
  hts_pos_t b_iref = iref - b->core.pos;

  int a_ret =
      cigar_iref2iseq_set(&a_cigar, a_cigar_max, &a_icig, &a_iseq, &a_iref);
  if (a_ret < 0)
    // no overlap or error
    return a_ret < -1 ? -1 : 0;

  int b_ret =
      cigar_iref2iseq_set(&b_cigar, b_cigar_max, &b_icig, &b_iseq, &b_iref);
  if (b_ret < 0)
    // no overlap or error
    return b_ret < -1 ? -1 : 0;

  // Determine which seq is the one getting modified qualities.
  uint8_t amul, bmul;
  if (__ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(a))) & 1) {
    amul = 1;
    bmul = 0;
  } else {
    amul = 0;
    bmul = 1;
  }

  // Loop over the overlapping region nulling qualities in either
  // seq a or b.
  int err = 0;
  while (1) {
    // Step to next matching reference position in a and b
    while (a_ret >= 0 && a_iref >= 0 && a_iref < iref - a->core.pos)
      a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max, &a_icig, &a_iseq,
                                   &a_iref);
    if (a_ret < 0) { // done
      err = a_ret < -1 ? -1 : 0;
      break;
    }

    while (b_ret >= 0 && b_iref >= 0 && b_iref < iref - b->core.pos)
      b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig, &b_iseq,
                                   &b_iref);
    if (b_ret < 0) { // done
      err = b_ret < -1 ? -1 : 0;
      break;
    }

    if (iref < a_iref + a->core.pos)
      iref = a_iref + a->core.pos;

    if (iref < b_iref + b->core.pos)
      iref = b_iref + b->core.pos;

    iref++;

    // If A or B has a deletion then we catch up the other to this point.
    // We also amend quality values using the same rules for mismatch.
    if (a_iref + a->core.pos != b_iref + b->core.pos) {
      if (a_iref + a->core.pos < b_iref + b->core.pos &&
          b_cigar > bam_get_cigar(b) && bam_cigar_op(b_cigar[-1]) == BAM_CDEL) {
        // Del in B means it's moved on further than A
        do {
          a_qual[a_iseq] = amul ? a_qual[a_iseq] * 0.8 : 0;
          a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max, &a_icig, &a_iseq,
                                       &a_iref);
          if (a_ret < 0)
            return -(a_ret < -1); // 0 or -1
        } while (a_iref + a->core.pos < b_iref + b->core.pos);
      } else if (a_cigar > bam_get_cigar(a) &&
                 bam_cigar_op(a_cigar[-1]) == BAM_CDEL) {
        // Del in A means it's moved on further than B
        do {
          b_qual[b_iseq] = bmul ? b_qual[b_iseq] * 0.8 : 0;
          b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig, &b_iseq,
                                       &b_iref);
          if (b_ret < 0)
            return -(b_ret < -1); // 0 or -1
        } while (b_iref + b->core.pos < a_iref + a->core.pos);
      } else {
        // Anything else, eg ref-skip, we don't support here
        continue;
      }
    }

    // fprintf(stderr, "a_cig=%ld,%ld b_cig=%ld,%ld iref=%ld "
    //         "a_iref=%ld b_iref=%ld a_iseq=%ld b_iseq=%ld\n",
    //         a_cigar-bam_get_cigar(a), a_icig,
    //         b_cigar-bam_get_cigar(b), b_icig,
    //         iref, a_iref+a->core.pos+1, b_iref+b->core.pos+1,
    //         a_iseq, b_iseq);

    if (a_iseq > a->core.l_qseq || b_iseq > b->core.l_qseq)
      // Fell off end of sequence, bad CIGAR?
      return -1;

    // We're finally at the same ref base in both a and b.
    // Check if the bases match (confident) or mismatch
    // (not so confident).
    if (bam_seqi(a_seq, a_iseq) == bam_seqi(b_seq, b_iseq)) {
      // We are very confident about this base.  Use sum of quals
      int qual = a_qual[a_iseq] + b_qual[b_iseq];
      a_qual[a_iseq] = amul * (qual > 200 ? 200 : qual);
      b_qual[b_iseq] = bmul * (qual > 200 ? 200 : qual);
      ;
    } else {
      // Not so confident about anymore given the mismatch.
      // Reduce qual for lowest quality base.
      if (a_qual[a_iseq] > b_qual[b_iseq]) {
        // A highest qual base; keep
        a_qual[a_iseq] = 0.8 * a_qual[a_iseq];
        b_qual[b_iseq] = 0;
      } else if (a_qual[a_iseq] < b_qual[b_iseq]) {
        // B highest qual base; keep
        b_qual[b_iseq] = 0.8 * b_qual[b_iseq];
        a_qual[a_iseq] = 0;
      } else {
        // Both equal, so pick randomly
        a_qual[a_iseq] = amul * 0.8 * a_qual[a_iseq];
        b_qual[b_iseq] = bmul * 0.8 * b_qual[b_iseq];
      }
    }
  }

  return err;
}
