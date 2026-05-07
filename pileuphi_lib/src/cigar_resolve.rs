use crate::alignment::{PileupAlignment, CIGAR_STATE_UNINIT};
use log::error;
use rust_htslib::bam::record::Cigar;

#[inline(always)]
/// This is a port of htslib's cigar_resolver2 from sam.c. I didn't try to reinvent the wheel by
/// implementing it in highly-idiomatic Rust; the algorithm is delicate enough as is. Shouldn't be
/// changed unless you have a very good reason.
pub fn resolve_cigar(plp: &mut PileupAlignment, pos: i64) {
    let cs = &mut plp.cstate;
    let cig = &cs.cig;
    let ncig = cig.len();
    let mut k: usize = 0;
    let mut op: Cigar;

    //////////////////////////////////////////////////////////////
    // PHASE 1: SEEK TO CIGAR OPERATOR CONTAINING QUERY COORDINATE
    //////////////////////////////////////////////////////////////

    // never processed
    if cs.icig == CIGAR_STATE_UNINIT {
        plp.qpos = 0;

        if ncig == 1 {
            match cig[0] {
                Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                    cs.icig = 0;
                    cs.bam_pos = plp.rec.pos();
                    cs.iseq = 0;
                }
                _ => (),
            }
        } else {
            cs.icig = 0;
            cs.bam_pos = plp.rec.pos();
            cs.iseq = 0;

            for idx in 0..ncig {
                k = idx;
                match cig[k] {
                    Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_) => break,

                    Cigar::Ins(l) | Cigar::SoftClip(l) => {
                        cs.iseq += l as usize;
                    }

                    // pad and hardclip, they don't consume anything, so we move on.
                    _ => (),
                }
            }

            assert!(k < ncig);
            cs.icig = k;
        }
    } else {
        op = cig[cs.icig];
        // the position we want is not in this operator, but one downstream.
        if pos - cs.bam_pos >= i64::from(op.len()) {
            assert!(cs.icig < ncig);

            // Now we peek the next operator after this one. If it is reference consuming, we jump
            // to it.
            let op2 = cig[cs.icig + 1];
            match op2 {
                Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                    // if the old cigar we choose to move past is a read or reference consuming,
                    // update the indexes accordingly.
                    let op = cig[cs.icig];
                    match op {
                        Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                            cs.iseq += l as usize;
                        }

                        _ => (),
                    };

                    cs.bam_pos += op.len() as i64;
                    cs.icig += 1;
                }

                // the next operator is not reference consuming, so we move on until we find one.
                _ => {
                    // 1. update state by moving past current cigar operator.
                    let op = cig[cs.icig];
                    match op {
                        Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                            cs.iseq += l as usize;
                        }

                        _ => (),
                    }

                    cs.bam_pos += op.len() as i64;

                    // search all following cigar strings for a ref consumer
                    for idx in (cs.icig + 1)..ncig {
                        k = idx;

                        let next_op = cig[k];
                        match next_op {
                            Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                                break
                            } // found it!

                            // didn't find it, but need to up query pos...
                            Cigar::Ins(l) | Cigar::SoftClip(l) => cs.iseq += l as usize,
                            _ => (),
                        }
                    }

                    cs.icig = k;
                }
            }

            assert!(cs.icig < ncig);
        }
    }

    ////////////////////////////////////////////////////
    // PHASE 2: EXTRACTING PILEUP INFORMATION FROM CIGAR
    ////////////////////////////////////////////////////

    // At this stage, we have hit the cigar operator that contains our queried position. Now look
    // at the alignment and extract needed info.

    op = cig[cs.icig];
    // TODO: RESET OPERATORS HERE
    plp.indel = 0;
    plp.refskip = false;
    plp.del = false;

    // our position is right at the edge of an operation, so peek the next one.
    if (cs.bam_pos + op.len() as i64) as usize - 1 == pos as usize && cs.icig + 1 < ncig {
        let op2 = cig[cs.icig + 1];

        match op2 {
            Cigar::Del(l) if !matches!(op, Cigar::Del(_)) => {
                plp.indel = -(l as i32);
                for idx in (cs.icig + 2)..ncig {
                    k = idx;

                    // check that this is correct
                    if matches!(cs.cig[k], Cigar::Del(_)) {
                        plp.indel -= cs.cig[k].len() as i32;
                    } else {
                        break;
                    }
                }
            }

            Cigar::Ins(l) => {
                plp.indel = l as i32;
                for idx in (cs.icig + 2)..ncig {
                    k = idx;

                    match cig[k] {
                        Cigar::Ins(l) => {
                            plp.indel += l as i32;
                        }

                        Cigar::Pad(_) => (),
                        _ => break,
                    };
                }
            }

            Cigar::Pad(_) if cs.icig + 2 < ncig => {
                let mut next_op: Cigar;
                let mut total_insertion_length = 0;
                for idx in (cs.icig + 2)..ncig {
                    k = idx;
                    next_op = cig[k];
                    match next_op {
                        Cigar::Ins(l) => {
                            total_insertion_length += l;
                        }

                        Cigar::Pad(_) => (),
                        _ => break,
                    }
                }

                if total_insertion_length > 0 {
                    plp.indel = total_insertion_length as i32;
                }
            }

            // CHECK THAT THIS IS CORRECT
            _ => (),
        }
    }

    match op {
        Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
            plp.qpos = cs.iseq + (pos - cs.bam_pos) as usize;
        }

        Cigar::Del(_) | Cigar::RefSkip(_) => {
            plp.del = true;
            plp.refskip = matches!(op, Cigar::RefSkip(_));
            plp.qpos = cs.iseq;
        }

        _ => {
            error!("Bad cigar resolution");
        }
    }

    plp.cigar_index = plp.cstate.icig;
    plp.head = pos == plp.rec.pos();
    plp.tail = pos == plp.rec.pos() + plp.cstate.read_len_from_cigar - 1;
}
