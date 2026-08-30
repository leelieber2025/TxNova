//! 合同 K — canonical splice dinucleotides on the transcript strand.

use crate::error::Result;
use crate::fasta::FastaIndex;
use crate::gtf::Transcript;

pub fn is_canonical(donor: &[u8; 2], acceptor: &[u8; 2]) -> bool {
    matches!((donor, acceptor), (b"GT", b"AG") | (b"GC", b"AG") | (b"AT", b"AC"))
}

/// For an unstranded transcript (GTF strand `.`), decide whether the plus
/// or the reverse-complement reading of its introns is the real one.
///
/// The library carries no strand information, so `t.strand` cannot be
/// trusted; the intron donor/acceptor motif can. Score every intron in
/// both directions and keep whichever direction has more canonical
/// (GT-AG/GC-AG/AT-AC) hits, ties and no-evidence loci defaulting to plus.
/// Both the canonical-splice gate (below) and the spliced-sequence /
/// ORF scan (`orf::splice_seq`) call this so a locus gets one consistent
/// strand call rather than being silently read as plus by each in turn.
pub fn infer_unstranded_plus(t: &Transcript, fa: &FastaIndex) -> Result<bool> {
    let introns = crate::coverage::transcript_introns(t);
    let mut plus_can = 0usize;
    let mut minus_can = 0usize;
    for &(intron_g_start, intron_g_end) in &introns {
        if intron_g_end < intron_g_start + 3 {
            continue;
        }
        let fwd = (
            fa.dinuc_tx(&t.chrom, intron_g_start, true)?,
            fa.dinuc_tx(&t.chrom, intron_g_end - 1, true)?,
        );
        if is_canonical(&fwd.0, &fwd.1) {
            plus_can += 1;
        }
        let rev = (
            fa.dinuc_tx(&t.chrom, intron_g_end - 1, false)?,
            fa.dinuc_tx(&t.chrom, intron_g_start, false)?,
        );
        if is_canonical(&rev.0, &rev.1) {
            minus_can += 1;
        }
    }
    Ok(plus_can >= minus_can)
}

/// Per-intron donor/acceptor on the transcript strand + canonical fraction.
///
/// Intron set is `coverage::transcript_introns` so `donors` and
/// `junction_support` stay the same length and order (合同 K).
pub fn splice_features(
    t: &Transcript,
    fa: &FastaIndex,
) -> Result<(String, String, f64, usize, Option<String>)> {
    let n_exons = t.n_exons();
    if n_exons < 2 {
        // No intron to score. Gate skips canonical splice when n_introns=0.
        return Ok((String::new(), String::new(), 0.0, 0, None));
    }
    let introns = crate::coverage::transcript_introns(t);
    // Abutting / overlapping exons happen in stringtie --merge. Do not
    // fail the whole run; the locus is marked and dropped at the gate.
    let structure_error = if introns.len() != n_exons - 1 {
        Some(format!(
            "transcript {} has {} exons but {} introns (overlapping or abutting exons)",
            t.transcript_id,
            n_exons,
            introns.len()
        ))
    } else {
        None
    };
    // Unstranded loci (strand '.') must not be silently scored as plus:
    // pick the better-supported reading once for the whole transcript.
    let plus = if t.strand == '-' {
        false
    } else if t.strand == '+' {
        true
    } else {
        infer_unstranded_plus(t, fa)?
    };
    let mut donors = Vec::new();
    let mut acceptors = Vec::new();
    let mut n_can = 0usize;
    for &(intron_g_start, intron_g_end) in &introns {
        // Donor + acceptor are 2 bp each; shorter than 4 bp cannot be scored.
        if intron_g_end < intron_g_start + 3 {
            donors.push("NN".into());
            acceptors.push("NN".into());
            continue;
        }
        let (donor, acceptor) = if plus {
            (
                fa.dinuc_tx(&t.chrom, intron_g_start, true)?,
                fa.dinuc_tx(&t.chrom, intron_g_end - 1, true)?,
            )
        } else {
            (
                fa.dinuc_tx(&t.chrom, intron_g_end - 1, false)?,
                fa.dinuc_tx(&t.chrom, intron_g_start, false)?,
            )
        };
        if is_canonical(&donor, &acceptor) {
            n_can += 1;
        }
        donors.push(format!("{}{}", donor[0] as char, donor[1] as char));
        acceptors.push(format!("{}{}", acceptor[0] as char, acceptor[1] as char));
    }
    let frac = if donors.is_empty() {
        0.0
    } else {
        n_can as f64 / donors.len() as f64
    };
    Ok((
        donors.join(","),
        acceptors.join(","),
        frac,
        donors.len(),
        structure_error,
    ))
}

#[cfg(test)]
mod tests {
    use super::is_canonical;

    #[test]
    fn canonical_includes_atac() {
        assert!(is_canonical(b"GT", b"AG"));
        assert!(is_canonical(b"GC", b"AG"));
        assert!(is_canonical(b"AT", b"AC"));
        assert!(!is_canonical(b"AT", b"AG"));
    }
}
