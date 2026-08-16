//! 合同 K — canonical splice dinucleotides on the transcript strand.

use crate::error::Result;
use crate::fasta::FastaIndex;
use crate::gtf::Transcript;

pub fn is_canonical(donor: &[u8; 2], acceptor: &[u8; 2]) -> bool {
    matches!((donor, acceptor), (b"GT", b"AG") | (b"GC", b"AG"))
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
    let plus = t.strand != '-';
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
