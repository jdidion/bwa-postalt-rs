extern crate bio_types;
use crate::utils;

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::rc::Rc;

use bio_types::genome::AbstractInterval;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{Read, Reader, Record};

use utils::{CigarStats, Strings};

pub trait Interval {
    fn start(&self) -> &u32;
    fn end(&self) -> &u32;
}

pub struct AltInterval {
    alt_start: u32,
    alt_end: u32,
    pub alt_len: u32,
    pub ref_name: Rc<String>,
    ref_end: u32,
    pub is_reverse: bool,
    pub pos: u32,
    cigar: Vec<Cigar>,
}

impl AltInterval {
    fn create<'a>(
        alt_read: &'a Record,
        alt_contig: Rc<String>,
        alt_cigar: Vec<Cigar>,
        cigar_stats: &'a CigarStats,
    ) -> AltInterval {
        let pos: u32 = alt_read.pos() as u32;
        let alt_len: u32 = cigar_stats.match_len + cigar_stats.ins_len;
        let ref_len: u32 = cigar_stats.match_len + cigar_stats.del_len + cigar_stats.skip_len;
        let clip_len: u32 = cigar_stats.hard_len + cigar_stats.soft_len;
        let start_idx = if alt_read.is_reverse() {
            alt_cigar.len() - 1
        } else {
            0
        };
        let alt_start = match alt_cigar[start_idx] {
            Cigar::SoftClip(l) => l,
            _ => 0,
        };
        AltInterval {
            alt_start,
            alt_end: alt_start + alt_len,
            alt_len: alt_len + clip_len,
            ref_name: alt_contig,
            ref_end: pos + ref_len,
            is_reverse: alt_read.is_reverse(),
            pos: pos - 1,
            cigar: alt_cigar,
        }
    }

    /// Returns the pos on REF given `pos` on ALT
    pub fn get_ref_pos(&self, pos: u32) -> Option<u32> {
        let mut x: u32 = 0;
        let mut y: u32 = 0;
        for i in 0..self.cigar.len() {
            match self.cigar[i] {
                Cigar::Match(l) => {
                    if y <= pos && pos < y + l {
                        return Some(x + (pos - y));
                    }
                    x += l;
                    y += l;
                }
                Cigar::Del(l) => x += l,
                Cigar::Ins(l) => {
                    if y <= pos && pos < y + l {
                        return Some(x);
                    }
                    y += l;
                }
                Cigar::SoftClip(l) | Cigar::HardClip(l) => {
                    if y <= pos && pos < y + l {
                        return None;
                    }
                    y += l;
                }
                _ => (),
            }
        }
        None
    }
}

impl Interval for AltInterval {
    fn start(&self) -> &u32 {
        &self.alt_start
    }

    fn end(&self) -> &u32 {
        &self.alt_end
    }
}

pub struct RefInterval {
    pub alt_name: Rc<String>,
    pub ref_pos: u32,
    pub ref_end: u32,
}

impl RefInterval {
    fn create<'a>(
        alt_read: &'a Record,
        alt_qname: Rc<String>,
        cigar_stats: &'a CigarStats,
    ) -> RefInterval {
        RefInterval {
            alt_name: alt_qname,
            ref_pos: alt_read.pos() as u32,
            ref_end: (alt_read.pos() as u32)
                + cigar_stats.match_len
                + cigar_stats.del_len
                + cigar_stats.skip_len,
        }
    }
}

impl Interval for RefInterval {
    fn start(&self) -> &u32 {
        &self.ref_pos
    }

    fn end(&self) -> &u32 {
        &self.ref_end
    }
}

/// number of bits to shift the genomic position to get the index bin
const BIN_SIZE_BITS: u8 = 13;

pub struct IntervalIndex<T: Interval> {
    ivls: Vec<T>,
    index: HashMap<u32, u32>,
    max_bin: u32,
}

impl<T: Interval> IntervalIndex<T> {
    fn create(mut ivls: Vec<T>) -> IntervalIndex<T> {
        // Sort the intervals by start position, using end position as tie-breaker. Use
        // sort-unstable since it's still possible to have ties and we don't care about the order
        // of intervals with the same start and end position.
        ivls.sort_unstable_by(|a: &T, b: &T| match a.start().cmp(b.start()) {
            Ordering::Equal => a.end().cmp(b.end()),
            ord => ord,
        });
        // build the index
        let mut index: HashMap<u32, u32> = HashMap::new();
        let mut max_bin: u32 = 0;
        ivls.iter().enumerate().for_each(|(i, ivl)| {
            let start_bin = ivl.start() >> BIN_SIZE_BITS;
            if !index.contains_key(&start_bin) {
                let end_bin = (ivl.end() - 1) >> BIN_SIZE_BITS;
                if start_bin != end_bin {
                    (start_bin..=end_bin).for_each(|j| {
                        index.insert(j, i as u32);
                    });
                } else {
                    index.insert(start_bin, i as u32);
                }
                if max_bin < end_bin {
                    max_bin = end_bin;
                }
            }
        });
        IntervalIndex {
            ivls,
            index,
            max_bin,
        }
    }

    pub fn find_overlaps(&self, start: &u32, end: &u32) -> Option<&[T]> {
        let start_bin: u32 = start >> BIN_SIZE_BITS;
        if start_bin > self.max_bin {
            None
        } else {
            let start_offset: usize = self.index.get(&start_bin).map(|i| *i as usize).unwrap_or(
                (0..=(((end - 1) >> BIN_SIZE_BITS) - 1))
                    .rev()
                    .find(|k| self.index.contains_key(k))
                    .map(|i| i as usize)
                    .unwrap_or(0),
            );
            let end_offset = (start_offset..self.ivls.len())
                .find(|i| {
                    let ivl: &T = &self.ivls[*i];
                    ivl.start() >= end || ivl.end() <= start
                })
                .unwrap_or(self.ivls.len());
            if start_offset == end_offset {
                None
            } else {
                Some(&self.ivls[start_offset..end_offset])
            }
        }
    }
}

pub struct AltToRef {
    pub alt_contigs: HashSet<Rc<String>>,
    pub hla_contig: Rc<String>,
    pub hla_counts: HashMap<String, u32>,
    pub alt_indexes: HashMap<Rc<String>, IntervalIndex<AltInterval>>,
    pub ref_indexes: HashMap<Rc<String>, IntervalIndex<RefInterval>>,
}

impl AltToRef {
    // process ALT-to-REF alignments
    pub fn parse(path: &PathBuf, contigs: &mut Strings) -> AltToRef {
        let mut alt = Reader::from_path(path).unwrap();
        // set of ALT contig names
        let mut alt_contigs: HashSet<Rc<String>> = HashSet::new();
        // the reference chromosome on which the HLA is located
        let mut hla_contig: Option<Rc<String>> = None;
        // counts for each HLA allele
        let mut hla_counts: HashMap<String, u32> = HashMap::new();
        // intervals that align to the reference for each ALT contig
        let mut alt_intervals: HashMap<Rc<String>, Vec<AltInterval>> = HashMap::new();
        // intervals that have an ALT alignment for each REF contig
        let mut ref_intervals: HashMap<Rc<String>, Vec<RefInterval>> = HashMap::new();
        for alt_read in alt
            .records()
            .map(|read| read.expect("Failure parsing alt file"))
        {
            let alt_name = contigs.intern(alt_read.qname());
            // add read names to set for later lookup
            alt_contigs.insert(alt_name.clone());
            // don't consider the record further if it is unmapped
            if alt_read.tid() == -1 || alt_read.is_unmapped() {
                continue;
            }
            let ref_name = contigs.intern_str(alt_read.contig());
            // if this is an HLA contig, check that the reference chromosomes are the same
            // and increment the count for the HLA allele
            if let Some(hla) = utils::get_hla_allele(alt_name.as_str()) {
                match &hla_contig {
                    Some(value) if *value != ref_name => {
                        panic!("HLA contigs map to more than one reference chromosome");
                    }
                    None => {
                        hla_contig = Some(ref_name.clone());
                    }
                    _ => (),
                }
                match hla_counts.get_mut(hla) {
                    Some(count) => *count += 1,
                    None => {
                        hla_counts.insert(hla.to_owned(), 1);
                    }
                }
            }
            let mut alt_cigar = alt_read.cigar().to_vec();
            let alt_cigar_stats = CigarStats::parse_cigar(&mut alt_cigar, true);
            let alt_interval =
                AltInterval::create(&alt_read, ref_name.clone(), alt_cigar, &alt_cigar_stats);
            match alt_intervals.entry(alt_name.clone()) {
                Entry::Vacant(e) => {
                    e.insert(vec![alt_interval]);
                }
                Entry::Occupied(mut e) => {
                    e.get_mut().push(alt_interval);
                }
            }
            let ref_interval = RefInterval::create(&alt_read, alt_name.clone(), &alt_cigar_stats);
            match ref_intervals.entry(ref_name.clone()) {
                Entry::Vacant(e) => {
                    e.insert(vec![ref_interval]);
                }
                Entry::Occupied(mut e) => {
                    e.get_mut().push(ref_interval);
                }
            }
        }
        // create indexes for intervals
        let alt_indexes: HashMap<Rc<String>, IntervalIndex<AltInterval>> = alt_intervals
            .into_iter()
            .map(|(contig, intervals)| (contig, IntervalIndex::create(intervals)))
            .collect();
        let ref_indexes: HashMap<Rc<String>, IntervalIndex<RefInterval>> = ref_intervals
            .into_iter()
            .map(|(contig, intervals)| (contig, IntervalIndex::create(intervals)))
            .collect();
        AltToRef {
            alt_contigs,
            hla_contig: hla_contig.unwrap(),
            hla_counts,
            alt_indexes,
            ref_indexes,
        }
    }
}
