use crate::cli::{Opt, StructOpt};

mod cli;
mod intervals;
mod utils;

use std::cell::Cell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;
use std::path::Path;
use std::rc::Rc;
use std::slice::IterMut;

use bio_types::genome::AbstractInterval;
use bio_types::sequence::SequenceRead;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{Header, Read, Reader, Record, Writer};

use cli::HtsFormat;
use intervals::{AltToRef, IntervalIndex, RefInterval};
use utils::{complement, get_hla_allele, revcomp, CigarStats, Strings};

struct Buffers {
    // buffer for all reads associated with a given query and segment number
    read_buf: Vec<Record>,
    // set of all HLA contigs with hits
    hla_hits: HashSet<Rc<String>>,
    ref_indexes: HashMap<Rc<String>, IntervalIndex<RefInterval>>,
    // prefix for HLA hit files
    hla_file_prefix: Option<String>,
    // files where HLA hits are written
    hla_files: Option<HashMap<Rc<String>, File>>,
}

impl Buffers {
    fn create(
        prefix: Option<String>,
        ref_indexes: HashMap<Rc<String>, IntervalIndex<RefInterval>>,
    ) -> Buffers {
        Buffers {
            read_buf: Vec::new(),
            hla_hits: HashSet::new(),
            ref_indexes,
            hla_files: if prefix.is_some() {
                Some(HashMap::new())
            } else {
                None
            },
            hla_file_prefix: prefix,
        }
    }

    fn push_read(&mut self, read: Record) {
        self.read_buf.push(read);
    }

    fn append_reads(&mut self, reads: &mut Vec<Record>) {
        self.read_buf.append(reads)
    }

    fn last_read(&self) -> Option<&Record> {
        self.read_buf.last()
    }

    fn read_iter_mut(&mut self) -> IterMut<'_, Record> {
        self.read_buf.iter_mut()
    }

    // adds HLA alleles for any intervals that overlap the given range
    fn collect_hla_hits(&mut self, contig: &Rc<String>, start: &u32, end: &u32) {
        if let Some(ref idx) = self.ref_indexes.get(contig) {
            if let Some(overlaps) = idx.find_overlaps(start, end) {
                for ivl in overlaps {
                    if !self.hla_hits.contains(&ivl.alt_name)
                        && get_hla_allele(&ivl.alt_name).is_some()
                    {
                        (&mut self.hla_hits).insert(ivl.alt_name.clone());
                    }
                }
            }
        }
    }

    fn write_to_hla_files(&mut self, entry: &[u8]) {
        if let Some(ref prefix) = self.hla_file_prefix {
            for contig in &self.hla_hits {
                let hla_file = self
                    .hla_files
                    .as_mut()
                    .unwrap()
                    .entry(contig.clone())
                    .or_insert_with(|| {
                        let file_name = format!("{}.{}.fq", prefix, contig.clone());
                        let file_path = Path::new(file_name.as_str());
                        match File::create(&file_path) {
                            Ok(file) => file,
                            Err(why) => panic!("couldn't open {}: {}", file_path.display(), why),
                        }
                    });
                hla_file.write_all(entry).unwrap();
            }
        }
    }

    // flushes the buffer if it contains segments of a different read than the current one
    fn flush(&mut self, out: &mut Writer) {
        if !self.read_buf.is_empty() {
            self.read_buf.iter().for_each(|r| out.write(r).unwrap());
            if let Some(ref hla_files) = self.hla_files {
                let first_read = self.read_buf.first().unwrap();
                let entry = format!(
                    "@{}/{}{}\n{}\n+\n{}\n",
                    std::str::from_utf8(first_read.qname()).unwrap(),
                    first_read.flags() >> 6 & 3,
                    if first_read.is_reverse() { "-" } else { "+" },
                    String::from_utf8(first_read.seq().as_bytes()).unwrap(),
                    first_read
                        .qual()
                        .iter()
                        .map(|q| char::try_from(u32::from(q + 33)).unwrap())
                        .collect::<String>(),
                );
                self.write_to_hla_files(entry.as_bytes());
            }
            self.read_buf.clear();
        }
        self.hla_hits.clear();
    }
}

struct Lifted {
    contig: Rc<String>,
    start: u32,
    end: u32,
    is_reverse: bool,
}

#[derive(Default)]
struct Hit {
    contig: Rc<String>,
    start: u32,
    end: u32,
    query_len: u32,
    is_reverse: bool,
    cigar: Vec<Cigar>,
    nm: u32,
    score: u32,
    index: usize,
    lifted: Option<Vec<Lifted>>,
    group: Option<u32>,
    ordering_contig: Rc<String>,
    ordering_start: u32,
    ordering_end: u32,
}

impl Hit {
    fn create(
        contig: Rc<String>,
        start: u32,
        end: u32,
        cigar: Vec<Cigar>,
        cigar_stats: &CigarStats,
        is_reverse: bool,
        nm: u32,
        index: usize,
    ) -> Hit {
        Hit {
            contig: contig.clone(),
            start,
            end,
            query_len: cigar_stats.match_len + cigar_stats.ins_len + cigar_stats.soft_len,
            is_reverse,
            cigar,
            nm,
            score: cigar_stats.calc_score(nm),
            index,
            lifted: None,
            group: None,
            ordering_contig: contig.clone(),
            ordering_start: start,
            ordering_end: end,
        }
    }

    fn set_lifted(&mut self, lifted: Vec<Lifted>) {
        if !lifted.is_empty() {
            let l0 = lifted.first().unwrap();
            self.ordering_contig = l0.contig.clone();
            self.ordering_start = l0.start;
            self.ordering_end = l0.end;
            self.lifted = Some(lifted);
        }
    }

    fn lifted_str(&self) -> Option<String> {
        if self.lifted.is_some() {
            Some(
                self.lifted
                    .as_ref()
                    .unwrap()
                    .iter()
                    .map(|l| {
                        format!(
                            "{},{},{};{}",
                            l.contig,
                            l.start,
                            l.end,
                            if l.is_reverse { "-" } else { "+" }
                        )
                    })
                    .collect(),
            )
        } else {
            None
        }
    }
}

impl PartialEq for Hit {
    fn eq(&self, other: &Self) -> bool {
        self.ordering_contig == other.ordering_contig
            && self.ordering_start == other.ordering_start
            && self.ordering_end == other.ordering_end
    }
}

impl Eq for Hit {}

impl Ord for Hit {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.ordering_contig == other.ordering_contig {
            self.ordering_start.cmp(&other.ordering_start)
        } else {
            self.ordering_contig.cmp(&other.ordering_contig)
        }
    }
}

impl PartialOrd for Hit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn main() {
    let opt: Opt = Opt::from_args();
    let mut contigs: Strings = Strings::new();
    // process alt-to-ref
    let alt_to_ref = AltToRef::parse(&opt.alt, &mut contigs);
    // process alignments
    let mut sam: Reader = opt
        .sam
        .as_ref()
        .map_or(Reader::from_stdin(), |path| Reader::from_path(path))
        .unwrap();
    let header = Header::from_template(sam.header());
    let format = HtsFormat::convert(if let Some(fmt) = opt.output_format {
        fmt
    } else {
        opt.output
            .as_ref()
            .or(opt.sam.as_ref())
            .and_then(|filename| {
                filename
                    .extension()
                    .map(|ext| HtsFormat::from_string_ignore_case(ext.to_str().unwrap()))
            })
            .unwrap_or(HtsFormat::Sam)
    });
    let mut out = if let Some(ref output) = opt.output {
        Writer::from_path(output, &header, format).unwrap()
    } else {
        Writer::from_stdout(&header, format).unwrap()
    };
    // buffers for holding reads for the same query
    let mut bufs = Buffers::create(opt.prefix, alt_to_ref.ref_indexes);
    // we need to access the header later (within the context of the mutable borrow of sam)
    let sam_header = sam.header().to_owned();
    for mut sam_read in sam
        .records()
        .map(|read| read.expect("Failure parsing sam file"))
    {
        // flush the buffer if it contains segments of a different read than the current one
        if bufs
            .last_read()
            .map(|r| {
                sam_read.qname() != r.qname()
                    || !(sam_read.is_first_in_template() == r.is_first_in_template()
                        && sam_read.is_last_in_template() == r.is_last_in_template())
            })
            .unwrap_or(false)
        {
            bufs.flush(&mut out);
        }
        // skip unmapped reads
        if sam_read.is_unmapped() {
            bufs.push_read(sam_read);
            continue;
        }
        // parse cigar
        let sam_cigar_stats = CigarStats::parse_cigar(&mut sam_read.cigar().to_vec(), false);
        // collect hits to HLA contig if any
        let sam_contig = contigs.intern_str(sam_read.contig());
        let sam_hit_start: u32 = u32::try_from(sam_read.pos()).unwrap();
        let sam_hit_end: u32 = sam_hit_start
            + sam_cigar_stats.match_len
            + sam_cigar_stats.del_len
            + sam_cigar_stats.skip_len;
        bufs.collect_hla_hits(&sam_contig, &sam_hit_start, &sam_hit_end);
        if sam_cigar_stats.hard_len > 0 {
            // cannot handle hard-clipped alignments
            bufs.push_read(sam_read);
            continue;
        }
        // parse NM (number of mutations) tag if it exists
        let sam_nm = sam_read
            .aux(b"NM")
            .map(|nm| if let Aux::U32(v) = nm { v } else { 0 })
            .unwrap_or(0);
        // parse XA (alternative hits) tag if it exits
        let sam_xa: Option<Vec<&str>> = sam_read.aux(b"XA").ok().and_then(|hits_str| {
            if let Aux::String(s) = hits_str {
                Some(s.split(';').filter(|s| !s.is_empty()).collect())
            } else {
                None
            }
        });
        // allocate a vector for the primary and alternative hits
        let mut hits = Vec::with_capacity(sam_xa.as_ref().map(|v| v.len()).unwrap_or(0) + 1);
        // add the primary hit
        hits[0] = Hit::create(
            sam_contig.clone(),
            sam_hit_start,
            sam_hit_end,
            sam_read.cigar().to_vec(),
            &sam_cigar_stats,
            sam_read.is_reverse(),
            sam_nm,
            0,
        );
        // add the alternative hits
        if let Some(xa_hits) = sam_xa {
            xa_hits.iter().enumerate().for_each(|(i, hit)| {
                let mut parts = hit.splitn(4, ',');
                let hit_contig = contigs.intern_str(parts.next().unwrap());
                let (hit_dir, hit_pos_str) = parts.next().unwrap().split_at(1);
                let hit_pos = hit_pos_str.parse().unwrap();
                let hit_cigar = CigarString::try_from(parts.next().unwrap())
                    .unwrap()
                    .to_vec();
                let hit_nm = parts.next().unwrap().parse().unwrap();
                let hit_cigar_stats = CigarStats::parse_cigar(&mut hit_cigar.to_vec(), false);
                hits[i + 1] = Hit::create(
                    hit_contig,
                    hit_pos,
                    hit_pos
                        + hit_cigar_stats.match_len
                        + hit_cigar_stats.del_len
                        + hit_cigar_stats.skip_len,
                    hit_cigar,
                    &hit_cigar_stats,
                    hit_dir == "-",
                    hit_nm,
                    i + 1,
                );
            });
        }
        // continue if none of the hits are to an ALT contig
        let alt_contigs = &alt_to_ref.alt_contigs;
        if !hits
            .iter()
            .any(|hit| alt_contigs.contains(hit.contig.as_ref()))
        {
            bufs.push_read(sam_read);
            continue;
        }
        // lift mapping positions to the primary assembly
        for hit in hits.iter_mut() {
            if !alt_contigs.contains(&hit.contig) {
                continue;
            }
            let alt_ivls = alt_to_ref
                .alt_indexes
                .get(&hit.contig)
                .unwrap()
                .find_overlaps(&hit.start, &hit.end);
            if alt_ivls.is_none() {
                continue;
            }
            hit.set_lifted(
                alt_ivls
                    .unwrap()
                    .iter()
                    .map(|ivl| {
                        let (ref_start, ref_end) = if ivl.is_reverse {
                            (
                                ivl.get_ref_pos(ivl.alt_len - hit.end),
                                ivl.get_ref_pos(ivl.alt_len - hit.start - 1).or(Some(0)),
                            )
                        } else {
                            (
                                ivl.get_ref_pos(hit.start),
                                ivl.get_ref_pos(hit.end - 1).or(Some(0)),
                            )
                        };
                        if ref_start.is_some() && ref_end.is_some() {
                            Some(Lifted {
                                contig: ivl.ref_name.clone(),
                                start: ref_start.unwrap() + ivl.pos,
                                end: ref_end.unwrap() + ivl.pos,
                                is_reverse: hit.is_reverse != ivl.is_reverse,
                            })
                        } else {
                            None
                        }
                    })
                    .filter_map(|x| x)
                    .collect(),
            );
        }
        // Group hits based on the lifted position on non-ALT sequences.
        // Find the index and group of the reported hit.
        // Find the size of the reported group.
        let mut reported_group: Option<u32> = None;
        let mut reported_index: Option<usize> = None;
        let mut num_reported: u32 = 0;
        if hits.len() > 1 {
            hits.sort_unstable();
            let mut cur_contig: Option<Rc<String>> = None;
            let mut max_end: u32 = 0;
            let mut cur_group: Option<u32> = None;
            for (i, hit) in hits.iter_mut().enumerate() {
                if cur_contig.as_ref().map_or(true, |c| c != &hit.contig) {
                    cur_contig = Some(hit.contig.clone());
                    cur_group = if let Some(g) = cur_group {
                        Some(g + 1)
                    } else {
                        Some(0)
                    };
                    max_end = 0;
                } else if hit.start > max_end {
                    cur_group = if let Some(g) = cur_group {
                        Some(g + 1)
                    } else {
                        Some(0)
                    };
                    max_end = hit.end;
                } else if hit.end > max_end {
                    max_end = hit.end;
                }
                hit.group = cur_group;
                if hit.index == 0 {
                    reported_group = cur_group;
                    reported_index = Some(i);
                }
            }
            num_reported = u32::try_from(
                hits.iter()
                    .filter(|hit| hit.group == reported_group)
                    .count(),
            )
            .unwrap();
        } else {
            let primary_hit = hits.first_mut().unwrap();
            if alt_contigs.contains(&primary_hit.contig) {
                bufs.push_read(sam_read);
                continue;
            }
            primary_hit.group = Some(0);
            reported_group = Some(0);
            reported_index = Some(0);
            num_reported = 1;
        }
        // re-estimate MAPQ
        let new_mapq = if num_reported > 1 {
            // first get all scores in each group
            let mut group_score = hits
                .iter()
                .map(|h| (h.group.unwrap(), h.score))
                .collect::<Vec<(u32, u32)>>();
            // sorting puts the highest score in the last position for each group
            group_score.sort();
            // collapsing to HashMap retains only the last value for each key
            // TODO: probably not good to rely on this behavior - replace with group_by when it's stable
            let group_max_map = group_score.into_iter().collect::<HashMap<_, _>>();
            let mut group_max = Vec::from_iter(group_max_map.iter());
            if group_max.len() > 1 {
                // sort groups by max score with the highest score first
                group_max.sort_by(|a, b| b.1.cmp(a.1));
            }
            let mapq: u8 = if group_max[0].0 != &reported_group.unwrap() {
                0
            } else if group_max.len() == 1 {
                60
            } else {
                u8::try_from(std::cmp::min(6 * (group_max[0].1 - group_max[1].1), 60)).unwrap()
            };
            match (
                alt_to_ref.alt_indexes.contains_key(&sam_contig),
                mapq.cmp(&sam_read.mapq()),
            ) {
                (true, Ordering::Less) => Some(mapq),
                (false, Ordering::Greater) => Some(mapq),
                _ => None,
            }
        } else {
            None
        };
        // collect HLA hits for the reported hit
        if hits[reported_index.unwrap()].ordering_contig == alt_to_ref.hla_contig {
            let mut start: Option<u32> = None;
            let mut end: Option<u32> = None;
            for hit in &hits {
                if hit.group.unwrap() == reported_group.unwrap() {
                    let hit_start = hit.ordering_start;
                    let hit_end = hit.ordering_end;
                    if start.is_none() || start.unwrap() > hit_start {
                        start = Some(hit_start);
                    }
                    if end.is_none() || end.unwrap() < hit_end {
                        end = Some(hit_end);
                    }
                }
            }
            if start.is_some() && end.is_some() {
                bufs.collect_hla_hits(&alt_to_ref.hla_contig, &start.unwrap(), &end.unwrap());
            }
        }
        // adjust the mapQ of the primary hits
        let primary_lifted = hits[reported_index.unwrap()].lifted.as_ref();
        if primary_lifted.is_none() || primary_lifted.unwrap().len() <= 1 {
            for buf_hit in bufs.read_iter_mut() {
                let is_overlap = if primary_lifted.is_some() {
                    let buf_hit_contig =
                        Rc::new(std::str::from_utf8(buf_hit.qname()).unwrap().to_owned());
                    let p = primary_lifted.unwrap().first().unwrap();
                    if p.contig != buf_hit_contig {
                        false
                    } else if p.is_reverse != buf_hit.is_reverse() {
                        false
                    } else {
                        let start = u32::try_from(buf_hit.pos() - 1).unwrap();
                        if start >= p.end {
                            false
                        } else {
                            let cigar_len: u32 = hits[reported_index.unwrap()]
                                .cigar
                                .iter()
                                .map(|c| match c {
                                    Cigar::Match(l) | Cigar::Del(l) | Cigar::RefSkip(l) => *l,
                                    _ => 0,
                                })
                                .sum();
                            p.start < (start + cigar_len)
                        }
                    }
                } else {
                    false
                };
                // parse the OM (original MAPQ) and PA (primary assembly) tags
                let om = buf_hit.aux(b"OM").ok().and_then(|om_aux| {
                    if let Aux::U8(om) = om_aux {
                        Some(om)
                    } else {
                        None
                    }
                });
                let pa = buf_hit
                    .aux(b"PA")
                    .ok()
                    .and_then(|pa_aux| {
                        if let Aux::Float(pa) = pa_aux {
                            Some(pa)
                        } else {
                            None
                        }
                    })
                    .unwrap_or(10.);
                if is_overlap {
                    let final_mapq = new_mapq.unwrap_or(sam_read.mapq());
                    match om {
                        Some(o) if o > 0 && o < final_mapq => buf_hit.set_mapq(o),
                        Some(_) => buf_hit.set_mapq(final_mapq),
                        None if buf_hit.mapq() > final_mapq => buf_hit.set_mapq(final_mapq),
                        _ => (),
                    };
                } else if pa < opt.min_pa_ratio {
                    if om.is_none() {
                        buf_hit.push_aux(b"OM", Aux::U8(buf_hit.mapq()));
                    }
                    buf_hit.set_mapq(0);
                }
            }
        }
        // update the read and add it to the buffer
        if num_reported > 1 {
            sam_read.push_aux(b"OM", Aux::U8(sam_read.mapq())).unwrap();
        }
        hits[reported_index.unwrap()]
            .lifted_str()
            .into_iter()
            .for_each(|s| sam_read.push_aux(b"LT", Aux::String(s.as_str())).unwrap());
        sam_read.set_mapq(new_mapq.unwrap());
        // create records from the hits generated from the XA tag
        let mut xa_hit_buf: Option<Vec<Record>> = None;
        if !hits.is_empty() {
            xa_hit_buf = Some(Vec::with_capacity(hits.len()));
            for (i, hit) in hits.iter().enumerate() {
                if i == reported_index.unwrap()
                    || hit.group != reported_group
                    || !alt_to_ref.alt_indexes.contains_key(&hit.contig)
                {
                    continue;
                }
                let mut hit_read = Record::new();
                if hit.is_reverse == hits[reported_index.unwrap()].is_reverse {
                    hit_read.set(
                        sam_read.qname(),
                        Some(&CigarString::from(hit.cigar.clone())),
                        sam_read.seq().as_bytes().as_slice(),
                        sam_read.qual(),
                    );
                    hit_read.set_flags(sam_read.flags() | 0x800);
                } else {
                    let mut qual_reversed = sam_read.qual().to_owned();
                    qual_reversed.reverse();
                    hit_read.set(
                        sam_read.qname(),
                        Some(&CigarString::from(hit.cigar.clone())),
                        revcomp(sam_read.seq().as_bytes()).as_slice(),
                        qual_reversed.as_slice(),
                    );
                    hit_read.set_flags((sam_read.flags() ^ 0x10) | 0x800);
                }
                hit_read.set_tid(
                    i32::try_from(
                        sam_header
                            .tid(hit.contig.as_bytes())
                            .expect(format!("contig {} not found in header", hit.contig).as_str()),
                    )
                    .unwrap(),
                );
                hit_read.set_pos(i64::from(hit.start));
                hit_read.set_mapq(new_mapq.unwrap());
                if sam_read.mtid() == sam_read.tid() && sam_contig != hit.contig {
                    hit_read.set_mtid(sam_read.tid());
                } else {
                    hit_read.set_mtid(sam_read.mtid());
                }
                hit_read.set_mpos(sam_read.mpos());
                hit_read.set_insert_size(sam_read.insert_size());
                hit_read.push_aux(b"NM", Aux::U32(hit.nm));
                hit.lifted_str().iter().for_each(|s| {
                    hit_read.push_aux(b"LT", Aux::String(s));
                });
                if let Some(rg) = sam_read.aux(b"RG").ok() {
                    hit_read.push_aux(b"RG", rg);
                }
                bufs.push_read(hit_read);
            }
        }
        // push the primary read and any XA-derived reads
        bufs.push_read(sam_read);
        if xa_hit_buf.is_some() {
            bufs.append_reads(&mut xa_hit_buf.unwrap())
        }
    }
    bufs.flush(&mut out);
}
