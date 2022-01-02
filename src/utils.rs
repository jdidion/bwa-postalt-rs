use std::borrow::Borrow;
use std::collections::HashMap;
use std::rc::Rc;

use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::bam::record::Cigar;

/// Data structure for interning commonly used strings, e.g. contig names.
pub struct Strings {
    strings: HashMap<Vec<u8>, Rc<String>>,
}

impl Strings {
    pub fn new() -> Strings {
        Strings {
            strings: HashMap::new(),
        }
    }

    pub fn intern(&mut self, bytes: &[u8]) -> Rc<String> {
        if self.strings.contains_key(bytes) {
            self.strings.get(bytes).unwrap().clone()
        } else {
            let key = bytes.to_owned();
            let value = Rc::new(std::str::from_utf8(bytes).unwrap().to_owned());
            let retval = value.clone();
            self.strings.insert(key, value);
            retval
        }
    }

    pub fn intern_str(&mut self, s: &str) -> Rc<String> {
        self.intern(s.as_bytes())
    }
}

// borrowed from rust-bio ---
lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}

lazy_static! {
    // regular expression for matching HLA contig names
    static ref HLA_RE: Regex = {
        Regex::new(r"^(HLA-[^\s*]+)\*\d+").unwrap()
    };
}

pub fn get_hla_allele(contig: &str) -> Option<&str> {
    HLA_RE
        .captures(contig)
        .map(|caps| caps.get(1).unwrap().as_str())
}

#[derive(Default)]
pub struct CigarStats {
    pub match_len: u32,
    pub ins_len: u32,
    pub ins_num: u32,
    pub del_len: u32,
    pub del_num: u32,
    pub skip_len: u32,
    pub hard_len: u32,
    pub soft_len: u32,
}

struct ScoreOpts {
    a: u32,
    b: u32,
    o: u32,
    e: u32,
    x: f32,
}

static SCORE_OPTS: ScoreOpts = ScoreOpts {
    a: 1,
    b: 4,
    o: 6,
    e: 1,
    x: 0.499,
};

impl CigarStats {
    pub fn parse_cigar(cigar: &mut Vec<Cigar>, replace_hard: bool) -> CigarStats {
        let mut stats = CigarStats {
            ..Default::default()
        };
        for i in 0..cigar.len() {
            match cigar[i] {
                Cigar::Match(l) => stats.match_len += l,
                Cigar::Ins(l) => {
                    stats.ins_len += l;
                    stats.ins_num += 1;
                }
                Cigar::Del(l) => {
                    stats.del_len += l;
                    stats.del_num += 1;
                }
                Cigar::RefSkip(l) => {
                    stats.skip_len += l;
                }
                Cigar::HardClip(l) => {
                    stats.hard_len += l;
                    if replace_hard {
                        // convert hard clip to soft clip
                        cigar[i] = Cigar::SoftClip(l);
                    }
                }
                Cigar::SoftClip(l) => {
                    stats.soft_len += l;
                }
                _ => (),
            }
        }
        stats
    }

    pub fn calc_score(&self, nm: u32) -> u32 {
        let score: f32 = ((SCORE_OPTS.a * self.match_len)
            - ((SCORE_OPTS.a + SCORE_OPTS.b) * (nm - self.del_len - self.ins_len))
            - (SCORE_OPTS.o * (self.del_num + self.ins_num))
            - (SCORE_OPTS.e * (self.del_len + self.ins_len))) as f32
            / (SCORE_OPTS.a as f32 + SCORE_OPTS.x);
        score.floor() as u32
    }
}
