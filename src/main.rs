use anyhow::{anyhow, bail, ensure, Context, Result};
use bed_utils::bed::{io::Reader, BEDLike, BroadPeak, MergeBed, NarrowPeak};
use bed_utils::extsort::ExternalSorterBuilder;
use clap::{Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::{collections::HashMap, path::PathBuf};

#[derive(Clone, Debug, ValueEnum)]
enum MergeMode {
    /// Existing narrowPeak/ATAC-style behavior: resize each peak to summit +/- half-width.
    SummitWindow,
    /// broadPeak behavior: preserve original intervals and keep regions supported by N samples.
    BroadConsensus,
}

#[derive(Parser, Debug)]
#[command(author, version, about = "Merge narrowPeak/BED peak files")]
struct Cli {
    /// BED or narrowPeak files to merge
    #[arg(value_name = "BED", num_args = 1.., value_hint = clap::ValueHint::FilePath)]
    bed_files: Vec<PathBuf>,

    /// Reference genome FASTA file used to derive chromosome sizes
    #[arg(long = "fasta", value_name = "FASTA", value_hint = clap::ValueHint::FilePath)]
    genome_fasta: Option<PathBuf>,

    /// Precomputed chromosome sizes (TSV or JSON)
    #[arg(long = "chrom-sizes", value_name = "FILE", value_hint = clap::ValueHint::FilePath)]
    chrom_sizes: Option<PathBuf>,

    /// Output BED path for merged peaks
    #[arg(short, long, value_name = "BED", value_hint = clap::ValueHint::FilePath)]
    output: PathBuf,

    /// Half window size (bp) used to expand summits before merging
    #[arg(long, default_value_t = 250)]
    half_width: u64,

    /// Merge strategy. Use summit-window for narrowPeak/ATAC and broad-consensus for broadPeak.
    #[arg(long, value_enum, default_value_t = MergeMode::SummitWindow)]
    merge_mode: MergeMode,

    /// Absolute support threshold for broad-consensus mode. Overrides --min-support-frac.
    #[arg(long)]
    min_support: Option<usize>,

    /// Fraction of input files that must cover a broad-consensus interval.
    #[arg(long, default_value_t = 0.2)]
    min_support_frac: f64,

    /// Maximum gap (bp) to bridge between adjacent broad-consensus intervals.
    #[arg(long, default_value_t = 1000)]
    max_gap: u64,

    /// Minimum width (bp) retained in broad-consensus output after gap bridging.
    #[arg(long, default_value_t = 1000)]
    min_width: u64,

    #[arg(long, default_value_t = true)]
    normalize: bool,
    #[arg(long, default_value_t = 0.0)]
    score_threshold: f64,
    #[arg(long, default_value_t = 0)]
    overlap_threshold: usize,
}

/// Normalize the peak scores to "score per million" (SPM).
/// Consumes the peaks and returns a new Vec<NarrowPeak> with scores recalculated.

fn spm(mut peaks: Vec<NarrowPeak>) -> Result<Vec<NarrowPeak>> {
    let total_signal: f64 = peaks.iter().filter_map(|p| p.p_value).sum();

    if total_signal == 0.0 {
        // Prevent division by zero; just return peaks as-is.
        return Ok(peaks);
    }
    let factor = 1_000_000.0 / total_signal;
    for peak in &mut peaks {
        if let Some(score) = peak.p_value.as_mut() {
            *score *= factor;
        }
    }

    Ok(peaks)
}

fn format_optional_string(value: Option<String>) -> String {
    value.unwrap_or_else(|| ".".to_string())
}

fn format_optional_value<T>(value: Option<T>) -> String
where
    T: ToString,
{
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| ".".to_string())
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusDomain {
    chrom: String,
    start: u64,
    end: u64,
    support: usize,
}

fn resolve_min_support(explicit: Option<usize>, fraction: f64, n_sources: usize) -> Result<usize> {
    ensure!(n_sources > 0, "No input peak files were provided");
    if let Some(min_support) = explicit {
        ensure!(min_support > 0, "--min-support must be greater than zero");
        ensure!(
            min_support <= n_sources,
            "--min-support ({}) cannot exceed the number of input files ({})",
            min_support,
            n_sources
        );
        return Ok(min_support);
    }

    ensure!(
        fraction > 0.0 && fraction <= 1.0,
        "--min-support-frac must be in the interval (0, 1]"
    );
    let threshold = ((fraction * n_sources as f64).ceil() as usize).max(2);
    Ok(threshold.min(n_sources))
}

pub fn merge_peaks_narrowpeak(
    peaks: HashMap<String, Vec<NarrowPeak>>,
    chrom_sizes: HashMap<String, u64>,
    half_width: u64,
    score_threshold: f64,
    overlap_threshold: usize,
) -> Result<Vec<NarrowPeak>> {
    let peak_list: Vec<_> = peaks.into_iter().collect::<Vec<_>>();
    let chrom_sizes = chrom_sizes.into_iter().collect();
    let merged_peaks: Vec<_> = merge_peaks(
        peak_list.iter().flat_map(|x| x.1.clone()),
        half_width,
        score_threshold,
        overlap_threshold,
    )
    .flatten()
    .map(|x| clip_peak(x, &chrom_sizes))
    .collect();

    Ok(merged_peaks)
}

pub fn clip_peak(mut peak: NarrowPeak, chrom_sizes: &HashMap<String, u64>) -> NarrowPeak {
    let chr = peak.chrom();
    let max_len = *chrom_sizes
        .get(chr)
        .expect(&format!("Size missing for chromosome: {}", chr));
    let new_start = peak.start().max(0).min(max_len);
    let new_end = peak.end().min(max_len);
    peak.set_start(new_start);
    peak.set_end(new_end);
    peak.peak = (new_start + peak.peak).min(new_end) - new_start;
    peak
}

pub fn merge_peaks_broad_consensus(
    peaks: HashMap<String, Vec<BroadPeak>>,
    chrom_sizes: HashMap<String, u64>,
    min_support: usize,
    max_gap: u64,
    min_width: u64,
) -> Result<Vec<ConsensusDomain>> {
    ensure!(min_support > 0, "min_support must be greater than zero");

    let mut events_by_chrom: HashMap<String, Vec<(u64, i64)>> = HashMap::new();
    for peak in peaks.into_values().flatten() {
        let Some(chrom_size) = chrom_sizes.get(peak.chrom()) else {
            bail!("Size missing for chromosome: {}", peak.chrom());
        };
        let start = peak.start().min(*chrom_size);
        let end = peak.end().min(*chrom_size);
        if start >= end {
            continue;
        }
        let events = events_by_chrom.entry(peak.chrom().to_string()).or_default();
        events.push((start, 1));
        events.push((end, -1));
    }

    let mut segments = Vec::new();
    for (chrom, mut events) in events_by_chrom {
        events.sort_unstable_by_key(|(pos, _)| *pos);
        let mut current_support = 0i64;
        let mut previous_pos: Option<u64> = None;
        let mut idx = 0usize;

        while idx < events.len() {
            let pos = events[idx].0;
            if let Some(start) = previous_pos {
                if start < pos && current_support as usize >= min_support {
                    segments.push(ConsensusDomain {
                        chrom: chrom.clone(),
                        start,
                        end: pos,
                        support: current_support as usize,
                    });
                }
            }

            let mut delta = 0i64;
            while idx < events.len() && events[idx].0 == pos {
                delta += events[idx].1;
                idx += 1;
            }
            current_support += delta;
            previous_pos = Some(pos);
        }
    }

    segments.sort_unstable_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
    });

    let mut merged: Vec<ConsensusDomain> = Vec::new();
    for segment in segments {
        if let Some(last) = merged.last_mut() {
            if last.chrom == segment.chrom && segment.start <= last.end.saturating_add(max_gap) {
                last.end = last.end.max(segment.end);
                last.support = last.support.max(segment.support);
                continue;
            }
        }
        merged.push(segment);
    }

    Ok(merged
        .into_iter()
        .filter(|domain| domain.end.saturating_sub(domain.start) >= min_width)
        .collect())
}

pub fn merge_peaks<I>(
    peaks: I,
    half_window_size: u64,
    score_threshold: f64,
    overlap_threshold: usize,
) -> impl Iterator<Item = Vec<NarrowPeak>>
where
    I: Iterator<Item = NarrowPeak>,
{
    fn iterative_merge(
        mut peaks: Vec<NarrowPeak>,
        score_threshold: f64,
        overlap_threshold: usize,
    ) -> Vec<NarrowPeak> {
        let mut result = Vec::new();
        while !peaks.is_empty() {
            let best_peak = peaks
                .iter()
                .max_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap())
                .unwrap()
                .clone();
            let previous_size = peaks.len();
            peaks = peaks
                .into_iter()
                .filter(|x| x.n_overlap(&best_peak) == 0)
                .collect();
            let latter_size = previous_size - peaks.len() - 1; // Remove self from the count
            if latter_size >= overlap_threshold
                && best_peak.p_value.unwrap_or(0.0) >= score_threshold
            {
                result.push(best_peak);
            }
        }
        result
    }

    let input = peaks.map(move |mut x| {
        let summit = x.start() + x.peak;
        x.start = summit.saturating_sub(half_window_size);
        x.end = summit + half_window_size + 1;
        x.peak = summit - x.start;
        x
    });
    ExternalSorterBuilder::new()
        .with_compression(2)
        .build()
        .unwrap()
        .sort_by(input, BEDLike::compare)
        .unwrap()
        .map(|x| x.unwrap())
        .merge_sorted_bed_with(move |peaks| {
            iterative_merge(peaks, score_threshold, overlap_threshold)
        })
}

fn read_bed(path: &PathBuf) -> Result<Vec<NarrowPeak>> {
    let bed_file_open = File::open(path)?;
    let mut bed_reader = Reader::new(bed_file_open, None);

    let mut narrow_peaks = Vec::new();

    for bed_result in bed_reader.records::<NarrowPeak>() {
        let bed_record = bed_result?;
        // Convert BED<10> to NarrowPeak using all ten fields
        narrow_peaks.push(bed_record);
    }

    Ok(narrow_peaks)
}

fn read_broad_bed(path: &PathBuf) -> Result<Vec<BroadPeak>> {
    let bed_file_open = File::open(path)?;
    let mut bed_reader = Reader::new(bed_file_open, None);

    let mut broad_peaks = Vec::new();
    for bed_result in bed_reader.records::<BroadPeak>() {
        broad_peaks.push(bed_result?);
    }

    Ok(broad_peaks)
}

#[derive(Serialize, Deserialize, Debug)]
struct ChromSizes {
    genome: String,
    chromosomes: HashMap<String, u64>,
}

fn write_bed(path: &PathBuf, peaks: &[NarrowPeak]) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to create output file: {:?}", path))?;
    for peak in peaks {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}",
            peak.chrom(),
            peak.start(),
            peak.end(),
            format_optional_string(peak.name.clone()),
            format_optional_value(peak.score),
            format_optional_string(peak.strand.map(|strand| strand.to_string())),
            peak.signal_value,
            format_optional_value(peak.p_value),
            format_optional_value(peak.q_value),
            peak.peak
        )?;
    }
    Ok(())
}

fn write_consensus_bed(path: &PathBuf, domains: &[ConsensusDomain]) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to create output file: {:?}", path))?;
    for (idx, domain) in domains.iter().enumerate() {
        let score = domain.support.min(1000);
        writeln!(
            file,
            "{}\t{}\t{}\tconsensus_domain_{}_support_{}\t{}\t.\t{}\t-1\t-1",
            domain.chrom,
            domain.start,
            domain.end,
            idx + 1,
            domain.support,
            score,
            domain.support
        )?;
    }
    Ok(())
}

fn read_fasta_chrom_sizes(fasta_path: &PathBuf) -> Result<HashMap<String, u64>> {
    let file = File::open(fasta_path)
        .with_context(|| format!("Failed to open FASTA file: {:?}", fasta_path))?;
    let mut reader = fasta_reader(file)?;

    let mut chrom_sizes = HashMap::new();
    let mut current_chrom: Option<String> = None;
    let mut current_length = 0u64;
    let mut buffer = Vec::new();

    loop {
        buffer.clear();
        let bytes_read = reader.read_until(b'\n', &mut buffer)?;
        if bytes_read == 0 {
            break;
        }

        if buffer.starts_with(b">") {
            if let Some(chrom) = current_chrom.take() {
                chrom_sizes.insert(chrom, current_length);
            }

            let header = String::from_utf8_lossy(&buffer[1..]).trim().to_string();
            let chrom = header
                .split_whitespace()
                .next()
                .unwrap_or(header.as_str())
                .to_string();
            current_chrom = Some(chrom);
            current_length = 0;
        } else {
            let seq_len = buffer
                .iter()
                .filter(|b| !matches!(b, b'\n' | b'\r'))
                .count() as u64;
            current_length += seq_len;
        }
    }

    if let Some(chrom) = current_chrom {
        chrom_sizes.insert(chrom, current_length);
    }

    Ok(chrom_sizes)
}

fn fasta_reader(file: File) -> Result<Box<dyn BufRead>> {
    let mut file = file;
    let mut magic = [0u8; 2];
    let bytes_read = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    let reader: Box<dyn BufRead> = if bytes_read == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(reader)
}

fn read_chrom_sizes_file(path: &PathBuf) -> Result<HashMap<String, u64>> {
    let raw = fs::read_to_string(path)
        .with_context(|| format!("Failed to read chrom sizes file: {:?}", path))?;
    let trimmed = raw.trim();
    ensure!(!trimmed.is_empty(), "Chrom sizes file {:?} is empty", path);

    if let Ok(wrapper) = serde_json::from_str::<ChromSizes>(trimmed) {
        return Ok(wrapper.chromosomes);
    }

    if let Ok(map) = serde_json::from_str::<HashMap<String, u64>>(trimmed) {
        return Ok(map);
    }

    let mut sizes = HashMap::new();
    for (idx, line) in raw.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let chrom = fields
            .next()
            .ok_or_else(|| anyhow!("Missing chromosome name on line {}", idx + 1))?;
        let size_str = fields
            .next()
            .ok_or_else(|| anyhow!("Missing size for chromosome {} on line {}", chrom, idx + 1))?;
        ensure!(
            fields.next().is_none(),
            "Too many columns on line {}: {}",
            idx + 1,
            line
        );
        let size = size_str.parse::<u64>().with_context(|| {
            format!(
                "Failed to parse size '{}' for chromosome {} on line {}",
                size_str,
                chrom,
                idx + 1
            )
        })?;

        sizes.insert(chrom.to_string(), size);
    }

    if sizes.is_empty() {
        bail!(
            "Unable to parse chromosome sizes from {:?}; supported formats are TSV ('chr size') or JSON",
            path
        );
    }

    Ok(sizes)
}

fn load_chrom_sizes(cli: &Cli) -> Result<HashMap<String, u64>> {
    if let Some(path) = &cli.chrom_sizes {
        read_chrom_sizes_file(path)
    } else if let Some(fasta) = &cli.genome_fasta {
        read_fasta_chrom_sizes(fasta)
            .with_context(|| format!("Failed to read FASTA file: {:?}", fasta))
    } else {
        bail!("Please provide either --chrom-sizes or --fasta");
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    ensure!(
        !cli.bed_files.is_empty(),
        "Please provide at least one BED file"
    );
    ensure!(
        cli.chrom_sizes.is_some() || cli.genome_fasta.is_some(),
        "Provide --chrom-sizes or --fasta to supply chromosome sizes"
    );

    let chrom_sizes = load_chrom_sizes(&cli)?;

    match cli.merge_mode {
        MergeMode::SummitWindow => {
            let mut peaks_by_source = HashMap::new();
            for bed_path in &cli.bed_files {
                let peaks = read_bed(bed_path)?;
                peaks_by_source.insert(bed_path.display().to_string(), peaks);
            }
            if cli.normalize {
                for peaks in peaks_by_source.values_mut() {
                    *peaks = spm(peaks.clone())?;
                }
            }
            let merged = merge_peaks_narrowpeak(
                peaks_by_source,
                chrom_sizes,
                cli.half_width,
                cli.score_threshold,
                cli.overlap_threshold,
            )?;

            write_bed(&cli.output, &merged)?;
            eprintln!("Merged {} peaks written to {:?}", merged.len(), cli.output);
        }
        MergeMode::BroadConsensus => {
            let min_support =
                resolve_min_support(cli.min_support, cli.min_support_frac, cli.bed_files.len())?;
            let mut peaks_by_source = HashMap::new();
            for bed_path in &cli.bed_files {
                let peaks = read_broad_bed(bed_path)?;
                peaks_by_source.insert(bed_path.display().to_string(), peaks);
            }
            let merged = merge_peaks_broad_consensus(
                peaks_by_source,
                chrom_sizes,
                min_support,
                cli.max_gap,
                cli.min_width,
            )?;

            write_consensus_bed(&cli.output, &merged)?;
            eprintln!(
                "Merged {} broad consensus domains written to {:?} with min_support={}",
                merged.len(),
                cli.output,
                min_support
            );
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    // test spm
    #[test]
    fn test_spm() {
        let peak1 =
            NarrowPeak::from_str("chr1\t100\t200\tname1\t0\t.\t100\t100\t100\t100").unwrap();
        let peak2 =
            NarrowPeak::from_str("chr1\t100\t200\tname2\t0\t.\t100\t100\t100\t100").unwrap();
        let peaks = vec![peak1, peak2];
        let normalized = spm(peaks).unwrap();
        assert_eq!(normalized.len(), 2);
        assert_eq!(normalized[0].p_value, Some(500_000.0));
        assert_eq!(normalized[1].p_value, Some(500_000.0));
    }

    #[test]
    fn test_merge_peaks() {
        let peaks = vec![
            NarrowPeak::from_str("chr1\t100\t200\tname1\t0\t.\t100\t100\t100\t100").unwrap(),
            NarrowPeak::from_str("chr1\t105\t205\tname2\t0\t.\t100\t200\t100\t100").unwrap(),
            NarrowPeak::from_str("chr1\t110\t210\tname3\t0\t.\t100\t300\t100\t100").unwrap(),
            NarrowPeak::from_str("chr1\t115\t215\tname4\t0\t.\t100\t100\t100\t100").unwrap(),
            NarrowPeak::from_str("chr2\t115\t215\tname5\t0\t.\t100\t100\t100\t100").unwrap(),
        ];
        let merged: Vec<NarrowPeak> = merge_peaks(peaks.into_iter(), 250, 0.0, 1)
            .flatten()
            .collect();
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].p_value, Some(300.0));
    }

    #[test]
    fn test_resolve_min_support() {
        assert_eq!(resolve_min_support(None, 0.2, 98).unwrap(), 20);
        assert_eq!(resolve_min_support(None, 0.2, 64).unwrap(), 13);
        assert_eq!(resolve_min_support(None, 0.2, 1).unwrap(), 1);
        assert_eq!(resolve_min_support(Some(3), 0.2, 10).unwrap(), 3);
    }

    #[test]
    fn test_broad_consensus_basic() {
        let sample1 = vec![
            BroadPeak::from_str("chr1\t100\t300\tpeak1\t0\t.\t10\t20\t30").unwrap(),
            BroadPeak::from_str("chr1\t500\t800\tpeak2\t0\t.\t10\t20\t30").unwrap(),
        ];
        let sample2 = vec![
            BroadPeak::from_str("chr1\t200\t400\tpeak3\t0\t.\t10\t20\t30").unwrap(),
            BroadPeak::from_str("chr1\t700\t900\tpeak4\t0\t.\t10\t20\t30").unwrap(),
        ];
        let sample3 = vec![BroadPeak::from_str("chr1\t250\t350\tpeak5\t0\t.\t10\t20\t30").unwrap()];
        let mut peaks = HashMap::new();
        peaks.insert("sample1".to_string(), sample1);
        peaks.insert("sample2".to_string(), sample2);
        peaks.insert("sample3".to_string(), sample3);
        let chrom_sizes = HashMap::from([("chr1".to_string(), 1000)]);

        let merged = merge_peaks_broad_consensus(peaks, chrom_sizes, 2, 0, 1).unwrap();

        assert_eq!(
            merged,
            vec![
                ConsensusDomain {
                    chrom: "chr1".to_string(),
                    start: 200,
                    end: 350,
                    support: 3,
                },
                ConsensusDomain {
                    chrom: "chr1".to_string(),
                    start: 700,
                    end: 800,
                    support: 2,
                },
            ]
        );
    }

    #[test]
    fn test_broad_consensus_gap_bridge_and_min_width() {
        let sample1 = vec![
            BroadPeak::from_str("chr1\t100\t200\tpeak1\t0\t.\t10\t20\t30").unwrap(),
            BroadPeak::from_str("chr1\t260\t400\tpeak2\t0\t.\t10\t20\t30").unwrap(),
        ];
        let sample2 = vec![
            BroadPeak::from_str("chr1\t120\t220\tpeak3\t0\t.\t10\t20\t30").unwrap(),
            BroadPeak::from_str("chr1\t280\t420\tpeak4\t0\t.\t10\t20\t30").unwrap(),
        ];
        let mut peaks = HashMap::new();
        peaks.insert("sample1".to_string(), sample1);
        peaks.insert("sample2".to_string(), sample2);
        let chrom_sizes = HashMap::from([("chr1".to_string(), 1000)]);

        let merged = merge_peaks_broad_consensus(peaks, chrom_sizes, 2, 100, 250).unwrap();

        assert_eq!(
            merged,
            vec![ConsensusDomain {
                chrom: "chr1".to_string(),
                start: 120,
                end: 400,
                support: 2,
            }]
        );
    }
}
