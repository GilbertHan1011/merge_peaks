use std::{collections::HashMap, path::PathBuf};
use anyhow::{anyhow, bail, ensure, Context, Result};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use clap::Parser;
use bed_utils::bed::{io::Reader, BEDLike, MergeBed, NarrowPeak};
use bed_utils::extsort::ExternalSorterBuilder;
use serde::{Deserialize, Serialize};
use flate2::read::MultiGzDecoder;

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

    #[arg(long, default_value_t = true)]
    normalize: bool,
}

/// Calculate the total signal of a set of peaks (sum of scores).
/// Treats missing scores as zero.
fn calculate_total_signal(peaks: &[NarrowPeak]) -> u64 {
    peaks
        .iter()
        .filter_map(|peak| peak.p_value)
        .map(|score| score)
        .sum()
}

/// Normalize the peak scores to "score per million" (SPM).
/// Consumes the peaks and returns a new Vec<NarrowPeak> with scores recalculated.

fn spm(mut peaks: Vec<NarrowPeak>) -> Result<Vec<NarrowPeak>> {
    let total_signal = calculate_total_signal(&peaks);

    if total_signal == 0 {
        // Prevent division by zero; just return peaks as-is.
        return Ok(peaks);
    }

    for peak in peaks.iter_mut() {
        let normalized_score = (peak.p_value / total_signal * 1_000_000).unwrap_or(0.0);

        peak.p_value = normalized_score;
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
    value.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string())
}

pub fn merge_peaks_narrowpeak(
    peaks: HashMap<String, Vec<NarrowPeak>>,
    chrom_sizes: HashMap<String, u64>,
    half_width: u64,
) -> Result<Vec<NarrowPeak>> {
    let peak_list: Vec<_> = peaks
        .into_iter()
        .collect::<Vec<_>>();
    let chrom_sizes = chrom_sizes.into_iter().collect();
    let merged_peaks: Vec<_> = merge_peaks(
            peak_list.iter().flat_map(|x| x.1.clone()), 
            half_width)
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


pub fn merge_peaks<I>(peaks: I, half_window_size: u64) -> impl Iterator<Item = Vec<NarrowPeak>>
where
    I: Iterator<Item = NarrowPeak>,
{
    fn iterative_merge(mut peaks: Vec<NarrowPeak>,score_threshold: u8,overlap_threshold: u8) -> Vec<NarrowPeak> {
        let mut result = Vec::new();
        while !peaks.is_empty() {
            let best_peak = peaks.iter()
                .max_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap()).unwrap()
                .clone();
            previous_size = peaks.len();
            peaks = peaks.into_iter().filter(|x| x.n_overlap(&best_peak) == 0).collect();
            latter_size = previous_size - peaks.len() - 1; // Remove self from the count
            if latter_size >= overlap_threshold && best_peak.p_value >= score_threshold {
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
        .build().unwrap()
        .sort_by(input, BEDLike::compare).unwrap()
        .map(|x| x.unwrap())
        .merge_sorted_bed_with(|peaks| iterative_merge(peaks, score_threshold, overlap_threshold))
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

#[derive(Serialize, Deserialize, Debug)]
struct ChromSizes {
    genome: String,
    chromosomes: HashMap<String, u64>,
}

fn write_bed(path: &PathBuf, peaks: &[NarrowPeak]) -> Result<()> {
    let mut file = File::create(path)
        .with_context(|| format!("Failed to create output file: {:?}", path))?;
    for peak in peaks {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}",
            peak.chrom(),
            peak.start(),
            peak.end(),
            format_optional_string(peak.name.clone()),
            format_optional_value(peak.score),
            format_optional_string(peak
                .strand
                .map(|strand| strand.to_string())),
            peak.signal_value,
            format_optional_value(peak.p_value),
            format_optional_value(peak.q_value),
            peak.peak
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
    ensure!(
        !trimmed.is_empty(),
        "Chrom sizes file {:?} is empty",
        path
    );

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
    let merged = merge_peaks_narrowpeak(peaks_by_source, chrom_sizes, cli.half_width)?;

    write_bed(&cli.output, &merged)?;
    eprintln!("Merged {} peaks written to {:?}", merged.len(), cli.output);

    Ok(())
}
