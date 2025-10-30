use clap::Parser;
use flate2::read::GzDecoder;
use log::{debug, error, info, warn};
use noodles::bgzf;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path as FilePath;

#[inline]
fn base_to_code(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

#[inline]
fn base_to_index(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

#[inline]
fn fill_counts(counts: &mut [usize; 5], window: &[u8]) {
    counts.fill(0);
    for &base in window {
        counts[base_to_index(base)] += 1;
    }
}

// Function to calculate percentiles from sorted data
fn calculate_percentile(sorted_data: &[f64], percentile: f64) -> f64 {
    if sorted_data.is_empty() {
        return 0.0;
    }

    let index = (percentile / 100.0) * (sorted_data.len() - 1) as f64;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;

    if lower == upper {
        sorted_data[lower]
    } else {
        let weight = index - lower as f64;
        sorted_data[lower] * (1.0 - weight) + sorted_data[upper] * weight
    }
}

// Function to calculate IQR and derive threshold
fn calculate_iqr_threshold(complexity_values: &[f64], iqr_multiplier: f64) -> (f64, f64, f64, f64) {
    if complexity_values.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let mut sorted = complexity_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let q1 = calculate_percentile(&sorted, 25.0);
    let q3 = calculate_percentile(&sorted, 75.0);
    let iqr = q3 - q1;

    // Calculate threshold using Q1 - multiplier * IQR
    let threshold = q1 - iqr_multiplier * iqr;

    (q1, q3, iqr, threshold)
}

// Shannon entropy complexity function
pub fn shannon_entropy_complexity(seq: &[u8], window_size: usize, step: usize) -> Vec<f64> {
    if seq.is_empty() || window_size == 0 {
        return Vec::new();
    }

    let seq_len = seq.len();
    if seq_len <= window_size {
        let mut counts = [0usize; 5];
        fill_counts(&mut counts, seq);
        return vec![entropy_from_counts(&counts, seq_len)];
    }

    let step = step.max(1);
    let mut results = Vec::with_capacity((seq_len.saturating_sub(window_size)) / step + 1);

    let mut counts = [0usize; 5];
    fill_counts(&mut counts, &seq[..window_size]);
    results.push(entropy_from_counts(&counts, window_size));

    if step >= window_size {
        let mut start = step;
        while start + window_size <= seq_len {
            fill_counts(&mut counts, &seq[start..start + window_size]);
            results.push(entropy_from_counts(&counts, window_size));
            start += step;
        }
        return results;
    }

    let mut start = 0;
    let mut end = window_size;
    while end + step <= seq_len {
        for &base in &seq[start..start + step] {
            counts[base_to_index(base)] -= 1;
        }
        for &base in &seq[end..end + step] {
            counts[base_to_index(base)] += 1;
        }
        start += step;
        end += step;
        results.push(entropy_from_counts(&counts, window_size));
    }

    results
}

fn entropy_from_counts(counts: &[usize; 5], window_len: usize) -> f64 {
    if window_len == 0 {
        return 0.0;
    }

    let mut entropy = 0.0;
    let total = window_len as f64;
    for &count in counts {
        if count != 0 {
            let p = count as f64 / total;
            entropy -= p * p.log2();
        }
    }
    entropy
}

// Linguistic complexity function (produces sliding window results)
fn linguistic_complexity(seq: &[u8], k: u8, w: usize) -> Vec<f64> {
    let n = seq.len();
    assert!(k <= 31); // Assuming reasonable k-mer size
    assert!(usize::from(k) > 0 && usize::from(k) < w && w <= n);

    // Compute k-mers
    let k = usize::from(k);
    let mut kmers: Vec<u64> = Vec::with_capacity(n - k + 1);
    let mut kmer = 0u64;
    for &base in &seq[..k] {
        kmer = (kmer << 2) | u64::from(base_to_code(base));
    }
    kmers.push(kmer);
    let mask = (1u64 << (2 * k)) - 1;
    for &base in &seq[k..] {
        kmer = ((kmer << 2) | u64::from(base_to_code(base))) & mask;
        kmers.push(kmer);
    }

    let theoretical_max = std::cmp::min(w + 1 - k, 1 << (2 * k));
    let mult = 1.0 / theoretical_max as f64;

    let mut counts: HashMap<u64, u16> = HashMap::new();
    let mut unique = 0;

    // Initialize first window
    for &kmer in &kmers[..(w - k + 1)] {
        let c = counts.entry(kmer).or_insert(0);
        if *c == 0 {
            unique += 1;
        }
        *c += 1;
    }

    let mut res = Vec::with_capacity(n - w + 1);
    res.push(unique as f64 * mult);

    // Slide window
    for i in 0..(kmers.len() - (w - k + 1)) {
        let lag_kmer = kmers[i];
        let new_kmer = kmers[i + (w - k + 1)];

        // Remove outgoing k-mer
        if let Some(count) = counts.get_mut(&lag_kmer) {
            *count -= 1;
            if *count == 0 {
                unique -= 1;
            }
        }

        // Add incoming k-mer
        let c = counts.entry(new_kmer).or_insert(0);
        if *c == 0 {
            unique += 1;
        }
        *c += 1;

        res.push(unique as f64 * mult);
    }
    res
}

#[derive(Clone)]
struct Node {
    sequence: String,
    length: usize,
}

#[derive(Clone)]
struct Path {
    name: String,
    nodes: Vec<(String, bool)>, // (node_id, is_forward)
}

struct GFA {
    nodes: HashMap<String, Node>,
    paths: Vec<Path>,
    edges: Vec<String>, // Store edge lines as-is
}

fn parse_gfa(filename: &FilePath) -> std::io::Result<GFA> {
    let mut file = File::open(filename)?;

    // Read magic bytes once to detect compression
    let mut magic_bytes = [0u8; 18];
    let mut read = 0usize;
    while read < magic_bytes.len() {
        let bytes = file.read(&mut magic_bytes[read..])?;
        if bytes == 0 {
            break;
        }
        read += bytes;
    }
    file.seek(SeekFrom::Start(0))?;

    // Select an appropriate reader (plain, gzip, or bgzip) without reopening the file.
    let reader: Box<dyn BufRead> = if read >= 2 && magic_bytes[0] == 0x1f && magic_bytes[1] == 0x8b
    {
        let is_bgzip = read >= 14
            && magic_bytes[2] == 0x08
            && magic_bytes[12] == b'B'
            && magic_bytes[13] == b'C';
        if is_bgzip {
            debug!("Using bgzip reader for bgzip compressed file");
            let bgzf_reader = bgzf::Reader::new(file);
            Box::new(BufReader::new(bgzf_reader))
        } else {
            debug!("Using standard gzip reader for gzip compressed file");
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        }
    } else {
        debug!("Using uncompressed file reader");
        Box::new(BufReader::new(file))
    };

    // Parse segments, paths, and edges eagerly so later stages can work in-memory.
    let mut nodes = HashMap::new();
    let mut paths = Vec::new();
    let mut edges = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        match parts[0] {
            "S" => {
                // Segment/Node line
                let id = parts[1].to_string();
                let sequence = parts[2].to_string();
                let length = sequence.len();
                nodes.insert(id, Node { sequence, length });
            }
            "P" => {
                // Path line
                let name = parts[1].to_string();
                let path_str = parts[2];
                let mut path_nodes = Vec::new();

                for node_str in path_str.split(',') {
                    let is_forward = node_str.ends_with('+');
                    let node_id = node_str[..node_str.len() - 1].to_string();
                    path_nodes.push((node_id, is_forward));
                }

                paths.push(Path {
                    name,
                    nodes: path_nodes,
                });
            }
            "L" => {
                // Edge line - store as-is
                edges.push(line);
            }
            _ => {}
        }
    }

    Ok(GFA {
        nodes,
        paths,
        edges,
    })
}

fn reconstruct_path_sequence(path: &Path, nodes: &HashMap<String, Node>) -> Vec<u8> {
    let mut sequence = Vec::new();

    for (node_id, is_forward) in &path.nodes {
        if let Some(node) = nodes.get(node_id) {
            let node_seq = node.sequence.as_bytes();
            if *is_forward {
                sequence.extend_from_slice(node_seq);
            } else {
                // Reverse complement
                for &base in node_seq.iter().rev() {
                    let rc = match base {
                        b'A' | b'a' => b'T',
                        b'T' | b't' => b'A',
                        b'C' | b'c' => b'G',
                        b'G' | b'g' => b'C',
                        _ => base,
                    };
                    sequence.push(rc);
                }
            }
        }
    }

    sequence
}

fn find_low_complexity_regions(
    complexity: &[f64],
    threshold: f64,
    window_size: usize,
    step_size: usize,
    merge_threshold: usize,
) -> (Vec<(usize, usize)>, Vec<(usize, usize, f64)>) {
    // Single pass merge: expand active region while scores stay below threshold, tracking mean on the fly.
    let mut merged_regions = Vec::new();
    let mut merged_windows = Vec::new();

    let mut current_start = 0usize;
    let mut current_end = 0usize;
    let mut sum = 0.0;
    let mut count = 0usize;
    let mut active = false;

    for (i, &score) in complexity.iter().enumerate() {
        if score >= threshold {
            continue;
        }

        let start = if step_size == window_size {
            i
        } else {
            i * step_size
        };
        let end = start + window_size;

        if active && start <= current_end + merge_threshold {
            current_end = current_end.max(end);
            sum += score;
            count += 1;
        } else {
            if active {
                merged_regions.push((current_start, current_end));
                merged_windows.push((current_start, current_end, sum / count as f64));
            }
            active = true;
            current_start = start;
            current_end = end;
            sum = score;
            count = 1;
        }
    }

    if !active {
        return (Vec::new(), Vec::new());
    }

    merged_regions.push((current_start, current_end));
    merged_windows.push((current_start, current_end, sum / count as f64));

    (merged_regions, merged_windows)
}

fn map_regions_to_nodes(
    path: &Path,
    nodes: &HashMap<String, Node>,
    regions: &[(usize, usize)],
) -> HashSet<String> {
    let mut marked_nodes = HashSet::new();
    let mut current_pos = 0;

    for (node_id, _is_forward) in &path.nodes {
        if let Some(node) = nodes.get(node_id) {
            let node_start = current_pos;
            let node_end = current_pos + node.length;

            // Check if this node overlaps with any low-complexity region
            for &(region_start, region_end) in regions {
                if !(region_end < node_start || region_start >= node_end) {
                    marked_nodes.insert(node_id.clone());
                    break;
                }
            }

            current_pos = node_end;
        }
    }

    marked_nodes
}

fn map_complexity_to_nodes(
    path: &Path,
    nodes: &HashMap<String, Node>,
    complexity_values: &[f64],
    window_size: usize,
    step_size: usize,
) -> HashMap<String, Vec<f64>> {
    let mut node_complexities: HashMap<String, Vec<f64>> = HashMap::new();
    let mut positions = Vec::with_capacity(path.nodes.len());
    let mut node_refs = Vec::with_capacity(path.nodes.len());
    let mut current_pos = 0usize;

    for (node_id, _is_forward) in &path.nodes {
        if let Some(node) = nodes.get(node_id) {
            let start = current_pos;
            let end = current_pos + node.length;
            positions.push((start, end));
            node_refs.push(node_id);
            current_pos = end;
        }
    }

    if positions.is_empty() || complexity_values.is_empty() {
        return node_complexities;
    }

    let mut values_per_node: Vec<Vec<f64>> = vec![Vec::new(); positions.len()];
    let mut leading_idx = 0usize;

    for (idx, &value) in complexity_values.iter().enumerate() {
        let window_start = if step_size == window_size {
            idx
        } else {
            idx * step_size
        };
        let window_end = window_start + window_size;

        while leading_idx < positions.len() && positions[leading_idx].1 <= window_start {
            leading_idx += 1;
        }

        let mut node_idx = leading_idx;
        while node_idx < positions.len() {
            let (node_start, node_end) = positions[node_idx];
            if node_start >= window_end {
                break;
            }
            if node_end > window_start {
                values_per_node[node_idx].push(value);
            }
            node_idx += 1;
        }
    }

    for (node_id, values) in node_refs.into_iter().zip(values_per_node.into_iter()) {
        if !values.is_empty() {
            node_complexities
                .entry(node_id.clone())
                .or_insert_with(Vec::new)
                .extend(values);
        }
    }

    node_complexities
}

fn annotate_gfa(
    gfa: &GFA,
    marked_nodes: &HashSet<String>,
    path_regions: &HashMap<String, Vec<(usize, usize)>>,
    output_file: &FilePath,
) -> std::io::Result<()> {
    let file = File::create(output_file)?;

    // Check if output should be compressed based on file extension
    let should_compress = output_file
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz" || ext == "bgz")
        .unwrap_or(false);

    let mut writer: Box<dyn Write> = if should_compress {
        Box::new(bgzf::Writer::new(file))
    } else {
        Box::new(file)
    };

    // Write nodes with annotations
    for (id, node) in &gfa.nodes {
        write!(writer, "S\t{}\t{}", id, node.sequence)?;
        if marked_nodes.contains(id) {
            write!(writer, "\tLC:i:1\tCL:z:red")?;
        }
        writeln!(writer)?;
    }

    // Write edges
    for edge in &gfa.edges {
        writeln!(writer, "{}", edge)?;
    }

    // Write paths with region annotations in comments
    for path in &gfa.paths {
        write!(writer, "P\t{}\t", path.name)?;
        let path_str: Vec<String> = path
            .nodes
            .iter()
            .map(|(id, forward)| format!("{}{}", id, if *forward { '+' } else { '-' }))
            .collect();
        write!(writer, "{}", path_str.join(","))?;

        // Add low-complexity regions as optional tags
        if let Some(regions) = path_regions.get(&path.name) {
            if !regions.is_empty() {
                write!(writer, "\tLR:Z:")?;
                let region_strs: Vec<String> = regions
                    .iter()
                    .map(|(s, e)| format!("{}-{}", s, e))
                    .collect();
                write!(writer, "{}", region_strs.join(","))?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input GFA file
    #[arg(short = 'i', long = "input-gfa")]
    input_gfa: String,

    /// Window size for complexity calculation
    #[arg(short = 'w', long = "window-size")]
    window_size: usize,

    /// Threshold for low-complexity regions (number or "auto")
    #[arg(short = 't', long = "threshold")]
    threshold: String,

    /// IQR multiplier for automatic threshold (default: 1.5, use with threshold "auto")
    #[arg(long = "iqr-multiplier", default_value = "1.5")]
    iqr_multiplier: f64,

    /// Complexity measure type: "linguistic" or "entropy" (default: "linguistic")
    #[arg(long, default_value = "linguistic")]
    complexity: String,

    /// K-mer size (used with linguistic complexity)
    #[arg(
        short = 'k',
        long = "k-mer",
        default_value_t = 16,
        value_parser = clap::value_parser!(u8).range(1..=31),
        conflicts_with = "step_size"
    )]
    k: u8,

    /// Step size for sliding window (used with entropy complexity)
    #[arg(
        short = 's',
        long = "step-size",
        default_value = "50",
        conflicts_with = "k"
    )]
    step_size: usize,

    /// Distance threshold for merging close ranges
    #[arg(short = 'd', long = "distance", default_value = "100")]
    merge_threshold: usize,

    /// Output annotated GFA file
    #[arg(short = 'o', long = "output-gfa")]
    output_gfa: Option<String>,

    /// Output BED file with low-complexity ranges
    #[arg(short = 'b', long = "bed")]
    bed: Option<String>,

    /// Output CSV file for Bandage node coloring (Node,Colour format)
    #[arg(short = 'c', long = "csv")]
    csv: Option<String>,

    /// Output boolean mask file: 1 if node is not annotated, 0 if annotated
    #[arg(short = 'm', long = "mask")]
    mask: Option<String>,

    /// Output weights file: node_id and its associated complexity/entropy weight
    #[arg(long = "weights")]
    weights: Option<String>,

    /// Verbosity level (0: errors only, 1: errors and info, 2: debug, 3: trace)
    #[arg(short = 'v', long = "verbose", default_value = "1")]
    verbose: u8,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    // Initialize logger with verbosity level
    let log_level = match args.verbose {
        0 => "error",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    // Check that at least one output format is specified
    if args.output_gfa.is_none()
        && args.bed.is_none()
        && args.csv.is_none()
        && args.mask.is_none()
        && args.weights.is_none()
    {
        error!(
            "At least one output format must be specified (--output-gfa, --bed, --csv, --mask, or --weights)"
        );
        std::process::exit(1);
    }

    // Validate complexity type
    if args.complexity != "linguistic" && args.complexity != "entropy" {
        error!("complexity must be either 'linguistic' or 'entropy'");
        std::process::exit(1);
    }

    // Validate threshold argument
    let is_auto_threshold = args.threshold == "auto";
    if !is_auto_threshold {
        if let Err(_) = args.threshold.parse::<f64>() {
            error!("Threshold must be a number or 'auto'");
            std::process::exit(1);
        }
    }

    let input_file = FilePath::new(&args.input_gfa);
    let window_size = args.window_size;

    info!("Parsing GFA file...");
    let gfa = parse_gfa(input_file)?;

    let mut all_marked_nodes = HashSet::new();
    let mut path_regions = HashMap::new();
    let mut path_windows = HashMap::new();
    let mut all_node_weights: HashMap<String, Vec<f64>> = HashMap::new();

    // Two-pass processing
    let mut all_complexity_values = Vec::new();
    let mut path_complexity_data: HashMap<String, Vec<f64>> = HashMap::new();

    info!("Analyzing {} paths...", gfa.paths.len());

    // First pass: compute complexity for all paths
    for path in &gfa.paths {
        debug!("Processing path: {}", path.name);

        // Reconstruct path sequence
        let sequence = reconstruct_path_sequence(path, &gfa.nodes);

        if sequence.len() <= window_size {
            warn!("  Path too short for window size, skipping");
            continue;
        }

        // Compute complexity based on selected method
        let complexity = match args.complexity.as_str() {
            "linguistic" => linguistic_complexity(&sequence, args.k, window_size),
            "entropy" => shannon_entropy_complexity(&sequence, window_size, args.step_size),
            _ => unreachable!(), // Already validated above
        };

        // Store complexity values for this path and collect all values
        path_complexity_data.insert(path.name.clone(), complexity.clone());
        all_complexity_values.extend(complexity.iter().copied());
    }

    // Determine threshold (either automatic or manual)
    let threshold = if is_auto_threshold {
        // Calculate automatic threshold
        if all_complexity_values.is_empty() {
            error!("No complexity values computed, cannot determine automatic threshold");
            std::process::exit(1);
        }

        let (q1, q3, iqr, auto_threshold) =
            calculate_iqr_threshold(&all_complexity_values, args.iqr_multiplier);

        info!("Automatic threshold calculation:");
        info!(
            "  Q1, Q3, IQR (Q3 - Q1), IQR multiplier: {:.4}, {:.4}, {:.4}, {:.4}",
            q1, q3, iqr, args.iqr_multiplier
        );
        info!(
            "  Using threshold: {:.4} ({:.4} - {:.4} * {:.4})",
            auto_threshold, q1, args.iqr_multiplier, iqr
        );

        auto_threshold
    } else {
        let manual_threshold = args.threshold.parse::<f64>().unwrap();
        info!("Using manual threshold: {:.4}", manual_threshold);
        manual_threshold
    };

    info!("Identifying low-complexity regions...");

    // Second pass: identify low-complexity regions and collect node weights
    for path in &gfa.paths {
        if let Some(complexity) = path_complexity_data.get(&path.name) {
            // For linguistic complexity, step_size is same as window_size (no overlap)
            // For entropy complexity, use provided step_size
            let step_size = if args.complexity == "linguistic" {
                window_size
            } else {
                args.step_size
            };

            // Collect node weights for this path
            let node_complexities =
                map_complexity_to_nodes(path, &gfa.nodes, complexity, window_size, step_size);

            // Add to global node weights collection
            for (node_id, complexities) in node_complexities {
                all_node_weights
                    .entry(node_id)
                    .or_insert_with(Vec::new)
                    .extend(complexities);
            }

            let (regions, windows) = find_low_complexity_regions(
                complexity,
                threshold,
                window_size,
                step_size,
                args.merge_threshold,
            );

            if !regions.is_empty() {
                debug!(
                    "  Path {}: Found {} low-complexity regions",
                    path.name,
                    regions.len()
                );

                // Map to nodes
                let marked = map_regions_to_nodes(path, &gfa.nodes, &regions);
                all_marked_nodes.extend(marked);

                path_regions.insert(path.name.clone(), regions);
                path_windows.insert(path.name.clone(), windows);
            }
        }
    }

    info!("Marked {} nodes as low-complexity", all_marked_nodes.len());

    if let Some(output_gfa_path) = &args.output_gfa {
        info!("Writing annotated GFA...");
        let output_file = FilePath::new(output_gfa_path);
        annotate_gfa(&gfa, &all_marked_nodes, &path_regions, output_file)?;
        info!("Annotated GFA written to: {}", output_gfa_path);
    }

    if let Some(bed_file) = &args.bed {
        info!("Writing BED file with low-complexity regions...");
        write_bed_file(&path_windows, bed_file, &gfa.paths)?;
        info!("BED file written to: {}", bed_file);
    }

    if let Some(csv_file) = &args.csv {
        info!("Writing CSV file for Bandage node coloring...");
        write_csv_file(&all_marked_nodes, csv_file)?;
        info!("CSV file written to: {}", csv_file);
    }

    if let Some(mask_file) = &args.mask {
        info!("Writing boolean mask file...");
        write_mask_file(&gfa.nodes, &all_marked_nodes, mask_file)?;
        info!("Mask file written to: {}", mask_file);
    }

    if let Some(weights_file) = &args.weights {
        info!("Writing weights file...");
        write_weights_file(&gfa.nodes, &all_node_weights, weights_file)?;
        info!("Weights file written to: {}", weights_file);
    }

    info!("Done!");
    Ok(())
}

fn write_bed_file(
    path_windows: &HashMap<String, Vec<(usize, usize, f64)>>,
    bed_file: &str,
    paths: &[Path],
) -> std::io::Result<()> {
    let mut file = File::create(bed_file)?;

    for path in paths {
        if let Some(windows) = path_windows.get(&path.name) {
            for &(start, end, complexity) in windows.iter() {
                // Unified format: chrom start end complexity score strand
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t0\t+",
                    path.name, start, end, complexity
                )?;
            }
        }
    }

    Ok(())
}

fn write_csv_file(marked_nodes: &HashSet<String>, csv_file: &str) -> std::io::Result<()> {
    let mut file = File::create(csv_file)?;

    // Write CSV header
    writeln!(file, "Node,Colour")?;

    // Write low-complexity nodes in red
    for node_id in marked_nodes {
        writeln!(file, "{},red", node_id)?;
    }

    Ok(())
}

fn write_mask_file(
    nodes: &HashMap<String, Node>,
    marked_nodes: &HashSet<String>,
    mask_file: &str,
) -> std::io::Result<()> {
    let mut file = File::create(mask_file)?;

    // Get sorted node IDs (assuming numeric IDs for proper ordering)
    let mut node_ids: Vec<&String> = nodes.keys().collect();
    node_ids.sort_by(|a, b| {
        // Try to parse as numbers first, fall back to string comparison
        match (a.parse::<i32>(), b.parse::<i32>()) {
            (Ok(a_num), Ok(b_num)) => a_num.cmp(&b_num),
            _ => a.cmp(b),
        }
    });

    // Write mask: 1 if NOT annotated (not low-complexity), 0 if annotated (low-complexity)
    for node_id in node_ids {
        let mask_value = if marked_nodes.contains(node_id) { 0 } else { 1 };
        writeln!(file, "{}", mask_value)?;
    }

    Ok(())
}

fn write_weights_file(
    nodes: &HashMap<String, Node>,
    node_weights: &HashMap<String, Vec<f64>>,
    weights_file: &str,
) -> std::io::Result<()> {
    let mut file = File::create(weights_file)?;

    // Get sorted node IDs (assuming numeric IDs for proper ordering)
    let mut node_ids: Vec<&String> = nodes.keys().collect();
    node_ids.sort_by(|a, b| {
        // Try to parse as numbers first, fall back to string comparison
        match (a.parse::<i32>(), b.parse::<i32>()) {
            (Ok(a_num), Ok(b_num)) => a_num.cmp(&b_num),
            _ => a.cmp(b),
        }
    });

    // Write node weights (one per line, row N for node ID N)
    for node_id in node_ids {
        let weight = if let Some(complexities) = node_weights.get(node_id) {
            if complexities.is_empty() {
                0.0
            } else {
                // Calculate mean complexity/entropy as the weight
                complexities.iter().sum::<f64>() / complexities.len() as f64
            }
        } else {
            // No complexity data for this node (e.g., not covered by any path)
            0.0
        };

        writeln!(file, "{}", weight)?;
    }

    Ok(())
}
