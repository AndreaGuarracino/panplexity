use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, Read};
use std::path::Path as FilePath;
use clap::Parser;
use noodles::bgzf;

// Shannon entropy complexity function
pub fn shannon_entropy_complexity(seq: &[u8], window_size: usize, step: usize) -> Vec<f64> {
    if seq.is_empty() || window_size == 0 {
        return Vec::new();
    }
    
    let mut results = Vec::new();
    let seq_len = seq.len();
    
    for start in (0..=seq_len.saturating_sub(window_size)).step_by(step) {
        let end = (start + window_size).min(seq_len);
        let window_seq = &seq[start..end];
        
        // Calculate Shannon entropy for this window
        let entropy = shannon_entropy(window_seq);
        results.push(entropy);
    }
    
    results
}

fn shannon_entropy(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    
    let mut freqs = HashMap::new();
    for &base in seq {
        let normalized_base = match base {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            _ => b'N', // Handle N's and other chars
        };
        *freqs.entry(normalized_base).or_insert(0) += 1;
    }
    
    let seq_len = seq.len() as f64;
    let mut entropy = 0.0;
    
    for &count in freqs.values() {
        if count > 0 {
            let p = count as f64 / seq_len;
            entropy -= p * p.log2();
        }
    }
    
    entropy
}

// Linguistic complexity function (with minor fixes for compilation)
pub fn linguistic_complexity(seq: &[u8], k: u8, w: usize) -> Vec<f64> {
    let n = seq.len();
    assert!(k <= 16); // Assuming reasonable k-mer size
    assert!(usize::from(k) < w && w < n);

    // Compute k-mers (simplified version - you'd use your actual k-mer library)
    let mut kmers: Vec<u64> = Vec::with_capacity(n - usize::from(k) + 1);
    for i in 0..=(n - usize::from(k)) {
        let mut kmer = 0u64;
        for j in 0..usize::from(k) {
            let base = match seq[i + j] {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0, // Handle N's and other chars
            };
            kmer = (kmer << 2) | base;
        }
        kmers.push(kmer);
    }

    let k = usize::from(k);
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
    let file = File::open(filename)?;
    
    // Check if file is gzip/bgzip compressed by reading magic bytes
    let mut magic_bytes = [0u8; 3];
    let mut file_clone = File::open(filename)?;
    file_clone.read_exact(&mut magic_bytes)?;
    
    let reader: Box<dyn BufRead> = if magic_bytes == [0x1f, 0x8b, 0x08] {
        // bgzip/gzip compressed
        let bgzf_reader = bgzf::Reader::new(file);
        Box::new(BufReader::new(bgzf_reader))
    } else {
        // Uncompressed
        Box::new(BufReader::new(file))
    };
    
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
                    let node_id = node_str[..node_str.len()-1].to_string();
                    path_nodes.push((node_id, is_forward));
                }
                
                paths.push(Path { name, nodes: path_nodes });
            }
            "L" => {
                // Edge line - store as-is
                edges.push(line);
            }
            _ => {}
        }
    }
    
    Ok(GFA { nodes, paths, edges })
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
    with_entropy: bool,
) -> (Vec<(usize, usize)>, Vec<(usize, usize, f64)>) {
    let mut windows_with_entropy = Vec::new();
    
    // Collect all windows below threshold with their complexity values
    for (i, &score) in complexity.iter().enumerate() {
        if score < threshold {
            // For entropy complexity with overlapping windows, use step_size scaling
            // For linguistic complexity, each index represents sequence position
            let start = if step_size == window_size {
                // Non-overlapping windows (linguistic) - direct sequence position
                i
            } else {
                // Overlapping windows (entropy) - use step_size scaling
                i * step_size
            };
            let end = start + window_size;
            windows_with_entropy.push((start, end, score));
        }
    }
    
    if windows_with_entropy.is_empty() {
        return (Vec::new(), Vec::new());
    }
    
    // Sort by start position
    windows_with_entropy.sort_by_key(|w| w.0);
    
    let mut merged_regions = Vec::new();
    let mut merged_windows = Vec::new();
    let mut current_start = windows_with_entropy[0].0;
    let mut current_end = windows_with_entropy[0].1;
    let mut entropies = vec![windows_with_entropy[0].2];
    
    for &(start, end, entropy) in windows_with_entropy.iter().skip(1) {
        // Merge if overlapping or within merge_threshold distance
        if start <= current_end + merge_threshold {
            current_end = current_end.max(end);
            entropies.push(entropy);
        } else {
            // Finalize current merged region
            merged_regions.push((current_start, current_end));
            
            if with_entropy {
                let mean_entropy = entropies.iter().sum::<f64>() / entropies.len() as f64;
                merged_windows.push((current_start, current_end, mean_entropy));
            }
            
            // Start new region
            current_start = start;
            current_end = end;
            entropies = vec![entropy];
        }
    }
    
    // Add final region
    merged_regions.push((current_start, current_end));
    if with_entropy {
        let mean_entropy = entropies.iter().sum::<f64>() / entropies.len() as f64;
        merged_windows.push((current_start, current_end, mean_entropy));
    }
    
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

fn annotate_gfa(
    gfa: &GFA,
    marked_nodes: &HashSet<String>,
    path_regions: &HashMap<String, Vec<(usize, usize)>>,
    output_file: &FilePath,
) -> std::io::Result<()> {
    let file = File::create(output_file)?;
    
    // Check if output should be compressed based on file extension
    let should_compress = output_file.extension()
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
        let path_str: Vec<String> = path.nodes
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
    #[arg(short = 'g', long = "gfa")]
    input_gfa: String,
    
    /// Output annotated GFA file
    #[arg(short = 'o', long = "output")]
    output_gfa: Option<String>,
    
    /// K-mer size (used with linguistic complexity)
    #[arg(short, long)]
    k: Option<u8>,
    
    /// Window size for complexity calculation
    #[arg(short, long)]
    window_size: usize,
    
    /// Threshold for low-complexity regions
    #[arg(short, long)]
    threshold: f64,
    
    /// Complexity measure type: "linguistic" or "entropy" (default: "linguistic")
    #[arg(long, default_value = "linguistic")]
    complexity_type: String,
    
    /// Step size for sliding window (used with entropy complexity, default: window_size/2)
    #[arg(long)]
    step_size: Option<usize>,
    
    /// BED file to emit low-complexity ranges with info
    #[arg(short, long)]
    bed: Option<String>,
    
    /// CSV file for Bandage node coloring (Node,Colour format)
    #[arg(short = 'c', long = "csv")]
    csv: Option<String>,
    
    /// Distance threshold for merging close ranges in base pairs (default: 100)
    #[arg(short = 'd', long = "distance", default_value = "100")]
    merge_threshold: usize,
    
    /// Boolean mask file: 1 if node is not annotated, 0 if annotated
    #[arg(short = 'm', long = "mask")]
    mask: Option<String>,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    
    // Check that at least one output format is specified
    if args.output_gfa.is_none() && args.bed.is_none() && args.csv.is_none() && args.mask.is_none() {
        eprintln!("Error: At least one output format must be specified (--output, --bed, --csv, or --mask)");
        std::process::exit(1);
    }
    
    // Validate complexity type
    if args.complexity_type != "linguistic" && args.complexity_type != "entropy" {
        eprintln!("Error: complexity_type must be either 'linguistic' or 'entropy'");
        std::process::exit(1);
    }
    
    // Validate k parameter for linguistic complexity
    if args.complexity_type == "linguistic" && args.k.is_none() {
        eprintln!("Error: k parameter is required when using linguistic complexity");
        std::process::exit(1);
    }
    
    // Validate step_size parameter - only for entropy complexity
    if args.complexity_type == "linguistic" && args.step_size.is_some() {
        eprintln!("Error: step_size parameter cannot be used with linguistic complexity");
        std::process::exit(1);
    }
    
    // Validate k parameter - only for linguistic complexity
    if args.complexity_type == "entropy" && args.k.is_some() {
        eprintln!("Error: k parameter cannot be used with entropy complexity");
        std::process::exit(1);
    }
    
    let input_file = FilePath::new(&args.input_gfa);
    let window_size = args.window_size;
    let threshold = args.threshold;
    // For linguistic complexity, step_size is same as window_size (no overlap)
    // For entropy complexity, use provided step_size or default to window_size/2
    let step_size = if args.complexity_type == "linguistic" {
        window_size
    } else {
        args.step_size.unwrap_or(window_size / 2)
    };
    
    println!("Parsing GFA file...");
    let gfa = parse_gfa(input_file)?;
    
    let mut all_marked_nodes = HashSet::new();
    let mut path_regions = HashMap::new();
    let mut path_windows = HashMap::new();
    
    println!("Analyzing {} paths...", gfa.paths.len());
    for path in &gfa.paths {
        println!("Processing path: {}", path.name);
        
        // Reconstruct path sequence
        let sequence = reconstruct_path_sequence(path, &gfa.nodes);
        
        if sequence.len() <= window_size {
            println!("  Path too short for window size, skipping");
            continue;
        }
        
        // Compute complexity based on selected method
        let complexity = match args.complexity_type.as_str() {
            "linguistic" => {
                let k = args.k.unwrap(); // Already validated above
                linguistic_complexity(&sequence, k, window_size)
            },
            "entropy" => {
                shannon_entropy_complexity(&sequence, window_size, step_size)
            },
            _ => unreachable!(), // Already validated above
        };
        
        // Find low-complexity regions with optional complexity averaging
        let with_scores = true; // Always compute scores for unified output
        let (regions, windows) = find_low_complexity_regions(
            &complexity, 
            threshold, 
            window_size, 
            step_size,
            args.merge_threshold,
            with_scores
        );
        
        if !regions.is_empty() {
            println!("  Found {} low-complexity regions", regions.len());
            
            // Map to nodes
            let marked = map_regions_to_nodes(path, &gfa.nodes, &regions);
            all_marked_nodes.extend(marked);
            
            path_regions.insert(path.name.clone(), regions);
            path_windows.insert(path.name.clone(), windows);
        }
    }
    
    println!("Marked {} nodes as low-complexity", all_marked_nodes.len());
    
    if let Some(output_gfa_path) = &args.output_gfa {
        println!("Writing annotated GFA...");
        let output_file = FilePath::new(output_gfa_path);
        annotate_gfa(&gfa, &all_marked_nodes, &path_regions, output_file)?;
        println!("Annotated GFA written to: {}", output_gfa_path);
    }
    
    if let Some(bed_file) = &args.bed {
        let with_entropy = args.complexity_type == "entropy";
        if with_entropy {
            println!("Writing BED file with low-complexity windows...");
        } else {
            println!("Writing BED file with low-complexity regions...");
        }
        write_bed_file(&path_regions, &path_windows, bed_file, &gfa.paths, with_entropy)?;
        println!("BED file written to: {}", bed_file);
    }
    
    if let Some(csv_file) = &args.csv {
        println!("Writing CSV file for Bandage node coloring...");
        write_csv_file(&all_marked_nodes, csv_file)?;
        println!("CSV file written to: {}", csv_file);
    }
    
    if let Some(mask_file) = &args.mask {
        println!("Writing boolean mask file...");
        write_mask_file(&gfa.nodes, &all_marked_nodes, mask_file)?;
        println!("Mask file written to: {}", mask_file);
    }
    
    println!("Done!");
    Ok(())
}

fn write_bed_file(
    path_regions: &HashMap<String, Vec<(usize, usize)>>,
    path_windows: &HashMap<String, Vec<(usize, usize, f64)>>,
    bed_file: &str,
    paths: &[Path],
    with_entropy: bool,
) -> std::io::Result<()> {
    let mut file = File::create(bed_file)?;
    
    for path in paths {
        if with_entropy {
            // Write with entropy values
            if let Some(windows) = path_windows.get(&path.name) {
                for &(start, end, entropy) in windows.iter() {
                    writeln!(
                        file,
                        "{}\t{}\t{}\t{:.4}",
                        path.name, start, end, entropy
                    )?;
                }
            }
        } else {
            // Write with generic region labels
            if let Some(regions) = path_regions.get(&path.name) {
                for (i, &(start, end)) in regions.iter().enumerate() {
                    writeln!(
                        file,
                        "{}\t{}\t{}\tlow_complexity_region_{}\t0\t+",
                        path.name, start, end, i + 1
                    )?;
                }
            }
        }
    }
    
    Ok(())
}

fn write_csv_file(
    marked_nodes: &HashSet<String>,
    csv_file: &str,
) -> std::io::Result<()> {
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