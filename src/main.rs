use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path as FilePath;
use clap::Parser;

// Your complexity function (with minor fixes for compilation)
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
    let reader = BufReader::new(file);
    
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
    merge_threshold: usize,
) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut in_region = false;
    let mut start = 0;
    
    for (i, &score) in complexity.iter().enumerate() {
        if score < threshold {
            if !in_region {
                start = i;
                in_region = true;
            }
        } else if in_region {
            // End of low-complexity region
            regions.push((start, i + window_size - 1));
            in_region = false;
        }
    }
    
    // Handle region extending to end
    if in_region {
        regions.push((start, complexity.len() + window_size - 1));
    }
    
    // Merge regions within merge_threshold distance
    if regions.is_empty() {
        return regions;
    }
    
    regions.sort_by_key(|r| r.0);
    let mut merged = vec![regions[0]];
    
    for region in regions.iter().skip(1) {
        let last = merged.last_mut().unwrap();
        // Merge if regions overlap or are within merge_threshold distance
        if region.0 <= last.1 + merge_threshold {
            last.1 = last.1.max(region.1);
        } else {
            merged.push(*region);
        }
    }
    
    merged
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
    let mut file = File::create(output_file)?;
    
    // Write nodes with annotations
    for (id, node) in &gfa.nodes {
        write!(file, "S\t{}\t{}", id, node.sequence)?;
        if marked_nodes.contains(id) {
            write!(file, "\tLC:i:1\tCL:z:red")?;
        }
        writeln!(file)?;
    }
    
    // Write edges
    for edge in &gfa.edges {
        writeln!(file, "{}", edge)?;
    }
    
    // Write paths with region annotations in comments
    for path in &gfa.paths {
        write!(file, "P\t{}\t", path.name)?;
        let path_str: Vec<String> = path.nodes
            .iter()
            .map(|(id, forward)| format!("{}{}", id, if *forward { '+' } else { '-' }))
            .collect();
        write!(file, "{}", path_str.join(","))?;
        
        // Add low-complexity regions as optional tags
        if let Some(regions) = path_regions.get(&path.name) {
            if !regions.is_empty() {
                write!(file, "\tLR:Z:")?;
                let region_strs: Vec<String> = regions
                    .iter()
                    .map(|(s, e)| format!("{}-{}", s, e))
                    .collect();
                write!(file, "{}", region_strs.join(","))?;
            }
        }
        writeln!(file)?;
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
    
    /// K-mer size
    #[arg(short, long)]
    k: u8,
    
    /// Window size for complexity calculation
    #[arg(short, long)]
    window_size: usize,
    
    /// Threshold for low-complexity regions
    #[arg(short, long)]
    threshold: f64,
    
    /// BED file to emit low-complexity ranges with info
    #[arg(long)]
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
    
    let input_file = FilePath::new(&args.input_gfa);
    let k = args.k;
    let window_size = args.window_size;
    let threshold = args.threshold;
    
    println!("Parsing GFA file...");
    let gfa = parse_gfa(input_file)?;
    
    let mut all_marked_nodes = HashSet::new();
    let mut path_regions = HashMap::new();
    
    println!("Analyzing {} paths...", gfa.paths.len());
    for path in &gfa.paths {
        println!("Processing path: {}", path.name);
        
        // Reconstruct path sequence
        let sequence = reconstruct_path_sequence(path, &gfa.nodes);
        
        if sequence.len() <= window_size {
            println!("  Path too short for window size, skipping");
            continue;
        }
        
        // Compute complexity
        let complexity = linguistic_complexity(&sequence, k, window_size);
        
        // Find low-complexity regions
        let regions = find_low_complexity_regions(&complexity, threshold, window_size, args.merge_threshold);
        
        if !regions.is_empty() {
            println!("  Found {} low-complexity regions", regions.len());
            
            // Map to nodes
            let marked = map_regions_to_nodes(path, &gfa.nodes, &regions);
            all_marked_nodes.extend(marked);
            
            path_regions.insert(path.name.clone(), regions);
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
        println!("Writing BED file with low-complexity regions...");
        write_bed_file(&path_regions, bed_file, &gfa.paths)?;
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
    bed_file: &str,
    paths: &[Path],
) -> std::io::Result<()> {
    let mut file = File::create(bed_file)?;
    
    for path in paths {
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