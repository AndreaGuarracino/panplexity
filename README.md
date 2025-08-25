# Panplexity

A tool for identifying and annotating low-complexity regions in pangenome graphs (GFA format).

## Overview

Panplexity analyzes sequences in GFA files to identify low-complexity regions using linguistic complexity measures based on k-mer diversity. It can output results in multiple formats for downstream analysis and visualization.

## Installation

```bash
git clone https://github.com/AndreaGuarracino/panplexity.git
cd panplexity
cargo build --release
```

## Usage

```bash
./target/release/panplexity -g input.gfa -k 3 -w 100 -t 0.5 [OPTIONS]
```

### Required Parameters

- `-g, --gfa <FILE>`: Input GFA file
- `-k, --k <SIZE>`: K-mer size for complexity calculation
- `-w, --window-size <SIZE>`: Window size for sliding window analysis
- `-t, --threshold <VALUE>`: Threshold for identifying low-complexity regions

### Output Options (at least one required)

- `-o, --output <FILE>`: Output annotated GFA file with LC tags
- `--bed <FILE>`: BED file with low-complexity region coordinates
- `-c, --csv <FILE>`: CSV file for Bandage node coloring (Node,Colour format)
- `-m, --mask <FILE>`: Boolean mask file (1=not annotated, 0=annotated)

### Optional Parameters

- `-d, --distance <BP>`: Distance threshold for merging close ranges (default: 100)

## Example

```bash
# Basic analysis with annotated GFA output
./target/release/panplexity -g graph.gfa -k 3 -w 50 -t 0.3 -o annotated.gfa

# Multiple output formats
./target/release/panplexity -g graph.gfa -k 3 -w 50 -t 0.3 \
    -o annotated.gfa --bed regions.bed -c colors.csv -m mask.txt

# Merge close regions within 200bp
./target/release/panplexity -g graph.gfa -k 3 -w 50 -t 0.3 \
    -o annotated.gfa -d 200
```

## Output Formats

- **Annotated GFA**: Original GFA with `LC:i:1` and `CL:z:red` tags on low-complexity nodes
- **BED**: Tab-separated format with region coordinates for each path
- **CSV**: Bandage-compatible node coloring file
- **Mask**: Boolean mask with one line per node (sorted by node ID)
