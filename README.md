# Panplexity

Identify low-complexity regions in pangenome graphs using linguistic complexity or Shannon entropy.

## Installation

```bash
git clone https://github.com/AndreaGuarracino/panplexity.git
cd panplexity
cargo build --release
```

## Usage

```bash
./target/release/panplexity -g input.gfa -w 100 -t 0.5 [OPTIONS]
```

**Required:**
- `-g`: Input GFA file (supports `.gfa`, `.gfa.gz`, `.gfa.bgz`)
- `-w`: Window size  
- `-t`: Low-complexity threshold

**Complexity methods:**
- `--complexity-type linguistic` (default): K-mer diversity, requires `-k`
- `--complexity-type entropy`: Shannon entropy of bases, optional `--step-size`

**Outputs (at least one):**
- `-o`: Annotated GFA
- `--bed`: BED file with regions
- `-c`: CSV for Bandage coloring
- `-m`: Boolean mask

## Examples

```bash
# Linguistic complexity (default)
./target/release/panplexity -g graph.gfa -k 3 -w 50 -t 0.3 -o output.gfa

# Shannon entropy
./target/release/panplexity -g graph.gfa -w 50 -t 1.5 --complexity-type entropy -o output.gfa

# Multiple outputs
./target/release/panplexity -g input.gfa.gz -k 3 -w 50 -t 0.3 \
    -o output.gfa.gz --bed regions.bed -c colors.csv
```