# Panplexity

Find low-complexity regions in pangenome graphs using linguistic complexity or Shannon entropy.

## Installation

```bash
git clone https://github.com/AndreaGuarracino/panplexity.git
cd panplexity
cargo build --release
```

## Usage

```bash
./target/release/panplexity -i input.gfa -w 100 -t 0.5 [OPTIONS]
```

### Required
- `-i/--input-gfa`: GFA file (`.gfa`, `.gfa.gz`, `.gfa.bgz`)
- `-w/--window-size`: Window size
- `-t/--threshold`: Complexity threshold

### Parameters
- `--complexity`: Complexity type: "linguistic" or "entropy" (default: "linguistic")
- `-k/--k-mer`: K-mer size for linguistic complexity (default: 16)
- `-s/--step-size`: Step size for entropy sliding window (default: 50)
- `-d/--distance`: Distance for merging nearby regions (default: 100)


### Output options (at least one required)
- `-o/--output-gfa`: Annotated GFA with nodes marked with `LC:i:1` and `CL:z:red` tags
- `-b/--bed`: BED file with low-complexity ranges and scores
- `-c/--csv`: CSV file for Bandage node coloring (Node,Colour format)
- `-m/--mask`: Boolean mask file (0=low-complexity, 1=normal)

## Examples

```bash
# Linguistic complexity with BED output
./target/release/panplexity -i input.gfa -w 100 -t 0.5 -b output.bed

# Shannon entropy with annotated GFA output
./target/release/panplexity -i input.gfa -w 50 -t 1.5 --complexity entropy -o annotated.gfa

# Multiple output formats
./target/release/panplexity -i input.gfa -w 200 -t 0.3 -b regions.bed -c bandage.csv -m mask.txt
```

## Output Formats

- **BED**: `chrom start end complexity_score 0 +`
- **CSV**: `Node,Colour` format for Bandage visualization
- **Mask**: One value per node (0=low-complexity, 1=normal)
- **GFA**: Original GFA with `LC:i:1` and `CL:z:red` tags on low-complexity nodes