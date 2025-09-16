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
# Manual threshold
panplexity -i input.gfa -t 0.9 [OPTIONS]

# Automatic threshold
panplexity -i input.gfa -t auto [OPTIONS]
```

### Key options
- `-i/--input-gfa`: GFA file (`.gfa`, `.gfa.gz`, `.gfa.bgz`)
- `-w/--window-size`: Window size for complexity calculation
- `-t/--threshold`: Complexity threshold (number or "auto")
- `--iqr-multiplier`: IQR multiplier for auto-threshold (default: 1.5)
- `--complexity`: "linguistic" (default) or "entropy"

### Output formats (choose one or more)

- `-o/--output-gfa`: Annotated GFA
- `-b/--bed`: BED file with regions
- `-c/--csv`: Bandage coloring file
- `-m/--mask`: Node boolean mask
- `--weights`: Node complexity weights

## Examples

```bash
# Linguistic complexity with BED output and automatic threshold
panplexity -i input.gfa -w 100 -t auto -b regions.bed

# Shannon entropy with annotated GFA output and stricter threshold
panplexity -i input.gfa -w 100 --complexity entropy -t auto --iqr-multiplier 3.0 -o output.gfa

# Multiple output formats
panplexity -i input.gfa -w 100 -t 0.9 -b regions.bed -c bandage.csv --weights weights.txt
```

## Output Formats

- **GFA**: Original GFA with `LC:i:1` and `CL:z:red` tags on low-complexity nodes
- **BED**: `chrom start end complexity_score 0 +`
- **CSV**: `Node,Colour` format for Bandage visualization
- **Mask**: One value per node (0=low-complexity, 1=normal)
- **Weights**: One weight per node
