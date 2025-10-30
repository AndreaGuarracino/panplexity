# Panplexity

Find low-complexity regions in pangenome graphs using linguistic complexity or Shannon entropy.

## Installation

```bash
git clone https://github.com/AndreaGuarracino/panplexity.git
cd panplexity
cargo build --release
```

## Usage

```
panplexity -i input.gfa [OPTIONS]
```

### Key options
- `-i/--input-gfa`: GFA file (`.gfa`, `.gfa.gz`, `.gfa.bgz`)
- `-w/--window-size`: Window size for complexity calculation (default: 100)
- `-t/--threshold`: Complexity threshold: number or "auto" (default)
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
# Linguistic complexity with BED output (defaults to window size 100, auto threshold)
panplexity -i input.gfa -b regions.bed

# Shannon entropy with annotated GFA output and stricter auto multiplier
panplexity -i input.gfa --complexity entropy --iqr-multiplier 3.0 -o output.gfa

# Multiple output formats with manual threshold
panplexity -i input.gfa -t 0.9 -b regions.bed -c bandage.csv --weights weights.txt
```

## Output Formats

- **GFA**: Original GFA with `LC:i:1` and `CL:z:red` tags on low-complexity nodes
- **BED**: `chrom start end complexity_score 0 +`
- **CSV**: `Node,Colour` format for Bandage visualization
- **Mask**: One value per node (0=low-complexity, 1=normal)
- **Weights**: One weight per node
