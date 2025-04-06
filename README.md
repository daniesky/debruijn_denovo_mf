# 🧬 De Bruijn Motif Discovery

This project implements a **De Bruijn graph-based algorithm** for **de novo motif discovery** in DNA sequences.

Results include **motif scores** and optional **motif logos** for visual inspection.

---

## Input

- A standard **FASTA** file containing DNA sequences.

---

## Usage Options

### Option 1: Run with Python (Recommended for full customization)

Use the `main.py` interface to access all parameters.

```bash
python src/main.py data/sequences.fa \
    --limit 10000 \
    --k 5 \
    --gaps True \
    --inst_limit 0.05 \
    --save_logos 10 \
    --print_motifs 10 \
    --overlap_factor 0.25 \
    --threshold 0.25 \
    --scoring_limit 10
```

### Option 2: Run via Makefile (For quick dev runs)

```bash
make setup           # Setup environment and install dependencies
make run             # Run motif discovery with default test parameters
make clean           # Clean environment and generated logos
make help            # List all Makefile targets
```

> ⚠️ `make` is intended for development/testing. Use direct Python calls for access to all parameters during experiments.

---

## 🔍 Parameters

| Argument         | Type    | Description                                                                 |
|------------------|---------|-----------------------------------------------------------------------------|
| `fasta_file`     | str     | Path to the input FASTA file                                                |
| `--limit`        | int     | Maximum number of sequences to process                                      |
| `--k`            | int     | Length of k-mers for De Bruijn graph                                        |
| `--gaps`         | bool    | Allow ambiguous (IUPAC) motif discovery                                     |
| `--inst_limit`   | float   | Threshold on motif instance frequency                                       |
| `--save_logos`   | int     | Number of top motifs to visualize as logos                                  |
| `--print_motifs` | int     | Number of top motifs to print                                               |
| `--overlap_factor`| float  | Minimum overlap ratio required during path traversal                        |
| `--threshold`    | float   | Weight threshold (as fraction of max) for edge traversal                    |
| `--scoring_limit`| int     | Limit for number of motifs to score (impacts runtime)                       |

---

## 📜 Requirements

Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
```