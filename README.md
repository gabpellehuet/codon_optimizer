# Codon Optimizer

A small Python tool to reverse-translate a protein sequence into DNA and optimize
the resulting codons for a given organism's codon usage table.

## Features

- Reverse-translates an amino acid sequence using a user-supplied codon usage table
- Or takes an existing nucleotide sequence and optimizes it directly
- Detects and repairs:
  - Homopolymer runs of 4+ identical nucleotides (e.g. `AAAA`, `GGGG`)
  - Identical adjacent codons
  - Identical codons separated by exactly one codon (skip-1 repeats)
- Upgrades codons to the most frequent synonym above a configurable frequency threshold
- Rich, color-coded terminal output (falls back to plain text if `rich` isn't installed)

## Usage

Edit the `CODON_TABLE_TEXT`, `AMINO_ACID_SEQUENCE` (or `INPUT_NT_SEQUENCE`), and
`FREQ_THRESHOLD` constants near the bottom of `reverse_translation.py`, then run:

```bash
python reverse_translation.py
```

## Requirements

- Python 3.10+
- `rich` (optional, for colored output): `pip install rich`

## License

MIT
