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
- GUI: detects and eliminates restriction enzyme recognition sites (any of the ~620
  commercially available enzymes) via the same synonymous-codon substitution approach,
  as a separate step after the main optimization

## Usage

### Option A — Run from source (Python, any OS including Linux/macOS)

Edit the `CODON_TABLE_TEXT`, `AMINO_ACID_SEQUENCE` (or `INPUT_NT_SEQUENCE`), and
`FREQ_THRESHOLD` constants near the bottom of `reverse_translation.py`, then run:

```bash
python reverse_translation.py
```

Requirements: Python 3.10+, and optionally `rich` for colored terminal output
(`pip install rich`). No GUI dependencies needed for this path.

### Option B — Windows GUI app (no Python required)

Download the latest `CodonOptimizer.exe` from the
[Releases page](https://github.com/gabpellehuet/codon_optimizer/releases) and run it
directly — no installation, no Python needed.

If you'd rather run the GUI from source instead:

```bash
pip install -r requirements.txt
python gui.py
```

The GUI itself only needs `tkinter` (ships with the standard Python installer on
Windows/macOS; on Linux it may need a separate package, e.g. `sudo apt install
python3-tk`), but the restriction-site tab depends on `biopython` — see
`requirements.txt`. The GUI is intended for Windows users; Linux users will
typically prefer Option A, which has no dependencies at all.

The GUI includes an "Organism" dropdown with codon usage tables bundled locally
(no internet connection needed) for *E. coli* K12, *S. cerevisiae*, *H. sapiens*,
*M. musculus*, *P. pastoris*, and *B. subtilis*, sourced from the
[Kazusa Codon Usage Database](https://www.kazusa.or.jp/codon/). You can also load
a custom table from a file, or use "Fetch from Kazusa online..." to search and
download a table for any organism in the live database directly (requires an
internet connection).

### Restriction site detection/elimination (GUI)

The "Restriction Sites" tab lets you select one or more restriction enzymes
(alphabetically listed, with a filter box) and:
- **Show restriction sites**: lists where each selected enzyme's recognition
  site occurs in the current optimized sequence (uses Biopython's
  `Bio.Restriction`).
- **Eliminate restriction sites**: removes them via synonymous codon
  substitution — same frequency-threshold logic as the main optimizer, with a
  fallback override (and matching warning) if no synonym above the threshold
  can remove a given site. Runs as a separate step on top of the main
  Optimize result, and its changes/warnings are added to the existing
  Changes/Warnings tabs.

## Building the Windows executable yourself

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name CodonOptimizer gui.py
```

The executable is written to `dist/CodonOptimizer.exe`.

## License

MIT
