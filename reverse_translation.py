"""
Codon Optimizer for Plasmid Construction
=========================================
- Reverse-translates an amino acid sequence using a user-supplied codon table
- OR takes a nucleotide sequence directly and optimizes it
- Detects and repairs:
    1. Homopolymer runs of >= 4 identical nucleotides (AAAA, CCCC, GGGG, TTTT)
    2. Identical adjacent codons (position i == i+1), only when the fix does
       not introduce a homopolymer run
- Replacements use the next-most-frequent synonym above a frequency threshold
- Rich terminal output with colour-coded annotations (requires `rich`)
  Falls back to plain text if `rich` is not installed.

Usage
-----
  python codon_optimizer.py

Edit CODON_TABLE_TEXT, AMINO_ACID_SEQUENCE (or INPUT_NT_SEQUENCE), and
FREQ_THRESHOLD at the bottom of this file.
"""

from __future__ import annotations
import re
from collections import defaultdict

# ─────────────────────────────────────────────────────────────────────────────
# Codon-table parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_codon_table(text: str) -> dict[str, list[tuple[str, float]]]:
    """
    Parse a codon table in the format:
        AmAcid  Codon  Number  /1000  Fraction
        Gly     GGG    50527   11.12  0.15

    Returns {AA_three_letter: [(codon, fraction), ...]} sorted by fraction DESC.
    """
    table: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith(';'):
            continue
        if re.match(r'AmAcid|Amino', line, re.I):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        aa, codon, fraction = parts[0], parts[1].upper(), parts[4]
        if not re.match(r'^[ACGT]{3}$', codon):
            continue
        try:
            table[aa].append((codon, float(fraction)))
        except ValueError:
            continue
    for aa in table:
        table[aa].sort(key=lambda x: x[1], reverse=True)
    return dict(table)


# ─────────────────────────────────────────────────────────────────────────────
# Amino acid helpers
# ─────────────────────────────────────────────────────────────────────────────

ONE_TO_THREE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'End',
}


def _parse_aa_input(seq: str) -> list[str]:
    """Accept 1-letter or 3-letter AA sequences (space-separated or concatenated)."""
    seq = seq.strip()
    parts = re.split(r'[\s,]+', seq)
    if len(parts) > 1 and all(len(p) == 3 for p in parts if p):
        return [p.capitalize() for p in parts if p]
    if len(seq) % 3 == 0 and seq[:3].capitalize() in ONE_TO_THREE.values():
        return [seq[i:i+3].capitalize() for i in range(0, len(seq), 3)]
    result = []
    for ch in seq.upper():
        if ch in ONE_TO_THREE:
            result.append(ONE_TO_THREE[ch])
        elif ch.isspace():
            continue
        else:
            raise ValueError(f"Unknown amino acid character: '{ch}'")
    return result


def _aa_for_codon(codon: str, codon_table: dict) -> str:
    for aa, choices in codon_table.items():
        for c, _ in choices:
            if c == codon:
                return aa
    return '?'


def reverse_translate(aa_sequence: str,
                      codon_table: dict,
                      freq_threshold: float = 0.10) -> list[str]:
    """Reverse-translate a protein sequence using the most frequent valid codon."""
    aas = _parse_aa_input(aa_sequence)
    codons = []
    for aa in aas:
        choices = codon_table.get(aa, [])
        if not choices:
            raise ValueError(f"Amino acid '{aa}' not found in codon table.")
        above = [(c, f) for c, f in choices if f >= freq_threshold]
        codons.append((above[0][0] if above else choices[0][0]))
    return codons


# ─────────────────────────────────────────────────────────────────────────────
# Violation detection
# ─────────────────────────────────────────────────────────────────────────────

def find_nt_repeats(codons: list[str], min_len: int = 4) -> list[tuple[str, list[int]]]:
    """
    Find nucleotide substrings of length >= min_len appearing >= min_len times.
    Returns unique substrings (longer ones subsume shorter), sorted by length DESC.
    """
    seq = ''.join(codons)
    seen: dict[str, list[int]] = {}
    for length in range(min_len, len(seq) // 2 + 1):
        for i in range(len(seq) - length + 1):
            s = seq[i:i+length]
            seen.setdefault(s, []).append(i)
    results = [(s, pos) for s, pos in seen.items() if len(pos) >= min_len]
    results.sort(key=lambda x: len(x[0]), reverse=True)
    # keep only longest non-subsumed patterns
    kept, kept_strs = [], set()
    for s, pos in results:
        if any(s in longer for longer in kept_strs):
            continue
        kept.append((s, pos))
        kept_strs.add(s)
    return kept


def find_homopolymer_runs(codons: list[str], min_len: int = 4) -> list[tuple[str, int, int]]:
    """Return (nucleotide, start_pos, length) for each run of >= min_len identical nt."""
    seq = ''.join(codons)
    runs = []
    i = 0
    while i < len(seq):
        j = i + 1
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            runs.append((seq[i], i, j - i))
        i = j
    return runs


def find_identical_adjacent(codons: list[str]) -> list[tuple[int, int, str]]:
    """[(i, i+1, codon)] for identical adjacent codon pairs."""
    return [(i, i+1, codons[i])
            for i in range(len(codons) - 1)
            if codons[i] == codons[i+1]]


def find_identical_skip1(codons: list[str]) -> list[tuple[int, int, str]]:
    """[(i, i+2, codon)] for identical codons separated by exactly one codon."""
    return [(i, i+2, codons[i])
            for i in range(len(codons) - 2)
            if codons[i] == codons[i+2]]


# ─────────────────────────────────────────────────────────────────────────────
# Synonym lookup
# ─────────────────────────────────────────────────────────────────────────────

def get_synonyms(codon: str, codon_table: dict, freq_threshold: float = 0.0) -> list[str]:
    """
    Synonyms for codon excluding itself, sorted by fraction DESC.
    Pass freq_threshold=0 to get every synonym regardless of frequency.
    """
    for aa, choices in codon_table.items():
        for c, _ in choices:
            if c == codon:
                return [c2 for c2, f2 in choices if c2 != codon and f2 >= freq_threshold]
    return []


# ─────────────────────────────────────────────────────────────────────────────
# Core optimizer
# ─────────────────────────────────────────────────────────────────────────────

def _run_optimization(codons: list[str], codon_table: dict,
                      freq_threshold: float, max_iterations: int
                      ) -> tuple[list[str], list, list, dict]:
    """
    Multi-phase optimizer.

    Phases 1-3 (homopolymer runs, adjacent identical codons, skip-1 identical
    codons) first try to fix each violation using only synonyms at or above
    freq_threshold. If none of those can fix it, they fall back to the full
    synonym pool (ignoring the threshold) as a last resort, so a violation is
    never left unfixed just because every fixing codon happens to be rare.
    Such fallback substitutions are tagged "(threshold override)" in the
    change log. If no synonym at all (even ignoring frequency) can fix it,
    the violation is warned and skipped.

    Frequency upgrade
        Each remaining codon is upgraded to its most frequent synonym
        (>= freq_threshold) that keeps the sequence free of runs and
        adjacent/skip-1 identical pairs.
    """
    import random
    rng = random.Random(42)
    codons = list(codons)
    changes: list[tuple[int, str, str, str]] = []
    warnings: list[str] = []

    def _n_homo(seq_list: list[str]) -> int:
        return len(find_homopolymer_runs(seq_list))

    # ── Phase 1: Homopolymer runs ─────────────────────────────────────────
    stuck_starts: set[int] = set()
    for _ in range(max_iterations):
        runs = [(nt, s, l) for nt, s, l in find_homopolymer_runs(codons)
                if s not in stuck_starts]
        if not runs:
            break

        improved = False
        for nt, start, run_len in runs:
            end = start + run_len
            cand_idxs = sorted({p // 3 for p in range(start, end)
                                 if 0 <= p // 3 < len(codons)})
            cur_n = _n_homo(codons)

            fix_idx = fix_alt = fix_pair = None
            used_fallback = False
            for threshold_try, is_fallback in [(freq_threshold, False), (0.0, True)]:
                candidates = [(ci, alt)
                              for ci in cand_idxs
                              for alt in get_synonyms(codons[ci], codon_table, threshold_try)]
                rng.shuffle(candidates)

                best_idx, best_alt, best_n = None, None, cur_n
                for ci, alt in candidates:
                    trial = codons[:]
                    trial[ci] = alt
                    n = _n_homo(trial)
                    if n < best_n:
                        best_n, best_idx, best_alt = n, ci, alt

                if best_idx is not None:
                    fix_idx, fix_alt, used_fallback = best_idx, best_alt, is_fallback
                    break

                # Single change insufficient — try all pairs of substitutions
                best_pair, best_pair_n = None, cur_n
                for ci1, alt1 in candidates:
                    trial1 = codons[:]
                    trial1[ci1] = alt1
                    for ci2, alt2 in candidates:
                        if ci2 == ci1:
                            continue
                        trial2 = trial1[:]
                        trial2[ci2] = alt2
                        n = _n_homo(trial2)
                        if n < best_pair_n:
                            best_pair_n, best_pair = n, (ci1, alt1, ci2, alt2)

                if best_pair is not None:
                    fix_pair, used_fallback = best_pair, is_fallback
                    break

            suffix = " (threshold override)" if used_fallback else ""
            if fix_idx is not None:
                old = codons[fix_idx]
                codons[fix_idx] = fix_alt
                changes.append((fix_idx, old, fix_alt,
                                 f"homopolymer '{nt * run_len}' @{start}{suffix}"))
                if used_fallback:
                    warnings.append(
                        f"Homopolymer '{nt * run_len}' @{start}: no synonym >= "
                        f"freq_threshold {freq_threshold} could fix this; used "
                        f"{fix_alt} below threshold as a last resort."
                    )
                improved = True
                break
            elif fix_pair is not None:
                ci1, alt1, ci2, alt2 = fix_pair
                for ci, alt in [(ci1, alt1), (ci2, alt2)]:
                    old = codons[ci]
                    codons[ci] = alt
                    changes.append((ci, old, alt,
                                     f"homopolymer '{nt * run_len}' @{start} (pair){suffix}"))
                if used_fallback:
                    warnings.append(
                        f"Homopolymer '{nt * run_len}' @{start}: no synonym pair >= "
                        f"freq_threshold {freq_threshold} could fix this; used a pair "
                        f"below threshold as a last resort."
                    )
                improved = True
                break
            else:
                stuck_starts.add(start)
                warnings.append(
                    f"Homopolymer '{nt * run_len}' @{start}: "
                    f"no synonymous substitution reduces runs."
                )

        if not improved:
            break

    # ── Phase 2: Adjacent identical codons ───────────────────────────────
    n_homo = _n_homo(codons)
    for _ in range(max_iterations):
        adjs = find_identical_adjacent(codons)
        if not adjs:
            break

        cur_v = len(adjs)
        improved = False
        for i, j, c in adjs:
            best_idx = best_alt = None
            used_fallback = False
            for threshold_try, is_fallback in [(freq_threshold, False), (0.0, True)]:
                best_v = cur_v
                cand_idx, cand_alt = None, None
                for cand in [i, j]:
                    for alt in get_synonyms(codons[cand], codon_table, threshold_try):
                        trial = codons[:]
                        trial[cand] = alt
                        if _n_homo(trial) > n_homo:
                            continue
                        v = len(find_identical_adjacent(trial))
                        if v < best_v:
                            best_v, cand_idx, cand_alt = v, cand, alt
                if cand_idx is not None:
                    best_idx, best_alt, used_fallback = cand_idx, cand_alt, is_fallback
                    break
            if best_idx is not None:
                suffix = " (threshold override)" if used_fallback else ""
                old = codons[best_idx]
                codons[best_idx] = best_alt
                changes.append((best_idx, old, best_alt,
                                 f"adjacent ({i}&{j}: {c}){suffix}"))
                if used_fallback:
                    warnings.append(
                        f"Codons {i}&{j} ({c}): no synonym >= freq_threshold "
                        f"{freq_threshold} could fix this adjacent pair; used "
                        f"{best_alt} below threshold as a last resort."
                    )
                improved = True
                break

        if not improved:
            for i, j, c in adjs:
                warnings.append(
                    f"Codons {i}&{j} ({c}): adjacent — "
                    f"no synonym fixes without introducing homopolymer run."
                )
            break

    # ── Phase 3: Skip-1 identical codons ─────────────────────────────────
    n_homo  = _n_homo(codons)
    n_adj   = len(find_identical_adjacent(codons))
    for _ in range(max_iterations):
        skips = find_identical_skip1(codons)
        if not skips:
            break

        cur_v = len(skips)
        improved = False
        for i, j, c in skips:
            best_idx = best_alt = None
            used_fallback = False
            for threshold_try, is_fallback in [(freq_threshold, False), (0.0, True)]:
                best_v = cur_v
                cand_idx, cand_alt = None, None
                for cand in [i, j]:
                    for alt in get_synonyms(codons[cand], codon_table, threshold_try):
                        trial = codons[:]
                        trial[cand] = alt
                        if (_n_homo(trial) > n_homo or
                                len(find_identical_adjacent(trial)) > n_adj):
                            continue
                        v = len(find_identical_skip1(trial))
                        if v < best_v:
                            best_v, cand_idx, cand_alt = v, cand, alt
                if cand_idx is not None:
                    best_idx, best_alt, used_fallback = cand_idx, cand_alt, is_fallback
                    break
            if best_idx is not None:
                suffix = " (threshold override)" if used_fallback else ""
                old = codons[best_idx]
                codons[best_idx] = best_alt
                changes.append((best_idx, old, best_alt,
                                 f"skip-1 ({i}&{j}: {c}){suffix}"))
                if used_fallback:
                    warnings.append(
                        f"Codons {i}&{j} ({c}): no synonym >= freq_threshold "
                        f"{freq_threshold} could fix this skip-1 repeat; used "
                        f"{best_alt} below threshold as a last resort."
                    )
                improved = True
                break

        if not improved:
            for i, j, c in skips:
                warnings.append(
                    f"Codons {i}&{j} ({c}): skip-1 — "
                    f"no synonym fixes without introducing homopolymer run or adjacent pair."
                )
            break

    # ── Frequency upgrade ─────────────────────────────────────────────────
    n_homo_f = _n_homo(codons)
    n_adj_f  = len(find_identical_adjacent(codons))
    n_skip_f = len(find_identical_skip1(codons))
    for ci in range(len(codons)):
        cur_freq = next((f for _, ch in codon_table.items()
                         for c, f in ch if c == codons[ci]), 0.0)
        best_codon, best_freq = codons[ci], cur_freq
        for alt in get_synonyms(codons[ci], codon_table, freq_threshold):
            trial = codons[:]
            trial[ci] = alt
            if (_n_homo(trial) <= n_homo_f and
                    len(find_identical_adjacent(trial)) <= n_adj_f and
                    len(find_identical_skip1(trial)) <= n_skip_f):
                alt_freq = next((f for _, ch in codon_table.items()
                                 for c, f in ch if c == alt), 0.0)
                if alt_freq > best_freq:
                    best_freq, best_codon = alt_freq, alt
        if best_codon != codons[ci]:
            old = codons[ci]
            codons[ci] = best_codon
            changes.append((ci, old, best_codon, "freq upgrade"))

    remaining = {
        'homopolymer_runs':   find_homopolymer_runs(codons),
        'identical_adjacent': find_identical_adjacent(codons),
        'identical_skip1':    find_identical_skip1(codons),
    }
    return codons, changes, warnings, remaining


def optimize_from_aa(aa_sequence: str, codon_table_text: str,
                     freq_threshold: float = 0.10,
                     max_iterations: int = 150) -> dict:
    """Reverse-translate then optimize."""
    ct = parse_codon_table(codon_table_text)
    initial = reverse_translate(aa_sequence, ct, freq_threshold)
    final, changes, warnings, remaining = _run_optimization(
        initial, ct, freq_threshold, max_iterations)
    return dict(codon_table=ct, initial_codons=initial, final_codons=final,
                changes=changes, warnings=warnings, remaining_issues=remaining)


def optimize_from_nt(nt_sequence: str, codon_table_text: str,
                     freq_threshold: float = 0.10,
                     max_iterations: int = 150) -> dict:
    """Optimize an existing nucleotide sequence (read codon-by-codon from position 0)."""
    ct = parse_codon_table(codon_table_text)
    nt = re.sub(r'\s+', '', nt_sequence).upper()
    if len(nt) % 3 != 0:
        print(f"WARNING: sequence length {len(nt)} is not a multiple of 3; "
              f"last {len(nt) % 3} nt ignored.")
    initial = [nt[i:i+3] for i in range(0, len(nt) - len(nt) % 3, 3)]
    warn_pre = []
    for c in initial:
        if _aa_for_codon(c, ct) == '?':
            warn_pre.append(f"Codon '{c}' not found in table — will keep as-is.")
    final, changes, warnings, remaining = _run_optimization(
        initial, ct, freq_threshold, max_iterations)
    return dict(codon_table=ct, initial_codons=initial, final_codons=final,
                changes=changes, warnings=warn_pre + warnings,
                remaining_issues=remaining)


# ─────────────────────────────────────────────────────────────────────────────
# Display
# ─────────────────────────────────────────────────────────────────────────────

def display_results(result: dict, freq_threshold: float = 0.10):
    """Use rich if available, else plain text."""
    try:
        import rich  # noqa: F401
        _display_rich(result, freq_threshold)
    except ImportError:
        _display_plain(result, freq_threshold)


def _changed_and_remaining_sets(result: dict):
    changed = {idx for idx, *_ in result['changes']}
    issues = result['remaining_issues']
    rem_codons: set[int] = set()
    for _, start, run_len in issues.get('homopolymer_runs', []):
        for pos in range(start, start + run_len):
            rem_codons.add(pos // 3)
    for i, j, _ in issues.get('identical_adjacent', []):
        rem_codons |= {i, j}
    for i, j, _ in issues.get('identical_skip1', []):
        rem_codons |= {i, j}
    return changed, rem_codons


def _display_rich(result: dict, freq_threshold: float):
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text
    from rich import box
    from rich.panel import Panel

    console = Console()
    ct = result['codon_table']
    initial = result['initial_codons']
    final   = result['final_codons']
    changes = result['changes']
    warnings = result['warnings']
    issues  = result['remaining_issues']
    changed, rem_codons = _changed_and_remaining_sets(result)

    console.rule("[bold cyan]CODON OPTIMIZER RESULTS[/bold cyan]")

    # Initial sequence
    console.print()
    init_text = Text(' '.join(initial))
    console.print(Panel(init_text,
                        title="[bold yellow]Initial sequence[/bold yellow]",
                        border_style="yellow"))
    console.print(f"[dim]{''.join(initial)}[/dim]\n")

    # Change log
    if changes:
        t = Table(title="[bold]Codon Changes[/bold]", box=box.SIMPLE_HEAD)
        t.add_column("#",         style="dim",   width=4)
        t.add_column("Codon #",   justify="right")
        t.add_column("AA",        width=6)
        t.add_column("Original",  style="red bold")
        t.add_column("→ Replaced", style="green bold")
        t.add_column("Reason",    style="cyan")
        for n, (idx, old, new, reason) in enumerate(changes, 1):
            aa = _aa_for_codon(old, ct)
            t.add_row(str(n), str(idx), aa, old, new, reason)
        console.print(t)
    else:
        console.print("[green]No changes needed.[/green]")

    # Annotated sequence
    console.print()
    console.rule("[bold green]Annotated Final Sequence[/bold green]")
    console.print("[dim]Legend: [bold green]green[/bold green]=changed  "
                  "[bold red]red[/bold red]=remaining issue  "
                  "[bold magenta]magenta[/bold magenta]=changed but still flagged[/dim]\n")

    codon_line = Text()
    nt_line    = Text()
    for idx, codon in enumerate(final):
        in_changed = idx in changed
        in_rem     = idx in rem_codons
        if in_changed and in_rem:
            style = "bold magenta"
        elif in_changed:
            style = "bold green"
        elif in_rem:
            style = "bold red underline"
        else:
            style = "white"
        codon_line.append(f" {codon}", style=style)
        nt_line.append(codon, style=style)

    console.print(codon_line)
    console.print()
    console.print(Panel(nt_line,
                        title="[bold green]Final nucleotide sequence[/bold green]",
                        border_style="green"))

    # Remaining issues
    console.print()
    has_issues = any(len(v) > 0 for v in issues.values())
    if has_issues:
        console.rule("[bold red]Remaining Issues[/bold red]")
        for i, j, c in issues.get('identical_adjacent', []):
            aa = _aa_for_codon(c, ct)
            console.print(f"  [red]Adjacent:[/red]  codons {i} & {j}: "
                          f"[bold]{c}[/bold] | [bold]{c}[/bold]  ({aa})")
        for i, j, c in issues.get('identical_skip1', []):
            mid = final[i+1]
            aa = _aa_for_codon(c, ct)
            console.print(f"  [red]Skip-1:[/red]    codons {i} & {j}: "
                          f"[bold]{c}[/bold] _ {mid} _ [bold]{c}[/bold]  ({aa})")
        for nt, start, run_len in issues.get('homopolymer_runs', []):
            ci = sorted(set(range(start // 3, (start + run_len - 1) // 3 + 1)))
            console.print(f"  [red]Homopolymer:[/red] '{nt * run_len}' @{start} "
                          f"(codon ~{ci})")
    else:
        console.print(Panel("[bold green]✓ All violations resolved![/bold green]",
                            border_style="green"))

    # Warnings
    if warnings:
        console.print()
        console.rule("[bold yellow]Warnings[/bold yellow]")
        for w in warnings:
            console.print(f"  [yellow]⚠[/yellow]  {w}")

    console.print()
    console.rule()
    final_nt = ''.join(final)
    console.print(f"[bold]Changes:[/bold] {len(changes)}  |  "
                  f"[bold]Warnings:[/bold] {len(warnings)}  |  "
                  f"[bold]Threshold:[/bold] {freq_threshold}")
    console.print(f"[bold]Final:[/bold]   {final_nt}")
    console.print(f"[bold]Length:[/bold]  {len(final_nt)} nt  ({len(final)} codons)")


def _display_plain(result: dict, freq_threshold: float):
    ct      = result['codon_table']
    initial = result['initial_codons']
    final   = result['final_codons']
    changes = result['changes']
    warnings = result['warnings']
    issues  = result['remaining_issues']
    changed, rem_codons = _changed_and_remaining_sets(result)

    SEP = "=" * 72
    print(SEP)
    print("  CODON OPTIMIZER RESULTS")
    print(SEP)
    print("\nINITIAL SEQUENCE:")
    print(' '.join(initial))
    print(''.join(initial))

    if changes:
        print(f"\nCODON CHANGES ({len(changes)}):")
        print(f"  {'#':>3}  {'Idx':>4}  {'AA':<6}  {'Original':<8}  {'New':<8}  Reason")
        print("  " + "-" * 64)
        for n, (idx, old, new, reason) in enumerate(changes, 1):
            aa = _aa_for_codon(old, ct)
            print(f"  {n:>3}  {idx:>4}  {aa:<6}  {old:<8}  {new:<8}  {reason}")
    else:
        print("\nNo changes needed.")

    print("\nFINAL SEQUENCE (changed=*codon*, remaining issue=[codon]):")
    annotated = []
    for i, c in enumerate(final):
        if i in changed and i in rem_codons:
            annotated.append(f"!{c}!")
        elif i in changed:
            annotated.append(f"*{c}*")
        elif i in rem_codons:
            annotated.append(f"[{c}]")
        else:
            annotated.append(c)
    print(' '.join(annotated))
    print("\nFINAL SEQUENCE (plain):")
    print(''.join(final))

    has_issues = any(len(v) > 0 for v in issues.values())
    if has_issues:
        print("\nREMAINING ISSUES:")
        for i, j, c in issues.get('identical_adjacent', []):
            print(f"  Adjacent:  codons {i}&{j}: {c}|{c}  ({_aa_for_codon(c,ct)})")
        for i, j, c in issues.get('identical_skip1', []):
            mid = final[i+1]
            print(f"  Skip-1:    codons {i}&{j}: {c}_{mid}_{c}  ({_aa_for_codon(c,ct)})")
        for nt, start, run_len in issues.get('homopolymer_runs', []):
            ci = sorted(set(range(start // 3, (start + run_len - 1) // 3 + 1)))
            print(f"  Homopolymer: '{nt * run_len}' @{start} (codon ~{ci})")
    else:
        print("\nAll violations resolved!")

    if warnings:
        print("\nWARNINGS:")
        for w in warnings:
            print(f"  WARNING: {w}")

    final_nt = ''.join(final)
    print(f"\n{SEP}")
    print(f"Changes: {len(changes)}  |  Warnings: {len(warnings)}  |  Threshold: {freq_threshold}")
    print(f"Final:   {final_nt}")
    print(f"Length:  {len(final_nt)} nt  ({len(final)} codons)")
    print(SEP)


# ─────────────────────────────────────────────────────────────────────────────
# ██  USER CONFIGURATION  ─────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────

CODON_TABLE_TEXT = """
AmAcid  Codon      Number    /1000     Fraction

Gly     GGG     50527.00     11.12      0.15
Gly     GGA     39036.00      8.59      0.12
Gly     GGT    114185.00     25.14      0.34
Gly     GGC    130043.00     28.63      0.39

Glu     GAG     83804.00     18.45      0.32
Glu     GAA    179460.00     39.51      0.68
Asp     GAT    146794.00     32.32      0.63
Asp     GAC     87759.00     19.32      0.37

Val     GTG    115687.00     25.47      0.36
Val     GTA     51020.00     11.23      0.16
Val     GTT     86572.00     19.06      0.27
Val     GTC     67356.00     14.83      0.21

Ala     GCG    146264.00     32.20      0.34
Ala     GCA     93390.00     20.56      0.22
Ala     GCT     73677.00     16.22      0.17
Ala     GCC    113412.00     24.97      0.27

Arg     AGG      7423.00      1.63      0.03
Arg     AGA     12345.00      2.72      0.05
Ser     AGT     41544.00      9.15      0.15
Ser     AGC     70867.00     15.60      0.26

Lys     AAG     51685.00     11.38      0.25
Lys     AAA    156169.00     34.38      0.75
Asn     AAT     84846.00     18.68      0.46
Asn     AAC     98018.00     21.58      0.54

Met     ATG    123604.00     27.21      1.00
Ile     ATA     24233.00      5.34      0.09
Ile     ATT    135873.00     29.92      0.50
Ile     ATC    111878.00     24.63      0.41

Thr     ACG     63696.00     14.02      0.26
Thr     ACA     35995.00      7.93      0.15
Thr     ACT     43256.00      9.52      0.18
Thr     ACC    103121.00     22.70      0.42

Trp     TGG     65630.00     14.45      1.00
End     TGA      4428.00      0.97      0.30
Cys     TGT     23461.00      5.17      0.45
Cys     TGC     28747.00      6.33      0.55

End     TAG      1172.00      0.26      0.08
End     TAA      9006.00      1.98      0.62
Tyr     TAT     75774.00     16.68      0.58
Tyr     TAC     55847.00     12.30      0.42

Leu     TTG     60322.00     13.28      0.13
Leu     TTA     62823.00     13.83      0.13
Phe     TTT    100128.00     22.05      0.57
Phe     TTC     74885.00     16.49      0.43

Ser     TCG     39546.00      8.71      0.15
Ser     TCA     35837.00      7.89      0.13
Ser     TCT     42367.00      9.33      0.16
Ser     TCC     40365.00      8.89      0.15

Arg     CGG     25751.00      5.67      0.10
Arg     CGA     16607.00      3.66      0.07
Arg     CGT     93997.00     20.70      0.37
Arg     CGC     96053.00     21.15      0.38

Gln     CAG    130898.00     28.82      0.66
Gln     CAA     67129.00     14.78      0.34
His     CAT     57585.00     12.68      0.57
His     CAC     43743.00      9.63      0.43

Leu     CTG    231373.00     50.94      0.49
Leu     CTA     18067.00      3.98      0.04
Leu     CTT     51442.00     11.33      0.11
Leu     CTC     48147.00     10.60      0.10

Pro     CCG    101467.00     22.34      0.51
Pro     CCA     38663.00      8.51      0.20
Pro     CCT     32678.00      7.19      0.17
Pro     CCC     24383.00      5.37      0.12
"""

# ── Input: choose ONE of the two options below ────────────────────────────────

# Option A: amino acid sequence (1-letter or 3-letter codes)
# Leave as "" to use Option B instead.
AMINO_ACID_SEQUENCE = "LAQAKALVDAYLAAIRAHDGAAAAALMAQLTALDPNDVLAAIDAKIADPSVSPGVKRQLTLIKTLIADNKPIKAVEAAAQAVREGRRDALPEALAIIDEIL"

# Option B: existing nucleotide sequence to check and fix
# Leave as "" to use Option A instead.
INPUT_NT_SEQUENCE = (
    ""
)

# Minimum codon usage fraction; synonyms below this are never used as replacements
FREQ_THRESHOLD = 0.10

# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if not AMINO_ACID_SEQUENCE and not INPUT_NT_SEQUENCE:
        raise ValueError("Set either AMINO_ACID_SEQUENCE or INPUT_NT_SEQUENCE.")

    if INPUT_NT_SEQUENCE:
        result = optimize_from_nt(INPUT_NT_SEQUENCE, CODON_TABLE_TEXT,
                                  freq_threshold=FREQ_THRESHOLD)
    else:
        result = optimize_from_aa(AMINO_ACID_SEQUENCE, CODON_TABLE_TEXT,
                                  freq_threshold=FREQ_THRESHOLD)

    display_results(result, freq_threshold=FREQ_THRESHOLD)
