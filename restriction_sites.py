"""Restriction site detection and elimination.

A separate post-processing step on top of the main codon optimization: find
where selected restriction enzymes' recognition sites occur in a sequence,
and optionally eliminate them using the same synonymous-codon substitution
approach (frequency-threshold-first, with fallback override) used for
homopolymer/repeat fixing in reverse_translation.py.

Detection uses Biopython's Bio.Restriction.
"""

from __future__ import annotations
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, CommOnly

from reverse_translation import get_synonyms


def list_enzyme_names() -> list[str]:
    """Alphabetically sorted names of commercially available restriction enzymes."""
    return sorted(str(e) for e in CommOnly)


def find_sites(codons: list[str], enzyme_names: list[str]) -> dict[str, list[tuple[int, int]]]:
    """Find recognition site spans for the given enzymes in the sequence.

    Returns {enzyme_name: [(start0, end0), ...]}, 0-indexed nucleotide
    positions, end exclusive. Only enzymes with at least one site are included.
    """
    if not enzyme_names:
        return {}
    seq = Seq("".join(codons))
    batch = RestrictionBatch(enzyme_names)
    result = batch.search(seq)
    sites: dict[str, list[tuple[int, int]]] = {}
    for enzyme, cut_positions in result.items():
        spans = []
        for cut_pos in cut_positions:
            start0 = cut_pos - 1 - enzyme.fst5
            spans.append((start0, start0 + enzyme.size))
        if spans:
            sites[str(enzyme)] = spans
    return sites


def eliminate_sites(codons: list[str], codon_table: dict, enzyme_names: list[str],
                     freq_threshold: float = 0.10, max_iterations: int = 200
                     ) -> tuple[list[str], list[tuple[int, str, str, str]], list[str],
                                dict[str, list[tuple[int, int]]]]:
    """Remove recognition sites for the given enzymes via codon substitution.

    Same threshold-first-then-fallback strategy as the pattern-fixing phases
    in reverse_translation.py: tries synonyms >= freq_threshold first; only
    falls back to the full synonym pool if no threshold-respecting synonym
    can remove the site, tagging such substitutions "(threshold override)"
    with a matching warning. A site that no synonym (even ignoring frequency)
    can remove is reported as a warning and left in place.

    Returns (codons, changes, warnings, remaining_sites).
    """
    codons = list(codons)
    changes: list[tuple[int, str, str, str]] = []
    warnings: list[str] = []
    stuck: set[tuple[str, int]] = set()

    for _ in range(max_iterations):
        sites = find_sites(codons, enzyme_names)
        all_spans = [(name, s, e) for name, spans in sites.items()
                     for s, e in spans if (name, s) not in stuck]
        if not all_spans:
            break

        enzyme_name, start, end = all_spans[0]
        cand_idxs = sorted({p // 3 for p in range(start, end)
                             if 0 <= p // 3 < len(codons)})
        cur_total = sum(len(v) for v in sites.values())

        fixed = False
        for threshold_try, is_fallback in [(freq_threshold, False), (0.0, True)]:
            candidates = [(ci, alt) for ci in cand_idxs
                          for alt in get_synonyms(codons[ci], codon_table, threshold_try)]
            best_ci, best_alt, best_total = None, None, cur_total
            for ci, alt in candidates:
                trial = codons[:]
                trial[ci] = alt
                n_total = sum(len(v) for v in find_sites(trial, enzyme_names).values())
                if n_total < best_total:
                    best_total, best_ci, best_alt = n_total, ci, alt

            if best_ci is not None:
                old = codons[best_ci]
                codons[best_ci] = best_alt
                suffix = " (threshold override)" if is_fallback else ""
                changes.append((best_ci, old, best_alt,
                                 f"restriction site {enzyme_name} @{start}{suffix}"))
                if is_fallback:
                    warnings.append(
                        f"Restriction site {enzyme_name} @{start}: no synonym >= "
                        f"freq_threshold {freq_threshold} could remove it; used "
                        f"{best_alt} below threshold as a last resort."
                    )
                fixed = True
                break

        if not fixed:
            stuck.add((enzyme_name, start))
            warnings.append(
                f"Restriction site {enzyme_name} @{start}: no synonymous "
                f"substitution removes this site."
            )

    remaining = find_sites(codons, enzyme_names)
    return codons, changes, warnings, remaining
