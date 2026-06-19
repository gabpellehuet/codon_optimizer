"""Optional online lookup against the Kazusa Codon Usage Database.

This is a convenience for organisms not in the bundled codon_tables.py.
It requires an internet connection and depends on kazusa.or.jp staying up and
keeping its current HTML format -- for fully offline/reliable use, prefer the
bundled tables in codon_tables.py or a local file.
"""

from __future__ import annotations
import re
import urllib.parse
import urllib.request

BASE_URL = "https://www.kazusa.or.jp/codon/cgi-bin"
TIMEOUT = 10

ONE_TO_THREE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'End',
}

SEARCH_RE = re.compile(
    r'<A HREF="[^"]*?species=([\w.]+)">\s*<I>(.*?)</I>\s*\[(\w+)\]:\s*(\d+)</A>',
    re.S,
)
ENTRY_RE = re.compile(r'([UCAG]{3})\s+([A-Z*])\s+([\d.]+)\s+([\d.]+)\s+\(\s*(\d+)\s*\)')


def search_organisms(query: str, limit: int = 50) -> list[tuple[str, str]]:
    """Search Kazusa by organism name (case-insensitive substring match).

    Returns a list of (species_id, display_label) tuples.
    Raises urllib.error.URLError / socket.timeout on connection problems.
    """
    url = f"{BASE_URL}/spsearch.cgi?species={urllib.parse.quote(query)}&c=i"
    with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
        html = resp.read().decode("iso-8859-1")
    results = []
    for species_id, name, dataset, n_cds in SEARCH_RE.findall(html):
        label = f"{name.strip()} [{dataset}] ({n_cds} CDS)"
        results.append((species_id, label))
    return results[:limit]


def fetch_codon_table(species_id: str) -> str:
    """Fetch and parse a codon usage table for one species ID.

    Returns text in the same format as CODON_TABLE_TEXT (DNA triplets,
    3-letter amino acid codes), ready for parse_codon_table().
    """
    url = f"{BASE_URL}/showcodon.cgi?species={species_id}&aa=1&style=N"
    with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
        html = resp.read().decode("iso-8859-1")

    pre_match = re.search(r'<PRE>(.*?)</PRE>', html, re.S)
    if not pre_match:
        raise ValueError("Could not find a codon usage table in the Kazusa response.")

    entries = ENTRY_RE.findall(pre_match.group(1))
    if not entries:
        raise ValueError("No codon entries found for this organism.")

    lines = ["AmAcid  Codon      Number    /1000     Fraction", ""]
    for codon_rna, aa1, fraction, per1000, count in entries:
        codon_dna = codon_rna.replace('U', 'T')
        aa3 = ONE_TO_THREE[aa1]
        lines.append(f"{aa3:<8}{codon_dna}     {count:>8}.00     {per1000:>6}      {fraction}")
    return "\n".join(lines)
