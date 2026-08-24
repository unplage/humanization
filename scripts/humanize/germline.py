"""Germline database handling and framework selection.

Uses the NCBI IgBLAST germline V/J FASTA files (human_gl_V.fasta,
human_gl_J.fasta) as the authoritative germline source. Sequences are
numbered with the portable Kabat engine, so matching is done position-wise
in Kabat space (no alignment needed).

Selection strategy (industry standard "hybrid" approach):
  1. Rank human V genes by FRAMEWORK identity (FR1+FR2+FR3, CDRs excluded).
  2. Among the top framework matches, prefer the gene with the HIGHEST CDR
     identity (minimizes required back-mutations and structural perturbation).
  3. J genes ranked by identity and FR4 compatibility.
"""

from __future__ import annotations

import os
import re
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .numbering import NumberedChain, number_heavy, number_light


GERMLINE_URLS = {
    # NCBI IgBLAST release directory (IMGT-derived human germline sets)
    "v": "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/database/",
    "j": "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/database/",
}
# Actual file names inside the IgBLAST tarball's database/ directory
IGHV_FILE = "human_gl_V.fasta"
IGHD_FILE = "human_gl_D.fasta"
IGLJ_FILE = "human_gl_J.fasta"
IGKV_IGLV_FILE = "human_gl_L.fasta"   # contains both kappa and lambda V
IGKJ_IGLJ_FILE = "human_gl_J.fasta"


@dataclass
class GermlineGene:
    gene_id: str        # e.g. IGHV1-69*01, IGKV1-39*01, IGHJ4*01
    chain_type: str     # "H" or "L"
    kind: str           # "V" or "J"
    sequence: str
    numbered: Optional[NumberedChain] = None

    @property
    def family(self) -> str:
        m = re.match(r"IG[HL][KV]?\d+", self.gene_id)
        return m.group(0) if m else self.gene_id


@dataclass
class GermlineDB:
    v_genes: List[GermlineGene]
    j_genes: List[GermlineGene]
    source_dir: str = ""

    def v_for(self, chain_type: str) -> List[GermlineGene]:
        return [g for g in self.v_genes if g.chain_type == chain_type]

    def j_for(self, chain_type: str) -> List[GermlineGene]:
        return [g for g in self.j_genes if g.chain_type == chain_type]

    def human(self, chain_type: str) -> List[GermlineGene]:
        """Only genes whose id contains 'IGHV'/'IGKV'/'IGLV' (human NCBI set)."""
        pat = "IGHV" if chain_type == "H" else "IG[KL]V"
        return [g for g in self.v_genes if g.chain_type == chain_type and re.match(pat, g.gene_id)]


# ---------------------------------------------------------------------------
# Parsing / loading
# ---------------------------------------------------------------------------

def _parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    records = []
    name, lines = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(lines)))
            name = line[1:].split()[0].strip()
            lines = []
        else:
            lines.append(line)
    if name is not None:
        records.append((name, "".join(lines)))
    return records


def _classify_gene(gene_id: str, chain_type: str, kind: str, seq: str) -> GermlineGene:
    g = GermlineGene(gene_id=gene_id, chain_type=chain_type, kind=kind, sequence=seq)
    try:
        if chain_type == "H":
            g.numbered = number_heavy(seq)
        else:
            g.numbered = number_light(seq)
    except ValueError:
        g.numbered = None
    return g


def load_germline_db(db_dir: str) -> GermlineDB:
    """Load the germline dataset. Priority:
      1. NCBI IgBLAST FASTA files in db_dir (human_gl_*.fasta)
      2. bundled Kabat-space JSON (data/germline/human_germline_kabat.json)
    Returns an empty-ish GermlineDB on total failure (caller falls back)."""
    # NCBI FASTA source
    vh_path = os.path.join(db_dir, IGHV_FILE)
    vl_path = os.path.join(db_dir, IGKV_IGLV_FILE)
    if os.path.exists(vh_path) or os.path.exists(vl_path):
        return _load_igblast_fasta(db_dir)
    # bundled JSON source
    bundled = os.path.join(os.path.dirname(os.path.dirname(__file__)), "..", "data", "germline", "human_germline_kabat.json")
    if os.path.exists(bundled):
        return _load_bundled_json(bundled)
    raise FileNotFoundError(
        f"no germline data found in {db_dir} (NCBI FASTA) nor at {bundled} (bundled)"
    )


def _load_bundled_json(path: str) -> GermlineDB:
    import json
    with open(path) as fh:
        blob = json.load(fh)
    genes: List[GermlineGene] = []
    for gene, pmap in blob["V"].items():
        ctype = "H" if gene.startswith("IGHV") else "L"
        seq = "".join(pmap.values())
        g = GermlineGene(gene, ctype, "V", seq)
        g.numbered = _chain_from_map(seq, pmap, ctype)
        genes.append(g)
    for gene, pmap in blob["J"].items():
        ctype = "H" if gene.startswith("IGHJ") else "L"
        seq = "".join(pmap.values())
        g = GermlineGene(gene, ctype, "J", seq)
        g.numbered = _chain_from_map(seq, pmap, ctype)
        genes.append(g)
    return GermlineDB(
        v_genes=[g for g in genes if g.kind == "V"],
        j_genes=[g for g in genes if g.kind == "J"],
        source_dir="bundled:" + path,
    )


def _chain_from_map(seq: str, pmap: dict, ctype: str) -> Optional[NumberedChain]:
    """Reconstruct a NumberedChain from an explicit {pos: aa} map."""
    from .numbering import NumberedChain, NumberedResidue, _cdr_segments
    res = []
    for i, (pos, aa) in enumerate(pmap.items()):
        region = _region_of_pos(pos)
        res.append(NumberedResidue(pos, aa, region, i))
    chain = NumberedChain(ctype, seq, res)
    chain.cdrs = _cdr_segments(chain)
    return chain


def _region_of_pos(pos: str) -> str:
    """Infer region from a Kabat position label (approx; used for germlines
    where the exact region walk is not available)."""
    from .numbering import NumberedResidue
    chain = pos[0]
    num = int("".join(c for c in pos if c.isdigit()))
    if chain == "H":
        if num <= 30:
            return "FR1"
        if num <= 35:
            return "CDR1"
        if num <= 49:
            return "FR2"
        if num <= 65:
            return "CDR2"
        if num <= 92:
            return "FR3"
        if num <= 102:
            return "CDR3"
        return "FR4"
    else:
        if num <= 23:
            return "FR1"
        if num <= 34:
            return "CDR1"
        if num <= 49:
            return "FR2"
        if num <= 56:
            return "CDR2"
        if num <= 88:
            return "FR3"
        if num <= 97:
            return "CDR3"
        return "FR4"


def _load_igblast_fasta(db_dir: str) -> GermlineDB:
    genes: List[GermlineGene] = []
    paths = [
        (os.path.join(db_dir, IGHV_FILE), r"IGHV", "H", "V"),
        (os.path.join(db_dir, IGKV_IGLV_FILE), r"IG[KL]V", "L", "V"),
        (os.path.join(db_dir, IGLJ_FILE), r"IGHJ", "H", "J"),
        (os.path.join(db_dir, IGLJ_FILE), r"IG[KL]J", "L", "J"),
    ]
    for path, pattern, ctype, kind in paths:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for gene_id, seq in _parse_fasta_text(fh.read()):
                if re.match(pattern, gene_id):
                    genes.append(_classify_gene(gene_id, ctype, kind, seq))
    if not genes:
        raise FileNotFoundError(
            f"no germline FASTA files found in {db_dir} "
            "(expected human_gl_V.fasta / human_gl_L.fasta / human_gl_J.fasta; "
            "run `humanize setup germline` to download)"
        )
    return GermlineDB(
        v_genes=[g for g in genes if g.kind == "V"],
        j_genes=[g for g in genes if g.kind == "J"],
        source_dir=db_dir,
    )


def download_germline_db(db_dir: str, url_base: Optional[str] = None) -> GermlineDB:
    """Download the germline FASTA files from NCBI IgBLAST release tarball.

    The human V/J germline FASTA files are distributed inside the IgBLAST
    binary tarball. We download and extract only the needed files.
    """
    os.makedirs(db_dir, exist_ok=True)
    tar = os.path.join(db_dir, "ncbi-igblast.tar.gz")
    latest = "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/LATEST/"
    # resolve the actual tarball name
    with urllib.request.urlopen(latest, timeout=60) as r:
        html = r.read().decode()
    m = re.search(r'href="(ncbi-igblast-[\d.]+-x64-linux\.tar\.gz)"', html)
    if not m:
        raise RuntimeError(f"could not locate IgBLAST tarball at {latest}")
    url = latest + m.group(1)
    print(f"[germline] downloading {url} ...")
    urllib.request.urlretrieve(url, tar)
    import tarfile
    with tarfile.open(tar) as tf:
        wanted = {IGHV_FILE, IGKV_IGLV_FILE, IGLJ_FILE}
        for member in tf.getmembers():
            if member.name.endswith("/database/" + IGHV_FILE) or \
               member.name.endswith("/database/" + IGKV_IGLV_FILE) or \
               member.name.endswith("/database/" + IGLJ_FILE):
                fname = os.path.basename(member.name)
                src = tf.extractfile(member)
                if src:
                    with open(os.path.join(db_dir, fname), "wb") as out:
                        out.write(src.read())
                    print(f"[germline] extracted {fname}")
    os.remove(tar)
    return load_germline_db(db_dir)


# ---------------------------------------------------------------------------
# Position-wise comparison (Kabat space)
# ---------------------------------------------------------------------------

FR_REGIONS = ("FR1", "FR2", "FR3")
CDR_REGIONS = ("CDR1", "CDR2")


def compare_to_germline(query: NumberedChain, gene: GermlineGene) -> Dict[str, float]:
    """Per-position comparison between query and a germline V gene.
    Returns identity fractions over FR, CDR1-2, and all."""
    if gene.numbered is None:
        return {"fr_identity": 0.0, "cdr_identity": 0.0, "all_identity": 0.0, "n_fr": 0, "n_cdr": 0}
    q, g = query.posmap(), gene.numbered.posmap()
    common_fr = [p for p in q if p in g and query.region_of(p) in FR_REGIONS]
    common_cdr = [p for p in q if p in g and query.region_of(p) in CDR_REGIONS]
    fr_id = sum(1 for p in common_fr if q[p] == g[p]) / len(common_fr) if common_fr else 0.0
    cdr_id = sum(1 for p in common_cdr if q[p] == g[p]) / len(common_cdr) if common_cdr else 0.0
    common = common_fr + common_cdr
    all_id = sum(1 for p in common if q[p] == g[p]) / len(common) if common else 0.0
    return {
        "fr_identity": round(fr_id, 4),
        "cdr_identity": round(cdr_id, 4),
        "all_identity": round(all_id, 4),
        "n_fr": len(common_fr),
        "n_cdr": len(common_cdr),
    }


def score_j_match(query: NumberedChain, jgene: GermlineGene) -> Tuple[float, int]:
    """Identity over the J region (FR4), which must be grafted from the
    human J gene. Position 103 (H) / 98 (L) onward."""
    if jgene.numbered is None:
        return 0.0, 0
    anchor_num = 103 if query.chain_type == "H" else 98
    q = query.posmap()
    g = jgene.numbered.posmap()
    n = 0
    same = 0
    for p in sorted(q.keys(), key=lambda x: (x[0], int("".join(c for c in x if c.isdigit())))):
        if int("".join(c for c in p if c.isdigit())) < anchor_num:
            continue
        if p in g and q[p] and g[p]:
            n += 1
            same += 1 if q[p] == g[p] else 0
    return (same / n if n else 0.0), n


# ---------------------------------------------------------------------------
# Germline selection
# ---------------------------------------------------------------------------

@dataclass
class GermlineChoice:
    v_gene: Optional[GermlineGene] = None
    j_gene: Optional[GermlineGene] = None
    scores: Dict[str, float] = field(default_factory=dict)
    runner_up: Optional[GermlineGene] = None
    alternatives: List[Tuple[GermlineGene, Dict[str, float]]] = field(default_factory=list)


def choose_germlines(
    query: NumberedChain,
    db: GermlineDB,
    n_alternatives: int = 5,
    j_pick: str = "identity",
) -> GermlineChoice:
    """Select the best human germline V and J genes for a query chain.

    Strategy (hybrid, per Olimpieri et al. 2015 / industry practice):
      * require FR identity >= 60% (below that, grafting is not viable);
      * primary ranking by FR identity;
      * among the top-FR hits, pick the gene that also maximizes CDR1+2
        identity (fewer back-mutations / less structural perturbation).
    """
    ctype = query.chain_type
    v_genes = db.human(ctype)
    scored = []
    for g in v_genes:
        s = compare_to_germline(query, g)
        scored.append((g, s))
    scored.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))
    viable = [t for t in scored if t[1]["fr_identity"] >= 0.60]
    if not viable:
        viable = scored
    top_fr = viable[:max(1, int(len(viable) * 0.3) or 1)]
    # within the top-FR set, pick max CDR identity (minimize back-mutations)
    best = max(top_fr, key=lambda t: (t[1]["cdr_identity"], t[1]["fr_identity"]))
    v_gene, v_scores = best
    runner_up = None
    if len(top_fr) > 1:
        alt = [t for t in top_fr if t[0] is not v_gene]
        if alt:
            runner_up = alt[0][0]

    # J gene: identity first (canonical J for the chosen V), else by length
    j_genes = db.j_for(ctype)
    j_scored = []
    for jg in j_genes:
        idn, n = score_j_match(query, jg)
        j_scored.append((jg, idn, n))
    j_scored.sort(key=lambda t: (-t[1], -t[2]))
    j_best = j_scored[0][0] if j_scored else None

    return GermlineChoice(
        v_gene=v_gene,
        j_gene=j_best,
        scores=v_scores,
        runner_up=runner_up,
        alternatives=[
            (g, s) for g, s in scored[:n_alternatives]
        ],
    )
