"""Anchor-based Kabat numbering engine (portable, zero-dependency).

Implements *Kabat* numbering of antibody variable domains using conserved
anchor residues that are structurally invariant:

    VH:  Cys22 (FR1) . Trp36 (FR2) . Cys92 (FR3) . Trp103 (FR4)
    VL:  Cys23 (FR1) . Trp35 (FR2) . Cys88 (FR3) . Phe98 (FR4)

Framework lengths between anchors are fixed by the scheme:

    VH:  FR1 1-30, CDR1 31-35(+35A..), FR2 36-49, CDR2 50-65(+52A-D),
         FR3 66-92(+82A-C), CDR3 95-102(+100A-K), FR4 103-113
    VL:  FR1 1-23, CDR1 24-34(+27A-E), FR2 35-49, CDR2 50-56,
         FR3 57-88(+82A-C), CDR3 89-97, FR4 98-107

Anchor candidates are validated by (a) conserved-motif context and (b) the
total number of residues spanned, so exotic loop lengths and species quirks
(e.g. mouse FR3 ending 'YFC', VHH long CDR3) are handled robustly. The engine
is scheme-accurate for standard V domains; ANARCI mode (exact) is used on the
server for cross-validation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


@dataclass
class NumberedResidue:
    pos: str            # e.g. "H31", "H100A", "L27E"
    aa: str
    region: str
    index: int          # 0-based index into the original sequence

    def __repr__(self) -> str:
        return f"{self.pos}:{self.aa}"


@dataclass
class NumberedChain:
    chain_type: str     # "H" or "L"
    sequence: str
    residues: List[NumberedResidue]
    species_hint: str = "unknown"
    warnings: List[str] = field(default_factory=list)
    cdrs: Dict[str, Tuple[str, str]] = field(default_factory=dict)

    def residue(self, pos: str) -> Optional[NumberedResidue]:
        for r in self.residues:
            if r.pos == pos:
                return r
        return None

    def posmap(self) -> Dict[str, str]:
        return {r.pos: r.aa.upper() for r in self.residues}

    def region_of(self, pos: str) -> Optional[str]:
        r = self.residue(pos)
        return r.region if r else None

    def seq_range(self, start_pos: str, end_pos: str) -> str:
        start = self.residue(start_pos)
        end = self.residue(end_pos)
        if start is None or end is None or end.index < start.index:
            return ""
        return self.sequence[start.index : end.index + 1]


def _letters(n: int) -> List[str]:
    """Return n insertion letters A, B, ... Z, AA, ..."""
    out = []
    for i in range(n):
        a = i // 26
        b = i % 26
        out.append(chr(65 + b) if a == 0 else chr(64 + a) + chr(65 + b))
    return out


# ---------------------------------------------------------------------------
# Heavy chain
# ---------------------------------------------------------------------------

def _find_fr2_trp(s: str, start: int, end: int, chain_type: str) -> Optional[int]:
    """Find the FR2-start Trp (H: W36, L: W35). Prefers a W whose downstream
    context matches the FR2-start consensus (H: 'W[VFLI][RQ]' e.g. WVRQ/WVKQ/
    WFRQ; L: 'WYQ') and whose CDR1 length is canonical."""
    if chain_type == "H":
        ctx_pat = re.compile(r"W[VFLIY][RQ]")
        len_ok = lambda n: 4 <= n <= 16
        prefer = 7
    else:
        ctx_pat = re.compile(r"WY[QFL]")
        len_ok = lambda n: 6 <= n <= 17
        prefer = 11
    best, best_score = None, -1e9
    for m in re.finditer("W", s[start:end]):
        i = start + m.start()
        cdr1_len = i - start
        score = 0
        if ctx_pat.match(s[i : i + 3]):
            score += 10
        if len_ok(cdr1_len):
            score += 5 - abs(cdr1_len - prefer)
        if score > best_score:
            best, best_score = i, score
    return best


def number_heavy(seq: str, species_hint: str = "unknown") -> NumberedChain:
    s = seq.upper().strip()
    n = len(s)
    warnings: List[str] = []
    res: List[NumberedResidue] = []

    # ---- FR1 (1-30), Cys22 anchor ----
    c_idx = None
    for i in range(18, min(28, n)):
        if s[i] == "C":
            c_idx = i
            break
    if c_idx is None:
        raise ValueError(f"[VH] no conserved Cys22 anchor found in first 28 residues: {s[:28]}")
    if c_idx != 21:
        warnings.append(f"[VH] Cys22 anchor at index {c_idx + 1} (expected 22); FR1 length adjusted")
    fr1_len = c_idx + 9                  # Cys pos 22 -> FR1 ends at pos 30
    if fr1_len > n:
        raise ValueError("[VH] FR1 extends beyond sequence end")
    for i, aa in enumerate(s[:fr1_len]):
        res.append(NumberedResidue(f"H{i + 1}", aa, "FR1", i))

    # ---- CDR1 (31-35 + 35A..) and FR2 anchor Trp36 ----
    w_idx = _find_fr2_trp(s, fr1_len, min(fr1_len + 22, n), "H")
    if w_idx is None:
        raise ValueError("[VH] no Trp36 (FR2 start) anchor found")
    cdr1_len = w_idx - fr1_len
    if not (4 <= cdr1_len <= 16):
        warnings.append(f"[VH] unusual CDR1 length {cdr1_len} (expected 5-12)")
    cdr1_seq = s[fr1_len:w_idx]
    if len(cdr1_seq) <= 5:
        for j, aa in enumerate(cdr1_seq):
            res.append(NumberedResidue(f"H{31 + j}", aa, "CDR1", fr1_len + j))
    else:
        for j, aa in enumerate(cdr1_seq):
            if j < 5:
                res.append(NumberedResidue(f"H{31 + j}", aa, "CDR1", fr1_len + j))
            else:
                res.append(NumberedResidue(f"H35{_letters(j - 4)[-1]}", aa, "CDR1", fr1_len + j))

    # ---- FR2 (36-49, 13-14 residues) / CDR2 (50-65, 16 slots) / FR3 (66-92) ----
    # all anchored on Cys92; fixed slot count from W36 (inclusive) to C92
    # (inclusive) = 14 + 16 + 27 = 57; insertions 52A-D and 82A-C allowed.
    c92_idx = _find_c92(s, w_idx, n)
    if c92_idx is None:
        raise ValueError("[VH] no Cys92 (FR3 end) anchor found")
    seg = s[w_idx : c92_idx + 1]
    seg_len = len(seg)
    # FR2 length: 13 for VH3/VH5 (position 49 empty, "WVRQAPGKGLEWV") unless
    # camelid VHH (hallmark H37F/Y => full 14-residue FR2); 14 otherwise.
    fr2_pat = seg[7:9] if len(seg) > 9 else "  "
    vhh_37 = seg[1] in ("F", "Y") if len(seg) > 1 else False
    fr2_len = 13 if fr2_pat == "KG" and not vhh_37 else 14
    expected = fr2_len + 16 + 27
    diff = seg_len - expected
    if diff > 7 or diff < -2:
        warnings.append(f"[VH] FR2+CDR2+FR3 span {seg_len} vs {expected} expected (diff {diff:+d})")
    if diff < 0:
        fr2_len += diff              # rare: extra gap in FR2
    rem = seg_len - fr2_len
    # FR3 length: 30 by default (VH1/VH3/VH4 82A-C insertions), unless the
    # camelid VHH signature is present (no insertion, "MDSL" loop) or the
    # 'MNSL' VH3-family motif forces a 30-residue FR3.
    vhh_score = int(seg[1] in ("F", "Y")) + int(len(seg) > 8 and seg[8] in ("E", "Q")) \
        + int(len(seg) > 9 and seg[9] == "R") + int(len(seg) > 11 and seg[11] in ("G", "F", "S"))
    if "MNSL" in seg[-30:]:
        fr3_len = 30
    elif vhh_score >= 3:
        fr3_len = 29
    else:
        fr3_len = 30
    cdr2_len = rem - fr3_len
    if not (16 <= cdr2_len <= 20):
        if cdr2_len > 20:
            warnings.append(f"[VH] CDR2 length {cdr2_len} > 20 (unusual H2 loop)")
            cdr2_len = 20
            fr3_len = rem - cdr2_len
        else:
            warnings.append(f"[VH] CDR2 length {cdr2_len} < 16 (unusual; check FR2/FR3 split)")
            cdr2_len = max(16, cdr2_len)
            fr3_len = rem - cdr2_len
    cdr2_ins = max(0, cdr2_len - 16)
    fr3_ins = max(0, fr3_len - 27)
    if fr3_ins > 3:
        warnings.append(f"[VH] {fr3_ins} FR3 insertions (82A..) - unusual")
    if cdr2_ins > 4:
        warnings.append(f"[VH] {cdr2_ins} CDR2 insertions (52A..) - unusual")

    # FR2
    for j in range(fr2_len):
        res.append(NumberedResidue(f"H{36 + j}", seg[j], "FR2", w_idx + j))
    # CDR2
    off = fr2_len
    ins = max(0, cdr2_len - 16)
    for j in range(cdr2_len):
        if j < 3:
            pos = f"H{50 + j}"
        elif j < 3 + ins:
            pos = f"H52{_letters(j - 2)[-1]}"
        else:
            pos = f"H{50 + j - ins}"
        res.append(NumberedResidue(pos, seg[off + j], "CDR2", w_idx + off + j))
    # FR3
    off += cdr2_len
    for j in range(fr3_len):
        if j < 17:
            pos = f"H{66 + j}"
        elif j < 17 + fr3_ins:
            pos = f"H82{_letters(j - 16)[-1]}"
        else:
            pos = f"H{66 + j - fr3_ins}"
        res.append(NumberedResidue(pos, seg[off + j], "FR3", w_idx + off + j))

    # ---- CDR3 (95-102 + 100A..) and FR4 (103-113), W103 anchor ----
    j_start = c92_idx + 1
    w103_idx = _find_j_anchor(s, j_start, min(j_start + 35, n), "H")
    if w103_idx is None:
        raise ValueError("[VH] no Trp103 (J region start) anchor found")
    cdr3_seq = s[j_start:w103_idx]
    cdr3_len = len(cdr3_seq)
    ins = max(0, cdr3_len - 8)
    for j, aa in enumerate(cdr3_seq):
        if j < 6:
            pos = f"H{95 + j}"
        elif j < 6 + ins:
            pos = f"H100{_letters(j - 5)[-1]}"
        else:
            pos = f"H{101 + j - 6 - ins}"
        res.append(NumberedResidue(pos, aa, "CDR3", j_start + j))
    fr4 = s[w103_idx : w103_idx + 11]
    if len(fr4) < 11:
        warnings.append(f"[VH] FR4 shorter than expected ({len(fr4)}/11) - truncated J region?")
    for j, aa in enumerate(fr4):
        res.append(NumberedResidue(f"H{103 + j}", aa, "FR4", w103_idx + j))
    tail_end = w103_idx + 11
    if tail_end < n:
        warnings.append(f"[VH] {n - tail_end} trailing residue(s) after FR4 ignored: {s[tail_end:]}")

    chain = NumberedChain("H", s, res, species_hint, warnings)
    chain.cdrs = _cdr_segments(chain)
    return chain


def _find_c92(s: str, w_idx: int, n: int) -> Optional[int]:
    """Locate Cys92: prefer 'Y[YF]C' motif at span 57 (+insertions) after W36."""
    best = None
    for m in re.finditer(r"Y[YF]?C", s[w_idx : min(w_idx + 130, n)]):
        idx = w_idx + m.start() + len(m.group(0)) - 1
        rel = idx - w_idx
        if 57 - 3 <= rel <= 57 + 8:
            if best is None:
                best = idx
            elif abs(rel - 57) < abs(best - w_idx - 57):
                best = idx
    if best is not None:
        return best
    for m in re.finditer("C", s[w_idx : min(w_idx + 130, n)]):
        idx = w_idx + m.start()
        if 57 - 3 <= idx - w_idx <= 57 + 8:
            return idx
    return None


def _find_j_anchor(s: str, start: int, end: int, chain_type: str) -> Optional[int]:
    """Find the J-region anchor (H: W of 'WGXG..', L: F of 'FGQG..').
    Uses the LAST plausible match so internal CDR3 motifs are not mistaken
    for the J region."""
    if chain_type == "H":
        pats = [r"W[GRASLVQFY]G", r"W[GRASLVQFY][QGS]"]
    else:
        pats = [r"FG[GQSTPE][GTKQE]", r"F[GSTPE][GQSTPEK]"]
    for pat in pats:
        ms = list(re.finditer(pat, s[start:end]))
        if ms:
            return start + ms[-1].start()
    # fallback: the anchor is the first F of the J region (<=3 residues past
    # the CDR3 end), typically right after the C-terminal Cys of the FR3
    for i in range(start, min(start + 4, end)):
        if s[i] == "F":
            return i
    return None


# ---------------------------------------------------------------------------
# Light chain (kappa / lambda)
# ---------------------------------------------------------------------------

def number_light(seq: str, species_hint: str = "unknown") -> NumberedChain:
    s = seq.upper().strip()
    n = len(s)
    warnings: List[str] = []
    res: List[NumberedResidue] = []

    # ---- FR1 (1-23), Cys23 anchor ----
    c_idx = None
    for i in range(20, min(27, n)):
        if s[i] == "C":
            c_idx = i
            break
    if c_idx is None:
        raise ValueError(f"[VL] no conserved Cys23 anchor found in first 27 residues: {s[:27]}")
    if c_idx != 22:
        warnings.append(f"[VL] Cys23 anchor at index {c_idx + 1} (expected 23); FR1 length adjusted")
    fr1_len = c_idx + 1
    for i, aa in enumerate(s[:fr1_len]):
        res.append(NumberedResidue(f"L{i + 1}", aa, "FR1", i))

    # ---- CDR1 (24-34 + 27A..) and FR2 anchor Trp35 ----
    w_idx = _find_fr2_trp(s, fr1_len, min(fr1_len + 24, n), "L")
    if w_idx is None:
        raise ValueError("[VL] no Trp35 (FR2 start) anchor found")
    cdr1_len = w_idx - fr1_len
    if not (6 <= cdr1_len <= 17):
        warnings.append(f"[VL] unusual CDR1 length {cdr1_len} (expected 11-17)")
    cdr1_seq = s[fr1_len:w_idx]
    cdr1_len = len(cdr1_seq)
    ins = max(0, cdr1_len - 11)
    for j, aa in enumerate(cdr1_seq):
        if j < 4:
            pos = f"L{24 + j}"
        elif j < 4 + ins:
            pos = f"L27{_letters(j - 3)[-1]}"
        else:
            pos = f"L{28 + j - 4 - ins}"
        res.append(NumberedResidue(pos, aa, "CDR1", fr1_len + j))

    # ---- FR2 (35-49, 15) / CDR2 (50-56, 7) / FR3 (57-88, 32) ----
    # fixed span from W35 (inclusive) to C88 (inclusive) = 15+7+32 = 54
    c88_idx = _find_c88(s, w_idx, n)
    if c88_idx is None:
        raise ValueError("[VL] no Cys88 (FR3 end) anchor found")
    seg = s[w_idx : c88_idx + 1]
    seg_len = len(seg)
    diff = seg_len - 54
    if diff > 4 or diff < -2:
        warnings.append(f"[VL] FR2+CDR2+FR3 span {seg_len} vs 54 expected (diff {diff:+d})")
    fr2_len = 15
    if diff < 0:
        fr2_len = 15 + diff
    rem = seg_len - fr2_len
    cdr2_len = min(7, rem)
    fr3_len = rem - cdr2_len
    fr3_ins = max(0, fr3_len - 32)
    if fr3_ins > 3:
        warnings.append(f"[VL] {fr3_ins} FR3 insertions (82A..) - unusual")
    for j in range(fr2_len):
        res.append(NumberedResidue(f"L{35 + j}", seg[j], "FR2", w_idx + j))
    off = fr2_len
    for j in range(cdr2_len):
        res.append(NumberedResidue(f"L{50 + j}", seg[off + j], "CDR2", w_idx + off + j))
    off += cdr2_len
    for j in range(fr3_len):
        if j < 26:
            pos = f"L{57 + j}"
        elif j < 26 + fr3_ins:
            pos = f"L82{_letters(j - 25)[-1]}"
        else:
            pos = f"L{57 + j - fr3_ins}"
        res.append(NumberedResidue(pos, seg[off + j], "FR3", w_idx + off + j))

    # ---- CDR3 (89-97) and FR4 (98-107), Phe98 anchor ----
    j_start = c88_idx + 1
    f98_idx = _find_j_anchor(s, j_start, min(j_start + 16, n), "L")
    if f98_idx is None:
        raise ValueError("[VL] no Phe98 (J region start) anchor found")
    cdr3_seq = s[j_start:f98_idx]
    for j, aa in enumerate(cdr3_seq):
        if j < 9:
            res.append(NumberedResidue(f"L{89 + j}", aa, "CDR3", j_start + j))
        else:
            res.append(NumberedResidue(f"L95{_letters(j - 7)[-1]}", aa, "CDR3", j_start + j))
    fr4 = s[f98_idx : f98_idx + 10]
    if len(fr4) < 10:
        warnings.append(f"[VL] FR4 shorter than expected ({len(fr4)}/10) - truncated J region?")
    for j, aa in enumerate(fr4):
        res.append(NumberedResidue(f"L{98 + j}", aa, "FR4", f98_idx + j))
    tail_end = f98_idx + 10
    if tail_end < n:
        warnings.append(f"[VL] {n - tail_end} trailing residue(s) after FR4 ignored: {s[tail_end:]}")

    chain = NumberedChain("L", s, res, species_hint, warnings)
    chain.cdrs = _cdr_segments(chain)
    return chain


def _find_c88(s: str, w_idx: int, n: int) -> Optional[int]:
    """Locate Cys88: prefer 'Y[YF]C' motif at span 54 (+insertions) after W35."""
    best = None
    for m in re.finditer(r"Y[YF]?C", s[w_idx : min(w_idx + 90, n)]):
        idx = w_idx + m.start() + len(m.group(0)) - 1
        rel = idx - w_idx
        if 54 - 2 <= rel <= 54 + 5:
            if best is None or abs(rel - 54) < abs(best - w_idx - 54):
                best = idx
    if best is not None:
        return best
    for m in re.finditer("C", s[w_idx : min(w_idx + 90, n)]):
        idx = w_idx + m.start()
        if 54 - 2 <= idx - w_idx <= 54 + 5:
            return idx
    return None


# ---------------------------------------------------------------------------
# CDR segment definitions per scheme (position ranges on the Kabat backbone)
# ---------------------------------------------------------------------------

KABAT_CDRS = {
    "H": {"CDR1": ("H31", "H35"), "CDR2": ("H50", "H65"), "CDR3": ("H95", "H102")},
    "L": {"CDR1": ("L24", "L34"), "CDR2": ("L50", "L56"), "CDR3": ("L89", "L97")},
}
CHOTHIA_CDRS = {
    "H": {"CDR1": ("H26", "H35"), "CDR2": ("H52", "H56"), "CDR3": ("H95", "H102")},
    "L": {"CDR1": ("L24", "L34"), "CDR2": ("L50", "L56"), "CDR3": ("L89", "L97")},
}
ABM_CDRS = {
    "H": {"CDR1": ("H26", "H35"), "CDR2": ("H50", "H58"), "CDR3": ("H95", "H102")},
    "L": {"CDR1": ("L24", "L34"), "CDR2": ("L50", "L56"), "CDR3": ("L89", "L97")},
}
IMGT_CDRS = {
    "H": {"CDR1": ("H26", "H34"), "CDR2": ("H56", "H65"), "CDR3": ("H95", "H102")},
    "L": {"CDR1": ("L24", "L34"), "CDR2": ("L50", "L56"), "CDR3": ("L89", "L97")},
}
CDR_SCHEMES = {"kabat": KABAT_CDRS, "chothia": CHOTHIA_CDRS, "abm": ABM_CDRS, "imgt": IMGT_CDRS}


def _cdr_segments(chain: NumberedChain, scheme: str = "kabat") -> Dict[str, Tuple[str, str]]:
    """Return {CDR1: (start_pos, end_pos), ...}; end positions include
    insertion codes (e.g. H35B, H100D)."""
    cdr_def = CDR_SCHEMES[scheme][chain.chain_type]
    out = {}
    for name, (start, _end) in cdr_def.items():
        start_res = chain.residue(start)
        if start_res is None:
            continue
        region = start_res.region
        last = start_res
        idx = chain.residues.index(start_res)
        while idx + 1 < len(chain.residues) and chain.residues[idx + 1].region == region:
            idx += 1
            last = chain.residues[idx]
        out[name] = (start, last.pos)
    return out


def cdr_regions(chain: NumberedChain) -> Dict[str, str]:
    return {k: f"{v[0]}..{v[1]}" for k, v in chain.cdrs.items()}


# ---------------------------------------------------------------------------
# VHH detection (single-domain camelid)
# ---------------------------------------------------------------------------

VHH_HALLMARK = {
    37: ("F", "Y"),     # Kabat 37: VH3 V -> VHH F/Y
    44: ("E", "Q"),     # 44: G -> E/Q
    45: ("R",),         # 45: L -> R
    47: ("G", "F", "S"),  # 47: W -> G/F/S
}


def is_vhh(chain: NumberedChain) -> Tuple[bool, float, List[str]]:
    """Detect camelid single-domain (VHH) hallmark in FR2 (Kabat 37/44/45/47).
    Returns (is_vhh, matched_count, matched_positions)."""
    matched = []
    for pos, allowed in VHH_HALLMARK.items():
        r = chain.residue(f"H{pos}")
        if r is not None and r.aa.upper() in allowed:
            matched.append(f"H{pos}{r.aa}")
    return (len(matched) >= 3, len(matched), matched)
