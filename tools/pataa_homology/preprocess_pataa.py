#!/usr/bin/env python3
"""PATAA patent protein DB -> antibody CDR reference database.

Extracts antibody V-domain sequences from the GenBank Patent Division BLAST
database (PATAA), detects Ig domains with IgBLAST (igblastp, human germline),
numbers them with ANARCI (IMGT) and builds a CDR JSON database compatible with
the `tools/cdr_homology` query workflow.

Pipeline stages (all resumable; delete a stage output to re-run it):

  dump     blastdbcmd -entry all -target_only  -> pataa.fasta
  split    assign unique synthetic ids (p00000001...) -> chunks/ + index.tsv
  detect   parallel igblastp (human_V, outfmt 7); keep max V bitscore >= THRESH
  extract  collect positive sequences -> positives.fasta
  number   parallel ANARCI (IMGT) -> cdr_records.jsonl
  repair   re-extract CDRs from stored numbering, then build (no ANARCI)
  build    dedup + patent metadata -> data/pataa_cdrs.json + pataa_stats.json

Usage:
  python3 tools/pataa_homology/preprocess_pataa.py all
  python3 tools/pataa_homology/preprocess_pataa.py detect --workers 64
  python3 tools/pataa_homology/preprocess_pataa.py repair   # after CDR fixes

Requires: IgBLAST 1.22.0 (data/igblast/ or IGBLAST_HOME), ANARCI + HMMER,
BioPython. The PATAA BLAST DB is auto-extracted from database/pataa.tar.gz.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else default


def _tool_path(configured: Path, name: str) -> Path:
    """Use the configured path, else fall back to PATH."""
    if configured.exists():
        return configured
    found = shutil.which(name)
    return Path(found) if found else configured


IGB = _env_path("IGBLAST_HOME",
                REPO / "data" / "igblast" / "ncbi-igblast-1.22.0")
IGDATA = _env_path("IGDATA", IGB)
BLASTDBCMD = _tool_path(_env_path("BLASTDBCMD", IGB / "bin" / "blastdbcmd"),
                        "blastdbcmd")
IGBLASTP = _tool_path(_env_path("IGBLASTP", IGB / "bin" / "igblastp"),
                      "igblastp")
GERMLINE_V = _env_path("GERMLINE_V",
                       IGDATA / "internal_data" / "human" / "human_V")

DEFAULT_DB = REPO / "tools" / "pataa_homology" / "database" / "pataa"
WORK = REPO / "data" / "pataa_work"
CHUNKS = WORK / "chunks"
IGOUT = WORK / "igout"
ANOUT = WORK / "anarci"

FASTA = WORK / "pataa.fasta"
INDEX_TSV = WORK / "index.tsv"
POSITIVES_TSV = WORK / "positives.tsv"
POSITIVES_FASTA = WORK / "positives.fasta"
CDR_RECORDS = WORK / "cdr_records.jsonl"
OUT_JSON = REPO / "data" / "pataa_cdrs.json"
OUT_STATS = REPO / "data" / "pataa_stats.json"

THRESHOLD = 40.0          # min max-V bitscore to call a sequence antibody-like
CHUNK_SIZE = 20000        # sequences per detection/split chunk
NUMBER_BATCH = 2000       # sequences per ANARCI worker task
# minimal CDR lengths: drop numbering artifacts / partial domains
MIN_CDR = {"CDR1": 4, "CDR2": 2, "CDR3": 4}
_IMGT_CDR_RANGES = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}
# A V-domain CDR3 is rarely longer than 40 aa. Longer "CDR3"s come from
# fusion constructs (toxin/cytokine + GGGGS linker) where ANARCI numbers the
# C-terminal partner as 111A..111Z / 112Z..112A insertions.
MAX_CDR3 = 40


def log(msg: str):
    print(f"[pataa] {msg}", flush=True)


def _require_tool(path: Path, name: str):
    if path.exists():
        return
    raise RuntimeError(
        f"{name} not found: {path}\n"
        "Install IgBLAST 1.22.0 under data/igblast/ or set IGBLAST_HOME "
        "(overrides: BLASTDBCMD, IGBLASTP, GERMLINE_V, IGDATA).")


def _require_germline_v():
    path = Path(str(GERMLINE_V) + ".nhr")
    if path.exists():
        return
    raise RuntimeError(
        f"IgBLAST human V germline DB not found: {path}\n"
        "Set IGBLAST_HOME or GERMLINE_V to the human_V BLAST DB prefix.")


def ensure_pataa_db(db: Path = DEFAULT_DB) -> Path:
    """Return an extracted PATAA BLAST DB prefix, extracting pataa.tar.gz."""
    if db.with_suffix(".pin").exists():
        return db
    archive = db.parent / "pataa.tar.gz"
    if not archive.exists():
        raise RuntimeError(
            f"PATAA BLAST DB not found: {db}.p*\n"
            f"Copy pataa.tar.gz into {db.parent}/ (it is auto-extracted), "
            "or pass 'dump --db <prefix>'.")
    db.parent.mkdir(parents=True, exist_ok=True)
    log(f"dump: extracting {archive} -> {db.parent}")
    root = db.parent.resolve()
    with tarfile.open(archive, "r:gz") as tf:
        for member in tf.getmembers():
            target = (root / member.name).resolve()
            if target != root and root not in target.parents:
                raise RuntimeError(
                    f"unsafe path in {archive.name}: {member.name}")
        tf.extractall(root)
    if not db.with_suffix(".pin").exists():
        raise RuntimeError(
            f"{archive} extracted but {db}.pin is still missing; "
            "check the archive layout")
    return db


def iter_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].rstrip("\n")
                buf = []
            else:
                line = line.strip()
                if line:
                    buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def parse_defline(defline: str) -> dict:
    out = {"accession": "", "patents": [], "species": ""}
    m = re.match(r"\s*([A-Za-z0-9_.]+)", defline)
    if m:
        out["accession"] = m.group(1).split(".")[0]
    for pm in re.finditer(
        r"Sequence\s+\d+\s+from\s+patent\s+([A-Z]{2}\s*\d+(?:\s*B\d)?)", defline
    ):
        p = re.sub(r"\s+", " ", pm.group(1)).strip()
        if p not in out["patents"]:
            out["patents"].append(p)
    for sm in re.findall(r"\[([^\]]+)\]", defline):
        s = sm.strip()
        if s and s.lower() not in ("", "unidentified"):
            out["species"] = s
            break
    return out


# ---------------------------------------------------------------------------
# stage: dump
# ---------------------------------------------------------------------------

def stage_dump(db: Path):
    if FASTA.exists() and FASTA.stat().st_size > 0:
        log(f"dump: exists ({FASTA.stat().st_size/1e9:.2f} GB), skip")
        return
    _require_tool(BLASTDBCMD, "blastdbcmd")
    db = ensure_pataa_db(db)
    WORK.mkdir(parents=True, exist_ok=True)
    log(f"dump: {db} -> {FASTA}")
    t0 = time.time()
    tmp = FASTA.with_name(FASTA.name + ".tmp")
    try:
        with open(tmp, "w") as out:
            r = subprocess.run(
                [str(BLASTDBCMD), "-db", str(db),
                 "-entry", "all", "-target_only"],
                stdout=out, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"blastdbcmd failed: {r.stderr[:300]}")
        os.replace(tmp, FASTA)
    finally:
        if tmp.exists():
            tmp.unlink()
    log(f"dump: done in {time.time()-t0:.1f}s ({FASTA.stat().st_size/1e9:.2f} GB)")


# ---------------------------------------------------------------------------
# stage: split
# ---------------------------------------------------------------------------

def stage_split():
    if (CHUNKS / ".done").exists() and INDEX_TSV.exists():
        log(f"split: exists ({len(list(CHUNKS.glob('chunk_*.fa')))} chunks), skip")
        return
    CHUNKS.mkdir(parents=True, exist_ok=True)
    log("split: assigning unique ids and writing chunks")
    t0 = time.time()
    idx = 0
    chunk_i = -1
    fh = None
    index = open(INDEX_TSV, "w")
    try:
        for defline, seq in iter_fasta(FASTA):
            if fh is None or (idx % CHUNK_SIZE == 0):
                if fh is not None:
                    fh.close()
                chunk_i += 1
                fh = open(CHUNKS / f"chunk_{chunk_i:05d}.fa", "w")
            qid = f"p{idx:09d}"
            fh.write(f">{qid}\n{seq}\n")
            index.write(f"{qid}\t{defline}\n")
            idx += 1
        if fh is not None:
            fh.close()
    finally:
        index.close()
    (CHUNKS / ".done").write_text(f"{idx}\n")
    log(f"split: {idx} sequences, {chunk_i+1} chunks in {time.time()-t0:.1f}s")


# ---------------------------------------------------------------------------
# stage: detect
# ---------------------------------------------------------------------------

def _detect_one(arg):
    """Run igblastp on one chunk; write outfmt7 + positives TSV."""
    chunk_str, igout_str, threshold = arg
    chunk = Path(chunk_str)
    igout = Path(igout_str)
    name = chunk.stem
    out7 = igout / f"{name}.ig7"
    pos = igout / f"{name}.pos.tsv"
    tmp7 = out7.with_name(out7.name + ".tmp")
    tmppos = pos.with_name(pos.name + ".tmp")
    env = dict(os.environ, IGDATA=str(IGDATA))
    try:
        with open(tmp7, "w") as o:
            r = subprocess.run(
                [str(IGBLASTP), "-query", str(chunk),
                 "-germline_db_V", str(GERMLINE_V), "-organism", "human",
                 "-outfmt", "7"],
                stdout=o, stderr=subprocess.DEVNULL, env=env)
        if r.returncode != 0:
            return (name, -1, 0)
        n_total, n_pos = 0, 0
        qid = None
        best = {}
        top = {}
        with open(tmp7) as fh, open(tmppos, "w") as ph:
            for line in fh:
                if line.startswith("# Query: "):
                    if qid is not None and best.get(qid, -1.0) >= threshold:
                        ph.write(f"{qid}\t{best[qid]:.1f}\t"
                                 f"{top[qid][0]}\t{top[qid][1]:.1f}\n")
                        n_pos += 1
                    qid = line[len("# Query: "):].strip()
                    n_total += 1
                    best[qid] = -1.0
                    top[qid] = ("", 0.0)
                elif line and line[0] == "V" and "\t" in line:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 14 and f[0] == "V":
                        try:
                            bs = float(f[-1]); pid = float(f[3])
                        except ValueError:
                            continue
                        if bs > best.get(qid, -1.0):
                            best[qid] = bs
                            top[qid] = (f[2], pid)
            if qid is not None and best.get(qid, -1.0) >= threshold:
                ph.write(f"{qid}\t{best[qid]:.1f}\t"
                         f"{top[qid][0]}\t{top[qid][1]:.1f}\n")
                n_pos += 1
        os.replace(tmp7, out7)
        os.replace(tmppos, pos)
        return (name, n_total, n_pos)
    finally:
        for tmp in (tmp7, tmppos):
            if tmp.exists():
                tmp.unlink()


def stage_detect(workers: int, threshold: float):
    _require_tool(IGBLASTP, "igblastp")
    _require_germline_v()
    IGOUT.mkdir(parents=True, exist_ok=True)
    chunks = sorted(CHUNKS.glob("chunk_*.fa"))
    if not chunks:
        raise RuntimeError("detect: no chunks; run 'split' first")
    todo = [c for c in chunks if not (IGOUT / f"{c.stem}.pos.tsv").exists()]
    log(f"detect: {len(todo)}/{len(chunks)} chunks pending, {workers} workers, "
        f"bitscore >= {threshold}")
    if not todo:
        return
    t0 = time.time()
    done = npos = ntot = 0
    failed = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_detect_one, (str(c), str(IGOUT), threshold)): c
                for c in todo}
        for f in as_completed(futs):
            name, tot, pos = f.result()
            done += 1
            if tot < 0:
                failed.append(name)
            else:
                ntot += tot
                npos += pos
            if done % 20 == 0 or done == len(todo):
                el = time.time() - t0
                log(f"  detect: {done}/{len(todo)} chunks, {ntot} queries, "
                    f"{npos} positives, {el:.0f}s")
    if failed:
        raise RuntimeError(
            f"detect: igblastp failed on {len(failed)} chunk(s): "
            f"{', '.join(failed[:5])}{' ...' if len(failed) > 5 else ''}")
    log(f"detect: done, {ntot} queries, {npos} positives in {time.time()-t0:.1f}s")


# ---------------------------------------------------------------------------
# stage: extract
# ---------------------------------------------------------------------------

def stage_extract():
    if POSITIVES_FASTA.exists() and POSITIVES_TSV.exists():
        log("extract: exists, skip")
        return
    chunks = {p.stem for p in CHUNKS.glob("chunk_*.fa")}
    detected = {p.name[:-len(".pos.tsv")] for p in IGOUT.glob("*.pos.tsv")}
    missing = chunks - detected
    if missing:
        raise RuntimeError(
            f"extract: {len(missing)} chunk(s) not detected; run 'detect' "
            f"first (e.g. {sorted(missing)[0]})")
    log("extract: collecting positives")
    positives = {}
    for p in sorted(IGOUT.glob("*.pos.tsv")):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    positives[f[0]] = (float(f[1]), f[2], float(f[3]))
    if not positives:
        raise RuntimeError(
            "extract: no positive sequences; check igblastp/germline DB "
            "and the detection threshold")
    tmp_tsv = POSITIVES_TSV.with_name(POSITIVES_TSV.name + ".tmp")
    tmp_fa = POSITIVES_FASTA.with_name(POSITIVES_FASTA.name + ".tmp")
    try:
        with open(tmp_tsv, "w") as fh:
            for qid, (bs, gl, pid) in positives.items():
                fh.write(f"{qid}\t{bs:.1f}\t{gl}\t{pid:.1f}\n")
        n = 0
        with open(tmp_fa, "w") as out:
            # positive ids are the synthetic ids used in the chunks
            for chunk in sorted(CHUNKS.glob("chunk_*.fa")):
                for defline, seq in iter_fasta(chunk):
                    qid = defline.split()[0]
                    if qid in positives:
                        out.write(f">{qid}\n{seq}\n")
                        n += 1
        os.replace(tmp_tsv, POSITIVES_TSV)
        os.replace(tmp_fa, POSITIVES_FASTA)
    finally:
        for tmp in (tmp_tsv, tmp_fa):
            if tmp.exists():
                tmp.unlink()
    log(f"extract: {len(positives)} positive ids, {n} sequences -> {POSITIVES_FASTA}")


# ---------------------------------------------------------------------------
# stage: number (ANARCI)
# ---------------------------------------------------------------------------

def _extract_cdrs_from_posmap(posmap: dict, prefix: str) -> dict:
    """Extract IMGT CDR1/CDR2/CDR3 in ANARCI sequence order.

    ``posmap`` must preserve ANARCI's insertion order (Python dicts keep
    insertion order, and ``json.loads`` restores it). IMGT insertion labels
    are **not** numerically sortable: CDR3 uses 111, 111A, 111B, ...
    then 112C, 112B, 112A, 112; CDR2 can use 60A..60D then 61E..61A; CDR1
    can use 32A..32C then 33D..33A. Walking positions numerically with
    ascending letters transposes these loops, so we iterate the posmap
    as-is and only use the number for region assignment.
    """
    cdrs = {cdr: "" for cdr in _IMGT_CDR_RANGES}
    for label, aa in posmap.items():
        if not aa or aa == "-" or not label.startswith(prefix):
            continue
        digits = "".join(c for c in label[1:] if c.isdigit())
        if not digits:
            continue
        pos_num = int(digits)
        for cdr, (start, end) in _IMGT_CDR_RANGES.items():
            if start <= pos_num <= end:
                cdrs[cdr] += aa
                break
    return cdrs


def _number_batch(batch):
    """batch: list[(qid, seq)] -> list[record dict]."""
    from anarci import anarci
    records = []
    for qid, seq in batch:
        try:
            res = anarci([(qid, seq)], scheme="imgt", output=False)
        except Exception:
            continue
        try:
            numbering = res[0][0][0][0]
            details = res[1][0][0] if res[1] and res[1][0] else {}
        except (IndexError, TypeError):
            continue
        if not numbering:
            continue
        chain_type = (details or {}).get("chain_type", "") or ""
        prefix = "H" if chain_type.upper() == "H" else "L"
        posmap = {}
        for (pos_num, ins_flag), aa in numbering:
            if aa == "-":
                continue
            label = (f"{prefix}{pos_num}{ins_flag}"
                     if ins_flag and ins_flag != " " else f"{prefix}{pos_num}")
            posmap[label] = aa.upper()
        records.append({
            "qid": qid,
            "chain": "VH" if prefix == "H" else "VL",
            "sequence": seq,
            "numbered": posmap,
            "cdrs": _extract_cdrs_from_posmap(posmap, prefix),
        })
    return records


def stage_number(workers: int):
    if CDR_RECORDS.exists() and CDR_RECORDS.stat().st_size > 0:
        log("number: exists, skip")
        return
    if not POSITIVES_FASTA.exists() or POSITIVES_FASTA.stat().st_size == 0:
        raise RuntimeError(
            "number: positives.fasta missing/empty; run 'extract' first")
    try:
        from anarci import anarci  # noqa: F401
    except Exception as exc:
        raise RuntimeError(
            f"number: ANARCI unavailable ({exc}); install with "
            "'pip install anarci' and HMMER (hmmscan/hmmsearch on PATH)")
    if not (shutil.which("hmmscan") or shutil.which("hmmsearch")):
        raise RuntimeError(
            "number: HMMER not found; install hmmer 3.x "
            "(hmmscan/hmmsearch must be on PATH)")
    log("number: running ANARCI")
    seqs = list(iter_fasta(POSITIVES_FASTA))
    log(f"number: {len(seqs)} positive sequences")
    batches = [seqs[i:i+NUMBER_BATCH] for i in range(0, len(seqs), NUMBER_BATCH)]
    t0 = time.time()
    done = 0
    nrec = 0
    tmp = CDR_RECORDS.with_name(CDR_RECORDS.name + ".tmp")
    try:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_number_batch, b) for b in batches]
            with open(tmp, "w") as out:
                for f in as_completed(futs):
                    recs = f.result()
                    for rec in recs:
                        out.write(json.dumps(rec) + "\n")
                    nrec += len(recs)
                    done += 1
                    if done % 20 == 0 or done == len(batches):
                        log(f"  number: {done}/{len(batches)} batches, "
                            f"{nrec} records, {time.time()-t0:.0f}s")
        if nrec == 0:
            raise RuntimeError("number: ANARCI produced no records")
        os.replace(tmp, CDR_RECORDS)
    finally:
        if tmp.exists():
            tmp.unlink()
    if nrec < len(seqs):
        log(f"number: WARNING {len(seqs)-nrec} sequences failed ANARCI")
    log(f"number: {nrec} CDR records in {time.time()-t0:.1f}s")


# ---------------------------------------------------------------------------
# stage: repair (re-extract CDRs from stored numbering, then build)
# ---------------------------------------------------------------------------

def stage_repair():
    """Re-extract CDRs from stored ``numbered`` posmaps, then rebuild.

    ``stage_build`` trusts the ``cdrs`` already stored in
    ``cdr_records.jsonl``, and ``stage_number`` skips when that file exists.
    After a CDR-extraction fix, run this stage instead of re-numbering:
    it streams the JSONL, recomputes ``cdrs`` from ``numbered`` in ANARCI
    sequence order, atomically replaces the JSONL, and calls ``stage_build``.
    No ANARCI or source BLAST database is required.
    """
    if not CDR_RECORDS.exists() or CDR_RECORDS.stat().st_size == 0:
        raise RuntimeError(
            "repair: cdr_records.jsonl missing/empty; run 'number' first")
    tmp = CDR_RECORDS.with_name(CDR_RECORDS.name + ".repair.tmp")
    log(f"repair: re-extracting CDRs from {CDR_RECORDS}")
    t0 = time.time()
    n = n_changed = n_no_numbered = 0
    try:
        with open(CDR_RECORDS) as src, open(tmp, "w") as out:
            for line in src:
                if not line.strip():
                    continue
                rec = json.loads(line)
                numbered = rec.get("numbered")
                if not numbered:
                    n_no_numbered += 1
                else:
                    prefix = "H" if rec.get("chain") == "VH" else "L"
                    new_cdrs = _extract_cdrs_from_posmap(numbered, prefix)
                    if new_cdrs != rec.get("cdrs"):
                        n_changed += 1
                    rec["cdrs"] = new_cdrs
                out.write(json.dumps(rec) + "\n")
                n += 1
                if n % 20000 == 0:
                    log(f"  repair: {n} records, {n_changed} changed, "
                        f"{time.time()-t0:.0f}s")
        os.replace(tmp, CDR_RECORDS)
    finally:
        if tmp.exists():
            tmp.unlink()
    log(f"repair: {n} records, {n_changed} CDR sets changed"
        + (f", {n_no_numbered} without numbered (kept as-is)"
           if n_no_numbered else "")
        + f", {time.time()-t0:.0f}s")
    stage_build()


# ---------------------------------------------------------------------------
# stage: build
# ---------------------------------------------------------------------------

def stage_build():
    if not CDR_RECORDS.exists() or CDR_RECORDS.stat().st_size == 0:
        raise RuntimeError(
            "build: cdr_records.jsonl missing/empty; run 'number' "
            "(or 'repair' for an existing JSONL)")
    log("build: merging")
    # metadata from detection
    meta = {}
    with open(POSITIVES_TSV) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                meta[f[0]] = {"bitscore": float(f[1]), "v_germline": f[2],
                              "v_identity": float(f[3])}
    # deflines for positives only
    deflines = {}
    need = set(meta)
    with open(INDEX_TSV) as fh:
        for line in fh:
            qid, _, dl = line.rstrip("\n").partition("\t")
            if qid in need:
                deflines[qid] = dl
    # cdr records
    records = []
    skipped_short = 0
    skipped_long = 0
    with open(CDR_RECORDS) as fh:
        for line in fh:
            rec = json.loads(line)
            qid = rec["qid"]
            m = meta.get(qid)
            if not m:
                continue
            if any(len(rec["cdrs"].get(c, "")) < MIN_CDR[c] for c in MIN_CDR):
                skipped_short += 1
                continue
            if len(rec["cdrs"].get("CDR3", "")) > MAX_CDR3:
                skipped_long += 1
                continue
            info = parse_defline(deflines.get(qid, ""))
            records.append({
                "qid": qid,
                "accession": info["accession"],
                "patents": info["patents"],
                "species": info["species"],
                "sequence": rec.get("sequence", ""),
                "numbered": rec["numbered"],
                "chain": rec["chain"],
                "cdrs": rec["cdrs"],
                **m,
            })
    if not records:
        raise RuntimeError("build: no records passed the CDR filters")
    # map qid -> original defline for dedup key by sequence
    # We don't retain raw sequence in records; reconstruct is unnecessary — dedup
    # by the CDR+numbering signature instead (equivalent for homology DB).
    seen = {}
    deduped = []
    for r in records:
        key = (r["chain"], r["cdrs"].get("CDR1", ""), r["cdrs"].get("CDR2", ""),
               r["cdrs"].get("CDR3", ""))
        if key in seen:
            prev = seen[key]
            for p in r["patents"]:
                if p not in prev["patents"] and len(prev["patents"]) < 20:
                    prev["patents"].append(p)
            continue
        seen[key] = r
        deduped.append(r)
    deduped.sort(key=lambda x: -x["bitscore"])
    # emit cdr_homology-compatible records
    out = []
    for r in deduped:
        chain_key = r["chain"]
        out.append({
            "name": r["accession"] or r["qid"],
            "source": "PATAA",
            "species": r["species"],
            "patents": r["patents"],
            "v_germline": r["v_germline"],
            "v_identity": r["v_identity"],
            "v_bitscore": r["bitscore"],
            "format": "",
            "genetics": "",
            "target": "",
            "vd_lc": chain_key,
            "highest_trial": "",
            "chains": {chain_key: {"sequence": r.get("sequence", ""), **r["cdrs"]}},
        })
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    tmp_json = OUT_JSON.with_name(OUT_JSON.name + ".tmp")
    tmp_stats = OUT_STATS.with_name(OUT_STATS.name + ".tmp")
    stats = {
        "database_name": "PATAA_CDR",
        "source": "GenBank Patent Division (PATAA) BLAST v5",
        "detection": "igblastp human germline, max V bitscore >= %.0f" % THRESHOLD,
        "numbering": "ANARCI IMGT",
        "total_positive_records": len(records),
        "skipped_short_cdr": skipped_short,
        "skipped_long_cdr3": skipped_long,
        "total_after_dedup": len(out),
        "chain_distribution": {
            c: sum(1 for r in out if r["vd_lc"] == c) for c in ("VH", "VL")},
        "cdr_length_summary": {},
        "top_germlines": {},
    }
    for ck in ("VH", "VL"):
        for cdr in ("CDR1", "CDR2", "CDR3"):
            lens = [len(r["chains"][ck][cdr]) for r in out
                    if r["vd_lc"] == ck and r["chains"][ck].get(cdr)]
            if lens:
                stats["cdr_length_summary"][f"{ck}_{cdr}"] = {
                    "count": len(lens), "min": min(lens), "max": max(lens),
                    "mean": round(sum(lens)/len(lens), 1)}
    gl = {}
    for r in out:
        g = r["v_germline"]
        if g:
            gl[g] = gl.get(g, 0) + 1
    stats["top_germlines"] = dict(sorted(gl.items(), key=lambda x: -x[1])[:30])
    try:
        with open(tmp_json, "w") as fh:
            json.dump(out, fh, ensure_ascii=False)
        with open(tmp_stats, "w") as fh:
            json.dump(stats, fh, indent=2, ensure_ascii=False)
        os.replace(tmp_json, OUT_JSON)
        os.replace(tmp_stats, OUT_STATS)
    finally:
        for tmp in (tmp_json, tmp_stats):
            if tmp.exists():
                tmp.unlink()
    log(f"build: {len(out)} entries -> {OUT_JSON}")
    log(f"build: stats -> {OUT_STATS}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", nargs="?", default="all",
                    choices=["all", "dump", "split", "detect", "extract",
                             "number", "repair", "build"])
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    args = ap.parse_args()

    globals()["THRESHOLD"] = args.threshold
    if args.stage in ("all", "dump"):
        stage_dump(Path(args.db))
    if args.stage in ("all", "split"):
        stage_split()
    if args.stage in ("all", "detect"):
        stage_detect(args.workers, args.threshold)
    if args.stage in ("all", "extract"):
        stage_extract()
    if args.stage in ("all", "number"):
        stage_number(args.workers)
    if args.stage == "repair":
        stage_repair()
    if args.stage in ("all", "build"):
        stage_build()


if __name__ == "__main__":
    main()
